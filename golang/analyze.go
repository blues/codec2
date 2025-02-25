package codec2

import (
	"fmt"
	"math"
)

// grid spacing for LSP root searches
const LSP_DELTA1 = 0.01

// dftSpeech implements the C dft_speech() function from sine.c.
// It builds a real input vector from the analysis buffer Sn using window w,
// computes the FFT via the FFT interface (stored in fftFwdCfg),
// and prints debug information.
// (Note: Do not swap halves—the raw FFT output already matches the C version.)
func dftSpeech(c2const *C2Const, fftFwdCfg interface{}, Sw []COMP, Sn []float64, w []float64) {
	mPitch := c2const.MPitch // full analysis window length (e.g., 320)
	nw := c2const.Nw         // nonzero window length (e.g., 279)
	in := make([]float64, FFTSize)
	// Zero initialize input.
	for i := 0; i < FFTSize; i++ {
		in[i] = 0.0
	}
	// Copy second half of windowed signal into the beginning of FFT input.
	// (Indices 0 .. nw/2 - 1)
	for i := 0; i < nw/2; i++ {
		in[i] = Sn[i+(mPitch/2)] * w[i+(mPitch/2)]
	}
	// Copy first half of windowed signal into the end of FFT input.
	// (Indices FFTSize - nw/2 .. FFTSize - 1)
	for i := 0; i < nw/2; i++ {
		in[FFTSize-nw/2+i] = Sn[i+(mPitch/2)-(nw/2)] * w[i+(mPitch/2)-(nw/2)]
	}
	// Compute FFT.
	fftInst := fftFwdCfg.(FFT)
	fftRes := fftInst.Forward(in)
	// Copy the raw FFT result into Sw.
	for i, c := range fftRes {
		Sw[i] = COMP{Real: real(c), Imag: imag(c)}
	}
}

// nlp implements the C nlp() function from nlp.c.
func nlp(nlpState interface{}, Sn []float64, n int, Sw []COMP, W []float64, prevF0 float64) (float64, float64) {
	nlpObj := nlpState.(*NLP)
	m := nlpObj.m
	// Square the latest n samples.
	for i := m - n; i < m; i++ {
		nlpObj.sq[i] = Sn[i] * Sn[i]
	}
	// Apply notch filter at DC.
	for i := m - n; i < m; i++ {
		notch := nlpObj.sq[i] - nlpObj.mem_x
		notch += COEFF * nlpObj.mem_y
		nlpObj.mem_x = nlpObj.sq[i]
		nlpObj.mem_y = notch
		nlpObj.sq[i] = notch + 1.0 // add a small constant as in C code
	}
	// FIR filtering (decimation FIR filter).
	for i := m - n; i < m; i++ {
		for j := 0; j < NLP_NTAP-1; j++ {
			nlpObj.mem_fir[j] = nlpObj.mem_fir[j+1]
		}
		nlpObj.mem_fir[NLP_NTAP-1] = nlpObj.sq[i]
		nlpObj.sq[i] = 0.0
		for j := 0; j < NLP_NTAP; j++ {
			nlpObj.sq[i] += nlpObj.mem_fir[j] * nlpFir[j]
		}
	}
	// Prepare Fw array.
	Fw := make([]COMP, PE_FFT_SIZE)
	for i := 0; i < PE_FFT_SIZE; i++ {
		Fw[i] = COMP{0, 0}
	}
	// Decimate and window.
	for i := 0; i < m/DEC; i++ {
		Fw[i].Real = nlpObj.sq[i*DEC] * nlpObj.w[i]
	}
	// FFT of decimated data.
	temp := make([]float64, PE_FFT_SIZE)
	for i := 0; i < PE_FFT_SIZE; i++ {
		temp[i] = Fw[i].Real
	}
	fftInst := nlpObj.fft_cfg.(FFT)
	fftRes := fftInst.Forward(temp)
	for i, c := range fftRes {
		// Compute magnitude squared.
		Fw[i].Real = real(c)*real(c) + imag(c)*imag(c)
		Fw[i].Imag = 0
	}
	// Global peak search.
	min_bin := int(float64(PE_FFT_SIZE*DEC) / float64(nlpObj.pmax))
	gmax := 0.0
	gmax_bin := min_bin
	for i := int(float64(PE_FFT_SIZE*DEC) / float64(nlpObj.pmax)); i <= int(float64(PE_FFT_SIZE*DEC)/float64(nlpObj.pmin)); i++ {
		if Fw[i].Real > gmax {
			gmax = Fw[i].Real
			gmax_bin = i
		}
	}
	bestF0 := postProcessSubMultiples(Fw, float64(nlpObj.pmin), float64(nlpObj.pmax), gmax, gmax_bin, &prevF0)
	// Shift the square buffer to make room for new samples.
	for i := 0; i < m-n; i++ {
		nlpObj.sq[i] = nlpObj.sq[i+n]
	}
	return float64(nlpObj.Fs) / bestF0, bestF0
}

// postProcessSubMultiples replicates the C post_process_sub_multiples() function.
func postProcessSubMultiples(Fw []COMP, pmin, pmax float64, gmax float64, gmax_bin int, prevF0 *float64) float64 {
	mult := 2
	min_bin := int(float64(PE_FFT_SIZE*DEC) / pmax)
	cmax_bin := gmax_bin
	prev_f0_bin := int(*prevF0 * float64(PE_FFT_SIZE*DEC) / float64(SampleRate))
	for float64(gmax_bin)/float64(mult) >= float64(min_bin) {
		b := gmax_bin / mult
		bmin := int(0.8 * float64(b))
		bmax := int(1.2 * float64(b))
		if bmin < min_bin {
			bmin = min_bin
		}
		var thresh float64
		if prev_f0_bin > bmin && prev_f0_bin < bmax {
			thresh = CNLP * 0.5 * gmax
		} else {
			thresh = CNLP * gmax
		}
		lmax := 0.0
		lmax_bin := bmin
		for b2 := bmin; b2 <= bmax; b2++ {
			if Fw[b2].Real > lmax {
				lmax = Fw[b2].Real
				lmax_bin = b2
			}
		}
		if lmax > thresh {
			if lmax_bin-1 >= 0 && lmax_bin+1 < len(Fw) {
				if Fw[lmax_bin-1].Real < lmax && Fw[lmax_bin+1].Real < lmax {
					cmax_bin = lmax_bin
				}
			}
		}
		mult++
	}
	bestF0 := float64(cmax_bin) * float64(SampleRate) / float64(PE_FFT_SIZE*DEC)
	return bestF0
}

// hsPitchRefinement implements harmonic sum pitch refinement.
func hsPitchRefinement(model *Model, Sw []COMP, pmin, pmax, pstep float64) {
	var E, Emax float64
	var b int
	Emax = 0.0
	Wom := model.Wo
	r := TWO_PI / FFTSize
	oneOnR := 1.0 / r
	for p := pmin; p <= pmax; p += pstep {
		E = 0.0
		Wo := TWO_PI / p
		var bFloat float64 = Wo * oneOnR
		var currentBFloat float64 = bFloat
		for m := 1; m <= model.L; m++ {
			b = int(currentBFloat + 0.5)
			E += Sw[b].Real*Sw[b].Real + Sw[b].Imag*Sw[b].Imag
			currentBFloat += bFloat
		}
		if E > Emax {
			Emax = E
			Wom = Wo
		}
	}
	model.Wo = Wom
}

// twoStagePitchRefinement transliterates the C two_stage_pitch_refinement() function.
func twoStagePitchRefinement(c2const *C2Const, model *Model, Sw []COMP) {
	var pminVal, pmaxVal, pstep float64
	pmaxVal = TWO_PI/model.Wo + 5
	pminVal = TWO_PI/model.Wo - 5
	pstep = 1.0
	hsPitchRefinement(model, Sw, pminVal, pmaxVal, pstep)
	pmaxVal = TWO_PI/model.Wo + 1
	pminVal = TWO_PI/model.Wo - 1
	pstep = 0.25
	hsPitchRefinement(model, Sw, pminVal, pmaxVal, pstep)
	if model.Wo < TWO_PI/float64(c2const.PMax) {
		model.Wo = TWO_PI / float64(c2const.PMax)
	}
	if model.Wo > TWO_PI/float64(c2const.PMin) {
		model.Wo = TWO_PI / float64(c2const.PMin)
	}
	model.L = int(math.Floor(math.Pi / model.Wo))
	if model.Wo*float64(model.L) >= 0.95*math.Pi {
		model.L--
	}
}

// estimateAmplitudes transliterates estimate_amplitudes() from nlp.c.
func estimateAmplitudes(model *Model, Sw []COMP, W []float64, estPhase int) {
	var i, m, am, bm int
	var den float64
	r := TWO_PI / FFTSize
	oneOnR := 1.0 / r
	for m = 1; m <= model.L; m++ {
		den = 0.0
		am = int(((float64(m) - 0.5) * model.Wo * oneOnR) + 0.5)
		bm = int(((float64(m) + 0.5) * model.Wo * oneOnR) + 0.5)
		for i = am; i < bm; i++ {
			den += Sw[i].Real*Sw[i].Real + Sw[i].Imag*Sw[i].Imag
		}
		model.A[m] = math.Sqrt(den)
		if estPhase != 0 {
			b := int((float64(m) * model.Wo / r) + 0.5)
			if b < 0 {
				b = 0
			}
			if b >= len(Sw) {
				b = len(Sw) - 1
			}
			model.Phi[m] = math.Atan2(Sw[b].Imag, Sw[b].Real)
		}
	}
}

// estVoicingMbe transliterates est_voicing_mbe() from nlp.c.
func estVoicingMbe(c2const *C2Const, model *Model, Sw []COMP, W []float64) float64 {
	var l, al, bl, m int
	var Am COMP
	var offset int
	var den float64
	errorAcc := 1e-4
	sig := 1e-4
	for l = 1; l <= int(float64(model.L)*1000.0/(float64(c2const.Fs)/2)); l++ {
		sig += model.A[l] * model.A[l]
	}
	Wo := model.Wo
	for l = 1; l <= int(float64(model.L)*1000.0/(float64(c2const.Fs)/2)); l++ {
		Am = COMP{0, 0}
		den = 0.0
		al = int(math.Ceil((float64(l)-0.5)*Wo*float64(FFTSize)/TWO_PI + 0.5))
		bl = int(math.Ceil((float64(l)+0.5)*Wo*float64(FFTSize)/TWO_PI + 0.5))
		offset = int(float64(FFTSize)/2 - float64(l)*Wo*float64(FFTSize)/TWO_PI + 0.5)
		for m = al; m < bl; m++ {
			Am.Real += Sw[m].Real * W[offset+m]
			Am.Imag += Sw[m].Imag * W[offset+m]
			den += W[offset+m] * W[offset+m]
		}
		Am.Real /= den
		Am.Imag /= den
		for m = al; m < bl; m++ {
			EwReal := Sw[m].Real - Am.Real*W[offset+m]
			EwImag := Sw[m].Imag - Am.Imag*W[offset+m]
			errorAcc += EwReal*EwReal + EwImag*EwImag
		}
	}
	snr := 10.0 * math.Log10(sig/errorAcc)
	if snr > V_THRESH {
		model.Voiced = true
	} else {
		model.Voiced = false
	}
	l2000 := int(float64(model.L) * 2000.0 / (float64(c2const.Fs) / 2))
	l4000 := int(float64(model.L) * 4000.0 / (float64(c2const.Fs) / 2))
	elow := 1e-4
	ehigh := 1e-4
	for l = 1; l <= l2000; l++ {
		elow += model.A[l] * model.A[l]
	}
	for l = l2000; l <= l4000; l++ {
		ehigh += model.A[l] * model.A[l]
	}
	eratio := 10.0 * math.Log10(elow/ehigh)
	if !model.Voiced && eratio > 10.0 {
		model.Voiced = true
	}
	if model.Voiced {
		if eratio < -10.0 {
			model.Voiced = false
		}
		sixty := 60.0 * TWO_PI / float64(c2const.Fs)
		if eratio < -4.0 && model.Wo <= sixty {
			model.Voiced = false
		}
	}
	return snr
}

// speechToUQLSPS is a direct transliteration of the C speech_to_uq_lsps() function.
// It computes the windowed signal, its energy, then performs autocorrelation,
// applies Levinson–Durbin to obtain LPC coefficients, applies bandwidth expansion,
// converts the LPCs to LSPs via LpcToLsp, and returns the LPC energy E along with
// the computed LSP frequencies and LPC coefficients.
// Parameters:
//
//	sn:     input analysis buffer (length mPitch)
//	w:      analysis window (length mPitch)
//	mPitch: length of the analysis buffer (e.g. 320)
//	order:  LPC order (e.g. LpcOrder, typically 10)
//
// Returns:
//
//	energy: computed LPC energy E,
//	lsp:    computed LSP frequencies (in radians),
//	ak:     LPC coefficients (length order+1, with ak[0]==1)
func speechToUQLSPS(sn, w []float64, mPitch, order int) (energy float64, lsp []float64, ak []float64) {

	// Compute windowed signal Wn and energy.
	Wn := make([]float64, mPitch)
	energy = 0.0
	for i := 0; i < mPitch; i++ {
		Wn[i] = sn[i] * w[i]
		energy += Wn[i] * Wn[i]
	}

	// Trap the zero-energy case.
	if energy == 0.0 {
		lsp = make([]float64, order)
		for i := 0; i < order; i++ {
			lsp[i] = (math.Pi / float64(order)) * float64(i)
		}
		ak = make([]float64, order+1)
		ak[0] = 1.0
		return energy, lsp, ak
	}

	// Compute autocorrelation R[0..order] of Wn.
	R := make([]float64, order+1)
	for i := 0; i <= order; i++ {
		R[i] = 0.0
		for j := 0; j < mPitch-i; j++ {
			R[i] += Wn[j] * Wn[j+i]
		}
	}

	// Compute LPC coefficients via Levinson–Durbin.
	ak = levinsonDurbin(R, order)

	// Compute LPC energy: E = sum_{i=0}^{order} ak[i]*R[i].
	E := 0.0
	for i := 0; i <= order; i++ {
		E += ak[i] * R[i]
	}

	// Apply bandwidth expansion: for i=0..order, ak[i] *= (0.994)^i.
	for i := 0; i <= order; i++ {
		ak[i] *= math.Pow(0.994, float64(i))
	}

	lsp = make([]float64, order)

	// Convert LPCs to LSPs.
	roots := LpcToLsp(ak, lsp, order, 5, LSP_DELTA1)
	if roots != order {
		// If root finding fails, use benign (evenly spaced) LSP values.
		for i := 0; i < order; i++ {
			lsp[i] = (math.Pi / float64(order)) * float64(i)
		}
	}

	return E, lsp, ak
}

func levinsonDurbin(R []float64, order int) []float64 {
	// Allocate ak with length order+1; ak[0] is 1.0.
	ak := make([]float64, order+1)
	err := R[0]
	ak[0] = 1.0
	for i := 1; i <= order; i++ {
		sum := 0.0
		for j := 1; j < i; j++ {
			sum += ak[j] * R[i-j]
		}
		// Reflection coefficient with negative sign.
		k := -(R[i] + sum) / err
		// Make a temporary copy of the current coefficients.
		prevA := make([]float64, i+1)
		copy(prevA, ak[:i+1])
		for j := 1; j < i; j++ {
			ak[j] = prevA[j] + k*prevA[i-j]
		}
		ak[i] = k
		err *= (1 - k*k)
	}
	return ak

}

// analyzeOneFrame replicates the C function analyse_one_frame().
// It only updates the analysis buffer, computes the FFT,
// estimates pitch (setting Wo and L), refines the model parameters,
// and then estimates amplitudes and voicing.
func (c *Codec2) analyzeOneFrame(speech []float64, model *Model) error {
	// Get local variables from c2const.
	nSamp := c.c2const.NSamp   // e.g., 80 samples
	mPitch := c.c2const.MPitch // e.g., 320 samples
	if len(c.Sn) < mPitch {
		return fmt.Errorf("analysis buffer length (%d) is less than MPitch (%d)", len(c.Sn), mPitch)
	}

	// --- Buffer Management ---
	// Shift older samples in c.Sn.
	for i := 0; i < mPitch-nSamp; i++ {
		c.Sn[i] = c.Sn[i+nSamp]
	}
	// Copy the new speech samples into the tail.
	for i := 0; i < nSamp; i++ {
		c.Sn[i+mPitch-nSamp] = speech[i]
	}

	// --- FFT Computation ---
	// Allocate an FFT output array of length FFTSize.
	Sw := make([]COMP, FFTSize)
	// Call our dftSpeech function (analogous to dft_speech in C).
	dftSpeech(&c.c2const, c.fftFwdCfg, Sw, c.Sn, c.w)

	// --- Pitch Estimation ---
	var pitch float64
	pitch, c.prevF0Enc = nlp(c.nlp, c.Sn, nSamp, Sw, c.W, c.prevF0Enc)
	model.Wo = TWO_PI / pitch
	model.L = int(math.Floor(math.Pi / model.Wo))

	// --- Model Parameter Refinement ---
	twoStagePitchRefinement(&c.c2const, model, Sw)

	// --- Amplitude and Voicing Estimation ---
	flag := 0
	if c.fmlfeat != nil {
		flag = 1
	}
	estimateAmplitudes(model, Sw, c.W, flag)
	estVoicingMbe(&c.c2const, model, Sw, c.W)

	// (Note: There is no energy or LSP conversion here, as in the C version.)
	return nil
}

// makeAnalysisWindow computes the analysis window exactly like the C code.
// It allocates a window of length MPitch. Only the central nw samples (where nw is
// defined in C2Const) are nonzero and are computed using the Hamming formula.
// (You may adjust the normalization factor below so that the sum of the window
// matches the C version.)
func makeAnalysisWindow(c2const *C2Const) []float64 {
	mPitch := c2const.MPitch
	nw := c2const.Nw
	w := make([]float64, mPitch)
	// Zero out the parts outside the central window.
	for i := 0; i < mPitch/2-nw/2; i++ {
		w[i] = 0.0
	}
	start := mPitch/2 - nw/2
	sum := 0.0
	for i := start; i < start+nw; i++ {
		j := i - start
		w[i] = 0.5 - 0.5*math.Cos(2*math.Pi*float64(j)/float64(nw-1))
		sum += w[i] * w[i]
	}
	for i := start + nw; i < mPitch; i++ {
		w[i] = 0.0
	}
	// Compute normalization factor. In C: norm = 1/sqrt(sum*FFT_ENC)
	// You might adjust this factor if needed.
	norm := 1.0 / math.Sqrt(sum*float64(FFTSize))
	for i := 0; i < mPitch; i++ {
		w[i] *= norm
	}
	return w
}
