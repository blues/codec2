package codec2

import (
	"math"
)

// synthesizeOneFrame synthesizes a subframe of speech from the given model,
// performing phase synthesis, post filtering, inverse FFT synthesis, gain adjustment,
// and ear protection. This is a direct transliteration of the C function
// synthesise_one_frame() from codec2.c.
// The output speech is returned as a slice of int16.
func (c *Codec2) synthesizeOneFrame(model *Model, speech []int16, Aw []COMP, gain float64) {

	// number of samples per subframe (e.g., 80)
	nsam := c.nsam

	// Phase synthesis: if in 700C mode, use Aw directly; otherwise, use a temporary H.
	// (The C code distinguishes between 700C and others; here we assume 700C is handled elsewhere.)
	H := make([]COMP, MAX_AMP+1)
	samplePhase(model, H, Aw)
	phaseSynthZeroOrder(nsam, model, &c.ex_phase, H)

	// Apply postfilter to update model parameters and background noise estimate.
	postfilter(model, &c.bg_est)

	// Synthesize the time-domain speech: inverse FFT synthesis.
	synthesise(nsam, c.fftInvCfg, c.Sn_, model, c.Pn, 1)

	// Multiply synthesized signal by the gain.
	for i := 0; i < nsam; i++ {
		c.Sn_[i] *= gain
	}

	// Apply ear protection.
	earProtection(c.Sn_, nsam)

	// Convert synthesized float samples to int16 PCM output.
	for i := 0; i < nsam; i++ {
		s := c.Sn_[i]
		if s > 32767.0 {
			speech[i] = 32767
		} else if s < -32767.0 {
			speech[i] = -32767
		} else {
			speech[i] = int16(s)
		}
	}
}

// synthesise is a faithful transliteration of the C function synthesise().
// It reconstructs the full spectrum from the half-spectrum, calls the inverse FFT,
// and then performs overlap–add with the synthesis window Pn.
// Parameters:
//
//	nSamp: number of samples per subframe (e.g., 80)
//	fftInvCfg: the inverse FFT configuration (implements FFT interface)
//	Sn_: synthesis buffer (must be at least 2*nSamp in length)
//	model: model parameters (which include Wo, L, A[] and Phi[])
//	Pn: synthesis (trapezoidal) window (length = 2*nSamp)
//	shift: if nonzero, the new frame overwrites Sn_; if zero, it is added.
func synthesise(nSamp int, fftInvCfg FFT, Sn_ []float64, model *Model, Pn []float64, shift int) {
	N := FFTSize             // FFTSize (e.g., 512)
	half := N/2 + 1          // length of half-spectrum
	Sw := make([]COMP, half) // allocate half-spectrum
	// Update memories
	if shift != 0 {
		for i := 0; i < nSamp-1; i++ {
			Sn_[i] = Sn_[i+nSamp]
		}
		Sn_[nSamp-1] = 0.0
	}
	// Zero initialize Sw.
	for i := 0; i < half; i++ {
		Sw[i] = COMP{Real: 0.0, Imag: 0.0}
	}
	// Compute r = TWO_PI / FFTSize.
	r := TWO_PI / float64(N)

	// For each harmonic m = 1 to model.L, determine the corresponding DFT bin.
	// (Note: in C, b = (int)(m*model->Wo/r + 0.5);)
	for m := 1; m <= model.L; m++ {
		bin := int(float64(m)*model.Wo/r + 0.5)
		if bin >= half {
			bin = half - 1
		}
		// In the C code, the synthesis filter uses the conjugate of A[b].
		// Here, we assume that the analysis filter provided in Aw has been computed
		// so that we can re-create a harmonic with amplitude model.A[m] and phase model.Phi[m].
		Sw[bin] = COMP{
			Real: model.A[m] * math.Cos(model.Phi[m]),
			Imag: model.A[m] * math.Sin(model.Phi[m]),
		}
	}

	// Reconstruct the full complex spectrum.
	fullSpectrum := make([]complex128, N)
	// DC and Nyquist terms.
	fullSpectrum[0] = complex(Sw[0].Real, Sw[0].Imag)
	fullSpectrum[N/2] = complex(Sw[N/2].Real, Sw[N/2].Imag)
	// For k = 1 to N/2 - 1, set X[N-k] = conjugate(X[k]).
	for k := 1; k < N/2; k++ {
		fullSpectrum[k] = complex(Sw[k].Real, Sw[k].Imag)
		fullSpectrum[N-k] = complex(Sw[k].Real, -Sw[k].Imag)
	}

	// Call the inverse FFT on the full spectrum.
	// The Inverse() method should now return a slice of length N.
	sw := fftInvCfg.Inverse(fullSpectrum)

	// In our FFT implementation, the IFFT scales by 1/FFTSize.
	// To match the C version (which applies no scaling), multiply by FFTSize.
	scaleFactor := float64(N)
	for i := range sw {
		sw[i] *= scaleFactor
	}

	// Overlap-add:
	// For i = 0 to nSamp-1, add the tail of the IFFT result (indices N - nSamp ... N-1)
	// multiplied by the first half of the synthesis window Pn[0...nSamp-1] to Sn_.
	var FFTI_FACTOR = 1.0
	if false { // if NOT USE_KISS_FFT in C
		FFTI_FACTOR = float64(N)
	}
	for i := 0; i < nSamp-1; i++ {
		Sn_[i] += sw[N-nSamp+1+i] * Pn[i] * FFTI_FACTOR
	}

	// For the remaining nSamp samples (i = nSamp-1 to 2*nSamp-1):
	if shift != 0 {
		j := 0
		for i := nSamp - 1; i < 2*nSamp; i++ {
			if j >= len(sw) {
				break
			}
			Sn_[i] = sw[j] * Pn[i] * FFTI_FACTOR
			j++
		}
	} else {
		j := 0
		for i := nSamp - 1; i < 2*nSamp; i++ {
			if j >= len(sw) {
				break
			}
			Sn_[i] += sw[j] * Pn[i] * FFTI_FACTOR
			j++
		}
	}
}

// makeSynthesisWindow generates the trapezoidal (Parzen) synthesis window.
// This is a direct transliteration of the C function make_synthesis_window().
// It returns a slice of length 2*NSamp (i.e. 20ms worth of samples) where NSamp is
// the number of samples per 10ms frame.
func makeSynthesisWindow(c2const *C2Const) []float64 {
	// n_samp is the number of samples per 10ms subframe.
	nsamp := c2const.NSamp
	// tw is the trapezoidal overlap (in samples)
	tw := c2const.Tw
	total := 2 * nsamp // synthesis window length

	Pn := make([]float64, total)

	// First segment: indices 0 to (nsamp/2 - tw - 1)
	for i := 0; i < nsamp/2-tw; i++ {
		Pn[i] = 0.0
	}

	// Second segment: indices nsamp/2 - tw to nsamp/2 + tw - 1
	win := 0.0
	increment := 1.0 / (2.0 * float64(tw))
	for i := nsamp/2 - tw; i < nsamp/2+tw; i++ {
		Pn[i] = win
		win += increment
	}

	// Third segment: indices nsamp/2 + tw to 3*nsamp/2 - tw - 1
	for i := nsamp/2 + tw; i < 3*nsamp/2-tw; i++ {
		Pn[i] = 1.0
	}

	// Fourth segment: indices 3*nsamp/2 - tw to 3*nsamp/2 + tw - 1
	win = 1.0
	decrement := 1.0 / (2.0 * float64(tw))
	for i := 3*nsamp/2 - tw; i < 3*nsamp/2+tw; i++ {
		Pn[i] = win
		win -= decrement
	}

	// Fifth segment: indices 3*nsamp/2 + tw to 2*nsamp - 1
	for i := 3*nsamp/2 + tw; i < total; i++ {
		Pn[i] = 0.0
	}

	return Pn
}
