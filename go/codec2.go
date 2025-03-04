package codec2

import (
	"fmt"
	"math"
)

// NewCodec2 creates a new Codec2 instance for 2400bps mode.  We only support
// 8000 and 16000 sample rates.
func NewCodec2() (*Codec2, error) {
	Nw := 279 // For 8000 Hz.
	if SampleRate != 8000 {
		Nw = 511 // actually a bit shorter in time but lets us maintain constant FFT size
	}
	c := &Codec2{
		mode: Mode2400,
		nsam: SamplesPerSubFrame,
		nbit: BitsPerFrame,
		c2const: C2Const{
			Fs:     SampleRate,
			NSamp:  SamplesPerSubFrame,
			MaxAmp: int(math.Floor(float64(SampleRate) * PMaxS / 2)),
			PMin:   int(math.Floor(float64(SampleRate) * PMinS)),
			PMax:   int(math.Floor(float64(SampleRate) * PMaxS)),
			MPitch: int(math.Floor(float64(SampleRate) * MPitchS)),
			WoMin:  TWO_PI / math.Floor(float64(SampleRate)*PMaxS),
			WoMax:  TWO_PI / math.Floor(float64(SampleRate)*PMinS),
			Nw:     Nw,
			Tw:     int(float64(SampleRate) * TWS),
		},
		prevF0Enc: 1 / PMaxS,
	}

	// Allocate Hamming window for analysis: length = MPitch.
	c.w = makeAnalysisWindow(&c.c2const)

	// Allocate analysis buffer (length = MPitch).
	c.Sn = make([]float64, c.c2const.MPitch)
	for i := 0; i < c.c2const.MPitch; i++ {
		c.Sn[i] = 1.0
	}
	// Allocate synthesis buffer (length = SamplesPerFrame).
	c.Sn_ = make([]float64, SamplesPerFrame)

	// Initialize FFT configuration using your FFT package.
	c.fftFwdCfg = NewFFT(FFTSize)

	// Allocate and (ideally) compute the frequency-domain window vector.
	// (For now, we fill W with 1.0.)
	c.W = make([]float64, FFTSize)
	for i := 0; i < FFTSize; i++ {
		c.W[i] = 1.0
	}

	// Remove duplicate allocation of c.w here.
	// (Previously, a new window of length SamplesPerSubFrame was allocated and overwrote the correct window.)

	// Initialize previous LSPs.
	c.prevLsps = make([]float64, LpcOrder)
	for i := 0; i < LpcOrder; i++ {
		c.prevLsps[i] = float64(i) * math.Pi / float64(LpcOrder+1)
	}

	// Initialize current model.
	c.model = NewModel()
	c.model.Wo = math.Pi / 100.0
	c.model.L = int(math.Floor(math.Pi / (math.Pi / 100.0)))

	// Initialize quantizer state.
	c.xq_enc = []float64{0.0, 0.0}
	c.xq_dec = []float64{0.0, 0.0}

	// Set previous model state.
	c.prev_model = NewModel()
	c.prev_model.Wo = TWO_PI / float64(c.c2const.PMax)
	c.prev_model.L = int(PI / c.prev_model.Wo)
	c.prev_model.Voiced = false
	c.prevE = 1

	// Set post filter parameters.
	c.lpc_pf = true
	c.bass_boost = true
	c.beta = 0.2
	c.gamma = 0.5

	c.prev_Wo = c.model.Wo

	// Initialize the NLP state.
	// This function mimics the C nlp_create() routine.
	c.nlp = nlpCreate(&c.c2const)

	// Initialize the inverse FFT configuration.
	c.fftInvCfg = NewFFT(FFTSize) // same FFT size as forward

	// Initialize synthesis buffer Pn: allocate 2*SamplesPerSubFrame samples.
	c.Pn = makeSynthesisWindow(&c.c2const)

	// Initialize the excitation phase and background estimate.
	c.ex_phase = 0.0
	c.bg_est = 0.0

	return c, nil
}

// NewCodecModels creates a new Codec2 instance for 2400bps mode.
func NewModel() (model Model) {
	model = Model{
		Voiced: true,
		E:      1.0,
		A:      make([]float64, MAX_AMP+1),
		Phi:    make([]float64, MAX_AMP+1),
		LSPs:   make([]float64, LpcOrder),
	}
	for i := 0; i < LpcOrder; i++ {
		model.LSPs[i] = float64(i+1) * math.Pi / float64(LpcOrder+1)
	}
	return
}

// nlpCreate creates and initializes an NLP state based on the provided C2Const.
// It mimics the behavior of the C nlp_create() function.
func nlpCreate(c2const *C2Const) *NLP {
	m := c2const.MPitch
	// Allocate NLP state.
	nlpObj := &NLP{
		Fs:      c2const.Fs,
		m:       m,
		sq:      make([]float64, m),
		mem_x:   0.0,
		mem_y:   0.0,
		mem_fir: make([]float64, NLP_NTAP),
		fft_cfg: NewFFT(PE_FFT_SIZE), // Using our FFT package.
		pmin:    c2const.PMin,
		pmax:    c2const.PMax,
	}
	// Allocate decimated analysis window of length m/DEC.
	n := m / DEC
	nlpObj.w = make([]float64, n)
	for i := 0; i < n; i++ {
		// The C code uses: 0.5 - 0.5*cos(2*pi*i/(n-1))
		nlpObj.w[i] = 0.5 - 0.5*math.Cos(2*math.Pi*float64(i)/float64(n-1))
	}
	return nlpObj
}

// Encode is a direct transliteration of codec2_encode_2400() from codec2.c.
func (c *Codec2) Encode(pcm []int16) ([]byte, error) {
	if len(pcm) != SamplesPerFrame {
		return nil, fmt.Errorf("invalid PCM frame length: got %d, want %d", len(pcm), SamplesPerFrame)
	}

	// Convert PCM to float64.
	speech := make([]float64, SamplesPerFrame)
	for i := 0; i < SamplesPerFrame; i++ {
		speech[i] = float64(pcm[i])
	}

	// Create output bit buffer.
	bits := make([]byte, BytesPerFrame)
	var nbit uint = 0

	var model Model
	model.A = make([]float64, MAX_AMP+1)
	model.Phi = make([]float64, MAX_AMP+1)
	model.LSPs = make([]float64, LpcOrder)

	// --- First 10ms analysis frame: only voicing is needed.
	// This calls analyse_one_frame, which does pitch, amplitudes, voicing.
	err := c.analyzeOneFrame(speech[0:SamplesPerSubFrame], &model)
	if err != nil {
		return nil, err
	}
	packBits(bits, &nbit, boolToInt(model.Voiced), 1)

	// --- Second 10ms analysis frame: full analysis.
	err = c.analyzeOneFrame(speech[SamplesPerSubFrame:SamplesPerFrame], &model)
	if err != nil {
		return nil, err
	}
	packBits(bits, &nbit, boolToInt(model.Voiced), 1)

	// --- Joint quantisation of Wo and energy.
	// Use the analysis buffer (c.Sn) and window (c.w) along with mPitch and LPC order.
	energy, lsps, _ := speechToUQLSPS(c.Sn, c.w, c.c2const.MPitch, LpcOrder)
	model.E = energy
	woE := encodeWoE(&model, model.E, c.xq_enc)
	packBits(bits, &nbit, woE, uint(WoEBits))

	// --- Encode LSPs using scalar quantisation.
	lspIndexes := make([]int, LSPScalarIndexes)
	encodeLSPScalar(lspIndexes, lsps, LpcOrder)
	for i := 0; i < LSPScalarIndexes; i++ {
		packBits(bits, &nbit, lspIndexes[i], uint(lspBits(i)))
	}

	// Pack spare bits (2 bits set to 0).
	packBits(bits, &nbit, 0, 2)

	if nbit != uint(BitsPerFrame) {
		return nil, fmt.Errorf("bit packing error: got %d bits, expected %d", nbit, BitsPerFrame)
	}

	return bits, nil
}

// Decode converts a 48-bit frame (6 bytes) into a 160-sample PCM frame.
// Decode decodes a 2400bps frame (48 bits = 6 bytes) into a 20ms (160-sample) PCM frame.
// This function is a direct replacement for the original Decode() in codec2.go.
func (c *Codec2) Decode(bits []byte) ([]int16, error) {

	// Verify that the input bit slice is exactly 6 bytes.
	if len(bits) != BytesPerFrame {
		return nil, fmt.Errorf("invalid bits buffer length: got %d, want %d", len(bits), BytesPerFrame)
	}
	var nbit uint = 0

	// --- Unpack the bitstream ---
	// Unpack voicing bits for the two 10ms subframes.
	v0, err := unpackBits(bits, &nbit, 1)
	if err != nil {
		return nil, err
	}
	v1, err := unpackBits(bits, &nbit, 1)
	if err != nil {
		return nil, err
	}

	// Unpack the joint Wo/E index (using WoEBits bits).
	woE, err := unpackBits(bits, &nbit, uint(WoEBits))
	if err != nil {
		return nil, fmt.Errorf("failed to unpack Wo/E index: %v", err)
	}

	// --- Build models for the two 10ms subframes ---
	// Create two MODEL objects (for the two 10ms frames).
	var model [2]Model
	model[0] = NewModel()
	model[1] = NewModel()

	// The voicing flags come directly from the unpacked bits.
	model[0].Voiced = (v0 != 0)
	model[1].Voiced = (v1 != 0)

	// Decode the joint Wo/E index into model[1].
	// This call updates model[1].E and the decoder state c.xq_dec.
	decodeWoE(&c.c2const, &model[1], &model[1].E, c.xq_dec, woE)

	// Unpack the LSP scalar indexes.
	var lsps [2][]float64
	lsps[0] = make([]float64, LpcOrder)
	lsps[1] = make([]float64, LpcOrder)
	lspIndexes := make([]int, LSPScalarIndexes)
	for i := 0; i < LSPScalarIndexes; i++ {
		idx, err := unpackBits(bits, &nbit, uint(lspBits(i)))
		if err != nil {
			return nil, fmt.Errorf("failed to unpack LSP index %d: %v", i, err)
		}
		lspIndexes[i] = idx
	}

	// Decode the LSPs for model[1] from the scalar indexes.
	decodeLSPScalar(lsps[1], lspIndexes, LpcOrder)
	// Check and, if necessary, adjust LSP order.
	checkLspOrder(lsps[1], LpcOrder)
	// Apply bandwidth expansion to ensure a minimum separation.
	bwExpandLsps(lsps[1], LpcOrder, 50.0, 100.0)

	// --- Interpolation ---
	// Interpolate parameters between the previous frame (stored in c.prev_model, c.prevE, c.prevLsps)
	// and the current model[1].

	// Interpolate fundamental frequency (Wo).
	interp_Wo(&model[0], &c.prev_model, &model[1], c.c2const.WoMin)
	// Interpolate energy.
	model[0].E = interp_energy(c.prevE, model[1].E)
	// Interpolate LSPs (using a weight of 0.5).
	interpolate_lsp_ver2(lsps[0], c.prevLsps, lsps[1], 0.5, LpcOrder)
	// (Note: The voicing flag for the interpolated frame is implicitly voiced.)

	// --- LPC Conversion and Synthesis ---
	var ak [2][]float64
	ak[0] = make([]float64, LpcOrder+1)
	ak[1] = make([]float64, LpcOrder+1)
	speech := make([]int16, SamplesPerFrame)
	for i := 0; i < 2; i++ {
		// Convert LSPs to LPC coefficients for each 10ms subframe.
		lspToLpc(lsps[i], ak[i], LpcOrder)
		// Allocate a buffer for the synthesis filter output (frequency domain).
		Aw := make([]COMP, FFTSize)
		var snr float64
		// Synthesize the first 10ms subframe from the interpolated model.
		aksToM2(c.fftFwdCfg, ak[i], LpcOrder, &model[i], model[i].E, &snr, 0, 0,
			c.lpc_pf, c.bass_boost, c.beta, c.gamma, Aw)
		applyLpcCorrection(&model[i])
		c.synthesizeOneFrame(&model[i], speech[c.nsam*i:], Aw, 1.0)
	}

	// --- Update state for next frame ---
	c.prev_model = model[1] // Save model[1] as previous frame model.
	c.prevE = model[1].E    // Save energy.
	for i := 0; i < LpcOrder; i++ {
		c.prevLsps[i] = lsps[1][i]
	}

	return speech, nil
}
