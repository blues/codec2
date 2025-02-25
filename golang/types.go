package codec2

import "math"

// Basic constants.
const (
	PI     = math.Pi
	TWO_PI = 2.0 * math.Pi
)

// FFT and NLP related constants.
const (
	FFTSize     = 512
	PE_FFT_SIZE = 512  // DFT size for pitch estimation (from nlp.c)
	DEC         = 5    // Decimation factor (from nlp.c)
	CNLP        = 0.3  // Post processor constant (from nlp.c)
	V_THRESH    = 6.0  // Voicing threshold (from nlp.c)
	NLP_NTAP    = 48   // Decimation FIR filter order (from nlp.c)
	COEFF       = 0.95 // Notch filter parameter (from nlp.c)
	MAX_AMP     = 160  // Maximum number of harmonics (from nlp.c)
)

// nlpFir: 48-tap 600Hz low-pass FIR filter coefficients (from nlp.c)
var nlpFir = []float64{
	-1.0818124e-03, -1.1008344e-03, -9.2768838e-04, -4.2289438e-04,
	5.5034190e-04, 2.0029849e-03, 3.7058509e-03, 5.1449415e-03,
	5.5924666e-03, 4.3036754e-03, 8.0284511e-04, -4.8204610e-03,
	-1.1705810e-02, -1.8199275e-02, -2.2065282e-02, -2.0920610e-02,
	-1.2808831e-02, 3.2204775e-03, 2.6683811e-02, 5.5520624e-02,
	8.6305944e-02, 1.1480192e-01, 1.3674206e-01, 1.4867556e-01,
	1.4867556e-01, 1.3674206e-01, 1.1480192e-01, 8.6305944e-02,
	5.5520624e-02, 2.6683811e-02, 3.2204775e-03, -1.2808831e-02,
	-2.0920610e-02, -2.2065282e-02, -1.8199275e-02, -1.1705810e-02,
	-4.8204610e-03, 8.0284511e-04, 4.3036754e-03, 5.5924666e-03,
	5.1449415e-03, 3.7058509e-03, 2.0029849e-03, 5.5034190e-04,
	-4.2289438e-04, -9.2768838e-04, -1.1008344e-03, -1.0818124e-03,
}

// PCMBuffer defines PCM samples as int16.
type PCMBuffer []int16

// Frame and quantization constants.
var (
	MPitchS            = 0.0400 // pitch analysis window in s
	PMinS              = 0.0025 // minimum pitch period in s
	PMaxS              = 0.0200 // maximum pitch period in s
	TWS                = 0.0050 // trapezoidal synthesis window in s
	FrameLengthSecs    = 0.01   // internal proc frame length in secs
	SampleRate         = 8000
	SamplesPerSubFrame = int(math.Round(float64(SampleRate) * FrameLengthSecs))
	SamplesPerFrame    = 2 * SamplesPerSubFrame
	BitsPerFrame       = 48
	BytesPerFrame      = 6

	LpcOrder = 10

	WoBits           = 7
	WoEBits          = 8
	LSPScalarIndexes = 10

	Mode2400 = 0x01
)

// COMP represents a complex number.
type COMP struct {
	Real float64
	Imag float64
}

// Model represents the codec2 model parameters for one frame.
type Model struct {
	Wo     float64   // Fundamental frequency estimate in radians.
	L      int       // Number of harmonics.
	A      []float64 // Amplitudes for harmonics (indices 1..MaxAmp; index 0 unused).
	Phi    []float64 // Phases for harmonics (indices 1..MaxAmp; index 0 unused).
	Voiced bool      // Voiced flag.
	LSPs   []float64 // LSPs (length = LpcOrder).
	E      float64   // Energy.
}

// C2Const holds constants calculated at run time.
type C2Const struct {
	Fs     int     // Sample rate.
	NSamp  int     // Number of samples per 10ms frame.
	MaxAmp int     // Maximum number of harmonics.
	MPitch int     // Pitch estimation window size in samples.
	PMin   int     // Minimum pitch period in samples.
	PMax   int     // Maximum pitch period in samples.
	WoMin  float64 // Minimum fundamental frequency in radians.
	WoMax  float64 // Maximum fundamental frequency in radians.
	Nw     int     // Analysis window size in samples.
	Tw     int     // Trapezoidal synthesis window overlap.
}

// NLP holds the state for the nonlinear pitch estimator.
type NLP struct {
	Fs      int         // Sample rate.
	m       int         // Analysis window length.
	sq      []float64   // Squared speech samples.
	mem_x   float64     // Notch filter memory.
	mem_y   float64     // Notch filter memory.
	mem_fir []float64   // FIR filter memory (length = NLP_NTAP).
	fft_cfg interface{} // FFT configuration (should be of type FFT)
	w       []float64   // Analysis window for decimated signal (length = m/DEC)
	pmin    int         // Minimum pitch period in samples.
	pmax    int         // Maximum pitch period in samples.
}

// LPC
const (
	LPC_MAX_N = 512 // Maximum number of samples in frame.
	ALPHA     = 1.0
	BETA      = 0.94
)

// Our primary codec context
type Codec2 struct {
	mode     int       // Only 2400bps mode supported.
	nsam     int       // Number of speech samples per subframe (e.g., 80).
	nbit     int       // Number of bits per frame (48 for 2400bps).
	c2const  C2Const   // Run-time constants.
	w        []float64 // Hamming window (length = MPitch)
	Sn       []float64 // Analysis buffer (length = MPitch)
	Sn_      []float64 // Synthesis buffer (length ≥ SamplesPerFrame)
	model    Model     // Current frame's model.
	prevLsps []float64 // Previous frame's LSPs (length = LpcOrder)
	prevE    float64   // Previous frame's LPC energy.
	xq_enc   []float64 // Encoder quantizer state (length = 2).
	xq_dec   []float64 // Decoder quantizer state (length = 2).

	// Additional fields.
	prev_Wo    float64 // Previous frame's Wo.
	prev_model Model   // Previous frame's model.
	lpc_pf     bool    // LPC post filter flag.
	bass_boost bool    // Bass boost flag.
	beta       float64 // LPC post filter beta.
	gamma      float64 // LPC post filter gamma.

	// Fields for analysis.
	fftFwdCfg FFT         // FFT configuration for forward FFT.
	fftInvCfg FFT         // NEW: FFT configuration for inverse FFT.
	W         []float64   // Precomputed frequency-domain window vector.
	prevF0Enc float64     // Previous pitch estimate.
	nlp       interface{} // NLP state (should be *NLP)
	fmlfeat   interface{} // ML feature output handle (if any)

	// NEW: Fields required for synthesis:
	ex_phase float64   // Excitation model phase track.
	bg_est   float64   // Background noise estimate for post filter.
	Pn       []float64 // Trapezoidal synthesis window (length = 2 * SamplesPerSubFrame)
}
