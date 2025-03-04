package codec2

import (
	"math"
)

// -----------------------------------------------------------------------------
// Codebook type definitions
// -----------------------------------------------------------------------------

// LspCodebook mirrors the C struct lsp_codebook.
type LspCodebook struct {
	K     int       // Dimension of each vector.
	Log2M int       // Number of bits in the index.
	M     int       // Number of codebook entries.
	CB    []float64 // The codebook data.
}

// GeCodebook mirrors the GE codebook structure.
type GeCodebook struct {
	K     int          // Dimension (should be 2)
	Log2M int          // Number of bits (8)
	M     int          // Number of entries (256)
	CB    [][2]float64 // The codebook data (256 entries)
}

// -----------------------------------------------------------------------------
// LSP codebooks (from lsp_cb.c)
// -----------------------------------------------------------------------------

var codes0 = []float64{
	225, 250, 275, 300, 325, 350, 375, 400,
	425, 450, 475, 500, 525, 550, 575, 600,
}
var codes1 = []float64{
	325, 350, 375, 400, 425, 450, 475, 500,
	525, 550, 575, 600, 625, 650, 675, 700,
}
var codes2 = []float64{
	500, 550, 600, 650, 700, 750, 800, 850,
	900, 950, 1000, 1050, 1100, 1150, 1200, 1250,
}
var codes3 = []float64{
	700, 800, 900, 1000, 1100, 1200, 1300, 1400,
	1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200,
}
var codes4 = []float64{
	950, 1050, 1150, 1250, 1350, 1450, 1550, 1650,
	1750, 1850, 1950, 2050, 2150, 2250, 2350, 2450,
}
var codes5 = []float64{
	1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800,
	1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600,
}
var codes6 = []float64{
	1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200,
	2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000,
}
var codes7 = []float64{
	2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000,
}
var codes8 = []float64{
	2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200,
}
var codes9 = []float64{
	2900, 3100, 3300, 3500,
}

// lspCb is the LSP scalar codebook array used by encode_lsps_scalar().
var lspCb = []LspCodebook{
	{K: 1, Log2M: 4, M: 16, CB: codes0},
	{K: 1, Log2M: 4, M: 16, CB: codes1},
	{K: 1, Log2M: 4, M: 16, CB: codes2},
	{K: 1, Log2M: 4, M: 16, CB: codes3},
	{K: 1, Log2M: 4, M: 16, CB: codes4},
	{K: 1, Log2M: 4, M: 16, CB: codes5},
	{K: 1, Log2M: 4, M: 16, CB: codes6},
	{K: 1, Log2M: 3, M: 8, CB: codes7},
	{K: 1, Log2M: 3, M: 8, CB: codes8},
	{K: 1, Log2M: 2, M: 4, CB: codes9},
}

// -----------------------------------------------------------------------------
// GE codebook (from ge_cb.c)
// -----------------------------------------------------------------------------

var geCodebookData = [][2]float64{
	{2.71, 12.0184},
	{0.04675, -2.73881},
	{0.120993, 8.38895},
	{-1.58028, -0.892307},
	{1.19307, -1.91561},
	{0.187101, -3.27679},
	{0.332251, -7.66455},
	{-1.47944, 31.2461},
	{1.52761, 27.7095},
	{-0.524379, 5.25012},
	{0.55333, 7.4388},
	{-0.843451, -1.95299},
	{2.26389, 8.61029},
	{0.143143, 2.36549},
	{0.616506, 1.28427},
	{-1.71133, 22.0967},
	{1.00813, 17.3965},
	{-0.106718, 1.41891},
	{-0.136246, 14.2736},
	{-1.70909, -20.5319},
	{1.65787, -3.39107},
	{0.138049, -4.95785},
	{0.536729, -1.94375},
	{0.196307, 36.8519},
	{1.27248, 22.5565},
	{-0.670219, -1.90604},
	{0.382092, 6.40113},
	{-0.756911, -4.90102},
	{1.82931, 4.6138},
	{0.318794, 0.73683},
	{0.612815, -2.07505},
	{-0.410151, 24.7871},
	{1.77602, 13.1909},
	{0.106457, -0.104492},
	{0.192206, 10.1838},
	{-1.82442, -7.71565},
	{0.931346, 4.34835},
	{0.308813, -4.086},
	{0.397143, -11.8089},
	{-0.048715, 41.2273},
	{0.877342, 35.8503},
	{-0.759794, 0.476634},
	{0.978593, 7.67467},
	{-1.19506, 3.03883},
	{2.63989, -3.41106},
	{0.191127, 3.60351},
	{0.402932, 1.0843},
	{-2.15202, 18.1076},
	{1.5468, 8.32271},
	{-0.143089, -4.07592},
	{-0.150142, 5.86674},
	{-1.40844, -3.2507},
	{1.56615, -10.4132},
	{0.178171, -10.2267},
	{0.362164, -0.028556},
	{-0.070125, 24.3907},
	{0.594752, 17.4828},
	{-0.28698, -6.90407},
	{0.464818, 10.2055},
	{-1.00684, -14.3572},
	{2.32957, -3.69161},
	{0.335745, 2.40714},
	{1.01966, -3.15565},
	{-1.25945, 7.9919},
	{2.38369, 19.6806},
	{-0.094947, -2.41374},
	{0.20933, 6.66477},
	{-2.22103, 1.37986},
	{1.29239, 2.04633},
	{0.243626, -0.890741},
	{0.428773, -7.19366},
	{-1.11374, 41.3414},
	{2.6098, 31.1405},
	{-0.446468, 2.53419},
	{0.490104, 4.62757},
	{-1.11723, -3.24174},
	{1.79156, 8.41493},
	{0.156012, 0.183336},
	{0.532447, 3.15455},
	{-0.764484, 18.514},
	{0.952395, 11.7713},
	{-0.332567, 0.346987},
	{0.202165, 14.7168},
	{-2.12924, -15.559},
	{1.35358, -1.92679},
	{-0.010963, -16.3364},
	{0.399053, -2.79057},
	{0.750657, 31.1483},
	{0.655743, 24.4819},
	{-0.45321, -0.735879},
	{0.2869, 6.5467},
	{-0.715673, -12.3578},
	{1.54849, 3.87217},
	{0.271874, 0.802339},
	{0.502073, -4.85485},
	{-0.497037, 17.7619},
	{1.19116, 13.9544},
	{0.01563, 1.33157},
	{0.341867, 8.93537},
	{-2.31601, -5.39506},
	{0.75861, 1.9645},
	{0.24132, -3.23769},
	{0.267151, -11.2344},
	{-0.273126, 32.6248},
	{1.75352, 40.432},
	{-0.784011, 3.04576},
	{0.705987, 5.66118},
	{-1.3864, 1.35356},
	{2.37646, 1.67485},
	{0.242973, 4.73218},
	{0.491227, 0.354061},
	{-1.60676, 8.65895},
	{1.16711, 5.9871},
	{-0.137601, -12.0417},
	{-0.251375, 10.3972},
	{-1.43151, -8.90411},
	{0.98828, -13.209},
	{0.261484, -6.35497},
	{0.395932, -0.702529},
	{0.283704, 26.8996},
	{0.420959, 15.4418},
	{-0.355804, -13.7278},
	{0.527372, 12.3985},
	{-1.16956, -15.9985},
	{1.90669, -5.81605},
	{0.354492, 3.85157},
	{0.82576, -4.16264},
	{-0.49019, 13.0572},
	{2.25577, 13.5264},
	{-0.004956, -3.23713},
	{0.026709, 7.86645},
	{-1.81037, -0.451183},
	{1.08383, -0.18362},
	{0.135836, -2.26658},
	{0.375812, -5.51225},
	{-1.96644, 38.6829},
	{1.97799, 24.5655},
	{-0.704656, 6.35881},
	{0.480786, 7.05175},
	{-0.976417, -2.42273},
	{2.50215, 6.75935},
	{0.083588, 3.2588},
	{0.543629, 0.910013},
	{-1.23196, 23.0915},
	{0.785492, 14.807},
	{-0.213554, 1.688},
	{0.004748, 18.1718},
	{-1.54719, -16.1168},
	{1.50104, -3.28114},
	{0.080133, -4.63472},
	{0.476592, -2.18093},
	{0.44247, 40.304},
	{1.07277, 27.592},
	{-0.594738, -4.16681},
	{0.42248, 7.61609},
	{-0.927521, -7.27441},
	{1.99162, 1.29636},
	{0.291307, 2.39878},
	{0.721081, -1.95062},
	{-0.804256, 24.9295},
	{1.64839, 19.1197},
	{0.060852, -0.590639},
	{0.266085, 9.10325},
	{-1.9574, -2.88461},
	{1.11693, 2.6724},
	{0.35458, -2.74854},
	{0.330733, -14.1561},
	{-0.527851, 39.5756},
	{0.991152, 43.195},
	{-0.589619, 1.26919},
	{0.787401, 8.73071},
	{-1.0138, 1.02507},
	{2.8254, 1.89538},
	{0.24089, 2.74557},
	{0.427195, 2.54446},
	{-1.95311, 12.244},
	{1.44862, 12.0607},
	{-0.210492, -3.37906},
	{-0.056713, 10.204},
	{-1.65237, -5.10274},
	{1.29475, -12.2708},
	{0.111608, -8.67592},
	{0.326634, -1.16763},
	{0.021781, 31.1258},
	{0.455335, 21.4684},
	{-0.37544, -3.37121},
	{0.39362, 11.302},
	{-0.851456, -19.4149},
	{2.10703, -2.22886},
	{0.373233, 1.92406},
	{0.884438, -1.72058},
	{-0.975127, 9.84013},
	{2.0033, 17.3954},
	{-0.036915, -1.11137},
	{0.148456, 5.39997},
	{-1.91441, 4.77382},
	{1.44791, 0.537122},
	{0.194979, -1.03818},
	{0.495771, -9.95502},
	{-1.05899, 32.9471},
	{2.01122, 32.4544},
	{-0.30965, 4.71911},
	{0.436082, 4.63552},
	{-1.23711, -1.25428},
	{2.02274, 9.42834},
	{0.190342, 1.46077},
	{0.479017, 2.48479},
	{-1.07848, 16.2217},
	{1.20764, 9.65421},
	{-0.258087, -1.67236},
	{0.071852, 13.416},
	{-1.87723, -16.072},
	{1.28957, -4.87118},
	{0.067713, -13.4427},
	{0.435551, -4.1655},
	{0.46614, 30.5895},
	{0.904895, 21.598},
	{-0.518369, -2.53205},
	{0.337363, 5.63726},
	{-0.554975, -17.4005},
	{1.69188, 1.14574},
	{0.227934, 0.889297},
	{0.587303, -5.72973},
	{-0.262133, 18.6666},
	{1.39505, 17.0029},
	{-0.01909, 4.30838},
	{0.304235, 12.6699},
	{-2.07406, -6.46084},
	{0.920546, 1.21296},
	{0.284927, -1.78547},
	{0.209724, -16.024},
	{-0.636067, 31.5768},
	{1.34989, 34.6775},
	{-0.971625, 5.30086},
	{0.590249, 4.44971},
	{-1.56787, 3.60239},
	{2.1455, 4.51666},
	{0.296022, 4.12017},
	{0.445299, 0.868772},
	{-1.44193, 14.1284},
	{1.35575, 6.0074},
	{-0.012814, -7.49657},
	{-0.43, 8.50012},
	{-1.20469, -7.11326},
	{1.10102, -6.83682},
	{0.196463, -6.234},
	{0.436747, -1.12979},
	{0.141052, 22.8549},
	{0.290821, 18.8114},
	{-0.529536, -7.73251},
	{0.63428, 10.7898},
	{-1.33472, -20.3258},
	{1.81564, -1.90332},
	{0.394778, 3.79758},
	{0.732682, -8.18382},
	{-0.741244, 11.7683},
}

// Create a global GeCodebook instance.
var geCb = GeCodebook{
	K:     2,
	Log2M: 8,
	M:     256,
	CB:    geCodebookData,
}

// -----------------------------------------------------------------------------
// Quantization functions
// -----------------------------------------------------------------------------

// quantize is a direct transliteration of the C quantise() function.
func quantize(cb []float64, vec []float64, w []float64, k int, m int, se *float64) int {
	bestIndex := 0
	bestErr := 1.0e32
	for j := 0; j < m; j++ {
		var e float64
		for i := 0; i < k; i++ {
			diff := cb[j*k+i] - vec[i]
			e += diff * w[i] * diff * w[i]
		}
		if e < bestErr {
			bestErr = e
			bestIndex = j
		}
	}
	*se += bestErr
	return bestIndex
}

// --- Joint Wo/E Quantisation ---

// geCoeff is defined identically to the C version.
var geCoeff = []float64{0.8, 0.9}

// encodeWoE is a direct transliteration of C's encode_WoE().
// It updates the state vector xq (passed as a pointer) and returns the chosen index.
func encodeWoE(model *Model, e float64, xq []float64) int {
	var x [2]float64
	var err [2]float64
	var nb_entries = geCb.M
	var ndim = geCb.K
	codebook1 := geCb.CB

	if e < 0.0 {
		e = 0.0
	}
	x[0] = math.Log10((model.Wo/math.Pi)*4000.0/50.0) / math.Log10(2.0)
	x[1] = 10.0 * math.Log10(1e-4+e)

	// Compute weights using computeWeights2 (which now exactly squares the weights).
	w := make([]float64, len(xq))
	computeWeights2(x[:], xq, w)

	for i := 0; i < 2; i++ {
		err[i] = x[i] - geCoeff[i]*xq[i]
	}
	// Find the nearest codebook vector from geCb.
	n1 := findNearestWeighted(codebook1, nb_entries, err, w, ndim)
	// Update state vector xq.
	for i := 0; i < ndim; i++ {
		xq[i] = geCoeff[i]*xq[i] + codebook1[n1][i]
	}
	return n1
}

// findNearestWeighted finds the index of the codebook entry (in a codebook
// represented as a slice of [2]float64 vectors) that minimizes the weighted
// squared distance from x. The weight vector w and the dimension ndim (which
// should be 2) are used in the computation.
func findNearestWeighted(codebook [][2]float64, nbEntries int, x [2]float64, w []float64, ndim int) int {
	minDist := 1e15
	nearest := 0
	for i := 0; i < nbEntries; i++ {
		dist := 0.0
		for j := 0; j < ndim; j++ {
			diff := x[j] - codebook[i][j]
			dist += w[j] * diff * diff
		}
		if dist < minDist {
			minDist = dist
			nearest = i
		}
	}
	return nearest
}

// decodeWoE is a direct transliteration of C's decode_WoE().
// It updates the state vector xq and sets model.Wo and *e.
func decodeWoE(c2const *C2Const, model *Model, e *float64, xq []float64, index int) {
	for i := 0; i < len(xq); i++ {
		xq[i] = geCoeff[i]*xq[i] + geCb.CB[index][i]
	}
	model.Wo = math.Pow(2.0, xq[0]) * (math.Pi * 50.0) / 4000.0
	if model.Wo > c2const.WoMax {
		model.Wo = c2const.WoMax
	}
	if model.Wo < c2const.WoMin {
		model.Wo = c2const.WoMin
	}
	model.L = int(math.Floor(math.Pi / model.Wo))
	*e = math.Pow(10.0, xq[1]/10.0)
}

// computeWeights2 exactly matches the C compute_weights2: it computes
// weights from x (current vector) and xp (previous state) and then squares them.
func computeWeights2(x, xp, w []float64) {
	w[0] = 30.0
	w[1] = 1.0
	if x[1] < 0 {
		w[0] *= 0.6
		w[1] *= 0.3
	}
	if x[1] < -10 {
		w[0] *= 0.3
		w[1] *= 0.3
	}
	if math.Abs(x[0]-xp[0]) < 0.2 {
		w[0] *= 2.0
		w[1] *= 1.5
	} else if math.Abs(x[0]-xp[0]) > 0.5 {
		w[0] *= 0.5
	}
	if x[1] < xp[1]-10 {
		w[1] *= 0.5
	}
	if x[1] < xp[1]-20 {
		w[1] *= 0.5
	}
	w[0] = w[0] * w[0]
	w[1] = w[1] * w[1]
}

// -----------------------------------------------------------------------------
// LSP Scalar Quantization (difference-based)
// -----------------------------------------------------------------------------

// Get lsp bits
func lspBits(i int) int {
	return lspCb[i].Log2M
}

// encodeLSPScalar replicates the C function encode_lspds_scalar().
// It converts LSPs (in radians) to Hz, computes differences, quantizes each difference,
// and stores the resulting indexes.
func encodeLSPScalar(indexes []int, lsp []float64, order int) {
	var se float64 = 0.0
	lspHz := make([]float64, order)
	// Convert each LSP from radians to Hz
	for i := 0; i < order; i++ {
		lspHz[i] = (4000.0 / math.Pi) * lsp[i]
	}
	// For each coefficient, quantize using the corresponding scalar codebook
	for i := 0; i < order; i++ {
		k := lspCb[i].K
		m := lspCb[i].M
		cb := lspCb[i].CB
		// Quantize the single value (slice of length 1) using a constant weight of 1.0
		indexes[i] = quantize(cb, lspHz[i:i+1], []float64{1.0}, k, m, &se)
	}
}

// decodeLSPScalar replicates the C function decode_lspds_scalar().
// It recovers quantized LSPs (in Hz) from the difference codebook indexes and converts them back to radians.
func decodeLSPScalar(lsps []float64, indexes []int, order int) {
	lspHz := make([]float64, order)
	for i := 0; i < order; i++ {
		k := lspCb[i].K
		cb := lspCb[i].CB
		idx := indexes[i]
		if idx < 0 {
			idx = 0
		}
		if idx >= len(cb)/k {
			idx = (len(cb) / k) - 1
		}
		lspHz[i] = cb[idx*k]
	}
	for i := 0; i < order; i++ {
		lsps[i] = (math.Pi / 4000.0) * lspHz[i]
	}
}

// aksToM2 transforms the LPC coefficients (ak) into a spectral amplitude estimate
// for each harmonic, updating the model's A values. It also computes a signal-to-noise ratio.
// This function is a direct transliteration of the C function aks_to_M2.
func aksToM2(fftFwdCfg FFT, ak []float64, order int, model *Model, E float64, snr *float64, dump int, simPf int, pf bool, bassBoost bool, beta, gamma float64, Aw []COMP) {
	// r = 2*pi / FFTSize
	r := TWO_PI / float64(FFTSize)

	// Create temporary buffer "a" of length FFTSize and zero it.
	a := make([]float64, FFTSize)
	for i := 0; i < FFTSize; i++ {
		a[i] = 0.0
	}
	// Copy LPC coefficients into a[0..order]
	for i := 0; i <= order; i++ {
		a[i] = ak[i]
	}
	// Compute FFT of a to get A(exp(jw)). We assume fftFwdCfg implements our FFT interface.
	fftRes := fftFwdCfg.Forward(a)
	// Store the FFT result into Aw (assume Aw has length FFTSize).
	for i, c := range fftRes {
		Aw[i] = COMP{Real: real(c), Imag: imag(c)}
	}

	// Compute power spectrum Pw (only first FFTSize/2 samples)
	Pw := make([]float64, FFTSize/2)
	for i := 0; i < FFTSize/2; i++ {
		mag2 := Aw[i].Real*Aw[i].Real + Aw[i].Imag*Aw[i].Imag + 1e-6
		Pw[i] = 1.0 / mag2
	}

	// If post-filtering is enabled, call lpcPostFilter (which you must have implemented elsewhere)
	if pf {
		lpcPostFilter(fftFwdCfg, Pw, ak, order, dump, beta, gamma, bassBoost, E)
	} else {
		for i := 0; i < FFTSize/2; i++ {
			Pw[i] *= E
		}
	}

	// Now, compute the harmonic amplitude estimates.
	signal := 1e-30
	noise := 1e-32
	for m := 1; m <= model.L; m++ {
		am := int(((float64(m)-0.5)*model.Wo/r + 0.5))
		bm := int(((float64(m)+0.5)*model.Wo/r + 0.5))
		if bm > FFTSize/2 {
			bm = FFTSize / 2
		}
		Em := 0.0
		for i := am; i < bm; i++ {
			Em += Pw[i]
		}
		Am := math.Sqrt(Em)
		signal += model.A[m] * model.A[m]
		diff := model.A[m] - Am
		noise += diff * diff

		if simPf != 0 {
			if Am > model.A[m] {
				Am *= 0.7
			}
			if Am < model.A[m] {
				Am *= 1.4
			}
		}
		model.A[m] = Am
	}
	*snr = 10.0 * math.Log10(signal/noise)
}

// applyLpcCorrection applies a simple LPC correction: if the fundamental frequency Wo
// is below (pi*150/4000), then it scales the first harmonic amplitude.
func applyLpcCorrection(model *Model) {
	if model.Wo < (math.Pi * 150.0 / 4000.0) {
		model.A[1] *= 0.032
	}
}

// lpcPostFilter applies a post filter to the LPC synthesis filter power spectrum Pw,
// which suppresses the inter-formant energy. This is a direct transliteration of
// the C function lpc_post_filter() from quantise.c.
func lpcPostFilter(fftFwdCfg FFT, Pw []float64, ak []float64, order int, dump int, beta, gamma float64, bassBoost bool, E float64) {
	// Assume FFTSize is defined as FFTSize (corresponding to C FFT_ENC).
	// Allocate an input buffer x of length FFTSize and zero it.
	x := make([]float64, FFTSize)
	for i := 0; i < FFTSize; i++ {
		x[i] = 0.0
	}
	// Set up the weighting filter input: x[0] = ak[0], then for i=1..order, x[i] = ak[i]*coeff,
	// where coeff is initially gamma and is multiplied by gamma after each iteration.
	x[0] = ak[0]
	coeff := gamma
	for i := 1; i <= order; i++ {
		x[i] = ak[i] * coeff
		coeff *= gamma
	}

	// Compute the FFT of x using the provided forward FFT configuration.
	// The C code uses codec2_fftr(), which returns FFTSize/2+1 complex numbers.
	WwRaw := fftFwdCfg.Forward(x) // []complex128; length should be FFTSize/2+1

	// Convert WwRaw to our internal []COMP slice.
	numBins := FFTSize / 2 // We'll use bins 0 .. FFTSize/2 - 1
	Ww := make([]COMP, numBins)
	for i := 0; i < numBins; i++ {
		Ww[i] = COMP{
			Real: real(WwRaw[i]),
			Imag: imag(WwRaw[i]),
		}
	}

	// Square the FFT output so that for each bin, Ww[i].Real = (real^2 + imag^2)
	for i := 0; i < numBins; i++ {
		Ww[i].Real = Ww[i].Real*Ww[i].Real + Ww[i].Imag*Ww[i].Imag
		Ww[i].Imag = 0.0
	}

	// Measure the energy of Pw before post filtering.
	eBefore := 1e-4
	for i := 0; i < numBins; i++ {
		eBefore += Pw[i]
	}

	// Determine the combined filter R = sqrt(Ww * Pw), and compute max and min of R.
	Rw := make([]float64, numBins)
	maxRw := 0.0
	minRw := 1e32
	for i := 0; i < numBins; i++ {
		Rw[i] = math.Sqrt(Ww[i].Real * Pw[i])
		if Rw[i] > maxRw {
			maxRw = Rw[i]
		}
		if Rw[i] < minRw {
			minRw = Rw[i]
		}
	}

	// Apply post filtering: for each bin, multiply Pw[i] by (Rw[i]^beta)^2.
	eAfter := 1e-4
	for i := 0; i < numBins; i++ {
		Pfw := math.Pow(Rw[i], beta)
		Pw[i] = Pw[i] * Pfw * Pfw
		eAfter += Pw[i]
	}

	// Compute gain = (energy before)/(energy after), then scale gain by E.
	gain := eBefore / eAfter
	gain *= E
	for i := 0; i < numBins; i++ {
		Pw[i] *= gain
	}

	// If bass boost is enabled, add 3dB to the first 1 kHz.
	if bassBoost {
		limit := FFTSize / 8
		if limit > numBins {
			limit = numBins
		}
		for i := 0; i < limit; i++ {
			Pw[i] *= 1.4 * 1.4
		}
	}
}
