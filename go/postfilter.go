package codec2

import "math"

const (
	BG_THRESH       = 40.0
	BG_BETA         = 0.1
	BG_MARGIN       = 6.0
	CODEC2_RAND_MAX = 32767
)

// postfilter implements the C function postfilter.
// It computes the average energy (in dB) across harmonics and, if below a threshold
// and the frame is unvoiced, updates the background noise estimate.
// Then, for voiced frames, it checks each harmonic amplitude against a threshold and
// if below, randomizes its phase.
func postfilter(model *Model, bgEst *float64) {
	e := 1e-12
	for m := 1; m <= model.L; m++ {
		e += model.A[m] * model.A[m]
	}
	e = 10.0 * math.Log10(e/float64(model.L))
	if e < BG_THRESH && !model.Voiced {
		*bgEst = *bgEst*(1.0-BG_BETA) + e*BG_BETA
	}
	uv := 0
	thresh := math.Pow(10.0, (*bgEst+BG_MARGIN)/20.0)
	if model.Voiced {
		for m := 1; m <= model.L; m++ {
			if model.A[m] < thresh {
				model.Phi[m] = (TWO_PI / float64(CODEC2_RAND_MAX)) * float64(codec2Rand())
				uv++
			}
		}
	}
	// (Optional debug dump omitted.)
}
