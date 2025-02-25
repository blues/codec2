package codec2

import (
	"math"
)

// phaseSynthZeroOrder implements the C function phase_synth_zero_order.
// It updates the fundamental excitation phase, then for each harmonic m (1..L)
// it generates an excitation (using a random phase for unvoiced frames),
// filters that excitation with the provided filterPhase (typically computed
// from the LPC synthesis filter) and then sets model.Phi[m] to the phase of
// the filtered excitation.
func phaseSynthZeroOrder(nSamp int, model *Model, exPhase *float64, filterPhase []COMP) {
	// Update excitation phase:
	*exPhase += model.Wo * float64(nSamp)
	*exPhase -= TWO_PI * math.Floor((*exPhase/TWO_PI)+0.5)
	// For each harmonic m=1..model.L:
	for m := 1; m <= model.L; m++ {
		var Ex COMP
		if model.Voiced {
			Ex.Real = math.Cos(*exPhase * float64(m))
			Ex.Imag = math.Sin(*exPhase * float64(m))
		} else {
			// Generate a random phase between 0 and 2π.
			phi := TWO_PI * float64(codec2Rand()) / float64(CODEC2_RAND_MAX)
			Ex.Real = math.Cos(phi)
			Ex.Imag = math.Sin(phi)
		}
		// Multiply the filter phase (filterPhase[m]) by the excitation Ex.
		var A_ COMP
		A_.Real = filterPhase[m].Real*Ex.Real - filterPhase[m].Imag*Ex.Imag
		A_.Imag = filterPhase[m].Imag*Ex.Real + filterPhase[m].Real*Ex.Imag
		// Set model phase for harmonic m:
		newPhi := math.Atan2(A_.Imag, A_.Real+1e-12)
		model.Phi[m] = newPhi
	}
}

// samplePhase is a direct transliteration of the C function sample_phase().
// It computes, for m = 1..model.L, the harmonic phase values H[m] as the conjugate of A[b],
// where b = int(m * model.Wo / (TWO_PI/FFTSize) + 0.5).
// To mimic the C behavior, model.L is clamped to not exceed the allocated slice H.
func samplePhase(model *Model, H []COMP, A []COMP) {
	// Compute r = TWO_PI / FFT_ENC, where FFT_ENC is FFTSize.
	r := TWO_PI / float64(FFTSize)

	// Sample phase at harmonics: for m = 1 .. model.L.
	for m := 1; m <= model.L; m++ {
		b := int(float64(m)*model.Wo/r + 0.5)
		// Clamp b to a valid index.
		if b < 0 {
			b = 0
		} else if b >= len(A) {
			b = len(A) - 1
		}
		// Inline conjugation: set H[m] to the conjugate of A[b].
		H[m] = COMP{
			Real: A[b].Real,
			Imag: -A[b].Imag,
		}
	}
}

// randInt returns a random integer between 0 and CODEC2_RAND_MAX-1.
// Note that this is a direct transliteration of the C function codec2_rand(), so that
// our numbers are identical to the C version.  It is not a good random number generator.
var nextRand uint64 = 1

func codec2Rand() int {
	nextRand = nextRand*1103515245 + 12345
	return int((nextRand / 65536) % 32768)
}
