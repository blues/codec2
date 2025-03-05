package codec2

// earProtection replicates the C function ear_protection().
// It finds the maximum sample in the buffer and, if above a threshold,
// scales the entire buffer down.
func earProtection(inOut []float64, n int) {

	// find maximum sample in frame

	maxSample := 0.0
	for i := 0; i < n; i++ {
		if inOut[i] > maxSample {
			maxSample = inOut[i]
		}
	}

	// determine how far above set point

	over := maxSample / 30000.0

	// If we are x B over set point we reduce level by 2x dB, this
	// attenuates major excursions in amplitude (likely to be caused
	// by bit errors) more than smaller ones
	if over > 1.0 {
		gain := 1.0 / (over * over);
		for i := 0; i<n; i++ {
			inOut[i] *= gain
		}
	}

}

func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}
