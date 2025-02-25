package codec2

// earProtection replicates the C function ear_protection().
// It finds the maximum sample in the buffer and, if above a threshold,
// scales the entire buffer down.
func earProtection(inOut []float64, n int) {
	maxSample := 0.0
	for i := 0; i < n; i++ {
		if abs(inOut[i]) > maxSample {
			maxSample = abs(inOut[i])
		}
	}
	threshold := 30000.0
	if maxSample > threshold {
		gain := 1.0 / (maxSample * maxSample)
		for i := 0; i < n; i++ {
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
