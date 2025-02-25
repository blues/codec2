package codec2

import (
	"math"
)

// interp_Wo interpolates the fundamental frequency (Wo) between frames.
func interp_Wo(interp *Model, prev *Model, next *Model, woMin float64) {
	interp_Wo2(interp, prev, next, 0.5, woMin)
}

// interp_Wo2 interpolates the fundamental frequency (Wo) between frames.
func interp_Wo2(interp *Model, prev *Model, next *Model, weight float64, woMin float64) {
	if interp.Voiced && !prev.Voiced && !next.Voiced {
		interp.Voiced = false
	}
	if interp.Voiced {
		if prev.Voiced && next.Voiced {
			interp.Wo = (1.0-weight)*prev.Wo + weight*next.Wo
		}
		if !prev.Voiced && next.Voiced {
			interp.Wo = next.Wo
		}
		if prev.Voiced && !next.Voiced {
			interp.Wo = prev.Wo
		}
	} else {
		interp.Wo = woMin
	}
	interp.L = int(PI / interp.Wo)
}

// interp_energy interpolates energy between frames.
func interp_energy(prevE, nextE float64) float64 {
	// Geometric mean (linear interpolation in log domain).
	return math.Sqrt(prevE * nextE)
}

// interpolate_lsp_ver2 interpolates LSP parameters between frames.
func interpolate_lsp_ver2(interp []float64, prev []float64, next []float64, weight float64, order int) {
	for i := 0; i < order; i++ {
		interp[i] = (1.0-weight)*prev[i] + weight*next[i]
	}
}
