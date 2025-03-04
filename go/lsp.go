package codec2

import "math"

// LpcToLsp converts LPC coefficients a (length order+1, with a[0]==1)
// into LSP frequencies (in radians) stored in lsp (length at least order).
// nb is the number of sub-intervals (e.g. 5) and delta is the grid spacing (e.g. LSP_DELTA1).
// It returns the number of roots found.
func LpcToLsp(a []float64, lsp []float64, order int, nb int, delta float64) int {
	var psuml, psumr, psumm, tempXr, xl, xr, xm float64
	var tempSumr float64
	var flag bool
	var k int
	roots := 0
	m := order / 2

	// Allocate arrays for P and Q (length order+1).
	P := make([]float64, order+1)
	Q := make([]float64, order+1)

	// Initialize P[0] and Q[0] to 1.
	P[0] = 1.0
	Q[0] = 1.0

	// Compute the coefficients for P and Q.
	for i := 1; i <= m; i++ {
		P[i] = a[i] + a[order+1-i] - P[i-1]
		Q[i] = a[i] - a[order+1-i] + Q[i-1]
	}

	// Multiply the first m coefficients by 2.
	for i := 0; i < m; i++ {
		P[i] *= 2.0
		Q[i] *= 2.0
	}

	// Initialize xl and xr.
	xl = 1.0
	xr = 0.0

	// For each expected LSP (j = 0 .. order-1) alternate between P and Q.
	for j := 0; j < order; j++ {
		var poly []float64
		if j%2 == 1 {
			poly = Q
		} else {
			poly = P
		}
		// IMPORTANT: pass only the first m+1 coefficients.
		psuml = chebPolyEva(poly[:m+1], xl, order)
		flag = true
		for flag && (xr >= -1.0) {
			xr = xl - delta
			psumr = chebPolyEva(poly[:m+1], xr, order)
			tempSumr = psumr
			tempXr = xr
			if (psumr*psuml < 0.0) || (psumr == 0.0) {
				roots++
				psumm = psuml
				for k = 0; k <= nb; k++ {
					xm = (xl + xr) / 2.0
					psumm = chebPolyEva(poly[:m+1], xm, order)
					if psumm*psuml > 0.0 {
						psuml = psumm
						xl = xm
					} else {
						/* psumr = psumm */
						xr = xm
					}
				}
				lsp[j] = xm
				xl = xm
				flag = false
			} else {
				psuml = tempSumr
				xl = tempXr
			}
		}
	}

	// Convert the computed LSP values from the cosine domain to radians.
	for i := 0; i < order; i++ {
		if lsp[i] > 1.0 {
			lsp[i] = 1.0
		} else if lsp[i] < -1.0 {
			lsp[i] = -1.0
		}
		lsp[i] = math.Acos(lsp[i])
	}

	return roots
}

// chebPolyEva evaluates a Chebyshev polynomial series for the given coefficients
// at x. The slice coef should have length order/2+1. It mimics the C function
// cheb_poly_eva.
func chebPolyEva(coef []float64, x float64, order int) float64 {
	n := order/2 + 1
	T := make([]float64, n)
	T[0] = 1.0
	if n > 1 {
		T[1] = x
	}
	for i := 2; i < n; i++ {
		T[i] = 2*x*T[i-1] - T[i-2]
	}
	sum := 0.0
	for i := 0; i < n; i++ {
		// The C code uses coef[(order/2) - i]
		sum += coef[(order/2)-i] * T[i]
	}
	return sum
}

// checkLspOrder checks that the LSP values are in ascending order.
// If not, it “swaps” (adjusts) them slightly and restarts the check.
// Returns the number of swaps performed.
func checkLspOrder(lsp []float64, order int) int {
	swaps := 0
	// Loop from index 1 to order-1.
	// (Note: order is assumed to be len(lsp))
	for i := 1; i < order; i++ {
		if lsp[i] < lsp[i-1] {
			swaps++
			tmp := lsp[i-1]
			lsp[i-1] = lsp[i] - 0.1
			lsp[i] = tmp + 0.1
			// Restart the loop (set i to 0 so next iteration i becomes 1)
			i = 0
		}
	}
	return swaps
}

// bwExpandLsps applies bandwidth expansion to the LSPs.
// It forces a minimum separation between consecutive LSP values.
// minSepLow is used for the first four LSPs and minSepHigh for the rest.
// The LSPs are assumed to be in radians.
func bwExpandLsps(lsp []float64, order int, minSepLow, minSepHigh float64) {
	factor := math.Pi / 4000.0
	// For the first 4 LSPs:
	for i := 1; i < 4 && i < order; i++ {
		if lsp[i]-lsp[i-1] < minSepLow*factor {
			lsp[i] = lsp[i-1] + minSepLow*factor
		}
	}
	// For LSPs from index 4 up to order-1:
	for i := 4; i < order; i++ {
		if lsp[i]-lsp[i-1] < minSepHigh*factor {
			lsp[i] = lsp[i-1] + minSepHigh*factor
		}
	}
}

// lspToLpc converts LSP frequencies (in radians) to LPC coefficients.
// It is a direct transliteration of the C function lsp_to_lpc().
// 'lsp' should be a slice of length 'order' containing LSP frequencies in radians.
// 'ak' must be a slice of length order+1; the resulting LPC coefficients are stored in ak.
func lspToLpc(lsp []float64, ak []float64, order int) {
	// Create a temporary slice to hold cos(lsp[i]) values.
	xfreq := make([]float64, order)
	for i := 0; i < order; i++ {
		xfreq[i] = math.Cos(lsp[i])
	}

	// Allocate a temporary buffer Wp of length (order*4 + 2), zero-initialized.
	sizeWp := order*4 + 2
	Wp := make([]float64, sizeWp)
	for i := 0; i < sizeWp; i++ {
		Wp[i] = 0.0
	}

	// Set xin1 and xin2 to 1.0.
	xin1 := 1.0
	xin2 := 1.0

	// We will compute order+1 LPC coefficients (ak[0]...ak[order]).
	// Outer loop: for j from 0 to order (inclusive)
	// Note: In the C code, the pointer arithmetic is used; here we simulate it.
	// Let N = order/2 (integer division)
	N := order / 2
	// Outer loop: j=0 ... order
	for j := 0; j <= order; j++ {
		// For each j, perform inner loop over i = 0 ... N-1
		for i := 0; i < N; i++ {
			// Compute the starting index for block i in Wp.
			idx := i * 4
			// Compute xout1 and xout2:
			xout1 := xin1 - 2.0*xfreq[2*i]*Wp[idx] + Wp[idx+1]
			xout2 := xin2 - 2.0*xfreq[2*i+1]*Wp[idx+2] + Wp[idx+3]
			// Shift: copy current block to next positions.
			Wp[idx+1] = Wp[idx]
			Wp[idx+3] = Wp[idx+2]
			Wp[idx] = xin1
			Wp[idx+2] = xin2
			// Update xin1 and xin2.
			xin1 = xout1
			xin2 = xout2
		}
		// After inner loop, let n4Index be the index of the last block's fourth element.
		if N > 0 {
			n4Index := (N-1)*4 + 3
			// Compute xout1 and xout2 using the elements following the last block.
			// In the C code: xout1 = xin1 + *(n4+1), xout2 = xin2 - *(n4+2)
			// Here, we assume that Wp[n4Index+1] and Wp[n4Index+2] exist.
			xout1 := xin1 + Wp[n4Index+1]
			xout2 := xin2 - Wp[n4Index+2]
			// Compute LPC coefficient for this iteration.
			ak[j] = 0.5 * (xout1 + xout2)
			// Update the temporary buffer.
			Wp[n4Index+1] = xin1
			Wp[n4Index+2] = xin2
		} else {
			// If order is 0, set ak[0] = 1.0.
			ak[j] = 1.0
		}
		// Reset xin1 and xin2 for next iteration.
		xin1 = 0.0
		xin2 = 0.0
	}
}
