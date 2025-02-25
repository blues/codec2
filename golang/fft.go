package codec2

import "github.com/mjibson/go-dsp/fft"

// FFT is our FFT interface.
type FFT interface {
	Forward(in []float64) []complex128
	Inverse(in []complex128) []float64
}

// defaultFFT implements FFT using go-dsp/fft.
type defaultFFT struct {
	size int
}

// NewFFT creates a new FFT instance for the given size.
func NewFFT(size int) FFT {
	return &defaultFFT{size: size}
}

// Forward returns the FFT of a real-valued input.
func (f *defaultFFT) Forward(in []float64) []complex128 {
	return fft.FFTReal(in)
}

// Inverse returns the inverse FFT of a complex-valued input.
// NOTE: Unlike the standard go-dsp/fft.IFFT, we do not apply any scaling here.
// This is because our C version (kiss_fft) returns unscaled data.
func (f *defaultFFT) Inverse(in []complex128) []float64 {
	complexOut := fft.IFFT(in)
	realOut := make([]float64, len(complexOut))
	// Do NOT scale by 1/len(in) so that later in synthesise() we can
	// multiply by FFTSize to mimic the C behavior.
	for i, v := range complexOut {
		realOut[i] = real(v)
	}
	return realOut
}
