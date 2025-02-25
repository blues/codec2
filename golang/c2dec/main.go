package main

import (
	"encoding/binary"
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/blues/codec2/golang"
)

func usage() {
	fmt.Fprintf(os.Stderr, "usage: %s InputBitFile OutputRawFile\n", os.Args[0])
	fmt.Fprintf(os.Stderr, "e.g. (headerless)    %s input.bin output.raw\n", os.Args[0])
	fmt.Fprintf(os.Stderr, "e.g. (with header)   %s input.c2 output.raw\n", os.Args[0])
	os.Exit(1)
}

func main() {
	if len(os.Args) != 3 {
		usage()
	}

	inputFile := os.Args[1]
	outputFile := os.Args[2]

	// Handle stdin/stdout
	var input []byte
	var err error

	if inputFile == "-" {
		input, err = io.ReadAll(os.Stdin)
	} else {
		input, err = os.ReadFile(inputFile)
	}
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error reading input: %v\n", err)
		os.Exit(1)
	}

	// Check for C2 file header
	var data []byte
	if len(input) >= 7 && codec2.IsC2Header(input) {
		mode := codec2.GetC2Mode(input)
		if mode != codec2.Mode2400 {
			fmt.Fprintf(os.Stderr, "Error: input file mode %d not supported\n", mode)
			fmt.Fprintf(os.Stderr, "header: %X\n", input[0:7])
			os.Exit(1)
		}
		data = input[7:] // Skip header
	} else {
		data = input
	}

	// Create codec instance
	codec, err := codec2.NewCodec2()
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating codec: %v\n", err)
		os.Exit(1)
	}

	// Process input in frames
	frameSize := 6 // 48 bits = 6 bytes per frame
	var decoded []byte

	// Decode frames
	for i := 0; i < len(data); i += frameSize {
		end := i + frameSize
		if end > len(data) {
			break // Don't process partial frames
		}

		frame := data[i:end]
		pcm, err := codec.Decode(frame)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error decoding frame %d: %v\n", i/frameSize+1, err)
			os.Exit(1)
		}

		// Convert PCM samples to bytes
		decodedFrame := make([]byte, codec2.SamplesPerFrame*2)
		for j := 0; j < codec2.SamplesPerFrame; j++ {
			binary.LittleEndian.PutUint16(decodedFrame[j*2:], uint16(pcm[j]))
		}
		decoded = append(decoded, decodedFrame...)

		// Print frame number to stderr if reading from stdin
		if inputFile == "-" {
			fmt.Fprintf(os.Stderr, "Frame: %d\r", i/frameSize+1)
		}
	}

	// Write output
	var out io.Writer
	if outputFile == "-" {
		// Don't write to stdout if running tests
		if flag.Lookup("test.v") != nil {
			return
		}
		out = os.Stdout
	} else {
		f, err := os.Create(outputFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
			os.Exit(1)
		}
		defer f.Close()
		out = f
	}

	_, err = out.Write(decoded)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error writing output: %v\n", err)
		os.Exit(1)
	}
}
