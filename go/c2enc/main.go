package main

import (
	"encoding/binary"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"

	"github.com/blues/codec2/go"
)

func usage() {
	fmt.Fprintf(os.Stderr, "usage: %s InputRawspeechFile OutputBitFile\n", os.Args[0])
	fmt.Fprintf(os.Stderr, "e.g. (headerless)    %s input.raw output.bin\n", os.Args[0])
	fmt.Fprintf(os.Stderr, "e.g. (with header)   %s input.raw output.c2\n", os.Args[0])
	os.Exit(1)
}

func main() {
	if len(os.Args) != 3 {
		usage()
	}

	inputFile := os.Args[1]
	outputFile := os.Args[2]

	// Create codec instance
	codec, err := codec2.NewCodec2()
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating codec: %v\n", err)
		os.Exit(1)
	}

	// Open input file or use stdin
	var fin *os.File
	if inputFile == "-" {
		fin = os.Stdin
	} else if fin, err = os.Open(inputFile); err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input speech file: %s: %s\n", inputFile, err.Error())
		os.Exit(1)
	}
	defer fin.Close()

	// Open output file or use stdout
	var fout *os.File
	if outputFile == "-" {
		fout = os.Stdout
	} else if fout, err = os.Create(outputFile); err != nil {
		fmt.Fprintf(os.Stderr, "Error opening output compressed bit file: %s: %s\n", outputFile, err.Error())
		os.Exit(1)
	}
	defer fout.Close()

	// Add c2 header if output file has .c2 extension
	if strings.ToLower(filepath.Ext(outputFile)) == ".c2" {
		header := codec2.NewHeader(codec2.Mode2400Byte)
		if err := binary.Write(fout, binary.LittleEndian, header); err != nil {
			fmt.Fprintf(os.Stderr, "Error writing header: %v\n", err)
			os.Exit(1)
		}
	}

	// Process input in frames
	frame := make([]byte, codec2.SamplesPerFrame*2) // 160 samples * 2 bytes per sample
	frameCount := 0

	for {
		frameCount++
		// Read a frame of raw PCM samples
		_, err := io.ReadFull(fin, frame)
		if err == io.EOF || err == io.ErrUnexpectedEOF {
			break // Don't process partial frames
		}
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error reading input: %v\n", err)
			os.Exit(1)
		}

		// Convert bytes to PCM samples
		pcm := make(codec2.PCMBuffer, codec2.SamplesPerFrame)
		for j := 0; j < codec2.SamplesPerFrame; j++ {
			pcm[j] = int16(binary.LittleEndian.Uint16(frame[j*2:]))
		}

		// Encode frame to compressed format
		encodedFrame, err := codec.Encode(pcm)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error encoding frame: %v\n", err)
			os.Exit(1)
		}
		_, err = fout.Write(encodedFrame)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error writing output: %v\n", err)
			os.Exit(1)
		}

		// Print frame number to stderr if reading from stdin
		if fin == os.Stdin {
			fmt.Fprintf(os.Stderr, "Frame: %d\r", frameCount)
		}

		// Flush stdout for pipeline usage
		if fout == os.Stdout {
			fout.Sync()
		}
	}
}
