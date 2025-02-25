package codec2

import "fmt"

const (
	WordSize   = 8   // Size of an unsigned char in bits.
	IndexMask  = 0x7 // Mask to pick the bit index within a byte.
	ShiftRight = 3   // Right-shift amount to convert a bit index to a byte index.
)

// packBits packs a bit field into a bit string using Gray code conversion.
// It is a wrapper around packNaturalOrGray with gray enabled.
func packBits(bitArray []byte, bitIndex *uint, field int, fieldWidth uint) {
	packNaturalOrGray(bitArray, bitIndex, field, fieldWidth, 1)
}

// packNaturalOrGray packs a bit field into a bit string.
// Parameters:
//
//	bitArray: the output byte slice.
//	bitIndex: pointer to the current bit index (in bits).
//	field: the bit field to be packed (as an int).
//	fieldWidth: the width of the field in bits.
//	gray: non-zero if the field should be converted to Gray code.
func packNaturalOrGray(bitArray []byte, bitIndex *uint, field int, fieldWidth uint, gray uint) {
	// Convert the field to Gray code if required.
	if gray != 0 {
		field = (field >> 1) ^ field
	}
	// Loop until the entire fieldWidth has been packed.
	for fieldWidth != 0 {
		bI := *bitIndex
		bitsLeft := WordSize - (bI & IndexMask)
		var sliceWidth uint
		if bitsLeft < fieldWidth {
			sliceWidth = bitsLeft
		} else {
			sliceWidth = fieldWidth
		}
		wordIndex := bI >> ShiftRight
		// Extract the top sliceWidth bits from field and shift into position.
		slice := byte((field >> int(fieldWidth-sliceWidth)) << (bitsLeft - sliceWidth))
		bitArray[wordIndex] |= slice
		*bitIndex += sliceWidth
		fieldWidth -= sliceWidth
	}
}

// unpackBits unpacks a bit field from a bit string using Gray code conversion.
// It is a wrapper around unpackNaturalOrGray with gray enabled.
func unpackBits(bitArray []byte, bitIndex *uint, fieldWidth uint) (int, error) {
	return unpackNaturalOrGray(bitArray, bitIndex, fieldWidth, 1)
}

// unpackNaturalOrGray unpacks a bit field from a bit string.
// Parameters:
//
//	bitArray: the input byte slice.
//	bitIndex: pointer to the current bit index (in bits).
//	fieldWidth: the width of the field in bits.
//	gray: non-zero if the field is encoded in Gray code.
//
// Returns the unpacked field as an int.
func unpackNaturalOrGray(bitArray []byte, bitIndex *uint, fieldWidth uint, gray uint) (int, error) {
	var field uint = 0
	// Loop until the full fieldWidth is extracted.
	for fieldWidth != 0 {
		bI := *bitIndex
		bitsLeft := WordSize - (bI & IndexMask)
		var sliceWidth uint
		if bitsLeft < fieldWidth {
			sliceWidth = bitsLeft
		} else {
			sliceWidth = fieldWidth
		}
		wordIndex := bI >> ShiftRight
		if int(wordIndex) >= len(bitArray) {
			return 0, fmt.Errorf("buffer underrun: word index %d out of range", wordIndex)
		}
		// Extract sliceWidth bits from the current byte.
		value := (uint(bitArray[wordIndex]) >> (bitsLeft - sliceWidth)) & ((1 << sliceWidth) - 1)
		field |= value << (fieldWidth - sliceWidth)
		*bitIndex += sliceWidth
		fieldWidth -= sliceWidth
	}
	var t uint
	if gray != 0 {
		// Convert from Gray code to binary.
		t = field ^ (field >> 8)
		t ^= (t >> 4)
		t ^= (t >> 2)
		t ^= (t >> 1)
	} else {
		t = field
	}
	return int(t), nil
}

// boolToInt converts a boolean to an int (1 if true, 0 if false).
func boolToInt(b bool) int {
	if b {
		return 1
	}
	return 0
}
