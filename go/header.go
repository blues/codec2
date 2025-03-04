package codec2

// C2 file header constants
var (
	MagicByte1   = byte(0xc0)
	MagicByte2   = byte(0xde)
	MagicByte3   = byte(0xc2)
	VersionMajor = byte(0x01)
	VersionMinor = byte(0x00)
	Mode2400Byte = byte(Mode2400)
)

var Magic = []byte{MagicByte1, MagicByte2, MagicByte3}

// Header represents the Codec2 file format header
type Header struct {
	Magic        [3]byte
	VersionMajor byte
	VersionMinor byte
	Mode         byte
	Flags        byte
}

// NewC2Header creates a new C2 file header for 2400bps mode
func NewC2Header() []byte {
	return []byte{
		MagicByte1, // Magic bytes
		MagicByte2,
		MagicByte3,
		VersionMajor, // Version
		VersionMinor,
		Mode2400Byte, // Mode (2400bps)
		0x00,         // Flags
	}
}

// IsC2Header checks if the given bytes are a valid C2 file header
func IsC2Header(data []byte) bool {
	if len(data) < 7 {
		return false
	}
	return data[0] == MagicByte1 &&
		data[1] == MagicByte2 &&
		data[2] == MagicByte3
}

// IsValid checks if the header has valid magic bytes
func IsValid(h *Header) bool {
	return h.Magic[0] == Magic[0] && h.Magic[1] == Magic[1] && h.Magic[2] == Magic[2]
}

// GetC2Mode returns the mode from a C2 file header
func GetC2Mode(data []byte) int {
	if !IsC2Header(data) {
		return -1
	}
	return int(data[5])
}

// NewHeader creates a new header with the given mode
func NewHeader(mode byte) Header {
	h := Header{
		VersionMajor: 1,
		VersionMinor: 0,
		Mode:         mode,
	}
	copy(h.Magic[:], Magic)
	return h
}
