# Codec 2 README

Codec 2 is an open source (LGPL 2.1) low bit rate speech codec: http://rowetel.com/codec2.html derived directly from David Rowe's C implementation:
https://github.com/drowe67/codec2
https://www.rowetel.com/

It is extremely kind and generous of David Rowe to contribute this to the community, and we are extremely thankful for his skill and dedication.

The implementation in the repo herein is a subset of his implementation.  There is a singlular purpose of this fork, which is
a) to have the smallest possible encoder/decoder on a low-power microcontroller with a PDM mic and a DFSDM subsystem so as to convert the PDM to PCM.  This would be wrapped in a "push to talk" semantic.
b) to encode the PCM audio at as low a bitrate as possible so that voice could be transmitted over a highly-constrained transport such as cellular or satellite
c) to then, on the server, use Whisper to convert the voice to text
d) to communicate or process the text
e) to send the text reply back to the device
f) in a future implementation, to send audio data back to the device for an audio reply rather than text.

The reason 2400bps was selected was that, after trial and error, it was observed that this is the lowest practical CODEC2 bitrate that preserves the ability to use Whisper accurately.
