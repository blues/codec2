#!/bin/bash
set -x

# build
cd ..
cd c/c2enc
./build.sh
cd ../..
cd c/c2dec
./build.sh
cd ../..
cd go/c2enc
./build.sh
cd ../..
cd go/c2dec
./build.sh
cd ../..
cd test

# convert wav to raw PCM
ffmpeg -y -loglevel quiet -i 1-wav/male.wav -f s16le -ac 1 -ar 8000 2-raw/male.raw
ffmpeg -y -loglevel quiet -i 1-wav/female.wav -f s16le -ac 1 -ar 8000 2-raw/female.raw

# encode the raw using codec2
../c/c2enc/c2enc 2-raw/male.raw 3-c2/male-c.c2
../c/c2enc/c2enc 2-raw/female.raw 3-c2/female-c.c2
../go/c2enc/c2enc 2-raw/male.raw 3-c2/male-go.c2
../go/c2enc/c2enc 2-raw/female.raw 3-c2/female-go.c2

# decode back to raw using codec2
../c/c2dec/c2dec 3-c2/male-c.c2 4-raw/male-c.raw
../c/c2dec/c2dec 3-c2/female-c.c2 4-raw/female-c.raw
../go/c2dec/c2dec 3-c2/male-go.c2 4-raw/male-go.raw
../go/c2dec/c2dec 3-c2/female-go.c2 4-raw/female-go.raw

# convert back to wav so we can look at them using audacity
cp 1-wav/male.wav 5-wav
cp 1-wav/female.wav 5-wav
ffmpeg -y -loglevel quiet -f s16le -ar 8000 -ac 1 -i 4-raw/male-c.raw 5-wav/male-c.wav
ffmpeg -y -loglevel quiet -f s16le -ar 8000 -ac 1 -i 4-raw/female-c.raw 5-wav/female-c.wav
ffmpeg -y -loglevel quiet -f s16le -ar 8000 -ac 1 -i 4-raw/male-go.raw 5-wav/male-go.wav
ffmpeg -y -loglevel quiet -f s16le -ar 8000 -ac 1 -i 4-raw/female-go.raw 5-wav/female-go.wav




