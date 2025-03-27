#!/bin/bash

gcc -o c2enc -I. -I.. \
	../c2_codec2.c \
	../c2_kiss_fft.c \
	../c2_kiss_fftr.c \
	../c2_lpc.c \
	../c2_lsp.c \
	../c2_nlp.c \
	../c2_sine.c \
	../c2_phase.c \
	../c2_quantise.c \
	../c2_postfilter.c \
	../c2_interp.c \
	../c2_pack.c \
	../c2_fft.c \
	../c2_mbest.c \
	../c2_lsp_cb.c \
	../c2_ge_cb.c \
	-x c main.c.skip

	
