/*---------------------------------------------------------------------------*\

  FILE........: sine.h
  AUTHOR......: David Rowe
  DATE CREATED: 1/11/94

  Header file for sinusoidal analysis and synthesis functions.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2009 David Rowe

  All rights reserved.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License version 2.1, as
  published by the Free Software Foundation.  This program is
  distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __SINE__
#define __SINE__

#include "c2_fft.h"
#include "c2_comp.h"
#include "c2_defines.h"

C2CONST c2const_create(int Fs, float framelength_ms);

void make_analysis_window(C2CONST *c2const, codec2_fft_cfg fft_fwd_cfg,
                          float w[], float W[]);
float hpf(float x, float states[]);
void dft_speech(C2CONST *c2const, codec2_fft_cfg fft_fwd_cfg, COMP Sw[],
                float Sn[], float w[]);
void two_stage_pitch_refinement(C2CONST *c2const, MODEL *model, COMP Sw[]);
void estimate_amplitudes(MODEL *model, COMP Sw[], float W[], int est_phase);
float est_voicing_mbe(C2CONST *c2const, MODEL *model, COMP Sw[], float W[]);
void make_synthesis_window(C2CONST *c2const, float Pn[]);
void synthesise(int n_samp, codec2_fftr_cfg fftr_inv_cfg, float Sn_[],
                MODEL *model, float Pn[], int shift);

#define CODEC2_RAND_MAX 32767
int codec2_rand(void);

#endif
