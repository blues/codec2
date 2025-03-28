/*---------------------------------------------------------------------------*\

  FILE........: codec2.c
  AUTHOR......: David Rowe
  DATE CREATED: 21/8/2010

  Codec2 fully quantised encoder and decoder functions.  If you want use
  codec2, the codec2_xxx functions are for you.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2010 David Rowe

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

#include "c2_codec2.h"

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "c2_bpf.h"
#include "c2_bpfb.h"
#include "c2_fft.h"
#include "c2_internal.h"
#include "c2_alloc.h"
#include "c2_defines.h"
#include "c2_interp.h"
#include "c2_lpc.h"
#include "c2_lsp.h"
#include "c2_machdep.h"
#include "c2_nlp.h"
#include "c2_phase.h"
#include "c2_postfilter.h"
#include "c2_quantise.h"
#include "c2_sine.h"

/*---------------------------------------------------------------------------* \

                             FUNCTION HEADERS

\*---------------------------------------------------------------------------*/

void analyse_one_frame(struct CODEC2 *c2, MODEL *model, short speech[]);
void synthesise_one_frame(struct CODEC2 *c2, short speech[], MODEL *model,
                          COMP Aw[], float gain);
void codec2_encode_2400(struct CODEC2 *c2, unsigned char *bits, short speech[]);
void codec2_decode_2400(struct CODEC2 *c2, short speech[],
                        const unsigned char *bits);
static void ear_protection(float in_out[], int n);

/*---------------------------------------------------------------------------*\

                                FUNCTIONS

\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\

  FUNCTION....: codec2_create
  AUTHOR......: David Rowe
  DATE CREATED: 21/8/2010

  Create and initialise an instance of the codec.  Returns a pointer
  to the codec states or NULL on failure.  One set of states is
  sufficient for a full duuplex codec (i.e. an encoder and decoder).
  You don't need separate states for encoders and decoders.  See
  c2enc.c and c2dec.c for examples.

\*---------------------------------------------------------------------------*/

struct CODEC2 *codec2_create(int mode)
{
    struct CODEC2 *c2;
    int i, l;

    // ALL POSSIBLE MODES MUST BE CHECKED HERE!
    // we test if the desired mode is enabled at compile time
    // and return NULL if not

    if (false == (CODEC2_MODE_ACTIVE(CODEC2_MODE_2400, mode))) {
        return NULL;
    }

    c2 = (struct CODEC2 *)MALLOC(sizeof(struct CODEC2));
    if (c2 == NULL) {
        return NULL;
    }

    c2->mode = mode;

    /* store constants in a few places for convenience */

    c2->c2const = c2const_create(8000, N_S);
    c2->Fs = c2->c2const.Fs;
    int n_samp = c2->n_samp = c2->c2const.n_samp;
    int m_pitch = c2->m_pitch = c2->c2const.m_pitch;

    c2->Pn = (float *)MALLOC(2 * n_samp * sizeof(float));
    if (c2->Pn == NULL) {
        return NULL;
    }
    c2->Sn_ = (float *)MALLOC(2 * n_samp * sizeof(float));
    if (c2->Sn_ == NULL) {
        FREE(c2->Pn);
        return NULL;
    }
    c2->w = (float *)MALLOC(m_pitch * sizeof(float));
    if (c2->w == NULL) {
        FREE(c2->Pn);
        FREE(c2->Sn_);
        return NULL;
    }
    c2->Sn = (float *)MALLOC(m_pitch * sizeof(float));
    if (c2->Sn == NULL) {
        FREE(c2->Pn);
        FREE(c2->Sn_);
        FREE(c2->w);
        return NULL;
    }

    for (i = 0; i < m_pitch; i++) {
        c2->Sn[i] = 1.0;
    }
    c2->hpf_states[0] = c2->hpf_states[1] = 0.0;
    for (i = 0; i < 2 * n_samp; i++) {
        c2->Sn_[i] = 0;
    }
    c2->fft_fwd_cfg = codec2_fft_alloc(FFT_ENC, 0, NULL, NULL);
    c2->fftr_fwd_cfg = codec2_fftr_alloc(FFT_ENC, 0, NULL, NULL);
    make_analysis_window(&c2->c2const, c2->fft_fwd_cfg, c2->w, c2->W);
    make_synthesis_window(&c2->c2const, c2->Pn);
    c2->fftr_inv_cfg = codec2_fftr_alloc(FFT_DEC, 1, NULL, NULL);
    c2->prev_f0_enc = 1 / P_MAX_S;
    c2->bg_est = 0.0;
    c2->ex_phase = 0.0;

    for (l = 1; l <= MAX_AMP; l++) {
        c2->prev_model_dec.A[l] = 0.0;
    }
    c2->prev_model_dec.Wo = TWO_PI / c2->c2const.p_max;
    c2->prev_model_dec.L = (int) (PI / c2->prev_model_dec.Wo);
    c2->prev_model_dec.voiced = 0;

    for (i = 0; i < LPC_ORD; i++) {
        c2->prev_lsps_dec[i] = i * PI / (LPC_ORD + 1);
    }
    c2->prev_e_dec = 1;

    c2->nlp = nlp_create(&c2->c2const);
    if (c2->nlp == NULL) {
        return NULL;
    }

    c2->lpc_pf = 1;
    c2->bass_boost = 1;
    c2->beta = LPCPF_BETA;
    c2->gamma = LPCPF_GAMMA;

    c2->xq_enc[0] = c2->xq_enc[1] = 0.0;
    c2->xq_dec[0] = c2->xq_dec[1] = 0.0;

    c2->smoothing = 0;
    c2->se = 0.0;
    c2->nse = 0;
    c2->user_rate_K_vec_no_mean_ = NULL;
    c2->post_filter_en = true;

    c2->bpf_buf = (float *)MALLOC(sizeof(float) * (BPF_N + 4 * c2->n_samp));
    assert(c2->bpf_buf != NULL);
    for (i = 0; i < BPF_N + 4 * c2->n_samp; i++) {
        c2->bpf_buf[i] = 0.0;
    }

    c2->softdec = NULL;
    c2->gray = 1;

    /* newamp1 initialisation */

    c2->fmlfeat = NULL;
    c2->fmlmodel = NULL;

    // make sure that one of the two decode function pointers is empty
    // for the encode function pointer this is not required since we always set it
    // to a meaningful value

    c2->decode = NULL;
    c2->decode_ber = NULL;

    if (CODEC2_MODE_ACTIVE(CODEC2_MODE_2400, c2->mode)) {
        c2->encode = codec2_encode_2400;
        c2->decode = codec2_decode_2400;
    }

    return c2;
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: codec2_destroy
  AUTHOR......: David Rowe
  DATE CREATED: 21/8/2010

  Destroy an instance of the codec.

\*---------------------------------------------------------------------------*/

void codec2_destroy(struct CODEC2 *c2)
{
    assert(c2 != NULL);
    FREE(c2->bpf_buf);
    nlp_destroy(c2->nlp);
    codec2_fft_free(c2->fft_fwd_cfg);
    codec2_fftr_free(c2->fftr_fwd_cfg);
    codec2_fftr_free(c2->fftr_inv_cfg);
    FREE(c2->Pn);
    FREE(c2->Sn);
    FREE(c2->w);
    FREE(c2->Sn_);
    FREE(c2);
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: codec2_bits_per_frame
  AUTHOR......: David Rowe
  DATE CREATED: Nov 14 2011

  Returns the number of bits per frame.

\*---------------------------------------------------------------------------*/

int codec2_bits_per_frame(struct CODEC2 *c2)
{
    if (CODEC2_MODE_ACTIVE(CODEC2_MODE_2400, c2->mode)) {
        return 48;
    }
    return 0; /* shouldn't get here */
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: codec2_bytes_per_frame
  DATE CREATED: April 2021

  Returns the number of bytes per frame.  Useful for allocated storage for
  codec2_encode()/codec2_decode().  Note the number of bits may not be a
  multiple of 8, therefore some bits in the last byte may be unused.

\*---------------------------------------------------------------------------*/

int codec2_bytes_per_frame(struct CODEC2 *c2)
{
    return (codec2_bits_per_frame(c2) + 7) / 8;
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: codec2_samples_per_frame
  AUTHOR......: David Rowe
  DATE CREATED: Nov 14 2011

  Returns the number of speech samples per frame.

\*---------------------------------------------------------------------------*/

int codec2_samples_per_frame(struct CODEC2 *c2)
{
    if (CODEC2_MODE_ACTIVE(CODEC2_MODE_2400, c2->mode)) {
        return 160;
    }
    return 0; /* shouldn't get here */
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: codec2_encode
  AUTHOR......: David Rowe
  DATE CREATED: Nov 14 2011

  Take an input buffer of speech samples, and compress them to a packed buffer
  of bytes.

\*---------------------------------------------------------------------------*/

void codec2_encode(struct CODEC2 *c2, unsigned char *bytes, short speech[])
{
    assert(c2 != NULL);
    assert(c2->encode != NULL);

    c2->encode(c2, bytes, speech);
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: codec2_decode
  AUTHOR......: David Rowe
  DATE CREATED: Nov 14 2011

  Take an input packed buffer of bytes, and decode them to a buffer of speech
  samples.

\*---------------------------------------------------------------------------*/

void codec2_decode(struct CODEC2 *c2, short speech[],
                   const unsigned char *bytes)
{
    codec2_decode_ber(c2, speech, bytes, 0.0);
}

void codec2_decode_ber(struct CODEC2 *c2, short speech[],
                       const unsigned char *bits, float ber_est)
{
    assert(c2 != NULL);
    assert(c2->decode != NULL || c2->decode_ber != NULL);

    if (c2->decode != NULL) {
        c2->decode(c2, speech, bits);
    } else {
        c2->decode_ber(c2, speech, bits, ber_est);
    }
}


/*---------------------------------------------------------------------------*\

  FUNCTION....: codec2_encode_2400
  AUTHOR......: David Rowe
  DATE CREATED: 21/8/2010

  Encodes 160 speech samples (20ms of speech) into 48 bits.

  The codec2 algorithm actually operates internally on 10ms (80
  sample) frames, so we run the encoding algorithm twice.  On the
  first frame we just send the voicing bit.  On the second frame we
  send all model parameters.

  The bit allocation is:

    Parameter                      bits/frame
    --------------------------------------
    Harmonic magnitudes (LSPs)     36
    Joint VQ of Energy and Wo       8
    Voicing (10ms update)           2
    Spare                           2
    TOTAL                          48

\*---------------------------------------------------------------------------*/

void codec2_encode_2400(struct CODEC2 *c2, unsigned char *bits,
                        short speech[])
{
    MODEL model;
    float ak[LPC_ORD + 1];
    float lsps[LPC_ORD];
    float e;
    int WoE_index;
    int lsp_indexes[LPC_ORD];
    int i;
    int spare = 0;
    unsigned int nbit = 0;

    assert(c2 != NULL);

    memset(bits, '\0', ((codec2_bits_per_frame(c2) + 7) / 8));

    /* first 10ms analysis frame - we just want voicing */

    analyse_one_frame(c2, &model, speech);
    pack(bits, &nbit, model.voiced, 1);

    /* second 10ms analysis frame */

    analyse_one_frame(c2, &model, &speech[c2->n_samp]);
    pack(bits, &nbit, model.voiced, 1);

    e = speech_to_uq_lsps(lsps, ak, c2->Sn, c2->w, c2->m_pitch, LPC_ORD);
    WoE_index = encode_WoE(&model, e, c2->xq_enc);
    pack(bits, &nbit, WoE_index, WO_E_BITS);

    encode_lsps_scalar(lsp_indexes, lsps, LPC_ORD);
    for (i = 0; i < LSP_SCALAR_INDEXES; i++) {
        pack(bits, &nbit, lsp_indexes[i], lsp_bits(i));
    }
    pack(bits, &nbit, spare, 2);

    assert(nbit == (unsigned)codec2_bits_per_frame(c2));
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: codec2_decode_2400
  AUTHOR......: David Rowe
  DATE CREATED: 21/8/2010

  Decodes frames of 48 bits into 160 samples (20ms) of speech.

\*---------------------------------------------------------------------------*/

void codec2_decode_2400(struct CODEC2 *c2, short speech[],
                        const unsigned char *bits)
{
    MODEL model[2];
    int lsp_indexes[LPC_ORD];
    float lsps[2][LPC_ORD];
    int WoE_index;
    float e[2];
    float snr;
    float ak[2][LPC_ORD + 1];
    int i, j;
    unsigned int nbit = 0;
    COMP Aw[FFT_ENC];

    assert(c2 != NULL);

    /* only need to zero these out due to (unused) snr calculation */

    for (i = 0; i < 2; i++)
        for (j = 1; j <= MAX_AMP; j++) {
            model[i].A[j] = 0.0;
        }

    /* unpack bits from channel ------------------------------------*/

    /* this will partially fill the model params for the 2 x 10ms
       frames */

    model[0].voiced = unpack(bits, &nbit, 1);
    model[1].voiced = unpack(bits, &nbit, 1);
    WoE_index = unpack(bits, &nbit, WO_E_BITS);
    decode_WoE(&c2->c2const, &model[1], &e[1], c2->xq_dec, WoE_index);

    for (i = 0; i < LSP_SCALAR_INDEXES; i++) {
        lsp_indexes[i] = unpack(bits, &nbit, lsp_bits(i));
    }
    decode_lsps_scalar(&lsps[1][0], lsp_indexes, LPC_ORD);
    check_lsp_order(&lsps[1][0], LPC_ORD);
    bw_expand_lsps(&lsps[1][0], LPC_ORD, 50.0, 100.0);

    /* interpolate ------------------------------------------------*/

    /* Wo and energy are sampled every 20ms, so we interpolate just 1
       10ms frame between 20ms samples */

    interp_Wo(&model[0], &c2->prev_model_dec, &model[1], c2->c2const.Wo_min);
    e[0] = interp_energy(c2->prev_e_dec, e[1]);

    /* LSPs are sampled every 20ms so we interpolate the frame in
       between, then recover spectral amplitudes */

    interpolate_lsp_ver2(&lsps[0][0], c2->prev_lsps_dec, &lsps[1][0], 0.5,
                         LPC_ORD);
    for (i = 0; i < 2; i++) {
        lsp_to_lpc(&lsps[i][0], &ak[i][0], LPC_ORD);
        aks_to_M2(c2->fftr_fwd_cfg, &ak[i][0], LPC_ORD, &model[i], e[i], &snr, 0, 0,
                  c2->lpc_pf, c2->bass_boost, c2->beta, c2->gamma, Aw);
        apply_lpc_correction(&model[i]);
        synthesise_one_frame(c2, &speech[c2->n_samp * i], &model[i], Aw, 1.0);

    }

    /* update memories for next frame ----------------------------*/

    c2->prev_model_dec = model[1];
    c2->prev_e_dec = e[1];
    for (i = 0; i < LPC_ORD; i++) {
        c2->prev_lsps_dec[i] = lsps[1][i];
    }
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: codec2_get_energy()
  AUTHOR......: Jeroen Vreeken
  DATE CREATED: 08/03/2016

  Extract energy value from an encoded frame.

\*---------------------------------------------------------------------------*/

float codec2_get_energy(struct CODEC2 *c2, const unsigned char *bits)
{
    assert(c2 != NULL);
    assert((CODEC2_MODE_ACTIVE(CODEC2_MODE_2400, c2->mode)));
    MODEL model;
    float xq_dec[2] = {0};
    int WoE_index;
    float e = 0.0f;
    unsigned int nbit;

    if (CODEC2_MODE_ACTIVE(CODEC2_MODE_2400, c2->mode)) {
        nbit = 1 + 1;
        WoE_index = unpack(bits, &nbit, WO_E_BITS);
        decode_WoE(&c2->c2const, &model, &e, xq_dec, WoE_index);
    }

    return e;
}

/*---------------------------------------------------------------------------* \

  FUNCTION....: synthesise_one_frame()
  AUTHOR......: David Rowe
  DATE CREATED: 23/8/2010

  Synthesise 80 speech samples (10ms) from model parameters.

\*---------------------------------------------------------------------------*/

void synthesise_one_frame(struct CODEC2 *c2, short speech[], MODEL *model,
                          COMP Aw[], float gain)
{
    int i;

    /* LPC based phase synthesis */
    COMP H[MAX_AMP + 1];
    sample_phase(model, H, Aw);
    phase_synth_zero_order(c2->n_samp, model, &c2->ex_phase, H);

    postfilter(model, &c2->bg_est);
    synthesise(c2->n_samp, c2->fftr_inv_cfg, c2->Sn_, model, c2->Pn, 1);

    for (i = 0; i < c2->n_samp; i++) {
        c2->Sn_[i] *= gain;
    }

    ear_protection(c2->Sn_, c2->n_samp);

    for (i = 0; i < c2->n_samp; i++) {
        if (c2->Sn_[i] > 32767.0) {
            speech[i] = 32767;
        } else if (c2->Sn_[i] < -32767.0) {
            speech[i] = -32767;
        } else {
            speech[i] = (int) c2->Sn_[i];
        }
    }
}

/*---------------------------------------------------------------------------* \

  FUNCTION....: analyse_one_frame()
  AUTHOR......: David Rowe
  DATE CREATED: 23/8/2010

  Extract sinusoidal model parameters from 80 speech samples (10ms of
  speech).

\*---------------------------------------------------------------------------*/

void analyse_one_frame(struct CODEC2 *c2, MODEL *model, short speech[])
{
    COMP Sw[FFT_ENC];
    float pitch;
    int i;
    int n_samp = c2->n_samp;
    int m_pitch = c2->m_pitch;

    /* Read input speech */

    for (i = 0; i < m_pitch - n_samp; i++) {
        c2->Sn[i] = c2->Sn[i + n_samp];
    }
    for (i = 0; i < n_samp; i++) {
        c2->Sn[i + m_pitch - n_samp] = speech[i];
    }

    dft_speech(&c2->c2const, c2->fft_fwd_cfg, Sw, c2->Sn, c2->w);

    /* Estimate pitch */
    nlp(c2->nlp, c2->Sn, n_samp, &pitch, Sw, c2->W, &c2->prev_f0_enc);
    model->Wo = TWO_PI / pitch;
    model->L = (int) (PI / model->Wo);

    /* estimate model parameters */
    two_stage_pitch_refinement(&c2->c2const, model, Sw);

    /* estimate phases when doing ML experiments */
    if (c2->fmlfeat != NULL) {
        estimate_amplitudes(model, Sw, c2->W, 1);
    } else {
        estimate_amplitudes(model, Sw, c2->W, 0);
    }
    est_voicing_mbe(&c2->c2const, model, Sw, c2->W);
}

/*---------------------------------------------------------------------------* \

  FUNCTION....: ear_protection()
  AUTHOR......: David Rowe
  DATE CREATED: Nov 7 2012

  Limits output level to protect ears when there are bit errors or the input
  is overdriven.  This doesn't correct or mask bit errors, just reduces the
  worst of their damage.

\*---------------------------------------------------------------------------*/

static void ear_protection(float in_out[], int n)
{
    float max_sample, over, gain;
    int i;

    /* find maximum sample in frame */

    max_sample = 0.0;
    for (i = 0; i < n; i++)
        if (in_out[i] > max_sample) {
            max_sample = in_out[i];
        }

    /* determine how far above set point */

    over = max_sample / 30000.0;

    /* If we are x dB over set point we reduce level by 2x dB, this
       attenuates major excursions in amplitude (likely to be caused
       by bit errors) more than smaller ones */

    if (over > 1.0) {
        gain = 1.0 / (over * over);
        for (i = 0; i < n; i++) {
            in_out[i] *= gain;
        }
    }
}

void codec2_set_lpc_post_filter(struct CODEC2 *c2, int enable, int bass_boost,
                                float beta, float gamma)
{
    assert((beta >= 0.0) && (beta <= 1.0));
    assert((gamma >= 0.0) && (gamma <= 1.0));
    c2->lpc_pf = enable;
    c2->bass_boost = bass_boost;
    c2->beta = beta;
    c2->gamma = gamma;
}

/*
   Allows optional stealing of one of the voicing bits for use as a
   spare bit, only 1300 & 1400 & 1600 bit/s supported for now.
   Experimental method of sending voice/data frames for FreeDV.
*/

int codec2_get_spare_bit_index(struct CODEC2 *c2)
{
    assert(c2 != NULL);
	(void) c2;
    return -1;
}

void codec2_set_natural_or_gray(struct CODEC2 *c2, int gray)
{
    assert(c2 != NULL);
    c2->gray = gray;
}

void codec2_set_softdec(struct CODEC2 *c2, float *softdec)
{
    assert(c2 != NULL);
    c2->softdec = softdec;
}

void codec2_open_mlfeat(struct CODEC2 *codec2_state, char *feat_fn,
                        char *model_fn)
{
    if ((codec2_state->fmlfeat = fopen(feat_fn, "wb")) == NULL) {
        fprintf(stderr, "error opening machine learning feature file: %s\n",
                feat_fn);
        exit(1);
    }
    if (model_fn) {
        if ((codec2_state->fmlmodel = fopen(model_fn, "wb")) == NULL) {
            fprintf(stderr, "error opening machine learning Codec 2 model file: %s\n",
                    feat_fn);
            exit(1);
        }
    }
}

float codec2_get_var(struct CODEC2 *codec2_state)
{
    if (codec2_state->nse) {
        return codec2_state->se / codec2_state->nse;
    } else {
        return 0;
    }
}

float *codec2_enable_user_ratek(struct CODEC2 *codec2_state, int *K)
{
    codec2_state->user_rate_K_vec_no_mean_ =
        (float *)MALLOC(sizeof(float) * NEWAMP1_K);
    *K = NEWAMP1_K;
    return codec2_state->user_rate_K_vec_no_mean_;
}

