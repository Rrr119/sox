/* libSoX effect: deess - dynamic de-esser with a band-pass side-chain
 *
 * A side-chain built from an RBJ biquad band-pass filter measures the
 * energy around the sibilance band.  When the side-chain level rises
 * above the threshold, a downward gain (limited to a maximum
 * attenuation) is applied to the full-band signal.  The envelope
 * follower and the gain are each smoothed with attack/release
 * one-pole integrators so that the gain changes are free of the
 * audible "pumping" of a hard-switched limiter.
 *
 * Copyright 2024 Chris Bagwell And SoX Contributors
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "sox_i.h"
#include <math.h>

typedef struct {
  /* Options (kept in the units given on the command line). */
  double threshold_db;    /* side-chain level at which gain reduction starts */
  double center_freq;     /* side-chain band-pass centre frequency (Hz) */
  double bandwidth;       /* side-chain band-pass bandwidth (Hz) */
  double attack;          /* envelope/gain attack time (seconds) */
  double release;         /* envelope/gain release time (seconds) */
  double max_atten_db;    /* maximum gain reduction (dB, positive number) */

  /* Biquad band-pass coefficients (normalised so that a0 == 1). */
  double b0, b1, b2, a1, a2;

  /* Per-channel state. */
  size_t num_channels;
  struct {
    double x1, x2;        /* band-pass filter input history */
    double y1, y2;        /* band-pass filter output history */
    double envelope;      /* smoothed side-chain level (linear amplitude) */
    double gain;          /* smoothed full-band gain (linear, <= 1) */
  } * channel;

  double rate;            /* sample rate the coefficients were built for */
} priv_t;

#define DEESS_USAGE \
  "[threshold-dB [center-freq [bandwidth [attack [release [max-atten-dB]]]]]]\n" \
  "\tdefaults: -40dB, 6500Hz, 5000Hz, 0.005s, 0.05s, 20dB"

/* (Re)build the RBJ constant-0dB-peak-gain band-pass coefficients.
 * Called from start() and again from flow() if the sample rate changes. */
static int setup(sox_effect_t * effp)
{
  priv_t * p = (priv_t *)effp->priv;
  double rate = effp->in_signal.rate;
  double w0, alpha, a0;

  if (p->center_freq >= rate / 2) {
    lsx_fail("center frequency must be less than half the sample-rate (Nyquist rate)");
    return SOX_EOF;
  }

  w0 = 2 * M_PI * p->center_freq / rate;
  /* RBJ band-pass, bandwidth specified in Hz: Q = centre / bandwidth */
  alpha = sin(w0) / (2 * p->center_freq / p->bandwidth);
  a0 = 1 + alpha;

  p->b0 =  alpha / a0;
  p->b1 =  0;
  p->b2 = -alpha / a0;
  p->a1 = -2 * cos(w0) / a0;
  p->a2 = (1 - alpha) / a0;

  p->rate = rate;
  return SOX_SUCCESS;
}

static int getopts(sox_effect_t * effp, int argc, char **argv)
{
  priv_t * p = (priv_t *)effp->priv;
  char dummy;
  char * end_ptr;

  p->threshold_db = -40;
  p->center_freq  = 6500;
  p->bandwidth    = 5000;
  p->attack       = 0.005;
  p->release      = 0.05;
  p->max_atten_db = 20;

  --argc, ++argv;
  if (argc > 6)
    return lsx_usage(effp);

  if (argc > 0 && sscanf(argv[0], "%lf %c", &p->threshold_db, &dummy) != 1)
    return lsx_usage(effp);
  if (argc > 1 &&
      ((p->center_freq = lsx_parse_frequency(argv[1], &end_ptr)) <= 0 ||
       *end_ptr))
    return lsx_usage(effp);
  if (argc > 2 &&
      ((p->bandwidth = lsx_parse_frequency(argv[2], &end_ptr)) <= 0 ||
       *end_ptr))
    return lsx_usage(effp);
  if (argc > 3 && (sscanf(argv[3], "%lf %c", &p->attack, &dummy) != 1 ||
      p->attack < 0))
    return lsx_usage(effp);
  if (argc > 4 && (sscanf(argv[4], "%lf %c", &p->release, &dummy) != 1 ||
      p->release < 0))
    return lsx_usage(effp);
  if (argc > 5 && (sscanf(argv[5], "%lf %c", &p->max_atten_db, &dummy) != 1 ||
      p->max_atten_db < 0))
    return lsx_usage(effp);

  if (p->threshold_db > 0) {
    lsx_fail("threshold must be at or below 0 dB (full-scale)");
    return SOX_EOF;
  }
  return SOX_SUCCESS;
}

static int start(sox_effect_t * effp)
{
  priv_t * p = (priv_t *)effp->priv;
  size_t i;

  if (effp->in_signal.rate <= 0)
    return SOX_EOF;
  if (setup(effp) != SOX_SUCCESS)
    return SOX_EOF;

  p->num_channels = effp->in_signal.channels;
  p->channel = lsx_calloc(p->num_channels, sizeof(*p->channel));
  for (i = 0; i < p->num_channels; ++i)
    p->channel[i].gain = 1;
  return SOX_SUCCESS;
}

static int flow(sox_effect_t * effp, const sox_sample_t *ibuf,
    sox_sample_t *obuf, size_t *isamp, size_t *osamp)
{
  priv_t * p = (priv_t *)effp->priv;
  size_t channels = effp->out_signal.channels;
  size_t frames = min(*isamp, *osamp) / channels;
  double const threshold = dB_to_linear(p->threshold_db);
  double const min_gain = dB_to_linear(-p->max_atten_db);
  double const attack_coef = p->attack ?
      exp(-1 / (p->attack * effp->in_signal.rate)) : 0;
  double const release_coef = p->release ?
      exp(-1 / (p->release * effp->in_signal.rate)) : 0;
  size_t frame, ch;

  /* Coefficients depend on the sample rate; rebuild it if that changed. */
  if (effp->in_signal.rate != p->rate && setup(effp) != SOX_SUCCESS)
    return SOX_EOF;

  for (frame = 0; frame < frames; ++frame) {
    for (ch = 0; ch < channels; ++ch) {
      double in = SOX_SAMPLE_TO_FLOAT_64BIT(ibuf[ch], effp->clips);
      double side, level, over_db, target, coef;

      /* 1. Side-chain: RBJ biquad band-pass around the sibilance band. */
      side = in * p->b0 + p->channel[ch].x1 * p->b1 +
          p->channel[ch].x2 * p->b2 - p->channel[ch].y1 * p->a1 -
          p->channel[ch].y2 * p->a2;
      p->channel[ch].x2 = p->channel[ch].x1;
      p->channel[ch].x1 = in;
      p->channel[ch].y2 = p->channel[ch].y1;
      p->channel[ch].y1 = side;

      /* 2. Peak envelope follower with separate attack/release times. */
      level = fabs(side);
      coef = level > p->channel[ch].envelope ?
          attack_coef : release_coef;
      p->channel[ch].envelope =
          coef * p->channel[ch].envelope + (1 - coef) * level;

      /* 3. Downward gain above threshold, capped at the maximum
       *    attenuation; unity gain below the threshold. */
      if (p->channel[ch].envelope > threshold)
        over_db = 20 * log10(p->channel[ch].envelope / threshold);
      else
        over_db = 0;
      target = dB_to_linear(-min(over_db, p->max_atten_db));

      /* 4. Smooth the gain itself so the dynamic attenuation cannot
       *    pump or zipper; recovery (release) is slower than ducking. */
      coef = target < p->channel[ch].gain ? attack_coef : release_coef;
      p->channel[ch].gain =
          coef * p->channel[ch].gain + (1 - coef) * target;
      if (p->channel[ch].gain < min_gain)
        p->channel[ch].gain = min_gain;

      *obuf++ = SOX_ROUND_CLIP_COUNT(ibuf[ch] * p->channel[ch].gain,
          effp->clips);
    }
    ibuf += channels;
  }

  *isamp = *osamp = frames * channels;
  return SOX_SUCCESS;
}

static int stop(sox_effect_t * effp)
{
  priv_t * p = (priv_t *)effp->priv;
  free(p->channel);
  p->channel = NULL;
  return SOX_SUCCESS;
}

sox_effect_handler_t const * lsx_deess_effect_fn(void)
{
  static sox_effect_handler_t handler = {
    "deess", DEESS_USAGE, SOX_EFF_MCHAN | SOX_EFF_GAIN,
    getopts, start, flow, NULL, stop, NULL, sizeof(priv_t)
  };
  return &handler;
}
