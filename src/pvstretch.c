/* libSoX effect: phase-vocoder time stretching with transient phase reset.
 *
 * Changes duration without changing pitch.  The STFT uses a Hann window and
 * advances the synthesis phase independently of the analysis phase.  A
 * transient detector resets those accumulators so percussive attacks stay
 * sharp instead of being smeared by overlapping windows.
 */

#include "sox_i.h"
#include "fft4g.h"
#include "fifo.h"

#include <assert.h>
#include <math.h>
#include <string.h>

#define DEFAULT_WINDOW_MS 46.0
#define DEFAULT_OVERLAP_PERCENT 75.0
#define DEFAULT_TRANSIENT_SENSITIVITY 0.5
#define MIN_FFT_SIZE 64
#define MAX_FFT_SIZE 32768

typedef struct {
  size_t channels;
  double factor;
  size_t fft_size;
  size_t bins;
  size_t analysis_hop;
  double synthesis_step;

  double * window;
  double * window_sq;
  double * fft_work;
  double * magnitude;
  double * combined;
  double * phase;
  double * prev_phase;
  double * sum_phase;
  double * theta;
  size_t * ref_chan;
  size_t * peak_index;
  double * dc;
  double * nyquist;

  float * ola;
  double * ola_norm;
  fifo_t input_fifo;
  fifo_t output_fifo;

  uint64_t samples_in;
  uint64_t samples_out;
  uint64_t frames;
  uint64_t out_base;
  uint64_t out_generated;
  uint64_t expected_out;
  sox_bool expected_known;
  sox_bool flushed;
  sox_bool have_previous_frame;

  double previous_hop_energy;
  unsigned transient_cooldown;
  double transient_sensitivity;
} pvstretch_t;

static double wrap_phase(double phase)
{
  return phase - 2 * M_PI * floor(phase / (2 * M_PI) + 0.5);
}

static size_t nearest_power_of_2(double x)
{
  size_t lower = MIN_FFT_SIZE, upper;

  while ((double)lower * 2 <= x && lower < MAX_FFT_SIZE)
    lower *= 2;
  upper = min(lower * 2, (size_t)MAX_FFT_SIZE);
  if (lower >= MAX_FFT_SIZE || x - lower <= (double)upper - x)
    return lower;
  return upper;
}

static uint64_t rounded_position(uint64_t frame, double step)
{
  return (uint64_t)floor(frame * step + 0.5);
}

static void pvstretch_delete(pvstretch_t * p)
{
  if (!p)
    return;
  free(p->window);
  free(p->window_sq);
  free(p->fft_work);
  free(p->magnitude);
  free(p->combined);
  free(p->phase);
  free(p->prev_phase);
  free(p->sum_phase);
  free(p->theta);
  free(p->ref_chan);
  free(p->peak_index);
  free(p->dc);
  free(p->nyquist);
  free(p->ola);
  free(p->ola_norm);
  fifo_delete(&p->input_fifo);
  fifo_delete(&p->output_fifo);
  free(p);
}

static pvstretch_t * pvstretch_create(size_t channels, double factor,
    double sample_rate, double window_ms, double overlap_percent,
    double sensitivity)
{
  pvstretch_t * p = lsx_calloc(1, sizeof(*p));
  size_t i;
  double requested = sample_rate * window_ms / 1000;
  double overlap = overlap_percent / 100;
  double max_analysis_hop;

  p->channels = channels;
  p->factor = factor;
  p->transient_sensitivity = sensitivity;
  p->fft_size = nearest_power_of_2(requested);
  p->bins = p->fft_size / 2 + 1;
  overlap = range_limit(overlap, .5, .9);
  p->analysis_hop = (size_t)floor(p->fft_size * (1 - overlap) + 0.5);
  p->analysis_hop = max(p->analysis_hop, (size_t)1);

  max_analysis_hop = factor * p->fft_size / 2;
  if (p->analysis_hop > max_analysis_hop)
    p->analysis_hop = max((size_t)floor(max_analysis_hop), (size_t)1);
  p->synthesis_step = p->analysis_hop / factor;

  p->window = lsx_malloc(p->fft_size * sizeof(*p->window));
  p->window_sq = lsx_malloc(p->fft_size * sizeof(*p->window_sq));
  for (i = 0; i < p->fft_size; ++i) {
    p->window[i] = .5 - .5 * cos(2 * M_PI * i / p->fft_size);
    p->window_sq[i] = p->window[i] * p->window[i];
  }

  p->fft_work = lsx_calloc(channels * p->fft_size, sizeof(*p->fft_work));
  p->magnitude = lsx_calloc(channels * p->bins, sizeof(*p->magnitude));
  p->combined = lsx_calloc(p->bins, sizeof(*p->combined));
  p->phase = lsx_calloc(channels * p->bins, sizeof(*p->phase));
  p->prev_phase = lsx_calloc(channels * p->bins, sizeof(*p->prev_phase));
  p->sum_phase = lsx_calloc(p->bins, sizeof(*p->sum_phase));
  p->theta = lsx_calloc(p->bins, sizeof(*p->theta));
  p->ref_chan = lsx_calloc(p->bins, sizeof(*p->ref_chan));
  p->peak_index = lsx_calloc(p->bins, sizeof(*p->peak_index));
  p->dc = lsx_calloc(channels, sizeof(*p->dc));
  p->nyquist = lsx_calloc(channels, sizeof(*p->nyquist));
  p->ola = lsx_calloc(channels * p->fft_size, sizeof(*p->ola));
  p->ola_norm = lsx_calloc(p->fft_size, sizeof(*p->ola_norm));
  fifo_create(&p->input_fifo, channels * sizeof(float));
  fifo_create(&p->output_fifo, channels * sizeof(float));
  return p;
}

static sox_bool detect_transient(pvstretch_t * p, float const * frame)
{
  double energy = 0, strength, threshold;
  size_t i, n = p->analysis_hop * p->channels;

  for (i = 0; i < n; ++i)
    energy += frame[i] * frame[i];
  energy /= n;
  strength = (energy - p->previous_hop_energy) /
      (energy + p->previous_hop_energy + 1e-20);
  threshold = .45 * (1 - p->transient_sensitivity) + .02;

  if (p->transient_cooldown)
    --p->transient_cooldown;

  if (p->have_previous_frame && !p->transient_cooldown &&
      energy > 1e-8 && strength > threshold) {
    p->previous_hop_energy = energy;
    p->transient_cooldown = max((size_t)1,
        p->fft_size / (2 * p->analysis_hop));
    return sox_true;
  }
  p->previous_hop_energy = energy;
  p->have_previous_frame = sox_true;
  return sox_false;
}

static void choose_phase_reference(pvstretch_t * p)
{
  size_t k, c;

  for (k = 1; k < p->bins - 1; ++k) {
    size_t best = 0;
    double best_mag = p->magnitude[k];
    for (c = 1; c < p->channels; ++c) {
      double mag = p->magnitude[c * p->bins + k];
      if (mag > best_mag) {
        best = c;
        best_mag = mag;
      }
    }
    p->ref_chan[k] = best;
    p->combined[k] = best_mag;
  }
}

static void advance_phases(pvstretch_t * p, sox_bool reset)
{
  size_t k, c;

  for (k = 1; k < p->bins - 1; ++k) {
    size_t ref = p->ref_chan[k];
    double current = p->phase[ref * p->bins + k];
    if (reset || !p->frames) {
      p->sum_phase[k] = current;
    } else {
      double omega = -2 * M_PI * k * p->analysis_hop / p->fft_size;
      double delta = wrap_phase(current - p->prev_phase[ref * p->bins + k] - omega);
      p->sum_phase[k] += (omega + delta) / p->factor;
    }
  }
  for (c = 0; c < p->channels; ++c)
    memcpy(p->prev_phase + c * p->bins + 1,
        p->phase + c * p->bins + 1,
        (p->bins - 2) * sizeof(*p->prev_phase));
}

static void lock_phases(pvstretch_t * p)
{
  size_t k, npeaks = 0, pi = 0;
  double max_mag = 0;
  size_t lock_radius = max(p->fft_size / 8, (size_t)2);

  for (k = 1; k < p->bins - 1; ++k) {
    max_mag = max(max_mag, p->combined[k]);
  }
  for (k = 2; k + 1 < p->bins - 1; ++k) {
    if (p->combined[k] > p->combined[k - 1] &&
        p->combined[k] >= p->combined[k + 1] &&
        p->combined[k] > max_mag * 1e-4)
      p->peak_index[npeaks++] = k;
  }

  if (!npeaks) {
    for (k = 1; k < p->bins - 1; ++k) {
      size_t ref = p->ref_chan[k];
      p->theta[k] = wrap_phase(p->sum_phase[k] - p->phase[ref * p->bins + k]);
    }
    return;
  }

  for (k = 1; k < p->bins - 1; ++k) {
    size_t nearest;
    while (pi + 1 < npeaks &&
        k > (p->peak_index[pi] + p->peak_index[pi + 1]) / 2)
      ++pi;
    nearest = p->peak_index[pi];
    if ((k > nearest ? k - nearest : nearest - k) <= lock_radius) {
      size_t ref = p->ref_chan[nearest];
      p->theta[k] =
          wrap_phase(p->sum_phase[nearest] - p->phase[ref * p->bins + nearest]);
    } else {
      size_t ref = p->ref_chan[k];
      p->theta[k] = wrap_phase(p->sum_phase[k] - p->phase[ref * p->bins + k]);
    }
  }
}

static void pvstretch_process_frame(pvstretch_t * p)
{
  float * frame = fifo_read_ptr(&p->input_fifo);
  uint64_t next_out = rounded_position(p->frames + 1, p->synthesis_step);
  size_t hop = (size_t)(next_out - p->out_base);
  size_t i, c, k;
  sox_bool transient;

  assert(fifo_occupancy(&p->input_fifo) >= p->fft_size);
  assert(p->out_base == rounded_position(p->frames, p->synthesis_step));
  assert(hop > 0 && hop <= p->fft_size);

  transient = detect_transient(p, frame);

  for (c = 0; c < p->channels; ++c) {
    double * work = p->fft_work + c * p->fft_size;
    for (i = 0; i < p->fft_size; ++i)
      work[i] = frame[i * p->channels + c] * p->window[i];
    lsx_safe_rdft(p->fft_size, 1, work);
    p->dc[c] = work[0];
    p->nyquist[c] = work[1];
    p->magnitude[c * p->bins] = fabs(work[0]);
    p->phase[c * p->bins] = 0;
    p->magnitude[c * p->bins + p->bins - 1] = fabs(work[1]);
    p->phase[c * p->bins + p->bins - 1] = 0;
    for (k = 1; k < p->bins - 1; ++k) {
      double re = work[2 * k];
      double im = work[2 * k + 1];
      p->magnitude[c * p->bins + k] = hypot(re, im);
      p->phase[c * p->bins + k] = atan2(im, re);
    }
  }

  choose_phase_reference(p);
  advance_phases(p, transient);
  if (transient || !p->frames) {
    memset(p->theta + 1, 0, (p->bins - 2) * sizeof(*p->theta));
  } else {
    lock_phases(p);
  }

  for (c = 0; c < p->channels; ++c) {
    double * work = p->fft_work + c * p->fft_size;
    for (k = 1; k < p->bins - 1; ++k) {
      double mag = p->magnitude[c * p->bins + k];
      double phase = p->phase[c * p->bins + k] + p->theta[k];
      work[2 * k] = mag * cos(phase);
      work[2 * k + 1] = mag * sin(phase);
    }
    work[0] = p->dc[c];
    work[1] = p->nyquist[c];
    lsx_safe_rdft(p->fft_size, -1, work);
  }

  for (i = 0; i < p->fft_size; ++i) {
    double w = p->window[i];
    p->ola_norm[i] += p->window_sq[i];
    for (c = 0; c < p->channels; ++c) {
      double sample = p->fft_work[c * p->fft_size + i] / p->fft_size * w;
      p->ola[c * p->fft_size + i] += (float)sample;
    }
  }

  {
    float * out = fifo_write(&p->output_fifo, hop, NULL);
    for (i = 0; i < hop; ++i) {
      double norm = p->ola_norm[i];
      for (c = 0; c < p->channels; ++c) {
        double sample = norm > 1e-12 ? p->ola[c * p->fft_size + i] / norm : 0;
        out[i * p->channels + c] = (float)sample;
      }
    }
  }

  for (c = 0; c < p->channels; ++c) {
    memmove(p->ola + c * p->fft_size, p->ola + c * p->fft_size + hop,
        (p->fft_size - hop) * sizeof(*p->ola));
    memset(p->ola + c * p->fft_size + p->fft_size - hop, 0,
        hop * sizeof(*p->ola));
  }
  memmove(p->ola_norm, p->ola_norm + hop,
      (p->fft_size - hop) * sizeof(*p->ola_norm));
  memset(p->ola_norm + p->fft_size - hop, 0,
      hop * sizeof(*p->ola_norm));

  fifo_read(&p->input_fifo, p->analysis_hop, NULL);
  p->out_base = next_out;
  p->out_generated += hop;
  ++p->frames;
}

static void pvstretch_process_available(pvstretch_t * p)
{
  while (fifo_occupancy(&p->input_fifo) >= p->fft_size &&
      (!p->expected_known || p->out_generated < p->expected_out))
    pvstretch_process_frame(p);
}

static void pvstretch_flush(pvstretch_t * p)
{
  size_t pad, c;

  if (!p->expected_known) {
    p->expected_out = (uint64_t)floor(p->samples_in / p->factor + 0.5);
    p->expected_known = sox_true;
  }
  if (p->flushed)
    return;
  p->flushed = sox_true;

  while (p->out_generated < p->expected_out) {
    size_t occupancy = fifo_occupancy(&p->input_fifo);
    if (occupancy < p->fft_size) {
      pad = p->fft_size - occupancy;
      float * zeros = fifo_write(&p->input_fifo, pad, NULL);
      for (c = 0; c < pad * p->channels; ++c)
        zeros[c] = 0;
    }
    pvstretch_process_frame(p);
  }
  if (p->out_generated > p->expected_out) {
    size_t excess = (size_t)(p->out_generated - p->expected_out);
    excess = min(excess, fifo_occupancy(&p->output_fifo));
    fifo_trim_by(&p->output_fifo, excess);
  }
}

static float * pvstretch_input(pvstretch_t * p, float const * samples, size_t n)
{
  p->samples_in += n;
  return fifo_write(&p->input_fifo, n, samples);
}

static float const * pvstretch_output(pvstretch_t * p, float * samples, size_t * n)
{
  size_t available = fifo_occupancy(&p->output_fifo);

  if (p->expected_known && p->samples_out < p->expected_out)
    available = min(available, (size_t)(p->expected_out - p->samples_out));
  else if (p->expected_known)
    available = 0;
  *n = min(*n, available);
  p->samples_out += *n;
  return fifo_read(&p->output_fifo, *n, samples);
}

typedef struct {
  pvstretch_t * pvstretch;
  double factor;
  double window_ms;
  double overlap_percent;
  double sensitivity;
  sox_bool drained;
} priv_t;

static int getopts(sox_effect_t * effp, int argc, char ** argv)
{
  priv_t * p = (priv_t *)effp->priv;
  double factor = 1, window_ms = DEFAULT_WINDOW_MS;
  double overlap_percent = DEFAULT_OVERLAP_PERCENT;
  double sensitivity = DEFAULT_TRANSIENT_SENSITIVITY;
  char * end_ptr;

  --argc, ++argv;
  if (!argc)
    return lsx_usage(effp);
  factor = strtod(*argv, &end_ptr);
  if (end_ptr == *argv || *end_ptr || factor < .5 || factor > 2)
    return lsx_usage(effp);
  --argc, ++argv;
  if (argc) {
    window_ms = strtod(*argv, &end_ptr);
    if (end_ptr == *argv || *end_ptr || window_ms < 10 || window_ms > 200)
      return lsx_usage(effp);
    --argc, ++argv;
  }
  if (argc) {
    overlap_percent = strtod(*argv, &end_ptr);
    if (end_ptr == *argv || *end_ptr ||
        overlap_percent < 50 || overlap_percent > 90)
      return lsx_usage(effp);
    --argc, ++argv;
  }
  if (argc) {
    sensitivity = strtod(*argv, &end_ptr);
    if (end_ptr == *argv || *end_ptr || sensitivity < 0 || sensitivity > 1)
      return lsx_usage(effp);
    --argc, ++argv;
  }
  if (argc)
    return lsx_usage(effp);

  p->factor = factor;
  p->window_ms = window_ms;
  p->overlap_percent = overlap_percent;
  p->sensitivity = sensitivity;
  return SOX_SUCCESS;
}

static int start(sox_effect_t * effp)
{
  priv_t * p = (priv_t *)effp->priv;
  size_t channels = effp->in_signal.channels;

  if (p->factor == 1)
    return SOX_EFF_NULL;

  p->pvstretch = pvstretch_create(channels, p->factor,
      effp->in_signal.rate, p->window_ms, p->overlap_percent, p->sensitivity);
  p->pvstretch->expected_out = 0;
  if (effp->in_signal.length != SOX_UNKNOWN_LEN) {
    uint64_t in_length = effp->in_signal.length / channels;
    p->pvstretch->expected_out =
        (uint64_t)floor(in_length / p->factor + 0.5);
    p->pvstretch->expected_known = sox_true;
    effp->out_signal.length = p->pvstretch->expected_out * channels;
  } else {
    effp->out_signal.length = SOX_UNKNOWN_LEN;
  }
  lsx_report("factor=%g window=%gms (%" PRIuPTR " samples) "
      "overlap=%.1f%% (analysis hop %" PRIuPTR ") transient sensitivity=%g",
      p->factor, 1000. * p->pvstretch->fft_size / effp->in_signal.rate,
      p->pvstretch->fft_size,
      100. * (1 - (double)p->pvstretch->analysis_hop / p->pvstretch->fft_size),
      p->pvstretch->analysis_hop, p->sensitivity);
  return SOX_SUCCESS;
}

static int flow(sox_effect_t * effp, const sox_sample_t * ibuf,
    sox_sample_t * obuf, size_t * isamp, size_t * osamp)
{
  priv_t * p = (priv_t *)effp->priv;
  pvstretch_t * v = p->pvstretch;
  size_t channels = effp->in_signal.channels;
  size_t odone = *osamp / channels;
  size_t orequested = odone, i;
  float const * s;
  SOX_SAMPLE_LOCALS;

  s = pvstretch_output(v, NULL, &odone);
  for (i = 0; i < odone * channels; ++i)
    *obuf++ = SOX_FLOAT_32BIT_TO_SAMPLE(*s++, effp->clips);

  if (*isamp && odone < orequested) {
    size_t in_wide = *isamp / channels;
    float * t = pvstretch_input(v, NULL, in_wide);
    for (i = *isamp; i; --i)
      *t++ = SOX_SAMPLE_TO_FLOAT_32BIT(*ibuf++, effp->clips);
    pvstretch_process_available(v);

    *osamp = orequested * channels;
    {
      size_t more = *osamp / channels - odone;
      s = pvstretch_output(v, NULL, &more);
      for (i = 0; i < more * channels; ++i)
        obuf[i] = SOX_FLOAT_32BIT_TO_SAMPLE(*s++, effp->clips);
      odone += more;
    }
    *osamp = odone * channels;
  } else {
    *isamp = 0;
    *osamp = odone * channels;
  }
  return SOX_SUCCESS;
}

static int drain(sox_effect_t * effp, sox_sample_t * obuf, size_t * osamp)
{
  priv_t * p = (priv_t *)effp->priv;
  pvstretch_t * v = p->pvstretch;
  size_t channels = effp->in_signal.channels;
  size_t odone = *osamp / channels, i;
  float const * s;
  SOX_SAMPLE_LOCALS;

  if (!p->drained) {
    pvstretch_flush(v);
    p->drained = sox_true;
  }
  s = pvstretch_output(v, NULL, &odone);
  for (i = 0; i < odone * channels; ++i)
    *obuf++ = SOX_FLOAT_32BIT_TO_SAMPLE(*s++, effp->clips);
  *osamp = odone * channels;
  return !odone && v->expected_known && v->samples_out >= v->expected_out ?
      SOX_EOF : SOX_SUCCESS;
}

static int stop(sox_effect_t * effp)
{
  priv_t * p = (priv_t *)effp->priv;
  pvstretch_delete(p->pvstretch);
  p->pvstretch = NULL;
  return SOX_SUCCESS;
}

sox_effect_handler_t const * lsx_pvstretch_effect_fn(void)
{
  static sox_effect_handler_t handler = {
    "pvstretch",
    "factor [window-ms(46) [overlap-percent(75) [transient-sensitivity(0.5)]]]",
    SOX_EFF_MCHAN | SOX_EFF_LENGTH,
    getopts, start, flow, drain, stop, NULL, sizeof(priv_t)
  };
  return &handler;
}
