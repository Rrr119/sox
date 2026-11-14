/* Focused tests for the pvstretch phase-vocoder effect.
 *
 * The input is a deterministic drum pattern with fast attacks and a tonal
 * body.  Tests verify that stretched attacks land at the expected time and
 * that the tonal body keeps its frequency.
 */

#include "sox.h"

#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TEST_CHUNK 4096
#define HIT_COUNT 3
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

typedef struct {
  unsigned rate;
  unsigned channels;
  double factor;
} test_case_t;

static double hit_time(int hit)
{
  static const double times[HIT_COUNT] = {.35, .80, 1.25};
  return times[hit];
}

static unsigned noise_state = 0x12345678u;

static double noise_sample(void)
{
  noise_state ^= noise_state << 13;
  noise_state ^= noise_state >> 17;
  noise_state ^= noise_state << 5;
  return (double)(noise_state & 0xffff) / 32767.5 - 1;
}

static double drum_sample(double t, unsigned channel)
{
  double sample = .0005 * noise_sample() +
      .0008 * sin(2 * M_PI * (180 + 7 * channel) * t);
  int hit;

  for (hit = 0; hit < HIT_COUNT; ++hit) {
    double dt = t - hit_time(hit);
    if (dt >= 0 && dt < .24) {
      double attack = dt < .0015 ? dt / .0015 : 1;
      double envelope = attack * exp(-dt / .055);
      double body = .46 * sin(2 * M_PI * 440 * dt + .2 * channel) +
          .16 * sin(2 * M_PI * 1320 * dt + .4 * channel);
      double click = .18 * noise_sample() * exp(-dt / .006);
      sample += envelope * body + click;
    }
  }
  return sample;
}

static sox_sample_t * make_drums(unsigned rate, unsigned channels,
    size_t * total, double * duration)
{
  size_t i, c, n;
  sox_sample_t * samples;
  size_t clips = 0;
  SOX_SAMPLE_LOCALS;

  *duration = 1.8;
  n = (size_t)(rate * *duration + .5);
  *total = n * channels;
  samples = malloc(*total * sizeof(*samples));
  if (!samples) {
    fprintf(stderr, "out of memory\n");
    exit(1);
  }
  for (i = 0; i < n; ++i) {
    double t = (double)i / rate;
    for (c = 0; c < channels; ++c)
      samples[i * channels + c] =
          SOX_FLOAT_64BIT_TO_SAMPLE(drum_sample(t, c), clips);
  }
  if (clips)
    fprintf(stderr, "warning: input clipped %" PRIu64 " times\n",
        (uint64_t)clips);
  return samples;
}

static float * run_effect(const test_case_t * test, const sox_sample_t * input,
    size_t input_samples, size_t * output_samples)
{
  sox_effect_t * effect;
  sox_effect_handler_t const * handler = sox_find_effect("pvstretch");
  char factor[32], window[32], overlap[32], sensitivity[32];
  char * options[4];
  char * output = NULL;
  size_t output_cap = 0, output_used = 0, input_pos = 0;
  size_t expected = (size_t)floor(input_samples /
      test->channels / test->factor + .5) * test->channels;
  float * result;
  size_t i;
  int status;
  unsigned guard = 0;

  if (!handler) {
    fprintf(stderr, "pvstretch is not registered\n");
    exit(1);
  }
  effect = sox_create_effect(handler);
  snprintf(factor, sizeof(factor), "%g", test->factor);
  snprintf(window, sizeof(window), "%u", 32u);
  snprintf(overlap, sizeof(overlap), "%u", 75u);
  snprintf(sensitivity, sizeof(sensitivity), "%.1f", .6);
  options[0] = factor;
  options[1] = window;
  options[2] = overlap;
  options[3] = sensitivity;
  if (sox_effect_options(effect, 4, options) != SOX_SUCCESS) {
    fprintf(stderr, "could not set pvstretch options\n");
    exit(1);
  }

  memset(&effect->in_signal, 0, sizeof(effect->in_signal));
  effect->in_signal.rate = test->rate;
  effect->in_signal.channels = test->channels;
  effect->in_signal.precision = 32;
  effect->in_signal.length = input_samples;
  effect->out_signal = effect->in_signal;
  effect->flows = 1;
  effect->flow = 0;
  if (effect->handler.start(effect) != SOX_SUCCESS) {
    fprintf(stderr, "pvstretch failed to start\n");
    exit(1);
  }

  output_cap = expected + TEST_CHUNK;
  output = malloc(output_cap * sizeof(sox_sample_t));
  if (!output) {
    fprintf(stderr, "out of memory\n");
    exit(1);
  }

  while (input_pos < input_samples) {
    size_t isamp = MIN(TEST_CHUNK, input_samples - input_pos);
    size_t osamp = output_cap - output_used;
    if (!osamp)
      break;
    status = effect->handler.flow(effect, input + input_pos,
        (sox_sample_t *)output + output_used, &isamp, &osamp);
    if (status != SOX_SUCCESS) {
      fprintf(stderr, "pvstretch flow failed\n");
      exit(1);
    }
    input_pos += isamp;
    output_used += osamp;
    if (!isamp && !osamp && ++guard > 16) {
      fprintf(stderr, "pvstretch flow stalled\n");
      exit(1);
    }
    if (isamp || osamp)
      guard = 0;
  }

  for (;;) {
    size_t osamp = output_cap - output_used;
    if (!osamp)
      break;
    status = effect->handler.drain(effect,
        (sox_sample_t *)output + output_used, &osamp);
    output_used += osamp;
    if (status == SOX_EOF)
      break;
  }
  effect->handler.stop(effect);
  free(effect->priv);
  free(effect);

  *output_samples = output_used;
  result = malloc(output_used * sizeof(*result));
  if (!result) {
    fprintf(stderr, "out of memory\n");
    exit(1);
  }
  for (i = 0; i < output_used; ++i)
    result[i] = (float)SOX_SAMPLE_TO_FLOAT_64BIT(
        ((sox_sample_t *)output)[i], i);
  free(output);
  return result;
}

static double block_energy(const float * samples, size_t samples_total,
    unsigned channels, unsigned channel, size_t begin, size_t end)
{
  size_t i;
  double energy = 0;

  if (begin >= samples_total / channels)
    return 0;
  end = MIN(end, samples_total / channels);
  for (i = begin; i < end; ++i) {
    double value = samples[i * channels + channel];
    energy += value * value;
  }
  return energy / MAX((size_t)1, end - begin);
}

static double find_onset(const float * samples, size_t samples_total,
    unsigned channels, unsigned channel, unsigned rate, double around,
    double search)
{
  size_t frames = samples_total / channels;
  size_t win = MAX((size_t)1, (size_t)(rate * .002 + .5));
  size_t step = MAX((size_t)1, win / 2);
  size_t first = around > search ?
      (size_t)((around - search) * rate) : 0;
  size_t last = MIN(frames, (size_t)((around + search) * rate));
  size_t best = first;
  double best_rise = -1;
  size_t t;

  for (t = first; t < last; t += step) {
    double before = block_energy(samples, samples_total, channels, channel,
        t > win ? t - win : 0, t);
    double after = block_energy(samples, samples_total, channels, channel,
        t, t + win);
    double rise = after - before;
    if (rise > best_rise) {
      best_rise = rise;
      best = t;
    }
  }
  return (double)best / rate;
}

static double dominant_frequency(const float * samples, size_t samples_total,
    unsigned channels, unsigned channel, unsigned rate,
    double start, double fmin, double fmax)
{
  size_t frames = samples_total / channels;
  size_t begin = (size_t)(start * rate + .5);
  size_t n = MIN((size_t)(rate * .04 + .5), frames - MIN(begin, frames));
  size_t i;
  double best_f = fmin, best_power = -1;
  double f;

  n = MIN(n, (size_t)2048);
  n = MAX(n, (size_t)256);
  if (begin + n > frames)
    n = frames - begin;

  for (f = fmin; f <= fmax; f += 1) {
    double re = 0, im = 0;
    for (i = 0; i < n; ++i) {
      double w = .5 - .5 * cos(2 * M_PI * i / n);
      double x = samples[(begin + i) * channels + channel] * w;
      double phase = 2 * M_PI * f * i / rate;
      re += x * cos(phase);
      im -= x * sin(phase);
    }
    if (re * re + im * im > best_power) {
      best_power = re * re + im * im;
      best_f = f;
    }
  }
  return best_f;
}

static int run_case(const test_case_t * test)
{
  size_t input_samples, output_samples, i;
  double duration;
  sox_sample_t * input = make_drums(test->rate, test->channels,
      &input_samples, &duration);
  float * output;
  double input_onsets[HIT_COUNT], output_onsets[HIT_COUNT];
  double tolerance = .020;

  (void)duration;
  output_samples = 0;
  output = run_effect(test, input, input_samples, &output_samples);

  {
    float * input_float = malloc(input_samples * sizeof(*input_float));
    if (!input_float) {
      fprintf(stderr, "out of memory\n");
      exit(1);
    }
    for (i = 0; i < input_samples; ++i)
      input_float[i] = (float)SOX_SAMPLE_TO_FLOAT_64BIT(input[i], i);

    for (i = 0; i < HIT_COUNT; ++i) {
      double nominal = hit_time((int)i);
      input_onsets[i] = find_onset(input_float, input_samples,
          test->channels, 0, test->rate, nominal, .06);
      output_onsets[i] = find_onset(output, output_samples,
          test->channels, 0, test->rate, input_onsets[i] / test->factor,
          MAX(.15, .12 / test->factor));
      if (fabs(output_onsets[i] - input_onsets[i] / test->factor) > tolerance) {
        fprintf(stderr,
            "FAIL rate=%u ch=%u factor=%g hit=%" PRIuPTR
            " onset %.4fs expected %.4fs\n",
            test->rate, test->channels, test->factor, i,
            output_onsets[i], input_onsets[i] / test->factor);
        return 1;
      }
      if (test->channels > 1) {
        double right = find_onset(output, output_samples, test->channels, 1,
            test->rate, input_onsets[i] / test->factor,
            MAX(.15, .12 / test->factor));
        if (fabs(right - output_onsets[i]) > .010) {
          fprintf(stderr, "FAIL rate=%u factor=%g channel timing mismatch\n",
              test->rate, test->factor);
          return 1;
        }
      }
    }

    for (i = 0; i < HIT_COUNT; ++i) {
      double input_f = dominant_frequency(input_float, input_samples,
          test->channels, 0, test->rate, input_onsets[i] + .01, 360, 520);
      double output_f = dominant_frequency(output, output_samples,
          test->channels, 0, test->rate, output_onsets[i] + .04, 360, 520);
      if (fabs(input_f - 440) > 2 || fabs(output_f - 440) > 3 ||
          fabs(output_f - input_f) > 3) {
        double later_f = dominant_frequency(output, output_samples,
            test->channels, 0, test->rate, output_onsets[i] + .04, 360, 520);
        fprintf(stderr,
            "FAIL rate=%u ch=%u factor=%g hit=%" PRIuPTR
            " spectrum input %.0fHz output %.0fHz (later %.0fHz)\n",
            test->rate, test->channels, test->factor, i,
            input_f, output_f, later_f);
        return 1;
      }
    }
    free(input_float);
  }

  if ((size_t)floor(input_samples / test->channels / test->factor + .5) *
      test->channels != output_samples) {
    fprintf(stderr, "FAIL rate=%u factor=%g length %" PRIuPTR " expected %" PRIuPTR "\n",
        test->rate, test->factor, output_samples,
        (size_t)floor(input_samples / test->channels / test->factor + .5) *
        test->channels);
    free(input);
    free(output);
    return 1;
  }

  printf("ok     rate=%u channels=%u factor=%g\n",
      test->rate, test->channels, test->factor);
  free(input);
  free(output);
  return 0;
}

int main(void)
{
  static const test_case_t cases[] = {
    {48000, 2, 1.5},
    {44100, 1, .75},
    {8000, 2, 2.0},
    {96000, 1, .5}
  };
  size_t i;
  int failed = 0;

  if (sox_init() != SOX_SUCCESS) {
    fprintf(stderr, "sox_init failed\n");
    return 1;
  }
  for (i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i)
    failed |= run_case(&cases[i]);
  sox_quit();
  return failed ? 1 : 0;
}
