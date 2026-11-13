/* sox-quality: EBU R128 / BS.1770 quality report generator.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 */

#include "sox.h"

#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <dirent.h>
#include <sys/stat.h>
#endif

#define SCHEMA_VERSION 1
#define BIN_COUNT 30
#define BLOCK_BINS 4
#define DEFAULT_ABSOLUTE_GATE (-70.0)
#define DEFAULT_RELATIVE_GATE 10.0
#define TRUE_PEAK_FACTOR 4
#define TRUE_PEAK_TAPS 49
#define TRUE_PEAK_DELAY ((TRUE_PEAK_TAPS + TRUE_PEAK_FACTOR - 1) / TRUE_PEAK_FACTOR)
#define READ_FRAMES 4096

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
  int recursive;
  int quiet;
  int batch;
  char const * output;
  double absolute_gate;
  double relative_gate;
} options_t;

typedef struct {
  char ** items;
  size_t len;
  size_t cap;
} string_list_t;

typedef struct {
  double * items;
  size_t len;
  size_t cap;
} double_vector_t;

typedef struct {
  unsigned channels;
  unsigned long rate;
  size_t samples_per_bin;
  size_t block_frames;
  size_t short_term_frames;
  double a[5];
  double b[5];
  double * filter_state;
  double * bin_energy;
  double * current_bin;
  double * raw_bin_energy;
  double raw_current_bin;
  size_t bin_pos;
  size_t bins_closed;
  size_t frames_in_bin;
  uint64_t frames;
  uint64_t scalar_samples;
  uint64_t clipped_samples;
  sox_sample_t full_scale_positive;
  double dc_sum;
  double * channel_weight;
  double * true_peak_ring;
  double * true_peak_coeff;
  size_t true_peak_pos;
  double true_peak_linear;
  double short_term_max_lufs;
  double_vector_t block_energies;
  double_vector_t raw_block_energies;
} analyzer_t;

typedef struct {
  char * path;
  int ok;
  char error[160];
  char format[64];
  char encoding[64];
  unsigned rate;
  unsigned channels;
  unsigned precision;
  unsigned bits_per_sample;
  uint64_t frames;
  uint64_t scalar_samples;
  uint64_t clipped_samples;
  double duration;
  double integrated_loudness;
  double short_term_max_loudness;
  double true_peak;
  double dc_offset;
  double dc_offset_dbfs;
  double noise_floor_dbfs;
} file_result_t;

typedef struct {
  uint64_t files_analyzed;
  uint64_t files_failed;
  uint64_t files_with_clipping;
  uint64_t scalar_samples;
  uint64_t clipped_samples;
  double duration;
  double dc_sum;
  double short_term_max_loudness;
  double true_peak;
  double_vector_t block_energies;
  double_vector_t raw_block_energies;
} summary_t;

static int compare_doubles(void const * a, void const * b);

static void * xcalloc(size_t count, size_t size)
{
  void * p = calloc(count, size);
  if (!p) {
    fprintf(stderr, "sox-quality: out of memory\n");
    exit(1);
  }
  return p;
}

static char * dup_string(char const * s)
{
  size_t len = strlen(s) + 1;
  char * result = xcalloc(len, 1);
  memcpy(result, s, len);
  return result;
}

static void string_list_push(string_list_t * list, char const * s)
{
  if (list->len == list->cap) {
    list->cap = list->cap ? 2 * list->cap : 16;
    list->items = realloc(list->items, list->cap * sizeof(*list->items));
    if (!list->items) {
      fprintf(stderr, "sox-quality: out of memory\n");
      exit(1);
    }
  }
  list->items[list->len++] = dup_string(s);
}

static void string_list_free(string_list_t * list)
{
  size_t i;
  for (i = 0; i < list->len; ++i)
    free(list->items[i]);
  free(list->items);
}

static int compare_strings(void const * a, void const * b)
{
  return strcmp(*(char * const *)a, *(char * const *)b);
}

static void string_list_sort(string_list_t * list)
{
  qsort(list->items, list->len, sizeof(*list->items), compare_strings);
}

static void double_vector_push(double_vector_t * vector, double value)
{
  if (vector->len == vector->cap) {
    vector->cap = vector->cap ? 2 * vector->cap : 256;
    vector->items = realloc(vector->items, vector->cap * sizeof(*vector->items));
    if (!vector->items) {
      fprintf(stderr, "sox-quality: out of memory\n");
      exit(1);
    }
  }
  vector->items[vector->len++] = value;
}

static void double_vector_free(double_vector_t * vector)
{
  free(vector->items);
  vector->items = NULL;
  vector->len = vector->cap = 0;
}

static double energy_to_loudness(double energy)
{
  return energy > 0? -0.691 + 10 * log10(energy) : -HUGE_VAL;
}

static double energy_to_dbfs(double energy)
{
  return energy > 0? 10 * log10(energy) : -HUGE_VAL;
}

static double linear_to_dbfs(double value)
{
  return value > 0? 20 * log10(value) : -HUGE_VAL;
}

static double percentile_energy(double_vector_t const * vector, double fraction)
{
  double * copy;
  size_t index;
  double result;

  if (!vector->len)
    return -HUGE_VAL;
  copy = xcalloc(vector->len, sizeof(*copy));
  memcpy(copy, vector->items, vector->len * sizeof(*copy));
  qsort(copy, vector->len, sizeof(*copy), compare_doubles);
  index = (size_t)(fraction * (vector->len - 1));
  result = copy[index];
  free(copy);
  return result;
}

static double gated_loudness(double_vector_t const * blocks,
    double absolute_gate, double relative_gate)
{
  double absolute_sum = 0;
  double relative_sum = 0;
  double relative_threshold;
  size_t absolute_count = 0;
  size_t relative_count = 0;
  size_t i;

  for (i = 0; i < blocks->len; ++i) {
    double energy = blocks->items[i];
    if (energy_to_loudness(energy) >= absolute_gate) {
      absolute_sum += energy;
      ++absolute_count;
    }
  }
  if (!absolute_count)
    return -HUGE_VAL;

  relative_threshold = energy_to_loudness(absolute_sum / absolute_count) -
      relative_gate;
  for (i = 0; i < blocks->len; ++i) {
    double energy = blocks->items[i];
    if (energy_to_loudness(energy) >= absolute_gate &&
        energy_to_loudness(energy) >= relative_threshold) {
      relative_sum += energy;
      ++relative_count;
    }
  }
  return relative_count? energy_to_loudness(relative_sum / relative_count)
      : -HUGE_VAL;
}

static int compare_doubles(void const * a, void const * b)
{
  double x = *(double const *)a;
  double y = *(double const *)b;
  return x < y? -1 : x > y? 1 : 0;
}

static void init_k_weighting(analyzer_t * a)
{
  double const f0_shelf = 1681.974450955533;
  double const gain_shelf = 3.999843853973347;
  double const q_shelf = 0.7071752369554196;
  double const f0_highpass = 38.13547087602444;
  double const q_highpass = 0.5003270373238773;
  double k = tan(M_PI * f0_shelf / a->rate);
  double vh = pow(10, gain_shelf / 20);
  double vb = pow(vh, 0.4996667741545416);
  double a0 = 1 + k / q_shelf + k * k;
  double pb0 = (vh + vb * k / q_shelf + k * k) / a0;
  double pb1 = 2 * (k * k - vh) / a0;
  double pb2 = (vh - vb * k / q_shelf + k * k) / a0;
  double pa1 = 2 * (k * k - 1) / a0;
  double pa2 = (1 - k / q_shelf + k * k) / a0;
  double ra1, ra2;

  k = tan(M_PI * f0_highpass / a->rate);
  a0 = 1 + k / q_highpass + k * k;
  ra1 = 2 * (k * k - 1) / a0;
  ra2 = (1 - k / q_highpass + k * k) / a0;

  a->b[0] = pb0;
  a->b[1] = -2 * pb0 + pb1;
  a->b[2] = pb0 - 2 * pb1 + pb2;
  a->b[3] = pb1 - 2 * pb2;
  a->b[4] = pb2;
  a->a[0] = 1;
  a->a[1] = pa1 + ra1;
  a->a[2] = ra2 + pa1 * ra1 + pa2;
  a->a[3] = pa1 * ra2 + pa2 * ra1;
  a->a[4] = pa2 * ra2;
}

static void init_true_peak(analyzer_t * a)
{
  unsigned i;
  for (i = 0; i < TRUE_PEAK_TAPS; ++i) {
    double m = (double)i - (TRUE_PEAK_TAPS - 1) / 2.0;
    double c = fabs(m) < 1e-10? 1 :
        sin(m * M_PI / TRUE_PEAK_FACTOR) /
        (m * M_PI / TRUE_PEAK_FACTOR);
    a->true_peak_coeff[i] =
        c * 0.5 * (1 - cos(2 * M_PI * i / (TRUE_PEAK_TAPS - 1)));
  }
}

static void init_channel_weights(analyzer_t * a)
{
  unsigned i;
  for (i = 0; i < a->channels; ++i) {
    if (a->channels == 4)
      a->channel_weight[i] = i < 2? 1 : 1.41;
    else if (a->channels == 5)
      a->channel_weight[i] = i < 3? 1 : 1.41;
    else if (i < 3)
      a->channel_weight[i] = 1;
    else if (i == 3)
      a->channel_weight[i] = 0;
    else if (i < 6)
      a->channel_weight[i] = 1.41;
    else
      a->channel_weight[i] = 0;
  }
}

static sox_sample_t positive_full_scale(sox_encoding_t encoding,
    unsigned precision)
{
  if (encoding == SOX_ENCODING_FLOAT || encoding == SOX_ENCODING_FLOAT_TEXT ||
      encoding == SOX_ENCODING_WAVPACKF)
    return SOX_SAMPLE_MAX;
  if (precision == 0 || precision >= 32)
    return SOX_SAMPLE_MAX;
  return (sox_sample_t)(SOX_INT_MAX(precision) << (32 - precision));
}

static int analyzer_init(analyzer_t * a, unsigned long rate, unsigned channels,
    unsigned precision, sox_encoding_t encoding)
{
  memset(a, 0, sizeof(*a));
  if (!rate || !channels) {
    fprintf(stderr, "sox-quality: invalid sample rate or channel count\n");
    return -1;
  }
  a->rate = rate;
  a->channels = channels;
  a->samples_per_bin = (rate + 5) / 10;
  if (!a->samples_per_bin)
    a->samples_per_bin = 1;
  a->block_frames = BLOCK_BINS * a->samples_per_bin;
  a->short_term_frames = BIN_COUNT * a->samples_per_bin;
  a->full_scale_positive = positive_full_scale(encoding, precision);
  a->short_term_max_lufs = -HUGE_VAL;

  a->filter_state = xcalloc(channels * 4, sizeof(*a->filter_state));
  a->bin_energy = xcalloc(BIN_COUNT * channels, sizeof(*a->bin_energy));
  a->current_bin = xcalloc(channels, sizeof(*a->current_bin));
  a->raw_bin_energy = xcalloc(BIN_COUNT, sizeof(*a->raw_bin_energy));
  a->channel_weight = xcalloc(channels, sizeof(*a->channel_weight));
  a->true_peak_ring = xcalloc(channels * TRUE_PEAK_DELAY,
      sizeof(*a->true_peak_ring));
  a->true_peak_coeff = xcalloc(TRUE_PEAK_TAPS, sizeof(*a->true_peak_coeff));

  init_k_weighting(a);
  init_true_peak(a);
  init_channel_weights(a);
  return 0;
}

static void analyzer_free(analyzer_t * a)
{
  free(a->filter_state);
  free(a->bin_energy);
  free(a->current_bin);
  free(a->raw_bin_energy);
  free(a->channel_weight);
  free(a->true_peak_ring);
  free(a->true_peak_coeff);
  double_vector_free(&a->block_energies);
  double_vector_free(&a->raw_block_energies);
}

static double k_weight(analyzer_t * a, unsigned channel, double sample)
{
  double * state = a->filter_state + channel * 4;
  double v0 = sample - a->a[1] * state[0] - a->a[2] * state[1] -
      a->a[3] * state[2] - a->a[4] * state[3];
  double result = a->b[0] * v0 + a->b[1] * state[0] + a->b[2] * state[1] +
      a->b[3] * state[2] + a->b[4] * state[3];
  state[3] = state[2];
  state[2] = state[1];
  state[1] = state[0];
  state[0] = v0;
  return result;
}

static void update_true_peak(analyzer_t * a)
{
  unsigned channel;
  unsigned phase;
  unsigned tap;

  for (channel = 0; channel < a->channels; ++channel) {
    double * ring = a->true_peak_ring + channel * TRUE_PEAK_DELAY;
    for (phase = 0; phase < TRUE_PEAK_FACTOR; ++phase) {
      double value = 0;
      for (tap = phase; tap < TRUE_PEAK_TAPS; tap += TRUE_PEAK_FACTOR) {
        int index = (int)a->true_peak_pos - (int)(tap / TRUE_PEAK_FACTOR);
        if (index < 0)
          index += TRUE_PEAK_DELAY;
        value += ring[index] * a->true_peak_coeff[tap];
      }
      value = fabs(value);
      if (value > a->true_peak_linear)
        a->true_peak_linear = value;
    }
  }
}

static void close_bin(analyzer_t * a)
{
  unsigned channel;
  size_t i;

  for (channel = 0; channel < a->channels; ++channel) {
    a->bin_energy[a->bin_pos * a->channels + channel] =
        a->current_bin[channel];
    a->current_bin[channel] = 0;
  }
  a->raw_bin_energy[a->bin_pos] = a->raw_current_bin;
  a->raw_current_bin = 0;
  ++a->bins_closed;

  if (a->bins_closed >= BLOCK_BINS) {
    double energy = 0;
    double raw_energy = 0;
    for (i = 0; i < BLOCK_BINS; ++i) {
      size_t bin = (a->bin_pos + BIN_COUNT - i) % BIN_COUNT;
      for (channel = 0; channel < a->channels; ++channel)
        energy += a->channel_weight[channel] *
            a->bin_energy[bin * a->channels + channel];
      raw_energy += a->raw_bin_energy[bin];
    }
    double_vector_push(&a->block_energies, energy / a->block_frames);
    double_vector_push(&a->raw_block_energies,
        raw_energy / a->block_frames);
  }

  if (a->bins_closed >= BIN_COUNT) {
    double energy = 0;
    double loudness;
    for (i = 0; i < BIN_COUNT; ++i) {
      size_t bin = (a->bin_pos + BIN_COUNT - i) % BIN_COUNT;
      for (channel = 0; channel < a->channels; ++channel)
        energy += a->channel_weight[channel] *
            a->bin_energy[bin * a->channels + channel];
    }
    loudness = energy_to_loudness(energy / a->short_term_frames);
    if (loudness > a->short_term_max_lufs)
      a->short_term_max_lufs = loudness;
  }

  a->bin_pos = (a->bin_pos + 1) % BIN_COUNT;
}

static void analyzer_process(analyzer_t * a, sox_sample_t const * samples,
    size_t frames)
{
  size_t frame;
  unsigned channel;

  for (frame = 0; frame < frames; ++frame) {
    for (channel = 0; channel < a->channels; ++channel) {
      sox_sample_t native = samples[frame * a->channels + channel];
      double sample = SOX_SAMPLE_TO_FLOAT_64BIT(native,);
      double filtered = k_weight(a, channel, sample);
      double energy = filtered * filtered;

      if (native == SOX_SAMPLE_MIN || native == a->full_scale_positive)
        ++a->clipped_samples;
      if (fabs(sample) > a->true_peak_linear)
        a->true_peak_linear = fabs(sample);
      a->dc_sum += sample;
      ++a->scalar_samples;
      a->current_bin[channel] += energy;
      a->raw_current_bin += sample * sample;
      a->true_peak_ring[channel * TRUE_PEAK_DELAY + a->true_peak_pos] = sample;
    }
    update_true_peak(a);
    a->true_peak_pos = (a->true_peak_pos + 1) % TRUE_PEAK_DELAY;
    ++a->frames;
    if (++a->frames_in_bin == a->samples_per_bin) {
      a->frames_in_bin = 0;
      close_bin(a);
    }
  }
}

static void summary_init(summary_t * summary)
{
  memset(summary, 0, sizeof(*summary));
  summary->short_term_max_loudness = -HUGE_VAL;
  summary->true_peak = -HUGE_VAL;
}

static void summary_free(summary_t * summary)
{
  double_vector_free(&summary->block_energies);
  double_vector_free(&summary->raw_block_energies);
}

static void summary_add(summary_t * summary, analyzer_t const * a,
    file_result_t const * result)
{
  size_t i;

  ++summary->files_analyzed;
  summary->scalar_samples += result->scalar_samples;
  summary->clipped_samples += result->clipped_samples;
  summary->duration += result->duration;
  summary->dc_sum += a->dc_sum;
  if (result->clipped_samples)
    ++summary->files_with_clipping;
  if (result->short_term_max_loudness > summary->short_term_max_loudness)
    summary->short_term_max_loudness = result->short_term_max_loudness;
  if (result->true_peak > summary->true_peak)
    summary->true_peak = result->true_peak;
  for (i = 0; i < a->block_energies.len; ++i)
    double_vector_push(&summary->block_energies, a->block_energies.items[i]);
  for (i = 0; i < a->raw_block_energies.len; ++i)
    double_vector_push(&summary->raw_block_energies,
        a->raw_block_energies.items[i]);
}

static int analyze_file(char const * path, options_t const * options,
    file_result_t * result, summary_t * summary)
{
  sox_format_t * input = NULL;
  sox_sample_t * buffer = NULL;
  analyzer_t analyzer;
  size_t read;
  size_t frames;
  double noise_energy;
  int initialized = 0;
  int ok = 0;

  memset(result, 0, sizeof(*result));
  result->path = dup_string(path);
  result->integrated_loudness = -HUGE_VAL;
  result->short_term_max_loudness = -HUGE_VAL;
  result->true_peak = -HUGE_VAL;
  result->dc_offset_dbfs = -HUGE_VAL;
  result->noise_floor_dbfs = -HUGE_VAL;

  input = sox_open_read(path, NULL, NULL, NULL);
  if (!input) {
    snprintf(result->error, sizeof(result->error), "cannot open audio file");
    goto done;
  }
  if (!input->signal.rate || !input->signal.channels) {
    snprintf(result->error, sizeof(result->error),
        "audio format has unknown sample rate or channels");
    goto done;
  }

  result->rate = (unsigned)(input->signal.rate + 0.5);
  result->channels = input->signal.channels;
  result->precision = input->signal.precision;
  result->bits_per_sample = input->encoding.bits_per_sample;
  snprintf(result->format, sizeof(result->format), "%s",
      input->filetype? input->filetype : "unknown");
  if (input->encoding.encoding < SOX_ENCODINGS)
    snprintf(result->encoding, sizeof(result->encoding), "%s",
        sox_encodings_info[input->encoding.encoding].name);
  else
    snprintf(result->encoding, sizeof(result->encoding), "unknown");

  if (analyzer_init(&analyzer, result->rate, result->channels,
      result->precision, input->encoding.encoding) != 0) {
    snprintf(result->error, sizeof(result->error),
        "cannot initialize loudness analyzer");
    goto done;
  }
  initialized = 1;
  buffer = xcalloc((size_t)READ_FRAMES * result->channels, sizeof(*buffer));

  while ((read = sox_read(input, buffer,
      (size_t)READ_FRAMES * result->channels)) > 0) {
    frames = read / result->channels;
    if (!frames)
      break;
    analyzer_process(&analyzer, buffer, frames);
  }
  analyzer.clipped_samples += input->clips;

  result->frames = analyzer.frames;
  result->scalar_samples = analyzer.scalar_samples;
  result->clipped_samples = analyzer.clipped_samples;
  result->duration = result->rate? (double)analyzer.frames / result->rate : 0;
  result->integrated_loudness = gated_loudness(&analyzer.block_energies,
      options->absolute_gate, options->relative_gate);
  result->short_term_max_loudness = analyzer.short_term_max_lufs;
  result->true_peak = linear_to_dbfs(analyzer.true_peak_linear);
  result->dc_offset = result->scalar_samples?
      analyzer.dc_sum / result->scalar_samples : 0;
  result->dc_offset_dbfs = linear_to_dbfs(fabs(result->dc_offset));
  noise_energy = percentile_energy(&analyzer.raw_block_energies, .10);
  result->noise_floor_dbfs = energy_to_dbfs(noise_energy);

  if (sox_close(input) != SOX_SUCCESS) {
    input = NULL;
    snprintf(result->error, sizeof(result->error), "error closing audio file");
    goto done;
  }
  input = NULL;
  summary_add(summary, &analyzer, result);
  result->ok = 1;
  ok = 1;

done:
  if (!result->ok && !result->error[0])
    snprintf(result->error, sizeof(result->error), "audio analysis failed");
  free(buffer);
  if (initialized)
    analyzer_free(&analyzer);
  if (input)
    sox_close(input);
  if (!ok)
    ++summary->files_failed;
  return ok? 0 : -1;
}

static char * path_join(char const * directory, char const * name)
{
  size_t dlen = strlen(directory);
  size_t nlen = strlen(name);
  char * result = xcalloc(dlen + nlen + 2, 1);
  memcpy(result, directory, dlen);
  if (dlen && directory[dlen - 1] != '/' && directory[dlen - 1] != '\\')
    result[dlen++] = '/';
  memcpy(result + dlen, name, nlen + 1);
  return result;
}

static int path_is_directory(char const * path)
{
#ifdef _WIN32
  DWORD attributes = GetFileAttributesA(path);
  return attributes != INVALID_FILE_ATTRIBUTES &&
      (attributes & FILE_ATTRIBUTE_DIRECTORY);
#else
  struct stat info;
  return stat(path, &info) == 0 && S_ISDIR(info.st_mode);
#endif
}

#ifdef _WIN32
static void scan_directory(char const * path, options_t const * options,
    string_list_t * files)
{
  char * pattern = path_join(path, "*");
  WIN32_FIND_DATAA entry;
  HANDLE handle = FindFirstFileA(pattern, &entry);
  free(pattern);
  if (handle == INVALID_HANDLE_VALUE)
    return;
  do {
    char * child;
    if (!strcmp(entry.cFileName, ".") || !strcmp(entry.cFileName, ".."))
      continue;
    child = path_join(path, entry.cFileName);
    if (entry.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
      if (options->recursive)
        scan_directory(child, options, files);
    }
    else
      string_list_push(files, child);
    free(child);
  } while (FindNextFileA(handle, &entry));
  FindClose(handle);
}
#else
static void scan_directory(char const * path, options_t const * options,
    string_list_t * files)
{
  DIR * directory = opendir(path);
  struct dirent * entry;
  if (!directory) {
    fprintf(stderr, "sox-quality: cannot read directory `%s': %s\n",
        path, strerror(errno));
    return;
  }
  while ((entry = readdir(directory)) != NULL) {
    char * child;
    struct stat info;
    if (!strcmp(entry->d_name, ".") || !strcmp(entry->d_name, ".."))
      continue;
    child = path_join(path, entry->d_name);
    if (lstat(child, &info) == 0) {
      if (S_ISDIR(info.st_mode)) {
        if (options->recursive)
          scan_directory(child, options, files);
      }
      else if (S_ISREG(info.st_mode))
        string_list_push(files, child);
    }
    free(child);
  }
  closedir(directory);
}
#endif

static void collect_path(char const * path, options_t * options,
    string_list_t * files)
{
  if (path_is_directory(path)) {
    options->batch = 1;
    scan_directory(path, options, files);
  }
  else
    string_list_push(files, path);
}

static void print_json_string(FILE * out, char const * text)
{
  unsigned char const * p = (unsigned char const *)text;
  putc('"', out);
  for (; *p; ++p) {
    switch (*p) {
      case '"': fputs("\\\"", out); break;
      case '\\': fputs("\\\\", out); break;
      case '\b': fputs("\\b", out); break;
      case '\f': fputs("\\f", out); break;
      case '\n': fputs("\\n", out); break;
      case '\r': fputs("\\r", out); break;
      case '\t': fputs("\\t", out); break;
      default:
        if (*p < 0x20)
          fprintf(out, "\\u%04x", *p);
        else
          putc(*p, out);
    }
  }
  putc('"', out);
}

static void print_number(FILE * out, double value, int precision)
{
  if (isfinite(value)) {
    if (value == 0)
      value = 0;
    fprintf(out, "%.*f", precision, value);
  }
  else
    fputs("null", out);
}

static void print_metrics(FILE * out, file_result_t const * result,
    char const * indent, int trailing_comma)
{
  fprintf(out, "{\n");
  fprintf(out, "%s  \"integrated_loudness_lufs\": ", indent);
  print_number(out, result->integrated_loudness, 6);
  fprintf(out, ",\n%s  \"short_term_loudness_max_lufs\": ", indent);
  print_number(out, result->short_term_max_loudness, 6);
  fprintf(out, ",\n%s  \"true_peak_dbtp\": ", indent);
  print_number(out, result->true_peak, 6);
  fprintf(out, ",\n%s  \"clipped_samples\": %" PRIu64 ",\n",
      indent, result->clipped_samples);
  fprintf(out, "%s  \"dc_offset\": ", indent);
  print_number(out, result->dc_offset, 9);
  fprintf(out, ",\n%s  \"dc_offset_dbfs\": ", indent);
  print_number(out, result->dc_offset_dbfs, 6);
  fprintf(out, ",\n%s  \"noise_floor_dbfs\": ", indent);
  print_number(out, result->noise_floor_dbfs, 6);
  fprintf(out, "\n%s}%s\n", indent, trailing_comma? "," : "");
}

static void print_summary(FILE * out, summary_t const * summary,
    size_t files_total, options_t const * options)
{
  double integrated = gated_loudness(&summary->block_energies,
      options->absolute_gate, options->relative_gate);
  double noise = energy_to_dbfs(
      percentile_energy(&summary->raw_block_energies, .10));
  double dc = summary->scalar_samples?
      summary->dc_sum / summary->scalar_samples : 0;

  fprintf(out, "  \"summary\": {\n");
  fprintf(out, "    \"files_total\": %zu,\n", files_total);
  fprintf(out, "    \"files_analyzed\": %" PRIu64 ",\n",
      summary->files_analyzed);
  fprintf(out, "    \"files_failed\": %" PRIu64 ",\n", summary->files_failed);
  fprintf(out, "    \"files_with_clipping\": %" PRIu64 ",\n",
      summary->files_with_clipping);
  fprintf(out, "    \"duration_seconds\": ");
  print_number(out, summary->duration, 6);
  fprintf(out, ",\n    \"scalar_samples\": %" PRIu64 ",\n",
      summary->scalar_samples);
  fprintf(out, "    \"integrated_loudness_lufs\": ");
  print_number(out, integrated, 6);
  fprintf(out, ",\n    \"short_term_loudness_max_lufs\": ");
  print_number(out, summary->short_term_max_loudness, 6);
  fprintf(out, ",\n    \"true_peak_dbtp\": ");
  print_number(out, summary->true_peak, 6);
  fprintf(out, ",\n    \"clipped_samples\": %" PRIu64 ",\n",
      summary->clipped_samples);
  fprintf(out, "    \"dc_offset\": ");
  print_number(out, dc, 9);
  fprintf(out, ",\n    \"dc_offset_dbfs\": ");
  print_number(out, linear_to_dbfs(fabs(dc)), 6);
  fprintf(out, ",\n    \"noise_floor_dbfs\": ");
  print_number(out, noise, 6);
  fprintf(out, "\n  },\n");
}

static void print_method(FILE * out, options_t const * options)
{
  fprintf(out, "  \"method\": {\n");
  fprintf(out, "    \"standard\": \"ITU-R BS.1770-4 / EBU R128\",\n");
  fprintf(out, "    \"integrated_loudness\": {\n");
  fprintf(out, "      \"block_ms\": 400,\n");
  fprintf(out, "      \"overlap_percent\": 75,\n");
  fprintf(out, "      \"hop_ms\": 100,\n");
  fprintf(out, "      \"absolute_gate_lufs\": ");
  print_number(out, options->absolute_gate, 3);
  fprintf(out, ",\n      \"relative_gate_lu\": ");
  print_number(out, options->relative_gate, 3);
  fprintf(out, "\n    },\n");
  fprintf(out, "    \"short_term_loudness\": {\n");
  fprintf(out, "      \"window_ms\": 3000,\n");
  fprintf(out, "      \"hop_ms\": 100\n");
  fprintf(out, "    },\n");
  fprintf(out, "    \"true_peak\": {\n");
  fprintf(out, "      \"oversampling_factor\": 4,\n");
  fprintf(out, "      \"interpolation\": \"49-tap Hanning-windowed sinc\"\n");
  fprintf(out, "    },\n");
  fprintf(out, "    \"clipping\": \"sample at source positive or negative full-scale code\",\n");
  fprintf(out, "    \"dc_offset\": \"mean of all decoded samples, normalized to [-1, 1)\",\n");
  fprintf(out, "    \"noise_floor\": \"10th percentile of complete 400 ms unweighted RMS blocks, dBFS\"\n");
  fprintf(out, "  },\n");
}

static void print_report(FILE * out, options_t const * options,
    string_list_t const * files, file_result_t const * results,
    summary_t const * summary)
{
  time_t now = time(NULL);
  struct tm * utc = gmtime(&now);
  char timestamp[32] = "unknown";
  size_t i;

  if (utc)
    strftime(timestamp, sizeof(timestamp), "%Y-%m-%dT%H:%M:%SZ", utc);
  fprintf(out, "{\n");
  fprintf(out, "  \"schema_version\": %d,\n", SCHEMA_VERSION);
  fprintf(out, "  \"report\": \"sox-quality\",\n");
  fprintf(out, "  \"tool\": {\"name\": \"sox-quality\", \"version\": ");
  print_json_string(out, sox_version());
  fprintf(out, "},\n");
  fprintf(out, "  \"generated_at\": ");
  print_json_string(out, timestamp);
  fprintf(out, ",\n");
  fprintf(out, "  \"mode\": \"%s\",\n", options->batch? "batch" : "single");
  print_method(out, options);
  print_summary(out, summary, files->len, options);
  fprintf(out, "  \"files\": [\n");
  for (i = 0; i < files->len; ++i) {
    file_result_t const * result = &results[i];
    fprintf(out, "    {\n      \"path\": ");
    print_json_string(out, result->path);
    fprintf(out, ",\n      \"status\": \"%s\"", result->ok? "ok" : "error");
    if (!result->ok) {
      fprintf(out, ",\n      \"error\": ");
      print_json_string(out, result->error);
      fprintf(out, "\n    }%s\n", i + 1 < files->len? "," : "");
      continue;
    }
    fprintf(out, ",\n      \"format\": {\n");
    fprintf(out, "        \"type\": ");
    print_json_string(out, result->format);
    fprintf(out, ",\n        \"encoding\": ");
    print_json_string(out, result->encoding);
    fprintf(out, ",\n        \"sample_rate_hz\": %u,\n", result->rate);
    fprintf(out, "        \"channels\": %u,\n", result->channels);
    fprintf(out, "        \"bit_depth\": %u,\n", result->bits_per_sample);
    fprintf(out, "        \"precision\": %u,\n", result->precision);
    fprintf(out, "        \"frames\": %" PRIu64 ",\n", result->frames);
    fprintf(out, "        \"duration_seconds\": ");
    print_number(out, result->duration, 6);
    fprintf(out, "\n      },\n");
    fprintf(out, "      \"metrics\": ");
    print_metrics(out, result, "      ", 0);
    fprintf(out, "    }%s\n", i + 1 < files->len? "," : "");
  }
  fprintf(out, "  ]\n}\n");
}

static void usage(FILE * out)
{
  fprintf(out,
      "Usage: sox-quality [options] FILE|DIRECTORY...\n"
      "\n"
      "Generate a JSON quality report with EBU R128 loudness, true peak,\n"
      "clipping, DC offset and an estimated noise floor.\n"
      "\n"
      "Options:\n"
      "  -o, --output FILE       write JSON to FILE instead of stdout\n"
      "  -r, --recursive         recurse into directories (default)\n"
      "      --no-recursive      analyze only files directly in directories\n"
      "      --absolute-gate DB  absolute loudness gate (default -70 LUFS)\n"
      "      --relative-gate LU  relative loudness gate (default 10 LU)\n"
      "  -q, --quiet             suppress libSoX messages\n"
      "  -h, --help              show this help\n"
      "  -V, --version           show version\n");
}

static int parse_double(char const * text, double * value)
{
  char * end;
  errno = 0;
  *value = strtod(text, &end);
  return errno || end == text || *end;
}

static int parse_options(int argc, char ** argv, options_t * options,
    string_list_t * files)
{
  int i;
  int end_options = 0;
  memset(options, 0, sizeof(*options));
  options->recursive = 1;
  options->absolute_gate = DEFAULT_ABSOLUTE_GATE;
  options->relative_gate = DEFAULT_RELATIVE_GATE;

  for (i = 1; i < argc; ++i) {
    char const * arg = argv[i];
    if (!end_options && !strcmp(arg, "--")) {
      end_options = 1;
      continue;
    }
    if (!end_options && (!strcmp(arg, "-h") || !strcmp(arg, "--help"))) {
      usage(stdout);
      exit(0);
    }
    if (!end_options && (!strcmp(arg, "-V") || !strcmp(arg, "--version"))) {
      printf("sox-quality %s\n", sox_version());
      exit(0);
    }
    if (!end_options && (!strcmp(arg, "-r") || !strcmp(arg, "--recursive"))) {
      options->recursive = 1;
      continue;
    }
    if (!end_options && !strcmp(arg, "--no-recursive")) {
      options->recursive = 0;
      continue;
    }
    if (!end_options && (!strcmp(arg, "-q") || !strcmp(arg, "--quiet"))) {
      options->quiet = 1;
      continue;
    }
    if (!end_options && (!strcmp(arg, "-o") || !strcmp(arg, "--output"))) {
      if (++i == argc) {
        fprintf(stderr, "sox-quality: option `%s' requires an argument\n", arg);
        return -1;
      }
      options->output = argv[i];
      continue;
    }
    if (!end_options && !strncmp(arg, "--output=", 9)) {
      options->output = arg + 9;
      continue;
    }
    if (!end_options && !strcmp(arg, "--absolute-gate")) {
      if (++i == argc || parse_double(argv[i], &options->absolute_gate)) {
        fprintf(stderr, "sox-quality: invalid absolute gate\n");
        return -1;
      }
      continue;
    }
    if (!end_options && !strncmp(arg, "--absolute-gate=", 16)) {
      if (parse_double(arg + 16, &options->absolute_gate)) {
        fprintf(stderr, "sox-quality: invalid absolute gate\n");
        return -1;
      }
      continue;
    }
    if (!end_options && !strcmp(arg, "--relative-gate")) {
      if (++i == argc || parse_double(argv[i], &options->relative_gate)) {
        fprintf(stderr, "sox-quality: invalid relative gate\n");
        return -1;
      }
      continue;
    }
    if (!end_options && !strncmp(arg, "--relative-gate=", 16)) {
      if (parse_double(arg + 16, &options->relative_gate)) {
        fprintf(stderr, "sox-quality: invalid relative gate\n");
        return -1;
      }
      continue;
    }
    if (!end_options && arg[0] == '-' && arg[1]) {
      fprintf(stderr, "sox-quality: unknown option `%s'\n", arg);
      return -1;
    }
    collect_path(arg, options, files);
  }

  if (!files->len) {
    fprintf(stderr, "sox-quality: no input files\n");
    return -1;
  }
  if (files->len > 1)
    options->batch = 1;
  string_list_sort(files);
  if (options->absolute_gate >= 0 || options->absolute_gate < -120) {
    fprintf(stderr, "sox-quality: absolute gate must be in [-120, 0) LUFS\n");
    return -1;
  }
  if (options->relative_gate <= 0 || options->relative_gate > 40) {
    fprintf(stderr, "sox-quality: relative gate must be in (0, 40] LU\n");
    return -1;
  }
  return 0;
}

int main(int argc, char ** argv)
{
  options_t options;
  string_list_t files = {0};
  file_result_t * results;
  summary_t summary;
  FILE * out;
  size_t i;
  int status = 0;

  if (parse_options(argc, argv, &options, &files) != 0) {
    usage(stderr);
    string_list_free(&files);
    return 1;
  }
  if (sox_init() != SOX_SUCCESS) {
    fprintf(stderr, "sox-quality: cannot initialize libSoX\n");
    string_list_free(&files);
    return 1;
  }
  if (options.quiet)
    sox_globals.verbosity = 0;

  results = xcalloc(files.len, sizeof(*results));
  summary_init(&summary);
  for (i = 0; i < files.len; ++i) {
    if (analyze_file(files.items[i], &options, &results[i], &summary) != 0)
      status = 1;
  }

  out = stdout;
  if (options.output && strcmp(options.output, "-")) {
    out = fopen(options.output, "w");
    if (!out) {
      fprintf(stderr, "sox-quality: cannot open `%s': %s\n",
          options.output, strerror(errno));
      status = 1;
      goto cleanup;
    }
  }
  print_report(out, &options, &files, results, &summary);
  if (out != stdout && fclose(out) != 0) {
    fprintf(stderr, "sox-quality: error closing `%s'\n", options.output);
    status = 1;
  }

cleanup:
  for (i = 0; i < files.len; ++i)
    free(results[i].path);
  free(results);
  summary_free(&summary);
  string_list_free(&files);
  sox_quit();
  return status;
}
