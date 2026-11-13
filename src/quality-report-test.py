#!/usr/bin/env python3
"""Integration tests for sox-quality synthetic signals."""

import argparse
import json
import math
import os
import random
import struct
import subprocess
import sys
import tempfile
import wave


def quantize(value, bits):
    scale = 1 << (bits - 1)
    return max(-scale, min(scale - 1, int(round(value * scale))))


def write_wav(path, rate, bits, values):
    width = bits // 8
    packed = bytearray()
    for value in values:
        code = quantize(value, bits)
        if bits == 16:
            packed.extend(struct.pack("<h", code))
        elif bits == 24:
            packed.extend(struct.pack("<i", code)[0:3])
        elif bits == 32:
            packed.extend(struct.pack("<i", code))
        else:
            raise ValueError(bits)
    with wave.open(path, "wb") as output:
        output.setnchannels(1)
        output.setsampwidth(width)
        output.setframerate(rate)
        output.setcomptype("NONE", "not compressed")
        output.writeframes(packed)


def sine_values(rate, seconds, frequency, amplitude):
    return [amplitude * math.sin(2 * math.pi * frequency * n / rate)
            for n in range(int(rate * seconds))]


def clipped_sine_values(rate, seconds, frequency, amplitude):
    values = []
    clipped = 0
    low = -(1 << 15)
    high = (1 << 15) - 1
    for n in range(int(rate * seconds)):
        code = quantize(amplitude * math.sin(2 * math.pi * frequency * n / rate), 16)
        if code == low or code == high:
            clipped += 1
        values.append(code / 32768.0)
    return values, clipped


def run_report(binary, arguments):
    result = subprocess.run(
        [binary] + arguments,
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise AssertionError(
            "sox-quality failed (%d):\n%s\n%s"
            % (result.returncode, result.stdout, result.stderr)
        )
    return json.loads(result.stdout)


def assert_close(actual, expected, tolerance, label):
    if actual is None or not math.isfinite(actual):
        raise AssertionError("%s is not finite: %r" % (label, actual))
    if abs(actual - expected) > tolerance:
        raise AssertionError(
            "%s: expected %.6f +/- %.6f, got %.6f"
            % (label, expected, tolerance, actual)
        )


def get_file(report, path):
    for entry in report["files"]:
        if entry["path"] == path:
            return entry
    raise AssertionError("missing report entry for %s" % path)


def test_sine(binary, directory, path, rate, bits):
    write_wav(path, rate, bits, sine_values(rate, 4.0, 1000, 0.5))
    report = run_report(binary, [path])
    entry = get_file(report, path)
    metrics = entry["metrics"]

    assert_close(metrics["integrated_loudness_lufs"], -9.05, 0.30,
                 "%d-bit/%d Hz integrated loudness" % (bits, rate))
    assert_close(metrics["short_term_loudness_max_lufs"],
                 metrics["integrated_loudness_lufs"], 0.15,
                 "%d-bit/%d Hz short-term maximum" % (bits, rate))
    assert_close(metrics["true_peak_dbtp"], -6.02, 0.20,
                 "%d-bit/%d Hz true peak" % (bits, rate))
    assert_close(metrics["noise_floor_dbfs"], -9.03, 0.08,
                 "%d-bit/%d Hz noise floor" % (bits, rate))
    assert metrics["clipped_samples"] == 0
    assert abs(metrics["dc_offset"]) < 1e-5
    assert entry["format"]["bit_depth"] == bits
    assert entry["format"]["sample_rate_hz"] == rate


def test_white_noise(binary, directory):
    rate = 48000
    bits = 24
    path = os.path.join(directory, "white-noise-24.wav")
    generator = random.Random(20261113)
    values = [generator.uniform(-0.25, 0.25) for _ in range(rate * 4)]
    write_wav(path, rate, bits, values)
    rms = math.sqrt(sum(value * value for value in values) / len(values))
    expected_rms_dbfs = 20 * math.log10(rms)
    expected_peak_dbfs = 20 * math.log10(max(abs(value) for value in values))

    report = run_report(binary, [path])
    metrics = get_file(report, path)["metrics"]

    assert_close(metrics["noise_floor_dbfs"], expected_rms_dbfs, 0.60,
                 "white-noise noise floor")
    assert_close(metrics["integrated_loudness_lufs"],
                 expected_rms_dbfs + 3.10, 0.50,
                 "white-noise integrated loudness")
    assert_close(metrics["true_peak_dbtp"], expected_peak_dbfs + 5.40, 1.00,
                 "white-noise true peak")
    assert metrics["clipped_samples"] == 0
    assert abs(metrics["dc_offset"]) < 1e-3


def test_clipped_signal(binary, directory):
    rate = 48000
    path = os.path.join(directory, "clipped-16.wav")
    values, expected_clipped = clipped_sine_values(rate, 4.0, 1000, 1.5)
    write_wav(path, rate, 16, values)
    report = run_report(binary, [path])
    metrics = get_file(report, path)["metrics"]

    assert metrics["clipped_samples"] == expected_clipped
    assert metrics["clipped_samples"] > 0
    assert -0.25 <= metrics["true_peak_dbtp"] <= 0.10
    assert -2.5 <= metrics["integrated_loudness_lufs"] <= 0.0
    assert abs(metrics["dc_offset"]) < 1e-4


def test_batch(binary, directory, sine_path, noise_path, clipped_path):
    batch_directory = os.path.join(directory, "batch")
    os.mkdir(batch_directory)
    paths = [
        os.path.join(batch_directory, "sine.wav"),
        os.path.join(batch_directory, "noise.wav"),
        os.path.join(batch_directory, "clipped.wav"),
    ]
    for source, destination in zip(
            [sine_path, noise_path, clipped_path], paths):
        with open(source, "rb") as source_file, open(destination, "wb") as target:
            target.write(source_file.read())

    report = run_report(binary, [batch_directory])
    assert report["mode"] == "batch"
    summary = report["summary"]
    file_metrics = [entry["metrics"] for entry in report["files"]]

    assert summary["files_total"] == 3
    assert summary["files_analyzed"] == 3
    assert summary["files_failed"] == 0
    assert summary["files_with_clipping"] == 1
    assert_close(summary["duration_seconds"], 12.0, 0.01, "batch duration")
    assert summary["clipped_samples"] == sum(
        metrics["clipped_samples"] for metrics in file_metrics)
    assert_close(summary["true_peak_dbtp"], max(
        metrics["true_peak_dbtp"] for metrics in file_metrics), 1e-4,
        "batch true peak")
    assert_close(summary["short_term_loudness_max_lufs"], max(
        metrics["short_term_loudness_max_lufs"] for metrics in file_metrics),
        1e-4, "batch short-term maximum")
    assert min(metrics["integrated_loudness_lufs"] for metrics in file_metrics) - 0.1 \
        <= summary["integrated_loudness_lufs"] \
        <= max(metrics["integrated_loudness_lufs"] for metrics in file_metrics) + 0.1
    assert min(metrics["noise_floor_dbfs"] for metrics in file_metrics) - 0.1 \
        <= summary["noise_floor_dbfs"] \
        <= max(metrics["noise_floor_dbfs"] for metrics in file_metrics) + 0.1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bindir", default=".")
    parser.add_argument("--srcdir", default=".")
    args = parser.parse_args()
    binary = os.path.abspath(os.path.join(args.bindir, "sox-quality"))
    if not os.path.exists(binary):
        raise SystemExit("missing binary: %s" % binary)

    with tempfile.TemporaryDirectory(prefix="sox-quality-test-") as directory:
        sine16 = os.path.join(directory, "sine-16-44100.wav")
        sine24 = os.path.join(directory, "sine-24-48000.wav")
        sine32 = os.path.join(directory, "sine-32-96000.wav")
        test_sine(binary, directory, sine16, 44100, 16)
        test_sine(binary, directory, sine24, 48000, 24)
        test_sine(binary, directory, sine32, 96000, 32)
        test_white_noise(binary, directory)
        test_clipped_signal(binary, directory)
        test_batch(
            binary,
            directory,
            sine24,
            os.path.join(directory, "white-noise-24.wav"),
            os.path.join(directory, "clipped-16.wav"),
        )
    print("ok quality-report: sine, white noise, clipping, batch")


if __name__ == "__main__":
    main()
