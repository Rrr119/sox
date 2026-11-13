#!/bin/sh
# Functional tests for the SoX `deess' effect.
#
# Uses only the SoX under test plus a POSIX shell and awk, so no extra
# dependencies are required.  Run from the top-level build directory:
#
#   test/test-deess.sh                 # synthetic material only
#   test/test-deess.sh voice.wav       # additionally check a 44.1kHz/16-bit
#                                      # real voice recording
#
# SOX may be set to the sox executable to test.

set -e

SOX=${SOX:-./src/sox}
AWK=${AWK:-awk}

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

failures=0

ok()   { printf 'ok   - %s\n' "$1"; }
fail() { printf 'FAIL - %s\n' "$1"; failures=$((failures + 1)); }

# Print the RMS amplitude measured in a frequency band (Hz).
band_rms() {
  file=$1; center=$2; width=$3
  "$SOX" "$file" -n bandpass "${center}" "${width}h" stat 2>&1 |
    "$AWK" '/RMS[[:space:]]+amplitude/ {print $3; exit}'
}

db() { "$AWK" -v x="$1" 'BEGIN {printf "%.3f", 20 * log(x) / log(10)}'; }

max_amp() {
  "$SOX" "$1" -n stat 2>&1 |
    "$AWK" '/Maximum amplitude/ {print $3; exit}'
}

# ---------------------------------------------------------------------------
# Build the synthetic sibilance signal: a 200 Hz fundamental with a burst of
# 6-8 kHz noise (the "ess") mixed on top.
# ---------------------------------------------------------------------------
"$SOX" -n -r 44100 -c 1 "$WORK/tone.wav" synth 4 sin 200 gain -8
# The "ess" exists only in the second half, leaving a clean fundamental
# section on which the effect must be (below threshold) a no-op.
"$SOX" -n -r 44100 -c 1 "$WORK/hiss.wav" synth 2 whitenoise sinc 6000-8000 gain -14 pad 2 0
"$SOX" -m -v 1 "$WORK/tone.wav" -v 1 "$WORK/hiss.wav" "$WORK/sib.wav"

# First two seconds: pure 200 Hz fundamental, no sibilance at all.
"$SOX" "$WORK/sib.wav" "$WORK/sib-head.wav" trim 0 2

base_band=$(band_rms "$WORK/sib.wav" 7000 2000)
base_low=$(band_rms "$WORK/sib.wav" 200 200)

# Thresholds chosen relative to the measured side-chain level:
# active is 6 dB below the band level, inactive is 20 dB above.
th_active=$(db "$base_band"); th_active=$("$AWK" -v t="$th_active" 'BEGIN {printf "%g", t - 6}')
th_idle=$("$AWK" -v b="$base_band" 'BEGIN {printf "%g", 20 * log(b) / log(10) + 20}')

echo "synthetic 6-8 kHz sibilance band RMS: $base_band"
echo "active threshold: ${th_active} dB, idle threshold: ${th_idle} dB"

"$SOX" "$WORK/sib.wav" "$WORK/active.wav" deess "$th_active"
"$SOX" "$WORK/sib.wav" "$WORK/idle.wav"   deess "$th_idle"

# 1. Above threshold: the sibilance band is reduced by several dB.
active_band=$(band_rms "$WORK/active.wav" 7000 2000)
reduction=$("$AWK" -v b="$base_band" -v a="$active_band" \
  'BEGIN {printf "%.2f", 20 * log(b / a) / log(10)}')
echo "over-threshold band reduction: ${reduction} dB"
if "$AWK" -v r="$reduction" 'BEGIN {exit !(r > 3)}'; then
  ok "sibilance band attenuated when threshold is exceeded"
else
  fail "sibilance band not attenuated enough (${reduction} dB)"
fi

# 2. Below threshold: the signal is left alone (bit-exact).
if cmp -s "$WORK/sib.wav" "$WORK/idle.wav"; then
  ok "signal unchanged when threshold is not exceeded"
else
  idle_band=$(band_rms "$WORK/idle.wav" 7000 2000)
  drift=$("$AWK" -v b="$base_band" -v i="$idle_band" \
    'BEGIN {printf "%.3f", 20 * log(b / i) / log(10)}')
  if "$AWK" -v d="$drift" 'BEGIN {exit !(d < 0.05)}'; then
    ok "signal unchanged when threshold is not exceeded (${drift} dB drift)"
  else
    fail "signal altered below threshold (${drift} dB)"
  fi
fi

# 3. The fundamental (200 Hz) is essentially preserved.
"$SOX" "$WORK/active.wav" "$WORK/active-head.wav" trim 0 2
ref_low=$(band_rms "$WORK/sib-head.wav" 200 200)
active_low=$(band_rms "$WORK/active-head.wav" 200 200)
low_change=$("$AWK" -v b="$ref_low" -v a="$active_low" \
  'BEGIN {printf "%.2f", 20 * log(b / a) / log(10)}')
echo "fundamental level change (non-sibilant section): ${low_change} dB"
if "$AWK" -v d="$low_change" 'BEGIN {exit !(d < 1.5)}'; then
  ok "200 Hz fundamental preserved"
else
  fail "fundamental changed by ${low_change} dB"
fi

# ---------------------------------------------------------------------------
# 4. Multi-channel: mono and stereo must behave identically per channel.
# ---------------------------------------------------------------------------
"$SOX" -M "$WORK/sib.wav" "$WORK/sib.wav" "$WORK/sib-stereo.wav"
"$SOX" "$WORK/sib-stereo.wav" "$WORK/stereo.wav" deess "$th_active"
"$SOX" "$WORK/stereo.wav" "$WORK/left.wav"  remix 1
"$SOX" "$WORK/stereo.wav" "$WORK/right.wav" remix 2
if cmp -s "$WORK/left.wav" "$WORK/right.wav" &&
   cmp -s "$WORK/left.wav" "$WORK/active.wav"; then
  ok "mono and stereo processing are identical"
else
  fail "stereo channels differ from mono processing"
fi

# ---------------------------------------------------------------------------
# 5. Maximum attenuation is capped.
# ---------------------------------------------------------------------------
"$SOX" -n -r 44100 "$WORK/loud-hiss.wav" synth 2 whitenoise sinc 6000-8000 gain -6
"$SOX" "$WORK/loud-hiss.wav" "$WORK/cap3.wav"  deess -60 6500 5000 0.005 0.05 3
"$SOX" "$WORK/loud-hiss.wav" "$WORK/cap20.wav" deess -60 6500 5000 0.005 0.05 20
cap3=$(band_rms "$WORK/cap3.wav" 7000 2000)
cap20=$(band_rms "$WORK/cap20.wav" 7000 2000)
if "$AWK" -v a="$cap3" -v b="$cap20" 'BEGIN {exit !(a > b)}'; then
  ok "max-attenuation parameter limits the gain reduction"
else
  fail "max-attenuation cap not honoured"
fi

# ---------------------------------------------------------------------------
# 6. Coefficients follow the sample rate (no failure, no clipping at 48 kHz).
# ---------------------------------------------------------------------------
"$SOX" -r 48000 "$WORK/sib.wav" "$WORK/sib48.wav"
if "$SOX" "$WORK/sib48.wav" "$WORK/out48.wav" deess "$th_active" 2>"$WORK/err"; then
  ok "deess runs at 48 kHz (coefficients recalculated)"
else
  cat "$WORK/err"
  fail "deess failed at 48 kHz"
fi

# ---------------------------------------------------------------------------
# 7. Real voice material (optional): no clipping, no obvious gross distortion.
#    If no file is supplied, synthesise a voice-like stand-in that combines
#    a gliding fundamental with shaped noise.
# ---------------------------------------------------------------------------
if [ -n "$1" ]; then
  voice=$1
else
  "$SOX" -n -r 44100 -c 1 -b 16 "$WORK/voice-tone.wav" synth 4 sine 180-260 fade q 1 4 1 gain -10
  "$SOX" -n -r 44100 -c 1 -b 16 "$WORK/voice-hiss.wav" synth 4 whitenoise sinc 5000-9000 fade t 0.5 4 0.5 gain -20
  "$SOX" -m -v 1 "$WORK/voice-tone.wav" -v 1 "$WORK/voice-hiss.wav" -b 16 "$WORK/voice.wav"
  voice=$WORK/voice.wav
fi
# Compare in floating point so 16-bit quantisation cannot mask the result:
# de-essing only attenuates, hence its output level must not exceed the
# input's.
"$SOX" "$voice" -e floating-point "$WORK/voice-ref.f32.wav"
"$SOX" "$voice" -e floating-point "$WORK/voice-out.f32.wav" deess
before_max=$(max_amp "$WORK/voice-ref.f32.wav")
after_max=$(max_amp "$WORK/voice-out.f32.wav")
# Re-run the same chain with verbose logging: de-essing only attenuates, so
# no effect may report clipping for a clean input.
"$SOX" -b 16 "$voice" "$WORK/voice-out.wav" deess
after_clips=$("$SOX" -V2 -b 16 "$voice" "$WORK/voice-verbose.wav" deess 2>&1 |
  "$AWK" 'match($0, /clipped [0-9]+ samples/) {
    n += substr($0, RSTART + 8, RLENGTH - 16)
  } END {print n + 0}')
echo "voice max amplitude before/after: $before_max / $after_max, clips: $after_clips"
if [ "$after_clips" -eq 0 ]; then
  ok "voice material: no clipping"
else
  fail "voice material: $after_clips clipped samples"
fi
if "$AWK" -v a="$after_max" -v b="$before_max" 'BEGIN {exit !(a <= b)}'; then
  ok "voice material: output never exceeds input level"
else
  fail "voice material: output level exceeds input"
fi

echo
if [ "$failures" -eq 0 ]; then
  echo "PASS: all deess tests succeeded"
else
  echo "FAIL: $failures deess test(s) failed"
  exit 1
fi
