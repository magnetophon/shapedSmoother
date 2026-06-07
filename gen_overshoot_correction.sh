#!/usr/bin/env bash
# Generate overshoot_correction.lib — a LUT-based overshoot correction
# multiplier, measured in two passes:
#   Pass 1: coarse full-range sweep to find the boundary where overshoot
#           becomes negligible (the ~0.01 floor).
#   Pass 2: fine sweep zoomed into just the active region, so the full
#           NATT×NSHAPE grid resolution covers only the range that matters.
#
# The lib exports range constants so the plugin knows where the table
# applies.  Values outside the table's range get 1.0 (no correction).
#
# Usage:
#   gen_overshoot_correction.sh PLUGIN.dsp [NATT] [NSHAPE] [SR] [OUT.lib]
#     PLUGIN.dsp  the shapedSmoother core
#     NATT        grid points per axis (final table)           (default 32)
#     NSHAPE      grid points per axis (final table)           (default 32)
#     SR          target sample rate                           (default 48000)
#     OUT.lib     output library path
#   FLOOR_MULT    env var: overshoot threshold relative to the measured floor
#                 for boundary detection (default 2.0: active = >2× floor)
set -euo pipefail

plugin=${1:?usage: gen_overshoot_correction.sh PLUGIN.dsp [NATT] [NSHAPE] [SR] [OUT.lib]}
natt=${2:-32}
nshape=${3:-32}
sr=${4:-48000}
plugin_dir=$(cd "$(dirname "$plugin")" && pwd)
plugin_base=$(basename "$plugin")
out=${5:-"$plugin_dir/overshoot_correction.lib"}
floor_mult=${FLOOR_MULT:-2.0}

att_min_ms="0.046"
att_max_ms="50.0"

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# ── Shared: generate a sweep DSP ─────────────────────────────────────────
# Args: $1=output.dsp $2=natt $3=nshape $4=att_min_ms $5=att_max_ms
gen_sweep_dsp() {
    local dsp=$1 na=$2 ns=$3 amin=$4 amax=$5
    local max_att_samp=$(( sr * ${amax%.*} / 1000 + sr / 1000 ))  # generous ceiling
    local stl=$(( max_att_samp > 2400 ? max_att_samp * 2 : 4800 ))
    local mea=$stl
    cat >"$dsp" <<FAUST
import("stdfaust.lib");
ss = library("$plugin_base");

NATT       = $na;
NSHAPE     = $ns;
SETTLE     = $stl;
MEASURE    = $mea;
BLOCK      = SETTLE + MEASURE;
TARGET_SR  = $sr.0;

ATT_MIN_MS = $amin;
ATT_MAX_MS = $amax;
REL_MS     = 0.5;

MAX_HOLD   = 0.05;
MAX_HOLD_SAMP = int(MAX_HOLD * TARGET_SR);

blockIdx   = int(ba.time) / BLOCK;
posInBlock = int(ba.time) % BLOCK;
attIdx     = min(NATT  - 1, blockIdx / NSHAPE);
shapeIdx   = min(NSHAPE - 1, blockIdx - attIdx * NSHAPE);

currentAttMs = ATT_MIN_MS
             * pow(ATT_MAX_MS / ATT_MIN_MS,
                   float(attIdx) / float(max(1, NATT - 1)));
currentShape = float(shapeIdx) / float(max(1, NSHAPE - 1));

stepInput = 1 - (posInBlock >= SETTLE);

paramSmoother(attMsArg, shapeSlArg, x) =
    lookaheadX, delayedX : env ~ (_, _, _) : (_, !, !)
with {
    attP  = attMsArg / 1000.0;
    relP  = REL_MS   / 1000.0;
    attSamp = attP * TARGET_SR : max(1) : min(MAX_HOLD_SAMP);
    lookaheadX = x : ba.slidingMin(int(attSamp) + 1, 1 + MAX_HOLD_SAMP);
    delayedX   = x @ int(attSamp);

    aShape = ss.shapeMap(shapeSlArg);
    rShape = ss.shapeMap(0.0);
    aICS = 1.0 / ss.curveScale(aShape);
    rICS = 1.0 / ss.curveScale(rShape);
    aZero = ss.cheapCurveBase(aShape, 0);
    rZero = ss.cheapCurveBase(rShape, 0);
    aMDB = ss.maxDerivativeBaseAttack(aShape);
    rMDB = ss.maxDerivativeBaseAttack(rShape);
    attStep = 1.0 / (attP * TARGET_SR : max(1));
    relStep = 1.0 / (relP * TARGET_SR : max(1));

    env(prev, prevPhase, prevTotalStep, lhx, dx) = result, newPhase, totalStep
    with {
        attacking = (prev > lhx) & (attP > 0);
        releasing = (prev < lhx) & (relP > 0);
        active    = attacking | releasing;
        shape = select2(releasing, aShape, rShape);
        ics   = aICS + releasing * (rICS - aICS);
        zv    = select2(releasing, aZero, rZero);
        mdv   = select2(releasing, aMDB, rMDB);
        stp   = select2(releasing, attStep, relStep);
        cbPrev   = (ss.cheapCurveBase(shape,
                       select2(releasing, 1 - prevPhase, prevPhase))
                    - zv) * ics;
        fracDone = select2(releasing, 1 - cbPrev, cbPrev);
        rawTS = (lhx - prev) + prevTotalStep * fracDone;
        needsRecompute = (lhx != lhx') | (prevTotalStep <= 0);
        totalStep = select2(releasing,
                        rawTS : min(prevTotalStep),
                        select2(needsRecompute, prevTotalStep, rawTS))
                    * active;
        prevSpeed  = prev - prev';
        speedRatio = prevSpeed
                   / (totalStep * ics * stp + (1 - active) * 1e-30);
        cr = max(0, min(speedRatio, mdv));
        gd(ph) = select2(releasing, cb, 1 - cb) * totalStep
        with {
            cb = (ss.cheapCurveBase(shape,
                     select2(releasing, 1 - ph, ph)) - zv) * ics;
        };
        proj = gd(select2(releasing,
                   ss.inverseDerivativeBottomAttack(shape, cr),
                   ss.inverseDerivativeTopRelease(shape, cr))) + prev;
        gmi  = proj > lhx;
        pms = select2(releasing,
                select2(gmi,
                    ss.inverseDerivativeBottomAttack(shape, cr),
                    ss.inverseDerivativeTopAttack(shape, cr)),
                select2(gmi,
                    ss.inverseDerivativeBottomRelease(shape, cr),
                    ss.inverseDerivativeTopRelease(shape, cr)));
        newPhase = (pms + stp) : min(1 - stp) : max(stp) * active;
        spd   = totalStep
              * ss.derivativeBaseRelease(shape,
                    select2(releasing, 1 - newPhase, newPhase))
              * ics;
        delta = spd * stp;
        result = min(prev + delta, dx);
    };
};

smootherOut = stepInput : paramSmoother(currentAttMs, currentShape);
blockMin = smootherOut : tracker ~ _
with {
    inMeas = posInBlock >= SETTLE;
    reset  = posInBlock == SETTLE;
    tracker(prev, x) = select2(reset, select2(inMeas, prev, min(prev, x)), x);
};
process = blockMin;
FAUST
}

# ── Shared: compile, run, extract ─────────────────────────────────────────
# Args: $1=name $2=dsp $3=natt $4=nshape $5=settle $6=measure
run_sweep() {
    local name=$1 dsp=$2 na=$3 ns=$4
    local total=$(( na * ns ))
    # read SETTLE and MEASURE from the DSP (they're literal integers)
    local stl=$(grep -oP 'SETTLE\s*=\s*\K[0-9]+' "$dsp" | head -1)
    local mea=$(grep -oP 'MEASURE\s*=\s*\K[0-9]+' "$dsp" | head -1)
    local blk=$(( stl + mea ))
    local samples=$(( total * blk ))

    echo "  [$name] compiling …" >&2
    local dir=$(dirname "$dsp")
    local base=$(basename "$dsp")
    (cd "$dir" && faust2plot -double -I "$plugin_dir" "$base") >/dev/null 2>&1
    local bin="${dsp%.dsp}"

    echo "  [$name] running ($total blocks × $blk samples) …" >&2
    "$bin" -n "$samples" >"${bin}.m"

    echo "  [$name] extracting …" >&2
    grep -E '^[[:space:]]*-?[0-9]' "${bin}.m" \
      | awk -v blk="$blk" -F';' \
            'NR % blk == 0 { gsub(/[ ]/,"",$1); print $1 }' \
            >"${bin}_mins.txt"

    local got=$(wc -l < "${bin}_mins.txt")
    if [ "$got" -ne "$total" ]; then
        echo "  ERROR: expected $total minima, got $got" >&2; exit 1
    fi
    echo "  [$name] $got minima extracted" >&2
}

# ══════════════════════════════════════════════════════════════════════════
#  PASS 1 — coarse full-range sweep to find the active boundary
# ══════════════════════════════════════════════════════════════════════════
echo "=== Pass 1: full-range boundary detection ===" >&2
nshape_p1=8    # few shape points — enough to detect overshoot at any shape
gen_sweep_dsp "$work/p1_sweep.dsp" "$natt" "$nshape_p1" "$att_min_ms" "$att_max_ms"
run_sweep "pass1" "$work/p1_sweep.dsp" "$natt" "$nshape_p1"

# Analyze pass 1: find the att boundary
boundary_att_ms=$(
NATT=$natt NSHAPE_P1=$nshape_p1 SR=$sr FLOOR_MULT=$floor_mult \
ATT_MIN_MS=$att_min_ms ATT_MAX_MS=$att_max_ms \
MINS="$work/p1_sweep_mins.txt" \
python3 - <<'PY1'
import os, math, re, sys

natt = int(os.environ["NATT"])
nshape = int(os.environ["NSHAPE_P1"])
sr = int(os.environ["SR"])
floor_mult = float(os.environ["FLOOR_MULT"])
att_min = float(os.environ["ATT_MIN_MS"])
att_max = float(os.environ["ATT_MAX_MS"])

_num = re.compile(r'^[+-]?\d')
raw = [float(l) for l in open(os.environ["MINS"])
       if l.strip() and _num.match(l.strip())]
assert len(raw) == natt * nshape

overshoot = [max(0.0, -v) for v in raw]

# Max overshoot per attIdx
max_per_att = []
for ia in range(natt):
    mx = max(overshoot[ia * nshape + js] for js in range(nshape))
    max_per_att.append(mx)

# Floor estimate: median of the last quarter of attIdx values
tail = sorted(max_per_att[-(natt // 4):])
floor_val = tail[len(tail) // 2]

threshold = floor_val * floor_mult
print(f"  floor estimate: {floor_val:.6f}, threshold (×{floor_mult}): {threshold:.6f}",
      file=sys.stderr)

# Last attIdx above threshold
boundary_idx = 0
for ia in range(natt):
    if max_per_att[ia] > threshold:
        boundary_idx = ia

# Convert to attMs (add 1 index of margin, capped)
boundary_idx = min(boundary_idx + 1, natt - 1)
att_norm = boundary_idx / max(1, natt - 1)
boundary_ms = att_min * (att_max / att_min) ** att_norm
boundary_samp = boundary_ms / 1000.0 * sr

print(f"  active boundary: attIdx={boundary_idx}, "
      f"attMs={boundary_ms:.4f}, attSamp={boundary_samp:.1f}",
      file=sys.stderr)

# Emit just the boundary ms (stdout → captured by bash)
print(f"{boundary_ms:.6f}")
PY1
)

echo "  boundary_att_ms=$boundary_att_ms" >&2

# ══════════════════════════════════════════════════════════════════════════
#  PASS 2 — fine sweep zoomed into the active region
# ══════════════════════════════════════════════════════════════════════════
echo "=== Pass 2: fine sweep in [$att_min_ms, $boundary_att_ms] ms ===" >&2
gen_sweep_dsp "$work/p2_sweep.dsp" "$natt" "$nshape" "$att_min_ms" "$boundary_att_ms"
run_sweep "pass2" "$work/p2_sweep.dsp" "$natt" "$nshape"

# ══════════════════════════════════════════════════════════════════════════
#  Emit the LUT lib
# ══════════════════════════════════════════════════════════════════════════
echo "=== Building lib ===" >&2
PLUGIN="$plugin_base" \
NATT=$natt NSHAPE=$nshape SR=$sr \
ATT_MIN_MS=$att_min_ms ATT_MAX_MS=$boundary_att_ms \
SLIDER_ATT_MAX_MS=$att_max_ms \
MINS="$work/p2_sweep_mins.txt" OUT="$out" \
python3 - <<'PY2'
import os, math, sys, re

natt    = int(os.environ["NATT"])
nshape  = int(os.environ["NSHAPE"])
sr      = int(os.environ["SR"])
att_min_ms = float(os.environ["ATT_MIN_MS"])
att_max_ms = float(os.environ["ATT_MAX_MS"])       # table boundary (zoomed)
slider_att_max_ms = float(os.environ["SLIDER_ATT_MAX_MS"])  # full slider max

att_min_samp = att_min_ms / 1000.0 * sr
att_max_samp = att_max_ms / 1000.0 * sr
log_ratio    = math.log(att_max_samp / att_min_samp)

_num = re.compile(r'^[+-]?\d')
raw = [float(l) for l in open(os.environ["MINS"])
       if l.strip() and _num.match(l.strip())]
assert len(raw) == natt * nshape, \
    f"expected {natt*nshape} minima, got {len(raw)}"

overshoot_raw = [max(0.0, -v) for v in raw]

# Build LUT: correction multiplier = 1/(1+overshoot)
lut = [1.0 / (1.0 + os) for os in overshoot_raw]

max_overshoot = max(overshoot_raw)
min_mult = min(lut)

print(f"  max overshoot:  {max_overshoot:.6f}", file=sys.stderr)
print(f"  mult range:     [{min_mult:.6f}, {max(lut):.6f}]", file=sys.stderr)
print(f"  table range:    att [{att_min_samp:.4f}, {att_max_samp:.1f}] samp, "
      f"shape [0, 1]", file=sys.stderr)
print(f"  table size:     {len(lut)} entries "
      f"({len(lut) * 4 / 1024:.1f} KB @ float32)", file=sys.stderr)

def fmt_waveform(values, per_line=8):
    lines = []
    for i in range(0, len(values), per_line):
        chunk = values[i:i+per_line]
        lines.append("        " + ", ".join(f"{v:.10e}" for v in chunk))
    return ",\n".join(lines)

# ── Write the lib ────────────────────────────────────────────────────────
L = []
L.append( "// overshoot_correction.lib — LUT-based overshoot correction multiplier.")
L.append( "// GENERATED by gen_overshoot_correction.sh (two-pass), do not hand-edit.")
L.append(f"// {natt}×{nshape} table ({len(lut)} entries), training SR: {sr} Hz")
L.append(f"// table covers att [{att_min_samp:.4f}, {att_max_samp:.1f}] samp "
         f"({att_min_ms} – {att_max_ms:.4f} ms)")
L.append(f"// max overshoot: {max_overshoot:.6f}, "
         f"mult range: [{min_mult:.6f}, 1.0]")
L.append(f"// source of truth: {os.environ['PLUGIN']} via library().")
L.append(f"// regenerate: gen_overshoot_correction.sh"
         f" {os.environ['PLUGIN']} {natt} {nshape} {sr}")
L.append( "//")
L.append( "// Usage:")
L.append( "//   import(\"overshoot_correction.lib\");")
L.append( "//   attMs_q = quantizeAttMs(attMs);")
L.append( "//   shape_q = quantizeShape(attackShapeSlider);")
L.append( "//   mult    = overshootCorrectionMult(attMs_q * 0.001 * ma.SR, shape_q);")
L.append( "//   // mult ∈ (0, 1]: 1.0 = no correction, <1 = scale down totalStep")
L.append( "//   // Returns 1.0 for att_samples above the table range.")
L.append( "// ============================================================================")
L.append( "")
L.append(f"overshootNATT   = {natt};")
L.append(f"overshootNSHAPE = {nshape};")
L.append( "")
L.append( "// ── Table range (the region where correction != 1.0) ─────────────────────")
L.append(f"overshootATT_MIN_MS   = {att_min_ms};")
L.append(f"overshootATT_MAX_MS   = {att_max_ms:.10f};")
L.append(f"overshootATT_MIN_SAMP = {att_min_samp!r};")
L.append(f"overshootATT_MAX_SAMP = {att_max_samp!r};")
L.append(f"overshootLOG_RATIO    = {log_ratio!r};")
L.append(f"overshootSHAPE_MIN    = 0.0;")
L.append(f"overshootSHAPE_MAX    = 1.0;")
L.append( "")
L.append( "// ── Full slider range (for the quantizers) ────────────────────────────────")
L.append(f"overshootSLIDER_ATT_MIN_MS = {att_min_ms};")
L.append(f"overshootSLIDER_ATT_MAX_MS = {slider_att_max_ms};")
L.append( "")
L.append( "// ── LUT lookup ───────────────────────────────────────────────────────────")
L.append( "overshootCorrectionMult(att_samples, shapeSlider) =")
L.append( "    select2(att_samples > overshootATT_MAX_SAMP, tableVal, 1.0)")
L.append( "with {")
L.append( "    tableVal = overshootTable, idx : rdtable;")
L.append( "    attIdx = int(0.5 + log(max(overshootATT_MIN_SAMP,")
L.append( "                              min(overshootATT_MAX_SAMP, att_samples))")
L.append( "                          / overshootATT_MIN_SAMP)")
L.append( "                      / overshootLOG_RATIO * (overshootNATT - 1))")
L.append( "           : max(0) : min(overshootNATT - 1);")
L.append( "    shapeIdx = int(0.5 + (shapeSlider : max(0) : min(1))")
L.append( "                        * (overshootNSHAPE - 1))")
L.append( "             : max(0) : min(overshootNSHAPE - 1);")
L.append( "    idx = attIdx * overshootNSHAPE + shapeIdx;")
L.append( "    overshootTable = waveform {")
L.append(fmt_waveform(lut))
L.append( "    };")
L.append( "};")
L.append( "")
L.append( "// ── Quantization helpers ──────────────────────────────────────────────────")
L.append( "// Snap attMs to the LUT grid. For att above the table range,")
L.append( "// returns the raw value unchanged (no correction needed there).")
L.append( "quantizeAttMs(attMs) = select2(attMs > overshootATT_MAX_MS,")
L.append( "    overshootATT_MIN_MS")
L.append( "        * pow(overshootATT_MAX_MS / overshootATT_MIN_MS,")
L.append( "              float(attIdx) / float(overshootNATT - 1)),")
L.append( "    attMs)")
L.append( "with {")
L.append( "    attIdx = int(0.5 + log(max(overshootATT_MIN_MS,")
L.append( "                              min(overshootATT_MAX_MS, attMs))")
L.append( "                          / overshootATT_MIN_MS)")
L.append( "                      / log(overshootATT_MAX_MS / overshootATT_MIN_MS)")
L.append( "                      * (overshootNATT - 1))")
L.append( "           : max(0) : min(overshootNATT - 1);")
L.append( "};")
L.append( "")
L.append( "quantizeShape(s) = float(shapeIdx) / float(overshootNSHAPE - 1)")
L.append( "with {")
L.append( "    shapeIdx = int(0.5 + (s : max(0) : min(1)) * (overshootNSHAPE - 1))")
L.append( "             : max(0) : min(overshootNSHAPE - 1);")
L.append( "};")

open(os.environ["OUT"], "w").write("\n".join(L) + "\n")
print(f"wrote {os.environ['OUT']}  ({natt}×{nshape} = {len(lut)} entries)",
      file=sys.stderr)

# ── Dump grid for inspection ──────────────────────────────────────────────
grid_path = os.environ["OUT"].rsplit(".", 1)[0] + "_grid.tsv"
with open(grid_path, "w") as gf:
    gf.write("attIdx\tshapeIdx\tattMs\tattSamples\tshapeSlider"
             "\tovershoot\tcorrection_mult\n")
    for ia in range(natt):
        att_ms = att_min_ms * (att_max_ms / att_min_ms) ** (
            ia / max(1, natt - 1))
        att_samp = att_ms / 1000.0 * sr
        for js in range(nshape):
            shape_sl = js / max(1, nshape - 1)
            k = ia * nshape + js
            gf.write(f"{ia}\t{js}\t{att_ms:.4f}\t{att_samp:.2f}"
                     f"\t{shape_sl:.4f}"
                     f"\t{overshoot_raw[k]:.8e}\t{lut[k]:.8e}\n")
print(f"wrote {grid_path}", file=sys.stderr)
PY2
