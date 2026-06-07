#!/usr/bin/env bash
# Generate overshoot_correction.lib — a LUT-based overshoot correction for the
# shaped-smoother's undershoot on a 1→0 step, plus quantization helpers so the
# plugin can snap att/shape to the grid values the table was measured at.
#
# Single source of truth: every curve function (shapeMap, cheapCurveBase, the
# env loop, etc.) is imported from the plugin via library(). The ONLY things
# restated here are the grid/block structure and the env loop wiring (which
# must be inlined so that att/shape come from ba.time instead of sliders).
#
# Why a LUT instead of a polynomial: the overshoot surface has discontinuities
# at integer-step boundaries (the number of discrete steps through the curve
# changes abruptly at certain (att, shape) combos).  A polynomial rings around
# these jumps; a table captures them exactly.  Quantizing the plugin's att and
# shape sliders to the grid makes the lookup exact — no interpolation error.
#
# Approach:
#   1. Generate a Faust DSP that sweeps (att, shape) in blocks via ba.time.
#      Uses TARGET_SR (not ma.SR) so the step granularity matches the user's
#      target sample rate regardless of what faust2plot compiles at.
#      Each block: SETTLE samples at input=1, then MEASURE samples at input=0.
#   2. Compile once with faust2plot -double; run once to dump all samples.
#   3. Extract the per-block minimum (last sample of each block).
#   4. Write overshoot_correction.lib with a waveform table + quantizers.
#
# Usage:
#   gen_overshoot_correction.sh PLUGIN.dsp [NATT] [NSHAPE] [SR] [OUT.lib]
#     PLUGIN.dsp  the shapedSmoother core
#     NATT        grid points along attack-time axis       (default 32)
#     NSHAPE      grid points along shape axis             (default 32)
#     SR          target sample rate                       (default 48000)
#     OUT.lib     output library path  (default <plugin_dir>/overshoot_correction.lib)
#   FLOOR_CUTOFF  env var: att_samples above which the LUT value is zeroed
#                 (default 100). Set to 0 to keep all raw values.
set -euo pipefail

plugin=${1:?usage: gen_overshoot_correction.sh PLUGIN.dsp [NATT] [NSHAPE] [SR] [OUT.lib]}
natt=${2:-32}
nshape=${3:-32}
sr=${4:-48000}
plugin_dir=$(cd "$(dirname "$plugin")" && pwd)
plugin_base=$(basename "$plugin")
out=${5:-"$plugin_dir/overshoot_correction.lib"}

att_min_ms="0.046"
att_max_ms="50.0"

max_att_samp=$((sr / 20))
settle=$((max_att_samp * 2))
measure=$((max_att_samp * 2))
block=$((settle + measure))
total_blocks=$((natt * nshape))
total_samples=$((total_blocks * block))

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

echo "target SR=$sr  block=${block} (settle=${settle}+measure=${measure})" >&2
echo "grid ${natt}×${nshape} = ${total_blocks} blocks, ${total_samples} total samples" >&2

# ── 1. Generate the sweep DSP ──────────────────────────────────────────────
cat >"$work/overshoot_sweep.dsp" <<FAUST
import("stdfaust.lib");
ss = library("$plugin_base");

NATT       = $natt;
NSHAPE     = $nshape;
SETTLE     = $settle;
MEASURE    = $measure;
BLOCK      = SETTLE + MEASURE;
TARGET_SR  = $sr.0;

ATT_MIN_MS = $att_min_ms;
ATT_MAX_MS = $att_max_ms;
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

# ── 2. Compile & run ──────────────────────────────────────────────────────
echo "compiling overshoot_sweep.dsp …" >&2
(cd "$work" && faust2plot -double -I "$plugin_dir" overshoot_sweep.dsp) >/dev/null 2>&1

echo "running sweep …" >&2
"$work/overshoot_sweep" -n "$total_samples" >"$work/dump.m"

# ── 3. Extract per-block minimum ──────────────────────────────────────────
grep -E '^[[:space:]]*-?[0-9]' "$work/dump.m" \
  | awk -v blk="$block" -F';' \
        'NR % blk == 0 { gsub(/[ ]/,"",$1); print $1 }' \
        >"$work/mins.txt"

got=$(wc -l < "$work/mins.txt")
echo "extracted $got block minima (expected $total_blocks)" >&2
if [ "$got" -ne "$total_blocks" ]; then
    echo "ERROR: expected $total_blocks minima, got $got — check faust2plot output" >&2
    exit 1
fi

# ── 4. Write the LUT lib ─────────────────────────────────────────────────
PLUGIN="$plugin_base" SRC="$plugin_dir/$plugin_base" \
NATT=$natt NSHAPE=$nshape SR=$sr \
MINS="$work/mins.txt" OUT="$out" \
ATT_MIN_MS=$att_min_ms ATT_MAX_MS=$att_max_ms \
FLOOR_CUTOFF=${FLOOR_CUTOFF:-100} \
python3 - <<'PY'
import os, math, sys, re

natt    = int(os.environ["NATT"])
nshape  = int(os.environ["NSHAPE"])
sr      = int(os.environ["SR"])
att_min_ms = float(os.environ["ATT_MIN_MS"])
att_max_ms = float(os.environ["ATT_MAX_MS"])
floor_cutoff = float(os.environ.get("FLOOR_CUTOFF", "100"))

att_min_samp = att_min_ms / 1000.0 * sr
att_max_samp = att_max_ms / 1000.0 * sr
log_ratio    = math.log(att_max_samp / att_min_samp)

_num = re.compile(r'^[+-]?\d')
raw = [float(l) for l in open(os.environ["MINS"])
       if l.strip() and _num.match(l.strip())]
assert len(raw) == natt * nshape, \
    f"expected {natt*nshape} minima, got {len(raw)}"

overshoot_raw = [max(0.0, -v) for v in raw]

att_norm   = [i / max(1, natt  - 1) for i in range(natt)]
shape_norm = [j / max(1, nshape - 1) for j in range(nshape)]

att_samp_per_idx = [att_min_samp * (att_max_samp / att_min_samp) ** an
                    for an in att_norm]

# Build the LUT: correction multiplier = 1/(1+overshoot).
# Multiply totalStep by this value so the curve lands at the target instead
# of overshooting.  Above the floor cutoff → 1.0 (no correction).
lut = []
for ia in range(natt):
    for js in range(nshape):
        raw_val = overshoot_raw[ia * nshape + js]
        if floor_cutoff > 0 and att_samp_per_idx[ia] > floor_cutoff:
            lut.append(1.0)
        else:
            lut.append(1.0 / (1.0 + raw_val))

max_overshoot = max(overshoot_raw)
min_mult = min(lut)
n_active = sum(1 for s in att_samp_per_idx if floor_cutoff <= 0 or s <= floor_cutoff) * nshape

print(f"  max overshoot measured: {max_overshoot:.6f}", file=sys.stderr)
print(f"  correction mult range: [{min_mult:.6f}, 1.0]", file=sys.stderr)
print(f"  floor cutoff:          {floor_cutoff:.0f} samples "
      f"({n_active} active entries, rest = 1.0)", file=sys.stderr)
print(f"  table size:            {len(lut)} entries "
      f"({len(lut) * 4 / 1024:.1f} KB @ float32)", file=sys.stderr)

# ── Format waveform values ────────────────────────────────────────────────
# Emit as comma-separated lines, ~8 values per line for readability
def fmt_waveform(values, per_line=8):
    lines = []
    for i in range(0, len(values), per_line):
        chunk = values[i:i+per_line]
        lines.append("        " + ", ".join(f"{v:.10e}" for v in chunk))
    return ",\n".join(lines)

# ── Write the lib ────────────────────────────────────────────────────────
fc_note = f"  floor cutoff {floor_cutoff:.0f} samp" if floor_cutoff > 0 else ""
fc_regen = f"FLOOR_CUTOFF={floor_cutoff:.0f} " if floor_cutoff > 0 else ""

L = []
L.append( "// overshoot_correction.lib — LUT-based overshoot correction multiplier.")
L.append( "// GENERATED by gen_overshoot_correction.sh, do not hand-edit.")
L.append(f"// {natt}×{nshape} table ({len(lut)} entries), "
         f"training SR: {sr} Hz.{fc_note}")
L.append(f"// att_samples range [{att_min_samp:.4f}, {att_max_samp:.1f}], "
         f"max overshoot: {max_overshoot:.6f}")
L.append(f"// correction mult = 1/(1+overshoot), range [{min_mult:.6f}, 1.0]")
L.append(f"// source of truth: {os.environ['PLUGIN']} via library().")
L.append(f"// regenerate: {fc_regen}gen_overshoot_correction.sh"
         f" {os.environ['PLUGIN']} {natt} {nshape} {sr}")
L.append( "//")
L.append( "// Usage:")
L.append( "//   import(\"overshoot_correction.lib\");")
L.append( "//   // Snap slider values to the grid the table was measured at:")
L.append( "//   attMs_q = quantizeAttMs(attMs);")
L.append( "//   shape_q = quantizeShape(attackShapeSlider);")
L.append( "//   // Get the correction multiplier (multiply totalStep by this):")
L.append( "//   mult = overshootCorrectionMult(attMs_q * 0.001 * ma.SR, shape_q);")
L.append( "//   // mult ∈ (0, 1]: 1.0 = no correction, <1 = scale down to prevent overshoot")
L.append( "// ============================================================================")
L.append( "")
L.append(f"overshootNATT   = {natt};")
L.append(f"overshootNSHAPE = {nshape};")
L.append(f"overshootATT_MIN_MS  = {att_min_ms};")
L.append(f"overshootATT_MAX_MS  = {att_max_ms};")
L.append(f"overshootATT_MIN_SAMP = {att_min_samp!r};")
L.append(f"overshootATT_MAX_SAMP = {att_max_samp!r};")
L.append(f"overshootLOG_RATIO    = {log_ratio!r};")
L.append( "")
L.append( "// ── LUT lookup ───────────────────────────────────────────────────────────")
L.append( "// Takes att_samples and shapeSlider, quantizes to grid, returns the")
L.append( "// correction multiplier: 1/(1+overshoot).  Multiply totalStep by this.")
L.append( "overshootCorrectionMult(att_samples, shapeSlider) =")
L.append( "    overshootTable, idx : rdtable")
L.append( "with {")
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
L.append( "// Snap attMs to the grid values the table was measured at.")
L.append( "// Use these to drive the smoother so the correction is exact.")
L.append( "quantizeAttMs(attMs) = overshootATT_MIN_MS")
L.append( "    * pow(overshootATT_MAX_MS / overshootATT_MIN_MS,")
L.append( "          float(attIdx) / float(overshootNATT - 1))")
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

# ── Dump the grid for inspection ──────────────────────────────────────────
grid_path = os.environ["OUT"].rsplit(".", 1)[0] + "_grid.tsv"
with open(grid_path, "w") as gf:
    gf.write("attIdx\tshapeIdx\tattMs\tattSamples\tshapeSlider"
             "\tovershoot_raw\tcorrection_mult\n")
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
PY
