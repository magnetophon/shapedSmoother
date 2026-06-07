#!/usr/bin/env bash
# Generate overshoot_correction.lib — a LUT-based overshoot correction
# multiplier for the shaped-smoother envelope follower.
#
# Two modes, auto-detected:
#
#   INITIAL (no lib exists):  Two-pass sweep — coarse full-range to find the
#       boundary, fine zoomed pass for the table — outputs the first lib.
#
#   ITERATIVE (lib exists):  Sweeps with the correction applied, measures the
#       residual overshoot, and tightens the LUT:
#           new_mult[k] = old_mult[k] / (1 + residual_os[k])
#       Repeats until max residual < CONV_THRESHOLD or MAX_ITERS reached.
#
# Usage:
#   gen_overshoot_correction.sh PLUGIN.dsp [NATT] [NSHAPE] [SR] [OUT.lib]
#   MAX_ITERS=5 CONV_THRESHOLD=0.001 ./gen_overshoot_correction.sh PLUGIN.dsp
#   FLOOR_MULT=2.0 ./gen_overshoot_correction.sh PLUGIN.dsp   # initial only
set -euo pipefail

plugin=${1:?usage: gen_overshoot_correction.sh PLUGIN.dsp [NATT] [NSHAPE] [SR] [OUT.lib]}
natt=${2:-32}
nshape=${3:-32}
sr=${4:-48000}
plugin_dir=$(cd "$(dirname "$plugin")" && pwd)
plugin_base=$(basename "$plugin")
out=${5:-"$plugin_dir/overshoot_correction.lib"}
out_dir=$(cd "$(dirname "$out")" && pwd)

att_min_ms="0.046"
att_max_ms="50.0"

work=$(mktemp -d)
if [ "${KEEP_WORK:-}" = "1" ]; then
    echo "work dir (kept): $work" >&2
else
    trap 'rm -rf "$work"' EXIT
fi

# ── Generate a sweep DSP ─────────────────────────────────────────────────
# $1=dsp $2=natt $3=nshape $4=att_min_ms $5=att_max_ms $6=correction_lib(opt)
gen_sweep_dsp() {
    local dsp=$1 na=$2 ns=$3 amin=$4 amax=$5 correction_lib=${6:-}
    local max_att_samp=$(( sr * ${amax%.*} / 1000 + sr / 1000 ))
    local stl=$(( max_att_samp > 2400 ? max_att_samp * 2 : 4800 ))
    local mea=$stl

    local oc_import=""
    local oc_mult="1.0"
    if [ -n "$correction_lib" ]; then
        oc_import="oc = library(\"$(basename "$correction_lib")\");"
        oc_mult="oc.overshootCorrectionMult(attMsArg * 0.001 * TARGET_SR, shapeSlArg)"
    fi

    cat >"$dsp" <<FAUST
import("stdfaust.lib");
ss = library("$plugin_base");
$oc_import

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
    lookaheadX, delayedX, corrMult : env ~ (_, _, _) : (_, !, !)
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

    corrMult = $oc_mult;

    env(prev, prevPhase, prevTotalStep, lhx, dx, cm) = result, newPhase, totalStep
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
        // correctedTS: scale the span for speed, velocity matching, and
        // projection, but keep uncorrected totalStep for feedback so the
        // correction doesn't compound across samples.
        correctedTS = totalStep * cm;
        prevSpeed  = prev - prev';
        speedRatio = prevSpeed
                   / (correctedTS * ics * stp + (1 - active) * 1e-30);
        cr = max(0, min(speedRatio, mdv));
        gd(ph) = select2(releasing, cb, 1 - cb) * correctedTS
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
        spd   = correctedTS
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

# ── Compile, run, extract ─────────────────────────────────────────────────
# $1=name $2=dsp $3=natt $4=nshape $5=extra_include(opt)
run_sweep() {
    local name=$1 dsp=$2 na=$3 ns=$4 extra_include=${5:-}
    local total=$(( na * ns ))
    local stl=$(grep -oP 'SETTLE\s*=\s*\K[0-9]+' "$dsp" | head -1)
    local mea=$(grep -oP 'MEASURE\s*=\s*\K[0-9]+' "$dsp" | head -1)
    local blk=$(( stl + mea ))
    local samples=$(( total * blk ))
    local dir=$(dirname "$dsp")
    local base=$(basename "$dsp")

    local extra_I=""
    [ -n "$extra_include" ] && extra_I="-I $extra_include"

    echo "  [$name] compiling …" >&2
    (cd "$dir" && faust2plot -double -I "$plugin_dir" $extra_I "$base") >/dev/null

    echo "  [$name] running ($total blocks × $blk) …" >&2
    "${dsp%.dsp}" -n "$samples" >"${dsp%.dsp}.m"

    echo "  [$name] extracting …" >&2
    grep -E '^[[:space:]]*-?[0-9]' "${dsp%.dsp}.m" \
      | awk -v blk="$blk" -F';' \
            'NR % blk == 0 { gsub(/[ ]/,"",$1); print $1 }' \
            >"${dsp%.dsp}_mins.txt"

    local got=$(wc -l < "${dsp%.dsp}_mins.txt")
    if [ "$got" -ne "$total" ]; then
        echo "  ERROR: expected $total minima, got $got" >&2; exit 1
    fi
    echo "  [$name] $got minima" >&2
}

# ── Write a lib from Python ──────────────────────────────────────────────
# Shared Python for writing the lib.  Expects env vars:
#   NATT NSHAPE SR ATT_MIN_MS ATT_MAX_MS PLUGIN MINS OUT
#   Optional: OLD_LIB (path to existing lib for iterative update)
write_lib() {
    python3 - <<'PYLIB'
import os, math, sys, re

natt    = int(os.environ["NATT"])
nshape  = int(os.environ["NSHAPE"])
sr      = int(os.environ["SR"])
att_min_ms = float(os.environ["ATT_MIN_MS"])
att_max_ms = float(os.environ["ATT_MAX_MS"])
old_lib = os.environ.get("OLD_LIB", "")

att_min_samp = att_min_ms / 1000.0 * sr
att_max_samp = att_max_ms / 1000.0 * sr
log_ratio    = math.log(att_max_samp / att_min_samp)

_num = re.compile(r'^[+-]?\d')
raw = [float(l) for l in open(os.environ["MINS"])
       if l.strip() and _num.match(l.strip())]
assert len(raw) == natt * nshape, \
    f"expected {natt*nshape} minima, got {len(raw)}"
overshoot_raw = [max(0.0, -v) for v in raw]

if old_lib and os.path.isfile(old_lib):
    # ── Iterative: read old LUT, compute new = old / (1 + residual) ──
    lib_text = open(old_lib).read()
    wf_match = re.search(r'waveform\s*\{([^}]+)\}', lib_text)
    old_lut = [float(x) for x in re.findall(
        r'[-+]?\d+\.\d+[eE][+-]?\d+', wf_match.group(1))]
    assert len(old_lut) == natt * nshape, \
        f"old LUT has {len(old_lut)} entries, expected {natt*nshape}"
    lut = [old / (1.0 + res) for old, res in zip(old_lut, overshoot_raw)]
    max_residual = max(overshoot_raw)
    print(f"  max residual overshoot: {max_residual:.6e}", file=sys.stderr)
else:
    # ── Initial: correction = 1/(1+overshoot) ──
    lut = [1.0 / (1.0 + os) for os in overshoot_raw]
    max_residual = max(overshoot_raw)
    print(f"  max overshoot: {max_residual:.6f}", file=sys.stderr)

min_mult = min(lut)
print(f"  mult range: [{min_mult:.6f}, {max(lut):.6f}]", file=sys.stderr)
print(f"  table: {natt}×{nshape} = {len(lut)} entries", file=sys.stderr)

# Show the worst residual points
worst = sorted(range(len(overshoot_raw)), key=lambda k: -overshoot_raw[k])[:10]
print(f"  worst residual points:", file=sys.stderr)
for rank, k in enumerate(worst):
    ia, js = k // nshape, k % nshape
    a_ms = att_min_ms * (att_max_ms / att_min_ms) ** (ia / max(1, natt-1))
    a_samp = a_ms / 1000.0 * sr
    shape = js / max(1, nshape-1)
    old_m = old_lut[k] if (old_lib and os.path.isfile(old_lib)) else 1.0/(1.0+overshoot_raw[k])
    print(f"    #{rank+1:2d}  att={a_samp:6.2f} samp  shape={shape:.3f}  "
          f"residual={overshoot_raw[k]:.4f}  old_mult={old_m:.4f}  new_mult={lut[k]:.4f}",
          file=sys.stderr)

def fmt_waveform(values, per_line=8):
    lines = []
    for i in range(0, len(values), per_line):
        chunk = values[i:i+per_line]
        lines.append("        " + ", ".join(f"{v:.10e}" for v in chunk))
    return ",\n".join(lines)

L = []
L.append( "// overshoot_correction.lib — LUT-based overshoot correction multiplier.")
L.append( "// GENERATED by gen_overshoot_correction.sh, do not hand-edit.")
L.append(f"// {natt}×{nshape} table ({len(lut)} entries), training SR: {sr} Hz")
L.append(f"// table covers att [{att_min_samp:.4f}, {att_max_samp:.1f}] samp "
         f"({att_min_ms} – {att_max_ms:.4f} ms)")
L.append(f"// mult range: [{min_mult:.6f}, 1.0]")
L.append(f"// source of truth: {os.environ['PLUGIN']} via library().")
L.append(f"// regenerate: gen_overshoot_correction.sh"
         f" {os.environ['PLUGIN']} {natt} {nshape} {sr}")
L.append( "//")
L.append( "// Usage:")
L.append( "//   import(\"overshoot_correction.lib\");")
L.append( "//   att_samples_q = quantizeAttSamples(attMs * 0.001 * ma.SR);")
L.append( "//   shape_q       = quantizeShape(attackShapeSlider);")
L.append( "//   mult          = overshootCorrectionMult(att_samples_q, shape_q);")
L.append( "//   att           = att_samples_q / ma.SR;")
L.append( "// ============================================================================")
L.append( "")
L.append(f"overshootNATT   = {natt};")
L.append(f"overshootNSHAPE = {nshape};")
L.append( "")
L.append(f"overshootATT_MIN_SAMP = {att_min_samp!r};")
L.append(f"overshootATT_MAX_SAMP = {att_max_samp!r};")
L.append(f"overshootLOG_RATIO    = {log_ratio!r};")
L.append(f"overshootSHAPE_MIN    = 0.0;")
L.append(f"overshootSHAPE_MAX    = 1.0;")
L.append( "")
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
L.append( "quantizeAttSamples(att_samples) = select2(att_samples > overshootATT_MAX_SAMP,")
L.append( "    overshootATT_MIN_SAMP")
L.append( "        * pow(overshootATT_MAX_SAMP / overshootATT_MIN_SAMP,")
L.append( "              float(attIdx) / float(overshootNATT - 1)),")
L.append( "    att_samples)")
L.append( "with {")
L.append( "    attIdx = int(0.5 + log(max(overshootATT_MIN_SAMP,")
L.append( "                              min(overshootATT_MAX_SAMP, att_samples))")
L.append( "                          / overshootATT_MIN_SAMP)")
L.append( "                      / overshootLOG_RATIO * (overshootNATT - 1))")
L.append( "           : max(0) : min(overshootNATT - 1);")
L.append( "};")
L.append( "")
L.append( "quantizeShape(s) = float(shapeIdx) / float(overshootNSHAPE - 1)")
L.append( "with {")
L.append( "    shapeIdx = int(0.5 + (s : max(0) : min(1)) * (overshootNSHAPE - 1))")
L.append( "             : max(0) : min(overshootNSHAPE - 1);")
L.append( "};")

open(os.environ["OUT"], "w").write("\n".join(L) + "\n")
print(f"wrote {os.environ['OUT']}", file=sys.stderr)

# Grid TSV
grid_path = os.environ["OUT"].rsplit(".", 1)[0] + "_grid.tsv"
with open(grid_path, "w") as gf:
    gf.write("attIdx\tshapeIdx\tattMs\tattSamples\tshapeSlider"
             "\tovershoot\tcorrection_mult\n")
    for ia in range(natt):
        a_ms = att_min_ms * (att_max_ms / att_min_ms) ** (
            ia / max(1, natt - 1))
        a_samp = a_ms / 1000.0 * sr
        for js in range(nshape):
            s = js / max(1, nshape - 1)
            k = ia * nshape + js
            gf.write(f"{ia}\t{js}\t{a_ms:.4f}\t{a_samp:.2f}\t{s:.4f}"
                     f"\t{overshoot_raw[k]:.8e}\t{lut[k]:.8e}\n")
print(f"wrote {grid_path}", file=sys.stderr)

# Report convergence status (stdout → captured by bash)
conv = float(os.environ.get("CONV_THRESHOLD", "0.001"))
if max_residual < conv:
    print("converged")
else:
    print("not_converged")
PYLIB
}

# ══════════════════════════════════════════════════════════════════════════
if [ -f "$out" ]; then
# ══════════════════════════════════════════════════════════════════════════
#  ITERATIVE MODE — lib exists, refine it
# ══════════════════════════════════════════════════════════════════════════

    max_iters=${MAX_ITERS:-5}
    conv_threshold=${CONV_THRESHOLD:-0.001}
    echo "=== Iterative mode (found $out) ===" >&2
    echo "    max_iters=$max_iters, conv_threshold=$conv_threshold" >&2

    # Parse existing lib's range
    lib_natt=$(grep -oP 'overshootNATT\s*=\s*\K[0-9]+' "$out")
    lib_nshape=$(grep -oP 'overshootNSHAPE\s*=\s*\K[0-9]+' "$out")
    lib_att_min_samp=$(grep -oP 'overshootATT_MIN_SAMP\s*=\s*\K[0-9.eE+-]+' "$out")
    lib_att_max_samp=$(grep -oP 'overshootATT_MAX_SAMP\s*=\s*\K[0-9.eE+-]+' "$out")

    # Convert to ms for the sweep DSP
    iter_att_min_ms=$(python3 -c "print(f'{$lib_att_min_samp / $sr * 1000:.10f}')")
    iter_att_max_ms=$(python3 -c "print(f'{$lib_att_max_samp / $sr * 1000:.10f}')")

    echo "    grid: ${lib_natt}×${lib_nshape}, att range: [$iter_att_min_ms, $iter_att_max_ms] ms" >&2

    for iter in $(seq 1 "$max_iters"); do
        echo "" >&2
        echo "── Iteration $iter/$max_iters ──" >&2

        local_dsp="$work/iter${iter}_sweep.dsp"
        gen_sweep_dsp "$local_dsp" \
            "$lib_natt" "$lib_nshape" "$iter_att_min_ms" "$iter_att_max_ms" "$out"
        run_sweep "iter$iter" "$local_dsp" \
            "$lib_natt" "$lib_nshape" "$out_dir"

        status=$(
            PLUGIN="$plugin_base" \
            NATT=$lib_natt NSHAPE=$lib_nshape SR=$sr \
            ATT_MIN_MS=$iter_att_min_ms ATT_MAX_MS=$iter_att_max_ms \
            MINS="${local_dsp%.dsp}_mins.txt" OUT="$out" \
            OLD_LIB="$out" CONV_THRESHOLD=$conv_threshold \
            write_lib
        )

        if [ "$status" = "converged" ]; then
            echo "" >&2
            echo "Converged after $iter iteration(s)." >&2
            break
        fi
    done

    if [ "$status" != "converged" ]; then
        echo "" >&2
        echo "Did not converge after $max_iters iterations (threshold=$conv_threshold)." >&2
        echo "Run again to continue refining, or raise CONV_THRESHOLD." >&2
    fi

else
# ══════════════════════════════════════════════════════════════════════════
#  INITIAL MODE — no lib, create from scratch (two-pass)
# ══════════════════════════════════════════════════════════════════════════

    floor_mult=${FLOOR_MULT:-2.0}
    echo "=== Initial mode (no existing lib) ===" >&2

    # Write a stub lib so the plugin (which imports it) can compile.
    # Identity values: mult=1.0, quantizers pass through.
    echo "  writing stub $out …" >&2
    cat >"$out" <<'STUB'
// overshoot_correction.lib — STUB (identity, no correction)
overshootNATT   = 2;
overshootNSHAPE = 2;
overshootATT_MIN_SAMP = 1.0;
overshootATT_MAX_SAMP = 2.0;
overshootLOG_RATIO    = 0.6931471805599453;
overshootSHAPE_MIN    = 0.0;
overshootSHAPE_MAX    = 1.0;
overshootCorrectionMult(att_samples, shapeSlider) = 1.0;
quantizeAttSamples(att_samples) = att_samples;
quantizeShape(s) = s;
STUB

    # ── Pass 1: full-range boundary detection ──
    echo "" >&2
    echo "── Pass 1: boundary detection ──" >&2
    nshape_p1=8
    gen_sweep_dsp "$work/p1_sweep.dsp" "$natt" "$nshape_p1" "$att_min_ms" "$att_max_ms"
    run_sweep "pass1" "$work/p1_sweep.dsp" "$natt" "$nshape_p1"

    boundary_att_ms=$(
    NATT=$natt NSHAPE_P1=$nshape_p1 SR=$sr FLOOR_MULT=$floor_mult \
    ATT_MIN_MS=$att_min_ms ATT_MAX_MS=$att_max_ms \
    MINS="$work/p1_sweep_mins.txt" \
    python3 - <<'PY_BOUNDARY'
import os, math, re, sys

natt = int(os.environ["NATT"])
nshape = int(os.environ["NSHAPE_P1"])
floor_mult = float(os.environ["FLOOR_MULT"])
att_min = float(os.environ["ATT_MIN_MS"])
att_max = float(os.environ["ATT_MAX_MS"])

_num = re.compile(r'^[+-]?\d')
raw = [float(l) for l in open(os.environ["MINS"])
       if l.strip() and _num.match(l.strip())]
assert len(raw) == natt * nshape
overshoot = [max(0.0, -v) for v in raw]

max_per_att = [max(overshoot[ia*nshape+js] for js in range(nshape))
               for ia in range(natt)]

tail = sorted(max_per_att[-(natt // 4):])
floor_val = tail[len(tail) // 2]
threshold = floor_val * floor_mult
print(f"  floor: {floor_val:.6f}, threshold (×{floor_mult}): {threshold:.6f}",
      file=sys.stderr)

boundary_idx = 0
for ia in range(natt):
    if max_per_att[ia] > threshold:
        boundary_idx = ia
boundary_idx = min(boundary_idx + 1, natt - 1)
att_norm = boundary_idx / max(1, natt - 1)
boundary_ms = att_min * (att_max / att_min) ** att_norm
sr = int(os.environ["SR"])
print(f"  boundary: attIdx={boundary_idx}, "
      f"{boundary_ms:.4f} ms ({boundary_ms/1000*sr:.1f} samp)", file=sys.stderr)
print(f"{boundary_ms:.10f}")
PY_BOUNDARY
    )

    echo "  boundary_att_ms=$boundary_att_ms" >&2

    # ── Pass 2: fine sweep in active region ──
    echo "" >&2
    echo "── Pass 2: fine sweep [$att_min_ms, $boundary_att_ms] ms ──" >&2
    gen_sweep_dsp "$work/p2_sweep.dsp" "$natt" "$nshape" "$att_min_ms" "$boundary_att_ms"
    run_sweep "pass2" "$work/p2_sweep.dsp" "$natt" "$nshape"

    # ── Write initial lib ──
    echo "" >&2
    echo "── Writing initial lib ──" >&2
    PLUGIN="$plugin_base" \
    NATT=$natt NSHAPE=$nshape SR=$sr \
    ATT_MIN_MS=$att_min_ms ATT_MAX_MS=$boundary_att_ms \
    MINS="$work/p2_sweep_mins.txt" OUT="$out" \
    CONV_THRESHOLD=999 \
    write_lib >/dev/null

    echo "" >&2
    echo "Initial lib written. Run again to iterate." >&2
fi
