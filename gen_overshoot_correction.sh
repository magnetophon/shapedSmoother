#!/usr/bin/env bash
# Generate overshoot_correction.lib — a baked 2D correction surface for the
# shaped-smoother's undershoot on a 1→0 step.
#
# Single source of truth: every curve function (shapeMap, cheapCurveBase, the
# env loop, etc.) is imported from the plugin via library(). The ONLY things
# restated here are the grid/block structure and the env loop wiring (which
# must be inlined so that att/shape come from ba.time instead of sliders).
#
# The correction surface is parameterized by (att_samples, shapeSlider), NOT
# by milliseconds, because the overshoot is a numerical-integration error
# that depends on the number of discrete steps through the curve, not on
# wall-clock time.  The lib takes att_samples directly (the caller already
# has this as att*ma.SR:max(1)) and does log normalization internally.
#
# Approach:
#   1. Generate a Faust DSP that sweeps (att, shape) in blocks via ba.time.
#      Uses TARGET_SR (not ma.SR) so the step granularity matches the user's
#      target sample rate regardless of what faust2plot compiles at.
#      Each block: SETTLE samples at input=1 (let the smoother converge to 1),
#      then MEASURE samples at input=0 (step down; track the running minimum).
#   2. Compile once with faust2plot -double; run once to dump all samples.
#   3. Extract the per-block minimum (last sample of each block) with awk.
#   4. Fit a 2D polynomial in (u,v) = (2·attNorm−1, 2·shapeNorm−1)
#      where attNorm = log(attSamp/minSamp)/log(maxSamp/minSamp) (log-scaled).
#   5. Write overshoot_correction.lib with a Horner-form evaluator.
#
# Usage:
#   gen_overshoot_correction.sh PLUGIN.dsp [DEG_A] [DEG_S] [NATT] [NSHAPE] [SR] [OUT.lib]
#     PLUGIN.dsp  the shapedSmoother core (must export shapeMap, cheapCurveBase,
#                 curveScale, derivativeBaseRelease, maxDerivativeBaseAttack,
#                 peakPhaseRelease, inverseDerivative{Top,Bottom}{Attack,Release})
#     DEG_A       polynomial degree in attack-time axis   (default 6)
#     DEG_S       polynomial degree in shape axis          (default 6)
#     NATT        grid points along attack-time axis       (default 32)
#     NSHAPE      grid points along shape axis             (default 32)
#     SR          target sample rate                       (default 48000)
#     OUT.lib     output library path  (default <plugin_dir>/overshoot_correction.lib)
#   FLOOR_CUTOFF  env var: att_samples above which the fit target is zeroed
#                 (default 100). The measured overshoot at long attacks is a
#                 ~0.01 floor that the polynomial shouldn't waste DOF on.
set -euo pipefail

plugin=${1:?usage: gen_overshoot_correction.sh PLUGIN.dsp [DEG_A] [DEG_S] [NATT] [NSHAPE] [SR] [OUT.lib]}
deg_a=${2:-6}
deg_s=${3:-6}
natt=${4:-32}
nshape=${5:-32}
sr=${6:-48000}
plugin_dir=$(cd "$(dirname "$plugin")" && pwd)
plugin_base=$(basename "$plugin")
out=${7:-"$plugin_dir/overshoot_correction.lib"}

# Attack-time slider range (ms) — must match the plugin
att_min_ms="0.046"
att_max_ms="50.0"

# Block geometry scales with the target sample rate.
# max att_samples = maxHold(0.05 s) × SR.  SETTLE and MEASURE each get 2×
# that so the sliding-min window, delay line, and release all clear.
max_att_samp=$((sr / 20))            # 0.05 * sr, integer
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
#
# TARGET_SR replaces ma.SR everywhere so the env loop's step granularity
# matches the user's target rate, regardless of faust2plot's internal rate.
cat >"$work/overshoot_sweep.dsp" <<FAUST
import("stdfaust.lib");
ss = library("$plugin_base");

// ── Configuration (substituted by the shell) ─────────────────────────────
NATT       = $natt;
NSHAPE     = $nshape;
SETTLE     = $settle;
MEASURE    = $measure;
BLOCK      = SETTLE + MEASURE;
TARGET_SR  = $sr.0;        // the sample rate we are measuring for

ATT_MIN_MS = $att_min_ms;
ATT_MAX_MS = $att_max_ms;
REL_MS     = 0.5;          // fast, gentle release for settling between blocks

MAX_HOLD   = 0.05;         // seconds — must match the plugin
MAX_HOLD_SAMP = int(MAX_HOLD * TARGET_SR);

// ── Block / grid indexing ────────────────────────────────────────────────
blockIdx   = int(ba.time) / BLOCK;
posInBlock = int(ba.time) % BLOCK;
attIdx     = min(NATT  - 1, blockIdx / NSHAPE);
shapeIdx   = min(NSHAPE - 1, blockIdx - attIdx * NSHAPE);

// Log-spaced attack (ms), linear shape (0–1)
currentAttMs = ATT_MIN_MS
             * pow(ATT_MAX_MS / ATT_MIN_MS,
                   float(attIdx) / float(max(1, NATT - 1)));
currentShape = float(shapeIdx) / float(max(1, NSHAPE - 1));

// Step signal: 1 during SETTLE, drops to 0 at the MEASURE boundary
stepInput = 1 - (posInBlock >= SETTLE);

// ── Parametric shaped smoother ───────────────────────────────────────────
// Inlines the env loop from the plugin, but wires att/shape from the grid
// and uses TARGET_SR instead of ma.SR.  Curve functions come from ss.*.
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

// ── Run & track ──────────────────────────────────────────────────────────
smootherOut = stepInput : paramSmoother(currentAttMs, currentShape);

// Running minimum during MEASURE; resets at the SETTLE→MEASURE boundary.
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
# grep keeps only numeric data lines (skip %----ChunkBoundary---- etc.);
# awk picks every BLOCK-th one (the last sample of each block, where the
# running-min tracker holds the final MEASURE-phase minimum).
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

# ── 4 & 5. Fit 2D polynomial & write the lib ─────────────────────────────
PLUGIN="$plugin_base" SRC="$plugin_dir/$plugin_base" \
DEG_A=$deg_a DEG_S=$deg_s NATT=$natt NSHAPE=$nshape SR=$sr \
MINS="$work/mins.txt" OUT="$out" \
ATT_MIN_MS=$att_min_ms ATT_MAX_MS=$att_max_ms \
FLOOR_CUTOFF=${FLOOR_CUTOFF:-100} \
python3 - <<'PY'
import os, math, sys, re

# ── read config ───────────────────────────────────────────────────────────
natt    = int(os.environ["NATT"])
nshape  = int(os.environ["NSHAPE"])
deg_a   = int(os.environ["DEG_A"])
deg_s   = int(os.environ["DEG_S"])
sr      = int(os.environ["SR"])
att_min_ms = float(os.environ["ATT_MIN_MS"])
att_max_ms = float(os.environ["ATT_MAX_MS"])

# Sample-domain bounds for the log normalization baked into the lib
att_min_samp = att_min_ms / 1000.0 * sr   # e.g. 2.208  @ 48 kHz
att_max_samp = att_max_ms / 1000.0 * sr   # e.g. 2400.0 @ 48 kHz
log_ratio    = math.log(att_max_samp / att_min_samp)
floor_cutoff = float(os.environ.get("FLOOR_CUTOFF", "100"))

# ── read measurements ────────────────────────────────────────────────────
_num = re.compile(r'^[+-]?\d')
raw = [float(l) for l in open(os.environ["MINS"])
       if l.strip() and _num.match(l.strip())]
assert len(raw) == natt * nshape, \
    f"expected {natt*nshape} minima, got {len(raw)}"

# overshoot magnitude: how far below 0 the output went (≥ 0)
overshoot = [max(0.0, -v) for v in raw]

# grid coordinates, both in [0, 1] (log-spaced attack, linear shape)
att_norm   = [i / max(1, natt  - 1) for i in range(natt)]
shape_norm = [j / max(1, nshape - 1) for j in range(nshape)]

# att_samples for each attIdx (needed for the cutoff)
att_samp_per_idx = [att_min_samp * (att_max_samp / att_min_samp) ** an
                    for an in att_norm]

S, Y, Y_raw = [], [], []   # coords, fit-target, raw measurement
for ia in range(natt):
    for js in range(nshape):
        S.append((att_norm[ia], shape_norm[js]))
        raw_val = overshoot[ia * nshape + js]
        Y_raw.append(raw_val)
        # Zero the fit target above the cutoff — the polynomial shouldn't
        # spend degrees of freedom on the ~0.01 floor there.
        Y.append(raw_val if att_samp_per_idx[ia] <= floor_cutoff else 0.0)

max_overshoot = max(Y_raw)
n_active = sum(1 for s in att_samp_per_idx if s <= floor_cutoff) * nshape
print(f"  max overshoot measured: {max_overshoot:.6e}", file=sys.stderr)
print(f"  floor cutoff: {floor_cutoff:.0f} samples "
      f"({n_active} active / {len(S)} total points)", file=sys.stderr)
if max_overshoot < 1e-15:
    print("  WARNING: no measurable overshoot — the correction surface will "
          "be zero everywhere.  The lib is still written (all-zero poly).",
          file=sys.stderr)

# ── 2D least-squares fit ─────────────────────────────────────────────────
# Tensor-product basis  u^i · v^j,  u = 2·attNorm−1,  v = 2·shapeNorm−1.
n_terms = (deg_a + 1) * (deg_s + 1)
n_pts   = len(S)

def basis(a_n, s_n):
    u = 2.0 * a_n - 1.0
    v = 2.0 * s_n - 1.0
    row = []
    for i in range(deg_a + 1):
        ui = u ** i
        for j in range(deg_s + 1):
            row.append(ui * (v ** j))
    return row

V = [basis(a, s) for a, s in S]

# Modified Gram-Schmidt QR
m = n_terms
Q = [row[:] for row in V]
R = [[0.0] * m for _ in range(m)]
for j in range(m):
    for i in range(j):
        R[i][j] = sum(Q[k][i] * V[k][j] for k in range(n_pts))
        for k in range(n_pts):
            Q[k][j] -= R[i][j] * Q[k][i]
    R[j][j] = math.sqrt(sum(Q[k][j] ** 2 for k in range(n_pts)))
    if R[j][j] < 1e-30:
        R[j][j] = 1e-30
    for k in range(n_pts):
        Q[k][j] /= R[j][j]

b = [sum(Q[k][i] * Y[k] for k in range(n_pts)) for i in range(m)]
c = [0.0] * m
for k in range(m - 1, -1, -1):
    c[k] = (b[k] - sum(R[k][j] * c[j] for j in range(k + 1, m))) / R[k][k]

# Evaluate fit, compute residual
clamp01 = lambda x: min(max(x, 0.0), 1.0)

def poly_eval(a_n, s_n):
    u = 2.0 * a_n - 1.0
    v = 2.0 * s_n - 1.0
    val, idx = 0.0, 0
    for i in range(deg_a + 1):
        ui = u ** i
        for j in range(deg_s + 1):
            val += c[idx] * ui * (v ** j)
            idx += 1
    return clamp01(val)

pv = [poly_eval(a, s) for a, s in S]
# Error in the active region only (att_samples <= floor_cutoff)
active_err = [abs(p - y) for p, y, (an, _) in zip(pv, Y, S)
              if att_min_samp * (att_max_samp / att_min_samp) ** an <= floor_cutoff]
max_err  = max(active_err) if active_err else 0.0
rms_err  = math.sqrt(sum(e*e for e in active_err) / max(1, len(active_err)))
max_coef = max(abs(ci) for ci in c)

print(f"  fit residual (active, ≤{floor_cutoff:.0f} samp): "
      f"max {max_err:.2e}  rms {rms_err:.2e}", file=sys.stderr)
print(f"  max |coef|:   {max_coef:.2e}", file=sys.stderr)
if max_coef > 100:
    print(f"  WARNING: large coefficients — the monomial basis may be poorly "
          f"conditioned at degree ({deg_a},{deg_s}). Consider lowering the "
          f"degrees.", file=sys.stderr)

# ── Build the Horner-form Faust expression ────────────────────────────────
def horner_v(i):
    """Horner form in v for the i-th row of coefficients."""
    expr = f"a{i}_{deg_s}"
    for j in range(deg_s - 1, -1, -1):
        expr = f"a{i}_{j} + v*({expr})"
    return expr

horner_u = f"row{deg_a}"
for i in range(deg_a - 1, -1, -1):
    horner_u = f"row{i} + u*({horner_u})"

# ── Write the lib ────────────────────────────────────────────────────────
L = []
L.append( "// overshoot_correction.lib — baked 2D overshoot correction surface.")
L.append( "// GENERATED by gen_overshoot_correction.sh, do not hand-edit.")
L.append(f"// degree ({deg_a},{deg_s}) tensor-product fit in (u,v)")
L.append(f"//   u = 2·log(attSamp/{att_min_samp:.4f})/log({att_max_samp:.1f}/{att_min_samp:.4f}) − 1")
L.append(f"//   v = 2·shapeSlider − 1")
L.append(f"// training SR: {sr} Hz  (att_samples range [{att_min_samp:.4f}, {att_max_samp:.1f}])")
L.append(f"// floor cutoff: {floor_cutoff:.0f} samples (overshoot zeroed above this for the fit)")
L.append(f"// max abs residual (≤{floor_cutoff:.0f} samp): {max_err:.2e};  max measured overshoot: {max_overshoot:.6e}")
L.append(f"// grid: {natt}×{nshape} (att × shape)")
L.append(f"// source of truth: {os.environ['PLUGIN']} (env loop, curve functions) via library().")
L.append(f"// regenerate: FLOOR_CUTOFF={floor_cutoff:.0f} gen_overshoot_correction.sh"
         f" {os.environ['PLUGIN']} {deg_a} {deg_s} {natt} {nshape} {sr}")
L.append( "//")
L.append( "// Usage:  import(\"overshoot_correction.lib\");")
L.append( "//         correction = overshootCorrection(att_samples, shapeSlider);")
L.append( "//   where att_samples = att * ma.SR : max(1)  (number of attack steps)")
L.append( "//         shapeSlider ∈ [0,1]  (the raw attack-shape slider value)")
L.append( "//   Returns the overshoot magnitude (≥ 0): how far below 0 the")
L.append( "//   envelope goes on a unit step with those parameters.")
L.append( "// ============================================================================")
L.append( "")
L.append( "overshootCorrection(att_samples, shapeSlider) =")
L.append( "    poly(att_samples, shapeSlider : max(0) : min(1)) : max(0) : min(1)")
L.append( "with {")
L.append(f"    ATT_MIN_SAMP = {att_min_samp!r};")
L.append(f"    ATT_MAX_SAMP = {att_max_samp!r};")
L.append(f"    LOG_RATIO    = {log_ratio!r};")
L.append( "")
L.append( "    poly(a_samp, s) = " + horner_u + " with {")
L.append( "        // log-normalize att_samples to [−1, 1], clamped to training range")
L.append( "        u = 2 * log(max(ATT_MIN_SAMP, min(ATT_MAX_SAMP, a_samp))")
L.append( "                    / ATT_MIN_SAMP)")
L.append( "              / LOG_RATIO - 1;")
L.append( "        v = 2*s - 1;")
for i in range(deg_a + 1):
    L.append(f"        row{i} = {horner_v(i)};")
L.append( "    };")
for i in range(deg_a + 1):
    for j in range(deg_s + 1):
        idx = i * (deg_s + 1) + j
        L.append(f"    a{i}_{j} = {c[idx]!r};")
L.append( "};")

open(os.environ["OUT"], "w").write("\n".join(L) + "\n")
print(f"wrote {os.environ['OUT']}  (degree ({deg_a},{deg_s})  {n_terms} terms)",
      file=sys.stderr)

# ── Dump the raw grid for inspection ──────────────────────────────────────
grid_path = os.environ["OUT"].rsplit(".", 1)[0] + "_grid.tsv"
with open(grid_path, "w") as gf:
    gf.write("attIdx\tshapeIdx\tattMs\tattSamples\tshapeSlider"
             "\tovershoot_raw\tfit_target\tpolyFit\n")
    for ia in range(natt):
        att_ms = att_min_ms * (att_max_ms / att_min_ms) ** (
            ia / max(1, natt - 1))
        att_samp = att_ms / 1000.0 * sr
        for js in range(nshape):
            shape_sl = js / max(1, nshape - 1)
            k = ia * nshape + js
            gf.write(f"{ia}\t{js}\t{att_ms:.4f}\t{att_samp:.2f}"
                     f"\t{shape_sl:.4f}"
                     f"\t{Y_raw[k]:.8e}\t{Y[k]:.8e}\t{pv[k]:.8e}\n")
print(f"wrote {grid_path}  (raw data + fit for inspection)", file=sys.stderr)
PY
