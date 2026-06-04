#!/usr/bin/env bash
# Generate auc_poly.lib (the baked AUC level-compensation multiplier) from the
# plugin's own aucLevelMultFormula.
#
# Single source of truth: the ONLY Faust input is your plugin .dsp. The fit
# target — aucLevelMultFormula, and the shapeMap/cheapCurveBase it builds on —
# is read verbatim via library(); nothing about it is re-stated here or in
# Python. The only Python is the polynomial fit (pure stdlib, no numpy).
#
# For the REPORTED error we also build a high-accuracy reference: the same curve
# (ss.cheapCurveBase, ss.shapeMap, taken from the plugin) re-integrated at REF_M
# midpoint steps. Only the integration loop is restated, purely to vary M, so
# the error is measured against the true AUC rather than against the plugin's
# own coarse-aucM samples (which carry the integration bias and would otherwise
# make the fit look better than it is).
#
# Usage: regen_auc_poly.sh PLUGIN.dsp [DEGREE] [N] [OUT.lib]
#   PLUGIN.dsp  plugin defining shapeMap, cheapCurveBase, aucLevelMultFormula
#   DEGREE      polynomial degree                    (default 8)
#   N           fit-grid resolution                  (default 257)
#   OUT.lib     output library path  (default <plugin_dir>/auc_poly.lib)
#   REF_M       env var: midpoint steps for the high-accuracy error reference
#               (default 256). Only the integration M is varied; cheapCurveBase
#               and shapeMap still come from the plugin, so the reference stays
#               source-of-truth. 256 puts the reference's own integration error
#               (~3e-9) far below any realistic fit floor.
#
# Run this BEFORE adding import("auc_poly.lib") to the plugin: at generation
# time the plugin must compile without the lib (att/rel unscaled, no reference
# to aucLevelMult). Once the lib exists, import it and scale att/rel.
set -euo pipefail

plugin=${1:?usage: regen_auc_poly.sh PLUGIN.dsp [DEGREE] [N] [OUT.lib]}
degree=${2:-8}
N=${3:-257}
ref_m=${REF_M:-256}
plugin_dir=$(cd "$(dirname "$plugin")" && pwd)
plugin_base=$(basename "$plugin")
out=${4:-"$plugin_dir/auc_poly.lib"}

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# Driver: just a tap on the plugin's formula, sampled by ba.time. No math here.
cat >"$work/auc_eval.dsp" <<EOF
import("stdfaust.lib");
ss = library("$plugin_base");
aucN = $N;
process = ss.aucLevelMultFormula(float(ba.time)/float(aucN-1)) : min(1);
EOF

# Build (double), dump N samples (octave matrix), keep the numeric column.
(cd "$work" && faust2plot -double -I "$plugin_dir" auc_eval.dsp) >/dev/null 2>&1
"$work/auc_eval" -n "$N" >"$work/dump.m"
awk -F';' '/^[ ]*[-0-9]/{gsub(/[ ]/,"",$1); print $1}' "$work/dump.m" >"$work/vals.txt"

# High-accuracy reference for the error report: the SAME formula, re-derived at
# REF_M midpoint steps. cheapCurveBase and shapeMap are the plugin's (source of
# truth); only the integration loop is restated, so we can integrate finer than
# the plugin's hardcoded aucM and measure the polynomial against the true AUC.
cat >"$work/auc_ref.dsp" <<EOF
import("stdfaust.lib");
ss = library("$plugin_base");
refN = $N;
refM = $ref_m;
rBaseI(c) = sum(k, refM, ss.cheapCurveBase(c, (k+0.5)/refM))/refM;
rAtt(c)   = 1-(rBaseI(c)-ss.cheapCurveBase(c, 0))/(ss.cheapCurveBase(c, 1)-ss.cheapCurveBase(c, 0));
rMin      = rAtt(ss.shapeMap(1.0));
process   = rMin/rAtt(ss.shapeMap(float(ba.time)/float(refN-1))) : min(1);
EOF
(cd "$work" && faust2plot -double -I "$plugin_dir" auc_ref.dsp) >/dev/null 2>&1
"$work/auc_ref" -n "$N" >"$work/refdump.m"
awk -F';' '/^[ ]*[-0-9]/{gsub(/[ ]/,"",$1); print $1}' "$work/refdump.m" >"$work/ref.txt"

# Fit (pure stdlib: modified Gram-Schmidt QR in u=2s-1) and WRITE the lib.
PLUGIN="$plugin_base" SRC="$plugin_dir/$plugin_base" DEGREE=$degree N=$N REF_M=$ref_m VALS="$work/vals.txt" REF="$work/ref.txt" OUT="$out" python3 - <<'PY'
import os, math, re, sys

def polyfit_u(s, y, deg):                      # least squares, coeffs ascending in u=2s-1
    n, m = len(s), deg + 1
    V = [[(2*x - 1)**j for j in range(m)] for x in s]
    Q = [row[:] for row in V]
    R = [[0.0]*m for _ in range(m)]
    for j in range(m):
        for i in range(j):
            R[i][j] = sum(Q[k][i]*V[k][j] for k in range(n))
            for k in range(n):
                Q[k][j] -= R[i][j]*Q[k][i]
        R[j][j] = math.sqrt(sum(Q[k][j]**2 for k in range(n)))
        for k in range(n):
            Q[k][j] /= R[j][j]
    b = [sum(Q[k][i]*y[k] for k in range(n)) for i in range(m)]
    c = [0.0]*m
    for k in range(m-1, -1, -1):
        c[k] = (b[k] - sum(R[k][j]*c[j] for j in range(k+1, m))) / R[k][k]
    return c

clamp = lambda v: min(max(v, 0.0), 1.0)        # mirror aucLevelMult's :max(0):min(1)

y    = [float(l) for l in open(os.environ["VALS"]) if l.strip()]   # plugin formula @ its own aucM (fit target)
yref = [float(l) for l in open(os.environ["REF"])  if l.strip()]   # same formula @ high REF_M (the truth)
N, deg, refM = int(os.environ["N"]), int(os.environ["DEGREE"]), int(os.environ["REF_M"])
assert len(y)    == N, f"expected {N} fit samples from Faust, got {len(y)}"
assert len(yref) == N, f"expected {N} reference samples from Faust, got {len(yref)}"
s = [i/(N-1) for i in range(N)]
c = polyfit_u(s, y, deg)                        # fit is UNCHANGED: still to the plugin's own formula
pv = [clamp(sum(ci*(2*x-1)**i for i, ci in enumerate(c))) for x in s]   # polynomial as shipped

fit_resid = max(abs(p - yi) for p, yi in zip(pv, y))      # poly vs the samples it was fit to (old, misleading metric)
true_err  = max(abs(p - ri) for p, ri in zip(pv, yref))   # poly vs high-REF_M truth  <- honest headline
int_bias  = max(abs(yi - ri) for yi, ri in zip(y, yref))  # plugin's aucM samples vs truth (the masked error)

horner = f"a{deg}"
for i in range(deg-1, -1, -1):
    horner = f"a{i} + u*({horner})"

# Pull the exact fit-target definitions (verbatim, in dependency order) from the
# plugin so this lib records the full code + constants it was generated from.
# curveScale is intentionally absent: aucLevelMultFormula does not use it.
chain = ["shapeMap", "cheapCurveBase", "aucM", "aucBaseIntegral",
         "aucAttack", "aucMinAttack", "aucLevelMultFormula"]
try:
    src = open(os.environ["SRC"]).read().splitlines()
except OSError:
    src = []
prov, missing = [], []
for nm in chain:
    hit = next((ln.rstrip() for ln in src if re.match(rf"{re.escape(nm)}\b\s*[=(]", ln)), None)
    prov.append(hit) if hit else missing.append(nm)
if missing:
    print(f"  WARNING: could not find {', '.join(missing)} in {os.environ['PLUGIN']} -- "
          f"embedded provenance is incomplete.", file=sys.stderr)
refm_prefix = "" if refM == 256 else f"REF_M={refM} "

L = []
L.append("// auc_poly.lib — baked AUC level-compensation multiplier. GENERATED, do not edit.")
L.append(f"// degree-{deg} fit of aucLevelMultFormula (in u=2s-1), max abs err {true_err:.2e} "
         f"(~{20*math.log10(1+true_err):.4f} dB) vs the curve integrated at M={refM}.")
L.append(f"// breakdown: fit residual vs the plugin's own aucM samples {fit_resid:.2e}; "
         f"integration bias of those samples {int_bias:.2e}.")
L.append(f"// source of truth: {os.environ['PLUGIN']} (aucLevelMultFormula) via library().")
L.append(f"// regenerate: {refm_prefix}regen_auc_poly.sh {os.environ['PLUGIN']} {deg} {N}")
L.append("//")
L.append("// == generation provenance (commented; the live code lives in the plugin) ====")
L.append(f"// constants: degree={deg}, N={N} (uniform grid on s in [0,1]), REF_M={refM}.")
L.append(f"// method: sample aucLevelMultFormula at the N grid points, then least-squares fit")
L.append(f"//   a degree-{deg} polynomial in u=2s-1. The integral inside is the midpoint rule.")
L.append(f"// fit-target definitions, copied verbatim from {os.environ['PLUGIN']}:")
for ln in prov:
    L.append(f"//   {ln}")
L.append(f"// error reference: the same curve re-integrated at M={refM} midpoint steps, used")
L.append(f"//   only to measure the error reported above (it does not affect the fit).")
L.append("// ============================================================================")
L.append("aucLevelMult(s) = poly(s:max(0):min(1)):max(0):min(1)")
L.append("with {")
L.append(f"    poly(x) = {horner} with {{ u = 2*x - 1; }};")
for i, ci in enumerate(c):
    L.append(f"    a{i} = {ci!r};")
L.append("};")
open(os.environ["OUT"], "w").write("\n".join(L) + "\n")
print(f"wrote {os.environ['OUT']}  (degree {deg})")
print(f"  true error vs M={refM} reference : {true_err:.2e}  (~{20*math.log10(1+true_err):.4f} dB)  <- baked into header")
print(f"  fit residual vs plugin samples   : {fit_resid:.2e}")
print(f"  integration bias (plugin aucM)   : {int_bias:.2e}")
maxcoef = max(abs(ci) for ci in c)
# The reported max is taken over the N fit nodes; if N is too small for the
# degree, that max understates the true sup error (worst point falls between
# nodes). ~8 samples per residual lobe (deg+1 lobes) keeps the metric honest.
if N < 8*(deg+1):
    print(f"  WARNING: N={N} is low for degree {deg} -- the reported error may understate the "
          f"true max (grid undersamples the residual). Use N >= {8*(deg+1)}.")
if int_bias > fit_resid:
    print(f"  NOTE: integration bias dominates -> raise aucM in {os.environ['PLUGIN']}; "
          f"it caps accuracy before the degree-{deg} fit does.")
elif maxcoef > 1.0:
    # monomial-in-u conditioning gives out ~degree 16-18 in double: coeffs sit
    # ~0.5 then balloon, and past the knee a higher degree fits WORSE.
    print(f"  NOTE: fit residual dominates, but max|coef|={maxcoef:.1f} (>1) -- the monomial-in-u basis "
          f"is out of conditioning at degree {deg}; higher degree fits WORSE. Switch to a Chebyshev "
          f"basis to go lower. aucM is fine.")
else:
    print(f"  NOTE: fit residual dominates -> raise DEGREE to tighten. (N={N} is already converged; "
          f"more N won't help and can slightly raise the sup error. aucM is fine.)")
PY
