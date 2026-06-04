#!/usr/bin/env bash
# Generate auc_poly.lib (the baked AUC level-compensation multiplier) from the
# plugin's own aucLevelMultFormula.
#
# Single source of truth: the ONLY Faust input is your plugin .dsp. shapeMap,
# cheapCurveBase AND the AUC reduction (aucLevelMultFormula) are all read from
# it via library(); nothing is re-stated here or in Python. The only Python is
# the polynomial fit (pure stdlib, no numpy).
#
# Usage: regen_auc_poly.sh PLUGIN.dsp [DEGREE] [N] [OUT.lib]
#   PLUGIN.dsp  plugin defining shapeMap, cheapCurveBase, aucLevelMultFormula
#   DEGREE      polynomial degree                    (default 8)
#   N           fit-grid resolution                  (default 257)
#   OUT.lib     output library path  (default <plugin_dir>/auc_poly.lib)
#
# Run this BEFORE adding import("auc_poly.lib") to the plugin: at generation
# time the plugin must compile without the lib (att/rel unscaled, no reference
# to aucLevelMult). Once the lib exists, import it and scale att/rel.
set -euo pipefail

plugin=${1:?usage: regen_auc_poly.sh PLUGIN.dsp [DEGREE] [N] [OUT.lib]}
degree=${2:-8}
N=${3:-257}
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

# Fit (pure stdlib: modified Gram-Schmidt QR in u=2s-1) and WRITE the lib.
PLUGIN="$plugin_base" DEGREE=$degree N=$N VALS="$work/vals.txt" OUT="$out" python3 - <<'PY'
import os, math

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

y = [float(line) for line in open(os.environ["VALS"]) if line.strip()]
N, deg = int(os.environ["N"]), int(os.environ["DEGREE"])
assert len(y) == N, f"expected {N} samples from Faust, got {len(y)}"
s = [i/(N-1) for i in range(N)]
c = polyfit_u(s, y, deg)
err = max(abs(sum(ci*(2*x-1)**i for i, ci in enumerate(c)) - yi) for x, yi in zip(s, y))

horner = f"a{deg}"
for i in range(deg-1, -1, -1):
    horner = f"a{i} + u*({horner})"

L = []
L.append("// auc_poly.lib — baked AUC level-compensation multiplier. GENERATED, do not edit.")
L.append(f"// degree-{deg} fit of aucLevelMultFormula (in u=2s-1), max abs err {err:.2e} "
         f"(~{20*math.log10(1+err):.4f} dB).")
L.append(f"// source of truth: {os.environ['PLUGIN']} (aucLevelMultFormula) via library().")
L.append(f"// regenerate: regen_auc_poly.sh {os.environ['PLUGIN']} {deg} {N}")
L.append("aucLevelMult(s) = poly(s:max(0):min(1)):max(0):min(1)")
L.append("with {")
L.append(f"    poly(x) = {horner} with {{ u = 2*x - 1; }};")
for i, ci in enumerate(c):
    L.append(f"    a{i} = {ci!r};")
L.append("};")
open(os.environ["OUT"], "w").write("\n".join(L) + "\n")
print(f"wrote {os.environ['OUT']}  (degree {deg}, max abs err {err:.2e} ~ {20*math.log10(1+err):.4f} dB)")
PY
