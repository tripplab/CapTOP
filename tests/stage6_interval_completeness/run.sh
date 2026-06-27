#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BIN="${CAPTOP_BIN:-$ROOT/build/captop}"
WORK="$ROOT/tests/stage6_interval_completeness/work"
mkdir -p "$WORK"
if ! "$BIN" --version | grep -q 'GUDHI support: enabled'; then echo "Stage 6 interval completeness tests require GUDHI; skipping."; exit 0; fi
run_case(){
  local case="$1" b0="$2" b1="$3" b2="$4"
  python3 "$ROOT/tests/stage6_interval_completeness/generate_mesh.py" "$case" "$WORK/$case.msh"
  rm -rf "$WORK/out_$case"
  "$BIN" persist "$WORK/$case.msh" --grid rectilinear --out "$WORK/out_$case" --overwrite --homology-dim 0,1,2 >/tmp/captop_stage6_${case}.txt
  grep -q "beta0(0) from intervals     : $b0" /tmp/captop_stage6_${case}.txt
  grep -q "beta1(0) from intervals     : $b1" /tmp/captop_stage6_${case}.txt
  grep -q "beta2(0) from intervals     : $b2" /tmp/captop_stage6_${case}.txt
  grep -q "Stage 5 consistency check   : PASS" /tmp/captop_stage6_${case}.txt
  test -f "$WORK/out_$case/diagram_dim0.csv" && test -f "$WORK/out_$case/diagram_dim1.csv" && test -f "$WORK/out_$case/diagram_dim2.csv"
  grep -q "^0,0,$b0," "$WORK/out_$case/betti_curve.csv"
}
run_case single 1 0 0
run_case two 2 0 0
run_case ring 1 1 0
run_case shell 1 0 1
