#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
mkdir -p build-stage5-fallback tests/stage5_fallback_reporting/generated
cmake -S . -B build-stage5-fallback -DCAPTOP_WITH_GUDHI=OFF >/dev/null
cmake --build build-stage5-fallback >/dev/null
BIN=build-stage5-fallback/captop
make_mesh(){ tests/make_cubic_gid_mesh.py "tests/stage5_fallback_reporting/generated/$1.msh" --nx "$2" --ny "$3" --nz "$4" --occupied "$5"; }
assert_contains(){ local f=$1 s=$2; if ! grep -Fq "$s" "$f"; then echo "missing expected string: $s" >&2; exit 1; fi; }
assert_not_contains(){ local f=$1 s=$2; if grep -Fq "$s" "$f"; then echo "forbidden string present: $s" >&2; exit 1; fi; }
make_mesh single_cube 1 1 1 all
make_mesh tunnel_3x3x1 3 3 1 0:0:0,1:0:0,2:0:0,0:1:0,2:1:0,0:2:0,1:2:0,2:2:0
make_mesh hollow_shell_3x3x3 3 3 3 0:0:0,1:0:0,2:0:0,0:1:0,1:1:0,2:1:0,0:2:0,1:2:0,2:2:0,0:0:1,1:0:1,2:0:1,0:1:1,2:1:1,0:2:1,1:2:1,2:2:1,0:0:2,1:0:2,2:0:2,0:1:2,1:1:2,2:1:2,0:2:2,1:2:2,2:2:2
out=tests/stage5_fallback_reporting/generated/single_out; rm -rf "$out"
$BIN betti tests/stage5_fallback_reporting/generated/single_cube.msh --out "$out" --overwrite > "$out.terminal.txt"
for s in "GUDHI support       : disabled at compile time" "GUDHI used          : no" "Fallback used       : yes" "Result authority    : diagnostic fallback, not GUDHI-authenticated" "Topology inferred from fallback diagnostics" "Ready for Stage 6 persistent homology: no"; do assert_contains "$out.terminal.txt" "$s"; done
for s in "Topology from GUDHI at filtration threshold 0" "Ready for Stage 6 persistent homology: yes" "VALID AND GUDHI-ANALYZED"; do assert_not_contains "$out.terminal.txt" "$s"; done
out=tests/stage5_fallback_reporting/generated/tunnel_out; rm -rf "$out"
$BIN betti tests/stage5_fallback_reporting/generated/tunnel_3x3x1.msh --out "$out" --overwrite --quiet
assert_contains "$out/captop_betti_report.txt" "H1 tunnels                : 1  [inferred from H0_closed + H2 - chi]"
tests/check_betti_summary.py "$out/captop_betti_summary.json" topology.H0=1 topology.H1=1 topology.H2=0 topology.chi=0
python3 - <<'PY'
import json
with open('tests/stage5_fallback_reporting/generated/tunnel_out/captop_betti_summary.json') as f: data=json.load(f)
assert data['topology']['H1_source']=='inferred_from_closed_H0_complement_H2_and_cellular_euler'
PY
out=tests/stage5_fallback_reporting/generated/shell_out; rm -rf "$out"
$BIN betti tests/stage5_fallback_reporting/generated/hollow_shell_3x3x3.msh --out "$out" --overwrite --quiet --write-json --write-csv
python3 - <<'PY'
import json
with open('tests/stage5_fallback_reporting/generated/shell_out/captop_betti_summary.json') as f: data=json.load(f)
expected={
('topology_engine','gudhi_compiled'): False, ('topology_engine','gudhi_used'): False, ('topology_engine','gudhi_success'): False,
('topology_engine','fallback_used'): True, ('topology_engine','topology_source'): 'fallback_diagnostic', ('topology_engine','topology_authoritative'): False,
('topology','H0_source'): 'closed_cube_26_neighbor_union_find', ('topology','H1_source'): 'inferred_from_closed_H0_complement_H2_and_cellular_euler', ('topology','H2_source'): 'complement_flood_fill',
('topology','H3_source'): 'assumed_zero_fallback', ('topology','chi_source'): 'cubical_cell_count', ('readiness','stage6_persistent_homology'): False}
for (section,key), value in expected.items(): assert data[section][key] == value, (section, key, data[section][key], value)
PY
python3 - <<'PY'
import csv
from pathlib import Path
p=Path('tests/stage5_fallback_reporting/generated/shell_out/captop_betti_summary.csv')
row=next(csv.DictReader(p.open()))
required=['gudhi_compiled','gudhi_used','gudhi_success','fallback_used','topology_source','topology_authoritative','H0_source','H1_source','H2','H3','chi']
missing=[c for c in required if c not in row]
assert not missing, missing
expected={'gudhi_compiled':'false','gudhi_used':'false','gudhi_success':'false','fallback_used':'true','topology_source':'fallback_diagnostic','topology_authoritative':'false'}
for k,v in expected.items(): assert row[k]==v, (k,row[k],v)
PY
if $BIN betti tests/stage5_fallback_reporting/generated/single_cube.msh --require-gudhi > /tmp/captop_req.out 2>/tmp/captop_req.err; then echo "--require-gudhi unexpectedly succeeded" >&2; exit 1; fi
assert_contains /tmp/captop_req.err "CAPTOP was built without GUDHI support"
assert_contains /tmp/captop_req.err "rebuild with -DCAPTOP_WITH_GUDHI=ON"
$BIN --version > /tmp/captop_version.out
assert_contains /tmp/captop_version.out "GUDHI support: disabled"
assert_contains /tmp/captop_version.out "Stage 5 fallback diagnostics: enabled"
echo "Stage 5 fallback reporting tests passed"
