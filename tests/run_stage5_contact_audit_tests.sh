#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
mkdir -p build-stage5 tests/stage5_contact_audit/generated
cmake -S . -B build-stage5 -DCAPTOP_WITH_GUDHI=OFF >/dev/null
cmake --build build-stage5 >/dev/null
BIN=build-stage5/captop
make_mesh(){ tests/make_cubic_gid_mesh.py "tests/stage5_contact_audit/generated/$1.msh" --nx "$2" --ny "$3" --nz "$4" --occupied "$5"; }
run_case(){ name=$1; out="tests/stage5_contact_audit/generated/${name}_out"; rm -rf "$out"; $BIN betti "tests/stage5_contact_audit/generated/$name.msh" --out "$out" --overwrite --quiet; }
check_json(){ python3 - "$@" <<'PY'
import json, pathlib, sys
path=pathlib.Path(sys.argv[1]); data=json.loads(path.read_text())
for spec in sys.argv[2:]:
    key, val = spec.split('=',1)
    cur=data
    for part in key.split('.'): cur=cur[part]
    exp = val if not (val.lstrip('-').isdigit()) else int(val)
    if cur != exp:
        raise SystemExit(f"{key}: expected {exp!r}, got {cur!r}")
PY
}
check_rows(){ python3 - "$@" <<'PY'
import csv, pathlib, sys
rows=list(csv.DictReader(open(sys.argv[1], newline='')))
if len(rows)!=int(sys.argv[2]): raise SystemExit(f"{sys.argv[1]} rows expected {sys.argv[2]} got {len(rows)}")
if len(sys.argv)>3 and rows and rows[0]['contact_type']!=sys.argv[3]: raise SystemExit(f"contact_type expected {sys.argv[3]} got {rows[0]['contact_type']}")
PY
}
make_mesh face2 2 1 1 0:0:0,1:0:0
run_case face2
check_json tests/stage5_contact_audit/generated/face2_out/captop_betti_summary.json topology.H0=1 contact_audit.H0_face=1 contact_audit.H0_closed=1 contact_audit.component_reduction=0 contact_audit.face_component_merger_count=0
check_rows tests/stage5_contact_audit/generated/face2_out/captop_contact_mergers.csv 0
make_mesh edge2 2 2 1 0:0:0,1:1:0
run_case edge2
check_json tests/stage5_contact_audit/generated/edge2_out/captop_betti_summary.json topology.H0=1 contact_audit.H0_face=2 contact_audit.H0_closed=1 contact_audit.component_reduction=1 contact_audit.face_component_merger_count=1 contact_audit.inter_face_component_edge_contact_count=1
check_rows tests/stage5_contact_audit/generated/edge2_out/captop_contact_mergers.csv 1 edge
make_mesh vertex2 2 2 2 0:0:0,1:1:1
run_case vertex2
check_json tests/stage5_contact_audit/generated/vertex2_out/captop_betti_summary.json topology.H0=1 contact_audit.H0_face=2 contact_audit.H0_closed=1 contact_audit.component_reduction=1 contact_audit.face_component_merger_count=1 contact_audit.inter_face_component_vertex_contact_count=1
check_rows tests/stage5_contact_audit/generated/vertex2_out/captop_contact_mergers.csv 1 vertex
make_mesh edge_chain 3 3 1 0:0:0,1:1:0,2:2:0
run_case edge_chain
check_json tests/stage5_contact_audit/generated/edge_chain_out/captop_betti_summary.json topology.H0=1 contact_audit.H0_face=3 contact_audit.H0_closed=1 contact_audit.component_reduction=2 contact_audit.face_component_merger_count=2
check_rows tests/stage5_contact_audit/generated/edge_chain_out/captop_contact_mergers.csv 2
make_mesh separated 3 1 1 0:0:0,2:0:0
run_case separated
check_json tests/stage5_contact_audit/generated/separated_out/captop_betti_summary.json topology.H0=2 contact_audit.H0_face=2 contact_audit.H0_closed=2 contact_audit.component_reduction=0 contact_audit.face_component_merger_count=0
check_rows tests/stage5_contact_audit/generated/separated_out/captop_contact_mergers.csv 0
make_mesh ring 3 3 1 0:0:0,1:0:0,2:0:0,0:1:0,2:1:0,0:2:0,1:2:0,2:2:0
run_case ring
check_json tests/stage5_contact_audit/generated/ring_out/captop_betti_summary.json topology.H0=1 topology.H1=1 topology.H2=0 contact_audit.H0_face=1 contact_audit.H0_closed=1 contact_audit.component_reduction=0
echo "Stage 5 contact audit tests passed"
