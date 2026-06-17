#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
mkdir -p build-stage5 tests/stage5/generated
cmake -S . -B build-stage5 -DCAPTOP_WITH_GUDHI=OFF >/dev/null
cmake --build build-stage5 >/dev/null
BIN=build-stage5/captop
make_mesh(){ tests/make_cubic_gid_mesh.py "tests/stage5/generated/$1.msh" --nx "$2" --ny "$3" --nz "$4" --occupied "$5" ${6:-}; }
run_case(){ name=$1; shift; out="tests/stage5/generated/${name}_out"; rm -rf "$out"; $BIN validate "tests/stage5/generated/$name.msh" "$@" >/dev/null; $BIN betti "tests/stage5/generated/$name.msh" --out "$out" --overwrite --quiet "$@"; }
make_mesh single 1 1 1 all
run_case single; tests/check_betti_summary.py tests/stage5/generated/single_out/captop_betti_summary.json topology.H0=1 topology.H1=0 topology.H2=0 topology.chi=1 topology.surface_faces=6
make_mesh block2 2 2 2 all
run_case block2; tests/check_betti_summary.py tests/stage5/generated/block2_out/captop_betti_summary.json topology.H0=1 topology.H1=0 topology.H2=0 topology.chi=1 topology.surface_faces=24
make_mesh sep 3 1 1 0:0:0,2:0:0
run_case sep; tests/check_betti_summary.py tests/stage5/generated/sep_out/captop_betti_summary.json topology.H0=2 topology.H1=0 topology.H2=0 topology.chi=2
make_mesh rect 3 2 1 all
run_case rect; tests/check_betti_summary.py tests/stage5/generated/rect_out/captop_betti_summary.json topology.H0=1 topology.H1=0 topology.H2=0 topology.chi=1
make_mesh ring 3 3 1 0:0:0,1:0:0,2:0:0,0:1:0,2:1:0,0:2:0,1:2:0,2:2:0
run_case ring; tests/check_betti_summary.py tests/stage5/generated/ring_out/captop_betti_summary.json topology.H0=1 topology.H1=1 topology.H2=0 topology.chi=0
make_mesh shell 3 3 3 0:0:0,1:0:0,2:0:0,0:1:0,1:1:0,2:1:0,0:2:0,1:2:0,2:2:0,0:0:1,1:0:1,2:0:1,0:1:1,2:1:1,0:2:1,1:2:1,2:2:1,0:0:2,1:0:2,2:0:2,0:1:2,1:1:2,2:1:2,0:2:2,1:2:2,2:2:2
run_case shell; tests/check_betti_summary.py tests/stage5/generated/shell_out/captop_betti_summary.json topology.H0=1 topology.H1=0 topology.H2=1 topology.chi=2 topology.surface_faces=60
make_mesh rlin 2 1 1 all '--spacing 1,2,5'
run_case rlin --grid rectilinear; tests/check_betti_summary.py tests/stage5/generated/rlin_out/captop_betti_summary.json topology.H0=1 topology.H1=0 topology.H2=0 topology.chi=1
echo "Stage 5 tests passed"
