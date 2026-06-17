#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
g++ -std=c++17 -O2 -Wall -Wextra -pedantic captop.cpp -o captop
./captop --help | grep -q -- "Convert options:"
./captop --help | grep -q -- "--filtration <policy>"
./captop convert --help | grep -q -- "--out <dir>"
./captop convert --help | grep -q -- "--threshold-op <op>"
./captop convert --help | grep -q -- "--cube-rel-tol"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
make_mesh(){ python3 - "$@" <<'PY'
from pathlib import Path
import sys
path=Path(sys.argv[1]); nx=int(sys.argv[2]); ny=int(sys.argv[3]); nz=int(sys.argv[4]); missing=set(sys.argv[5].split(';')) if sys.argv[5] else set(); mats=[int(x) for x in sys.argv[6].split(',')] if sys.argv[6] else []
nid={}; nodes=[]
def node(x,y,z):
 key=(x,y,z)
 if key not in nid: nid[key]=len(nid)+1; nodes.append(key)
 return nid[key]
elems=[]; eid=10; mi=0
for k in range(nz):
 for j in range(ny):
  for i in range(nx):
   if f'{i},{j},{k}' in missing: continue
   ids=[node(i,j,k),node(i+1,j,k),node(i+1,j+1,k),node(i,j+1,k),node(i,j,k+1),node(i+1,j,k+1),node(i+1,j+1,k+1),node(i,j+1,k+1)]
   mat = mats[mi] if mi < len(mats) else None; mi+=1
   elems.append((eid,ids,mat)); eid+=10
with path.open('w') as f:
 f.write('MESH "m" dimension 3 ElemType Hexahedra Nnode 8\nCoordinates\n')
 for idx,(x,y,z) in enumerate(nodes,1): f.write(f'{idx} {x} {y} {z}\n')
 f.write('End Coordinates\nElements\n')
 for eid,ids,mat in elems: f.write(' '.join(map(str,[eid]+ids+([] if mat is None else [mat])))+'\n')
 f.write('End Elements\n')
PY
}
make_mesh "$TMP/single.msh" 1 1 1 '' ''
./captop convert "$TMP/single.msh" --out "$TMP/out_single" >/dev/null
python3 tests/check_raw_arrays.py single "$TMP/out_single"
make_mesh "$TMP/block.msh" 2 2 2 '' ''
./captop convert "$TMP/block.msh" --out "$TMP/out_block" >/dev/null
python3 tests/check_raw_arrays.py block "$TMP/out_block"
make_mesh "$TMP/sparse.msh" 2 2 2 '1,1,1' ''
./captop convert "$TMP/sparse.msh" --out "$TMP/out_sparse" >/dev/null
python3 tests/check_raw_arrays.py sparse "$TMP/out_sparse"
make_mesh "$TMP/material.msh" 2 1 1 '' '5,9'
./captop convert "$TMP/material.msh" --filtration material --out "$TMP/out_mat" >/dev/null
python3 tests/check_raw_arrays.py material "$TMP/out_mat"
./captop convert "$TMP/material.msh" --filtration material --selected-material 9 --out "$TMP/out_mat9" >/dev/null
python3 tests/check_raw_arrays.py material9 "$TMP/out_mat9"
make_mesh "$TMP/scalar.msh" 2 1 1 '' ''
printf 'element_id,value\n10,1.5\n20,2.5\n' > "$TMP/values.csv"
./captop convert "$TMP/scalar.msh" --filtration scalar-file --scalar-file "$TMP/values.csv" --out "$TMP/out_scalar" >/dev/null
python3 tests/check_raw_arrays.py scalar "$TMP/out_scalar"
make_mesh "$TMP/threshold.msh" 3 1 1 '' ''
printf 'element_id,value\n10,0.1\n20,0.5\n30,0.9\n' > "$TMP/tvalues.csv"
./captop convert "$TMP/threshold.msh" --filtration binary-threshold --scalar-file "$TMP/tvalues.csv" --threshold 0.5 --threshold-op ge --out "$TMP/out_threshold" >/dev/null
python3 tests/check_raw_arrays.py threshold "$TMP/out_threshold"
printf 'element_id,value\n10,1.5\n' > "$TMP/missing.csv"
if ./captop convert "$TMP/scalar.msh" --filtration scalar-file --scalar-file "$TMP/missing.csv" --out "$TMP/out_missing" 2>"$TMP/err"; then exit 1; fi
grep -q 'missing scalar value for element 20' "$TMP/err"
printf 'element_id,value\n10,1.5\n10,2.0\n20,2.5\n' > "$TMP/dup.csv"
if ./captop convert "$TMP/scalar.msh" --filtration scalar-file --scalar-file "$TMP/dup.csv" --out "$TMP/out_dup" 2>"$TMP/err2"; then exit 1; fi
grep -q 'duplicate scalar value for element 10' "$TMP/err2"
if ./captop convert "$TMP/single.msh" --out "$TMP/out_single" 2>"$TMP/err3"; then exit 1; fi
grep -q 'output file already exists' "$TMP/err3"
./captop convert "$TMP/single.msh" --out "$TMP/out_single" --overwrite >/dev/null
echo 'Stage 4 tests passed'
