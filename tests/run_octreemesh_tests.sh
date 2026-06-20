#!/usr/bin/env bash
set -euo pipefail

g++ -std=c++17 -O2 -Wall -Wextra -pedantic captop.cpp -o captop
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

cat > "$TMP/single.octree" <<'EOF'
; Geometry file

{Nodes}
3 ; Dimension
8 ; Nodes count
; X1 X2 X3
0 0 0
1 0 0
1 1 0
0 1 0
0 0 1
1 0 1
1 1 1
0 1 1

{Mesh}
5 ; Element type (2=Triangle, 3=Quadrilateral, 4=Tetrahedra, 5=Hexahedra)
8 ; Nodes per element
1 ; Elements count
; Material Node1 Node2 ...
7 1 2 3 4 5 6 7 8
EOF

./captop validate "$TMP/single.octree" --grid strict-cube > "$TMP/validate.out" 2> "$TMP/validate.err"
rg -q 'Nodes               : 8' "$TMP/validate.out"
rg -q 'Hexahedra parsed    : 1' "$TMP/validate.out"
rg -q 'Parsing OctreeMesh nodes: 8' "$TMP/validate.err"
rg -q 'Finished OctreeMesh elements: 1' "$TMP/validate.err"

./captop convert "$TMP/single.octree" --filtration material --out "$TMP/out" > /dev/null 2> "$TMP/convert.err"
python3 - "$TMP/out" <<'PYRAW'
import pathlib, struct, sys
out=pathlib.Path(sys.argv[1])
vals=list(struct.unpack('<d', (out/'captop_cube_values_f64.raw').read_bytes()))
occ=list((out/'captop_occupied_u8.raw').read_bytes())
mat=list(struct.unpack('<q', (out/'captop_material_i64.raw').read_bytes()))
eid=list(struct.unpack('<q', (out/'captop_element_id_i64.raw').read_bytes()))
assert vals == [7.0], vals
assert occ == [1], occ
assert mat == [7], mat
assert eid == [1], eid
PYRAW
rg -q 'Parsing OctreeMesh elements: 1' "$TMP/convert.err"

cat > "$TMP/bad_extra_column.octree" <<'EOF'
{Nodes}
3
8
0 0 0
1 0 0
1 1 0
0 1 0
0 0 1
1 0 1
1 1 1
0 1 1
{Mesh}
5
8
1
7 1 2 3 4 5 6 7 8 9
EOF
if ./captop validate "$TMP/bad_extra_column.octree" > "$TMP/bad.out" 2> "$TMP/bad.err"; then
  echo 'expected extra OctreeMesh element column to fail' >&2
  exit 1
fi
rg -q 'must contain exactly material plus 8 node ids' "$TMP/bad.err"

cat > "$TMP/bad_type.octree" <<'EOF'
{Nodes}
3
0
{Mesh}
4
8
0
EOF
if ./captop validate "$TMP/bad_type.octree" > "$TMP/bad_type.out" 2> "$TMP/bad_type.err"; then
  echo 'expected unsupported OctreeMesh element type to fail' >&2
  exit 1
fi
rg -q 'only 5=Hexahedra is supported' "$TMP/bad_type.err"
