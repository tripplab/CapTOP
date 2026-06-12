#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT
BIN="${TMP_DIR}/captop"

g++ -std=c++17 -O2 -Wall -Wextra -pedantic "${ROOT_DIR}/captop.cpp" -o "${BIN}"

write_box_mesh() {
  local path="$1" dx="$2" dy="$3" dz="$4"
  cat > "${path}" <<EOF
mesh "box" dimension 3 elemtype Hexahedra nnode 8
coordinates
1 0 0 0
2 ${dx} 0 0
3 ${dx} ${dy} 0
4 0 ${dy} 0
5 0 0 ${dz}
6 ${dx} 0 ${dz}
7 ${dx} ${dy} ${dz}
8 0 ${dy} ${dz}
end coordinates
elements
1 1 2 3 4 5 6 7 8 1
end elements
EOF
}

write_box_mesh "${TMP_DIR}/exact_cube.msh" 1.0 1.0 1.0
write_box_mesh "${TMP_DIR}/rounded_cube.msh" 0.999 1.000 1.001
write_box_mesh "${TMP_DIR}/rectangular_voxel.msh" 1.0 2.0 1.0

cat > "${TMP_DIR}/non_axis_aligned.msh" <<'EOF'
mesh "skew" dimension 3 elemtype Hexahedra nnode 8
coordinates
1 0 0 0
2 1 0 0
3 1 1 0
4 0.2 1 0
5 0 0 1
6 1 0 1
7 1 1 1
8 0.2 1 1
end coordinates
elements
1 1 2 3 4 5 6 7 8 1
end elements
EOF

cat > "${TMP_DIR}/duplicate_cell.msh" <<'EOF'
mesh "duplicate" dimension 3 elemtype Hexahedra nnode 8
coordinates
1 0 0 0
2 1 0 0
3 1 1 0
4 0 1 0
5 0 0 1
6 1 0 1
7 1 1 1
8 0 1 1
end coordinates
elements
1 1 2 3 4 5 6 7 8 1
2 1 2 3 4 5 6 7 8 1
end elements
EOF

cat > "${TMP_DIR}/missing_node.msh" <<'EOF'
mesh "missing" dimension 3 elemtype Hexahedra nnode 8
coordinates
1 0 0 0
2 1 0 0
3 1 1 0
4 0 1 0
5 0 0 1
6 1 0 1
7 1 1 1
end coordinates
elements
1 1 2 3 4 5 6 7 8 1
end elements
EOF

expect_status() {
  local expected="$1"
  shift
  set +e
  "$@" >/dev/null 2>&1
  local status=$?
  set -e
  if [[ "${status}" -ne "${expected}" ]]; then
    echo "Expected exit ${expected}, got ${status}: $*" >&2
    exit 1
  fi
}

# Exact cube mesh: valid in all modes.
expect_status 0 "${BIN}" validate "${TMP_DIR}/exact_cube.msh" --grid strict-cube
expect_status 0 "${BIN}" validate "${TMP_DIR}/exact_cube.msh" --grid approximate-cube
expect_status 0 "${BIN}" validate "${TMP_DIR}/exact_cube.msh" --grid rectilinear

# Rounded cube mesh: strict invalid, approximate and rectilinear valid.
expect_status 2 "${BIN}" validate "${TMP_DIR}/rounded_cube.msh" --grid strict-cube
expect_status 0 "${BIN}" validate "${TMP_DIR}/rounded_cube.msh" --grid approximate-cube --cube-rel-tol 0.01
expect_status 0 "${BIN}" validate "${TMP_DIR}/rounded_cube.msh" --grid rectilinear

# Rectangular voxel mesh: only rectilinear valid.
expect_status 2 "${BIN}" validate "${TMP_DIR}/rectangular_voxel.msh" --grid strict-cube
expect_status 2 "${BIN}" validate "${TMP_DIR}/rectangular_voxel.msh" --grid approximate-cube
expect_status 0 "${BIN}" validate "${TMP_DIR}/rectangular_voxel.msh" --grid rectilinear

# Geometric/grid failures are invalid in all modes.
for mesh in non_axis_aligned duplicate_cell missing_node; do
  expect_status 2 "${BIN}" validate "${TMP_DIR}/${mesh}.msh" --grid strict-cube
  expect_status 2 "${BIN}" validate "${TMP_DIR}/${mesh}.msh" --grid approximate-cube
  expect_status 2 "${BIN}" validate "${TMP_DIR}/${mesh}.msh" --grid rectilinear
done

# CLI help should expose the revised modes and cube relative tolerance.
"${BIN}" validate --help | grep -q -- "--grid strict-cube"
"${BIN}" validate --help | grep -q -- "--grid approximate-cube"
"${BIN}" validate --help | grep -q -- "--cube-rel-tol"

echo "stage3 validation regression tests passed"
