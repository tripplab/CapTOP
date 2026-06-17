# CapTOP

by trippm@tripplab.com on june 2026

CapTOP is a small C++ command-line tool for cubical meshes. Its current
implementation reads GiD ASCII `.msh` files, extracts supported 3D 8-node
hexahedral elements, validates that those elements can be interpreted as a
bitmap-style cubical complex, and converts valid meshes into dense indexed
cubical bitmap files for downstream topology analysis.

The executable name used throughout this README is `captop`.

## What CapTOP does

CapTOP currently provides two operational commands:

```bash
captop validate <input.msh> [options]
captop convert <input.msh> [options]
```

During validation, CapTOP:

1. Parses GiD ASCII mesh blocks.
2. Accepts mesh blocks with `dimension = 3`, `elemtype = Hexahedra`, and
   `nnode = 8`.
3. Reads node coordinates and hexahedral element connectivity.
4. Verifies that each hexahedron references 8 distinct, existing nodes.
5. Checks that every supported element is axis-aligned and occupies exactly one
   cell in the selected grid interpretation.
6. Detects duplicate occupied cells.
7. Reports mesh bounds, inferred grid dimensions, origin, spacing, occupied
   cell count, face adjacencies, boundary faces, and material/layer counts when
   material IDs are present.

The validator exits with:

- `0` when the mesh is valid.
- `1` for command-line, input, or parse errors.
- `2` when parsing succeeds but validation fails.

The converter reuses validation and exits with:

- `0` when conversion succeeds.
- `1` for command-line, input, parse, or I/O setup errors.
- `2` when validation fails.
- `3` when validation succeeds but bitmap conversion fails.

## Grid validation modes

CapTOP supports three grid modes:

- `strict-cube` (default, legacy alias `strict`): every accepted element must be
  an exact equal-edge cube and all cubes must lie on one uniform cubic lattice.
- `approximate-cube`: elements must be axis-aligned and nearly cubic according
  to `--cube-rel-tol`, while grid mapping uses snapped rectilinear coordinate
  axes.
- `rectilinear`: accepted elements may be axis-aligned rectangular voxels on
  rectilinear coordinate axes, but each element must still occupy one interval
  along each axis.

Use `strict-cube` for synthetic or debugging meshes that should be a conventional
equal-spacing voxel/cube model. Use `approximate-cube` for GiD/octree outputs
with small coordinate-rounding artifacts. Use `rectilinear` for topology-first
analysis of axis-aligned meshes whose spacing can vary by coordinate interval.

## Supported input format

CapTOP expects GiD-style ASCII mesh content with `mesh`, `coordinates`, and
`elements` sections. Supported element blocks must be 3D 8-node hexahedra.
Element records are interpreted as:

```text
<element_id> <node_1> <node_2> <node_3> <node_4> <node_5> <node_6> <node_7> <node_8> [material_id]
```

The optional tenth value is treated as a material or layer identifier. Extra
columns after that value are ignored with a warning.

Unsupported mesh blocks are reported as warnings when encountered. By default,
an element record inside an unsupported block is a parse error. Pass
`--ignore-non-hexa` to skip those element records instead.

## Requirements

- A C++17 compiler such as `g++` or `clang++`.
- Standard C++ library headers only; there are no third-party runtime
  dependencies in the current codebase.
- Optional: `git` for cloning the repository.
- Optional: `micromamba` if you want an isolated build environment.

## Clone the repository

```bash
git clone https://github.com/tripplab/CapTOP.git
cd CapTOP
```

## Build directly with a system compiler

From the repository root:

```bash
g++ -std=c++17 -O2 -Wall -Wextra -pedantic captop.cpp -o captop
```

Check that the binary runs:

```bash
./captop --version
./captop --help
```

## Create a micromamba build environment

If your system does not already have a suitable C++ compiler, create an isolated
Conda-style environment with `micromamba`:

```bash
micromamba create -n captop -c conda-forge cxx-compiler make git
micromamba activate captop
```

Then clone and build inside that environment:

```bash
git clone https://github.com/tripplab/CapTOP.git
cd CapTOP
g++ -std=c++17 -O2 -Wall -Wextra -pedantic captop.cpp -o captop
./captop --version
```

When you are done working on CapTOP, leave the environment with:

```bash
micromamba deactivate
```

If your shell has not yet been initialized for micromamba activation, run the
shell hook recommended by your micromamba installation first, for example:

```bash
eval "$(micromamba shell hook --shell bash)"
micromamba activate captop
```

## Usage

Show version information:

```bash
./captop --version
```

Show help:

```bash
./captop --help
```

Validate a GiD mesh with the default strict cubic lattice checks:

```bash
./captop validate mesh.msh
```

Validate with explicit strict mode:

```bash
./captop validate mesh.msh --grid strict-cube
```

Validate approximately cubic cells using a 2% relative spread tolerance:

```bash
./captop validate mesh.msh --grid approximate-cube --cube-rel-tol 0.02
```

Validate using rectilinear axis-aligned cells:

```bash
./captop validate mesh.msh --grid rectilinear
```

Adjust coordinate/snapping tolerance:

```bash
./captop validate mesh.msh --tol 1e-8
```

Skip unsupported non-hexahedral element blocks:

```bash
./captop validate mesh.msh --ignore-non-hexa
```

Limit the number of printed validation errors:

```bash
./captop validate mesh.msh --max-errors 10
```

Options can be combined:

```bash
./captop validate mesh.msh --grid approximate-cube --tol 1e-9 --cube-rel-tol 0.01 --ignore-non-hexa --max-errors 50
```

## Bitmap conversion

Stage 4 adds dense GiD-to-bitmap conversion without computing Betti numbers,
persistent homology, barcodes, or any GUDHI-dependent topology. Conversion first
validates the mesh; invalid meshes are never converted.

The dense array order is documented in the metadata and is used by every raw
file:

```text
linear_index(i,j,k) = i + nx * (j + ny * k)
```

A missing cell inside the validated bounding box is represented explicitly with
`occupied = 0`, `cube_value = +inf`, `material = -1`, and `element_id = -1`.

Basic occupancy conversion:

```bash
./captop convert mesh.msh --grid strict --filtration occupancy --out out_occ
```

Material-valued conversion and single-material extraction:

```bash
./captop convert mesh.msh --filtration material --out out_material
./captop convert mesh.msh --filtration material --selected-material 3 --out out_material3
```

Element scalar CSV conversion and binary thresholding:

```bash
./captop convert mesh.msh --filtration scalar-file --scalar-file values.csv --out out_scalar
./captop convert mesh.msh --filtration binary-threshold --scalar-file values.csv --threshold 0.5 --threshold-op ge --out out_threshold
```

The scalar CSV must have a header with `element_id` and `value` columns. Every
occupied element must have exactly one finite scalar value; duplicate and missing
entries are conversion errors.

The converter writes:

```text
outdir/
  captop_grid_metadata.json
  captop_cube_values_f64.raw
  captop_occupied_u8.raw
  captop_material_i64.raw
  captop_element_id_i64.raw
  captop_conversion_report.txt
```

Existing output files are protected by default. Use `--overwrite` to replace
them. Dense allocation is limited by `--max-memory-gb` and can be overridden with
`--force`.

## Example validation report

A valid mesh prints a report similar to:

```text
============================================================
CAPTOP validation report
============================================================
Software version      : 0.1.0-stage4
Input file            : mesh.msh
Grid mode             : strict-cube
Coordinate tolerance  : 1e-08
Status under selected mode: VALID

GiD mesh structure
  Mesh blocks         : 1
  Supported hex blocks: 1
  Unsupported blocks  : 0
  Nodes               : 8
  Hexahedra parsed    : 1

Coordinate bounds
  min                 : 0 0 0
  max                 : 1 1 1

Cubical grid
  Grid dimensions     : 1 x 1 x 1
  Origin              : 0 0 0
  Spacing             : 1 1 1
  Occupied cells      : 1
  Face adjacencies    : 0
  Boundary faces      : 6

Materials/layers
  1 : 1 elements

Ready for Stage 4 bitmap conversion under selected mode: yes
============================================================
```

## Minimal mesh example

Save the following as `one_cube.msh`:

```text
mesh "one_cube" dimension 3 elemtype Hexahedra nnode 8
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
end elements
```

Build, validate, and convert it:

```bash
g++ -std=c++17 -O2 -Wall -Wextra -pedantic captop.cpp -o captop
./captop validate one_cube.msh
./captop convert one_cube.msh --filtration occupancy --out one_cube_out
```

## Development notes

The current repository is intentionally compact:

- `captop.cpp` contains the command-line interface, GiD parser, geometry checks,
  grid mapping, dense bitmap conversion, diagnostics, and report generation.
- `tests/run_stage3_validation_tests.sh` exercises exact-cube, approximate-cube,
  rectilinear, duplicate-cell, non-axis-aligned, and missing-node cases.
- `tests/run_stage4_tests.sh` builds CapTOP, creates synthetic GiD meshes, runs
  conversion modes, and checks metadata plus raw arrays with `tests/check_raw_arrays.py`.
- `README.md` documents how to build and run the tool.
- `LICENSE` contains the MIT license.

The version output notes that GUDHI integration is not enabled in the current
stage. The validation report includes `Ready for Stage 4 bitmap conversion under
selected mode`; the conversion report includes `Ready for Stage 5 topology
computation: yes` after the bitmap files have been written.

## License

CapTOP is distributed under the MIT License. See [LICENSE](LICENSE) for details.
