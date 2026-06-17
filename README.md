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
- `cmake` for the supported Stage 5 build.
- GUDHI development headers for `captop betti` when building with
  `CAPTOP_WITH_GUDHI=ON`.
- Optional: `git` for cloning the repository.
- Optional: `micromamba` if you want an isolated build environment with all
  compiler, CMake, and GUDHI dependencies.

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

If your system does not already have a suitable C++ compiler, CMake, and GUDHI
installation, create an isolated Conda-style environment with `micromamba` from
the `conda-forge` channel:

```bash
micromamba create -n captop -c conda-forge \
  cxx-compiler \
  cmake \
  make \
  ninja \
  git \
  gudhi
micromamba activate captop
```

If your shell has not yet been initialized for micromamba activation, run the
shell hook recommended by your micromamba installation first, for example:

```bash
eval "$(micromamba shell hook --shell bash)"
micromamba activate captop
```

Then clone and build inside that environment with Stage 5 GUDHI support enabled:

```bash
git clone https://github.com/tripplab/CapTOP.git
cd CapTOP
cmake -S . -B build \
  -DCAPTOP_WITH_GUDHI=ON \
  -DGUDHI_INCLUDE_DIR="$CONDA_PREFIX/include"
cmake --build build
./build/captop --version
./build/captop --help
```

The `gudhi` Conda package installs the required headers under
`$CONDA_PREFIX/include`, which is why the CMake command above passes that path as
`GUDHI_INCLUDE_DIR`. If you only need to work on parser, validation, or bitmap
conversion code and do not need the GUDHI-backed `betti` command, you can build
without GUDHI instead:

```bash
cmake -S . -B build-no-gudhi -DCAPTOP_WITH_GUDHI=OFF
cmake --build build-no-gudhi
./build-no-gudhi/captop --version
```

When you are done working on CapTOP, leave the environment with:

```bash
micromamba deactivate
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

## Stage 5 Betti topology

Stage 5 adds plain occupied-domain topology descriptors for validated cubical GiD
meshes:

```bash
captop betti <input.msh> [options]
```

The `betti` command parses the mesh, runs the same Stage 3 grid validator used
by `validate`, builds the Stage 4 dense occupancy bitmap in memory with the
stable index order

```text
linear_index(i,j,k) = i + nx * (j + ny * k)
```

and computes topology for the binary occupied voxel domain. CAPTOP Stage 5
uses the topology of the closed occupied cubical complex as the primary topology
convention. Under this convention, occupied voxels whose closed cubes touch by a
face, edge, or vertex may belong to the same connected component; the primary
connectivity diagnostic is therefore 26-neighbor closed-cube connectivity.
CAPTOP also reports 6-neighbor face-adjacency components as a mesh-quality
diagnostic. This diagnostic may be larger than GUDHI `H0` when separate
face-connected components touch at edges or vertices. Missing voxels inside the
bounding box are assigned `+inf` in the cubical filtration, while occupied
top-dimensional cubes are assigned `0.0`. When CAPTOP is built with GUDHI
support, Betti numbers are queried at threshold `from = 0.0`, `to = 0.0`, not
from the final all-included rectangular filtration.

Reported descriptors include:

- `H0` connected components under the primary closed-cube topology.
- `H1` independent tunnels.
- `H2` enclosed cavities.
- Euler characteristic `chi = H0 - H1 + H2`.
- Occupied, missing, and total bounding-box voxel counts.
- Bounding-box dimensions.
- Surface voxel count and exposed surface face count.
- Closed-cube `H0_closed` (26-neighbor) and face-adjacency `H0_face`
  (6-neighbor) union-find diagnostics.
- Contact audit counts for face, edge, and vertex contacts.
- Raw inter-component contact counts for occupied voxel pairs whose
  6-neighbor face-adjacency components differ.
- A deterministic merger audit explaining reductions from `H0_face` to
  `H0_closed` when an independent subset of edge/vertex contacts merges face
  components. A mesh can have many raw contacts, for example 30 raw
  inter-component contacts, but only 6 merger records if many contacts join
  components already merged by earlier contacts.
- Cubical-cell Euler characteristic, GUDHI Euler characteristic, and cross-check
  status.
- Complement flood-fill `H2` diagnostic.

Example:

```bash
captop betti mesh.msh --grid strict --field 2 --out captop_betti_out --overwrite
```

Useful Stage 5 options:

```text
--field <prime>          coefficient field, default 2
--out <directory>        output directory, default captop_betti_out
--overwrite              allow replacing existing output files
--max-memory-gb <value>  dense-grid memory limit, default 2.0
--force                  continue above the memory limit
--write-diagram          reserve raw interval outputs for debugging
--write-contact-audit    request contact audit CSV files
--write-json             request JSON summary output
--write-csv              request CSV summary output
--quiet                  suppress the terminal report
```

By default, successful Betti runs write:

```text
captop_betti_summary.json
captop_betti_summary.csv
captop_betti_report.txt
captop_contact_audit_summary.json
```

Detailed contact CSV files are written when `--write-contact-audit` is passed,
when `H0_face != H0_closed`, when GUDHI `H0` disagrees with `H0_closed`, or when
contact-audit consistency checks fail:

```text
captop_contact_mergers.csv   independent merger records
captop_contact_contacts.csv  all raw inter-component contacts
```

Glossary:

- `face-adjacency component`: a component under 6-neighbor voxel connectivity.
- `closed-cube component`: a component under closed-cube 26-neighbor
  connectivity.
- `raw contact`: an edge or vertex contact between two distinct face-adjacency
  components.
- `merger record`: a raw contact selected by the face-component union-find that
  actually reduces the component count.
- `component reduction`: `H0_face - H0_closed`, which must match the number of
  independent merger records in a valid contact audit.

Exit statuses for `betti` are:

- `0` success.
- `1` CLI, parse, build, GUDHI, or I/O error.
- `2` validation error.
- `3` conversion/grid construction error.
- `4` topology cross-check failure.

## CMake build

CAPTOP now includes a CMake build. Stage 5 GUDHI support is requested by
default:

```bash
cmake -S . -B build -DCAPTOP_WITH_GUDHI=ON -DGUDHI_INCLUDE_DIR=/path/to/gudhi/include
cmake --build build
./build/captop --version
```

If GUDHI headers are unavailable, configure with `-DCAPTOP_WITH_GUDHI=OFF` for
validator/converter development builds. In that configuration, `captop betti`
still runs the independent diagnostics, but the JSON summary marks GUDHI as not
enabled and the report includes a warning.

## Stage 5 tests

The Stage 5 test runner builds CAPTOP, generates synthetic GiD meshes, runs
`validate` and `betti`, and checks JSON topology values exactly:

```bash
tests/run_stage5_tests.sh
tests/run_stage5_contact_audit_tests.sh
```

The synthetic cases include a single cube, solid `2 x 2 x 2` block, two
separated cubes, asymmetric `3 x 2 x 1` block, one-voxel-thick tunnel ring,
closed hollow shell, and rectilinear unequal-spacing block.

## Stage 5 topology sources

CAPTOP Stage 5 records explicit provenance for every reported topology value so
fallback diagnostics are not confused with GUDHI-authenticated Betti numbers.

With GUDHI enabled (`-DCAPTOP_WITH_GUDHI=ON`), the headline `H0`, `H1`, `H2`,
and `H3` values are computed from GUDHI persistent Betti numbers at filtration
threshold 0. GUDHI `H0` is expected to match `H0_closed`, the 26-neighbor
closed-cube union-find count. `H0_face`, the 6-neighbor face-adjacency count, is
a diagnostic only and is not a fatal cross-check when it differs from GUDHI.
Complement flood-fill and cubical-cell Euler calculations remain available as
independent cross-checks.

If `H0_face` differs from `H0_closed`, CAPTOP writes a contact audit showing
which edge or vertex contacts merge face components under the closed-cube
topology. `captop_contact_mergers.csv` contains one row for each successful
face-component merger, while `captop_contact_contacts.csv` contains all
inter-face-component contact records.

Without GUDHI (`-DCAPTOP_WITH_GUDHI=OFF`), `captop betti` can still produce
fallback diagnostics:

- `H0` from closed-cube 26-neighbor union-find.
- `H0_face` from 6-neighbor face-adjacency union-find as a diagnostic.
- `H2` from complement flood-fill.
- `chi` from cubical cell counts.
- `H1` inferred from `H0_closed + H2 - chi`.
- `H3` assumed zero for the finite voxel-subset fallback path.

Fallback results are useful for diagnostics and sanity checks, but they are not a
substitute for the required GUDHI-backed Stage 5 deliverable. Reports generated
without GUDHI must not be described as GUDHI Betti results. Use
`captop betti <input.msh> --require-gudhi` to fail rather than accept
fallback-only analysis when authoritative GUDHI topology is required.
