# CapTOP

by trippm@tripplab.com on june 2026

CapTOP is a small C++ command-line validator for cubical meshes. Its current
implementation reads GiD ASCII `.msh` files, extracts supported 3D 8-node
hexahedral elements, and checks whether those elements can be interpreted as a
bitmap-style cubical complex for downstream topology analysis.

The executable name used throughout this README is `captop`.

## What CapTOP does

CapTOP currently provides one operational command:

```bash
captop validate <input.msh> [options]
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

## Grid validation modes

CapTOP supports two grid modes:

- `strict` (default): every accepted element must be a true cube with equal
  edge lengths, and all cubes must lie on one uniform cubic lattice.
- `rectilinear`: accepted elements may be axis-aligned rectangular voxels on
  rectilinear coordinate axes, but each element must still occupy one interval
  along each axis.

Use `strict` when your data should be a conventional equal-spacing voxel/cube
model. Use `rectilinear` for axis-aligned meshes whose spacing can vary by
coordinate interval.

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
./captop validate mesh.msh --grid strict
```

Validate using rectilinear axis-aligned cells:

```bash
./captop validate mesh.msh --grid rectilinear
```

Adjust floating-point tolerance:

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
./captop validate mesh.msh --grid rectilinear --tol 1e-9 --ignore-non-hexa --max-errors 50
```

## Example validation report

A valid mesh prints a report similar to:

```text
============================================================
CAPTOP validation report
============================================================
Software version      : 0.1.0-stage3
Input file            : mesh.msh
Grid mode             : strict
Tolerance             : 1e-08
Status                : VALID

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

Ready for Stage 4 bitmap conversion: yes
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

Build and validate it:

```bash
g++ -std=c++17 -O2 -Wall -Wextra -pedantic captop.cpp -o captop
./captop validate one_cube.msh
```

## Development notes

The current repository is intentionally compact:

- `captop.cpp` contains the command-line interface, GiD parser, geometry checks,
  grid mapping, and report generation.
- `README.md` documents how to build and run the tool.
- `LICENSE` contains the MIT license.

The version output notes that GUDHI integration is not enabled in the current
stage. The final report includes `Ready for Stage 4 bitmap conversion` to show
whether the parsed cells passed the compatibility checks expected before a later
bitmap/topology-analysis stage.

## License

CapTOP is distributed under the MIT License. See [LICENSE](LICENSE) for details.
