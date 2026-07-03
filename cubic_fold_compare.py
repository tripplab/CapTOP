#!/usr/bin/env python3
"""
cubic_fold_compare.py
by trippm  (cubical-complex counterpart to alpha_fold_compare_v02.py)

Compare CUBICAL persistence diagrams produced by CapTOP (>= 0.1.0-stage7) for
icosahedral capsid folds, the way alpha_fold_compare_v02.py compares the
alpha-complex diagrams of the atomic substructure.

WHAT THIS CONSUMES (the CapTOP producer/consumer contract)
----------------------------------------------------------
Each fold is one CapTOP `persist` output DIRECTORY containing at least:
    captop_diagram_manifest.json   - units, mode, comparability_key, per-dim files
    diagram_dim0.csv / _dim1 / _dim2 - birth,death,...,persistence,alive_at_zero
and, when CapTOP was run with `--write-grid`, additionally:
    captop_grid_metadata.json      - nx/ny/nz, spacing, origin, index_order
    captop_occupied_u8.raw         - dense occupancy (uint8)
    captop_cube_values_f64.raw     - the exact physical field GUDHI saw (float64 LE)

The grid files are what make this module self-validating and what will let a
later stage compute a +/-1-voxel selector-stability noise floor IN-PROCESS,
without re-running CapTOP.

WHY A SELF-CHECK IS STAGE 0 (and lives in the module, not a throwaway script)
-----------------------------------------------------------------------------
Every later stage rests on one atomic operation: turn a cube-value FIELD into a
set of per-dimension persistence diagrams via a GUDHI cubical complex
(`field_to_diagrams`). The self-check is exactly that operation applied to a
fold's OWN field, asserting it reproduces the diagram CapTOP already wrote
(Betti-at-0 vector and the finite-bar lifetime multiset). If it reproduces the
authoritative diagram, the index/axis order and units are proven, and the same
conversion can be trusted when applied to a PERTURBED field (the noise floor) or
when recomputing for any other reason. So the self-check is not scaffolding; it
is the acceptance test that locks down the one function the whole module reuses.

STAGE MAP (this file grows additively; only Stage 0 is active today)
--------------------------------------------------------------------
  Stage 0  [ACTIVE]  load_fold, field_to_diagrams, self_check
                     -> prove the grid contract on a single fold (e.g. h=12).
  Stage 1  [stub]    perturbation_floor : +/-1-voxel dilation/erosion -> SDT ->
                     field_to_diagrams -> bottleneck to the unperturbed diagram.
                     The cubical analogue of alpha_fold_compare control 2
                     (selector stability). Reuses field_to_diagrams verbatim.
  Stage 2  [stub]    compare_folds : per-dim bottleneck + 1-Wasserstein between
                     two folds, on the UNSIGNED-LINEAR axis, gated by matching
                     comparability_key AND h within a relative tolerance (--h-rtol,
                     default 1e-2), with optional common-h renormalization.
  Stage 3  [stub]    verdict : per-pair floor (max of the two folds' Stage-1
                     floors) -> "REAL fold difference" iff a dimension clears it.
  Stage 4  [stub]    capsid_vs_mesh_table : explicitly NON-metric structural
                     comparison (top-feature multiplicity bands + the
                     H2-count-vs-h divergence), since alpha (atoms, exact radii)
                     and cubical (voxelized) do not share a metric space and the
                     raw/denoised counts are resolution-dependent by ~1/h^2.

DESIGN DECISIONS ALREADY SETTLED (recorded so later stages honor them)
----------------------------------------------------------------------
  * fold-vs-fold REQUIRES h to match within a RELATIVE tolerance (--h-rtol,
    default 1e-2), a hard error past that: finite-bar counts range over orders of
    magnitude with resolution (19 -> 50246 for 1CWP F2_0 across h=12..1), so
    genuinely cross-h distances measure discretization. The tolerance (not exact
    equality) admits the ~4th-significant-figure h difference between two same-
    resolution folds voxelized independently; a 2x octree-level change is ~100%,
    far above any sane tolerance, so real resolution mismatches still hard-fail.
    Optional --renorm-h rescales all folds to a shared h to erase the residual
    sub-tolerance skew (off by default; below the noise floor).
  * default fold-vs-fold axis is UNSIGNED-LINEAR (Angstrom). The signed-square
    (Angstrom^2) axis is reserved for the capsid-vs-mesh shared-axis table; its
    huge dynamic range lets one dominant feature swamp the bottleneck.
  * the metric pre-filter is a LIGHT persistence floor (anti-swarm / numerics),
    not a cavity-isolating cut: no fixed or relative threshold recovers a stable
    cavity count, so we do not pretend one does.
  * comparability_key (policy|units|mode|sdt-sign/square) gates metric comparison;
    h is checked SEPARATELY on top of it (the key intentionally omits h so the
    capsid-vs-mesh table may compare across resolutions).

Self-contained by choice (no import from the alpha modules yet); shared helpers
will be factored out only once two stages actually need the same code.

Requires: numpy, gudhi (>= 3.5, with CubicalComplex). gudhi only needed for the
topology operations (self-check, and later floor); pure loading/reporting do not.
"""

import argparse
import json
import os
import sys
from dataclasses import dataclass, field
from typing import Optional

import numpy as np


# ===========================================================================
# Tunable defaults
# ===========================================================================
# Relative tolerance for the fold-vs-fold spacing (h) match. Two folds that are
# the SAME nominal resolution but voxelized independently (e.g. different fold
# orientations meshed separately by OctreeMesh) have h values that differ at the
# ~4th significant figure -- export/orientation quantization, NOT a resolution
# difference. Adjacent octree levels differ by 2x (100%), so a 1e-2 gate cleanly
# admits quantization noise (~1e-3 observed) while still rejecting a true level
# change with ~50x headroom. Override on the CLI with --h-rtol.
DEFAULT_H_RTOL = 1e-2


# ===========================================================================
# Loaded representation of one fold (one CapTOP persist directory)
# ===========================================================================
@dataclass
class FoldData:
    """Everything loaded from one CapTOP `persist` output directory.

    This is the object every stage passes around. Topology arrays (occupancy,
    field) are loaded lazily and only when grid files are present, so a fold with
    only diagrams (no --write-grid) is still usable for the pure-diagram stages
    (compare/verdict) but not for the self-check or the perturbation floor.
    """
    path: str
    label: str                              # short name for reports (dir basename)
    manifest: dict                          # captop_diagram_manifest.json
    diagrams: dict                          # dim -> (N,2) float array [birth, death] (finite only)
    essentials: dict                        # dim -> int count of infinite/essential bars
    diagram_raw: dict                       # dim -> structured ndarray of full CSV rows
    # --- comparability fields lifted from the manifest for quick gating ---
    comparability_key: Optional[str] = None
    units: Optional[str] = None
    h: Optional[float] = None               # uniform spacing (Angstrom)
    # --- grid (present only if CapTOP was run with --write-grid) ---
    has_grid: bool = False
    grid_meta: Optional[dict] = None        # captop_grid_metadata.json
    nx: Optional[int] = None
    ny: Optional[int] = None
    nz: Optional[int] = None
    _occupancy: Optional[np.ndarray] = field(default=None, repr=False)
    _field: Optional[np.ndarray] = field(default=None, repr=False)

    # ---- lazy loaders for the heavy arrays -------------------------------
    def occupancy(self) -> np.ndarray:
        """Dense occupancy as a flat uint8 array in CapTOP index order
        (i + nx*(j + ny*k)). Raises if --write-grid was not used."""
        if not self.has_grid:
            raise RuntimeError(
                f"fold '{self.label}' has no grid files; re-run CapTOP persist "
                f"with --write-grid to enable self-check / perturbation floor")
        if self._occupancy is None:
            p = os.path.join(self.path, "captop_occupied_u8.raw")
            occ = np.fromfile(p, dtype=np.uint8)
            self._expect_size(occ, p)
            self._occupancy = occ
        return self._occupancy

    def cube_field(self) -> np.ndarray:
        """The exact physical filtration field CapTOP fed to GUDHI, flat float64
        in CapTOP index order. Raises if --write-grid was not used."""
        if not self.has_grid:
            raise RuntimeError(
                f"fold '{self.label}' has no grid files; re-run CapTOP persist "
                f"with --write-grid to enable self-check / perturbation floor")
        if self._field is None:
            p = os.path.join(self.path, "captop_cube_values_f64.raw")
            fld = np.fromfile(p, dtype="<f8")        # little-endian per contract
            self._expect_size(fld, p)
            self._field = fld
        return self._field

    def _expect_size(self, arr, p):
        want = self.nx * self.ny * self.nz
        if arr.size != want:
            raise ValueError(
                f"{p}: read {arr.size} values but grid is "
                f"{self.nx}x{self.ny}x{self.nz} = {want}; "
                f"endianness or dtype mismatch?")


# ===========================================================================
# Stage 0a: loading a fold from a CapTOP persist directory
# ===========================================================================
def _read_diagram_csv(path):
    """Read a CapTOP diagram_dimN.csv. Returns (finite_pairs (N,2), n_essential,
    raw_structured). 'death' may be the literal 'inf'."""
    # columns: birth,death,birth_physical,death_physical,death_type,persistence,alive_at_zero
    births_phys, deaths_phys, dtypes, raw_rows = [], [], [], []
    with open(path) as fh:
        header = fh.readline().strip().split(",")
        col = {name: i for i, name in enumerate(header)}
        need = ("birth_physical", "death_physical", "death_type")
        for n in need:
            if n not in col:
                raise ValueError(f"{path}: missing column '{n}' (got {header})")
        for line in fh:
            line = line.strip()
            if not line:
                continue
            t = line.split(",")
            raw_rows.append(t)
            dt = t[col["death_type"]]
            dtypes.append(dt)
            bp = float(t[col["birth_physical"]])
            dp_s = t[col["death_physical"]]
            if dt == "finite" and dp_s not in ("inf", "+inf", "Inf"):
                births_phys.append(bp)
                deaths_phys.append(float(dp_s))
    finite = (np.column_stack([births_phys, deaths_phys]).astype(float)
              if births_phys else np.empty((0, 2)))
    n_essential = sum(1 for d in dtypes if d != "finite")
    return finite, n_essential, raw_rows


def load_fold(path):
    """Load one CapTOP persist directory into a FoldData. Pure I/O; no topology.

    Reads the manifest (required), the per-dimension diagram CSVs it lists, and -
    if the manifest advertises grid_files (i.e. CapTOP was run with --write-grid) -
    records the grid metadata so the heavy arrays can be lazily loaded later.
    """
    if not os.path.isdir(path):
        raise NotADirectoryError(f"not a directory: {path}")
    man_p = os.path.join(path, "captop_diagram_manifest.json")
    if not os.path.isfile(man_p):
        raise FileNotFoundError(
            f"{path}: no captop_diagram_manifest.json (is this a CapTOP persist "
            f"output dir from >= stage7?)")
    manifest = json.load(open(man_p))

    filt = manifest.get("filtration", {})
    units = filt.get("units")
    comparability_key = manifest.get("comparability_key")
    grid_blk = filt.get("grid") or manifest.get("grid") or {}
    h = grid_blk.get("spacing_h")
    if h is None:
        # spacing may live as a 3-vector in the grid_files block
        gf = manifest.get("grid_files") or {}
        gg = gf.get("grid") if isinstance(gf, dict) else None
        if gg and "spacing" in gg:
            h = gg["spacing"][0]

    # per-dimension diagrams
    diagrams, essentials, diagram_raw = {}, {}, {}
    dims_listed = manifest.get("dimensions", [])
    if dims_listed:
        items = [(d["dim"], d["file"]) for d in dims_listed]
    else:                                    # fallback: probe dim0..2
        items = [(d, f"diagram_dim{d}.csv") for d in (0, 1, 2)
                 if os.path.isfile(os.path.join(path, f"diagram_dim{d}.csv"))]
    for d, fname in items:
        fpath = os.path.join(path, fname)
        if not os.path.isfile(fpath):
            raise FileNotFoundError(f"{path}: manifest lists {fname} but it is absent")
        fin, ness, raw = _read_diagram_csv(fpath)
        diagrams[d] = fin
        essentials[d] = ness
        diagram_raw[d] = raw

    fd = FoldData(
        path=path, label=os.path.basename(os.path.normpath(path)),
        manifest=manifest, diagrams=diagrams, essentials=essentials,
        diagram_raw=diagram_raw, comparability_key=comparability_key,
        units=units, h=h)

    # grid files (optional)
    gridfiles = manifest.get("grid_files")
    if gridfiles and gridfiles.get("present"):
        gmeta_p = os.path.join(path, "captop_grid_metadata.json")
        if os.path.isfile(gmeta_p):
            gmeta = json.load(open(gmeta_p))
            g = gmeta.get("grid", {})
            fd.has_grid = True
            fd.grid_meta = gmeta
            fd.nx, fd.ny, fd.nz = g.get("nx"), g.get("ny"), g.get("nz")
            if fd.h is None:
                sp = g.get("spacing")
                if sp:
                    fd.h = sp[0]
        else:
            print(f"[warn] {path}: manifest advertises grid_files but "
                  f"captop_grid_metadata.json is missing; grid disabled",
                  file=sys.stderr)
    return fd


# ===========================================================================
# Stage 0b: the atomic operation - FIELD -> per-dimension diagrams via GUDHI
#           This is the one function self-check, perturbation-floor, and any
#           recompute all share. Locking its index/axis order here (proven by
#           the self-check) is what makes every later stage trustworthy.
# ===========================================================================
def field_to_diagrams(flat_field, nx, ny, nz, max_dim=2):
    """Build a GUDHI cubical complex from a FLAT cube-value field (CapTOP index
    order i + nx*(j + ny*k)) and return per-dimension finite diagrams + Betti@0.

    CapTOP constructs `Bitmap_cubical_complex(dimensions={nx,ny,nz}, field, true)`
    where `field` is the same flat array in i-fastest order. GUDHI-Python's
    CubicalComplex(dimensions=[nx,ny,nz], top_dimensional_cells=flat) consumes the
    flat array in the SAME convention (first dimension varies fastest), so passing
    CapTOP's flat array directly is the matching call. The self-check is what
    proves this; if it ever fails with H1/H2 looking swapped, the fix is a
    transpose HERE (and only here), not anywhere downstream.

    Returns
    -------
    dict with:
      'diagrams'  : {d: (N,2) finite [birth, death]}
      'essentials': {d: int}                       # infinite-death bars
      'betti0'    : {d: int}                        # persistent Betti at threshold 0
    All values are in the field's own units (physical, per the manifest).
    """
    import gudhi
    flat = np.ascontiguousarray(flat_field, dtype=float)
    if flat.size != nx * ny * nz:
        raise ValueError(f"field size {flat.size} != {nx}*{ny}*{nz}")
    cc = gudhi.CubicalComplex(dimensions=[nx, ny, nz],
                              top_dimensional_cells=flat)
    cc.compute_persistence(homology_coeff_field=2, min_persistence=0.0)

    diagrams = {d: [] for d in range(max_dim + 1)}
    essentials = {d: 0 for d in range(max_dim + 1)}
    for dim, (b, dth) in cc.persistence():
        if dim > max_dim:
            continue
        if dth == float("inf"):
            essentials[dim] += 1
        else:
            diagrams[dim].append((b, dth))
    diagrams = {d: (np.asarray(v, float).reshape(-1, 2) if v else np.empty((0, 2)))
                for d, v in diagrams.items()}

    # Betti at threshold 0 (the occupied/interface slice). persistent_betti_numbers
    # counts classes born <= 0 and dying > 0.
    betti0 = {}
    pb = cc.persistent_betti_numbers(0.0, 0.0)
    for d in range(max_dim + 1):
        betti0[d] = int(pb[d]) if d < len(pb) else 0
    return {"diagrams": diagrams, "essentials": essentials, "betti0": betti0}


# ===========================================================================
# Stage 0c: self-check - does the fold's own field reproduce its own diagram?
# ===========================================================================
def _lifetime_multiset(pairs, rel_tol=1e-4):
    """Collapse a set of birth/death pairs into (rounded-lifetime -> count) bands
    using a RELATIVE tolerance, so CapTOP and GUDHI-Python agreeing only to within
    a few ULPs still match. Mirrors the alpha-side relative-tolerance fix rather
    than exact float equality."""
    if len(pairs) == 0:
        return {}
    lives = np.sort((pairs[:, 1] - pairs[:, 0]).astype(float))[::-1]
    bands = {}
    i, n = 0, len(lives)
    while i < n:
        j = i + 1
        while j < n and abs(lives[j] - lives[i]) <= rel_tol * max(1.0, abs(lives[i])):
            j += 1
        key = round(float(np.mean(lives[i:j])), 3)
        bands[key] = bands.get(key, 0) + (j - i)
        i = j
    return bands


def _diagram_bands(pairs, rel_tol=1e-4):
    """Collapse (birth,death) pairs into bands keyed by ROUNDED (birth,death) using
    a relative tolerance, returning {(b_round,d_round): count}. Used to compare two
    finite diagrams as multisets tolerant to last-ULP differences."""
    if len(pairs) == 0:
        return {}
    P = np.asarray(pairs, float)
    order = np.lexsort((P[:, 1], P[:, 0]))
    P = P[order]
    bands = {}
    i, n = 0, len(P)
    while i < n:
        j = i + 1
        while (j < n
               and abs(P[j, 0] - P[i, 0]) <= rel_tol * max(1.0, abs(P[i, 0]))
               and abs(P[j, 1] - P[i, 1]) <= rel_tol * max(1.0, abs(P[i, 1]))):
            j += 1
        bb = round(float(np.mean(P[i:j, 0])), 3)
        dd = round(float(np.mean(P[i:j, 1])), 3)
        bands[(bb, dd)] = bands.get((bb, dd), 0) + (j - i)
        i = j
    return bands


def _bands_match(a, b, rel_tol):
    """Two (birth,death)->count band dicts match iff same total and every band in a
    has a counterpart in b within rel_tol on BOTH coordinates with equal count."""
    if sum(a.values()) != sum(b.values()):
        return False, f"bar counts differ: {sum(a.values())} vs {sum(b.values())}"
    bkeys = list(b.keys())
    used = [False] * len(bkeys)
    for (ba, da), ca in a.items():
        hit = False
        for idx, (bb, db) in enumerate(bkeys):
            if used[idx]:
                continue
            if (abs(ba - bb) <= rel_tol * max(1.0, abs(ba))
                    and abs(da - db) <= rel_tol * max(1.0, abs(da))
                    and b[(bb, db)] == ca):
                used[idx] = True
                hit = True
                break
        if not hit:
            return False, f"band (birth={ba}, death={da})x{ca} has no match"
    return True, ""


def self_check(fold, max_dim=2, rel_tol=1e-4, verbose=True):
    """Stage 0 gate (strict, complete-diagram equality). Reproduce the fold's
    diagram from its own grid field via field_to_diagrams and require it to equal,
    per dimension, the COMPLETE diagram CapTOP wrote:

      * finite bars match as a (birth_physical, death_physical) multiset
        (relative-tolerance banded), AND
      * the count of essential (infinite-death) bars matches.

    This is achievable only because CapTOP (>= stage7.2) records the TRUE finite
    deaths of threshold-zero-alive classes instead of synthetic (0, inf) tokens;
    against older outputs (schema version 1) the finite multisets will not match
    and the check correctly fails. Betti@0 agreement is also asserted as the
    primary proof the index/axis order and units are correct.
    """
    report = {"label": fold.label, "ok": False, "checks": []}

    def record(name, ok, detail=""):
        report["checks"].append({"name": name, "ok": ok, "detail": detail})
        if verbose:
            print(f"  [{'PASS' if ok else 'FAIL'}] {name}"
                  + (f"  {detail}" if detail else ""))

    if not fold.has_grid:
        record("grid present", False,
               "no grid files (re-run CapTOP persist with --write-grid)")
        return report

    # Warn (don't fail yet) if the producer predates faithful deaths.
    sch = fold.manifest.get("schema", {})
    if not sch.get("deaths_faithful", sch.get("version", 1) >= 2):
        record("producer records faithful deaths (schema >= 2)", False,
               "this CapTOP output predates stage7.2; alive-at-zero deaths are "
               "synthetic inf and the complete-diagram check cannot pass. Re-run "
               "CapTOP persist with the stage7.2 binary.")
        # continue anyway so the user sees the rest, but it will FAIL.

    # 1) recompute from the field
    try:
        field = fold.cube_field()
        out = field_to_diagrams(field, fold.nx, fold.ny, fold.nz, max_dim)
    except Exception as exc:
        record("recompute persistence from field", False, str(exc))
        return report
    record("recompute persistence from field", True,
           f"{fold.nx}x{fold.ny}x{fold.nz} = {field.size} cells")

    # 2) Betti@0 vector vs manifest (primary correctness proof)
    man_b0 = {d["dim"]: d.get("betti_at_zero")
              for d in fold.manifest.get("dimensions", [])}
    b0_ok = True
    for d in range(max_dim + 1):
        got = out["betti0"].get(d, 0)
        want = man_b0.get(d)
        ok = (want is None) or (got == want)
        b0_ok &= ok
        record(f"Betti@0 H{d}", ok, f"recomputed={got}  manifest={want}")

    # 3) COMPLETE diagram equality per dimension: finite multiset + essential count
    diag_ok = True
    for d in range(max_dim + 1):
        got_fin = _diagram_bands(out["diagrams"].get(d, np.empty((0, 2))), rel_tol)
        want_fin = _diagram_bands(fold.diagrams.get(d, np.empty((0, 2))), rel_tol)
        fin_ok, why = _bands_match(got_fin, want_fin, rel_tol)
        ess_got = out["essentials"].get(d, 0)
        ess_want = fold.essentials.get(d, 0)
        ess_ok = (ess_got == ess_want)
        ok = fin_ok and ess_ok
        diag_ok &= ok
        detail = (f"finite: recomputed {sum(got_fin.values())} vs CapTOP "
                  f"{sum(want_fin.values())}; essential: {ess_got} vs {ess_want}")
        if not fin_ok:
            detail += f"  [{why}]"
        record(f"H{d} complete diagram", ok, detail)

    report["ok"] = b0_ok and diag_ok
    return report


# ===========================================================================
# Stage 1+ : STUBS - named homes for the already-designed stages.
# ===========================================================================
def captop_sdt(occ_flat, nx, ny, nz, h, signed, square):
    """Reproduce CapTOP's signed distance transform from a FLAT occupancy array in
    CapTOP index order (i + nx*(j + ny*k)). Returns the flat float64 field in the
    same order.

    Convention (identical to captop.cpp build_sdt_field):
      * physical Euclidean distance (Angstrom) between cube centers, scaled by h;
      * NEGATIVE inside the occupied solid, POSITIVE in the void;
      * unsigned (signed=False): occupied cells -> 0, distance grows into the void;
      * signed  (signed=True) : occupied -> -(distance to nearest empty),
                                empty    -> +(distance to nearest occupied);
      * square  (square=True) : value = sign(d) * d^2  (Angstrom^2), applied last.

    scipy's distance_transform_edt(x) returns, at each cell, the distance to the
    nearest ZERO (background) cell of x. So distance-to-nearest-occupied = EDT(~occ)
    and distance-to-nearest-empty = EDT(occ). The (nz,ny,nx) reshape with C-order
    makes the last axis (i) vary fastest, matching CapTOP's flat order exactly
    (verified: reshape->flatten round-trips and r[k,j,i] == flat[i+nx*(j+ny*k)]).
    """
    from scipy import ndimage
    occ = np.asarray(occ_flat, bool).reshape(nz, ny, nx)     # [k, j, i]
    out = ndimage.distance_transform_edt(~occ) * h           # 0 on occupied, >0 void
    if not signed:
        field = out
    else:
        inn = ndimage.distance_transform_edt(occ) * h        # 0 on void, >0 inside
        field = np.where(occ, -inn, out)                     # negative inside
    if square:
        field = np.sign(field) * field * field
    return np.ascontiguousarray(field.reshape(-1), dtype=float)


def _sdt_flags_from_manifest(fold):
    """Read (signed, square, h) for the fold's SDT from its manifest. Raises if the
    fold is not an SDT filtration (the floor is only defined for SDT)."""
    filt = fold.manifest.get("filtration", {})
    sdt = filt.get("sdt")
    if not sdt or not sdt.get("enabled"):
        raise RuntimeError(
            f"fold '{fold.label}': perturbation floor is only defined for an SDT "
            f"filtration; manifest filtration.sdt is absent/disabled")
    signed = bool(sdt.get("signed"))
    square = bool(sdt.get("square"))
    if fold.h is None:
        raise RuntimeError(f"fold '{fold.label}': spacing h unknown; cannot scale SDT")
    return signed, square, float(fold.h)


def assert_sdt_reproduces_field(fold, rtol=1e-6, atol=1e-6):
    """Prove the Python SDT reproduces CapTOP's stored field on the UNPERTURBED
    occupancy, before any perturbed field is trusted. This is the Stage-1 analogue
    of the Stage-0 self-check: it guarantees the perturbation measures real boundary
    sensitivity, not a convention mismatch between captop.cpp and scipy.

    Returns (ok, max_abs_diff). Compares only where both are finite (the field is
    finite everywhere for a valid SDT, so that is the whole grid)."""
    signed, square, h = _sdt_flags_from_manifest(fold)
    occ = fold.occupancy().astype(bool)
    stored = fold.cube_field()
    recomputed = captop_sdt(occ, fold.nx, fold.ny, fold.nz, h, signed, square)
    finite = np.isfinite(stored) & np.isfinite(recomputed)
    if not finite.any():
        return False, np.inf
    diff = np.abs(stored[finite] - recomputed[finite])
    max_abs = float(diff.max())
    scale = np.maximum(np.abs(stored[finite]), np.abs(recomputed[finite]))
    ok = bool(np.all(diff <= atol + rtol * scale))
    return ok, max_abs


def _bottleneck(diag_a, diag_b):
    """Bottleneck distance between two finite diagrams (N,2 arrays). Empty-safe."""
    import gudhi
    A = diag_a if len(diag_a) else np.empty((0, 2))
    B = diag_b if len(diag_b) else np.empty((0, 2))
    return float(gudhi.bottleneck_distance(A, B))


def perturbation_floor(fold, max_dim=2, voxels=1, floor_axis="unsigned-linear",
                       verbose=True):
    """Stage 1. Selector-stability noise floor: the cubical analogue of
    alpha_fold_compare's control 2. Nudge the occupancy boundary by +/- `voxels`
    (binary dilation and erosion), recompute the SDT the way CapTOP does, rebuild
    the diagram via the Stage-0-proven field_to_diagrams, and bottleneck each
    perturbed diagram against the unperturbed one. The floor per dimension is the
    WORST (max) over the two perturbation directions.

    floor_axis selects the axis the floor is computed on, which MUST match the axis
    the eventual fold-vs-fold comparison uses:
      'unsigned-linear' (default) : d in Angstrom, unsigned  -> for fold-vs-fold
      'signed-linear'             : d in Angstrom, signed
      'signed-square'             : sign(d)*d^2 in Angstrom^2 -> capsid-vs-mesh axis
      'as-stored'                 : whatever the fold's own manifest says
    The occupancy is the same regardless of axis; only the SDT mapping differs, so
    the floor is derived from occupancy + the chosen axis, independent of how the
    stored .raw was squared.

    Requires grid files (--write-grid). Returns dict:
      {'floor': {d: bottleneck}, 'per_direction': {...}, 'axis': ..., 'voxels': ...,
       'sdt_check': (ok, max_abs_diff)}
    """
    from scipy import ndimage
    if not fold.has_grid:
        raise RuntimeError(
            f"fold '{fold.label}': perturbation floor needs grid files; re-run "
            f"CapTOP persist with --write-grid")

    # axis -> (signed, square); h from the fold
    _, _, h = _sdt_flags_from_manifest(fold)
    axis_map = {
        "unsigned-linear": (False, False),
        "signed-linear":   (True,  False),
        "signed-square":   (True,  True),
    }
    if floor_axis == "as-stored":
        signed, square, _ = _sdt_flags_from_manifest(fold)
    elif floor_axis in axis_map:
        signed, square = axis_map[floor_axis]
    else:
        raise ValueError(f"unknown floor_axis '{floor_axis}'")

    # 0) prove the SDT reproducer matches CapTOP on the unperturbed occupancy, in
    #    the fold's OWN stored axis (that is the only axis we can cross-check against
    #    the .raw). A convention bug would surface here before we trust perturbations.
    sdt_ok, sdt_maxdiff = assert_sdt_reproduces_field(fold)
    if verbose:
        print(f"  SDT reproducer vs stored field: "
              f"{'OK' if sdt_ok else 'MISMATCH'} (max |diff| = {sdt_maxdiff:.3g})")
    if not sdt_ok:
        raise RuntimeError(
            f"fold '{fold.label}': Python SDT does not reproduce CapTOP's stored "
            f"field (max |diff| = {sdt_maxdiff:.3g}); the perturbation floor would "
            f"be measuring a convention mismatch, not boundary sensitivity. Fix "
            f"captop_sdt to match captop.cpp build_sdt_field before trusting Stage 1.")

    occ = fold.occupancy().astype(bool).reshape(fold.nz, fold.ny, fold.nx)

    # unperturbed diagram IN THE CHOSEN AXIS, via the shared recompute path so the
    # floor's base diagram is bit-identical to what compare_folds uses for this fold.
    base = fold_diagrams_on_axis(fold, floor_axis, max_dim)

    struct = ndimage.generate_binary_structure(3, 1)   # 6-connectivity (faces)
    per_direction = {}
    # Two floors computed in the SAME pass (the perturbed diagrams are already in
    # hand, so measuring Wasserstein-to-self alongside bottleneck-to-self is just an
    # extra distance call, no extra persistence). 'floor' stays the bottleneck floor
    # for backward compatibility; 'floor_wass' is the Wasserstein floor.
    floor = {d: 0.0 for d in range(max_dim + 1)}
    floor_wass = {d: 0.0 for d in range(max_dim + 1)}
    for name, op in (("dilate", ndimage.binary_dilation),
                     ("erode", ndimage.binary_erosion)):
        occ_p = op(occ, structure=struct, iterations=voxels)
        # guard: perturbation must keep both phases (else SDT is undefined)
        n_occ = int(occ_p.sum()); n_tot = occ_p.size
        if n_occ == 0 or n_occ == n_tot:
            if verbose:
                print(f"  [skip] {name}: perturbation removed a phase "
                      f"(occupied={n_occ}/{n_tot})")
            per_direction[name] = None
            continue
        field_p = captop_sdt(occ_p.reshape(-1), fold.nx, fold.ny, fold.nz,
                             h, signed, square)
        pert = field_to_diagrams(field_p, fold.nx, fold.ny, fold.nz, max_dim)
        dists = {}; dists_w = {}
        for d in range(max_dim + 1):
            A = base["diagrams"].get(d, np.empty((0, 2)))
            B = pert["diagrams"].get(d, np.empty((0, 2)))
            bn = _bottleneck(A, B)
            wd = _wasserstein(A, B)
            dists[d] = bn
            dists_w[d] = wd
            floor[d] = max(floor[d], bn)
            if not np.isnan(wd):
                floor_wass[d] = max(floor_wass[d], wd)
        per_direction[name] = {"bottleneck": dists, "wasserstein": dists_w}
        if verbose:
            cells = "  ".join(f"H{d}={dists[d]:.4g}" for d in range(max_dim + 1))
            wcells = "  ".join(f"H{d}={dists_w[d]:.4g}" for d in range(max_dim + 1))
            print(f"  {name} (+/-{voxels} vox): bn[{cells}]  w[{wcells}]")

    return {"floor": floor, "floor_wass": floor_wass,
            "per_direction": per_direction, "axis": floor_axis,
            "voxels": voxels, "sdt_check": (sdt_ok, sdt_maxdiff)}


def _axis_to_flags(fold, axis):
    """Map an axis name to (signed, square, h) for a fold. 'as-stored' reads the
    fold's own manifest SDT flags; the named axes override them. h always comes
    from the fold (the occupancy grid is axis-independent)."""
    _, _, h = _sdt_flags_from_manifest(fold)
    named = {
        "unsigned-linear": (False, False),
        "signed-linear":   (True,  False),
        "signed-square":   (True,  True),
    }
    if axis == "as-stored":
        signed, square, _ = _sdt_flags_from_manifest(fold)
        return signed, square, h
    if axis in named:
        s, q = named[axis]
        return s, q, h
    raise ValueError(f"unknown axis '{axis}'")


def fold_diagrams_on_axis(fold, axis, max_dim=2):
    """Recompute a fold's per-dimension diagrams from its OCCUPANCY on the given
    axis, via the Stage-0-proven field_to_diagrams. THE single path by which any
    diagram enters a distance computation (both the cross-fold comparison and the
    perturbation floor route through here), so the compared diagrams and the floor
    are guaranteed to live on the identical axis by construction.

    Uses the CapTOP-exact SDT reproducer, which the Stage-1 self-assertion proves
    reproduces CapTOP's stored field bit-exactly on the fold's own axis."""
    if not fold.has_grid:
        raise RuntimeError(
            f"fold '{fold.label}': needs grid files (--write-grid) to recompute "
            f"diagrams on a chosen axis")
    signed, square, h = _axis_to_flags(fold, axis)
    occ = fold.occupancy().astype(bool)
    fld = captop_sdt(occ, fold.nx, fold.ny, fold.nz, h, signed, square)
    return field_to_diagrams(fld, fold.nx, fold.ny, fold.nz, max_dim)


def _wasserstein(diag_a, diag_b, order=1.0, internal_p=2.0):
    """1-Wasserstein between two finite diagrams (empty-safe). Returns NaN (with a
    one-time note) if the optimal-transport backend (POT, `import ot`) is missing,
    since Wasserstein is corroborating signal only - the verdict gates on bottleneck,
    which has no such dependency."""
    try:
        from gudhi.wasserstein import wasserstein_distance
    except (ImportError, ModuleNotFoundError):
        if not _wasserstein._warned:
            print("[note] Wasserstein unavailable (POT/`ot` not installed); "
                  "reporting bottleneck only. `pip install pot` to enable.",
                  file=sys.stderr)
            _wasserstein._warned = True
        return float("nan")
    A = diag_a if len(diag_a) else np.empty((0, 2))
    B = diag_b if len(diag_b) else np.empty((0, 2))
    try:
        return float(wasserstein_distance(A, B, order=order, internal_p=internal_p))
    except (ImportError, ModuleNotFoundError):
        if not _wasserstein._warned:
            print("[note] Wasserstein unavailable (POT/`ot` not installed); "
                  "reporting bottleneck only. `pip install pot` to enable.",
                  file=sys.stderr)
            _wasserstein._warned = True
        return float("nan")


_wasserstein._warned = False


def compute_fold_floor(fold, axis, max_dim=2, voxels=1, verbose=False, _cache={}):
    """Stage 1 floor for one fold on `axis`, COMPUTED ONCE and cached by
    (path, axis, voxels). Stage 3 combines per-fold floors per pair, so a fold that
    appears in many pairs pays for its floor only once. Returns a dict with both
    metric floors: {'bottleneck': {d: ...}, 'wasserstein': {d: ...}}."""
    key = (os.path.abspath(fold.path), axis, voxels, max_dim)
    if key in _cache:
        return _cache[key]
    res = perturbation_floor(fold, max_dim=max_dim, voxels=voxels,
                             floor_axis=axis, verbose=verbose)
    both = {"bottleneck": res["floor"], "wasserstein": res["floor_wass"]}
    _cache[key] = both
    return both


def _renorm_power(axis):
    """The exponent p such that a fold's diagram coordinates scale as h**p on the
    given axis: linear axes carry physical distance (h**1); the signed-square axis
    carries sign*distance**2 (h**2). A common-h renormalization multiplies a fold's
    diagram (and its floor, which is a bottleneck of such diagrams) by
    (h_ref/h_fold)**p."""
    return 2 if axis == "signed-square" else 1


def _renorm_factor(h_fold, h_ref, axis):
    """Scalar that maps a fold's diagram/floor from its own h_fold units to the
    shared h_ref units. Returns 1.0 when renormalization is off (h_ref is None) or
    h_fold is unusable, so callers can multiply unconditionally."""
    if h_ref is None or not h_fold or h_fold <= 0.0:
        return 1.0
    return (h_ref / h_fold) ** _renorm_power(axis)


def _scale_diagrams(diagrams, s):
    """Return a new {dim: (N,2) array} with birth/death multiplied by scalar s.
    s == 1.0 short-circuits to the original arrays (no copy)."""
    if s == 1.0:
        return diagrams
    return {d: (arr * s if len(arr) else arr) for d, arr in diagrams.items()}


def _check_comparable(fold_a, fold_b, h_rtol=DEFAULT_H_RTOL):
    """Hard gate before any cross-fold distance. Returns (ok, reason). Requires
    matching comparability_key AND h equal within a RELATIVE tolerance h_rtol (the
    key intentionally omits h, so h is checked separately; genuinely cross-h
    distances measure discretization, not folds). The tolerance exists because two
    folds at the same nominal resolution but voxelized independently differ in h at
    the ~4th significant figure (export/orientation quantization); a 2x octree-level
    change is ~100%, far above any sane h_rtol, so real resolution mismatches are
    still rejected."""
    if fold_a.comparability_key != fold_b.comparability_key:
        return False, (f"comparability_key mismatch: "
                       f"'{fold_a.comparability_key}' vs '{fold_b.comparability_key}' "
                       f"(different filtration/units/mode/sdt settings)")
    if fold_a.h is None or fold_b.h is None:
        return False, "spacing h unknown for one fold"
    denom = max(abs(fold_a.h), abs(fold_b.h), 1e-300)
    rel = abs(fold_a.h - fold_b.h) / denom
    if rel > h_rtol:
        return False, (f"resolution mismatch: h={fold_a.h:g} vs h={fold_b.h:g} "
                       f"(relative {rel:.3g} > --h-rtol {h_rtol:g}); fold-vs-fold "
                       f"requires matching h. If these are the same nominal "
                       f"resolution and the gap is export/orientation quantization, "
                       f"raise --h-rtol (a true octree-level change is ~1.0). "
                       f"Finite-bar counts scale with resolution, so genuinely "
                       f"cross-h distances measure the mesh, not the fold.")
    return True, ""


def compare_folds(fold_a, fold_b, axis="unsigned-linear", max_dim=2,
                  h_rtol=DEFAULT_H_RTOL, h_ref=None):
    """Stage 2. Per-dimension bottleneck + 1-Wasserstein between two folds, on
    `axis`. HARD-REQUIRES matching comparability_key AND h within h_rtol. Both
    folds' diagrams are recomputed on the SAME axis via fold_diagrams_on_axis, so
    the distances and the Stage-1 floor (also on `axis`) are guaranteed consistent.

    If h_ref is given (common-h renormalization on), each fold's diagram is scaled
    to the shared h_ref units by (h_ref/h_fold)**p before distances, removing the
    sub-tolerance scale skew between two independently-meshed folds. h_ref is None
    (off) by default; the skew it removes is ~1000x below the perturbation floor, so
    verdicts are unchanged either way. When on, the CALLER must scale each fold's
    floor by the same factor so distance and floor stay in the same units (verdict
    does this).

    Returns {d: {'bottleneck':bn, 'wasserstein':wd, 'npts_a':.., 'npts_b':..,
                 'ess_a':.., 'ess_b':..}}."""
    ok, why = _check_comparable(fold_a, fold_b, h_rtol=h_rtol)
    if not ok:
        raise RuntimeError(f"cannot compare '{fold_a.label}' vs '{fold_b.label}': {why}")
    da = fold_diagrams_on_axis(fold_a, axis, max_dim)
    db = fold_diagrams_on_axis(fold_b, axis, max_dim)
    sa = _renorm_factor(fold_a.h, h_ref, axis)
    sb = _renorm_factor(fold_b.h, h_ref, axis)
    diag_a = _scale_diagrams(da["diagrams"], sa)
    diag_b = _scale_diagrams(db["diagrams"], sb)
    out = {}
    for d in range(max_dim + 1):
        A = diag_a.get(d, np.empty((0, 2)))
        B = diag_b.get(d, np.empty((0, 2)))
        out[d] = {
            "bottleneck": _bottleneck(A, B),
            "wasserstein": _wasserstein(A, B),
            "npts_a": len(A), "npts_b": len(B),
            "ess_a": da["essentials"].get(d, 0), "ess_b": db["essentials"].get(d, 0),
        }
    return out


def verdict(folds, axis="unsigned-linear", max_dim=2, voxels=1, verbose=True,
            h_rtol=DEFAULT_H_RTOL, renorm_h=False, metric="bottleneck"):
    """Stage 3. For every pair of folds: per-dimension bottleneck/Wasserstein
    (Stage 2), a per-pair floor = max of the two folds' Stage-1 floors (each
    computed once and cached), and the verdict 'REAL fold difference' iff ANY
    dimension's distance clears its own pair floor (the permissive rule mirroring
    alpha_fold_compare). Per-dimension clears/noise is reported so the discriminating
    dimension is visible.

    metric selects which distance drives the verdict AND which floor it is gated
    against - they MUST be the same metric (a Wasserstein distance cannot be gated
    by a bottleneck floor):
      'bottleneck'  (default): worst-single-point distance vs bottleneck floor.
                    Conservative; on voxelized diagrams the single worst-matched
                    point is noisy, so signals tend to sit near the floor.
      'wasserstein': whole-diagram transport distance vs Wasserstein floor.
                    Integrates the entire cavity/channel population, so it is more
                    stable and (empirically) more discriminating for folds, at the
                    cost of being sensitive to the swarm of small bars.
    Both floors are computed in the same perturbation pass, so switching metric
    costs nothing extra.

    h_rtol is the relative spacing-match tolerance (see _check_comparable). If
    renorm_h is True, all folds are renormalized to a single shared h_ref (the mean
    of the folds' h) before distances: each fold's diagram AND its floor are scaled
    by (h_ref/h_fold)**p, keeping distance and floor in identical units.

    Returns a list of per-pair result dicts."""
    import itertools
    if metric not in ("bottleneck", "wasserstein"):
        raise ValueError(f"metric must be 'bottleneck' or 'wasserstein', got '{metric}'")
    # gate all pairs up front so we fail fast on any incompatibility
    labels = [f.label for f in folds]
    for a, b in itertools.combinations(folds, 2):
        ok, why = _check_comparable(a, b, h_rtol=h_rtol)
        if not ok:
            raise RuntimeError(f"cannot compare '{a.label}' vs '{b.label}': {why}")

    # common-h renormalization target: mean of the folds' spacings (only when on)
    h_ref = None
    if renorm_h:
        hs = [f.h for f in folds if f.h]
        h_ref = (sum(hs) / len(hs)) if hs else None
        if verbose and h_ref is not None:
            print(f"common-h renormalization ON: h_ref={h_ref:.6g} Angstrom "
                  f"(diagrams and floors scaled by (h_ref/h_fold)**p per axis)")

    # per-fold floors, computed once each (cached), then scaled to h_ref units.
    # Select the floor for the chosen metric.
    if verbose:
        print(f"per-fold selector-stability floors (metric={metric}, "
              f"axis={axis}, +/-{voxels} vox), computed once each:")
    floors = {}
    for f in folds:
        raw_both = compute_fold_floor(f, axis, max_dim, voxels, verbose=False)
        raw = raw_both[metric]
        s = _renorm_factor(f.h, h_ref, axis)
        floors[f.label] = {d: raw[d] * s for d in range(max_dim + 1)}
        if verbose:
            fl = floors[f.label]
            print(f"  {f.label:<20} "
                  + "  ".join(f"H{d}={fl[d]:.4g}" for d in range(max_dim + 1)))

    results = []
    for a, b in itertools.combinations(folds, 2):
        per = compare_folds(a, b, axis=axis, max_dim=max_dim,
                            h_rtol=h_rtol, h_ref=h_ref)
        pair_floor = {d: max(floors[a.label][d], floors[b.label][d])
                      for d in range(max_dim + 1)}
        clears = {d: per[d][metric] > pair_floor[d] for d in range(max_dim + 1)}
        is_real = any(clears.values())
        results.append({"a": a.label, "b": b.label, "per": per,
                        "pair_floor": pair_floor, "clears": clears,
                        "real": is_real, "axis": axis, "metric": metric})
    return results


def capsid_vs_mesh_table(alpha_diagram_dir, mesh_fold, max_dim=2):
    """Stage 4 [NOT IMPLEMENTED]. The headline capsid-vs-mesh comparison as an
    explicitly NON-metric structural table: top-feature multiplicity bands and the
    H2-count-vs-h divergence, NOT a bottleneck distance (alpha and cubical do not
    share a metric space and the counts are ~1/h^2 resolution-dependent)."""
    raise NotImplementedError(
        "Stage 4 (capsid_vs_mesh_table): structural (non-metric) comparison.")


# ===========================================================================
# CLI
# ===========================================================================
def _check_gudhi():
    try:
        import gudhi  # noqa: F401
    except ImportError:
        sys.exit("[error] gudhi not importable; activate the env you ran CapTOP "
                 "from (e.g. micromamba activate gudhi)")


def main():
    ap = argparse.ArgumentParser(
        prog="cubic_fold_compare.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
        epilog="examples:\n"
               "  # Stage 0: prove the grid contract on the small h=12 fold\n"
               "  %(prog)s --self-check captop_out/out_persist_sdt/1CWP_F2_0_h12\n"
               "\n"
               "  # Stage 2-3: compare two folds (same nominal resolution)\n"
               "  %(prog)s --compare DIR_A DIR_B --floor-axis unsigned-linear\n"
               "\n"
               "  # loosen the spacing match if a same-resolution pair is refused\n"
               "  %(prog)s --compare DIR_A DIR_B --h-rtol 2e-2\n"
               "\n"
               "  # erase the residual sub-tolerance h skew (shared-h units)\n"
               "  %(prog)s --compare DIR_A DIR_B --renorm-h\n")

    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--self-check", metavar="DIR",
                      help="STAGE 0 (active): reproduce a fold's diagram from its "
                           "own grid field (requires the fold was produced with "
                           "CapTOP persist --write-grid) and assert it matches the "
                           "CapTOP-written Betti@0 vector and finite-bar lifetimes. "
                           "This proves the field->complex->diagram conversion the "
                           "whole module reuses. Exit 0 on PASS, 1 on FAIL.")
    mode.add_argument("--compare", nargs="+", metavar="DIR",
                      help="STAGE 2-3 (active): per-dimension bottleneck/Wasserstein "
                           "between folds on --floor-axis (default unsigned-linear), "
                           "gated by matching comparability_key AND h within "
                           "--h-rtol (default 1e-2; admits export/orientation "
                           "quantization but still rejects a true octree-level "
                           "change), with a per-pair selector-stability floor (each "
                           "fold's Stage-1 floor computed once and cached) and a "
                           "'REAL iff any dimension clears its pair floor' verdict. "
                           "Optionally renormalize all folds to a shared h with "
                           "--renorm-h. Needs --write-grid output for every fold.")
    mode.add_argument("--floor", metavar="DIR",
                      help="STAGE 1 (active): compute the selector-stability noise "
                           "floor for one fold - dilate/erode the occupancy by +/-1 "
                           "voxel, recompute the SDT, rebuild the diagram, and "
                           "bottleneck each perturbation against the unperturbed "
                           "diagram (worst per dimension). Requires --write-grid "
                           "output. First asserts the Python SDT reproduces CapTOP's "
                           "stored field. This is the cubical analogue of "
                           "alpha_fold_compare control 2.")
    mode.add_argument("--info", metavar="DIR",
                      help="load a fold and print what was found (manifest units, "
                           "comparability_key, h, per-dim bar counts, whether grid "
                           "files are present). Pure I/O; no topology, no gudhi.")

    ap.add_argument("--max-dim", type=int, default=2,
                    help="highest homology dimension to handle (default 2)")
    ap.add_argument("--rel-tol", type=float, default=1e-4,
                    help="relative tolerance for matching finite-bar lifetimes in "
                         "the self-check (default 1e-4; mirrors the alpha-side "
                         "relative-tolerance band collapse)")
    ap.add_argument("--floor-axis", default="unsigned-linear",
                    choices=["unsigned-linear", "signed-linear", "signed-square",
                             "as-stored"],
                    help="axis on which to compute the perturbation floor; MUST "
                         "match the axis of the eventual fold-vs-fold comparison "
                         "(default unsigned-linear, the fold-vs-fold axis)")
    ap.add_argument("--floor-voxels", type=int, default=1,
                    help="perturbation size in voxels for --floor (default 1)")
    ap.add_argument("--h-rtol", type=float, default=DEFAULT_H_RTOL,
                    help="relative tolerance for the fold-vs-fold spacing (h) match "
                         "(default %(default)g). Two same-resolution folds voxelized "
                         "independently differ in h at the ~4th significant figure; "
                         "this admits that quantization while a true octree-level "
                         "change (~1.0 relative) is still rejected. Raise it if a "
                         "same-nominal-resolution pair is refused; lower it to be "
                         "stricter.")
    ap.add_argument("--renorm-h", action="store_true",
                    help="renormalize all --compare folds to a shared h (the mean of "
                         "their spacings) before distances: each fold's diagram and "
                         "floor are scaled by (h_ref/h_fold)**p, removing the "
                         "sub-tolerance scale skew between independently-meshed folds. "
                         "Off by default (the skew is ~1000x below the noise floor "
                         "and cannot change a verdict); the fold-vs-itself null "
                         "control stays exactly zero either way.")
    ap.add_argument("--metric", default="bottleneck",
                    choices=["bottleneck", "wasserstein"],
                    help="distance driving the --compare verdict, gated against the "
                         "matching-metric floor (default bottleneck). 'bottleneck' is "
                         "the conservative worst-single-point distance; 'wasserstein' "
                         "integrates the whole diagram and is empirically more "
                         "discriminating for folds. Both floors are computed in one "
                         "pass, so switching is free.")
    args = ap.parse_args()
    if args.info:
        fd = load_fold(args.info)
        print(f"fold '{fd.label}'  ({fd.path})")
        print(f"  units              : {fd.units}")
        print(f"  comparability_key  : {fd.comparability_key}")
        print(f"  spacing h          : {fd.h}")
        print(f"  grid files present : {fd.has_grid}"
              + (f"  ({fd.nx}x{fd.ny}x{fd.nz})" if fd.has_grid else ""))
        for d in sorted(fd.diagrams):
            print(f"  H{d}: {len(fd.diagrams[d])} finite bars, "
                  f"{fd.essentials[d]} essential")
        return

    # ---- --self-check : the Stage 0 gate (needs gudhi) ----
    if args.self_check:
        fd = load_fold(args.self_check)
        if not fd.has_grid:
            sys.exit(f"[error] fold '{fd.label}' has no grid files; re-run CapTOP "
                     f"persist with --write-grid to enable the self-check.")
        _check_gudhi()
        print(f"self-check fold '{fd.label}'  "
              f"(units={fd.units}, key={fd.comparability_key}, h={fd.h})")
        rep = self_check(fd, max_dim=args.max_dim, rel_tol=args.rel_tol,
                         verbose=True)
        print(f"\nself-check: {'PASS' if rep['ok'] else 'FAIL'}")
        if not rep["ok"]:
            print("  -> the field->diagram conversion does NOT reproduce CapTOP's "
                  "diagram. If H1/H2 look swapped, the index/axis order in "
                  "field_to_diagrams needs a transpose (fix there, only there). "
                  "Do not build later stages until this passes.", file=sys.stderr)
        sys.exit(0 if rep["ok"] else 1)

    # ---- --floor : Stage 1 perturbation floor (needs gudhi + scipy) ----
    if args.floor:
        _check_gudhi()
        fd = load_fold(args.floor)
        if not fd.has_grid:
            sys.exit(f"[error] fold '{fd.label}' has no grid files; re-run CapTOP "
                     f"persist with --write-grid to enable the perturbation floor.")
        print(f"perturbation floor for fold '{fd.label}'  "
              f"(axis={args.floor_axis}, +/-{args.floor_voxels} vox, "
              f"units={fd.units}, h={fd.h})")
        res = perturbation_floor(fd, max_dim=args.max_dim, voxels=args.floor_voxels,
                                 floor_axis=args.floor_axis, verbose=True)
        floor = res["floor"]; floor_w = res["floor_wass"]
        print("\n  selector-stability floor per dim (worst over dilate/erode):")
        print("    bottleneck : " + "   ".join(f"H{d}={floor[d]:.4g}"
                                                for d in range(args.max_dim + 1)))
        print("    wasserstein: " + "   ".join(f"H{d}={floor_w[d]:.4g}"
                                                for d in range(args.max_dim + 1)))
        print("\n  a cross-fold distance must EXCEED the matching-metric floor to "
              "count as real fold signal (Stage 3 verdict; --metric selects which).")
        return

    # ---- --compare : Stage 2-3 (needs gudhi + scipy) ----
    if args.compare:
        if len(args.compare) < 2:
            sys.exit("[error] --compare needs at least two fold directories")
        _check_gudhi()
        folds = [load_fold(p) for p in args.compare]
        for f in folds:
            if not f.has_grid:
                sys.exit(f"[error] fold '{f.label}' has no grid files; re-run "
                         f"CapTOP persist with --write-grid (needed for the "
                         f"axis-consistent recompute and the floor).")
        print(f"comparing {len(folds)} folds on axis={args.floor_axis}, "
              f"floor +/-{args.floor_voxels} vox, metric={args.metric}")
        print("=" * 72)
        try:
            results = verdict(folds, axis=args.floor_axis, max_dim=args.max_dim,
                              voxels=args.floor_voxels, verbose=True,
                              h_rtol=args.h_rtol, renorm_h=args.renorm_h,
                              metric=args.metric)
        except RuntimeError as exc:
            sys.exit(f"[error] {exc}")

        print("\n" + "=" * 72)
        print("cross-fold distances (bottleneck bn, 1-Wasserstein w) on "
              f"axis={args.floor_axis}")
        print("=" * 72)
        for r in results:
            cells = " | ".join(
                f"H{d}: bn={r['per'][d]['bottleneck']:.4g} "
                f"w={r['per'][d]['wasserstein']:.4g}"
                for d in range(args.max_dim + 1))
            print(f"  {r['a']} vs {r['b']}\n    {cells}")

        mlabel = args.metric
        print("\n" + "=" * 72)
        print(f"verdict (metric={mlabel}: REAL iff ANY dimension's {mlabel} "
              f"clears its pair floor)")
        print("=" * 72)
        for r in results:
            tag = ("REAL fold difference" if r["real"]
                   else "within noise floor (not significant)")
            print(f"  {r['a']} vs {r['b']}  ->  {tag}")
            print("      " + "   ".join(
                f"H{d}: {r['per'][d][mlabel]:.4g} vs "
                f"floor {r['pair_floor'][d]:.4g} "
                f"{'CLEARS' if r['clears'][d] else 'noise'}"
                for d in range(args.max_dim + 1)))
        return


if __name__ == "__main__":
    main()
    
