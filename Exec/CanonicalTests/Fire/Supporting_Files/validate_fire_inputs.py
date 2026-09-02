#!/usr/bin/env python3
"""Static validation of the fire canonical-test input files.

Checks each ``inputs_*`` file against the constraints ERF and AMReX enforce at
startup, so a broken case is caught here instead of aborting several minutes
into a run.  Nothing is executed; this only parses the input files.

Constraints checked
-------------------
AMReX (``AmrMesh::checkInput``)
  * ``amr.n_cell`` divisible by ``amr.blocking_factor``.  ERF resets that
    default to 1 in ``Source/main.cpp``, so this only bites a case that sets
    the blocking factor itself.

Fire prerequisites (``verify_fire_prerequisites``, Source/Fire)
  * box x/y lengths divisible by ``erf.fire.grid_ratio``
  * no z-decomposition: the z box length must equal ``amr.n_cell`` z, which
    means ``amr.max_grid_size_z`` has to be set whenever n_cell z exceeds the
    AMReX default (32 on CPU, 64 on GPU)
  * domain height greater than ``erf.fire.wind_ref_ht``
  * divisors read from the inputs file are positive:
    ``levelset.reinit_every``, ``spotting.spotting_interval``,
    ``levelset.cfl``, ``farsite.cfl_fire``

Ignition reachability
  * ``FireLayer::initialize`` asserts at least one cell is marked burned.  A
    circular ignition only marks a cell when its radius reaches a fire-cell
    centre, so ``ignition_r`` must be comparable to the fire cell size.

File references
  * fuel maps, polygon/polyline vertex lists, ignition schedules and terrain
    files must resolve relative to the case directory or the repository root.

The duplicate-key check applies to any ERF inputs file, so this is also worth
pointing at other canonical-test trees:

    python3 Fire/Supporting_Files/validate_fire_inputs.py Dust

Usage
-----
    python3 Supporting_Files/validate_fire_inputs.py [root]

Exits non-zero if any case fails.
"""

import os
import sys
import glob

# ERF sets its own default in Source/main.cpp ('blocking_factor = 1'), which is
# far looser than the AMReX default of 8. Use ERF's, or every case looks broken.
BLOCKING_FACTOR_DEFAULT = 1
# AMReX_AmrMesh.H: IntVect(128,128,32) on CPU, IntVect(128,128,64) on GPU.
# Use the CPU value so a case that only works on GPU is still reported.
MAX_GRID_SIZE_Z_DEFAULT = 32

FILE_KEYS = [
    "erf.fire.fuel_map.fuel_map_file",
    "erf.fire.ignition.polygon_file",
    "erf.fire.ignition.ignition_schedule_file",
    "erf.fire.terrain_file_name",
]


def parse(path):
    """Parse an ERF inputs file into {key: value-string}, stripping comments."""
    out = {}
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if "=" not in line:
                continue
            key, val = line.split("=", 1)
            out[key.strip()] = val.strip().strip('"')
    return out


def ints(d, key, n=None):
    if key not in d:
        return None
    try:
        vals = [int(float(v)) for v in d[key].split()]
    except ValueError:
        return None
    return vals[:n] if n else vals


def real(d, key, idx=0):
    if key not in d:
        return None
    try:
        return float(d[key].split()[idx])
    except (ValueError, IndexError):
        return None


def duplicate_keys(path):
    """Keys assigned more than once. AMReX aborts with 'Duplicate inputs
    detected in inputs file', so this is fatal regardless of the values."""
    seen, dups = {}, []
    with open(path, errors="replace") as fh:
        for n, line in enumerate(fh, 1):
            line = line.split("#", 1)[0].strip()
            if "=" not in line:
                continue
            key = line.split("=", 1)[0].strip()
            if key in seen:
                dups.append((key, seen[key], n))
            else:
                seen[key] = n
    return dups


def check(path, repo_root):
    d = parse(path)
    problems = []

    # Applies to every ERF inputs file, fire or not.
    for key, first, second in duplicate_keys(path):
        problems.append(
            f"'{key}' assigned twice (lines {first} and {second}); AMReX aborts "
            "with 'Duplicate inputs detected in inputs file'")

    if d.get("erf.fire.enable", "false").lower() not in ("true", "1"):
        return problems  # remaining checks are fire-specific

    n_cell = ints(d, "amr.n_cell", 3)
    if not n_cell or len(n_cell) != 3:
        problems.append("amr.n_cell missing or malformed")
        return problems
    nx, ny, nz = n_cell

    bf = ints(d, "amr.blocking_factor")
    bf = bf[0] if bf else BLOCKING_FACTOR_DEFAULT
    for name, n in (("x", nx), ("y", ny), ("z", nz)):
        if bf > 1 and n % bf:
            problems.append(
                f"amr.n_cell {name}={n} not divisible by blocking_factor {bf} "
                "(AMReX: 'domain size not divisible by blocking_factor')")

    mgs = ints(d, "amr.max_grid_size")
    mgs = mgs[0] if mgs else 128
    if bf > 1 and mgs % bf:
        problems.append(
            f"amr.max_grid_size={mgs} not divisible by blocking_factor {bf}")

    # Fire grid ratio must divide the box x/y lengths.
    C = ints(d, "erf.fire.grid_ratio")
    C = C[0] if C else 5
    if C < 1:
        problems.append(f"erf.fire.grid_ratio={C} must be >= 1")
    else:
        bx, by = min(mgs, nx), min(mgs, ny)
        if bx % C or by % C:
            problems.append(
                f"box {bx}x{by} not divisible by grid_ratio {C} "
                "(fire prerequisite: 'Box sizes not divisible by grid_ratio')")

    # No z-decomposition.
    mgsz = ints(d, "amr.max_grid_size_z")
    box_nz = mgsz[0] if mgsz else min(mgs, MAX_GRID_SIZE_Z_DEFAULT)
    if box_nz != nz:
        problems.append(
            f"z box length {box_nz} != n_cell z {nz}; set amr.max_grid_size_z = {nz} "
            "(fire prerequisite: 'Cannot decompose in z direction')")

    # Domain height vs wind reference height.
    hi = None
    if "geometry.prob_hi" in d:
        hi = real(d, "geometry.prob_hi", 2)
    elif "geometry.prob_extent" in d:
        lo = real(d, "geometry.prob_lo", 2) or 0.0
        ext = real(d, "geometry.prob_extent", 2)
        hi = lo + ext if ext is not None else None
    wref = real(d, "erf.fire.wind_ref_ht")
    wref = 6.1 if wref is None else wref
    if hi is not None and hi <= wref:
        problems.append(f"domain height {hi} <= wind_ref_ht {wref}")

    # Positive divisors.
    for key, minimum in (("erf.fire.levelset.reinit_every", 1),
                         ("erf.fire.spotting.spotting_interval", 1)):
        v = real(d, key)
        if v is not None and v < minimum:
            problems.append(f"{key}={v} must be >= {minimum} (used as a modulus divisor)")
    for key in ("erf.fire.levelset.cfl", "erf.fire.farsite.cfl_fire"):
        v = real(d, key)
        if v is not None and v <= 0.0:
            problems.append(f"{key}={v} must be > 0")

    # Ignition must actually mark a cell.
    extent_x = None
    if "geometry.prob_hi" in d and "geometry.prob_lo" in d:
        extent_x = real(d, "geometry.prob_hi", 0) - real(d, "geometry.prob_lo", 0)
    elif "geometry.prob_extent" in d:
        extent_x = real(d, "geometry.prob_extent", 0)
    uses_polygon = bool(d.get("erf.fire.ignition.polygon_file"))
    if extent_x and C >= 1 and not uses_polygon:
        fire_dx = extent_x / (nx * C)
        r = real(d, "erf.fire.ignition_r")
        r = 20.0 if r is None else r
        # Worst case the ignition centre sits half a cell from the nearest centre
        # in each direction, so the diagonal is the distance that must be covered.
        if r < 0.71 * fire_dx:
            problems.append(
                f"ignition_r={r} m is small next to the {fire_dx:.1f} m fire cell; "
                "may mark no cell and trip the 'no cells were marked as burned' assert")

    # Referenced files must resolve.
    case_dir = os.path.dirname(path)
    for key in FILE_KEYS:
        ref = d.get(key)
        if not ref:
            continue
        if not any(os.path.exists(p) for p in
                   (os.path.join(case_dir, ref), os.path.join(repo_root, ref), ref)):
            problems.append(f"{key} -> '{ref}' not found relative to the case directory")

    return problems


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..")
    root = os.path.abspath(root)
    files = sorted(f for f in
                   glob.glob(os.path.join(root, "**", "inputs_*"), recursive=True)
                   + glob.glob(os.path.join(root, "**", "inputs"), recursive=True)
                   if os.path.isfile(f))
    total_bad = 0
    for path in files:
        problems = check(path, root)
        if problems:
            total_bad += 1
            print(f"\nFAIL {os.path.relpath(path, root)}")
            for p in problems:
                print(f"      - {p}")
    print(f"\nchecked {len(files)} input files; {total_bad} with problems")
    return 1 if total_bad else 0


if __name__ == "__main__":
    sys.exit(main())
