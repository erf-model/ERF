#!/usr/bin/env python3
"""Check the face balance's refined level against a refined level that holds both buildings.

    check_ibseb_refined_levels.py --both <dir> --one <dir>

Each directory holds the step-0 face dumps of one run of
Tests/test_files/IBSEB_RefinedLevels: <dir>/faces/set.lev1.rank*.csv (level 1)
and <dir>/faces/set.rank*.csv (level 0). In --both level 1 holds the cube and
the tower; in --one it holds the cube only, so the tower reaches its rays
only through the column map of level 0 (IBFaceSet::add_outside_occluders()).

Asserts, each naming the defect it guards:
  1. the level-1 dump of --one exists and is not level 0's (the dumps of a
     refined level must not overwrite the level below);
  2. --one's level 1 has one building, --both's two, and both have the same
     cube faces (the same face list from the same cells);
  3. the tower shades the cube's east-facing wall in --both (more than a
     quarter of its faces), so the comparison below is not vacuous;
  4. in --one the cube's faces see the same shadow on every face as in
     --both: without the coarse column map the east-facing wall comes out
     sunlit;
  5. their view fractions agree within 5 of the 128 hemisphere rays (the
     deck's 16 x 8), the difference the coarser tower makes: level 0 resolves
     it at 20 m with a ramp cell on every side, which catches a few more rays.
     Without the coarse column map the cube's faces lose 19 rays to the sky.
"""
import argparse
import glob
import sys


def load(pattern):
    """The rows of every rank's dump as {column: list}, or None without files."""
    rows, hdr = [], None
    for fn in sorted(glob.glob(pattern)):
        with open(fn) as f:
            hdr = f.readline().strip().split(",")
            rows += [[float(v) for v in line.split(",")] for line in f if line.strip()]
    if not rows:
        return None
    return {h: [r[n] for r in rows] for n, h in enumerate(hdr)}


def faces(d, bid):
    return {(int(d["i"][n]), int(d["j"][n]), int(d["k"][n]), int(d["dir"][n]), int(d["side"][n])): n
            for n in range(len(d["bid"])) if int(d["bid"][n]) == bid}


def mean(v):
    return sum(v) / len(v) if v else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--both", required=True)
    ap.add_argument("--one", required=True)
    a = ap.parse_args()
    ok = True

    def check(name, cond, detail):
        nonlocal ok
        print(f"  {name}: {'PASS' if cond else 'FAIL'} ({detail})")
        ok &= bool(cond)

    both = load(a.both + "/faces/set.lev1.rank*.csv")
    one = load(a.one + "/faces/set.lev1.rank*.csv")
    one0 = load(a.one + "/faces/set.rank*.csv")
    check("1. level 1 dumps its own faces", one is not None and both is not None and one0 is not None
          and len(one["i"]) != len(one0["i"]),
          "missing" if one is None or one0 is None else f"level 1 {len(one['i'])} faces, level 0 {len(one0['i'])}")
    if one is None or both is None:
        return 1
    nb_one, nb_both = int(max(one["bid"])), int(max(both["bid"]))
    f1, fb = faces(one, 1), faces(both, 1)
    check("2. buildings and cube faces on level 1", nb_one == 1 and nb_both == 2 and set(f1) == set(fb),
          f"{nb_one} and {nb_both} buildings, {len(f1)} and {len(fb)} cube faces")
    common = sorted(set(f1) & set(fb))
    # The wall whose solid is at i - 1 has its fluid to the east: it faces the sun.
    east = [c for c in common if c[3] == 0 and c[4] == -1]
    sh_b = [both["shadow"][fb[c]] for c in east]
    check("3. the tower shades the cube's east-facing wall", len(east) > 0 and mean(sh_b) > 0.25,
          f"{len(east)} faces, shadow fraction {mean(sh_b):.3f}")
    s1 = [one["shadow"][f1[c]] for c in common]
    sb = [both["shadow"][fb[c]] for c in common]
    ndiff = sum(1 for u, v in zip(s1, sb) if u != v)
    check("4. the same shadow on every cube face", len(common) > 0 and ndiff == 0,
          f"{ndiff} of {len(common)} faces differ; east-facing wall "
          f"{mean([one['shadow'][f1[c]] for c in east]):.3f} vs {mean(sh_b):.3f}")
    worst = 0.0
    for q in ("f_sky", "f_ground", "f_bldg"):
        worst = max([worst] + [abs(one[q][f1[c]] - both[q][fb[c]]) for c in common])
    nrays = 16 * 8
    check("5. view fractions within 5 of the 128 rays", worst <= 5.0 / nrays + 1e-12,
          f"largest difference {worst:.4f} = {worst * nrays:.0f} rays")
    print("ALL PASS" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
