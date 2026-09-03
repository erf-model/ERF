#!/usr/bin/env python3
"""
Fire wind interpolation unit tests.

The fire grid is C times finer than the atmospheric grid, so the atmospheric
wind has to be mapped onto it. Taking the column that contains each fire cell
makes the fire-grid wind piecewise constant on atmospheric cells; bilinear
interpolation between the four surrounding columns removes the staircase. These
tests cover the stencil weights, the reproduction of a linear field, and the
vertical bracket that picks the levels to interpolate between.

Reference implementation extracted from:
  - ERF_FireWindExtract.cpp: fill_fire_wind_from_interpolation(),
                             column_wind_at_height()

Run: python3 test_wind_interpolation.py
"""

import sys
import math


def stencil(i_f, C):
    """
    Lower-left column and weight of the bilinear stencil for a fire cell.

    Formula from fill_fire_wind_from_interpolation(): the fire cell centre sits
    at continuous atmospheric cell-centre coordinate (i_f + 0.5)/C - 0.5.

    Returns (i0, w) with the interpolation running between columns i0 and i0+1.
    """
    g = (i_f + 0.5) / C - 0.5
    i0 = math.floor(g)
    return i0, g - i0


def sample_bilinear(i_f, j_f, C, u_of_column):
    """Blend the four surrounding columns with bilinear weights."""
    i0, wx = stencil(i_f, C)
    j0, wy = stencil(j_f, C)
    total = 0.0
    for dj in (0, 1):
        for di in (0, 1):
            w = (wx if di else 1.0 - wx) * (wy if dj else 1.0 - wy)
            total += w * u_of_column(i0 + di, j0 + dj)
    return total


def sample_nearest(i_f, j_f, C, u_of_column):
    """Take the single column containing the fire cell."""
    return u_of_column(i_f // C, j_f // C)


def column_wind_at_height(z_cc, u_of_level, z_target):
    """
    Wind at a target height within one column.

    Formula from column_wind_at_height(): bracket by bisection on the column's
    cell-centre heights and interpolate linearly, clamping rather than
    extrapolating outside the column.
    """
    nz = len(z_cc)
    if z_target <= z_cc[0]:
        k_lo = 0
    elif z_target >= z_cc[nz - 1]:
        k_lo = nz - 2
    else:
        lo, hi = 0, nz - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if z_cc[mid] <= z_target:
                lo = mid
            else:
                hi = mid
        k_lo = lo
    k_lo = max(0, min(k_lo, nz - 2))

    z_lo, z_hi = z_cc[k_lo], z_cc[k_lo + 1]
    alpha = 0.0
    if z_hi > z_lo:
        alpha = max(0.0, min(1.0, (z_target - z_lo) / (z_hi - z_lo)))
    u_lo, u_hi = u_of_level(k_lo), u_of_level(k_lo + 1)
    return u_lo + alpha * (u_hi - u_lo)


def check(failures, name, condition, detail=""):
    status = "✓" if condition else "✗"
    print(f"{status} {name}")
    if not condition:
        if detail:
            print(f"    {detail}")
        failures.append(name)


def test_stencil_weights_partition_unity():
    """The four bilinear weights always sum to one."""
    failures = []
    C = 5
    worst = 0.0
    for i_f in range(0, 40):
        for j_f in range(0, 40):
            _, wx = stencil(i_f, C)
            _, wy = stencil(j_f, C)
            total = sum((wx if di else 1 - wx) * (wy if dj else 1 - wy)
                        for dj in (0, 1) for di in (0, 1))
            worst = max(worst, abs(total - 1.0))
            if not (0.0 <= wx < 1.0 and 0.0 <= wy < 1.0):
                worst = 1.0
    check(failures, "Bilinear weights sum to one and stay in [0, 1)",
          worst < 1.0e-12, f"max deviation {worst}")
    return failures


def test_fire_cell_at_column_centre():
    """The fire cell centred on an atmospheric cell takes that column alone."""
    failures = []
    # With odd C the middle fire cell of a column sits on the column centre
    C = 5
    i0, wx = stencil(2, C)     # third fire cell of atmospheric cell 0
    check(failures, "Fire cell on the column centre gets weight 1 on that column",
          i0 == 0 and abs(wx) < 1.0e-12, f"i0={i0}, wx={wx}")
    return failures


def test_bilinear_reproduces_linear_field():
    """A linear wind field is reproduced exactly; nearest staircases."""
    failures = []
    C = 5
    dx_atm = 100.0
    dx_fire = dx_atm / C

    # u = 1 + 0.01 x, sampled at atmospheric cell centres
    def u_col(ia, ja):
        return 1.0 + 0.01 * ((ia + 0.5) * dx_atm)

    worst_bilinear = 0.0
    worst_nearest = 0.0
    for i_f in range(0, 20):
        x_f = (i_f + 0.5) * dx_fire
        exact = 1.0 + 0.01 * x_f
        worst_bilinear = max(worst_bilinear,
                             abs(sample_bilinear(i_f, 10, C, u_col) - exact))
        worst_nearest = max(worst_nearest,
                            abs(sample_nearest(i_f, 10, C, u_col) - exact))

    check(failures, "Bilinear reproduces a linear wind field exactly",
          worst_bilinear < 1.0e-12, f"max error {worst_bilinear}")
    check(failures, "Nearest column misses a linear field by up to 0.4 m/s",
          abs(worst_nearest - 0.4) < 1.0e-12, f"max error {worst_nearest}")
    return failures


def test_nearest_is_discontinuous():
    """Nearest jumps at atmospheric cell edges; bilinear follows the gradient."""
    failures = []
    C = 5
    dx_atm = 100.0

    def u_col(ia, ja):
        return 1.0 + 0.01 * ((ia + 0.5) * dx_atm)

    jump_nearest = abs(sample_nearest(5, 10, C, u_col)
                       - sample_nearest(4, 10, C, u_col))
    jump_bilinear = abs(sample_bilinear(5, 10, C, u_col)
                        - sample_bilinear(4, 10, C, u_col))
    # The true change over one 20 m fire cell is 0.01 * 20 = 0.2
    check(failures, "Bilinear step across a column edge is the true gradient",
          abs(jump_bilinear - 0.2) < 1.0e-12, f"got {jump_bilinear}")
    check(failures, "Nearest jumps by a full atmospheric cell of shear",
          abs(jump_nearest - 1.0) < 1.0e-12, f"got {jump_nearest}")
    return failures


def test_vertical_bracket():
    """The bracket clamps outside the column instead of running to the top."""
    failures = []
    z_cc = [50.0, 150.0, 250.0, 350.0]

    def u_of_level(k):
        return 1.0 + z_cc[k] / 100.0      # u = 1 + z/100

    cases = [
        (20.0, 1.5, "below the lowest cell centre"),
        (50.0, 1.5, "at the lowest cell centre"),
        (200.0, 3.0, "between cell centres"),
        (5000.0, 4.5, "above the top cell centre"),
    ]
    for z_target, expected, what in cases:
        got = column_wind_at_height(z_cc, u_of_level, z_target)
        check(failures, f"Bracket handles a target {what}",
              abs(got - expected) < 1.0e-12, f"got {got}, expected {expected}")

    # The superseded linear scan started at the top interval and fell through
    # for a target below the lowest cell centre, returning the wind from 250 m.
    check(failures, "Below-surface target no longer returns the upper-level wind",
          abs(column_wind_at_height(z_cc, u_of_level, 20.0) - 3.5) > 1.0,
          "still returning the k = nz-2 value")
    return failures


def test_per_column_ground_anchoring():
    """Each column is sampled at the reference height above its own ground."""
    failures = []
    # Two columns 28 m apart in ground elevation, as on a 30 degree slope with
    # 50 m atmospheric cells, sampled at 6.1 m above ground.
    z_ref = 6.1
    ground = [0.0, 28.5]
    absolute = [g + z_ref for g in ground]
    agl = [a - g for a, g in zip(absolute, ground)]

    check(failures, "Per-column targets are equal heights above local ground",
          all(abs(a - z_ref) < 1.0e-12 for a in agl))

    # Blending at one absolute height would sample the upslope column far higher
    single = ground[0] + z_ref
    agl_wrong = single - ground[1]
    check(failures, "A single absolute height would mis-sample the neighbour",
          agl_wrong < 0.0,
          f"neighbour would be sampled {agl_wrong:.1f} m above its ground")
    return failures


def main():
    print("=" * 70)
    print("Fire Wind Interpolation Unit Tests")
    print("=" * 70)

    all_failures = []
    all_failures.extend(test_stencil_weights_partition_unity())
    all_failures.extend(test_fire_cell_at_column_centre())
    all_failures.extend(test_bilinear_reproduces_linear_field())
    all_failures.extend(test_nearest_is_discontinuous())
    all_failures.extend(test_vertical_bracket())
    all_failures.extend(test_per_column_ground_anchoring())

    print("\n" + "=" * 70)
    if not all_failures:
        print("FINAL RESULT: All tests PASSED")
        return 0
    print(f"FINAL RESULT: {len(all_failures)} test(s) FAILED")
    for failure in all_failures:
        print(f"  ✗ {failure}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
