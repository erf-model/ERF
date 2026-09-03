#!/usr/bin/env python3
"""
Terrain projection unit tests.

The rate of spread produced by every ROS model is a speed along the terrain
surface. Both propagation paths advance the front in map coordinates, so a
surface step of ds up a slope of angle theta covers only ds*cos(theta)
horizontally. These tests cover the projection factor used by the level-set
gradient operator and by the FARSITE displacement, and the datum used for wind
extraction.

Reference implementations extracted from:
  - ERF_NumericalSchemes.H: godunov_norm_grad_phi()
  - ERF_FarsiteEllipse.H:   advance_farsite_one_step()
  - ERF_TerrainSlope.cpp:   compute_fire_surface_height()

Run: python3 test_terrain_projection.py
"""

import sys
import math


def surface_grad_norm(phi_x, phi_y, slope_x=0.0, slope_y=0.0):
    """
    Godunov gradient norm projected onto the terrain surface.

    Formula from ERF_NumericalSchemes.H, godunov_norm_grad_phi():
      |grad(phi)|_surface = |grad(phi)|^2 / sqrt(|grad(phi)|^2 + (s.grad(phi))^2)

    which equals |grad(phi)| / sqrt(1 + (s.n)^2) for the unit normal
    n = grad(phi)/|grad(phi)|.

    Args:
        phi_x, phi_y:     Level-set gradient components [1/m]
        slope_x, slope_y: Terrain slopes dz/dx, dz/dy [-]

    Returns:
        Projected gradient norm [1/m]
    """
    grad_sq = phi_x * phi_x + phi_y * phi_y
    s_dot_grad = phi_x * slope_x + phi_y * slope_y
    denom = math.sqrt(max(grad_sq + s_dot_grad * s_dot_grad, 1.0e-14))
    return grad_sq / denom


def projection_factor(nx, ny, slope_x, slope_y):
    """
    Analytic map projection factor cos(theta) for spread direction (nx, ny).

    tan(theta) = grad(z) . n, so cos(theta) = 1 / sqrt(1 + (s.n)^2).
    """
    s_n = slope_x * nx + slope_y * ny
    return 1.0 / math.sqrt(1.0 + s_n * s_n)


def farsite_map_displacement(ds, nx, ny, slope_x=0.0, slope_y=0.0):
    """
    Map-view displacement of the FARSITE front for a surface step ds.

    Formula from ERF_FarsiteEllipse.H, advance_farsite_one_step():
      ds_map = ds / sqrt(1 + (s.n)^2)
    """
    s_n = slope_x * nx + slope_y * ny
    return ds / math.sqrt(1.0 + s_n * s_n)


def atm_ground_elevation(z_nodes):
    """
    Ground elevation of an atmospheric cell.

    Formula from ERF_TerrainSlope.cpp, compute_fire_surface_height(): the mean
    of the four k = 0 nodes of the cell.

    Args:
        z_nodes: (z_ll, z_lr, z_ul, z_ur) nodal heights [m]
    """
    return 0.25 * sum(z_nodes)


def check(failures, name, condition, detail=""):
    status = "✓" if condition else "✗"
    print(f"{status} {name}")
    if not condition:
        if detail:
            print(f"    {detail}")
        failures.append(name)


def test_flat_terrain_is_identity():
    """Zero slope must leave the gradient norm untouched."""
    failures = []
    worst = 0.0
    for phi_x, phi_y in [(1.0, 0.0), (0.0, 1.0), (0.6, 0.8), (-0.3, 0.4)]:
        plain = math.sqrt(phi_x * phi_x + phi_y * phi_y)
        worst = max(worst, abs(surface_grad_norm(phi_x, phi_y) - plain))
    check(failures, "Flat terrain leaves |grad(phi)| unchanged", worst < 1.0e-12,
          f"max deviation {worst}")
    return failures


def test_upslope_projection():
    """Head-on slope gives exactly cos(theta)."""
    failures = []
    for deg in (15.0, 30.0, 45.0, 60.0):
        s = math.tan(math.radians(deg))
        # Front moving in +x, slope in +x
        ratio = surface_grad_norm(1.0, 0.0, s, 0.0) / surface_grad_norm(1.0, 0.0)
        expected = math.cos(math.radians(deg))
        check(failures, f"Upslope projection at {deg:.0f} deg equals cos(theta)",
              abs(ratio - expected) < 1.0e-12,
              f"got {ratio}, expected {expected}")
    return failures


def test_cross_slope_is_unprojected():
    """A front running along the contour is already horizontal."""
    failures = []
    s = math.tan(math.radians(60.0))
    # Front moving in +y, slope entirely in +x
    ratio = surface_grad_norm(0.0, 1.0, s, 0.0) / surface_grad_norm(0.0, 1.0)
    check(failures, "Cross-slope spread is not projected", abs(ratio - 1.0) < 1.0e-12,
          f"got {ratio}, expected 1.0")
    return failures


def test_oblique_matches_directional_formula():
    """Oblique spread follows the directional factor, not a component-wise one."""
    failures = []
    worst = 0.0
    inv_sqrt2 = 1.0 / math.sqrt(2.0)
    for deg in (30.0, 60.0):
        for nx, ny in [(inv_sqrt2, inv_sqrt2), (0.6, 0.8), (-inv_sqrt2, inv_sqrt2)]:
            s = math.tan(math.radians(deg))
            ratio = surface_grad_norm(nx, ny, s, 0.0) / surface_grad_norm(nx, ny)
            worst = max(worst, abs(ratio - projection_factor(nx, ny, s, 0.0)))
    check(failures, "Oblique spread matches the directional projection factor",
          worst < 1.0e-12, f"max deviation {worst}")
    return failures


def test_projection_never_amplifies():
    """The factor is a cosine: it can only slow map-view spread."""
    failures = []
    ok = True
    for sx in (0.0, 0.2, 0.5, 1.0, 2.0):
        for sy in (0.0, 0.3, 1.5):
            for nx, ny in [(1.0, 0.0), (0.0, 1.0), (0.6, 0.8)]:
                ratio = surface_grad_norm(nx, ny, sx, sy) / surface_grad_norm(nx, ny)
                if not (0.0 < ratio <= 1.0 + 1.0e-12):
                    ok = False
    check(failures, "Projection factor stays in (0, 1]", ok)
    return failures


def test_farsite_displacement_projection():
    """FARSITE map displacement carries the same factor as the level set."""
    failures = []
    worst = 0.0
    for deg in (0.0, 30.0, 60.0):
        s = math.tan(math.radians(deg))
        for nx, ny in [(1.0, 0.0), (0.0, 1.0), (0.6, 0.8)]:
            ds = 7.0
            ds_map = farsite_map_displacement(ds, nx, ny, s, 0.0)
            level_set = ds * (surface_grad_norm(nx, ny, s, 0.0)
                              / surface_grad_norm(nx, ny))
            worst = max(worst, abs(ds_map - level_set))
    check(failures, "FARSITE and level-set projections agree", worst < 1.0e-12,
          f"max deviation {worst}")
    return failures


def test_ground_elevation_datum():
    """Wind extraction measures from the ground, not the first cell centre."""
    failures = []

    # Planar terrain z = 0.1 x over an atmospheric cell spanning x in [100, 140]
    def z_of(x):
        return 0.1 * x

    dx = 40.0
    x_lo, x_hi = 100.0, 100.0 + dx
    nodes = (z_of(x_lo), z_of(x_hi), z_of(x_lo), z_of(x_hi))
    ground = atm_ground_elevation(nodes)
    exact = z_of(0.5 * (x_lo + x_hi))
    check(failures, "Four-node average recovers the terrain at the cell centre",
          abs(ground - exact) < 1.0e-12, f"got {ground}, expected {exact}")

    # The superseded datum was the first cell centre, half a cell above ground
    dz = 10.0
    z_ref = 6.1
    old_target = (ground + 0.5 * dz) + z_ref
    new_target = ground + z_ref
    check(failures, "Old datum overshot the requested height by dz/2",
          abs((old_target - new_target) - 0.5 * dz) < 1.0e-12,
          f"offset {old_target - new_target}, expected {0.5 * dz}")
    return failures


def main():
    print("=" * 70)
    print("Terrain Projection Unit Tests")
    print("=" * 70)

    all_failures = []

    all_failures.extend(test_flat_terrain_is_identity())
    all_failures.extend(test_upslope_projection())
    all_failures.extend(test_cross_slope_is_unprojected())
    all_failures.extend(test_oblique_matches_directional_formula())
    all_failures.extend(test_projection_never_amplifies())
    all_failures.extend(test_farsite_displacement_projection())
    all_failures.extend(test_ground_elevation_datum())

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
