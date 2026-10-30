#!/usr/bin/env python3
"""Check the fixed east-sun hill case and optionally plot the shadow and fluxes."""
import argparse
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def read(path):
    with Dataset(path) as ds:
        names = ("z_surf", "topographic_shadow", "terrain_direct_factor", "terrain_incidence",
                 "unshaded_sw_flux_dn", "terrain_sw_flux_dn", "sw_flux_dn", "cos_zenith_angle")
        fields = {name: np.array(ds[name][:]).squeeze() for name in names}
        dx, dy = np.array(ds["CellSize"][:])[:2]
    for name, value in fields.items():
        if value.ndim != 2 or not np.all(np.isfinite(value)) or np.any(value == -999):
            raise AssertionError(f"missing or invalid {name} in {path}")
    return fields, dx, dy


def check(fields, dx):
    h = fields["z_surf"]
    mask = fields["topographic_shadow"]
    assert np.all((mask == 0) | (mask == 1)), "shadow mask is not binary"
    assert 0 < mask.sum() < mask.size, "expected both shadow and sunlit ground"
    # Independent horizon oracle: compute the maximum elevation angle among
    # all mass points eastward within 3 km, including periodic copies. With the
    # cardinal sun, the maximum of each line segment's slope is at an endpoint.
    horizon = np.full_like(h, -np.inf)
    for distance in range(1, int(3000/dx) + 1):
        horizon = np.maximum(horizon, (np.roll(h, -distance, axis=1) - h)/(distance*dx))
    expected = horizon > np.tan(np.deg2rad(10)) + 1.e-8
    np.testing.assert_array_equal(mask, expected, err_msg="terrain horizon disagrees with ERF mask")
    yy, xx = np.where(mask == 1)
    assert np.all(xx < h.shape[1]//2), "east sun should cast a westward shadow"
    factor = fields["terrain_direct_factor"]
    assert np.all(factor[expected] == 0), "direct beam survived in a terrain shadow"
    np.testing.assert_allclose(fields["cos_zenith_angle"], np.sin(np.deg2rad(10)), rtol=1.e-6)
    np.testing.assert_allclose(fields["terrain_sw_flux_dn"], fields["sw_flux_dn"], atol=1.e-8,
                               err_msg="Noah-MP did not receive the corrected flux")
    raw, corrected = fields["unshaded_sw_flux_dn"], fields["terrain_sw_flux_dn"]
    assert np.all(raw > 0) and np.all(corrected >= 0)
    assert np.all(corrected[expected] < raw[expected]), "shadows did not reduce shortwave"
    assert np.all(corrected[expected] > 0), "diffuse sunlight should remain in shadows"
    assert np.any(corrected > raw), "expected slope enhancement on the sun-facing hill"
    print(f"PASS: {int(mask.sum())}/{mask.size} cells shadowed; Noah-MP receives corrected SW. "
          f"Shadowed SW range {corrected[expected].min():.3f}..{corrected[expected].max():.3f} W/m2.")


def plot(fields, dx, dy, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    h = fields["z_surf"]
    ny, nx = h.shape
    x, y = (np.arange(nx)+0.5)*dx/1000, (np.arange(ny)+0.5)*dy/1000
    extent = (0, nx*dx/1000, 0, ny*dy/1000)
    fig, axes = plt.subplots(2, 2, figsize=(12, 7), constrained_layout=True)
    datasets = [("Terrain height", h, "terrain", "m"),
                ("Terrain shadow (1 = blocked)", fields["topographic_shadow"], "Greys", ""),
                ("Uncorrected downwelling shortwave", fields["unshaded_sw_flux_dn"], "magma", "W/m²"),
                ("Downwelling shortwave supplied to Noah-MP", fields["terrain_sw_flux_dn"], "magma", "W/m²")]
    vmax = max(fields["terrain_sw_flux_dn"].max(), fields["unshaded_sw_flux_dn"].max())
    for ax, (title, values, cmap, units) in zip(axes.flat, datasets):
        kw = {"vmin": 0, "vmax": vmax} if units == "W/m²" else {}
        im = ax.imshow(values, origin="lower", extent=extent, cmap=cmap, **kw)
        ax.contour(x, y, h, levels=[25, 50, 100, 150, 200], colors="cyan", linewidths=0.6, alpha=0.7)
        ax.set(title=title, xlabel="x (km; east →)", ylabel="y (km)")
        fig.colorbar(im, ax=ax, label=units, shrink=0.8)
    fig.suptitle("Witch of Agnesi hill · sun 10° above eastern horizon", fontsize=15)
    fig.savefig(output, dpi=160)
    print(f"Wrote {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plotfile", type=Path, help="surface00001_d01.nc (or a later step)")
    parser.add_argument("--held", type=Path, help="a second output from the same radiation interval")
    parser.add_argument("--plot", type=Path, help="write a PNG visualization")
    args = parser.parse_args()
    fields, dx, dy = read(args.plotfile)
    check(fields, dx)
    if args.held:
        held, _, _ = read(args.held)
        for name in fields:
            np.testing.assert_array_equal(fields[name], held[name], err_msg=f"{name} changed during held radiation interval")
        print("PASS: surface radiation and shadow remain unchanged between radiation updates.")
    if args.plot:
        plot(fields, dx, dy, args.plot)


if __name__ == "__main__":
    main()
