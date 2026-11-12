#!/usr/bin/env python3
"""Check Golovin coalescence of super-droplets on a stretched mesh, layer by layer.

Usage: check_stretched_coalescence.py INPUTS PLOTFILE_0 PLOTFILE_T

INPUTS is the deck (erf.terrain_z_levels, geometry.prob_lo/prob_hi, amr.n_cell,
stop_time and super_droplets_moisture.initial_number_density are read from it), and
PLOTFILE_0 and PLOTFILE_T hold the super-droplets at the start and at stop_time. Without
advection every layer is a closed box, so from the particle data alone:

- every layer starts at the specified number density;
- every layer keeps its liquid volume fraction phi (sum of multiplicity times droplet
  volume over the layer volume): coalescence conserves water;
- the number density decays as exp(-b phi t), b = 1500 s^-1 for the Golovin kernel. The
  rate is checked over the thin layers together and over the thick layers together.

With the computational cell volume in place of the physical one, the thin (50 m) layers
start at twice the specified density and decay at half the rate, the thick (150 m) layers
at two thirds of it and 1.5 times the rate.
"""
import bisect
import math
import struct
import sys
from pathlib import Path

GOLOVIN_B = 1.5e3
RATE_TOL = 0.1


def tolerances(plotfile, name="super_droplets_moisture"):
    """Initial-density and water tolerances for the precision of the particle data."""
    with (Path(plotfile) / name / "Header").open() as stream:
        single = stream.readline().strip().endswith("_single")
    return (1.0e-5, 1.0e-5) if single else (1.0e-9, 1.0e-10)


def read_inputs(path):
    """key = values pairs of an inputs file, comments and continuation lines aside."""
    values = {}
    with open(path) as stream:
        for line in stream:
            line = line.split("#")[0].strip()
            if "=" in line:
                key, value = line.split("=", 1)
                values[key.strip()] = value.split()
    return values


def particles(plotfile, name="super_droplets_moisture"):
    """(z, multiplicity, radius) of every super-droplet in a plotfile."""
    directory = Path(plotfile) / name
    with (directory / "Header").open() as stream:
        version = stream.readline().strip()
        assert int(stream.readline()) == 3
        nr = int(stream.readline())
        names = [stream.readline().strip() for _ in range(nr)]
        ni = int(stream.readline())
        for _ in range(ni):
            stream.readline()
        assert int(stream.readline()) == 1  # integer id/cpu and attributes are present
        total = int(stream.readline())
        stream.readline()  # next particle id
        finest = int(stream.readline())
        grids = [int(stream.readline()) for _ in range(finest + 1)]
        chunks = [[tuple(map(int, stream.readline().split())) for _ in range(n)]
                  for n in grids]
    code = "d" if version.endswith("_double") else "f"
    width = struct.calcsize(code)
    xi = 3 + names.index("multiplicity")
    radius = 3 + names.index("radius")
    result = []
    for level, entries in enumerate(chunks):
        for which, count, offset in entries:
            if count == 0:
                continue
            with (directory / ("Level_%d" % level) / ("DATA_%05d" % which)).open("rb") as stream:
                stream.seek(offset + 4 * count * (2 + ni))
                data = struct.unpack("=" + code * count * (3 + nr),
                                     stream.read(width * count * (3 + nr)))
            for p in range(count):
                values = data[p * (3 + nr):(p + 1) * (3 + nr)]
                result.append((values[2], values[xi], values[radius]))
    assert len(result) == total, (plotfile, len(result), total)
    return result


def layers(plotfile, z_levels, area):
    """Number density and liquid volume fraction of each layer."""
    number = [0.0] * (len(z_levels) - 1)
    volume = [0.0] * (len(z_levels) - 1)
    for z, mult, radius in particles(plotfile):
        k = bisect.bisect_right(z_levels, z) - 1
        assert 0 <= k < len(number), ("particle outside the column", z)
        number[k] += mult
        volume[k] += mult * 4.0 / 3.0 * math.pi * radius**3
    dv = [area * (z_levels[k + 1] - z_levels[k]) for k in range(len(number))]
    return ([n / v for n, v in zip(number, dv)], [w / v for w, v in zip(volume, dv)], dv)


def main(argv):
    if len(argv) != 4:
        print(__doc__)
        return 2
    inputs = read_inputs(argv[1])
    z_levels = [float(z) for z in inputs["erf.terrain_z_levels"]]
    lo = [float(x) for x in inputs["geometry.prob_lo"]]
    hi = [float(x) for x in inputs["geometry.prob_hi"]]
    nz = int(inputs["amr.n_cell"][2])
    t_end = float(inputs["stop_time"][0])
    n_init = float(inputs["super_droplets_moisture.initial_number_density"][0])
    area = (hi[0] - lo[0]) * (hi[1] - lo[1])
    dz_comp = (hi[2] - lo[2]) / nz

    init_tol, phi_tol = tolerances(argv[2])
    n_0, phi_0, dv = layers(argv[2], z_levels, area)
    n_t, phi_t, _ = layers(argv[3], z_levels, area)

    failed = 0
    groups = {"thin": [0.0, 0.0, 0.0, 0.0], "thick": [0.0, 0.0, 0.0, 0.0]}
    print("layer  dz      n0/n_init     phi_t/phi_0-1   rate/(b phi)")
    for k in range(len(n_0)):
        dz = z_levels[k + 1] - z_levels[k]
        init = n_0[k] / n_init
        drift = phi_t[k] / phi_0[k] - 1.0
        rate = -math.log(n_t[k] / n_0[k]) / t_end / (GOLOVIN_B * phi_0[k])
        print("%5d  %6.1f  %.10f  %+.3e      %.4f" % (k, dz, init, drift, rate))
        if abs(init - 1.0) > init_tol:
            print("  FAIL: initial number density differs from the specified one")
            failed += 1
        if abs(drift) > phi_tol:
            print("  FAIL: liquid volume fraction not conserved")
            failed += 1
        group = groups["thin" if dz < dz_comp else "thick"]
        group[0] += n_0[k] * dv[k]
        group[1] += n_t[k] * dv[k]
        group[2] += phi_0[k] * dv[k]
        group[3] += dv[k]
    for name, (count_0, count_t, water, vol) in groups.items():
        if vol == 0.0:
            continue
        rate = -math.log(count_t / count_0) / t_end / (GOLOVIN_B * water / vol)
        print("%s layers: number decay rate / (b phi) = %.4f" % (name, rate))
        if abs(rate - 1.0) > RATE_TOL:
            print("  FAIL: decay rate differs from b phi by more than %g" % RATE_TOL)
            failed += 1
    print("PASS" if failed == 0 else "FAILED (%d checks)" % failed)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
