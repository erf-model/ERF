#!/usr/bin/env python3
"""Check sampled multiplicities of super-droplets with more than one varying component.

Usage: check_sampled_multiplicity.py INPUTS PLOTFILE

INPUTS is the deck, PLOTFILE holds the super-droplets at t = 0. Mode 0 (soil mass > 0) carries
two log-normal aerosols, a soil core and its NH4HSO4 coating; mode 1 carries an exponential
water mass on a constant aerosol. Each component is sampled uniformly in the logarithm and
weighted by its density, so the multiplicity-weighted mean and standard deviation of the
logarithm of each component match those of its truncated distribution, the two aerosols of
mode 0 are uncorrelated, and the multiplicities of each mode add up to its number density.
The cell-by-cell rescaling of the weights widens the sampled distributions by about 2% at the
128 super-droplets per cell of the deck, well inside the tolerances. With the weights of the two
aerosols summed instead of multiplied, a third of the cores and two thirds of the coatings are
log-uniform; with the water weight left out, the water mass is log-uniform.
"""
import math
import struct
import sys
from pathlib import Path

MEAN_TOL = 0.05   # absolute, in the logarithm
STD_TOL = 0.05    # relative
CORR_TOL = 0.05   # absolute


def num_tol(plotfile, name="super_droplets_moisture"):
    """Relative tolerance of the summed multiplicities for the precision of the particle data."""
    with (Path(plotfile) / name / "Header").open() as stream:
        return 1.0e-4 if stream.readline().strip().endswith("_single") else 1.0e-6


def read_inputs(path):
    """key = values pairs of an inputs file, comments aside."""
    values = {}
    with open(path) as stream:
        for line in stream:
            line = line.split("#")[0].strip()
            if "=" in line:
                key, value = line.split("=", 1)
                values[key.strip()] = [v.strip('"') for v in value.split()]
    return values


def particles(plotfile, fields, name="super_droplets_moisture"):
    """The named real attributes of every super-droplet in a plotfile."""
    directory = Path(plotfile) / name
    with (directory / "Header").open() as stream:
        version = stream.readline().strip()
        assert int(stream.readline()) == 3
        nr = int(stream.readline())
        names = [stream.readline().strip() for _ in range(nr)]
        ni = int(stream.readline())
        for _ in range(ni):
            stream.readline()
        assert int(stream.readline()) == 1
        total = int(stream.readline())
        stream.readline()
        finest = int(stream.readline())
        grids = [int(stream.readline()) for _ in range(finest + 1)]
        chunks = [[tuple(map(int, stream.readline().split())) for _ in range(n)] for n in grids]
    code = "d" if version.endswith("_double") else "f"
    width = struct.calcsize(code)
    index = [3 + names.index(f) for f in fields]
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
                result.append(tuple(values[i] for i in index))
    assert len(result) == total, (plotfile, len(result), total)
    return result


def moments(density, lo, hi, n=20000):
    """Mean and standard deviation of x for a density in x truncated to [lo, hi]."""
    h = (hi - lo) / n
    s0 = s1 = s2 = 0.0
    for i in range(n + 1):
        x = lo + i * h
        w = density(x) * (0.5 if i in (0, n) else 1.0)
        s0, s1, s2 = s0 + w, s1 + w * x, s2 + w * x * x
    mean = s1 / s0
    return mean, math.sqrt(max(s2 / s0 - mean * mean, 0.0))


def weighted(xs, ws):
    total = sum(ws)
    mean = sum(w * x for x, w in zip(xs, ws)) / total
    return mean, math.sqrt(max(sum(w * (x - mean) ** 2 for x, w in zip(xs, ws)) / total, 0.0))


def lognormal(rg, sg):
    s = math.log(sg)
    return lambda x: math.exp(-(x - math.log(rg)) ** 2 / (2.0 * s * s))


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        return 2
    inp = read_inputs(argv[1])
    key = "super_droplets_moisture."
    lo = [float(v) for v in inp["geometry.prob_lo"]]
    hi = [float(v) for v in inp["geometry.prob_hi"]]
    volume = (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2])
    # densities of soil and NH4HSO4 (setProperties_soil, setProperties_NH4HSO4)
    rho_soil, rho_bisulfate = 1220.0, 1780.0

    sds = particles(argv[2], ["multiplicity", "aerosol_mass_soil", "aerosol_mass_NH4HSO4",
                              "species_mass_H2O"])
    dust = [p for p in sds if p[1] > 0.0]
    water = [p for p in sds if p[1] == 0.0]
    if not dust or not water:
        print("FAIL: %d super-droplets in mode 0 and %d in mode 1" % (len(dust), len(water)))
        return 1
    radius = lambda m, rho: (3.0 * m / (4.0 * math.pi * rho)) ** (1.0 / 3.0)

    checks = []
    xs_core = [math.log(radius(p[1], rho_soil)) for p in dust]
    xs_coat = [math.log(radius(p[2], rho_bisulfate)) for p in dust]
    xs_water = [math.log(p[3]) for p in water]
    m0 = key + "0.initial_aerosol_"
    value = lambda k: float(inp[k][0])
    core = lognormal(value(m0 + "mean_radius_soil"), value(m0 + "geomstd_radius_soil"))
    coat = lognormal(value(m0 + "mean_radius_NH4HSO4"), value(m0 + "geomstd_radius_NH4HSO4"))
    m1 = key + "1.initial_species_"
    mean_m, min_m = value(m1 + "mean_mass_H2O"), value(m1 + "min_mass_H2O")
    delta = mean_m - min_m
    expo = lambda x: math.exp(x) / delta * math.exp(-(math.exp(x) - min_m) / delta)
    # sampling ranges of the deck, in the logarithm
    bounds = [(m0 + "min_radius_soil", m0 + "max_radius_soil"),
              (m0 + "min_radius_NH4HSO4", m0 + "max_radius_NH4HSO4"),
              (m1 + "min_mass_H2O", m1 + "max_mass_H2O")]
    for (label, xs, sds_, dens), (klo, khi) in zip([("dust core ln r", xs_core, dust, core),
                                                  ("dust coating ln r", xs_coat, dust, coat),
                                                  ("water ln m", xs_water, water, expo)], bounds):
        target = moments(dens, math.log(value(klo)), math.log(value(khi)))
        got = weighted(xs, [p[0] for p in sds_])
        checks.append((label + " mean", got[0], target[0], abs(got[0] - target[0]) <= MEAN_TOL))
        checks.append((label + " std", got[1], target[1], abs(got[1] / target[1] - 1.0) <= STD_TOL))
    w = [p[0] for p in dust]
    (mc, sc), (mt, st) = weighted(xs_core, w), weighted(xs_coat, w)
    corr = sum(wi * (a - mc) * (b - mt) for wi, a, b in zip(w, xs_core, xs_coat)) / (sum(w) * sc * st)
    checks.append(("core-coating correlation", corr, 0.0, abs(corr) <= CORR_TOL))
    for label, sds_, mode in [("dust number density", dust, 0), ("water number density", water, 1)]:
        n0 = float(inp[key + "%d.initial_number_density" % mode][0])
        got = sum(p[0] for p in sds_) / volume
        checks.append((label, got, n0, abs(got / n0 - 1.0) <= num_tol(argv[2])))

    failed = 0
    print("%-24s %14s %14s" % ("quantity", "particles", "expected"))
    for label, got, target, ok in checks:
        print("%-24s %14.6g %14.6g  %s" % (label, got, target, "ok" if ok else "FAIL"))
        failed += 0 if ok else 1
    print("%d super-droplets in mode 0, %d in mode 1" % (len(dust), len(water)))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
