"""Exercise authoritative vapor diffusion and particle condensate ownership."""

import argparse
import math
import shutil
import shlex
import subprocess
import sys
import struct
from pathlib import Path


def particle_budget(plotfile, header):
    """Independently integrate binary particle attributes, not mesh deposition."""
    directory = plotfile / "super_droplets_moisture"
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
    species = [3 + names.index("species_mass_" + name) for name in ("water", "agua")]
    masses = [0.0, 0.0]
    counts = [0] * (finest + 1)
    for level, entries in enumerate(chunks):
        for which, count, offset in entries:
            if count == 0:
                continue
            with (directory / ("Level_%d" % level) / ("DATA_%05d" % which)).open("rb") as stream:
                stream.seek(offset)
                integers = struct.unpack("=" + "i" * count * (2 + ni),
                                         stream.read(4 * count * (2 + ni)))
                data = struct.unpack("=" + code * count * (3 + nr),
                                     stream.read(width * count * (3 + nr)))
            for p in range(count):
                values = data[p * (3 + nr):(p + 1) * (3 + nr)]
                assert integers[p * (2 + ni) + 2] > 0, "inactive stored particle"
                # Temperature and cold-only attributes are not initialized in
                # this warm fixture. Check positions and represented reservoirs.
                assert all(math.isfinite(values[c]) for c in (0, 1, 2, xi, *species))
                assert values[xi] >= 1
                fine_index = tuple(math.floor(values[d] / (header["dx"][d] / 2))
                                   for d in range(3))
                owner = int(fine_index in header["fine_cells"])
                assert level == owner, ("particle is not on its finest level", level, values[:3])
                for s in range(2):
                    assert values[species[s]] >= 0
                    masses[s] += values[xi] * values[species[s]]
                counts[level] += 1
    assert sum(counts) == total
    return masses, counts


def composite(plotfile):
    """Read this flat test's one/two-level, ratio-two composite cell data."""
    header = _read_header(str(plotfile))
    assert header["finest"] in (0, 1)
    names = header["names"]
    fields = ["density", "rhoQ7", "rhoQ8", "rhoQ9", "rhoQ10", "Kmv"]
    fields += [
        "super_droplets_moisture_species_mass_density_" + species
        for species in ("water", "agua")
    ]
    components = [names.index(field) for field in fields]
    levels = []
    for level in range(header["finest"] + 1):
        directory = plotfile / ("Level_%d" % level)
        ncomp, _, boxes, fabs = _read_cell_h(str(directory))
        cells = {}
        for (lo, hi), (filename, offset) in zip(boxes, fabs):
            flo, fhi, data = _read_fab(str(directory / filename), offset, ncomp)
            nx, ny, nz = [fhi[d] - flo[d] + 1 for d in range(3)]
            for k in range(lo[2], hi[2] + 1):
                for j in range(lo[1], hi[1] + 1):
                    for i in range(lo[0], hi[0] + 1):
                        index = ((k - flo[2]) * ny + j - flo[1]) * nx + i - flo[0]
                        cells[i, j, k] = [
                            data[c * nx * ny * nz + index] for c in components
                        ]
        levels.append(cells)
    result = []
    header["fine_cells"] = frozenset(levels[1]) if len(levels) > 1 else frozenset()
    for level, cells in enumerate(levels):
        covered = (
            {tuple(n // 2 for n in index) for index in levels[level + 1]}
            if level < header["finest"] else set()
        )
        volume = math.prod(header["dx"]) / 8**level
        result.extend((volume, values) for index, values in cells.items()
                      if index not in covered)
    return header, result


def metrics(plotfile, tolerance):
    header, cells = composite(plotfile)
    mass = [0.0, 0.0]
    vapor_square = [0.0, 0.0]
    condensate = [0.0, 0.0]
    for volume, values in cells:
        assert all(math.isfinite(value) for value in values), plotfile
        rho = values[0]
        assert rho > 0
        for species in range(2):
            vapor = values[1 + 2 * species]
            projected = values[2 + 2 * species]
            particle = values[6 + species]
            assert vapor >= -tolerance * 0.005, (plotfile, vapor)
            assert projected >= 0 and particle >= 0, (plotfile, projected, particle)
            assert abs(projected - particle) <= tolerance * max(1.e-8, particle), (
                "condensate is not the particle projection", plotfile, projected, particle
            )
            mass[species] += volume * (vapor + particle)
            condensate[species] += volume * particle
            vapor_square[species] += volume * vapor * vapor / rho
    return header, mass, vapor_square, condensate, max(v[5] for _, v in cells)


def run(args, label, options, injected_volume=0):
    directory = Path(label)
    directory.mkdir(exist_ok=True)
    # Remove only this checker's named outputs, so reruns cannot reuse old plots.
    for plot in directory.glob("plt[0-9][0-9][0-9][0-9][0-9]"):
        if plot.is_dir():
            shutil.rmtree(plot)
    command = []
    if args.mpi_exec:
        command = shlex.split(args.mpi_exec)
        command += [args.mpi_numproc_flag, str(args.ranks)]
        command += shlex.split(args.mpi_preflags)
    command += [
        args.exe, args.input,
        "erf.input_sounding_file=" + str(Path(args.input).with_name("input_sounding")),
    ] + options
    with (directory / "simulation.log").open("w") as log:
        subprocess.run(command, cwd=directory, stdout=log, stderr=subprocess.STDOUT,
                       check=True, timeout=180)
    plots = sorted(directory.glob("plt[0-9][0-9][0-9][0-9][0-9]"))
    assert len(plots) >= 2, directory
    first = metrics(plots[0], args.tolerance)
    last = metrics(plots[-1], args.tolerance)
    assert abs(first[0]["time"]) < 1.e-12
    assert abs(last[0]["time"] - 0.08) < 1.e-6
    for plot in plots[1:]:
        current = metrics(plot, args.tolerance)
        for s, (initial, measured) in enumerate(zip(first[1], current[1])):
            expected = initial + injected_volume * 1.e6 * current[0]["time"] * (s + 1) * 1.e-12
            assert abs(measured - expected) <= args.tolerance * expected, (
                "species mass is not conserved (including injection)", plot, expected, measured
            )
    if injected_volume:
        for plot in plots:
            current = metrics(plot, args.tolerance)
            particles, counts = particle_budget(plot, current[0])
            for deposited, independent in zip(current[3], particles):
                assert abs(deposited - independent) <= args.tolerance * max(independent, 1.e-20), (
                    "mesh projection disagrees with independent particle mass", plot, deposited, independent
                )
            print(plot, "independent particle mass:", particles, "level counts:", counts, flush=True)
        if args.case == "injection":
            for s in range(2):
                expected = first[3][s] + injected_volume * 1.e6 * last[0]["time"] * (s + 1) * 1.e-12
                assert abs(last[3][s] - expected) <= args.tolerance * expected
        else:
            unexchanged = first[3][0] + injected_volume * 1.e6 * last[0]["time"] * 1.e-12
            assert last[3][0] < 0.9 * unexchanged, "phase change must be active"
    if args.case.startswith("regrid"):
        coverage = [composite(plot)[0]["fine_cells"] for plot in plots]
        assert coverage[0] and coverage[-1], "both refined regions must exist"
        assert any(not cells for cells in coverage), "fine level must disappear"
        assert coverage[0] != coverage[-1], "refinement coverage must actually change"
        for plot in plots[1:]:
            current = metrics(plot, args.tolerance)
            for initial, measured in zip(first[3], current[3]):
                assert abs(measured - initial) <= args.tolerance * initial, (
                    "redistribution changed particle species mass", plot, initial, measured
                )
        print("refined cell counts:", [len(cells) for cells in coverage], flush=True)
    assert last[4] > 1, ("test must activate turbulent diffusion", last[4])
    print(label, "mass:", last[1], "vapor L2:", last[2],
          "particle condensate:", last[3], "max Kmv:", last[4], flush=True)
    return first, last


def main(args):
    if args.case.startswith("injection"):
        for advection in ("false", "true"):
            for splitmerge in ("false", "true"):
                for crossing in (False, True):
                    zlo, zhi = (6, 10) if crossing else (2, 6)
                    options = [
                        "amr.max_level=1", "amr.ref_ratio=2",
                        "erf.refinement_indicators=box1", "erf.box1.max_level=1",
                        "erf.box1.in_box_lo=0 0 0", "erf.box1.in_box_hi=8 4 8",
                        "super_droplets_moisture.include_phase_change=" +
                        str(args.case == "injection_phase").lower(),
                        "super_droplets_moisture.include_advection=" + advection,
                        "super_droplets_moisture.split_merge_amr=" + splitmerge,
                        "super_droplets_moisture.num_injections=1",
                        "super_droplets_moisture.injection.distribution_type=uniform",
                        "super_droplets_moisture.injection.particle_box_lo=2 0 %d" % zlo,
                        "super_droplets_moisture.injection.particle_box_hi=6 4 %d" % zhi,
                        "super_droplets_moisture.injection.rate=1.e6",
                        "super_droplets_moisture.injection.particles_per_cell=1",
                        "super_droplets_moisture.injection.species_mean_mass_H2O=1.e-12",
                        "super_droplets_moisture.injection.species_mean_mass_water=1.e-12",
                        "super_droplets_moisture.injection.species_mean_mass_agua=2.e-12",
                        "super_droplets_moisture.injection.aerosol_distribution_type_NaCl=mass_constant",
                        "super_droplets_moisture.injection.aerosol_mean_mass_NaCl=1.e-19",
                    ]
                    if args.case == "injection_phase":
                        # Latent heating creates thermal gradients: unlike vapor,
                        # their diffusion is explicit and must respect its limit.
                        options.append("erf.Cs=2.0")
                    label = "adv_%s_split_%s_cross_%s" % (advection, splitmerge, crossing)
                    run(args, label, options, injected_volume=64)
        return
    options = []
    if args.case == "amr" or args.case.startswith("regrid"):
        options += [
            "amr.max_level=1", "amr.ref_ratio=2",
            "erf.refinement_indicators=box1", "erf.box1.max_level=1",
            "erf.box1.in_box_lo=0 0 0", "erf.box1.in_box_hi=8 4 8",
        ]
    if args.case.startswith("regrid"):
        options += [
            "erf.regrid_int=1", "erf.refinement_indicators=box1 box2",
            "erf.box1.end_time=0.019", "erf.box2.start_time=0.06",
            "erf.box2.max_level=1", "erf.box2.in_box_lo=0 0 8",
            "erf.box2.in_box_hi=8 4 16",
        ]
        if args.case == "regrid_splitmerge":
            options.append("super_droplets_moisture.split_merge_amr=true")
    if args.case == "phase":
        options.append("super_droplets_moisture.include_phase_change=true")
    initial, implicit = run(args, "implicit", options)
    if args.case == "phase":
        assert implicit[3][0] < 0.9 * initial[3][0], "phase change must be active"
    else:
        for before, after in zip(initial[2], implicit[2]):
            assert after < 0.85 * before, "the vapor slab must actually diffuse"
        for before, after in zip(initial[3], implicit[3]):
            assert abs(before - after) <= args.tolerance * before
    if args.case == "column":
        _, explicit = run(args, "explicit", [
            "erf.implicit_moisture_diffusion=false", "erf.fixed_dt=0.0005",
            "max_step=160", "erf.plot_int_1=160",
        ])
        # Coarse backward-Euler RK stages are deliberately dissipative; the
        # well-resolved explicit trajectory is a consistency check, not bit parity.
        for reference, measured in zip(explicit[2], implicit[2]):
            assert abs(measured - reference) < 0.08 * reference, (
                "implicit and resolved explicit trajectories disagree", reference, measured
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--reader-dir", required=True)
    parser.add_argument("--case", choices=("column", "amr", "phase", "regrid", "regrid_splitmerge",
                                         "injection", "injection_phase"), required=True)
    parser.add_argument("--mpi-exec", default="")
    parser.add_argument("--mpi-numproc-flag", default="-n")
    parser.add_argument("--mpi-preflags", default="")
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--tolerance", type=float, default=2.e-11)
    arguments = parser.parse_args()
    sys.path.insert(0, arguments.reader_dir)
    from erf_plotfile import _read_cell_h, _read_fab, _read_header
    main(arguments)
