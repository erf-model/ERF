"""Check local DSDs against two immobile, monodisperse populations of known multiplicity."""
import argparse
import math
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile


def check(directory, name, times, edges, bins):
    lines = (directory / f"{name}.dsd").read_text().splitlines()
    actual_edges = [float(v) for v in next(s for s in lines if s.startswith("# radius_edges_m:")).split()[2:]]
    rows = [[float(v) for v in s.split()] for s in lines if not s.startswith("#")]
    for actual, expected in [(actual_edges, edges), ([r[0] for r in rows], times), *[(r[1:], bins) for r in rows]]:
        assert len(actual) == len(expected) and all(math.isclose(a, e, rel_tol=2e-6) for a, e in zip(actual, expected)), (
            name, actual, expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--mpi-exec", default="")
    parser.add_argument("--mpi-numproc-flag", default="-n")
    parser.add_argument("--mpi-preflags", default="")
    parser.add_argument("--ranks", default="1")
    args = parser.parse_args()
    launch = shlex.split(args.mpi_exec)
    if launch:
        launch += [args.mpi_numproc_flag, args.ranks, *shlex.split(args.mpi_preflags)]
    root = Path(tempfile.mkdtemp(dir=Path.cwd(), prefix="probe-"))
    print(f"Probe regression directory: {root}", flush=True)  # kept on failure
    prefix = "super_droplets_moisture."
    config = ["num_initializations=2", "1.initial_species_mean_mass_H2O=8e-12",
              "1.initial_number_density=2e7", "dsd_probe_names=Cell Whole Small Empty", "dsd_probe_interval=1",
              "Cell.cell=2 2 3", "Cell.r_edges=0 7e-6 10e-6 15e-6", "Cell.write_mass=true",
              "Whole.box_lo=0 0 0", "Whole.box_hi=8 4 16", "Whole.rmin=0", "Whole.rmax=15e-6", "Whole.nbins=3",
              "Small.center=2.5 2.5 3.5", "Small.side_lengths=0.5", "Small.rmin=1e-6", "Small.rmax=16e-6",
              "Small.nbins=4", "Small.spacing=log", "Empty.box_lo=0 0 0", "Empty.box_hi=0.5 0.5 0.5",
              "Empty.r_edges=0 1e-4"]

    def run(directory, label, options, fail=None):
        directory.mkdir(exist_ok=True)
        log = directory / f"{label}.log"
        command = [*launch, args.exe, str(args.input), f"erf.input_sounding_file={args.input.with_name('input_sounding')}",
                   "amr.max_grid_size_x=4", "amr.blocking_factor=2", "erf.check_int=1", "erf.plot_int_1=-1",
                   *[prefix + s for s in config], *options]
        with log.open("w") as f:
            result = subprocess.run(command, cwd=directory, stdout=f, stderr=subprocess.STDOUT, timeout=180, check=False)
        if fail:
            assert result.returncode and fail in log.read_text(), log
        elif result.returncode:
            raise RuntimeError(f"{log}:\n{log.read_text()[-3000:]}")

    cell_edges = [0, 7e-6, 10e-6, 15e-6]
    cell_bins = [1e7, 0, 2e7, 1e-5, 0, 1.6e-4]
    flat = root / "flat"
    run(flat, "straight", ["max_step=3"])
    times = [0, .02, .04, .06]
    check(flat, "Cell", times, cell_edges, cell_bins)
    check(flat, "Whole", times, [0, 5e-6, 10e-6, 15e-6], [0, 1e7, 2e7])
    # One cell's particles in 1/8 m³: this must NOT be a cell-averaged concentration.
    check(flat, "Small", times, [1e-6, 2e-6, 4e-6, 8e-6, 16e-6], [0, 0, 8e7, 16e7])
    check(flat, "Empty", times, [0, 1e-4], [0])  # particle on the excluded upper face
    assert not list(flat.glob("*_g_lnR_*.txt")), "probe should not enable domain output"
    # Rewind to an off-cadence step. Old rows past the checkpoint must not survive.
    restart = ["erf.restart=chk00001", "max_step=2", prefix + "dsd_probe_interval=2", "erf.check_int=-1"]
    run(flat, "rewind", restart)
    check(flat, "Cell", [0, .04], cell_edges, cell_bins)
    run(flat, "on_cadence", [*restart, "erf.restart=chk00002", "max_step=3"])
    check(flat, "Cell", [0, .04], cell_edges, cell_bins)
    run(flat, "mismatch", [*restart, prefix + "Cell.r_edges=0 8e-6 10e-6 15e-6"], "incompatible probe header")

    # A refined level owns the particles at Cell; the box spans coarse and fine cells.
    refined = root / "refined"
    run(refined, "refined", ["max_step=0", "amr.max_level=1", "amr.ref_ratio=2", "erf.refinement_indicators=box1",
                             "erf.box1.max_level=1", "erf.box1.in_box_lo=0 0 0", "erf.box1.in_box_hi=8 4 8"])
    check(refined, "Cell", [0], cell_edges, cell_bins)
    check(refined, "Whole", [0], [0, 5e-6, 10e-6, 15e-6], [0, 1e7, 2e7])

    # First eight layers are 0.5 m deep, the next eight are 1.5 m. The particle
    # at computational z=3.5 is at physical z=1.75, and represents half as many drops.
    stretched = root / "stretched"
    levels = [i*.5 for i in range(9)] + [4+i*1.5 for i in range(1, 9)]
    run(stretched, "stretched", ["max_step=0", "erf.terrain_z_levels=" + " ".join(map(str, levels)),
                                prefix + "Small.center=2.5 2.5 1.75", prefix + "Small.side_lengths=0.5 0.5 0.25"])
    check(stretched, "Cell", [0], cell_edges, cell_bins)
    check(stretched, "Small", [0], [1e-6, 2e-6, 4e-6, 8e-6, 16e-6], [0, 0, 8e7, 16e7])
    run(root / "bad", "bad_edges", ["max_step=0", prefix + "Cell.r_edges=0 1e-6 1e-6"], "strictly increasing")
    shutil.rmtree(root)
    print("PASS: number/mass bins, spatial selection, AMR, stretched z, cadence, restart and rejected inputs")


if __name__ == "__main__":
    main()
