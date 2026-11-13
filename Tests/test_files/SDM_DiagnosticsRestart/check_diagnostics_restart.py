#!/usr/bin/env python3
"""Step labels and times of the super-droplet diagnostics, straight and restarted.

Runs a deck to --step-end with checkpoints, plotfiles and diagnostics every --step-chk
steps ("straight"), then again from its checkpoint at --step-chk ("restart"). Checks:

- each run writes the g(ln R) file of every diagnostics step from its own first step on,
  labelled with the step of the state it holds (the restart writes --step-chk at
  start-up, not 00000);
- the restarted run's g(ln R) files equal the straight run's (phase change makes the
  distribution change from step to step);
- every binned distribution plotfile (ERF_ENABLE_ML_UPHYS_DIAGNOSTICS builds) holds the
  time and step of the plotfile written at the same step in the same run.
"""
import argparse
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

NAME = "super_droplets_moisture"


def run(args, directory, options):
    if directory.exists():
        shutil.rmtree(directory)
    directory.mkdir()
    command = []
    if args.mpi_exec:
        command = shlex.split(args.mpi_exec) + [args.mpi_numproc_flag, str(args.ranks)]
        command += shlex.split(args.mpi_preflags)
    command += [args.exe, args.input,
                "erf.input_sounding_file=" + str(Path(args.input).with_name("input_sounding")),
                "max_step=%d" % args.step_end,
                "erf.check_file=chk", "erf.check_int=%d" % args.step_chk,
                "erf.plot_file_1=plt", "erf.plot_int_1=%d" % args.step_chk,
                "%s.diagnostics_interval=%d" % (NAME, args.step_chk),
                "%s.include_phase_change=true" % NAME] + options
    with (directory / "simulation.log").open("w") as log:
        subprocess.run(command, cwd=directory, stdout=log, stderr=subprocess.STDOUT,
                       check=True, timeout=300)


def header_time_and_step(plotfile):
    """Time and level-0 step from a plotfile Header."""
    lines = (plotfile / "Header").read_text().splitlines()
    nvars = int(lines[1])
    return float(lines[3 + nvars]), int(lines[9 + nvars].split()[0])


def g_lnR(directory, step):
    path = directory / ("%s_g_lnR_%05d.txt" % (NAME, step))
    return [[float(x) for x in line.split()] for line in path.read_text().splitlines()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--step-chk", type=int, default=2)
    parser.add_argument("--step-end", type=int, default=4)
    parser.add_argument("--mpi-exec", default="")
    parser.add_argument("--mpi-numproc-flag", default="-n")
    parser.add_argument("--mpi-preflags", default="")
    parser.add_argument("--ranks", type=int, default=1)
    args = parser.parse_args()

    straight = Path("straight").resolve()
    restart = Path("restart").resolve()
    run(args, straight, [])
    run(args, restart, ["erf.restart=%s/chk%05d" % (straight, args.step_chk)])

    failed = []
    steps = range(0, args.step_end + 1, args.step_chk)
    for directory, first in ((straight, 0), (restart, args.step_chk)):
        expected = ["%s_g_lnR_%05d.txt" % (NAME, s) for s in steps if s >= first]
        found = sorted(p.name for p in directory.glob(NAME + "_g_lnR_*.txt"))
        if found != expected:
            failed.append("%s: g(ln R) files %s, expected %s" % (directory.name, found, expected))
        for binned in sorted(directory.glob(NAME + "_binned_dsd_*_[0-9]*")):
            if not binned.is_dir():
                continue
            plotfile = directory / ("plt" + binned.name[-5:])
            got, want = header_time_and_step(binned), header_time_and_step(plotfile)
            print("%s/%s: time %.17g step %d; %s: time %.17g step %d"
                  % (directory.name, binned.name, got[0], got[1], plotfile.name, want[0], want[1]))
            if abs(got[0] - want[0]) > 1.e-6 * max(1.0, abs(want[0])) or got[1] != want[1]:
                failed.append("%s/%s: time and step differ from %s"
                              % (directory.name, binned.name, plotfile.name))

    for step in steps:
        if step < args.step_chk:
            continue
        try:
            a, b = g_lnR(straight, step), g_lnR(restart, step)
        except FileNotFoundError:
            continue
        same = len(a) == len(b) and all(
            len(ra) == len(rb) and all(abs(x - y) <= 1.e-5 * max(abs(x), abs(y), 1.e-300)
                                       for x, y in zip(ra, rb))
            for ra, rb in zip(a, b))
        if not same:
            failed.append("g(ln R) of step %d differs between the straight and restarted runs" % step)

    for message in failed:
        print("FAIL: " + message)
    print("PASS" if not failed else "FAILED (%d checks)" % len(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
