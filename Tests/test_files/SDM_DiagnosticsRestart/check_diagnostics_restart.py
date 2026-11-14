#!/usr/bin/env python3
"""Step labels and times of the super-droplet diagnostics, straight and restarted.

Runs a deck to --step-end with a checkpoint every step and plotfiles and diagnostics
every --interval steps ("straight"), then again from each checkpoint in --restarts.
File names use --digits digits (erf.file_name_digits). Checks:

- each run writes the g(ln R) file of every multiple of --interval from its own first
  step on, labelled with the step of the state it holds (a restart from a multiple
  writes that step at start-up, a restart from another step writes nothing then);
- the restarted runs' g(ln R) files equal the straight run's (phase change makes the
  distribution change from step to step);
- every binned distribution plotfile (ERF_ENABLE_ML_UPHYS_DIAGNOSTICS builds, which
  pass --expect-binned) holds the time and step of the plotfile written at the same
  step in the same run.
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
                "max_step=%d" % args.step_end, "erf.file_name_digits=%d" % args.digits,
                "erf.check_file=chk", "erf.check_int=1",
                "erf.plot_file_1=plt", "erf.plot_int_1=%d" % args.interval,
                "%s.diagnostics_interval=%d" % (NAME, args.interval),
                "%s.include_phase_change=true" % NAME] + options
    with (directory / "simulation.log").open("w") as log:
        subprocess.run(command, cwd=directory, stdout=log, stderr=subprocess.STDOUT,
                       check=True, timeout=300)


def header_time_and_step(plotfile):
    """Time and level-0 step from a plotfile Header."""
    lines = (plotfile / "Header").read_text().splitlines()
    nvars = int(lines[1])
    return float(lines[3 + nvars]), int(lines[9 + nvars].split()[0])


def g_lnR(directory, name):
    path = directory / name
    return [[float(x) for x in line.split()] for line in path.read_text().splitlines()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--interval", type=int, default=2)
    parser.add_argument("--restarts", type=int, nargs="+", default=[2, 3])
    parser.add_argument("--step-end", type=int, default=6)
    parser.add_argument("--digits", type=int, default=6)
    parser.add_argument("--expect-binned", action="store_true")
    parser.add_argument("--mpi-exec", default="")
    parser.add_argument("--mpi-numproc-flag", default="-n")
    parser.add_argument("--mpi-preflags", default="")
    parser.add_argument("--ranks", type=int, default=1)
    args = parser.parse_args()

    def label(step):
        return "%0*d" % (args.digits, step)

    def g_lnR_name(step):
        return "%s_g_lnR_%s.txt" % (NAME, label(step))

    straight = Path("straight").resolve()
    run(args, straight, [])
    restarts = []
    for step in args.restarts:
        directory = Path("restart_%d" % step).resolve()
        run(args, directory, ["erf.restart=%s/chk%s" % (straight, label(step))])
        restarts.append((directory, step))

    failed = []
    steps = range(0, args.step_end + 1, args.interval)
    num_binned = 0
    for directory, first in [(straight, 0)] + restarts:
        expected = [g_lnR_name(s) for s in steps if s >= first]
        found = sorted(p.name for p in directory.glob(NAME + "_g_lnR_*.txt"))
        if found != expected:
            failed.append("%s: g(ln R) files %s, expected %s" % (directory.name, found, expected))
        for binned in sorted(directory.glob(NAME + "_binned_dsd_*_[0-9]*")):
            if not binned.is_dir():
                continue
            num_binned += 1
            plotfile = directory / ("plt" + binned.name.rsplit("_", 1)[1])
            got, want = header_time_and_step(binned), header_time_and_step(plotfile)
            print("%s/%s: time %.17g step %d; %s: time %.17g step %d"
                  % (directory.name, binned.name, got[0], got[1], plotfile.name, want[0], want[1]))
            if abs(got[0] - want[0]) > 1.e-6 * max(1.0, abs(want[0])) or got[1] != want[1]:
                failed.append("%s/%s: time and step differ from %s"
                              % (directory.name, binned.name, plotfile.name))
    if args.expect_binned and num_binned == 0:
        failed.append("no binned distribution plotfiles, although the build writes them")

    for directory, first in restarts:
        for step in steps:
            if step < first:
                continue
            try:
                a, b = g_lnR(straight, g_lnR_name(step)), g_lnR(directory, g_lnR_name(step))
            except FileNotFoundError:
                continue
            same = len(a) == len(b) and all(
                len(ra) == len(rb) and all(abs(x - y) <= 1.e-5 * max(abs(x), abs(y), 1.e-300)
                                           for x, y in zip(ra, rb))
                for ra, rb in zip(a, b))
            if not same:
                failed.append("g(ln R) of step %d differs between the straight run and %s"
                              % (step, directory.name))

    for message in failed:
        print("FAIL: " + message)
    print("PASS" if not failed else "FAILED (%d checks)" % len(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
