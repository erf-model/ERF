#!/usr/bin/env python3
"""Check ERF inputs files before running them.

    python3 Exec/Tools/check_erf_inputs.py inputs [more inputs ...]
    python3 Exec/Tools/check_erf_inputs.py --ranks 64 --umax 30 Exec/ABL/inputs_*
    python3 Exec/Tools/check_erf_inputs.py --min-severity WARNING --strict Exec/RegTests/*/inputs*
    python3 Exec/Tools/check_erf_inputs.py --self-test

For each deck the script prints, from most to least severe:

    ERROR    the run would stop at start-up
    WARNING  the run starts but the setting is very likely not what was meant
    ADVICE   a setting that could be chosen better
    INFO     derived numbers (cell sizes, time-step estimates, box counts)

The exit code is 1 when any deck has an ERROR, or with --strict when any
deck has a WARNING.  See erf_input_checker/__init__.py for the design and
erf_input_checker/checks/ for what is checked.
"""

import argparse
import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from erf_input_checker import (Report, RunContext, Severity, parse_deck,  # noqa: E402
                               run_all_checks)
from erf_input_checker.context import Options  # noqa: E402
from erf_input_checker.deck import DeckError  # noqa: E402
from erf_input_checker.checker import registered_checkers  # noqa: E402


def check_one(path: str, opts: Options, only=None) -> Report:
    report = Report(path)
    try:
        deck = parse_deck(path)
    except (DeckError, OSError) as exc:
        from erf_input_checker.findings import Finding
        report.add(Finding(Severity.ERROR, "deck", "parse", str(exc), "Fix the line so ParmParse can read it."))
        return report
    ctx = RunContext(deck, opts)
    run_all_checks(deck, ctx, report, only=only)
    return report


_INCLUDE_RE = re.compile(r"^\s*FILE\s*=\s*\"?([^\s\"#]+)")


def included_by(path: str):
    """Name of a sibling file that includes ``path`` with FILE =, or None."""
    d = os.path.dirname(os.path.abspath(path))
    base = os.path.basename(path)
    try:
        names = os.listdir(d)
    except OSError:
        return None
    for fn in names:
        if fn == base or not fn.startswith("inputs") and not fn.endswith(".i"):
            continue
        full = os.path.join(d, fn)
        if not os.path.isfile(full):
            continue
        try:
            with open(full, errors="replace") as fh:
                for line in fh:
                    m = _INCLUDE_RE.match(line)
                    if m and os.path.basename(m.group(1)) == base:
                        return fn
        except OSError:
            continue
    return None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="*", help="inputs files (globs are expanded)")
    ap.add_argument("--ranks", type=int, help="MPI ranks the run will use; enables the box-count balance advice")
    ap.add_argument("--gpu", action="store_true", help="the run is on GPUs (box-size advice changes)")
    ap.add_argument("--umax", type=float, default=20.0, help="largest wind speed [m/s] for the CFL estimates (default 20)")
    ap.add_argument("--theta", type=float, default=300.0, help="potential temperature [K] for the sound speed (default 300)")
    ap.add_argument("--min-severity", default="INFO", choices=[s.name for s in Severity], help="hide findings below this level")
    ap.add_argument("--strict", action="store_true", help="exit 1 on WARNING as well as ERROR")
    ap.add_argument("--json", action="store_true", help="JSON output")
    ap.add_argument("--only", action="append", help="run only the named checker(s), e.g. --only grid --only timestep.compressible")
    ap.add_argument("--list-checkers", action="store_true")
    ap.add_argument("--docs", help="path to Docs/sphinx_doc/Inputs.rst (found automatically inside the repository)")
    ap.add_argument("--include-bases", action="store_true",
                    help="also check files that a sibling deck pulls in with FILE = (skipped by default, since they are not complete decks)")
    ap.add_argument("--self-test", action="store_true", help="run the built-in tests and exit")
    args = ap.parse_args(argv)

    if args.self_test:
        from erf_input_checker import selftest
        return selftest.main()
    if args.list_checkers:
        import erf_input_checker.checks  # noqa: F401
        for cls in registered_checkers():
            print("%-24s %s" % (cls.name, (cls.__doc__ or "").strip().splitlines()[0] if cls.__doc__ else ""))
        return 0
    if args.docs:
        from erf_input_checker.checks.keys import UnknownKeys
        UnknownKeys.docs_path = args.docs

    paths = []
    for pat in args.inputs:
        hits = sorted(glob.glob(pat))
        paths.extend(hits if hits else [pat])
    if not paths:
        ap.error("give at least one inputs file")

    opts = Options(ranks=args.ranks, gpu=args.gpu, umax=args.umax, theta=args.theta)
    worst = Severity.INFO
    minsev = Severity[args.min_severity]
    reports = []
    for path in paths:
        if os.path.isdir(path):
            continue
        parent = included_by(path)
        if parent and not args.include_bases:
            if not args.json:
                print("== %s\n-- skipped: included by %s (use --include-bases to check it alone)\n" % (path, parent))
            continue
        if "master" in os.path.basename(path) and not args.include_bases:
            if not args.json:
                print("== %s\n-- skipped: a *master* deck lists every key for documentation and is not meant to run (use --include-bases)\n" % path)
            continue
        rep = check_one(path, opts, only=args.only)
        reports.append(rep)
        if rep.worst is not None and rep.worst > worst:
            worst = rep.worst
        if args.json:
            print(rep.as_json())
        else:
            print(rep.format(minsev))
            print()
    if len(reports) > 1 and not args.json:
        print("== %d decks: %d with errors, %d with warnings" % (
            len(reports), sum(1 for r in reports if r.count(Severity.ERROR)),
            sum(1 for r in reports if r.count(Severity.WARNING))))
    if worst >= Severity.ERROR or (args.strict and worst >= Severity.WARNING):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
