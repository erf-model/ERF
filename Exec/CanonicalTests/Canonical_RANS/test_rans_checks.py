#!/usr/bin/env python3
"""Self-test of Report.check in rans_checks.py: the pass/fail logic itself.

    python3 test_rans_checks.py

Every case states what the check must decide; the script reports each
disagreement and exits non-zero. The range cases are the ones that motivated
this test: the comparison used to accept half a band width outside the band.

It also runs the checkers' command lines with malformed arguments: each must
print its usage and exit with status 2 before reading any plotfile, rather than
fail with a traceback (check_implicit_explicit_ke.py raised StopIteration when
--tol was the last argument).
"""
import os
import subprocess
import sys
import rans_checks as rc

HERE = os.path.dirname(os.path.abspath(__file__))

CLI_CASES = [
    # (checker relative to this directory, arguments)
    ("Convective_ABL_Flat/check_implicit_explicit_ke.py", []),
    ("Convective_ABL_Flat/check_implicit_explicit_ke.py", ["--tol"]),
    ("Convective_ABL_Flat/check_implicit_explicit_ke.py", ["a_plt", "b_plt", "--tol"]),
    ("Convective_ABL_Flat/check_implicit_explicit_ke.py", ["--tol", "abc", "a_plt", "b_plt"]),
    ("Convective_ABL_Flat/check_implicit_explicit_ke.py", ["--tol", "1e-4", "a_plt"]),
]

CASES = [
    # (name, value, target, tol, kind, expected pass)
    ("abs inside",       1.05, 1.0,        0.1, "abs",   True),
    ("abs outside",      1.20, 1.0,        0.1, "abs",   False),
    ("rel inside",       1.05, 1.0,        0.1, "rel",   True),
    ("rel outside",      1.20, 1.0,        0.1, "rel",   False),
    ("min at target",    8.0,  8.0,        0.0, "min",   True),
    ("min below",        7.9,  8.0,        0.0, "min",   False),
    ("max at target",    2.0,  2.0,        0.0, "max",   True),
    ("max above",        2.1,  2.0,        0.0, "max",   False),
    ("range low edge",   0.5,  (0.5, 2.0), 0.0, "range", True),
    ("range middle",     1.25, (0.5, 2.0), 0.0, "range", True),
    ("range high edge",  2.0,  (0.5, 2.0), 0.0, "range", True),
    ("range just above", 2.03, (0.5, 2.0), 0.0, "range", False),
    ("range just below", 0.47, (0.5, 2.0), 0.0, "range", False),
    ("range far above",  2.70, (0.5, 2.0), 0.0, "range", False),
]


def main():
    bad = 0
    for name, value, target, tol, kind, expected in CASES:
        rep = rc.Report()
        rep.check(name, value, target, tol, kind)
        got = (rep.failed == 0)
        if got != expected:
            print("  %s: value %g with %s target %s -> %s, expected %s"
                  % (name, value, kind, target, "pass" if got else "fail",
                     "pass" if expected else "fail"))
            bad += 1
    print("rans_checks self-test: %d of %d cases as expected"
          % (len(CASES) - bad, len(CASES)))

    cli_bad = 0
    for script, args in CLI_CASES:
        p = subprocess.run([sys.executable, "-B", os.path.join(HERE, script)] + args,
                           capture_output=True, text=True)
        if p.returncode != 2 or "Traceback" in p.stderr:
            print("  %s %s: exit %d, expected 2 with the usage%s"
                  % (script, " ".join(args), p.returncode,
                     " (traceback)" if "Traceback" in p.stderr else ""))
            cli_bad += 1
    print("checker command lines: %d of %d malformed calls rejected with the usage"
          % (len(CLI_CASES) - cli_bad, len(CLI_CASES)))
    sys.exit(1 if (bad or cli_bad) else 0)


if __name__ == "__main__":
    main()
