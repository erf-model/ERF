#!/usr/bin/env python3
"""Sensitivity check for check_theta_above_ground.py (inputs_theta with the flag off).

Usage: check_theta_flag_off.py [--smoke] <plotfile>

With erf.input_sounding_theta_above_ground = false the interior starts with the
inversion at physical heights, so the step-0 theta must differ from the
sounding's theta above the local ground: the median error must be at least
0.5 K and the maximum at least 4 K. If this fails, the comparison in
check_theta_above_ground.py could pass without the flag doing anything.
"""

import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_here, os.path.join(_here, "..")]
import rans_checks as rc  # noqa: E402
from check_theta_above_ground import measure, median, step0_plotfile  # noqa: E402


def main(argv):
    args = [a for a in argv[1:] if not a.startswith("--")]
    if len(args) != 1:
        print(__doc__)
        return 2
    m = measure(args[0])
    rep = rc.Report()
    print("step 0 from %s: %d cells" % (step0_plotfile(args[0]), len(m["err"])))
    rep.check("flag off, step 0: median |theta - theta(z above ground)| [K]", median(m["err"]), 0.5, 0.0, "min")
    rep.check("flag off, step 0: max |theta - theta(z above ground)| [K]", max(m["err"]), 4.0, 0.0, "min")
    rep.check("final plotfile: all fields finite", 1.0 if m["finite"] else 0.0, 1.0, 0.0)
    rep.dump()
    return 1 if rep.failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
