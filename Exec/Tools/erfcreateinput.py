#!/usr/bin/env python3
"""Write a commented ERF inputs file, and the dummy files it names, for a chosen configuration.

    python3 Exec/Tools/erfcreateinput.py --list
    python3 Exec/Tools/erfcreateinput.py --preset abl_convective_les --outdir my_cbl
    python3 Exec/Tools/erfcreateinput.py --setup abl --turbulence deardorff --surface most \\
        --roughness file --coriolis --forcing geostrophic --output plt,chk,profiles --outdir my_abl
    python3 Exec/Tools/erfcreateinput.py --setup bubble --dycore anelastic --outdir my_bubble
    python3 Exec/Tools/erfcreateinput.py --self-test

The chosen option of each group is active in the file; the alternatives are
commented lines beside it.  See erf_input_creator/README.md.

This script complements check_erf_inputs.py: it writes a new case, and the
checker then checks the case and the files it names before a run.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from erf_input_creator.generate import CaseError, build, command_line, summary, write  # noqa: E402
from erf_input_creator.options import OptionError, case_options, list_choices, parse, resolve  # noqa: E402


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    try:
        args = parse(argv)
    except OptionError as exc:
        print("erfcreateinput: %s" % exc, file=sys.stderr)
        return 2
    if args.self_test:
        from erf_input_creator import selftest
        return selftest.main(schema_path=args.schema)
    if args.list:
        print("\n".join(list_choices()))
        return 0
    try:
        spec = resolve(args)
        deck, files = build(spec, command_line(argv))
        write(spec, deck, files, case_options(args))
    except CaseError as exc:
        for line in str(exc).splitlines():
            print("erfcreateinput: " + line, file=sys.stderr)
        return 2
    print("\n".join(summary(spec, deck, files)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
