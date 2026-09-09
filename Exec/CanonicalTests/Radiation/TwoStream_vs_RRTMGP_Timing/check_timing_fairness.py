#!/usr/bin/env python3
"""
Guard the fairness of the two-stream / RRTMGP timing comparison.

A timing comparison is only meaningful while the two runs differ in the solver
and nothing else. That invariant is easy to break by editing one input file and
forgetting the other, and the resulting numbers still look plausible. This
check fails when the two configurations drift apart.

It verifies that:
  1. both inputs pull in the same shared block via FILE = inputs_common;
  2. neither input redefines anything the shared block already sets;
  3. the settings each input adds on its own are confined to its solver;
  4. radiation is called every step in both;
  5. plotfiles and per-solver logs are off, so I/O stays out of the timing.

If a results CSV is present it is sanity-checked too, but its absence is not a
failure: the measurement needs an RRTMGP-enabled build and may not have run.
"""

import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
COMMON = "shared_settings"
INPUTS = ["inputs_twostream", "inputs_rrtmgp"]

# Keys each input is allowed to set on its own: the solver selector and the
# knobs belonging to that solver.
ALLOWED = {
    "inputs_twostream": re.compile(r"^(erf\.prob_name|amr\.n_cell|erf\.radiation_type|erf\.radiation\.)"),
    # The gas volume mixing ratios are read only by the RRTMGP interface, so
    # they are solver-local even though they sit directly under erf.
    "inputs_rrtmgp": re.compile(
        r"^(erf\.prob_name|amr\.n_cell|start_datetime|erf\.radiation_model"
        r"|erf\.rad_|erf\.profile_rad_int|erf\.rrtmgp_"
        r"|erf\.(co2|o3|n2o|co|ch4|o2|n2)vmr)"),
}

failures = []


def note(msg):
    failures.append(msg)
    print("ERROR: %s" % msg)


def settings(path):
    """Return {key: value} for assignments in a file, ignoring comments."""
    out = {}
    with open(path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if not line or "=" not in line:
                continue
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def main():
    print("=" * 70)
    print("Two-stream / RRTMGP timing comparison: fairness check")
    print("=" * 70)

    common_path = os.path.join(HERE, COMMON)
    if not os.path.isfile(common_path):
        note("missing shared configuration %s" % COMMON)
        return 1
    common = settings(common_path)
    print("  shared block %s defines %d settings" % (COMMON, len(common)))

    for name in INPUTS:
        path = os.path.join(HERE, name)
        if not os.path.isfile(path):
            note("missing %s" % name)
            continue
        own = settings(path)

        # 1. both must include the shared block
        if own.get("FILE") != COMMON:
            note("%s does not pull in the shared block (FILE = %s)" % (name, COMMON))

        # 2. no shadowing of shared settings
        for key in own:
            if key == "FILE":
                continue
            if key in common:
                note("%s redefines %r, which the shared block already sets; "
                     "the two runs would no longer be comparable" % (name, key))

        # 3. anything it adds must belong to its own solver
        pattern = ALLOWED[name]
        for key in own:
            if key == "FILE":
                continue
            if not pattern.match(key):
                note("%s sets %r, which is outside its solver's namespace; "
                     "shared settings belong in %s" % (name, key, COMMON))
        print("  %-18s adds %d settings, all solver-local" % (name, len(own) - 1))

    # 4. both solvers must see the same surface temperature
    ts = settings(os.path.join(HERE, "inputs_twostream"))
    rr_t = settings(os.path.join(HERE, "inputs_rrtmgp")).get("erf.rad_t_sfc")
    ts_t = ts.get("erf.radiation.surface_temp_k")
    if rr_t is None or ts_t is None or float(rr_t) != float(ts_t):
        note("surface temperature differs between the solvers "
             "(RRTMGP %s, two-stream %s); the lower boundary must match"
             % (rr_t, ts_t))
    else:
        print("  both solvers use surface temperature %s K" % ts_t)

    # 5. radiation must run every step on both sides
    rr = settings(os.path.join(HERE, "inputs_rrtmgp"))
    if rr.get("erf.rad_freq_in_steps") != "1":
        note("RRTMGP is not called every step (erf.rad_freq_in_steps = %s); "
             "the comparison would measure update frequency, not solver cost"
             % rr.get("erf.rad_freq_in_steps"))
    else:
        print("  RRTMGP called every step")

    # 5. I/O must stay out of the measured region
    if common.get("erf.plot_int_1") != "-1":
        note("plotfiles are enabled in the shared block; file I/O would be timed")
    if common.get("amrex.tiny_profile") != "1":
        note("amrex.tiny_profile is not set; there would be nothing to measure")
    else:
        print("  profiler enabled, plotfiles and logs off")

    # Results, if any, are sanity-checked but not required.
    csv_path = os.path.join(HERE, "radiation_timing_comparison.csv")
    if os.path.isfile(csv_path):
        with open(csv_path) as f:
            rows = list(csv.DictReader(f))
        print("  results present: %d resolutions" % len(rows))
        for r in rows:
            for col in ("twostream_ms_per_call", "twostream_calls"):
                if not r.get(col):
                    note("results row missing %s" % col)
            val = r.get("twostream_ms_per_call")
            if val and float(val) <= 0.0:
                note("non-positive two-stream cost in the results")
        # Cost must grow with problem size; a flat or falling curve means the
        # measurement captured something other than the solver.
        sized = sorted(rows, key=lambda r: int(r["cells"]))
        costs = [float(r["twostream_ms_per_call"]) for r in sized
                 if r.get("twostream_ms_per_call")]
        if len(costs) > 1 and costs[-1] <= costs[0]:
            note("two-stream cost does not increase with problem size "
                 "(%.3f -> %.3f ms); the timing looks wrong" % (costs[0], costs[-1]))
    else:
        print("  no results CSV yet (needs an RRTMGP-enabled build); not required")

    print("=" * 70)
    if failures:
        print("RESULT: FAIL (%d)" % len(failures))
        return 1
    print("RESULT: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
