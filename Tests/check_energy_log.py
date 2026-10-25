#!/usr/bin/env python3
"""Check the first row of erf.energy_data_log against the plotfile of the same time.

Usage: check_energy_log.py <plotfile> <energy log> [--tol T]

ERF writes the domain means of the total mass rho_d (1 + q_t) and of the total energy

    rho_d [ (c_vd + c_vv q_v + c_l q_l + c_i q_i) T - L_v q_l - L_s q_i
            + (1 + q_t) (|u|^2 / 2 + g z) ]

with q_v, q_l (cloud water + rain), q_i the mixing ratios and z the cell-centre
height. This script computes both from the plotfile's density, qv, qc, qrain (if
plotted), temp and cell-centred velocities on a flat, single-level, uniform mesh,
and compares them with the log, which carries 6 significant digits. It stops if the
plotfile has no cloud water, since the liquid terms would then go unchecked.
Exit status 0 on pass, 1 on a mismatch, 2 if the case cannot test the liquid terms.
"""

import argparse
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "Exec", "CanonicalTests", "Canonical_RANS"))
import erf_plotfile  # noqa: E402

# ERF_Constants.H and ERF_MicrophysicsConstants.H
R_D, R_V = 287.0, 461.505
CP_D, CP_V, CP_L, C_I = 1004.5, 1859.0, 4200.0, 2106.0
L_V, L_F = 2.5e6, 3.337e5
G = 9.81


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("plotfile")
    ap.add_argument("log")
    ap.add_argument("--tol", type=float, default=5e-6)
    a = ap.parse_args(argv)

    hdr = erf_plotfile._read_header(a.plotfile)
    names = ["density", "qv", "qc", "temp", "x_velocity", "y_velocity", "z_velocity"]
    has_rain = "qrain" in hdr["names"]
    if has_rain:
        names.append("qrain")
    _, d = erf_plotfile.read_fields(a.plotfile, names)
    lo, hi, dx = hdr["lo"], hdr["hi"], hdr["dx"]
    n = [hi[i] - lo[i] + 1 for i in range(3)]
    zlo = hdr["prob_lo"][2]

    mass = energy = 0.0
    qc_max = 0.0
    for i in range(n[0]):
        for j in range(n[1]):
            for k in range(n[2]):
                rho = d["density"][i][j][k]
                qv = d["qv"][i][j][k]
                ql = d["qc"][i][j][k] + (d["qrain"][i][j][k] if has_rain else 0.0)
                qc_max = max(qc_max, d["qc"][i][j][k])
                qt = qv + ql
                T = d["temp"][i][j][k]
                ke = 0.5 * sum(d[c][i][j][k] ** 2 for c in ("x_velocity", "y_velocity", "z_velocity"))
                z = zlo + (k + 0.5) * dx[2]
                mass += rho * (1.0 + qt)
                energy += rho * ((CP_D - R_D + (CP_V - R_V) * qv + CP_L * ql) * T - L_V * ql
                                 + (1.0 + qt) * (ke + G * z))
    ncell = n[0] * n[1] * n[2]
    mass /= ncell
    energy /= ncell
    if qc_max <= 0.0:
        print("FAIL: no cloud water in the plotfile; the liquid terms would go unchecked")
        return 2

    with open(a.log) as fh:
        rows = [ln.split() for ln in fh if ln.strip()]
    if len(rows) != 2:
        print("FAIL: %s should hold a header and one row (it is appended to; remove it first)" % a.log)
        return 2
    row = [float(v) for v in rows[1]]   # the data row, after the header
    status = 0
    for name, ref, val in (("tot_mass", mass, row[1]), ("tot_energy", energy, row[2])):
        rel = abs(val - ref) / abs(ref)
        ok = rel <= a.tol
        print("%-10s log %.6e  from the plotfile %.6e  relative %.2e  %s"
              % (name, val, ref, rel, "ok" if ok else "FAIL"))
        if not ok:
            status = 1
    print("PASS" if status == 0 else "FAIL")
    return status


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
