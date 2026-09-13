#!/usr/bin/env python3
"""Check the phase 4 wall function and the heat budget.

    python3 check_sensible.py formula faces_diag 0.01 0.001 1004.5
    python3 check_sensible.py budget  run_couple.log plt_diag_00040 plt00040 0.125 1004.5

formula: on every face of the dump, H must equal
    rho c_p u* kappa (theta_skin - theta_air) / ln(delta / z0h),
    u* = kappa U_tan / ln(delta / z0),  theta_skin = T_skin theta_air / T_air,
    delta = half the cell size in the face direction (2.5 m here),
and be positive (the faces are warmer than the air).
budget: the heat the faces added over the run, sum over steps of
    H_total dt, against the difference of the air's internal energy
    sum(rho c_v T dV) between the coupled run and the diagnostic run at the
    same step. The two runs differ only in whether the face flux is applied,
    so the difference isolates the heating from the (much larger) background
    adjustment of the whole domain. The domain is a rigid closed box (periodic
    sides, slip top, adiabatic ground), so heat added raises the internal
    energy, c_v = c_p - R_d; the enthalpy rises by c_p/c_v times more.
"""
import sys, glob, re, os
import numpy as np

KAPPA = 0.41; DX = 5.0

def load(prefix):
    rows = []
    for fn in sorted(glob.glob(prefix + ".rank*.csv")):
        with open(fn) as f:
            hdr = f.readline().strip().split(",")
            for line in f: rows.append([float(v) for v in line.strip().split(",")])
    a = np.array(rows); return {h: a[:, n] for n, h in enumerate(hdr)}

def read_plotfile(path):
    hdr = open(os.path.join(path, "Header")).read().split("\n")
    nvar = int(hdr[1]); names = hdr[2:2 + nvar]; i = 2 + nvar + 1 + 2 + 3
    m = re.findall(r"\((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\)", hdr[i])[0]
    lo = [int(m[0]), int(m[1]), int(m[2])]; n = [int(m[3]) - lo[0] + 1, int(m[4]) - lo[1] + 1, int(m[5]) - lo[2] + 1]
    cell_h = open(os.path.join(path, "Level_0", "Cell_H")).read().split("\n")
    k = 0
    while not cell_h[k].startswith("("): k += 1
    nfab = int(cell_h[k].strip("(").split()[0]); k += 1
    boxes = []
    for _ in range(nfab):
        b = re.findall(r"\((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\)", cell_h[k])[0]; boxes.append([int(v) for v in b]); k += 1
    fabs = [(l.split()[1], int(l.split()[2])) for l in cell_h if l.startswith("FabOnDisk:")]
    data = {nm: np.zeros(n) for nm in names}
    for b, (fname, off) in zip(boxes, fabs):
        with open(os.path.join(path, "Level_0", fname), "rb") as f:
            f.seek(off); f.readline()
            nx, ny, nz = b[3] - b[0] + 1, b[4] - b[1] + 1, b[5] - b[2] + 1
            arr = np.frombuffer(f.read(8 * nx * ny * nz * nvar), dtype="<f8").reshape(nvar, nz, ny, nx)
        sl = (slice(b[0]-lo[0], b[3]-lo[0]+1), slice(b[1]-lo[1], b[4]-lo[1]+1), slice(b[2]-lo[2], b[5]-lo[2]+1))
        for c, nm in enumerate(names): data[nm][sl] = arr[c].transpose(2, 1, 0)
    return data

def check_formula(prefix, z0, z0h, cp):
    d = load(prefix); ok = True
    delta = 0.5 * DX
    ustar = KAPPA * d["U_tan"] / np.log(delta / z0)
    th_skin = d["T_skin"] * d["theta_air"] / d["T_air"]
    H = d["rho"] * cp * ustar * KAPPA * (th_skin - d["theta_air"]) / np.log(delta / z0h)
    def rep(name, cond, extra=""):
        nonlocal ok; print(f"  {name:52s} {'PASS' if cond else 'FAIL'} {extra}"); ok = ok and bool(cond)
    rep("u* = kappa U_tan / ln(delta/z0) on every face", np.allclose(d["ustar"], ustar, rtol=1e-7, atol=1e-12))
    rep("H = rho c_p u* kappa (theta_skin - theta_air)/ln(delta/z0h)", np.allclose(d["H"], H, rtol=1e-7, atol=1e-9))
    rep("H > 0 on every face (walls warmer than the air)", np.all(d["H"] > 0), f"H {d['H'].min():.1f}-{d['H'].max():.1f} W/m2")
    rep("H largest on the windward face (largest U_tan)", d["U_tan"][np.argmax(d["H"])] > 0.5 * d["U_tan"].max())
    return ok

def check_budget(log, plt_diag, plt_couple, dt, cp):
    H = [float(re.search(r"H_total_W=([-\d.eE+]+)", l).group(1)) for l in open(log) if "H_total_W=" in l]
    # One [IBSEB] line per step from post_timestep plus the one at initialisation;
    # the flux applied during step n is the one computed at its start.
    Q_faces = sum(H[:-1]) * dt
    a = read_plotfile(plt_diag); b = read_plotfile(plt_couple)
    V = DX ** 3; cv = cp - 287.0
    dE = cv * np.sum(b["density"] * b["temp"] - a["density"] * a["temp"]) * V
    rel = abs(dE - Q_faces) / Q_faces
    ok = rel < 0.10
    print(f"  heat from the faces {Q_faces/1e6:.3f} MJ, extra internal energy of the coupled air {dE/1e6:.3f} MJ, relative difference {rel*100:.1f}% -> {'PASS' if ok else 'FAIL'}")
    return ok

if __name__ == "__main__":
    if sys.argv[1] == "formula": ok = check_formula(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))
    else: ok = check_budget(sys.argv[2], sys.argv[3], sys.argv[4], float(sys.argv[5]), float(sys.argv[6]))
    sys.exit(0 if ok else 1)
