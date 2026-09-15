#!/usr/bin/env python3
"""Check the immersed-boundary awareness of the MRF and YSUNew schemes.

    python3 check_pbl_ib.py scheme plt_mrf_on00200 [plt_mrf_off00200 | none] run_mrf_off.log
    python3 check_pbl_ib.py flat   plt_flat_off00200 plt_flat_on00200

scheme: the run with the switch on must be finite everywhere with zero
  diffusivity in every solid cell, and the roof column's viscosity as a
  function of height above the roof must peak at the same local height as
  the open-ground column's as a function of height above the ground (the
  wall distance restarts at the roof); the first fluid cell above the roof
  must be within a factor five of the ground's first cell. The run with
  the switch off is reported: on this case it either fails (YSUNew drives
  the density negative at the second step) or fills the domain with NaN
  (MRF), because the ground surface layer evaluates u* and L on cells
  inside the building; that is the defect the switch removes.
flat: without a building the two plotfiles are bit-identical.
"""
import sys, os
import numpy as np
import yt

def load(pf, mask=True):
    ds = yt.load(pf); ad = ds.all_data()
    g = {"Kmv": ad["boxlib", "Kmv"].value}
    if mask: g["terrain_IB_mask"] = ad["boxlib", "terrain_IB_mask"].value
    for c in "xyz": g[c] = ad["index", c].in_units("code_length").value
    return g

def report(name, ok, detail):
    print(f"  {name}: {'PASS' if ok else 'FAIL'} ({detail})"); return ok

def column(g, x, y):
    m = (np.abs(g["x"] - x) < 1) & (np.abs(g["y"] - y) < 1)
    o = np.argsort(g["z"][m]); return g["z"][m][o], g["Kmv"][m][o]

def main():
    mode = sys.argv[1]; ok = True
    if mode == "flat":
        a, b = load(sys.argv[2], mask=False), load(sys.argv[3], mask=False)
        d = np.abs(a["Kmv"] - b["Kmv"]).max()
        ok &= report("no building: switch on and off bit-identical", d == 0.0, f"max |dK| {d:.1e}")
    else:
        on = load(sys.argv[2]); off_pf = sys.argv[3]; off_log = sys.argv[4]
        solid = on["terrain_IB_mask"] >= 0.5
        ok &= report("switch on: finite everywhere, zero in every solid cell",
                     np.isfinite(on["Kmv"]).all() and on["Kmv"][solid].max() == 0.0,
                     f"{int((~np.isfinite(on['Kmv'])).sum())} non-finite, max K inside the cube {on['Kmv'][solid].max():.1e}")
        zr, kr = column(on, 165, 165); zg, kg = column(on, 45, 45)
        fr = kr[zr > 40.0]; zfr = zr[zr > 40.0] - 40.0          # above the roof, local height
        # Shape of the profile against the local height: the roof column's,
        # normalised by its own maximum, must follow the ground column's over
        # the first six cells (the wall distance restarts at the roof). For a
        # flat profile (a neutral YSUNew) the peak position is meaningless and
        # only the shape is compared.
        n = 6; sr = fr[:n] / fr[:n].max(); sg = kg[:n] / kg[:n].max()
        flat = (kg[:n].max() / max(kg[:n].min(), 1e-12)) < 1.2
        ipk_r = np.argmax(fr[:n]); ipk_g = np.argmax(kg[:n])
        ok &= report("switch on: the roof column's profile against height above the roof follows the ground's against height above the ground",
                     np.abs(sr - sg).max() < 0.35 and (flat or abs(zfr[ipk_r] - zg[ipk_g]) <= 20.0),
                     f"shape difference {np.abs(sr - sg).max():.2f}; " + ("flat profile" if flat else f"peak {zfr[ipk_r]:.0f} m above the roof, {zg[ipk_g]:.0f} m above the ground"))
        r = fr[0] / max(kg[0], 1e-12)
        ok &= report("switch on: first cell above the roof within a factor five of the first cell above the ground", 0.2 < r < 5.0,
                     f"roof {fr[0]:.3f}, ground {kg[0]:.3f}")
        # the run without the switch, for the record
        if os.path.isdir(off_pf):
            off = load(off_pf); nn = int((~np.isfinite(off["Kmv"])).sum())
            print(f"  switch off: run completed with {nn} non-finite viscosity cells" + (" (unusable)" if nn else ""))
        else:
            msg = [l.strip() for l in open(off_log) if "Rho is negative" in l or "Abort" in l]
            print("  switch off: run failed: " + (msg[0][:80] if msg else "see the log"))
    print(f"{mode}: {'PASS' if ok else 'FAIL'}"); sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
