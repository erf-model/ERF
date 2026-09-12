#!/usr/bin/env python3
"""Check the partial-cell regtest from its final plotfile.

    python3 check_partial_cells.py plt19000

The run completing is the main check (the original selection traps). On
the plotfile: the fluid temperature stays within 0.5 K of the neutral
300 K everywhere (the checkerboard reached 10 K), the vertical velocity
stays below 2 m/s (it reached 5), and the cube still blocks the flow (a
wake with the streamwise velocity below a third of the inflow behind it).
"""
import sys
import numpy as np
import yt

def main():
    ds = yt.load(sys.argv[1]); ad = ds.all_data()
    th = ad["boxlib", "theta"].value; m = ad["boxlib", "terrain_IB_mask"].value
    u = ad["boxlib", "x_velocity"].value; w = ad["boxlib", "z_velocity"].value
    x = ad["index", "x"].in_units("code_length").value; y = ad["index", "y"].in_units("code_length").value
    z = ad["index", "z"].in_units("code_length").value
    fluid = m < 0.5
    ok = True
    def rep(name, good, detail):
        nonlocal ok; ok &= good; print(f"  {name}: {'PASS' if good else 'FAIL'} ({detail})")
    dth = np.abs(th[fluid] - 300.0).max()
    rep("fluid temperature within 0.5 K of neutral", dth < 0.5, f"max |theta - 300| {dth:.3f} K at t = {float(ds.current_time):.0f} s")
    wmax = np.abs(w[fluid]).max()
    rep("vertical velocity below 2 m/s", wmax < 2.0, f"max |w| {wmax:.2f} m/s")
    wake = fluid & (x > 180) & (x < 230) & (np.abs(y - 160) < 15) & (z < 40)
    rep("wake behind the cube", u[wake].min() < 1.0, f"min u in the wake {u[wake].min():.2f} m/s (inflow 3)")
    print("partial cells:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
