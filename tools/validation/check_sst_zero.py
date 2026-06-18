#!/usr/bin/env python3
"""SST=0 pre-launch gate. Verify NO water cell has SST<=0 (or t_surf=0) in a
wrfinput/wrflowinp/met_em file (all time levels). Exit 0 = clean, 1 = holes.

Usage: check_sst_zero.py <file> [file2 ...] [--mask <wrfinput_with_XLAND>]
Water = XLAND>1.5 if present, else LANDMASK<0.5. wrflowinp has no land mask, so
pass --mask pointing at the matching wrfinput (auto-used for any landmask-less file).
"""
import sys, numpy as np, netCDF4 as nc

_EXTERN_MASK = [None]  # set from --mask

def water_mask(d):
    if "XLAND" in d.variables:
        return np.asarray(d.variables["XLAND"][:]) > 1.5
    if "LANDMASK" in d.variables:
        return np.asarray(d.variables["LANDMASK"][:]) < 0.5
    if _EXTERN_MASK[0] is not None:        # landmask-less file (e.g. wrflowinp)
        return _EXTERN_MASK[0]
    return None

def check(path):
    d = nc.Dataset(path)
    if "SST" not in d.variables:
        print(f"  {path}: no SST var -> skip"); d.close(); return 0
    sst = np.asarray(d.variables["SST"][:])           # (time,y,x) or (y,x)
    if sst.ndim == 2: sst = sst[None]
    w = water_mask(d)
    if w is None:
        print(f"  {path}: no land mask (pass --mask wrfinput) -> cannot classify"); d.close(); return 2
    if w.ndim == 2: w = w[None]
    W = np.broadcast_to(w, sst.shape) if w.shape[0]==1 else w
    if W.shape != sst.shape:  # static mask vs time-varying SST
        W = np.broadcast_to(w[:1], sst.shape)
    holes = int((W & (sst <= 0)).sum())
    near  = int((W & (sst < 150) & (sst > 0)).sum())
    nt = sst.shape[0]
    status = "CLEAN" if holes==0 and near==0 else "HOLES"
    print(f"  {path}: times={nt} water SST<=0={holes} water 0<SST<150={near}  [{status}]")
    d.close()
    return 0 if (holes==0 and near==0) else 1

if __name__ == "__main__":
    args = sys.argv[1:]
    if "--mask" in args:
        mi = args.index("--mask"); mf = args[mi+1]
        md = nc.Dataset(mf)
        _EXTERN_MASK[0] = (np.asarray(md.variables["XLAND"][0]) > 1.5
                           if "XLAND" in md.variables
                           else np.asarray(md.variables["LANDMASK"][0]) < 0.5)
        md.close()
        args = args[:mi] + args[mi+2:]
    if not args:
        print(__doc__); sys.exit(2)
    print("=== SST=0 pre-launch gate ===")
    rc = 0
    for f in args:
        try: rc |= check(f)
        except Exception as e: print(f"  {f}: ERROR {e}"); rc |= 2
    print("RESULT:", "ALL CLEAN" if rc==0 else "HOLES/ERROR FOUND (do NOT run)")
    sys.exit(rc)
