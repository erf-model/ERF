#!/usr/bin/env python3
"""
Generate wrfinput_chisholmview_d02 from wrfinput_chisholmview_d01 for the
inputs_wps_ml_interp regression test (ERF issue #4060, surface-only init path).

The d02 is a 3x horizontal refinement of an interior subdomain of d01 with the
same eta levels, carrying only the variables ERF's surface-only read needs
(init_from_wrfinput with read_atmos_state=false), plus the base-state scalars.
Its start time is d01's plus DELAY_S seconds, so ERF creates level 1 mid-run
via MakeNewLevelFromCoarse -- the trigger for the surface-only path.

Refinement: piecewise-constant (np.repeat) for cell-centered fields so every
fine value equals its parent value; linear interpolation for horizontally
staggered fields (MAPFAC_U/V) so the nodes stay monotone.
"""
import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta

SRC = "wrfinput_d01"
DST = "wrfinput_d02"

RATIO = 3
IS, JS = 35, 35   # 0-based coarse cell of the fine patch lower-left corner
NCX, NCY = 30, 30   # coarse cells covered by the fine patch
DELAY_S = 25        # d02 start = d01 start + 25 s -> level 1 appears at step 5 (dt = 5)

NFX, NFY = RATIO * NCX, RATIO * NCY   # 90 x 90 fine cells

# Variables ERF reads on the surface-only path (ERF_InitFromWRFInput.cpp,
# read_atmos_state = false, LSM active), plus the base-state scalars
# (read_base_state_params_from_wrfinput) and the Times stamp.
KEEP = [
    "Times",
    "PH", "PHB",                       # terrain (z from geopotential)
    "ALB", "AL",                       # density (CheckForDensity requires it even surface-only)
    "PSFC", "MUB",
    "MAPFAC_U", "MAPFAC_V", "MAPFAC_M",
    "SST", "TSK", "LANDMASK",
    "C1H", "C2H", "RDNW",
    "XLAT", "XLONG",
    "IVGTYP", "ISLTYP",
    "TSLB", "SMOIS", "SH2O", "LAI",    # LSM fields (SLM::wrfinput_map)
    "ZS", "DZS", "VEGFRA", "TMN", "SHDMIN", "SHDMAX",
    "HGT",                             # not read today; kept for tooling
    "T00", "P00", "TLP", "TISO", "TLP_STRAT", "P_STRAT",
]

def refine_unstag(a, axis, start, ncoarse):
    """Slice [start, start+ncoarse) along axis and repeat each entry RATIO times."""
    sl = [slice(None)] * a.ndim
    sl[axis] = slice(start, start + ncoarse)
    return np.repeat(a[tuple(sl)], RATIO, axis=axis)

def refine_stag(a, axis, start, ncoarse):
    """Linear interpolation of coarse nodes [start, start+ncoarse] onto the
    RATIO*ncoarse+1 fine nodes that span the same interval."""
    sl = [slice(None)] * a.ndim
    sl[axis] = slice(start, start + ncoarse + 1)
    sub = a[tuple(sl)]
    xf = np.arange(RATIO * ncoarse + 1) / RATIO   # fine nodes in coarse node coords
    i0 = np.minimum(xf.astype(int), ncoarse - 1)
    w = xf - i0
    lo = np.take(sub, i0, axis=axis)
    hi = np.take(sub, i0 + 1, axis=axis)
    shape = [1] * a.ndim
    shape[axis] = len(xf)
    w = w.reshape(shape)
    out = lo * (1.0 - w) + hi * w
    return out.astype(a.dtype)

src = nc.Dataset(SRC)
dst = nc.Dataset(DST, "w", format="NETCDF4")

# --- dimensions ---
dim_size = {
    "west_east": NFX,        "south_north": NFY,
    "west_east_stag": NFX + 1, "south_north_stag": NFY + 1,
}
needed_dims = set()
for v in KEEP:
    needed_dims.update(src.variables[v].dimensions)
for d in needed_dims:
    n = dim_size.get(d, len(src.dimensions[d]))
    dst.createDimension(d, n)

# --- global attributes: copy, then override the nest bookkeeping ---
t0 = datetime.strptime(src.getncattr("START_DATE"), "%Y-%m-%d_%H:%M:%S")
t2 = (t0 + timedelta(seconds=DELAY_S)).strftime("%Y-%m-%d_%H:%M:%S")

attrs = {a: src.getncattr(a) for a in src.ncattrs()}
attrs.update({
    "START_DATE": t2, "SIMULATION_START_DATE": t2,
    "WEST-EAST_GRID_DIMENSION": np.int32(NFX + 1),
    "SOUTH-NORTH_GRID_DIMENSION": np.int32(NFY + 1),
    "DX": np.float32(src.getncattr("DX") / RATIO),
    "DY": np.float32(src.getncattr("DY") / RATIO),
    "I_PARENT_START": np.int32(IS + 1),   # 1-based, ERF subtracts 1
    "J_PARENT_START": np.int32(JS + 1),
    "PARENT_GRID_RATIO": np.int32(RATIO),
    "PARENT_ID": np.int32(1), "GRID_ID": np.int32(2),
    "WEST-EAST_PATCH_START_UNSTAG": np.int32(1),
    "WEST-EAST_PATCH_END_UNSTAG": np.int32(NFX),
    "WEST-EAST_PATCH_START_STAG": np.int32(1),
    "WEST-EAST_PATCH_END_STAG": np.int32(NFX + 1),
    "SOUTH-NORTH_PATCH_START_UNSTAG": np.int32(1),
    "SOUTH-NORTH_PATCH_END_UNSTAG": np.int32(NFY),
    "SOUTH-NORTH_PATCH_START_STAG": np.int32(1),
    "SOUTH-NORTH_PATCH_END_STAG": np.int32(NFY + 1),
})
dst.setncatts(attrs)

# --- variables ---
for name in KEEP:
    sv = src.variables[name]
    dv = dst.createVariable(name, sv.dtype, sv.dimensions)
    dv.setncatts({a: sv.getncattr(a) for a in sv.ncattrs()})
    if name == "Times":
        stamp = np.array(list(t2.ljust(len(src.dimensions["DateStrLen"]))), dtype="S1")
        dv[0, :] = stamp
        continue
    data = sv[:].data if np.ma.isMaskedArray(sv[:]) else sv[:]
    dims = sv.dimensions
    for axis, d in enumerate(dims):
        if d == "west_east":
            data = refine_unstag(data, axis, IS, NCX)
        elif d == "south_north":
            data = refine_unstag(data, axis, JS, NCY)
        elif d == "west_east_stag":
            data = refine_stag(data, axis, IS, NCX)
        elif d == "south_north_stag":
            data = refine_stag(data, axis, JS, NCY)
    dv[:] = data

dst.close()
src.close()

chk = nc.Dataset(DST)
print("wrote", DST)
print("dims:", {k: len(v) for k, v in chk.dimensions.items()})
for a in ["START_DATE", "DX", "I_PARENT_START", "J_PARENT_START",
          "PARENT_GRID_RATIO", "WEST-EAST_GRID_DIMENSION",
          "BOTTOM-TOP_GRID_DIMENSION", "USE_THETA_M"]:
    print(f"  {a} = {chk.getncattr(a)}")
print("Times:", chk.variables["Times"][:].tobytes().decode())
for v in ["TSK", "TSLB", "SMOIS", "PH", "PHB", "VEGFRA"]:
    arr = chk.variables[v][:]
    print(f"  {v}: shape {arr.shape}  min {arr.min():.4g}  max {arr.max():.4g}")
chk.close()
