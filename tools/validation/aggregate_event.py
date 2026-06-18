#!/usr/bin/env python3
"""Event-aggregate ERF output and observations onto a common footprint.

The cadences don't match (ERF 30 min, MRMS 2 min, NEXRAD ~5 min, ASOS ~5 min),
so instead of per-timestep collocation we compare EVENT AGGREGATES, which are
cadence-robust:

  - peak 10 m wind   : max over the event   (ERF wspd10max vs ASOS max gust/2-min)
  - total precip     : sum over the event   (ERF rain_accum final vs MRMS QPE sum)
  - max reflectivity : max over the event   (ERF max_reflectivity-over-plotfiles
                       or the running cmpref_max field, vs MRMS/NEXRAD composite max)

This script builds each aggregate and writes them to NetCDF on the ERF grid plus
a per-station CSV for ASOS. Regridding obs -> ERF grid is nearest-neighbour
(MRMS 1 km -> ERF 3/9 km is a mild coarsening; fine for event maxima).

Subcommands:
  erf      <plt2d_glob> <plt_glob>      -> erf_event.nc  (wspd10max, rain_accum, refl_max)
  mrms     <mrms_dir>                   -> mrms_event.nc  (refl_max, qpe_total)
  asos     <asos_dir>                   -> asos_event.csv (peak wind, peak gust, ...)
  nexrad   <nexrad_dir>                 -> nexrad_event.nc (per-site gridded refl_max)  [needs pyart]

This is a SCAFFOLD: the heavy regridding/exact-grid bits are marked TODO and
filled once the ERF run + obs are on disk. The aggregation math (max/sum over
the time axis) is the part that matters and is implemented.
"""
import sys, os, glob, csv

def erf_aggregate(plt2d_glob, plt_glob, out="erf_event.nc"):
    """Event aggregates from ERF plotfiles.

    wspd10max : already cumulative in the run -> just read the LAST plt2d.
    rain_accum: cumulative in the run        -> read the LAST plt (3D, k=0).
    refl_max  : if cmpref_max present in plt2d, read last; else take pixelwise
                max of max_reflectivity across all plt (3D -> k=0 slab).
    """
    import numpy as np, yt
    yt.set_log_level(50)
    p2 = sorted(glob.glob(plt2d_glob), key=_stepkey)
    p3 = sorted(glob.glob(plt_glob),  key=_stepkey)
    if not p2:
        print("no plt2d match", file=sys.stderr); return
    last2 = yt.load(p2[-1])
    names2 = [f[1] for f in last2.field_list if f[0]=="boxlib"]
    cg = last2.covering_grid(0, last2.domain_left_edge, last2.domain_dimensions)
    out_fields = {}
    if "wspd10max" in names2:
        out_fields["wspd10max"] = np.asarray(cg["boxlib","wspd10max"])[:,:,0]
    if "cmpref_max" in names2:
        out_fields["refl_max"] = np.asarray(cg["boxlib","cmpref_max"])[:,:,0]
    # rain_accum (final cumulative) from last 3D plt
    if p3:
        last3 = yt.load(p3[-1]); n3=[f[1] for f in last3.field_list if f[0]=="boxlib"]
        cg3 = last3.covering_grid(0, last3.domain_left_edge, last3.domain_dimensions)
        if "rain_accum" in n3:
            out_fields["rain_accum"] = np.asarray(cg3["boxlib","rain_accum"])[:,:,0]
        # fallback refl_max = max over plotfiles of max_reflectivity
        if "refl_max" not in out_fields and "max_reflectivity" in n3:
            acc=None
            for p in p3:
                ds=yt.load(p); g=ds.covering_grid(0,ds.domain_left_edge,ds.domain_dimensions)
                m=np.asarray(g["boxlib","max_reflectivity"])[:,:,0]
                acc=m if acc is None else np.maximum(acc,m)
            out_fields["refl_max"]=acc
    _write_nc(out, out_fields, "ERF event aggregates")
    for k,v in out_fields.items():
        print(f"ERF {k}: min/mean/max = {np.nanmin(v):.2f}/{np.nanmean(v):.2f}/{np.nanmax(v):.2f}")

def mrms_aggregate(mrms_dir, out="mrms_event.nc"):
    """Event aggregates from MRMS grib2: composite-reflectivity MAX and QPE SUM."""
    import numpy as np
    try:
        import xarray as xr
    except ImportError:
        print("need xarray+cfgrib for MRMS read (pip install xarray cfgrib)", file=sys.stderr); return
    refl_files = sorted(glob.glob(os.path.join(mrms_dir,"**","MergedReflectivityQCComposite*"),recursive=True))
    qpe_files  = sorted(glob.glob(os.path.join(mrms_dir,"**","RadarOnly_QPE_01H*"),recursive=True))
    print(f"MRMS: {len(refl_files)} refl, {len(qpe_files)} qpe files")
    refl_max=None; lat=lon=None
    for f in refl_files:
        ds=xr.open_dataset(f, engine="cfgrib", backend_kwargs={"indexpath":""})
        v=ds[list(ds.data_vars)[0]].values
        v=np.where(v<-90, np.nan, v)   # MRMS missing = -999/-99
        refl_max=v if refl_max is None else np.fmax(refl_max,v)
        if lat is None: lat=ds.latitude.values; lon=ds.longitude.values
        ds.close()
    qpe_total=None
    for f in qpe_files:   # 1-hour QPE; sum the hourly stamps for storm total
        ds=xr.open_dataset(f, engine="cfgrib", backend_kwargs={"indexpath":""})
        v=ds[list(ds.data_vars)[0]].values
        v=np.where(v<0,0.0,v)
        qpe_total=v if qpe_total is None else (qpe_total+v)
        ds.close()
    flds={}
    if refl_max is not None: flds["refl_max"]=refl_max
    if qpe_total is not None: flds["qpe_total"]=qpe_total
    if lat is not None: flds["lat"]=lat; flds["lon"]=lon
    _write_nc(out, flds, "MRMS event aggregates (refl max, QPE sum)")
    if refl_max is not None: print(f"MRMS refl_max: max={np.nanmax(refl_max):.1f} dBZ")
    if qpe_total is not None: print(f"MRMS qpe_total: max={np.nanmax(qpe_total):.1f} mm")

def asos_aggregate(asos_dir, out="asos_event.csv"):
    """Per-station event peaks from the IEM CSVs."""
    rows=[]
    KT2MS=0.514444
    for fn in sorted(glob.glob(os.path.join(asos_dir,"ASOS_*.csv"))):
        st=os.path.basename(fn)[5:-4]
        lat=lon=None; peak_wind=peak_gust=0.0; tmin=1e9; tmax=-1e9; pmax=0.0
        with open(fn) as f:
            for r in csv.DictReader(f):
                try: lat=float(r["lat"]); lon=float(r["lon"])
                except: pass
                for src,key in [("sknt","peak_wind"),("gust","peak_gust")]:
                    try:
                        v=float(r[src])*KT2MS
                        if key=="peak_wind": peak_wind=max(peak_wind,v)
                        else: peak_gust=max(peak_gust,v)
                    except: pass
                try:
                    t=(float(r["tmpf"])-32)*5/9; tmin=min(tmin,t); tmax=max(tmax,t)
                except: pass
                try: pmax=max(pmax,float(r["p01i"])*25.4)
                except: pass
        rows.append(dict(station=st,lat=lat,lon=lon,
                         peak_wind_ms=round(peak_wind,2),peak_gust_ms=round(peak_gust,2),
                         tmin_C=None if tmin>1e8 else round(tmin,1),
                         tmax_C=None if tmax<-1e8 else round(tmax,1),
                         max_1h_precip_mm=round(pmax,1)))
    with open(out,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"ASOS event peaks -> {out} ({len(rows)} stations)")
    for r in rows: print(f"  {r['station']}: gust={r['peak_gust_ms']} m/s wind={r['peak_wind_ms']} m/s")

def _stepkey(p):
    import re
    m=re.findall(r"(\d+)", os.path.basename(p)); return int(m[-1]) if m else 0

def _write_nc(out, fields, title):
    import netCDF4 as nc, numpy as np
    if not fields: print("no fields to write", file=sys.stderr); return
    d=nc.Dataset(out,"w"); d.title=title
    any2d=next(iter(fields.values()))
    ny,nx=any2d.shape
    d.createDimension("y",ny); d.createDimension("x",nx)
    for k,v in fields.items():
        var=d.createVariable(k,"f4",("y","x"),zlib=True)
        var[:]=np.asarray(v,dtype="f4")
    d.close(); print(f"wrote {out}: {list(fields.keys())}")

if __name__=="__main__":
    if len(sys.argv)<2:
        print(__doc__); sys.exit(1)
    cmd=sys.argv[1]
    if cmd=="erf":    erf_aggregate(sys.argv[2], sys.argv[3])
    elif cmd=="mrms": mrms_aggregate(sys.argv[2])
    elif cmd=="asos": asos_aggregate(sys.argv[2])
    else: print("unknown cmd",cmd); sys.exit(1)
