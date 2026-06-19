#!/usr/bin/env python3
"""Compare ERF event-aggregates against NEXRAD / MRMS / ASOS observations.

Everything is reduced to EVENT AGGREGATES (cadence-robust), then obs are put on
the ERF grid (gridded sources) or ERF is sampled at station points (ASOS):

  peak 10 m wind   ERF wspd10max      vs  ASOS peak gust / 2-min wind
  storm total rain ERF rain_accum     vs  MRMS RadarOnly_QPE_01H summed
  max reflectivity ERF cmpref_max     vs  MRMS composite max / NEXRAD composite max

The ERF grid geography comes from the plt2d lat_m/lon_m fields, so regridding is
just a nearest-neighbour (obs -> ERF cell) via scipy.cKDTree -- a mild coarsening
from MRMS 1 km / NEXRAD polar to ERF 3/9 km, which is appropriate for maxima.

Subcommands (run from a dir, pass paths):
  erf-grid   <plt2d_glob> <plt_glob>          -> erf_grid.npz   (lat,lon,wspd10max,rain_accum,refl_max)
  mrms-grid  <mrms_dir> <erf_grid.npz>        -> mrms_on_erf.npz (refl_max, qpe_total on ERF grid)
  nexrad-grid<nexrad_dir> <erf_grid.npz>      -> nexrad_on_erf.npz (per-site composite refl_max on ERF grid)
  asos-pts   <asos_event.csv> <erf_grid.npz>  -> asos_compare.csv (obs peak + ERF wspd10max at station cell)
  plot       <erf_grid.npz> [mrms_on_erf.npz] [nexrad_on_erf.npz] [asos_compare.csv] -> panels PNG

Heavy obs reads (cfgrib/pyart) only happen in the *-grid subcommands.
"""
import sys, os, glob, csv, re
import numpy as np

CLIP = 9   # boundary relaxation cells to mask in ERF fields

# ---------------------------------------------------------------- ERF side
def _stepkey(p):
    m = re.findall(r"(\d+)", os.path.basename(p)); return int(m[-1]) if m else 0

def erf_grid(plt2d_glob, plt_glob, out="erf_grid.npz"):
    import yt; yt.set_log_level(50)
    p2 = sorted(glob.glob(plt2d_glob), key=_stepkey)
    p3 = sorted(glob.glob(plt_glob),  key=_stepkey)
    if not p2: print("no plt2d", file=sys.stderr); return
    last2 = yt.load(p2[-1]); n2=[f[1] for f in last2.field_list if f[0]=="boxlib"]
    cg = last2.covering_grid(0, last2.domain_left_edge, last2.domain_dimensions)
    lat = np.asarray(cg["boxlib","lat_m"])[:,:,0]
    lon = np.asarray(cg["boxlib","lon_m"])[:,:,0]
    out_d = dict(lat=lat, lon=lon)
    if "wspd10max" in n2: out_d["wspd10max"]=np.asarray(cg["boxlib","wspd10max"])[:,:,0]
    if "cmpref_max" in n2: out_d["refl_max"]=np.asarray(cg["boxlib","cmpref_max"])[:,:,0]
    if p3:
        last3=yt.load(p3[-1]); n3=[f[1] for f in last3.field_list if f[0]=="boxlib"]
        cg3=last3.covering_grid(0,last3.domain_left_edge,last3.domain_dimensions)
        if "rain_accum" in n3: out_d["rain_accum"]=np.asarray(cg3["boxlib","rain_accum"])[:,:,0]
        if "refl_max" not in out_d and "max_reflectivity" in n3:
            acc=None
            for p in p3:
                ds=yt.load(p); g=ds.covering_grid(0,ds.domain_left_edge,ds.domain_dimensions)
                m=np.asarray(g["boxlib","max_reflectivity"])[:,:,0]
                acc=m if acc is None else np.maximum(acc,m)
            out_d["refl_max"]=acc
    # mask boundary ring
    mask=np.zeros_like(lat,dtype=bool); mask[CLIP:-CLIP,CLIP:-CLIP]=True
    out_d["interior"]=mask
    np.savez(out, **out_d)
    print(f"wrote {out}: {[k for k in out_d if k not in ('lat','lon','interior')]}")
    for k in out_d:
        if k in ("lat","lon","interior"): continue
        a=out_d[k][mask]; print(f"  ERF {k}: min/mean/max={np.nanmin(a):.2f}/{np.nanmean(a):.2f}/{np.nanmax(a):.2f}")

# ------------------------------------------------------- nearest-neighbour regrid
def _to_erf(src_lat, src_lon, src_val, erf_lat, erf_lon, agg="max"):
    """Bin source points into ERF cells via KDTree on the ERF grid; aggregate."""
    from scipy.spatial import cKDTree
    # ERF grid as target points (lon in -180..180)
    elon = np.where(erf_lon>180, erf_lon-360, erf_lon)
    tree = cKDTree(np.column_stack([erf_lat.ravel(), elon.ravel()]))
    slon = np.where(src_lon>180, src_lon-360, src_lon)
    pts = np.column_stack([src_lat.ravel(), slon.ravel()])
    val = src_val.ravel()
    good = np.isfinite(val)
    d,idx = tree.query(pts[good])
    out = np.full(erf_lat.size, np.nan)
    flat_idx = idx
    v = val[good]
    # aggregate by target cell
    order = np.argsort(flat_idx)
    fi = flat_idx[order]; vv=v[order]
    uniq, start = np.unique(fi, return_index=True)
    for u,s,e in zip(uniq, start, list(start[1:])+[len(fi)]):
        seg=vv[s:e]
        out[u]= np.nanmax(seg) if agg=="max" else np.nansum(seg)
    return out.reshape(erf_lat.shape)

def _subset_box(lat, lon, blat, blon, pad=0.5):
    lo = np.where(lon>180, lon-360, lon)
    m = (lat>=blat[0]-pad)&(lat<=blat[1]+pad)&(lo>=blon[0]-pad)&(lo<=blon[1]+pad)
    return m

def mrms_grid(mrms_dir, erf_npz, out="mrms_on_erf.npz"):
    import xarray as xr
    g=np.load(erf_npz); elat=g["lat"]; elon=g["lon"]
    blat=(elat.min(),elat.max()); blon=(np.where(elon>180,elon-360,elon).min(),
                                         np.where(elon>180,elon-360,elon).max())
    refl_files=sorted(f for f in glob.glob(os.path.join(mrms_dir,"**","*MergedReflectivityQCComposite*"),recursive=True)
                      if os.path.isfile(f) and f.endswith((".grib2",".grib2.gz")))
    qpe_files =sorted(f for f in glob.glob(os.path.join(mrms_dir,"**","*RadarOnly_QPE_01H*"),recursive=True)
                      if os.path.isfile(f) and f.endswith((".grib2",".grib2.gz")))
    print(f"MRMS: {len(refl_files)} refl, {len(qpe_files)} qpe")
    def _open(f):
        import gzip,shutil,tempfile
        if f.endswith(".gz"):
            tmp=tempfile.NamedTemporaryFile(suffix=".grib2",delete=False).name
            with gzip.open(f,"rb") as i, open(tmp,"wb") as o: shutil.copyfileobj(i,o)
            return tmp,True
        return f,False
    refl_max=None; sub=None; slat=slon=None
    for f in refl_files:
        fp,tmp=_open(f)
        try:
            ds=xr.open_dataset(fp,engine="cfgrib",backend_kwargs={"indexpath":""})
            v=ds[list(ds.data_vars)[0]].values; v=np.where(v<-90,np.nan,v)
            la=ds.latitude.values; lo=ds.longitude.values
            La,Lo=np.meshgrid(la,lo,indexing="ij") if la.ndim==1 else (la,lo)
            if sub is None:
                sub=_subset_box(La,Lo,blat,blon); slat=La[sub]; slon=Lo[sub]
            cur=v[sub]
            refl_max=cur if refl_max is None else np.fmax(refl_max,cur)
            ds.close()
        finally:
            if tmp: os.unlink(fp)
    # QPE: RadarOnly_QPE_01H is a 1-hour ROLLING accumulation refreshed every
    # ~2 min. Summing every file would multiply-count the same hour ~30x. Take
    # ONE file per top-of-hour (nearest stamp) and sum those -> true storm total.
    import re as _re, datetime as _dt
    def _qpe_time(p):
        b=os.path.basename(p)
        for tok in b.replace(".","_").split("_"):
            if len(tok)==15 and tok[8]=="-" and tok[:8].isdigit():
                try: return _dt.datetime.strptime(tok,"%Y%m%d-%H%M%S")
                except ValueError: return None
        return None
    by_hour={}
    for f in qpe_files:
        t=_qpe_time(f)
        if t is None: continue
        hr=t.replace(minute=0,second=0,microsecond=0)
        # keep the file closest to the top of the hour
        off=abs((t-hr).total_seconds()); off=min(off, 3600-off)
        if hr not in by_hour or off < by_hour[hr][0]:
            by_hour[hr]=(off,f)
    hourly_qpe=[v[1] for k,v in sorted(by_hour.items())]
    print(f"MRMS QPE: {len(qpe_files)} files -> {len(hourly_qpe)} hourly stamps summed")
    qpe_total=None
    for f in hourly_qpe:
        fp,tmp=_open(f)
        try:
            ds=xr.open_dataset(fp,engine="cfgrib",backend_kwargs={"indexpath":""})
            v=ds[list(ds.data_vars)[0]].values; v=np.where(v<0,0.0,v)
            qcur=v[sub]
            qpe_total=qcur if qpe_total is None else (qpe_total+qcur)
            ds.close()
        finally:
            if tmp: os.unlink(fp)
    out_d=dict(lat=elat,lon=elon)
    if refl_max is not None: out_d["refl_max"]=_to_erf(slat,slon,refl_max,elat,elon,"max")
    if qpe_total is not None: out_d["qpe_total"]=_to_erf(slat,slon,qpe_total,elat,elon,"max")
    np.savez(out,**out_d)
    print(f"wrote {out}")
    for k in ("refl_max","qpe_total"):
        if k in out_d: print(f"  MRMS->ERF {k}: max={np.nanmax(out_d[k]):.1f}")

def nexrad_grid(nexrad_dir, erf_npz, out="nexrad_on_erf.npz"):
    import pyart
    g=np.load(erf_npz); elat=g["lat"]; elon=g["lon"]
    blat=(elat.min(),elat.max())
    elon180=np.where(elon>180,elon-360,elon); blon=(elon180.min(),elon180.max())
    sites=sorted([d for d in os.listdir(nexrad_dir) if os.path.isdir(os.path.join(nexrad_dir,d))])
    print(f"NEXRAD sites: {sites}")
    refl_max=np.full(elat.shape,np.nan)
    for site in sites:
        vols=sorted(glob.glob(os.path.join(nexrad_dir,site,"*")))
        print(f"  {site}: {len(vols)} volumes")
        for vf in vols:
            try:
                radar=pyart.io.read_nexrad_archive(vf)
                # composite (column-max) reflectivity on the radar's gates
                sweep_refl=radar.fields["reflectivity"]["data"]
                lats,lons,_=radar.get_gate_lat_lon_alt(0)  # approx; per-sweep below
                # build composite over all sweeps by max at gate lat/lon
                allval=[]; allla=[]; alllo=[]
                for s in range(radar.nsweeps):
                    la,lo,_=radar.get_gate_lat_lon_alt(s)
                    d=radar.get_field(s,"reflectivity")
                    allval.append(np.asarray(d).ravel()); allla.append(la.ravel()); alllo.append(lo.ravel())
                val=np.ma.filled(np.concatenate(allval),np.nan)
                la=np.concatenate(allla); lo=np.concatenate(alllo)
                m=_subset_box(la,lo,blat,blon,pad=0.1)
                if m.sum()==0: continue
                gx=_to_erf(la[m],lo[m],val[m],elat,elon,"max")
                refl_max=np.fmax(refl_max,gx)
            except Exception as e:
                print(f"    skip {os.path.basename(vf)}: {e}", file=sys.stderr)
    np.savez(out, lat=elat, lon=elon, refl_max=refl_max)
    print(f"wrote {out}: NEXRAD refl_max max={np.nanmax(refl_max):.1f}")

# ---------------------------------------------------------------- ASOS points
def asos_pts(asos_event_csv, erf_npz, out="asos_compare.csv"):
    from scipy.spatial import cKDTree
    g=np.load(erf_npz); elat=g["lat"]; elon=np.where(g["lon"]>180,g["lon"]-360,g["lon"])
    have_w = "wspd10max" in g.files
    wmax = g["wspd10max"] if have_w else None
    tree=cKDTree(np.column_stack([elat.ravel(),elon.ravel()]))
    rows=[]
    with open(asos_event_csv) as f:
        for r in csv.DictReader(f):
            if not r.get("lat") or r["lat"]=="None": continue
            la=float(r["lat"]); lo=float(r["lon"])
            d,idx=tree.query([la,lo])
            erf_w = float(wmax.ravel()[idx]) if have_w else None
            rows.append(dict(station=r["station"],lat=la,lon=lo,
                             obs_peak_gust_ms=r.get("peak_gust_ms"),
                             obs_peak_wind_ms=r.get("peak_wind_ms"),
                             erf_wspd10max_ms=None if erf_w is None else round(erf_w,2)))
    with open(out,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"wrote {out} ({len(rows)} stations)")
    print(f"{'stn':5s} {'obs_gust':>9s} {'obs_wind':>9s} {'erf_w10max':>11s}")
    for r in rows:
        print(f"{r['station']:5s} {str(r['obs_peak_gust_ms']):>9s} {str(r['obs_peak_wind_ms']):>9s} {str(r['erf_wspd10max_ms']):>11s}")

# ---------------------------------------------------------------- plot
def plot(erf_npz, *extras):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    g=np.load(erf_npz); elat=g["lat"]; elon=np.where(g["lon"]>180,g["lon"]-360,g["lon"])
    mrms=nex=asos=None
    for e in extras:
        if e.endswith(".npz") and "mrms" in e: mrms=np.load(e)
        elif e.endswith(".npz") and "nexrad" in e: nex=np.load(e)
        elif e.endswith(".csv"): asos=e
    # reflectivity comparison row
    panels=[("ERF refl_max", g.get("refl_max"))]
    if mrms is not None and "refl_max" in mrms.files: panels.append(("MRMS refl_max", mrms["refl_max"]))
    if nex  is not None and "refl_max" in nex.files:  panels.append(("NEXRAD refl_max", nex["refl_max"]))
    n=len(panels)
    fig,axs=plt.subplots(1,n,figsize=(5*n,4.4),constrained_layout=True,squeeze=False)
    for ax,(title,arr) in zip(axs[0],panels):
        if arr is None: continue
        im=ax.pcolormesh(elon,elat,np.where(np.isfinite(arr),arr,np.nan),
                         vmin=-10,vmax=60,cmap="turbo",shading="auto")
        ax.set_title(title); fig.colorbar(im,ax=ax,label="dBZ",shrink=0.8)
    fig.suptitle("Event-max composite reflectivity: ERF vs radar")
    fig.savefig("compare_reflmax.png",dpi=120,bbox_inches="tight"); print("wrote compare_reflmax.png")
    # ASOS scatter
    if asos:
        import csv as _csv
        st=[];og=[];ew=[]
        for r in _csv.DictReader(open(asos)):
            if r["erf_wspd10max_ms"] in (None,"None",""): continue
            st.append(r["station"]); og.append(float(r["obs_peak_gust_ms"])); ew.append(float(r["erf_wspd10max_ms"]))
        fig2,ax=plt.subplots(figsize=(5.5,5.5))
        ax.scatter(og,ew)
        for s,x,y in zip(st,og,ew): ax.annotate(s,(x,y),fontsize=7)
        lim=max(max(og),max(ew))*1.1; ax.plot([0,lim],[0,lim],"k--",lw=1)
        ax.set_xlabel("ASOS peak gust [m/s]"); ax.set_ylabel("ERF wspd10max [m/s]")
        ax.set_title("Peak 10 m wind: ERF vs ASOS"); ax.set_xlim(0,lim); ax.set_ylim(0,lim)
        fig2.savefig("compare_asos_wind.png",dpi=120,bbox_inches="tight"); print("wrote compare_asos_wind.png")

def refl_snapshots(plt_glob, mrms_dir, erf_grid_npz, outdir="refl_snapshots",
                   start="2024-09-17_00:00:00", tol_min=5):
    """Time-matched ERF max_reflectivity vs MRMS composite, per ERF plotfile.

    For each ERF plt (max_reflectivity, a column-max composite = MRMS analog):
      - convert its current_time(s)+start -> absolute UTC
      - find the MRMS MergedReflectivityQCComposite frame nearest that UTC (<= tol_min)
      - regrid that ONE MRMS frame to the ERF grid (cKDTree NN)
      - write a 3-panel PNG (ERF | MRMS | ERF-MRMS diff) with coastlines
    Also writes refl_timeseries.csv of domain stats per matched time.
    """
    import yt, datetime as dt, glob as _glob, os as _os, csv as _csv
    import xarray as xr
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    import cartopy.crs as ccrs, cartopy.feature as cfeature
    from scipy.spatial import cKDTree
    yt.set_log_level(50)
    _os.makedirs(outdir, exist_ok=True)
    t0 = dt.datetime.strptime(start, "%Y-%m-%d_%H:%M:%S")

    g = np.load(erf_grid_npz)
    elat = g["lat"]; elon = np.where(g["lon"]>180, g["lon"]-360, g["lon"])
    interior = g["interior"] if "interior" in g.files else np.ones_like(elat,bool)
    tree = cKDTree(np.column_stack([elat.ravel(), elon.ravel()]))

    # index MRMS composite files by UTC
    mrms = sorted(f for f in _glob.glob(_os.path.join(mrms_dir,"**","*MergedReflectivityQCComposite*"),recursive=True)
                  if _os.path.isfile(f) and f.endswith((".grib2",".grib2.gz")))
    def _mt(p):
        b=_os.path.basename(p)
        for tok in b.replace(".","_").split("_"):
            if len(tok)==15 and tok[8]=="-" and tok[:8].isdigit():
                try: return dt.datetime.strptime(tok,"%Y%m%d-%H%M%S")
                except ValueError: return None
        return None
    mrms_t = [(_mt(f),f) for f in mrms]; mrms_t=[(t,f) for t,f in mrms_t if t]
    print(f"{len(mrms_t)} MRMS composite frames indexed")

    def _open(f):
        import gzip,shutil,tempfile
        if f.endswith(".gz"):
            tmp=tempfile.NamedTemporaryFile(suffix=".grib2",delete=False).name
            with gzip.open(f,"rb") as i, open(tmp,"wb") as o: shutil.copyfileobj(i,o)
            return tmp,True
        return f,False

    ext=[elon.min(),elon.max(),elat.min(),elat.max()]; proj=ccrs.PlateCarree()
    def deco(ax):
        ax.set_extent(ext,crs=proj); ax.add_feature(cfeature.COASTLINE,linewidth=0.6)
        ax.add_feature(cfeature.STATES,linewidth=0.25,edgecolor="gray")
    rows=[]
    plts=sorted(_glob.glob(plt_glob), key=lambda p:int(''.join(filter(str.isdigit,_os.path.basename(p)))))
    sub=None; slat=slon=None
    for p in plts:
        ds=yt.load(p); n=[f[1] for f in ds.field_list if f[0]=="boxlib"]
        if "max_reflectivity" not in n: continue
        utc=t0+dt.timedelta(seconds=float(ds.current_time))
        # nearest MRMS
        best=min(mrms_t, key=lambda tf: abs((tf[0]-utc).total_seconds()))
        dmin=abs((best[0]-utc).total_seconds())/60.0
        if dmin>tol_min:
            print(f"{_os.path.basename(p)} {utc} -> no MRMS within {tol_min}min (nearest {dmin:.1f})"); continue
        cg=ds.covering_grid(0,ds.domain_left_edge,ds.domain_dimensions)
        erf=np.asarray(cg["boxlib","max_reflectivity"])[:,:,0]
        erf=np.where(interior,erf,np.nan)
        # MRMS frame -> ERF grid
        fp,tmp=_open(best[1])
        try:
            md=xr.open_dataset(fp,engine="cfgrib",backend_kwargs={"indexpath":""})
            v=md[list(md.data_vars)[0]].values; v=np.where(v<-90,np.nan,v)
            la=md.latitude.values; lo=md.longitude.values
            La,Lo=(np.meshgrid(la,lo,indexing="ij") if la.ndim==1 else (la,lo))
            if sub is None:
                lo180=np.where(Lo>180,Lo-360,Lo)
                sub=((La>=elat.min()-0.5)&(La<=elat.max()+0.5)&
                     (lo180>=elon.min()-0.5)&(lo180<=elon.max()+0.5))
                slat=La[sub]; slon=np.where(Lo>180,Lo-360,Lo)[sub]
            d,idx=tree.query(np.column_stack([slat,slon]))
            mg=np.full(elat.size,np.nan); vv=v[sub]
            order=np.argsort(idx); fi=idx[order]; vs=vv[order]
            uniq,st=np.unique(fi,return_index=True)
            for u,s,e in zip(uniq,st,list(st[1:])+[len(fi)]):
                seg=vs[s:e]; mg[u]=np.nanmax(seg) if np.any(np.isfinite(seg)) else np.nan
            mrg=mg.reshape(elat.shape); mrg=np.where(interior,mrg,np.nan)
            md.close()
        finally:
            if tmp: _os.unlink(fp)
        # stats over cells where MRMS has echo (>5 dBZ)
        m=interior & np.isfinite(mrg) & np.isfinite(erf)
        thr=m & (mrg>5)
        bias=float(np.nanmean((erf-mrg)[thr])) if thr.sum() else float("nan")
        rows.append(dict(time=utc.isoformat(), plt=_os.path.basename(p),
                         mrms_dt_min=round(dmin,1),
                         erf_maxdbz=round(float(np.nanmax(erf)),1),
                         mrms_maxdbz=round(float(np.nanmax(mrg)),1),
                         bias_dbz_echo=round(bias,1)))
        # 3-panel
        fig,axs=plt.subplots(1,3,figsize=(17,5),subplot_kw={"projection":proj},constrained_layout=True)
        for ax in axs: deco(ax)
        im0=axs[0].pcolormesh(elon,elat,erf,vmin=-10,vmax=60,cmap="turbo",shading="auto",transform=proj)
        axs[0].set_title(f"ERF max_reflectivity\n{utc:%Y-%m-%d %H:%M}Z"); fig.colorbar(im0,ax=axs[0],shrink=0.7,label="dBZ")
        im1=axs[1].pcolormesh(elon,elat,mrg,vmin=-10,vmax=60,cmap="turbo",shading="auto",transform=proj)
        axs[1].set_title(f"MRMS composite\n{best[0]:%H:%M}Z ({dmin:.0f} min off)"); fig.colorbar(im1,ax=axs[1],shrink=0.7,label="dBZ")
        diff=erf-mrg
        im2=axs[2].pcolormesh(elon,elat,diff,vmin=-30,vmax=30,cmap="RdBu_r",shading="auto",transform=proj)
        axs[2].set_title("ERF - MRMS [dBZ]"); fig.colorbar(im2,ax=axs[2],shrink=0.7,label="dBZ")
        fn=_os.path.join(outdir,f"refl_{utc:%Y%m%d_%H%M}.png")
        fig.savefig(fn,dpi=110,bbox_inches="tight"); plt.close(fig)
        print(f"{_os.path.basename(p)} {utc:%m-%d %H:%M}Z  ERFmax={rows[-1]['erf_maxdbz']} MRMSmax={rows[-1]['mrms_maxdbz']} bias={rows[-1]['bias_dbz_echo']}")
    if rows:
        with open(_os.path.join(outdir,"refl_timeseries.csv"),"w",newline="") as f:
            w=_csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
        print(f"wrote {len(rows)} snapshots + refl_timeseries.csv to {outdir}/")

def wind_panels(erf_npz, asos_event_csv, out="compare_wind_panels.png", which="gust"):
    """Two-panel map: (1) ASOS station points colored by peak wind,
    (2) ERF event-max wspd10max map. Shared colorbar, CONUS coastline."""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    import cartopy.crs as ccrs, cartopy.feature as cfeature
    g=np.load(erf_npz)
    elat=g["lat"]; elon=np.where(g["lon"]>180,g["lon"]-360,g["lon"])
    w=g["wspd10max"].copy()
    if "interior" in g.files:  # mask boundary relaxation ring
        w=np.where(g["interior"], w, np.nan)
    # ASOS points
    col = "peak_gust_ms" if which=="gust" else "peak_wind_ms"
    st=[]; la=[]; lo=[]; val=[]
    for r in csv.DictReader(open(asos_event_csv)):
        if not r.get("lat") or r["lat"]=="None": continue
        try: v=float(r[col])
        except: continue
        st.append(r["station"]); la.append(float(r["lat"])); lo.append(float(r["lon"])); val.append(v)
    la=np.array(la); lo=np.array(lo); val=np.array(val)
    # shared colour range over both obs points and ERF interior
    vmax=np.nanpercentile(np.concatenate([val, w[np.isfinite(w)]]), 99)
    vmax=max(vmax, val.max()); vmin=0.0
    ext=[elon.min(), elon.max(), elat.min(), elat.max()]
    proj=ccrs.PlateCarree()
    fig,axs=plt.subplots(1,2,figsize=(15,6.2),subplot_kw={"projection":proj},
                         constrained_layout=True)
    def deco(ax):
        ax.set_extent(ext,crs=proj)
        ax.add_feature(cfeature.COASTLINE,linewidth=0.7)
        ax.add_feature(cfeature.STATES,linewidth=0.3,edgecolor="gray")
        ax.add_feature(cfeature.BORDERS,linewidth=0.4)
        gl=ax.gridlines(draw_labels=True,linewidth=0.2,alpha=0.3)
        gl.top_labels=False; gl.right_labels=False
    # Panel 1: ASOS points
    deco(axs[0])
    sc=axs[0].scatter(lo,la,c=val,cmap="turbo",vmin=vmin,vmax=vmax,
                      s=70,edgecolor="k",linewidth=0.5,transform=proj,zorder=5)
    # Label stations only when the network is sparse enough to read.
    if len(st) <= 30:
        for s,x,y in zip(st,lo,la): axs[0].annotate(s,(x,y),fontsize=7,
                          xytext=(3,3),textcoords="offset points",transform=proj,zorder=6)
    axs[0].set_title(f"ASOS observed peak {('gust' if which=='gust' else '2-min wind')} [m/s]")
    fig.colorbar(sc,ax=axs[0],shrink=0.8,label="m/s")
    # Panel 2: ERF wspd10max map
    deco(axs[1])
    im=axs[1].pcolormesh(elon,elat,w,cmap="turbo",vmin=vmin,vmax=vmax,
                         shading="auto",transform=proj)
    # overlay station markers (outline only) for spatial reference
    axs[1].scatter(lo,la,facecolors="none",edgecolor="k",s=80,linewidth=0.7,
                   transform=proj,zorder=5)
    axs[1].set_title("ERF event-max 10 m wind (wspd10max) [m/s]")
    fig.colorbar(im,ax=axs[1],shrink=0.8,label="m/s")
    fig.suptitle("Peak 10 m wind: ASOS stations vs ERF event maximum",fontsize=13)
    fig.savefig(out,dpi=130,bbox_inches="tight")
    print(f"wrote {out}  (ASOS {which} max={val.max():.1f}, ERF interior max={np.nanmax(w):.1f} m/s)")

if __name__=="__main__":
    if len(sys.argv)<2: print(__doc__); sys.exit(1)
    c=sys.argv[1]
    if   c=="erf-grid":    erf_grid(sys.argv[2],sys.argv[3])
    elif c=="mrms-grid":   mrms_grid(sys.argv[2],sys.argv[3])
    elif c=="nexrad-grid": nexrad_grid(sys.argv[2],sys.argv[3])
    elif c=="asos-pts":    asos_pts(sys.argv[2],sys.argv[3])
    elif c=="wind-panels": wind_panels(sys.argv[2],sys.argv[3],
                                       *( [sys.argv[4]] if len(sys.argv)>4 else [] ))
    elif c=="refl-snapshots": refl_snapshots(sys.argv[2],sys.argv[3],sys.argv[4],
                                       *( [sys.argv[5]] if len(sys.argv)>5 else [] ))
    elif c=="plot":        plot(*sys.argv[2:])
    else: print("unknown cmd",c); sys.exit(1)
