#!/usr/bin/env python3
"""Download ASOS/AWOS surface observations for the Sep-2024 cyclone validation.

Source: Iowa Environmental Mesonet (IEM) ASOS archive -- open, no auth, returns
1-min/standard METAR obs as CSV. We pull the coastal/island stations inside (or
just outside) the WFIP3 domain (offshore Long Island / southern New England).

Variables retrieved (IEM names):
  sknt  -- wind speed [knots] (10 m)         -> convert to m/s
  drct  -- wind direction [deg]
  gust  -- wind gust [knots]                 -> convert to m/s   (ASOS peak)
  tmpf  -- 2 m air temperature [F]           -> convert to degC/K
  dwpf  -- 2 m dew point [F]
  mslp  -- mean sea-level pressure [hPa]
  p01i  -- 1-hour precip [inch]              -> convert to mm
  alti  -- altimeter [inHg]

These map to the ERF surface diagnostics: wspd10/u10/v10 (sknt/drct),
wspd10max (gust / running max), t2 (tmpf), q2 (dwpf), surf_pres (mslp).

Output: one CSV per station in ./asos/  (idempotent: skips existing non-empty).
Usage:  python get_asos.py [outdir]
"""
import os, sys, io, csv, time
import requests

# Domain bounding box (ERF 9km/3km WFIP3 grid, center 40.83N, -72.03W) plus a
# small halo, so we capture EVERY ASOS/AWOS station inside the domain rather than
# a hand-picked subset.
DOMAIN = dict(latmin=37.0, latmax=44.4, lonmin=-78.0, lonmax=-66.2)

# IEM networks that can overlap the NE-US / offshore domain. We query each
# network's station GeoJSON and keep stations whose coords fall in DOMAIN.
CANDIDATE_NETWORKS = ["NY_ASOS","NJ_ASOS","CT_ASOS","RI_ASOS","MA_ASOS","PA_ASOS",
                      "NH_ASOS","VT_ASOS","DE_ASOS","MD_ASOS","ME_ASOS"]

def discover_stations():
    """Return [(id, network, name)] for every station inside DOMAIN."""
    out = []
    for net in CANDIDATE_NETWORKS:
        try:
            r = requests.get(f"https://mesonet.agron.iastate.edu/geojson/network/{net}.geojson",
                             timeout=90)
            r.raise_for_status()
            feats = r.json()["features"]
        except Exception as e:
            print(f"  WARN network {net}: {e}", flush=True); continue
        for f in feats:
            lon, lat = f["geometry"]["coordinates"]
            if (DOMAIN["latmin"] <= lat <= DOMAIN["latmax"] and
                DOMAIN["lonmin"] <= lon <= DOMAIN["lonmax"]):
                out.append((f["id"], net, f["properties"].get("sname","")))
    out.sort()
    return out

START = ("2024", "09", "17")
END   = ("2024", "09", "24")  # inclusive day; IEM end is exclusive-ish, use 25 to be safe
END_SAFE = ("2024", "09", "25")

IEM = ("https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py")

# IEM throttles bursts; be polite between stations.
INTER_STATION_SLEEP = 8

def fetch(station, network, outdir):
    fn = os.path.join(outdir, f"ASOS_{station}.csv")
    if os.path.exists(fn) and os.path.getsize(fn) > 500:
        print(f"skip {fn}", flush=True); return
    params = {
        "station": station,
        "network": network,
        "data": "sknt,drct,gust,tmpf,dwpf,mslp,p01i,alti",
        "year1": START[0], "month1": START[1], "day1": START[2],
        "year2": END_SAFE[0], "month2": END_SAFE[1], "day2": END_SAFE[2],
        "tz": "UTC", "format": "onlycomma", "latlon": "yes",
        "missing": "M", "trace": "0.0001", "direct": "no",
    }
    for attempt in range(5):
        try:
            r = requests.get(IEM, params=params, timeout=120)
            r.raise_for_status()
            n = max(0, r.text.count("\n") - 1)
            if n == 0:
                # often transient throttling -> back off and retry
                print(f"WARN {station}: 0 rows (attempt {attempt}), backing off", flush=True)
                time.sleep(10 + 5*attempt); continue
            with open(fn, "w") as f:
                f.write(r.text)
            print(f"GET  {fn}  ({n} rows)", flush=True)
            return
        except Exception as e:
            print(f"retry {station} ({attempt}): {e}", flush=True); time.sleep(5)
    print(f"FAIL {station}", flush=True)

def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "asos"
    os.makedirs(outdir, exist_ok=True)
    print(f"ASOS -> {outdir}  ({START[0]}-{START[1]}-{START[2]} .. {END[0]}-{END[1]}-{END[2]} UTC)")
    stations = discover_stations()
    print(f"Discovered {len(stations)} stations inside domain "
          f"[{DOMAIN['latmin']},{DOMAIN['latmax']}]x[{DOMAIN['lonmin']},{DOMAIN['lonmax']}]")
    for st, net, _ in stations:
        fetch(st, net, outdir)
        time.sleep(INTER_STATION_SLEEP)
    print("ASOS_DONE", flush=True)

if __name__ == "__main__":
    main()
