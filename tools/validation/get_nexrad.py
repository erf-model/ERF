#!/usr/bin/env python3
"""Download NEXRAD WSR-88D Level II volumes for the Sep-2024 cyclone.

Source: NOAA Open Data on AWS, bucket `noaa-nexrad-level2` (public, no auth).
HTTPS + S3 ListObjectsV2 REST (XML) only -- no boto3.

Key layout: <YYYY>/<MM>/<DD>/<SITE>/<SITE><YYYYMMDD>_<HHMMSS>_V06[.gz]
Radars covering the WFIP3 domain (offshore Long Island / S. New England):
  KOKX -- Upton NY (primary, in-domain)
  KBOX -- Boston MA
  KDIX -- Mt Holly NJ / Philadelphia
  KGYX -- Gray ME (northern edge)

Each volume is ~few MB and ~5-10 min cadence; we subsample to one per
DEFAULT_STRIDE_MIN minutes (default 30). Pass --all for every volume.
Gridding to a Cartesian grid (for ERF comparison) is a later step and needs
Py-ART (`pip install arm-pyart`); this script only fetches the raw L2 files.

Output: ./nexrad/<SITE>/<files>  (idempotent).
Usage:  python get_nexrad.py [outdir] [--all] [--sites KOKX,KBOX]
"""
import os, sys, time, datetime as dt
import xml.etree.ElementTree as ET
import requests

# NOTE: the primary `noaa-nexrad-level2` PDS bucket denies anonymous REST
# listing (403) and anonymous GET. The Unidata mirror `unidata-nexrad-level2`
# holds the same key layout and DOES allow anonymous list + GET, so we use it.
BUCKET = "unidata-nexrad-level2"
BASE   = f"https://{BUCKET}.s3.amazonaws.com"
DEFAULT_SITES = ["KOKX", "KBOX", "KDIX", "KGYX"]
START = dt.datetime(2024, 9, 17, 0, 0)
END   = dt.datetime(2024, 9, 24, 0, 0)
DEFAULT_STRIDE_MIN = 30

def list_keys(prefix):
    keys, token = [], None
    ns = "{http://s3.amazonaws.com/doc/2006-03-01/}"
    while True:
        params = {"list-type": "2", "prefix": prefix}
        if token: params["continuation-token"] = token
        r = requests.get(BASE, params=params, timeout=120)
        r.raise_for_status()
        root = ET.fromstring(r.text)
        for c in root.findall(f"{ns}Contents"):
            keys.append(c.find(f"{ns}Key").text)
        trunc = root.find(f"{ns}IsTruncated")
        if trunc is not None and trunc.text == "true":
            token = root.find(f"{ns}NextContinuationToken").text
        else:
            break
    return keys

def key_time(key):
    base = os.path.basename(key)            # e.g. KOKX20240917_001234_V06
    try:
        stamp = base[4:19]                  # YYYYMMDD_HHMMSS
        return dt.datetime.strptime(stamp, "%Y%m%d_%H%M%S")
    except (ValueError, IndexError):
        return None

def fetch(key, outdir):
    fn = os.path.join(outdir, os.path.basename(key))
    if os.path.exists(fn) and os.path.getsize(fn) > 1000:
        return False
    for attempt in range(4):
        try:
            r = requests.get(f"{BASE}/{key}", timeout=180)
            r.raise_for_status()
            with open(fn, "wb") as f:
                f.write(r.content)
            return True
        except Exception as e:
            print(f"retry {key} ({attempt}): {e}", flush=True); time.sleep(5)
    print(f"FAIL {key}", flush=True)
    return False

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    take_all = "--all" in sys.argv
    sites = DEFAULT_SITES
    for a in sys.argv:
        if a.startswith("--sites"):
            sites = a.split("=", 1)[-1].split(",") if "=" in a else sites
    outdir = args[0] if args else "nexrad"
    stride = dt.timedelta(minutes=(0 if take_all else DEFAULT_STRIDE_MIN))
    os.makedirs(outdir, exist_ok=True)
    print(f"NEXRAD -> {outdir}  {START} .. {END} UTC  sites={sites}  stride={stride or 'all'}")

    for site in sites:
        sdir = os.path.join(outdir, site)
        os.makedirs(sdir, exist_ok=True)
        got = 0
        d = START.date()
        next_t = START
        while d <= END.date():
            prefix = f"{d.strftime('%Y/%m/%d')}/{site}/"
            try:
                keys = list_keys(prefix)
            except Exception as e:
                print(f"WARN list {prefix}: {e}", flush=True); keys = []
            keys = sorted(k for k in keys if k.endswith(("_V06", "_V06.gz", "_V08")))
            for k in keys:
                t = key_time(k)
                if t is None or t < START or t > END:
                    continue
                if stride and t < next_t:
                    continue
                if fetch(k, sdir):
                    got += 1
                if stride:
                    next_t = t + stride
            d += dt.timedelta(days=1)
        print(f"  {site}: {got} volumes", flush=True)
    print("NEXRAD_DONE", flush=True)

if __name__ == "__main__":
    main()
