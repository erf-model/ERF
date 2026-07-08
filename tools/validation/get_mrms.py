#!/usr/bin/env python3
"""Download MRMS (Multi-Radar Multi-Sensor) products for the Sep-2024 cyclone.

Source: NOAA Open Data on AWS, bucket `noaa-mrms-pds` (public, no auth). We use
the bucket's HTTPS endpoint + the S3 ListObjectsV2 REST API (XML), so only
`requests` is needed -- no boto3.

Products (CONUS, 1 km grid), gzipped GRIB2 (.grib2.gz):
  - MergedReflectivityQCComposite_00.50  -- composite reflectivity  (vs ERF max_reflectivity)
  - PrecipRate_00.00                      -- instantaneous rate     (vs ERF rain rate)
  - RadarOnly_QPE_01H_00.00               -- 1-hour QPE accumulation (vs ERF rain_accum)

Cadence is 2 min; to keep volume sane we subsample to one file per
DEFAULT_STRIDE_MIN minutes (default 30) over the event window. Pass --all for
every file. Domain subsetting happens later at read time (GRIB2 is CONUS-wide;
~a few MB each composite).

Output: ./mrms/<PRODUCT>/<files>.grib2.gz  (idempotent).
Usage:  python get_mrms.py [outdir] [--all]
"""
import os, sys, time, datetime as dt
import xml.etree.ElementTree as ET
import requests

BUCKET = "noaa-mrms-pds"
BASE   = f"https://{BUCKET}.s3.amazonaws.com"
PRODUCTS = [
    "CONUS/MergedReflectivityQCComposite_00.50",
    "CONUS/PrecipRate_00.00",
    "CONUS/RadarOnly_QPE_01H_00.00",
]
START = dt.datetime(2024, 9, 17, 0, 0)
END   = dt.datetime(2024, 9, 24, 0, 0)
DEFAULT_STRIDE_MIN = 30

def list_keys(prefix):
    """List object keys under prefix via the S3 REST ListObjectsV2 XML API."""
    keys = []
    token = None
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
    """Parse YYYYMMDD-HHMMSS from an MRMS key; None if not parseable."""
    base = os.path.basename(key)
    for tok in base.replace(".", "_").split("_"):
        if len(tok) == 15 and tok[8] == "-" and tok[:8].isdigit():
            try:
                return dt.datetime.strptime(tok, "%Y%m%d-%H%M%S")
            except ValueError:
                return None
    return None

def fetch(key, outdir):
    rel = key.split("CONUS/", 1)[-1]
    fn = os.path.join(outdir, rel)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    if os.path.exists(fn) and os.path.getsize(fn) > 100:
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
    outdir = args[0] if args else "mrms"
    stride = dt.timedelta(minutes=(0 if take_all else DEFAULT_STRIDE_MIN))
    os.makedirs(outdir, exist_ok=True)
    print(f"MRMS -> {outdir}  {START} .. {END} UTC  stride={stride or 'all'}")

    # MRMS keys are organized by product/day: CONUS/<product>/<YYYYMMDD>/...
    for prod in PRODUCTS:
        got = 0
        d = START.date()
        next_t = START
        while d <= END.date():
            prefix = f"{prod}/{d.strftime('%Y%m%d')}/"
            try:
                keys = list_keys(prefix)
            except Exception as e:
                print(f"WARN list {prefix}: {e}", flush=True); keys = []
            keys = sorted(k for k in keys if k.endswith(".grib2.gz"))
            for k in keys:
                t = key_time(k)
                if t is None or t < START or t > END:
                    continue
                if stride and t < next_t:
                    continue
                if fetch(k, outdir):
                    got += 1
                if stride:
                    next_t = t + stride
            d += dt.timedelta(days=1)
        print(f"  {prod}: {got} files", flush=True)
    print("MRMS_DONE", flush=True)

if __name__ == "__main__":
    main()
