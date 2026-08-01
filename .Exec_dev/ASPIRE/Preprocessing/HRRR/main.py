import sys
import os
import argparse
from mpi4py import MPI
import glob
from pyproj import CRS, Transformer
from numpy import *
import time

from DownloadHRRR3DData import DownloadHRRR3DData
from ReadHRRR3DData import ReadHRRR3DData

def CreateLCCMapping(area):

    lat1 = area[2]
    lat2 = area[0]
    lon1 = area[1]
    lon2 = area[3]

    # Build CRS
    delta = lat2 - lat1
    lon0 = (lon1 + lon2) / 2
    lat0 = (lat1 + lat2) / 2

    lat_1 = lat1 + delta/6
    lat_2 = lat2 - delta/6

    lambert_conformal = (
        f"+proj=lcc +lat_1={lat_1:.6f} +lat_2={lat_2:.6f} "
        f"+lat_0={lat0:.6f} +lon_0={lon0:.6f} +datum=WGS84 +units=m +no_defs"
    )

    return lambert_conformal

if __name__ == "__main__":

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    start_time = time.time()

    if len(sys.argv) == 1:
        print("Usage: python3 WriteICFromERA5Data.py <input_filename> [--do_forecast=true]")
        sys.exit(1)

    parser = argparse.ArgumentParser(description="Download and process ERA5 data.")
    parser.add_argument("input_filename", help="Input filename, e.g. inputs_for_Laura")

    parser.add_argument(
        "--do_forecast", type=lambda x: x.lower() == "true", default=False,
        help="Set to true to download forecast data"
    )

    # Forecast args (only required if --do_forecast is set)
    parser.add_argument("--forecast_time_hours", type=int,
                        help="Forecast length in hours, e.g. 72")
    parser.add_argument("--interval_hours", type=int,
                        help="Forecast interval in hours, e.g. 6")
    parser.add_argument(
        "--init_hurricane_latlon",
        type=str,
        default=None,
        help="Optional initial hurricane lat,lon (e.g. --init_hurricane_latlon=25,-80)"
    )
    args = parser.parse_args()

    # Enforce requirement only if do_forecast is true
    if args.do_forecast:
        if args.forecast_time_hours is None or args.interval_hours is None:
            parser.error("--do_forecast requires --forecast_time_hours and --interval_hours")

    input_filename = args.input_filename
    do_forecast = args.do_forecast
    forecast_time_hours = args.forecast_time_hours
    interval_hours = args.interval_hours

    if rank == 0:
        print(f"Input file: {input_filename}")

    os.makedirs("HRRRFiles", exist_ok=True)
    os.makedirs("Output/VTK/3D/HRRRDomain", exist_ok=True)
    os.makedirs("Output/Thunderstorm", exist_ok=True)

    comm.Barrier();

    # Download 3d data over pressure levels
    if do_forecast:
        filenames, area = DownloadHRRR3DData(input_filename, forecast_time_hours, interval_hours)
        lambert_conformal = CreateLCCMapping(area)
        comm.Barrier();

        # Find all HRRR files
        filenames = sorted(glob.glob(os.path.join("HRRRFiles", "*.grib2")))

        # Distribute files across MPI ranks
        my_files = filenames[rank::size]

        for filename in my_files:
            print(f"[Rank {rank}] Processing file: {filename}")
            ReadHRRR3DData(filename, area, lambert_conformal)

    comm.Barrier()

    end_time = time.time()
    elapsed = end_time - start_time

    # Gather or print only on rank 0
    max_elapsed = comm.reduce(elapsed, op=MPI.MAX, root=0)  # max over all ranks

    if rank == 0:
        print(f"Total runtime (wall-clock, across ranks): {max_elapsed:.2f} seconds", flush=True)
