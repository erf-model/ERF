import os
import requests
from datetime import datetime, timedelta
from mpi4py import MPI


def read_user_input(filename):
    data = {}

    with open(filename, "r") as f:
        for line in f:
            if ":" in line:
                key, value = line.strip().split(":", 1)

                key = key.strip().lower()
                value = value.strip()

                if key == "area":
                    data[key] = [float(x) for x in value.split(",")]
                elif key == "time":
                    data[key] = value
                else:
                    data[key] = int(value)

    return data


def generate_timestamps(start_dt, hours=72, interval=1):
    """
    Generate timestamps beginning at start_dt.

    Example:
        start_dt = 2025-06-17 00:00
        hours = 72
        interval = 1

    Returns:
        2025-06-17 00:00
        2025-06-17 01:00
        ...
        2025-06-20 00:00
    """
    timestamps = []

    for i in range(0, hours, interval):
        timestamps.append(start_dt + timedelta(hours=i))

    return timestamps


def build_hrrr_url(dt):
    """
    Construct HRRR archive URL.

    Example:
        dt = 2025-06-17 00:00

        https://pando-rgw01.chpc.utah.edu/hrrr/prs/20250617/
        hrrr.t00z.wrfsfcf00.grib2
    """

    ymd = dt.strftime("%Y%m%d")
    hour = dt.strftime("%H")

    url = (
        f"https://pando-rgw01.chpc.utah.edu/hrrr/sfc/"
        f"{ymd}/hrrr.t{hour}z.wrfsfcf01.grib2"
    )

    return url


def download_file(url, output_filename, idx):
    if os.path.exists(output_filename):
        print(f"[{idx}] Skipping existing: {output_filename}")
        return

    print(f"[{idx}] Downloading {output_filename}")

    response = requests.get(url, stream=True, timeout=300)

    if response.status_code != 200:
        print(f"[{idx}] ERROR: {url}")
        print(f"[{idx}] HTTP status = {response.status_code}")
        return

    with open(output_filename, "wb") as f:
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            if chunk:
                f.write(chunk)

    print(f"[{idx}] Done: {output_filename}")


def DownloadHRRR2DData(inputs_file,
                     forecast_time=72,
                     interval=1):
    """
    Download HRRR analysis files for a sequence of hours.

    Parameters
    ----------
    inputs_file : str
        User input file.

    forecast_time : int
        Number of forecast hours to download.

    interval : int
        Hour spacing between files.

    Returns
    -------
    filenames : list
        Files downloaded by this MPI rank.

    area : list
        [north, west, south, east]
    """

    user_inputs = read_user_input(inputs_file)

    area = user_inputs["area"]

    start_time = datetime(
        user_inputs["year"],
        user_inputs["month"],
        user_inputs["day"],
        int(user_inputs["time"].split(":")[0])
    )

    timestamps = generate_timestamps(
        start_time,
        forecast_time,
        interval
    )

    # MPI setup
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    filenames = []

    max_download_ranks = 4
    active_ranks = min(size, max_download_ranks)

    for idx, dt in enumerate(timestamps):

        if rank >= active_ranks:
            continue

        if (idx % active_ranks) != rank:
            continue

        ymdh = dt.strftime("%Y%m%d_%H")

        fname = f"HRRRFiles/2D/hrrr_{ymdh}.grib2"

        url = build_hrrr_url(dt)

        filenames.append(fname)

        download_file(
            url=url,
            output_filename=fname,
            idx=idx
        )

    comm.Barrier()

    return filenames, area
