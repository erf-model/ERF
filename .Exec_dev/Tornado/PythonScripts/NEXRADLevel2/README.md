# NEXRAD Level II Radar Data Visualization

This repository contains Python routines to read and plot **NEXRAD Level II radar data**.  
The script reads `.gz` radar files, extracts radial velocity data, and visualizes it with a specified city marker on a planar map.

---

## Overview

The script:

- Reads **NEXRAD Level II** data using [Py-ART](https://arm-doe.github.io/pyart/).
- Plots **radial velocity (m/s)** for a specified sweep.
- Marks a given **city location** on the radar map.
- Saves all plots to an `Images/` subdirectory.
- Optionally, you can assemble all images into a movie using `ffmpeg`.

---

## Usage

Run the script as:

```bash
python3 ReadRadarData.py --data_folder=<folder_with_gz_files> --lat=<latitude_of_city> --lon=<longitude_of_city>
```

## Example

An example is below. The data downloaded is for the Joplin tornado on 22 May 2011.  
The data can be accessed from the Google Cloud Public Dataset:

[Google Cloud NEXRAD Level II Data](https://console.cloud.google.com/storage/browser/gcp-public-data-nexrad-l2)

Follow these steps:

```bash
1. Create a folder for the event
mkdir Joplin
cd Joplin

2. Download the NEXRAD Level II tar file
wget https://storage.googleapis.com/gcp-public-data-nexrad-l2/2013/05/20/KTLX/NWS_NEXRAD_NXL2DP_KTLX_20130520200000_20130520205959.tar

3. Extract the tar file
tar xvf NWS_NEXRAD_NXL2DP_KTLX_20130520200000_20130520205959.tar

4. Run the Python plotting script
python3 ReadRadarData.py --data_folder=./Joplin --lat=37.063 --lon=-94.513
```
The images are written into a directory `Joplin/Images`.
