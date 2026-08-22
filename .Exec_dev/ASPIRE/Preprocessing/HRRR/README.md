# HRRR data pre-processor

This directory contains the Python scripts to process [HRRR data](https://home.chpc.utah.edu/~u0553130/Brian_Blaylock/cgi-bin/hrrr_download.cgi) which gives data over the contiguous US (CONUS) with a horizontal resolution of 3 km.

# Steps

1. Give the year, month, day, start time and the geographical area to download the data in a text file. For example, to download the data on June 18, 2025 at 00:00 (24 hour format)
```
year: 2025
month: 06
day: 18
time: 00:00
area: 40,-102,33,-93
```
Note: The geographical area is specified as latitude maximum, longitude minimum, latitude minimum, longitude maximum.

2. Run the script.    

```
srun -n 32 python3 main.py <input_file> --do_forecast=true --forecast_time_hours=48 --interval_hours=1
```

where `input_file` is the file in Step 1 above. This uses 32 MPI ranks to download and process the HRRR data for a total of 48 hours with an interval of 1 hours.   

3. The output VTK files for visualization is written into  `Output/VTK/3D/HRRRDomain` for 3D HRRR data.

An example inputs file is provided in this directory.
