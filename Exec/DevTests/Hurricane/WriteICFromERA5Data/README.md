# Read weather data from ERA5 and write initial condition file for ERF

This directory contains the python scripts to create an initial condition for hurricane simulations in ERF from the ERA5 weather data.

1. Follow the steps here https://cds.climate.copernicus.eu/how-to-api to create a free account   
   to download ERA5 data

2. Give the year, month, day, time and the geographical area to download the data in a text file.  
For example, to download the data on August 26, 2020 at 00:00 (24 hour format)
```
year: 2020
month: 08
day: 26
time: 00:00
area: 50,-130,10,-50
```
Note: The geographical area is specified as latitude maximum, longitude minimum, latitude minimum, longitude maximum.

3. `python3 ReadERA5DataAndWriteERF_IC.py <input_file>`   
If any packages are missing, install them using `pip install <package>`

4. The output VTK files for visualization is written into a directory `Output`. The initial condition binary file (`*bin`) for ERF   
is also written into `Output`.

## Examples

Example inputs are given in the input file `input_for_era5_Laura` and `input_for_era5_Henri`. 

1. Run `python3 ReadERA5DataAndWriteERF_IC.py input_for_era5_Laura`.  
2. Visualize the VTK files in the `Output` directory in VisIt or ParaView.

