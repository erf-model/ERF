# Read weather data from ERA5 and write initial condition file for ERF

1. Follow the steps here https://cds.climate.copernicus.eu/how-to-api to create a free account   
   to download ERA5 data

2. Give the year, month, day, time and the geographical area to download the data in a input file as  
```
year: 2020
month: 08
day: 26
time: 00:00
area: 50,-130,10,-50
```
Note: The geogprahical area is specified as latitude maximum, logitude minimum, latitude minimum, longitude maximum.

3. `python3 ReadERA5DataAndWriteERF_IC.py <input_file>`   


An example input is given in the input file `input_for_era5`. Run `python3 ReadERA5DataAndWriteERF_IC.py input_for_era5`. 

