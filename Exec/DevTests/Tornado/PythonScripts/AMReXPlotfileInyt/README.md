# AMReX plotfile in python using yt

This directory has the python routines for plotting AMReX data in WRF style  
using yt. The command is 
```
python3 PlotAMReXFile.py <plot_file_name> --var=<variable_to_plot> --area=<longitude_min, longitude_max, latitude_min, latitude_max>`
```
For example,
```
python3 PlotAMReXFile.py plt00100 --var=max_reflectivity --area=-103.5,-80.5,29.8,46.18
```


