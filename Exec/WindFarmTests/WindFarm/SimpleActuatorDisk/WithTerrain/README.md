# How to run this example

1. Download the USGS terrain file for the region as a GeoTIFF file from the USGS Earth Explorer  
2. Read the USGS terrain file and write out an ERF-readable terrain file 
   ```
	python3 ReadTerrainUSGS.py
	```
Note that python module `rasterio` is needed. Usually, it can be installed as  
	```
		pip install rasterio
	```

On supercomputer clusters, the following might be needed  
```
module load python
conda create -n raster_env -c conda-forge python=3.11 rasterio
conda activate raster_env
```
