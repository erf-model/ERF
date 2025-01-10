This folder contains an example of the Simple actuator disk wind farm parametetrization  
for a wind farm with multiple wind turbines with terrain 

The steps for running this example are

1. Download the USGS terrain file for the region as a GeoTIFF (`.tif`) file from the USGS Earth Explorer
2. Read the USGS terrain file and write out an ERF-readable terrain file
```
python3 ReadTerrainUSGS.py
```
This reads in the `.tif` file and writes out `ERF_terrain_file.txt`.
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

3. `make -j8`
4. `mpirun -np 4 <exe> inputs_AWAKEN_SimpleAD_KingPlains_WithTerrain`
