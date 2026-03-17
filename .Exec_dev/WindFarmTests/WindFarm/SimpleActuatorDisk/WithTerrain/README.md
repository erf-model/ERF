## Wind Farm - Simple actuator disk with terrain

This folder contains an example of the Simple actuator disk wind farm parametetrization 
for a wind farm with multiple wind turbines with terrain. The steps for running this example are

1. Download the USGS terrain file for the region as a GeoTIFF (`.tif`) file from the USGS Earth Explorer.
   The file for this example can be obtained by the following command  
```
wget https://zenodo.org/record/14629890/files/n36_w098_1arc_v3.tif
```

2. Read the USGS terrain file and write out an ERF-readable terrain file
```
python3 ReadTerrainUSGS.py n36_w098_1arc_v3.tif domain_bounds.txt
```
the `domain_bounds.txt` contains the extents of the terrain domain as 
```
min longitude
max longitude 
min latitude
max latitude
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
