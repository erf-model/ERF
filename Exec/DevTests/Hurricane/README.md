
# Hurricane simulations 

This folder contains examples for hurricane simulations from real weather data. First, a Python pre-processing   
step is required to curate the weather data for ingestion into ERF. This step produces files for the initial   
condition, lateral forcing and surface boundary fluxes for ERF. And then the ERF simulation can be run.

1. Follow the steps in the erftools directory to generate the initial condition and boundary 
   condition files.  
   For ERA5 data see the README section [here](https://github.com/erf-model/erftools/tree/main/notebooks/era5).    
   For GFS data see the README section [here](https://github.com/erf-model/erftools/tree/main/notebooks/gfs).  

2. From the step above, copy the following directories to the directory to the directory where  
   ERF will be run - `Output/ERA5Data_3D` and `Output/ERA5Data_Surface`. These folders have files for the lateral 
   forcing and surface fluxes.  

3. The domain extents to be used in the `inputs` are in `Output/domain_extents.txt`

4. Copy the initial condition file (`hindcast_IC_filename` in the inputs file) to the ERF run directory. This can  
   be the first file in the `ERA5Data_3D`. It can also be a file that was generated at the same time as the first file in `ERA5Data_3D`. 
   For eg. The initial condition can be from ERA5 data, but the boundary data can be from GFS data. But the initial condition 
   and the very first file in `ERA5Data_3D` should correspond to the same date and geographical area.

5. Run ERF.
  
   
