
# Hurricane simulations 

This folder contains examples for hurricane simulations from real weather data.

1. Follow the steps in the erftools directory to generate the initial condition and boundary 
   condition files.  
   For ERA5 data see the README section [here](https://github.com/erf-model/erftools/tree/main/notebooks/era5).    
   For GFS data see the README section [here](https://github.com/erf-model/erftools/tree/main/notebooks/gfs).  

2. From the step above, copy the `Output/HindCastBoundaryDataDir` directory to the directory where ERF will be run.

3. The python script also outputs the domain size to be used in the `inputs` file for the ERF run.
   ```
    geometry.prob_lo  =  -2593434.0 -2065213.0 0.0
    geometry.prob_hi  =  2593434.0 2328015.0 25000.0
    ```


4. Copy the initial condition file to the ERF run directory. This can be the first file in the `HindCastBoundaryDataDir`.  
   It can also be a file that was generated at the same time as the first file in `HindCastBoundaryDataDir`. 
   for eg. The initial condition can be from ERA5 data, but the boundary data can be from GFS data. But the initial condition 
   and the very first file in `HindCastBoundaryDataDir` should correspond to the same date and geographical area.

5. Run ERF.
  
   
