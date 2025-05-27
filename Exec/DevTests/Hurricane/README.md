
# Hurricane simulations 

This folder contains examples for hurricane simulations from real weather data.

1. Generate the initial condition file (`*bin`) for ERF from ERA5 weather data. See the README section of 
[WriteICFromERA5Data](https://github.com/erf-model/ERF/tree/development/Exec/DevTests/Hurricane) for details.

2. `make -j8
3. Run the simulation with the inputs option 
```
erf.IC_file = <ERF_IC_file (*bin)>
```
to initialize the simulation with the ERA5 data.

## Example

1. TBD
