# Running Hurricane Simulations with WRF Files

This folder contains the input files along with the links to WRF files (`wrfinput`, `wrfbdy` and `wrflow`) for running different hurricane simulations.

## Running a Hurricane Simulation

1. Download the required WRF files for the hurricane from the corresponding links below.
2. Compile ERF by running:
   ```
   make -j8
   ```
   in the `Exec` directory.
3. Run ERF using the corresponding `inputs` file for hurricane from the `InputFiles` directory.


## WRF files

### Katrina
```
wget "https://zenodo.org/records/21083216/files/WRFFiles_HurricaneKatrina2005.zip?download=1" -O WRFFiles_HurricaneKatrina2005.zip
```

## Example

For **Hurricane Katrina**:

1. Download the WRF files
2. Use one of the `inputs` files in the `InputFiles/Katrina` directory to run ERF.
