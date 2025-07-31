# Initializing ERF from a NetCDF file.

This folder contains a proof-of-concept demonstration for a newer workflow scheme, 
originally developed by Timothy Sliwinski at CIRA/CSU/NOAA GSL:
1. Create initial conditions using python.
2. Export those values to a NetCDF file.
3. Initialize an ERF simluation using that NetCDF data, along with a standard inputs file.

## Background

ERF allows for multiple methods of initializing the state of a simulation;
see [Initialization Pathways](https://erf.readthedocs.io/en/latest/Initialization.html).
For this demonstration, we adapted the [Isentropic Vortex](https://github.com/erf-model/ERF/tree/development/Exec/DryRegTests/IsentropicVortex) 
problem. Here, the file `inputs` corresponds to the inputs file [inputs_advecting](https://github.com/erf-model/ERF/blob/development/Exec/DryRegTests/IsentropicVortex/inputs_advecting). It preserves the problem geometry and input parameters, but is modified to accept initialization from NetCDF (see instructions below). 
The resulting plotfiles for the initial time step, i.e. `plt00000`, for these two problems and the respective inputs are
expected to be the same up to the order of $10^{-8}$. 

Note that ERF accepts input from a NetCDF file to specify the *total state* at the initial time step.
By contrast, the ERF C++ API allows users to calculate custom perturbations that are applied to a base state
by overriding the method `ProblemBase::init_custom_pert()`. Therefore, in adapting a perturbative problem like this one, it is
necessary to calculate `base_state + perturbation` for our initialization data.

## Setup

### Tell ERF to use NetCDF for initialization

In order to tell ERF to use NetCDF data, we minimally need the following in our inputs
file:

```
erf.init_type = NCFile
erf.nc_init_file_0 = "initial_data.nc"
```

Note that `initial_data.nc` is the name of our NetCDF file, specified in our python
script when we export the data. If you use a different file name, adjust your inputs
file accordingly.

### Calculate initialization data in python

Creating our initialization data in python allows us to use the familiar [numpy](https://numpy.org/doc/) array APIs to create our grid(s) and populate them.
As is standard for ERF, we use Arakawa C-cells in this example, where conserved quantities, e.g. $\rho$, $\theta$, are stored at
cell centroids, and velocity values are staggered. Here, our `inputs` file specifies a grid with dimensions 48x48x4. Therefore,
we create numpy arrays for our conserved quantities with those dimensions, e.g.

```python
Rho = np.ndarray(n_cell, np.float64)
```

where `n_cell` is an array with the values `[48, 48, 4]`. Similarly, velocity components
are stored in staggered arrays, e.g.

```python
x_vel = np.ndarray((Nx_face, Ny_cell, Nz_cell), np.float64)
```

where `Nx_face` specifies the number of faces in the x-direction (in this case 49).

Next we calculate the relevant values for our data arrays iteratively, replicating 
the problem specification, but taking into account the fact that we are calculating
`base_state + perturbation`, as mentioned above.

### Export to NetCDF

We use the python package [netCDF4](https://unidata.github.io/netcdf4-python/). 
Please refer to the documentation for details on installation, setup, and use.

We use the method `createDimension` to give a name to the grid axis that we populate with cell values. For example:

```python
bottom_top_dim = outfile.createDimension("BottomTop", Nz_cell)
bottom_top_stag_dim = outfile.createDimension("BottomTopStag", Nz_face)
```

creates a dimension `bottom_top_dim` enumerating the number of cell-centered values in the z-direction, and `bottom_top_stag_dim` gives the number of face-centered values.


We then use `createVariable` to create named variables that will be stored as NetCDF fields. Note that ERF expects the naming convention seen in this example file, i.e.

| Quantity   | Name   |
|------------|--------|
| rho        | "RHO"  |
| rho_theta  | "T"    |
| rho_scalar | "SCAL" |
| x-component of velocity | "U" |
| y-component of velocity | "V" |
| z-component of velocity | "W" |

Lastly, we populate these variable fields in the NetCDF file by copying the data
from the numpy arrays into the variable fields at time step 0. For example,

```python
uwind_var[0, :, :, :] = x_vel
```

sets the x-component of velocity at time coordinate 0 and spatial coordinate `i, j, k`.

**Note**: NetCDF input is expected to use the convention
`t, z, y, x` for space-time coordinates. In this example, we populated our data arrays in `x, y, z` format for clarity, but this requires a call to `numpy.swapaxes()` to rearrange the data into the expected format.

### Create files `ERF_Prob.H` and `ERF_Prob.cpp`

As usual, it is necessary to create your own `class Problem` which inherits from `ProblemBase` and defines a constructor. In this case, though, no additional setup is required in the constructor override, as all of our problem setup has been taken care of in our NetCDF file.
