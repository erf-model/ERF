
.. _sec:HindCast:

Hindcasting using weather data
===============================

Introduction
-------------

ERF supports hindcasting simulations using weather data from ERA5 and GFS. There are automated Python tools
that download and process the weather data, and writes out visualizable output as well as the
initial and boundary condition files for running ERF simulations in a seamless manner.  Currently there are
examples for performing hurricane hindcasting using weather data.

Hurricane simulations
----------------------

This folder contains examples for hurricane simulations from real weather data.

1. Follow the steps in the ``erftools`` directory to generate the initial condition and boundary
   condition files.

   - For ERA5 data see the README section `here <https://github.com/erf-model/erftools/tree/main/notebooks/era5>`_.
   - For GFS data see the README section `here <https://github.com/erf-model/erftools/tree/main/notebooks/gfs>`_.

2. From the step above, copy the ``Output/ERA5Data_3D`` (``Output/GFSData_3D`` for GFS) directory to the directory
   where ERF will be run. This is the boundary data for the lateral forcing of large-scale meteorology. Set the
   inputs option ``erf.hindcast_boundary_data_dir = ERA5Data_3D``.

3. The Python script also outputs the domain size to be used in the ``inputs`` file for the ERF run.
   For example::

       geometry.prob_lo  =  -2593434.0 -2065213.0 0.0
       geometry.prob_hi  =   2593434.0  2328015.0 25000.0

4. Copy the initial condition file to the ERF run directory. This can be the first file in the
   ``hindcast_boundary_data_dir``. It can also be a file that was generated at the same time as the first file in
   ``hindcast_boundary_data_dir``. For example, the initial condition can be from ERA5 data, but the
   boundary data can be from GFS data. However, the initial condition and the very first file in
   ``hindcast_boundary_data_dir`` should correspond to the same date and geographical area.

5. Run ERF.

Inputs for hindcast simulations
-------------------------------

The following are the inputs required for hindcast simulations.

.. code-block:: cpp

    // If using era5 or gfs data to initialize
    // and for boundary forcing
    erf.init_type = "hindcast"
    // Initial condition filename -
    // obtained from running the python script
    // with the inputs file specifying the date and
    // geographical area
    erf.hindcast_IC_filename = "ERF_IC_2025_08_18_00_00_000.bin"

    // Boundary conditions
    geometry.is_periodic = 0 0 0
    xlo.type = "Outflow"
    xhi.type = "Outflow"
    ylo.type = "Outflow"
    yhi.type = "Outflow"
    zlo.type = "SlipWall"
    zhi.type = "SlipWall"

    // Lateral forcing with reanalysis/forecast data
    erf.hindcast_boundary_data_dir = "ERA5Data_3D"
    // Time interval in hours between the boundary
    // data files
    erf.hindcast_data_interval_in_hrs = 3.0
    erf.hindcast_lateral_forcing = true
    // Sponge strength
    erf.hindcast_lateral_sponge_strength = 0.3
    // Sponge length of 144 km (for eg.)
    erf.hindcast_lateral_sponge_length = 144000

    // Sponge damping for top of domain
    // to absorb reflections
    erf.hindcast_zhi_sponge_damping = true
    // Sponge strnegth
    erf.hindcast_zhi_sponge_strength = 0.3
    // Sponge length of 5 km
    erf.hindcast_zhi_sponge_length = 5000.0

    // Coriolis force
    erf.use_coriolis = true
    erf.coriolis_3d = true
