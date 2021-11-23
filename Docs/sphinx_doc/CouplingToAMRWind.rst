
 .. role:: cpp(code)
    :language: c++

 .. role:: fortran(code)
    :language: fortran

 .. _CouplingToAMRWind:

Coupling To AMRWind
=====================

.. note::

    Everything below is a work in progress

The simplest form of coupling is one-way, file-based coupling. With this approach, an
ERF simulation is run and outfile files are stored at regular time intervals. Once the
ERF simulation is complete, an AMR-Wind simulation begins, using the ERF output files
as time-dependent boundary conditions.

File-based coupling
-------------------

+----------------------------+------------------+------------------+-------------+
| Parameter                  | Definition       | Acceptable       | Default     |
|                            |                  | Values           |             |
+============================+==================+==================+=============+
| **erf.bdry_coupling_type** | type of output   | "File_1D" or     | "None"      |
|                            | for coupling     | "None"           |             |
|                            | to AMR-Wind      |                  |             |
+----------------------------+------------------+------------------+-------------+
| **erf.output_file**        | prefix for       | String           | “prof”      |
|                            | output files     |                  |             |
+----------------------------+------------------+------------------+-------------+
| **erf.output_int**         | how often (by    | Integer          | -1          |
|                            | level-0 time     | :math:`> 0`      |             |
|                            | steps) to output |                  |             |
|                            | plot files       |                  |             |
+----------------------------+------------------+------------------+-------------+
| **erf.output_per**         | how often (by    | Real :math:`> 0` | -1.0        |
|                            | simulation time) |                  |             |
|                            | to write output  |                  |             |
|                            | files            |                  |             |
+----------------------------+------------------+------------------+-------------+

*  You should specify either **erf.output_int** or **erf.output_per**, but not both.

*  The name of the output files will begin with "prof" unless an alternate name is specified in the inputs file.

1D File output
~~~~~~~~~~~~~~

To create the 1D profiles to be output, the ERF data is horizontally averaged to create
profiles that are a function of the vertical coordinate only.  Because ERF uses height-based
coordinates, the averaging operation is straightforward.

We note that the solution variables in ERF are

.. math::

    (\rho, \rho u, \rho v, \rho w, \rho \theta)

while those in AMR-Wind are

.. math::

    (u, v, T, p)

We convert the variables from conservative to primitive form before/after the horizontal averaging.

.. note::

    What do we do about the fact that AMR-Wind currently uses Boussinesq model which assumes constant density.

The format of the 1D output file is described in the AMR-Wind documentation here <Add link>


