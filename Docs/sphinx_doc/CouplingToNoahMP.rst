   .. role:: cpp(code)
      :language: c++

Coupling to Noah-MP
===================

Overview
--------

The NOAH Land Surface Model (LSM) is integrated with ERF to facilitate
interaction with the Noah-MP (Multi-Physics) land surface processes.

This documentation covers key components of this interface and its
associated data structures and routines, which are implemented in C++
and Fortran, and provides details on initialization and management of
data structures necessary for these processes.

Files Overview
--------------

-  **Source/LandSurfaceModel/NOAH/ERF_NOAH.H**: Contains the declaration
   of the NOAH class, which extends the NullSurf.

-  **Source/LandSurfaceModel/NOAH/ERF_NOAH.cpp**: Implements the
   initialization routine for the NOAH class.

-  **Submodules/NOAH-MP/drivers/hrldas/NoahmpIO.H**: Defines the C++
   `NoahmpIO_type` that is used to interface with Noah-MP implementations
   following similar structure as the underlying Fortran interface
   (https://dx.doi.org/10.5065/ew8g-yr95).

-  **Submodules/NOAH-MP/drivers/hrldas/NoahmpIO.cpp**: Contains the
   implementation of C++ routines interfacing with Fortran.

-  **Submodules/NOAH-MP/drivers/hrldas/NoahmpIO_fi.F90**: Fortran module
   responsible for managing mapping data between C++ and Fortran.

NOAH Class
----------

The NOAH class serves as the handler for initializing and managing the
data structures required for NOAH-MP operations. It inherits from the
`NullSurf` class. This class declares private variable `NoahmpIO_type
noahmpio`, that is passed to NoahMP routines similar to the Fortran
interface in the Noah-MP documentation
(https://dx.doi.org/10.5065/ew8g-yr95)

NoahmpIO_type Structure
-----------------------

This structure is key for handling the input and output operations for
NOAH-MP through C++ and Fortran interoperation. Contains various
variables for domain, memory, and tile configuration. Also, contains
arrays for geographic variables. At present this type exposes only a
select set of variables. More variables should be exposed as needed by
applications in ERF. The process of adding new variables is as follows:

#. In **Submodules/NOAH-MP/drivers/hrldas/NoahmpIO.H** add pointers to
   the desired variable and set their initialization for
   `NoahmpIO_type_fi` similar to implementation of `WSLAKEXY` and
   `XLAT`.

#. In **Submodules/NOAH-MP/drivers/hrldas/NoahmpIO.H** declare objects
   for Fortran-style multidimensional arrays for the same variables in
   `NoahmpIO_type` similar to implemnation of `NoahArray2D<double> XLAT`
   and `NoahArray2D<double> WSLAKEXY`.

#. In **Submodules/NOAH-MP/drivers/hrldas/NoahmpIO.cpp** cast the
   pointers from `NoahmpIO_type_fi` to multidimensional arrays in
   `NoahmpIO_type` within the implementation of `void
   NoahmpIOVarInitDefault(NoahmpIO_type* noahmpio)`.

Fortran Interoperability
------------------------

The connection between C++ and Fortran is managed through `NoahmpIO_fi`.
This module contains a mirroring of the C++ structure for NOAH-MP
input-output operations.

The following functions are used to operate on the `NoahmpIO_type` and
interface with their respective Fortran implementations:

-  `void NoahmpIOVarInitDefault(NoahmpIO_type* noahmpio)`: Initializes
   default variables of `NoahmpIO_type`. Create C pointer for Fortran
   data.

-  `void NoahmpInitMain(NoahmpIO_type* noahmpio)`: Main initialization
   function for the NOAH-MP operations in C++.

Usage
-----

To use the NOAH class and associated functions, ensure the correct
initialization sequence is followed within the simulation setup. The
interplay between C++ and Fortran necessitates careful memory and data
handling, which is crucial for ensuring performance and correctness in
simulations. The interface is designed to mimic the Fortran interface
from documentation(https://dx.doi.org/10.5065/ew8g-yr95), therefore
similar practices should be followed.
