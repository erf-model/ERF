   .. role:: cpp(code)
      :language: c++

Coupling to Noah-MP
===================

Overview
--------

The Noah-MP Land Surface Model (LSM) is integrated with ERF to facilitate
interaction with the Noah-MP (Multi-Physics) land surface processes.

This documentation covers key components of this interface and its
associated data structures and routines, which are implemented in C++
and Fortran, and provides details on initialization and management of
data structures necessary for these processes.

Building and Running with Noah-MP
---------------------------------
To build ERF with Noah-MP support, the ``NetCDF``, ``NetCDF Fortran``, and ``HDF5`` libraries are
required. Furthermore, ``ERF_ENABLE_NOAHMP=ON`` must be specified with CMake builds or ``USE_NOAHMP=TRUE``
and ``USE_NETCDF=TRUE`` must be specified with GNU Make. Once an executable has been generated, the
inputs file for the simulation must specify Noah-MP as the land surface model type:

.. code-block:: bash

    erf.land_surface_model = "NOAHMP"

Currently, Noah-MP may only be utilized for simulations that are initialized from a WRF input file
(``erf.init_type = "WRFInput"``). Additionally, two files are required to be in the run directory
for Noah-MP initialization: ``namelist.erf`` and ``NoahmpTable.TBL``. Sample files are provided for
the :download:`namelist.erf <namelist.erf>` and :download:`NoahmpTable.TBL <NoahmpTable.TBL>`.

Files Overview
--------------

-  **Source/LandSurfaceModel/Noah-MP/ERF_NOAHMP.H**: Contains the declaration
   of the NOAH class, which extends the NullSurf.

-  **Source/LandSurfaceModel/Noah-MP/ERF_NOAHMP.cpp**: Implements the
   initialization routine for the NOAH class.

-  **Submodules/Noah-MP/drivers/erf/NoahmpIO.H**: Defines the C++
   `NoahmpIO_type` that is used to interface with Noah-MP implementations
   following similar structure as the underlying Fortran interface
   (https://dx.doi.org/10.5065/ew8g-yr95).

-  **Submodules/Noah-MP/drivers/erf/NoahmpIO.cpp**: Contains the
   implementation of C++ routines interfacing with Fortran.

-  **Submodules/Noah-MP/drivers/erf/NoahmpIO_fi.F90**: Fortran module
   responsible for managing mapping data between C++ and Fortran.

NOAHMP Class
------------

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

#. In **Submodules/Noah-MP/drivers/erf/NoahmpIO.H** add pointers to
   the desired variable and set their initialization for
   `NoahmpIO_type_fi` similar to implementation of `WSLAKEXY` and
   `XLAT`.

#. In **Submodules/Noah-MP/drivers/erf/NoahmpIO.H** declare objects
   for Fortran-style multidimensional arrays for the same variables in
   `NoahmpIO_type` similar to implemnation of `NoahArray2D<double> XLAT`
   and `NoahArray2D<double> WSLAKEXY`.

#. In **Submodules/Noah-MP/drivers/erf/NoahmpIO.cpp** cast the
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

Generating Fortran–C++ Bindings using CodeScribe
================================================

**CodeScribe** (https://github.com/akashdhruv/CodeScribe) can be used to
automatically generate Fortran–C interoperability bindings and C++ interface code
for the Noah-MP land-surface model within ERF.

The following files can be generated or updated using CodeScribe:

-  **Source/LandSurfaceModel/Noah-MP/ERF_NOAHMP.cpp**
-  **Submodules/Noah-MP/drivers/erf/NoahmpIO.H**
-  **Submodules/Noah-MP/drivers/erf/NoahmpIO.cpp**
-  **Submodules/Noah-MP/drivers/erf/NoahmpIO_fi.F90**


Follow the instructions in the `CodeScribe repository <https://github.com/akashdhruv/CodeScribe>`_
to configure your LLM environment. You will need API access for your preferred
model (e.g., OpenAI, Argo, etc.). Tutorials are available at
`https://github.com/akashdhruv/codescribe-tutorial <https://github.com/akashdhruv/codescribe-tutorial>`_.

1. Edit the prompt file **prompts/noahmpio_update.toml** to specify which
   variables should be exposed to the C++ interface.

2. Run the following commands to generate or update bindings in **Submodules/Noah-MP/drivers/erf** directory:

.. code-block:: bash

   code-scribe update NoahmpIO.H NoahmpIO.cpp NoahmpIO_fi.F90 \
       -p prompts/noahmpio_update.toml \
       -q "Write a natural language prompt with variable names, dimensions, etc." \
       -m <openai|argo-gpt4o|...>

3. Run the following to generate or update bindings in **Source/LandSurfaceModel/Noah-MP** directory:

.. code-block:: bash

   code-scribe update ERF_NOAHMP.cpp \
       -p prompts/noahmpio_update.toml \
       -q "Write a natural language prompt with variable names, dimensions, etc." \
       -m <openai|argo-gpt4o|...>

You may need to manually edit **Submodules/Noah-MP/drivers/erf/NoahmpIOVarType.F90** to replace:

.. code-block:: fortran

   real(kind=kind_noahmp)

with:

.. code-block:: fortran

   real(kind=C_DOUBLE)

This ensures compatibility with the C++ side. Alternatively, CodeScribe can perform this update automatically (depending on your
model’s context length) using:

.. code-block:: bash

   code-scribe update NoahmpIOVarType.F90 \
       -p prompts/noahmpio_update.toml \
       -q "Write a natural language prompt with variable names, dimensions, etc." \
       -m <openai|argo-gpt4o|...>

If you want to control Noah-MP plot variables, you can update **Submodules/Noah-MP/drivers/erf/NoahmpWriteLandMod.F90** file:

.. code-block:: bash

   code-scribe update NoahmpWriteLandMod.F90 \
       -p prompts/noahmpwriteland_update.toml \
       -q "Write a natural language prompt with variable names, dimensions, etc." \
       -m <openai|argo-gpt4o|...>
