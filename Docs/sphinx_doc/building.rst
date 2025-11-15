.. _sec:building:

Building ERF
============

ERF can be built using either GNUMake or CMake. Here are the three most common build commands:

.. code:: shell

   # GNUMake with GNU compiler and MPI
   make COMP=gnu USE_MPI=TRUE

   # GNUMake with GPU support
   make COMP=gnu USE_MPI=TRUE USE_CUDA=TRUE

   # CMake configuration script
   cd Build && ./cmake.sh

Cleanup Commands
----------------

Remove build artifacts:

.. code:: shell

   # GNUMake
   make clean

   # CMake
   cd Build && make distclean

Detailed Build Documentation
-----------------------------

For comprehensive build system documentation, compiler options, HPC-specific configurations, and troubleshooting, see the :ref:`building:configuration` section, which includes:

* :ref:`sec:build:overview` - Comparison of GNUMake and CMake build systems, key variables, and workflows
* :ref:`sec:build:systems` - Detailed technical reference for all build system options and configurations  
* :ref:`sec:build:hpc` - Architectural overview of building on HPC platforms with automated Cray detection
* :ref:`sec:build:library` - Guide to configuring NetCDF, HDF5, and E3SM physics packages with parallel I/O
* :ref:`sec:build:troubleshooting` - Solutions for common build failures and debugging strategies
