.. _sec:building:

Building ERF
============

ERF supports two build systems: GNUMake and CMake. Most users should start with GNUMake. When on Cray HPC platforms like Perlmutter, CMake provides automated detection to simplify configuration. This page provides quick-start commands; comprehensive documentation follows in the sections below.

Quick Start Commands
--------------------

Here are the three most common build commands:

.. code-block:: bash

   # GNUMake with GNU compiler and MPI
   make COMP=gnu USE_MPI=TRUE

   # GNUMake with GPU support
   make COMP=gnu USE_MPI=TRUE USE_CUDA=TRUE

   # CMake configuration script
   cd Build && ./cmake.sh

Cleanup Commands
----------------

Remove build artifacts:

.. code-block:: bash

   # GNUMake
   make clean

   # CMake
   cd Build && make distclean

Detailed Build Documentation
-----------------------------

For comprehensive build system documentation, compiler options, HPC-specific configurations, and troubleshooting, see the sections below:

* :ref:`sec:build:overview` - Comparison of GNUMake and CMake build systems, key variables, and workflows
* :ref:`sec:build:systems` - Technical reference for all build system options and configurations  
* :ref:`sec:build:hpc` - Building on HPC platforms with automated Cray detection
* :ref:`sec:build:library` - Configuring NetCDF, HDF5, and E3SM physics packages with parallel I/O
* :ref:`sec:build:troubleshooting` - Solutions for common build failures and debugging strategies
