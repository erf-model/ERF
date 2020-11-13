 .. role:: cpp(code)
    :language: c++
 
.. _GettingStarted:

Getting Started
===============

Navigation
----------

The ERF directory structure is as shown below:

* **Source** - C++ source code

* **Util** - third party utilities

  * BLAS
  * LAPACK
  * VODE

* **Docs**   - ERF documentation 

  * sphinx_doc

* **Exec** - regression tests and various capability demonstrations
  
  * :ref:`ScalarAdvection -- scalar advection with constant velocity and pressure

Setting up a problem to run with ERF involves writing an input file and problem specific code in the 
directory where you would like to run the problem.
ERF is built using the AMReX build system which supports out-of-source builds but as configured in 
ERF requires a specific directory structure. 
Within each case directory in Exec, are the source files that specify the setup of that particular case. 
The user has to build each case by compiling source files using a GNUMakefile which also compiles and links together AMReX and ERF sources.
The source files contained in the case directory are treated preferentially and can override ERF/AMReX source files.  
A few key files that need to be supplied for (most) cases are:

**inputs** -- a text file containing parameters that are ready by the ParmParse capability in AMReX. These include things like number of time steps, grid size, output file frequency, which physics to include, etc. 
A list of available data in the ERF group can be found in ERF/Source/param_includes/erf_params.H

**prob.cpp** -- Routines called at:

  * Initialization (`amrex_probinit`) 
  * To set initial values on the grid (`pc_initdata`)
  * Problem teardown (`pc_prob_close`)

**prob.H** -- Something about prob.H

**prob_parm.H** -- Something about prob_parm.H

**GNUMakefile** -- Options to build profiling, debugging, MPI, OpenMP, and Compiler toolchain options are set here for compile time selection. The GNUMakefile includes the ``Make.ERF`` file from the `Exec` directory that contains build configuration common across the examples.


.. include:: building.rst


.. include:: InputFiles.rst


.. include:: tutorials.rst


.. include:: testing.rst
