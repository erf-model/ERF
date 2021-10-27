Energy Research and Forecasting (ERF): An atmospheric modeling code
----

`ERF` is built upon the `AMReX <https://amrex-codes.github.io/amrex/>`_ software framework
for massively parallel block-structured applications.

Test Status
~~~~~~~~~~~

=================  =============  ====================
Regression Tests    |regtests|     |regtest-coverage|
=================  =============  ====================

.. |regtests| image:: https://github.com/erf-model/ERF/actions/workflows/ci.yml/badge.svg?branch=development
.. |regtest-coverage| image:: https://codecov.io/gh/erf-model/ERF/branch/development/graph/badge.svg?token=74S5Q66M1M
    :target: https://codecov.io/gh/erf-model/ERF
.. |unittests| image:: https://github.com/erf-model/ERF/actions/workflows/ci.yml/badge.svg?branch=development


Getting Started
~~~~~~~~~~~~~~~

* To compile and run the `ERF` suite of codes, one needs a C++ compiler that supports the C++14 standard.  A hierarchical strategy for parallelism is supported, based on MPI + OpenMP or MPI + CUDA.  The codes work with all major MPI and OpenMP implementations.  The codes should build and run with no modifications to the `make` system if using a Linux system with the GNU compilers, version 4.9.3 and above.

To clone the source code of `ERF`:

1. One can have ERF use the default submodule for ``AMReX`` in its own repo by simply performing: ::

    export ERF_HOME=<location for ERF code>
    git clone --recursive git@github.com:ERF-model/ERF.git ${ERF_HOME}

    export AMREX_HOME=${ERF_HOME}/Submodules/AMReX # w.r.t. ERF_HOME

Note that the path ``AMREX_HOME`` is dependent on ``ERF_HOME``.

Alternatively, one can set environment variable ``AMREX_HOME`` to use ``AMReX`` repo from external locations independent of ``ERF_HOME``: ::

1. Set the environment variable ``AMREX_HOME`` and clone a copy of ``AMReX`` there: ::

    export AMREX_HOME=<location for AMReX>  # Can be anywhere
    git clone git@github.com:AMReX-Codes/amrex.git ${AMREX_HOME}

2. Set the environment variable ``ERF_HOME`` and clone a copy of ``ERF`` there. You should be placed in the ``development`` branch: ::

    export ERF_HOME=<location for ERF code>
    git clone git@github.com:ERF-model/ERF.git ${ERF_HOME}

Note that cloning using the format ``git@github.com:ERF-model/ERF.git`` requires that your public rsa keys are set in your github profile. See `here <https://docs.github.com/en/free-pro-team@latest/github/authenticating-to-github/adding-a-new-ssh-key-to-your-github-account>`_ how to do this.

If one doesn't have the rsa keys setup in github profile, one can clone using ``https://github.com/erf-model/ERF.git`` format.

To build the code and run a sample problem:

1. Move to an example build folder, build an executable, and run a test case: ::

    cd ${ERF_HOME}/Exec/ScalarAdvection
    make
    ./ERF3d.xxx.yyy.ex inputs_ex

   In order to make on multiple processors use ``make -j N`` instead of ``make``, where ``N`` is the number of processors to use for building and is typically 8 or 16. The executable name is of the format ``ERF3d.xxx.yyy.ex``.

* Notes:

   A. In the exec line above, ``xxx.yyy`` is a tag identifying your compiler and various build options, for example, ``gnu.MPI`` and will vary across pltaform.  (Note that GNU compilers must be at least 4.8.4, and MPI should be at least version 3).
   B. In addition to informative output to the terminal, periodic plot files are written in the run folder.  These may be viewed with CCSE's `Amrvis <https://ccse.lbl.gov/Downloads/downloadAmrvis.html>`_ or `Vis-It <http://vis.lbl.gov/NERSC/Software/visit/>`_ (!! FIX the URLs !!):

      1. In Vis-It, direct the File->Open dialogue to select the file named "Header" that is inside each plotfile folder.
      2. With Amrvis, ``amrvis3d plt00030``, for example.


Development model
~~~~~~~~~~~~~~~~~

See CONTRIBUTING.md for how to contribute to ERF development.

.. note::

   Github CI uses the ``CMake`` build system and ``CTest`` to test the core source files of ERF. If you are adding source files, you will need to add them to the list of source files in the ``CMake`` directory for the tests to pass. Make sure to add them to the GNU make makefiles as well.

Acknowledgment
~~~~~~~~~~~~~~

The development of the Energy Research and Forecasting (ERF) code is funded by the Wind Energy Technologies Office (WETO), part of the U.S. Department of Energy (DOE)'s Office of Energy Efficiency & Renewable Energy (EERE).

The developers of ERF acknowledge and thank the developers of the AMReX-based
`PeleC <https://github.com/AMReX-combustion/PeleC>`_ ,
`FHDeX <https://github.com/AMReX-FHD/FHDeX>`_ and
`AMR-Wind <https://github.com/Exawind/amr-wind>`_ codes.  In the spirit of open source code
development, the ERF project has ported sections of code from each of these projects rather
than writing them from scratch.
ERF is built on the `AMReX <https://github.com/AMReX-codes/AMReX>`_ library.

