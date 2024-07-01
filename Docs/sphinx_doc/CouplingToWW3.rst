
 .. role:: cpp(code)
    :language: c++

 .. _CouplingToWW3:

Coupling To WW3
===============

Coupling with WaveWatch III is currently a work in progress.
Currently, we have a one-way coupling between ERF and WaveWatch III (WW3), where WW3 sends ERF Hwave (significant wave height) and Lwave (mean wavelength) over a grid. The

One-way coupling WW3 to ERF
---------------------------

Values are used to compute the surface roughness z0 through a fixed-point iteration:

.. math::

    z0 = 1200.0 Hwave \left(\frac{Hwave}{Lwave}\right)^{4.5} + \frac{0.11 \mu}{u_*}

To run the coupled model:

.. code-block:: bash

    git clone --recursive git@github.com:erf-model/ERF
    cd ERF/Exec/ABL
    make -j4 USE_WW3_COUPLING=TRUE
    cd ../../Submodules/WW3
    ./model/bin/w3_setup model -c gnu -s Ifremer1
    cd regtests
    ./bin/run_cmake_test -C MPMD -n 2 -p mpirun -f -s PR1_MPI ../model ww3_tp2.2


