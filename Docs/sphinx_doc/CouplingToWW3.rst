
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

Modifications to the problem size and geometry, as well as other parameters can be done in the `inputs_mpmd` file. The `plt` files as well as relevant outputs can be viewed in the `regtests/ww3_tp2.2/work` directory.

Two-way coupling:
-----------------

Disclaimer: Two-way coupling is currently a work in progress. Two-way coupling involves sending the wind velocity and direction to WW3. We convert the x and y velocities from ERF to a wind speed (U_{wind}) and wind direction (\theta).

.. math::

    \text{Wind speed at reference height} \quad U_{wind} = \sqrt{u^2 + v^2}

.. math::

    \text{Wind direction at reference height} \quad \theta = \arctan{\frac{v}{u}}

Both U_{wind} and \theta are then used in the wind source term S_{in} in the ST6 subroutine in WW3. To run the model with two-way coupling:

.. code-block:: bash

    ./bin/run_cmake_test -C MPMD -n 2 -p mpirun -f -s PR1_ST6 ../model ww3_ts2

