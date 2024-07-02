
 .. role:: cpp(code)
    :language: c++

 .. _CouplingToWW3:

Coupling To WW3
===============

Coupling with WaveWatch III is currently a work in progress.
Currently, we have a one-way coupling between ERF and WaveWatch III (WW3), where WW3 sends ERF Hwave (significant wave height) and Lwave (mean wavelength) over a grid. 

One-way coupling WW3 to ERF
---------------------------

The values are used to compute the surface roughness :math:`\overline{z_{0}}` through a fixed-point iteration:

.. math::
  \overline{z_{0}} = 1200.0 Hwave (\frac{Hwave}{Lwave})^{4.5} + \frac{0.11 \mu}{u_*}

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

Disclaimer: Two-way coupling is currently a work in progress. Two-way coupling involves sending the wind velocity and direction to WW3. We convert the x and y velocities from ERF to a wind speed (:math:`\overline{U_{wind}}`) and wind direction (:math:`\theta`) from the reference height.

.. math::

  \overline{U_{wind}} = \sqrt{\overline{u^{2}} + \overline{v^2}}

.. math::

  \overline{\theta} = \mathrm{arctan}{\frac{v}{u}}

Both :math:`\overline{U_{wind}}` and :math:`\theta` are then used in the wind source term :math:`S_{in}` in the ST6 subroutine in WW3
