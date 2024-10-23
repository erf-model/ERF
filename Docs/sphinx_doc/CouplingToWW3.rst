
 .. role:: cpp(code)
    :language: c++

 .. _CouplingToWW3:

Coupling To WW3
===============

Coupling with WaveWatch III is currently a work in progress.
We have established two-way coupling between ERF and WaveWatch III (WW3),
in which WW3 sends ERF :math:`Hwave` (significant wave height) and :math:`Lwave` (mean wavelength) at the lower boundary,
and ERF sends WW3 the wind magnitude and direction:

.. math::

  U_{wind} = \sqrt{u^{2} + v^2}

.. math::

  \theta = \mathrm{arctan}{\frac{v}{u}}


The bidirectional send and receive currently implemented within the WaveWatch III (WW3) regression test ``ww3_ts3``. This regression test uses the wind velocity :math:`U_{wind}` and wind direction :math:`\theta` sent from ERF to compute an atmospheric source term :math:`S_in` in the ST6 subroutine in WW3. The values for :math:`Hwave` and :math:`Lwave` are then sent to ERF to compute the surface roughness parameter :math:`z_0` through a fixed point iteration:

.. math::
  z_{0} = 1200.0 Hwave (\frac{Hwave}{Lwave})^{4.5} + \frac{0.11 \mu}{u_*}

To run the coupled model:

.. code-block:: bash

    git clone --recursive git@github.com:erf-model/ERF
    cd ERF/Exec/DevTest/ABL_with_WW3
    make -j4 USE_WW3_COUPLING=TRUE
    cd ../../Submodules/WW3
    ./model/bin/w3_setup model -c gnu -s Ifremer1
    cd regtests
    ./bin/run_cmake_test -C MPMD -n 2 -p mpirun -f -s ST6_PR1_MPI -w work ../model ww3_ts3

Modifications to the problem size and geometry, as well as other parameters can be done in the `inputs_mpmd` file. The `plt` files as well as relevant outputs can be viewed in the `regtests/ww3_ts3/work` directory.

