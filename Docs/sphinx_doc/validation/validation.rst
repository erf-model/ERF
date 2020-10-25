
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

 
.. _Validation:


.. highlight:: rst

Validation of ERF
-------------------


The ERF validation plan is aimed at exercising and validating the ERF reacting flow capabilities. The following cases, described further on, are used for validation.

* Decay of homogeneous isotropic turbulence
* Non reacting Taylor-Green vortex breakdown
* Reacting Taylor-Green vortex breakdown

.. warning::

   This section is a work in progress. Several of these cases have yet
   to be performed and are noted as such.


Taylor-Green vortex breakdown
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This setup is one of the test problems outlined by `the High-Order CFD
workshop <https://www.grc.nasa.gov/hiocfd>`_. A complete description
of the problem can be found `at NASA HOCFDW website
<https://www.grc.nasa.gov/hiocfd/wp-content/uploads/sites/22/case_c3.3.pdf>`_
and the reference data is found `here
<https://www.grc.nasa.gov/wp-content/uploads/sites/22/C3.3_datafiles.zip>`_. More
details of the problem and methods used to obtain the reference data
can be found in `Bull and Jameson (2014) 7th AIAA Theoretical Fluid
Mechanics Conference (doi: 10.2514/6.2014-3210)` and `DeBonis (2013)
51st AIAA Aerospace Sciences Meeting (doi:10.2514/6.2013-382)`.

Building and running
####################

The Taylor-Green vortex breakdown case can be found in ``Exec/RegTests/TG``:

.. code-block:: bash

   $ make -j 16 DIM=3 USE_MPI=TRUE
   $ mpiexec -n 36 $EXECUTABLE inputs_3d amr.ncell=64 64 64

The user can run a convergence study by varying ``amr.ncell``.


Results
#######

As the resolution increases, there is good agreement between the Pele
data and reference data.

.. figure:: ./validation/tg/dissipation.png
   :align: center
   :figwidth: 40%

   Dissipation as a function of time. Solid red: :math:`32^3`, dashed green :math:`64^3`, dash-dotted blue: :math:`128^3`, dotted orange: :math:`256^3`, black squares: HOCFDW.

.. figure:: ./validation/tg/enstrophy.png
   :align: center
   :figwidth: 40%

   Enstrophy as a function of time. Solid red: :math:`32^3`, dashed green :math:`64^3`, dash-dotted blue: :math:`128^3`, dotted orange: :math:`256^3`, black squares: HOCFDW.

.. figure:: ./validation/tg/KE.png
   :align: center
   :figwidth: 40%

   Kinetic energy as a function of time. Solid red: :math:`32^3`, dashed green :math:`64^3`, dash-dotted blue: :math:`128^3`, dotted orange: :math:`256^3`, black: HOCFDW.

.. figure:: ./validation/tg/spectrum.png
   :align: center
   :figwidth: 40%

   Spectrum at :math:`t=9 t_c`. Solid red: :math:`32^3`, dashed green :math:`64^3`, dash-dotted blue: :math:`128^3`, dotted orange: :math:`256^3`, black: HOCFDW.

Reacting Taylor-Green vortex breakdown
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This test case is based on work by `Abdelsamie et al. (Mini-Symposium
on Verification and Validation of Combustion DNS, 17th Int. Conference
on Numerical Combustion, Aachen, Germany May 7, 2019` where a
Taylor-Green vortex setup used in non-reacting CFD is adapted to a
reacting flow configuration. Comparison of results from several
well-established codes such as Nek5000, DINO and YALES are provided in
the workshop documentation. We have performed the entire suite of
cases described in the workshop documentation and only present the
final 3D reacting case.

Good comparisons with the reference simulations were obtained in most
of the quantities of interest.

.. figure:: ./validation/rtg/ux.png
   :align: center
   :figwidth: 40%

   :math:`x`-velocity at :math:`t=5e-4 \tau`. Solid red: :math:`32^3`, dashed green :math:`64^3`, dash-dotted blue: :math:`128^3`, black: reference solution (DINO).

.. figure:: ./validation/rtg/yh2.png
   :align: center
   :figwidth: 40%

   :math:`Y_{H_2}` at :math:`t=5e-4 \tau`. Solid red: :math:`32^3`, dashed green :math:`64^3`, dash-dotted blue: :math:`128^3`, black: reference solution (DINO).

.. figure:: ./validation/rtg/hr.png
   :align: center
   :figwidth: 40%

   Heat release at :math:`t=5e-4 \tau`. Solid red: :math:`32^3`, dashed green :math:`64^3`, dash-dotted blue: :math:`128^3`, black: reference solution (DINO).

.. figure:: ./validation/rtg/tmax.png
   :align: center
   :figwidth: 40%

   Maximum temperature in the domain as a function of time. Solid red: :math:`32^3`, dashed green :math:`64^3`, dash-dotted blue: :math:`128^3`, black: reference solution (DINO).


.. note::

   We are not using the constant Lewis approximation that is
   prescribed in the workshop documentation. Instead we rely on
   transport coefficients resulting from PelePhysics. This may lead to
   discrepancies with the reference results.

.. note::

   Because of computational constraints, we have not been able to
   perform higher resolution simulations that may show better
   convergence.
