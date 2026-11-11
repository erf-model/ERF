
 .. role:: cpp(code)
    :language: c++

.. _ConstantsAndUnits:

Constants and Units
===================

Units
-----

The following units are used in ERF:

   +-----------------------+-----------------------+-----------------------+
   | Name                  | Units                 | Description           |
   +=======================+=======================+=======================+
   | :math:`t`             | :math:`s`             | time                  |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\rho`          | :math:`kg/m^3`        | density               |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\mathbf{u}`    | :math:`m/s`           | velocity              |
   +-----------------------+-----------------------+-----------------------+
   | :math:`p`             | :math:`Pa`            | pressure              |
   +-----------------------+-----------------------+-----------------------+
   | :math:`T`             | :math:`K`             | temperature           |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\theta`        | :math:`K`             | potential temperature |
   +-----------------------+-----------------------+-----------------------+


Constants
---------

The following are ERF's fixed reference values. They should not be
confused with height-dependent reference-state fields or runtime inputs.

.. list-table:: Selected physical constants
   :header-rows: 1
   :widths: 23 27 50

   * - Symbol
     - Value
     - Meaning
   * - :math:`R_d`
     - 287.0 J/(kg K)
     - Dry-air gas constant.
   * - :math:`R_v`
     - 461.505 J/(kg K)
     - Water-vapor gas constant.
   * - :math:`C_{p,d}`
     - 1004.5 J/(kg K)
     - Fixed reference dry-air heat capacity.
   * - :math:`P_{00}`
     - :math:`10^5` Pa
     - Fixed reference pressure for potential temperature and the EOS
       (the C++ constant is named ``p_0``); not the hydrostatic profile
       :math:`p_0(z)`.
   * - :math:`g`
     - 9.81 m/s\ :sup:`2`
     - Magnitude of gravitational acceleration.
   * - :math:`\Gamma`
     - 1.4
     - Fixed compressible EOS exponent, consistent with
       :math:`C_{p,d}/(C_{p,d}-R_d)` at the reference constants.

The runtime option ``erf.c_p`` defaults to :math:`C_{p,d}` but does not
change the fixed EOS exponent :math:`\Gamma`. If ``erf.c_p`` is changed,
thermodynamic conversions that use :math:`R_d/c_p` can be inconsistent
with inverse relations using the fixed exponent. See
:ref:`GoverningEquations` and :ref:`Buoyancy`.

These constants are defined in the file  :cpp:`Source/ERF_Constants.H`, which holds the
thermodynamic and dynamical constants used throughout ERF.  Two companion headers hold
the rest:

- :cpp:`Source/ERF_NumericalConstants.H` -- dimensionless numeric literals
  (:cpp:`zero`, :cpp:`one`, :cpp:`myhalf`, ...) and the mathematical constant
  :math:`\pi`.  These carry no physical meaning; they exist so the code stays
  precision-agnostic between double and single builds.  :cpp:`ERF_Constants.H`
  includes this header.

- :cpp:`Source/Microphysics/ERF_MicrophysicsConstants.H` -- constants used only by the
  moisture and cloud-physics code: hydrometeor densities, the temperature thresholds that
  partition condensate among the hydrometeor species, terminal fall-speed coefficients,
  autoconversion thresholds and collection efficiencies, size-distribution intercepts,
  and the latent heats of condensation, fusion and sublimation.  This directory is on the
  include path for every ERF build, so any file may include it.
