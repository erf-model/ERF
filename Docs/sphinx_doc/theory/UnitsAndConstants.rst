
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

   +-----------------------+-----------------------+--------------------------+
   | Name                  | Value                 | Description              |
   +=======================+=======================+==========================+
   | :math:`R_d`           | 287.0  J/kg/K         | gas constant for dry air |
   +-----------------------+-----------------------+--------------------------+
   | :math:`c_p`           | 1004.5 J/kg/K         | specific heat capacity   |
   +-----------------------+-----------------------+--------------------------+
   | :math:`p_0`           | 1.0e5  Pa             | reference pressure       |
   +-----------------------+-----------------------+--------------------------+
   | :math:`g`             | 9.81   m/s/s          | gravity                  |
   +-----------------------+-----------------------+--------------------------+
   | :math:`\gamma`        | 1.4                   | :math:`c_p /(c_p-R_d)`   |
   +-----------------------+-----------------------+--------------------------+

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
