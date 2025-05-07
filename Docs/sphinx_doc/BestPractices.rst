 .. role:: cpp(code)
    :language: c++

.. _BestPractices:

Best Practices
==============

Please note this section is a work in progress and aims to offer guidelines
but, in general, your mileage may vary.

Large-Eddy Simulations
----------------------

* Advection Scheme

  - A WRF-like configuration is generally robust, but may under-resolve
    turbulence features.

    .. code-block:: python

       erf.dycore_horiz_adv_type  = "Upwind_5th"
       erf.dycore_vert_adv_type   = "Upwind_3rd"
       erf.dryscal_horiz_adv_type = "Upwind_5th"
       erf.dryscal_vert_adv_type  = "Upwind_3rd"

  - Centered difference schemes will generally give some non-physical
    numerical noise, clearly visible in the free atmosphere, but may also
    better resolve turbulence features. With ``Centered_2nd``, the simulation
    may remain numerically stable but without any upwinding or numerical
    diffusion, these results should be carefully interpreted.

  - For higher-order central differencing alone (i.e., without any added
    upwinding), at least 5% numerical diffusion should be included to stabilize
    the solution; this was tested with ``Centered_6th``. Note that this does not
    necessarily kill the numerical noise and is only for numerical stability.
    These options are identical to WRF's ``diff_6th_opt`` (default: off) and
    ``diff_6th_factor`` (default: 12%) options.

    .. code-block:: python

       erf.use_NumDiff  = true
       erf.NumDiffCoeff = 0.05

* Time Integration

  - Split timestepping offers some computational cost savings but still does
    not allow you to run with an incompressible/anelastic time-step size in
    general.
  - The acoustic CFL should conservatively be less than or equal to 0.5, with
    4--6 fast timesteps (substeps) according to WRF best practices. If not
    explicitly specified (through ``erf.fixed_mri_dt_ratio`` or
    ``erf.fixed_fast_dt``, the number of substeps in ERF is chosen based on the
    same algorithm as WRF. If the user follows the recommendation that
    dt [s] ~ 6 dx [km],
    then 4 substeps will be used, giving an effective CFL of approximately 0.5.
    This meets the stability criteria from Wicker & Skamarock 2002 that, for a
    5th-order scheme, the CFL be less than 1.42/sqrt(3) = 0.820.

    .. code-block:: python

       erf.fixed_dt           = 0.06  # slow timestep

       # These are equivalent and result in a fixed fast timestep size
       #   if dx=10, speed of sound ~ 300 m/s
       erf.fixed_mri_dt_ratio = 4
       #   or
       #erf.fixed_fast_dt      = 0.015  # ==> CFL~0.45
       #   or, let ERF chose the fast timestep
       #erf.substepping_cfl    = 0.5

  - Following the WRF guidelines for dt is conservative. More aggressive time
    integration--i.e., larger time steps with more substeps--is possible. We
    note that ERF LESs with 10 or more fast timesteps have successfully been
    run but your mileage may vary.


Single-Column Model
-------------------

* An SCM is set up with a single cell in the lateral directions:

  .. code-block:: python

     geometry.prob_extent = 400  400  400
     amr.n_cell           =   1    1   64
     geometry.is_periodic =   1    1    0

* An SCM was successfully run with third-order advection in the horizontal and
  vertical.


2-D Cases
---------

* A 2-D planar domain can be configured as follows:

  .. code-block:: python

     geometry.prob_extent = 10000  100  1000
     amr.n_cell           =   100    1    20
     geometry.is_periodic =     0    0     0

     ylo.type = "SlipWall"
     yhi.type = "SlipWall"
