 .. role:: cpp(code)
    :language: c++

.. _GettingStarted:

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
    not allow you to run with an incompressible time-step size.
  - The acoustic CFL should be less than 0.5, with 4--6 fast timesteps
    (substeps) according to WRF best practices.

    .. code-block:: python

       erf.fixed_dt           = 0.06  # slow timestep

       # These are equivalent and result in a fixed fast timestep size
       #   if dx=10, speed of sound ~ 350 m/s
       erf.fixed_fast_dt      = 0.01  # ==> CFL~0.35
       #erf.fixed_mri_dt_ratio = 6

       # Alternatively, let ERF chose the fast timestep
       #erf.cfl                = 0.5

  - We note that ERF LESs with up to 10 fast timesteps have successfully been
    run but your mileage may vary.


Single-Column Model
-------------------

* Currently, ERF does not have the ability to run a true single-column model
  (SCM). The grid size in the lateral directions must have a minimum number of
  cells. This will give comparable results, e.g.:

  .. code-block:: python

     geometry.prob_extent = 400  400  400
     amr.n_cell           =   2    2   64
     geometry.is_periodic =   1    1    0

  When set up this way, the solution is not sensitive to horizontal problem
  extent.

* An SCM was successfully run with third-order advection in the horizontal and
  vertical.
