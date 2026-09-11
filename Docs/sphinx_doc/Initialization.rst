 .. role:: cpp(code)
    :language: c++

.. _sec:Initialization:

Initialization Pathways
========================

This section describes different ways for defining
or reading in the initial data for an ERF simulation.
When you run the ERF executable you must specify an inputs file which includes
the specification of which initialization pathway ERF will take.

Custom Initialization
----------------------------------

The deterministic rectangular ``Cloud Chamber`` initializer is documented in
:ref:`CloudChamber`.  Dry mode initializes a physical-temperature profile
without moisture state variables.  SatAdj mode additionally initializes water
vapor from the requested relative humidity and sets cloud water to zero before
the model begins its equilibrium partitioning.  Wall temperature and moisture
settings are a separate boundary-condition contract that controls subsequent
resolved wall-normal transfer.

When not reading the initial data as described in the section below,
the initialization in ERF has two steps: creation of the background state
and creation of optionally non-zero initial perturbations from the background state.

If **erf.init_type = Uniform** the user must provide values in the
inputs file,  **prob.rho_0** and  **prob.T_0**, to
specify the background density and temperature which will
be assumed to be constant in space throughout the domain.
Base state pressure is computed from the EOS.

If **erf.init_type = ConstantDensity** the user must provide
**prob.rho_0** to specify the background density which will
be assumed to be constant in space throughout the domain.
If gravity is set to be non-zero then the density will be vertically
integrated to generate the background pressure
as described in :ref:`sec:BaseState`.  Once density and pressure are
known, the base state potential temperature will be computed from the EOS.

If **erf.init_type = Isentropic**, the background state is computed by
iterating to find base state pressure and density which satisfy both the
dry EOS and hydrostatic equilibrium (HSE).

If **erf.init_type = MoistBaseState**, the background state is computed by
iterating to find base state pressure and density which satisfy both the
moist EOS and hydrostatic equilibrium (HSE).

If **erf.init_type = InputSounding**, then the thermodynamic profiles in the
provided **erf.input_sounding_file** are used to set initial conditions and the
base state depending on **erf.sounding_type**.
For an ``Ideal`` sounding (default), a stratified, hydrostatically balanced base
state is reconstructed from the 1-D input sounding data as described in
:ref:`sec:BaseState`. The stored base-state density is dry-air density
:math:`\rho_0`; base-state pressure is computed from
:math:`\rho_0\theta_0` and :math:`q_{v,0}` using the equation of state, while
hydrostatic balance uses the total moist base-state density
:math:`\rho_0(1 + q_{v,0})`. The initial dry-density,
:math:`\rho\theta`, and water-vapor fields match the base state. This
configuration corresponds to WRF's ideal.exe initialization.

If the sounding is ``Isentropic`` or ``DryIsentropic``, a set of
thermodynamically consistent, isentropic (constant :math:`\theta`) conditions
are determined. An initial pressure profile is calculated by integrating
:math:`p_{surf}` and :math:`\theta_{surf}` from the surface; the resulting
base-state pressure profile is related to :math:`\rho\theta_d` through the
equation of state. This reference profile may or may not take into account the
input water vapor mixing ratio profile, depending the choice of isentropic or
dry isentropic, respectively. Using the input profiles of potential temperature
and water vapor mixing ratio, we determine the initial dry air density. These
three quantities are used to initialize the solution fields. Note that the base
state corresponds to the integrated :math:`p(z)` and
:math:`\theta(z)=\theta_{surf}`.

If the sounding is ``ConstantDensity``, then the initial density field is
uniformly set to 1.0; the potential temperature (and water vapor mixing ratio)
field(s) are set to the sounding values.

For the ``Ideal``, ``Isentropic`` and ``DryIsentropic`` soundings, the profiles are
interpolated onto the cell-centered heights of the mesh and each column is then
rebalanced so that the base state is in discrete hydrostatic equilibrium on that mesh.
With terrain-fitted coordinates this step matters, because the terrain-following heights
differ from the nominal levels the 1-D sounding was integrated on.

This is done on **every level**, against that level's own heights, so a refined level is
hydrostatic on its own mesh rather than inheriting an interpolation of its parent's base
state.  A refined level whose terrain is read from a text file via
``erf.terrain_file_name`` re-reads that file at its own resolution, so it genuinely
resolves topography its parent does not.  A refined region may also cover only part of the
depth of the domain; see :ref:`subsec:base-state-refined` for why that does not change the
base state in the cells it does contain, and :ref:`subsec:partial-depth-refinement` for how
to specify such a region.

The same restrictions apply as for a refined region in any terrain run: the PBL models
(MYJ, MYNN2.5, MYNN-EDMF, YSU, MRF), the SHOC PBL model, and the column-integral derived
quantities (``helicity``, ``precipitable``, ``max_reflectivity``, ``mucape``) all need
entire columns, so they cannot be used on a level whose grids do not span the domain in z.

.. note::

   You can optionally replace only the velocity fields (``u``, ``v``, ``w``)
   by reading them from an existing checkpoint. This is useful when restarting
   with updated thermodynamics or a new base state while keeping a prior wind
   field.

   Add a line such as the following to your inputs file:

   .. code-block:: none

      erf.init_vels_from_checkpoint = chk00010

   The value should be the checkpoint directory name (relative to the run
   directory or an absolute path). When set, ERF reads the velocity fields from
   that checkpoint and uses the usual initialization pathway for all other fields.

In any of these cases, the user can specify any perturbations from the
base state by editing the routines that live in the ``Source/Prob`` directory
and are called in **Exec/ERF_Prob.cpp**

Initialization From Real Data
----------------------------------

There are three options for ingesting the full 3D initial data from a NetCDF file.
In these cases, no additional initial conditions must be supplied by the user but the
file **Exec/ERF_Prob.cpp** must still be present for the build.

* **erf.init_type = WRFInput**

  - In this case both the base state and the full state are read in from a NetCDF
    file generated by the WPS preprocessing system.  This case is designed for
    realistic atmospheric flow problems and allows for terrain-fitted coordinates,
    map factors and a non-zero base state.   In addition to the initial conditions,
    time-varying boundary conditions can be read from a NetCDF file
    specified using ``erf.nc_bdy_file``
    and time-varying sea surface temperatures can be read from a NetCDF file
    specified using ``erf.nc_low_file``.
    We also note that wrfinput files that have later start times than the base level
    can be read and the new level will be initialized once the start time of the finer
    level has been reached.
    A file may be supplied for each level -- ``erf.nc_init_file_0``,
    ``erf.nc_init_file_1``, and so on -- but only level 0 is required to have one;
    see :ref:`sec:nested-wrfinput` below.

* **erf.init_type = Metgrid**

  - In this case the problem is initialized with data
    contained in the first NetCDF file provided via ``erf.nc_init_file_0``.
    Lateral boundary conditions are derived from the sequence of NetCDF
    files provided via ``erf.nc_init_file_0``. The sequence of
    ``erf.nc_init_file_0`` should be output from the WRF Preprocessing
    System (WPS) listed chronologically starting with the earliest
    timestamp. A minimum of two files are required to derive lateral
    boundary conditions.

* **erf.init_type = NCFile**

  - In this case the base state defaults to zero and the full state is read in from
    a much simplified NetCDF file.  Right now, only the density, horizontal and
    vertical velocity components, potential temperature, and water vapor mixing ratio can be read in.
    This case is designed for idealized problems and does not allow
    for terrain-fitted coordinates or map factors.
    The variables in the NC file should have dimensions of (Time, bottom_top, south_north, west_east).
    Variable names include ``RHO``, ``U``, ``V``, ``W``, ``T``, and ``QV``.
    Optional HSE variables include ``RHO_HSE``, ``T_HSE``, and ``P_HSE``; the base state will be
    calculated if it is not specified.

.. _sec:nested-wrfinput:

Nested Initialization From WRF Input Files
------------------------------------------

With ``erf.init_type = WRFInput``, a refined level may be initialized either from its
own wrfinput file or by interpolation from its parent; it may be refined in the
horizontal only or in the vertical as well; and it may cover either the full depth of the
domain or only part of it.  The refinement ratio is given per direction with
``amr.ref_ratio_vect``, or in all three directions at once with ``amr.ref_ratio``::

          amr.ref_ratio_vect = 3 3 1     # refine horizontally only
          amr.ref_ratio_vect = 3 3 3     # refine in the vertical as well
          amr.ref_ratio      = 3         # the same thing, written the short way

.. _subsec:wrfinput-vertical-refinement:

Refining in the vertical
~~~~~~~~~~~~~~~~~~~~~~~~

A WRF nest is refined in the horizontal only: ``wrfinput_d02`` carries exactly the eta
levels of ``wrfinput_d01``, so every file in the set is written on level 0's vertical
grid no matter how finely it resolves the horizontal.  Asking for a refinement ratio
greater than 1 in z therefore means asking ERF to build a vertical grid that the files do
not have, and ERF does so as it reads each file:

-  Every field read from the file -- the geopotential ``PH`` and ``PHB``, the state, and
   the horizontal velocities -- is interpolated from the file's eta levels onto this
   level's, linearly in the vertical index, which is the coordinate the file is
   discretized on.  Each WRF layer is thereby split into :math:`r_z` sublayers of equal
   eta thickness.  Two-dimensional surface fields and the soil fields are left alone;
   the soil layers are not atmospheric levels.

-  Because the geopotential is refined in exactly the same way as the state, the nodal
   heights ERF then builds from ``PH + PHB`` and the state it holds remain collocated.
   The vertical remap onto the ERF mesh that follows therefore has essentially nothing
   left to move, rather than smoothing the profile a second time.

-  The base state is rebuilt on the refined heights by the same construction level 0
   uses -- the analytic WRF reference profile evaluated at this level's own
   cell-centered heights, then rebalanced into discrete hydrostatic equilibrium.  It is
   not interpolated from the coarse level.  See `The base state at each level`_ below.

With the default ``erf.avg_grid_faces_to_nodes = false``, ERF does not use the refined
WRF heights as its mesh; it reconstructs the nodal heights of the first two layers, takes
the first-layer thickness from them, and solves for the stretching factor that fills the
domain with this level's (now :math:`r_z` times as many) layers.  The fine level thus gets
a stretched grid whose first layer is :math:`1/r_z` as thick as the parent's, and the
state is remapped onto it from the refined WRF heights.  With
``erf.avg_grid_faces_to_nodes = true`` the mesh is the refined WRF mesh itself, so every
:math:`r_z`-th fine interface coincides exactly with a WRF interface.

Vertical refinement composes with everything else described here: the nest may still be
given its own file or be interpolated from its parent, and it may still stop below the
domain top.

Where the refined region goes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a file is given for the finer level, ERF reads the nest's position and size from that
file's ``I_PARENT_START``, ``J_PARENT_START`` and ``PARENT_GRID_RATIO`` attributes, and
by default refines exactly the region the file covers, over the full depth of the domain.

To refine only part of that region, add a refinement indicator.  The box it specifies
must be contained in the region covered by the file; ERF checks this and aborts with both
boxes printed if it is not.  The same mechanism is what lets the refined region stop below
the domain top::

          amr.max_level      = 1
          amr.ref_ratio_vect = 3 3 1

          erf.init_type      = WRFInput
          erf.nc_init_file_0 = "wrfinput_d01"
          erf.nc_init_file_1 = "wrfinput_d02"

          erf.refinement_indicators = box1
          erf.box1.max_level = 1
          erf.box1.in_box_lo =  70000.   90000.     0.
          erf.box1.in_box_hi = 110000.  150000.  6000.

Note that all three components must be given.  Supplying only x and y is the documented
way to ask for the full depth of the domain (see :ref:`subsec:full-depth-refinement`),
which is the opposite of what is wanted here.

If no file is given for the finer level, omit ``erf.nc_init_file_1`` and use a refinement
indicator alone.  The finer level's terrain and state are then interpolated from its
parent, and its base state is built as described below.  This works whether or not the
refined region reaches the domain top.

The base state at each level
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The WRF reference state is a closed-form function of height -- piecewise in log-pressure
across a troposphere, an isothermal layer and an optional stratosphere -- built from the
six parameters ``T00``, ``P00``, ``TLP``, ``TISO``, ``TLP_STRAT`` and ``P_STRAT``.  ERF
evaluates that profile at each level's own cell-centered heights and then rebalances each
column so that the result satisfies :math:`dp_0/dz = -\rho_0 g` discretely on that
level's mesh.  Every level runs this same construction, which is what makes the base
states of the levels consistent with one another.

Two properties matter when the refined region does not reach the domain top:

-  The rebalance integrates **upward from the bottom of each column**, so the value in a
   cell depends only on the cells at or below it.  A refined region that stops partway up
   therefore gets exactly the base state it would have had if it had been refined to the
   model top -- the values in the cells it does contain are unchanged, to the last bit.

-  The six parameters are read from the **level 0** file only, and every finer level uses
   those same values.  The whole base state is built from those six numbers, so levels
   that disagreed about them could not have base states that agree.  This matters in
   practice because idealized and hand-built wrfinput files often declare these variables
   but leave them zero-filled, in which case they fall back to the ERF defaults.  If a
   nested file carries values that differ from level 0's, ERF prints a warning naming
   both sets and uses level 0's.

The base states of two levels that each have their own file will still not be identical
to one another, because each is built on its own terrain and a nested file resolves finer
topography than its parent.  That difference is set by the difference in surface elevation
between the files and is carried up the column as an essentially constant hydrostatic
offset; it does not grow with height.  It is the same difference a full-depth nest already
has across its lateral faces.

Restrictions
~~~~~~~~~~~~

-  The refinement ratio in z may be greater than 1, in which case the wrfinput data is
   interpolated onto the finer vertical grid as it is read (see
   :ref:`subsec:wrfinput-vertical-refinement`).  Note that ``amr.ref_ratio`` applies one
   ratio in all three directions, so it always refines in z as well; use
   ``amr.ref_ratio_vect`` with 1 in the third slot to refine horizontally only.

-  Refining in the vertical does not add information the wrfinput file does not carry.
   The extra layers are interpolated from the file's eta levels, so what the finer
   vertical grid buys is resolution for the *solution* to develop on, not a better
   resolved initial condition.

-  The PBL models (MYJ, MYNN2.5, MYNN-EDMF, YSU, MRF), the SHOC PBL model, and the
   column-integral derived quantities (``helicity``, ``precipitable``,
   ``max_reflectivity``, ``mucape``) all need entire columns and cannot be used on a
   level whose grids do not span the domain in z.  Set ``erf.pbl_type = "None"`` for
   that level and leave those variables out of the plotfile variable list.

-  For a level that has its own wrfinput file and does not reach the domain top, the
   vertical remap of the **state** is incomplete in the cells near the top of the refined
   region.  ERF remaps the state from the WRF vertical coordinate onto its own
   terrain-following grid by searching the column of geopotential heights in the file;
   both that column and the state being interpolated are stored on the level's own grids,
   so above the top of the refined region there is nothing to interpolate from and those
   cells keep the values read from the file.  ERF prints a warning naming the level when
   this happens.  The base state is not affected.  A level whose state is interpolated
   from its parent is not affected either.

-  ``erf.box1.in_box_lo`` / ``in_box_hi`` in z are interpreted on the nominal (undeformed)
   vertical grid, not on the terrain-following heights the wrfinput file supplies.  Check
   the box ERF reports (``Saving in 'boxes at level'``) and the ``z_phys`` field in the
   plotfile to confirm the refined region reaches the height you intended.

See :ref:`subsec:partial-depth-refinement` for the mesh-refinement side of this, and
:ref:`MeshRefinement` for refinement in general.

TKE Initialization
--------------------

    When a turbulence closure that uses prognostic TKE is active, ERF initializes
    TKE at startup to ``erf.tke_min`` (default ``1.e-6 m^2/s^2`` but can be read from the inputs file).
    This includes Deardorff LES, k-equation RANS, MYJ, MYNN2.5, MYNN-EDMF, and SHOC.

Workflows
--------------------
For a summary of initialization strategies for real-data simulations, see the table below.

.. list-table:: Simulation Workflows
   :header-rows: 1

   * -
     - Large-Scale Data (reanalysis, HRRR, ...)
     - Intermediate Processing
     - Weather Simulation
   * - WRF
     - Manual download
     - WPS + ``real.exe``
     - ``wrf.exe``
   * - WRF --> ERF
     - Manual download
     - WPS + ``real.exe``
     - ``erf_exec`` (init from wrfinput)
   * -
     - Manual download
     - ``ndown.exe``
     - ``erf_exec`` (init from wrfinput)
   * - WPS --> ERF
     - Manual download
     - WPS
     - ``erf_exec`` (init from metgrid)
   * - E3SM --> ERF
     - ``run_e3sm``
     -  *Under development*
     - ``erf_exec``
   * - ERF standalone
     - Python tools
     - Python tools *(under development)*
     - ``erf_exec``
