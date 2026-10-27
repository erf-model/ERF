.. role:: cpp(code)
   :language: c++

.. _sec:FireCoupling:

Fire-Atmosphere Coupling
========================

A wildland fire and the atmosphere drive each other: the wind spreads the
fire, and the fire's heat drives the wind. ERF provides the atmosphere side of
that coupling: where the fire grid sits, what the fire is handed each step,
and how its heat, water vapour and smoke enter the atmosphere. The fire
itself, the model that spreads the front and computes the heat it releases,
is a *fire engine* that plugs into this coupling. Fire engines live outside
this repository, each in its own *engine directory*; a user builds ERF plus
one or more engines in their own build directory against the same ERF source
tree, and ERF itself needs no change to add an engine. The coupling layer is
the common part every engine shares.

Building with the fire coupling
-------------------------------

The coupling layer is compiled with the CMake option ``ERF_ENABLE_FIRE=ON``
(``USE_FIRE = TRUE`` with GNU make); it is off by default. With it on, the
state carries one more passive scalar, the smoke tracer (``smoke`` in the
plotfile), whether or not a fire runs. Every checkpoint, from any build,
records the number of passive scalars in a small file ``ScalarLayout``, so
that a build configured differently refuses to restart from it instead of
reading the components shifted; a checkpoint without the file is read as one
scalar, so older checkpoints restart in builds without the fire coupling.

The compile definition is ``ERF_USE_FIRE``, like the other optional modules.

:cpp:`erf.fire.enable = true` needs an engine compiled into the executable;
without one the run stops at start-up and says how to add one. A build
without the fire coupling stops on :cpp:`erf.fire.enable = true` as well.

Building ERF with a fire engine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An engine directory provides two build files next to its sources:

* ``Make.package`` for GNU make. ``Exec/Make.ERF`` includes it, with the
  directory on ``VPATH`` and the include path and ``FIRE_ENGINE_DIR`` set to
  it, so it lists its sources (``CEXE_sources += ...``) and any library it
  links (``LIBRARIES += ...``).
* ``ERFFireEngine.cmake`` for CMake, defining ``erf_fire_engine_add(target)``.
  ERF includes the file and calls the function with its object library, with
  ``ERF_FIRE_ENGINE_DIR`` set to the directory; it adds the sources
  (``target_sources``), include directories and libraries.

With GNU make, from any directory:

.. code-block:: make

   ERF_HOME   = /path/to/ERF
   FIRE_HOME  = /path/to/fire-engine      # one or more engine directories
   USE_FIRE   = TRUE
   include $(ERF_HOME)/Exec/Make.ERF

With CMake:

.. code-block:: console

   cmake -S /path/to/ERF -B build -DERF_ENABLE_FIRE=ON \
         -DERF_FIRE_ENGINE_DIRS="/path/to/fire-engine;/path/to/other-engine"

Several engines can be compiled into one executable; :cpp:`erf.fire.engine`
then picks one at run time (with a single engine it is the default).
``Tests/FireEngineExample`` is a complete engine directory, a prescribed
growing disk of heat that ERF's CI builds both ways.

Coupling layer and fire engines
-------------------------------

The code is split into two parts, so that more than one fire model can be
coupled to ERF the same way.

* **The coupling layer** (``Source/FireCoupling``, the ``FireCoupler`` class)
  is everything on the atmosphere side. It checks that the mesh can carry a
  fire grid, builds that grid on the level :cpp:`erf.fire.anchor_level`
  selects and hands the fire the atmospheric state at the start of each of
  its steps. With ``"lagged"`` or ``"synchronous"`` coupling it averages the
  heat the fire reports onto the atmospheric columns, holds it for one step,
  and injects it at every Runge-Kutta stage of the next step; with
  :cpp:`erf.fire.smoke_enable` it injects the smoke tracer the same way, also
  with ``"passive"`` coupling, which injects no heat. It saves the heat it
  holds in the checkpoint, so a restarted run injects the same heat in its
  first step as an uninterrupted one. It also provides helpers a fire engine
  may use: a wind sampler (``ERF_FireWindExtract.H``) and the ground height
  of each atmospheric column (``ERF_FireSurfaceHeight.H``).
* **The fire engine** (``FireEngine`` in ``ERF_FireEngine.H``) is the fire
  model itself: it advances the front on the fire grid and reports, after
  each step, the sensible and latent heat flux into the air on the fire grid
  [W/m²], averaged over its step. It may write its own checkpoint state and
  fire plotfile, list fields for the ERF plotfile, and offer a level set, a
  reference wind and the heat release to other surface models.
  :cpp:`erf.fire.engine` selects it among the engines compiled in.

The fire grid refines the columns of its level by :cpp:`erf.fire.grid_ratio`
in x and y. Each fire box is an atmospheric box refined by the ratio, owned
by the same rank, and every x and y box length on that level must be
divisible by the ratio; the run stops at start-up naming the first box that
is not, with the inputs to change. The level must keep its boxes for the whole
run (on a level above 0, :cpp:`erf.regrid_int` must be off), its boxes must
form one rectangle, no box may be split in z, and the bottom boundary must be
a surface layer unless the engine says it reads nothing from one
(``FireEngine::needs_surface_layer``; the example engine and LSFIRE do not).

The engine reads, at the start of each step, the face velocities of the fire
grid's level (before the dynamical core with ``"passive"`` and ``"lagged"``
coupling, after it with ``"synchronous"``), the cell-centre heights, the air
temperature and relative humidity of the lowest cell before the dynamical
core, and, when it asks for it, the surface precipitation of every column
from the microphysics accumulators (the run then needs a microphysics scheme
with precipitation). What it does with them, including the height it samples
the wind at and any correction it applies to that wind, is the engine's
choice.

Adding a fire engine
~~~~~~~~~~~~~~~~~~~~

``Source/FireCoupling/ERF_FireEngineTemplate.H`` is the starting point for a
new fire model. It declares every method a fire engine provides, each with a
comment on what ERF hands it, what ERF expects back and when it is called. It
holds declarations only: nothing in it has a body, no input selects it and
nothing creates one, so it adds no code to a run. ``ERF_FireCoupler.cpp``
checks at compile time that it declares every method ``FireEngine`` requires
(the pure virtual ones), and ``override`` checks the signature of every method
it declares; the optional methods are kept in it by hand.

To add an engine:

#. Copy ``Tests/FireEngineExample`` (or the template) into the new model's
   engine directory, outside ERF, rename the class, and write the method
   bodies there. Delete the optional methods the engine does not provide;
   their defaults then apply.
#. Read the model's inputs under its own prefix (``fire.<engine>.*`` is
   recommended), so they never mix with the coupling's keys
   (:cpp:`erf.fire.*`) or another engine's.
#. Register the engine in its ``.cpp``, at namespace scope:
   ``ERF_REGISTER_FIRE_ENGINE(MyEngine, "my_fire")``
   (``ERF_FireEngineRegistry.H``). The registration runs during static
   initialisation, so that ``.cpp`` must be compiled as source into the
   build, not hidden inside a prebuilt static library that nothing
   references; an engine shipping a prebuilt library keeps its registering
   wrapper as source. List the sources in the directory's ``Make.package``
   and ``ERFFireEngine.cmake``.
#. Run the ``FireCoupler`` unit tests (``Tests/Unit/Fire/ERF_GTestFireCoupler.cpp``,
   whose stand-in engine is the smallest engine that works) with the new
   engine, and a restart of a coupled run against an uninterrupted one: the
   atmosphere and fire plotfiles must agree at zero tolerance.

The engine owns its fire physics and its own state, including its checkpoint
and plotfile; the coupling layer does the rest, as described on this page.
An engine without restart support keeps the default ``read_checkpoint``,
which stops a restarted run with a message instead of continuing from the
engine's initial state.

The engine contract
~~~~~~~~~~~~~~~~~~~

**Fluxes.** ``fill_fluxes`` returns the sensible and latent heat flux into
the air, per unit ground area [W/m²], positive upward, on the fire grid.
Each is the time mean over exactly the interval of the last ``advance``,
:math:`[t_n, t_{n+1}]` with :math:`t_n` = ``forcing.time`` and
:math:`t_{n+1} - t_n` = ``forcing.dt``, however the engine steps inside it:
an engine that subcycles sums its substep fluxes times their steps and
divides by ``forcing.dt``. The coupling then injects that mean, unchanged, at
every Runge-Kutta stage of the next ERF step of the fire grid's level, over
:math:`[t_{n+1}, t_{n+2}]` (the one-step lag above), so the energy a column
receives over the run equals the energy the engine released, shifted by one
step.

**Water.** The latent flux is the energy of the water vapour the fire
releases. The coupling turns it into a vapour mass flux
:math:`F_q = Q_{lat} / L_v` [kg/m²/s], with :math:`L_v` the latent heat of
vaporisation of ``ERF_Constants.H``, and injects it into
:math:`\rho_d q_v` with the vertical profile of the heat. An engine that
computes a mass flux multiplies it by that same :math:`L_v`, so the water the
atmosphere receives is exactly the water the engine released.

**Memory.** Every MultiFab the coupling hands an engine and every one an
engine hands back is an ordinary AMReX MultiFab in the default arena, so on
a GPU build its data is on the device:

* the fields of ``FireForcing`` are temporaries of the ``advance`` call; an
  engine that needs them on the host copies them during the call
  (``fire_host_copy`` in ``ERF_FireUtils.H`` copies and synchronises the GPU
  stream);
* the engine owns every MultiFab it returns (fluxes, plot fields, getters),
  on the fire grid's ``BoxArray`` and ``DistributionMapping``, and keeps it
  alive and unchanged until its next ``advance``; a host-side engine copies
  its results into those device MultiFabs (``fire_host_copy`` the other way)
  before returning;
* an engine may hold managed or pinned memory internally, but synchronises
  the GPU stream before reading on the host what the device may still be
  writing;
* the coupling never frees, resizes or keeps across steps a pointer an
  engine returns.

**Shared helpers.** ``fill_fire_wind_from_interpolation``
(``ERF_FireWindExtract.H``) samples the 3D wind at a height above ground on
the fire grid; ``extract_fire_wind_k0`` gives the wind of the lowest cell of
each atmospheric column and its height above that column's ground, for an
engine driven by the lowest-cell wind; ``compute_fire_surface_height``
(``ERF_FireSurfaceHeight.H``) gives the ground of each fire cell's column.

Fire fields in the ERF plotfile
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An engine may list fields (``plot_var_names``, names starting with
``fire_``) that the deck can request in :cpp:`erf.plot_vars_1` and
:cpp:`erf.plot_vars_2`, like any ERF variable. The coupler writes the mean of
each over the fire cells under every atmospheric column into the lowest
cell of that column on the fire grid's level, and zero elsewhere. The engine's
own fire-resolution plotfile is separate (:cpp:`erf.fire.plot_int`,
:cpp:`erf.fire.plot_per`, which default to the 3D plotfile schedule).

Coupling modes
--------------

:cpp:`erf.fire.coupling_type` selects which wind the fire sees and whether its
heat returns to the atmosphere:

- ``"passive"``: the fire spreads on the wind from the start of the step and
  its heat flux is never injected. Use it for fast screening runs and for
  tests where the atmosphere should not respond.
- ``"lagged"`` (default): the fire spreads on the wind from the start of the
  step; the flux it produces is stored and injected during the next step's
  Runge-Kutta stages. This one-step explicit lag is the WRF-SFIRE approach
  (Mandel et al., 2011).
- ``"synchronous"``: the fire spreads on the wind *after* the dynamical core
  has advanced, so the spread at step :math:`n+1` already reflects the
  momentum response to the heat injected from step :math:`n`; the flux still
  enters at the next step.

The fire runs once per step of its level in all three modes.

Injection into the atmosphere
-----------------------------

The fire-grid fluxes are area-averaged onto the atmospheric columns of the
fire grid's level and distributed vertically with an exponential profile,

.. math::

   H(z) = \frac{Q}{c_p}\, e^{-z/\alpha_g}, \qquad
   \frac{\partial (\rho\theta)}{\partial t} = -\frac{\partial H}{\partial z},

where :math:`z` is height above the first cell centre of the column (the
profile's datum, half a cell above the terrain, so the lowest cell takes the
share :math:`1 - e^{-\Delta z_1/\alpha_g}` of the flux, with
:math:`\Delta z_1` the distance from the first cell centre to the second, and
nothing is lost below it) and
:math:`\alpha_g` = :cpp:`erf.fire.heat_flux_alfg` (default 45 m) is the
e-folding depth: 37% of the flux remains at 45 m and less than 1% at 225 m.
Each cell's tendency is the flux difference between the profile at its own
centre and at the next one, divided by the cell's own thickness (face to
face, from the nodal heights on a terrain-fitted or stretched grid; only a
degenerate cell of zero thickness is floored), so the
column receives :math:`(1 - e^{-z_{top}/\alpha_g})\,Q` on any vertical
spacing, where :math:`z_{top}` is the height of the top cell centre above the
first one plus the mean vertical spacing of the level (the profile closes one
mean spacing above the top centre). On a uniform grid, heights above the first
centre equal the heights of the cell faces above the ground, so each cell takes
the profile between its own faces; on a grid stretched by a constant ratio
:math:`r` the weights (apart from the top cell's) are those of the face-based profile with the e-folding
depth :math:`2\alpha_g/(1+r)`, which places slightly more of the heat in the
lowest cells (about 1 % of the flux moves down at :math:`r = 1.05`). The latent
flux is distributed the same way into the vapour equation (divided by the
latent heat of vaporisation) when :cpp:`erf.fire.inject_latent` is true and a
moisture model is active.

:cpp:`erf.fire.heat_profile = bottom_cell` replaces the profile: the whole
sensible and latent flux enters the lowest cell of each column, divided by
that cell's thickness, as in the LSFIRE coupling; :cpp:`erf.fire.heat_flux_alfg`
is then unused and :cpp:`erf.fire.heat_open_fraction` cannot be set. The
default, ``exponential``, is the profile above. Either way the heat enters
through the source term at every Runge-Kutta stage of the step after the
engine released it.

The tendency above is a :math:`\rho\theta` budget: heating :math:`q` per unit
volume gives :math:`\partial\theta/\partial t = q / (\rho c_p \Pi)` with the
Exner function :math:`\Pi = (p/p_0)^{R_d/c_p}`, so the enthalpy
:math:`c_p \Pi\, \partial(\rho\theta)/\partial t` the column receives is
:math:`\Pi` times the flux, 4.5 % short at 850 hPa.
:cpp:`erf.fire.heat_tendency_exner = true` divides the tendency by
:math:`\Pi` of each cell so that the enthalpy integrates to the flux; this
is the default. Setting it to false omits the factor, as WRF-SFIRE does. :cpp:`erf.fire.heat_tendency_density = true`
multiplies the tendency by the air density, so the column receives
:math:`\rho` times the flux (kept for comparison with runs made in that
form); set together, the two factors stack (:math:`\rho / \Pi`). With
:cpp:`erf.fire.fire_debug` the coupling prints the column-integrated heating
against the flux handed in, whose ratio is :math:`1 - e^{-z_{top}/\alpha_g}`
with neither factor (the factors multiply it by the profile-weighted mean of
:math:`\rho`, of :math:`1/\Pi`, or of :math:`\rho/\Pi` with both, essentially
their values near the surface). The Exner factor applies only to sensible
heating. The latent flux is converted to a water-mass flux by dividing by
the latent heat of vaporisation :math:`L_v` of ``ERF_Constants.H`` and enters the conserved water-vapour source
without Exner scaling. The optional legacy density multiplier still scales
both the sensible and water-vapour tendencies.

:cpp:`erf.fire.fire_atm_feedback` multiplies both fluxes before injection;
zero gives one-way coupling with the fire still responding to the wind, one
gives the full feedback. At zero the coupling writes nothing into the source,
so with ``source_mode = "overwrite"`` the other theta sources survive there,
while any positive value replaces them.

The cell source is rebuilt by the atmosphere at every Runge-Kutta stage and
the fire tendency is applied after that rebuild. :cpp:`erf.fire.source_mode`
selects how: ``"overwrite"`` (default) replaces the potential-temperature and
vapour slots, which keeps the injected energy independent of the number of
stages but discards any other source already in those slots, such as the
immersed-forcing scalar relaxation, radiative heating or Rayleigh damping;
``"add"`` accumulates the fire tendency into the freshly rebuilt source, once
per stage, so nothing is lost and nothing is double counted. Use ``"add"``
whenever another source writes into those slots: immersed forcing applied on
the slow step (``erf.immersed_forcing_substep = false``, the anelastic
default), radiative heating, Rayleigh damping. With immersed forcing on the
acoustic substeps (the compressible default) the scalar relaxation never
meets the slow-step source and the two modes are identical, as they are with
no other source. Because the default would silently drop such a source, a run
that leaves :cpp:`erf.fire.source_mode` unset while the fire injects heat and
the atmosphere has radiation, :cpp:`erf.four_stream_radiation`,
:cpp:`erf.rayleigh_damp_T`, custom theta forcing or subsidence (or custom
moisture forcing with latent heat on), or slow-step immersed forcing stops at
start-up and names them; setting the key to either value is the explicit
choice.

The heat held between steps carries five ghost columns, holding the
neighbouring box's value (or the edge column's outside a non-periodic face),
for a reader that samples a halo around each tile.

Levels
~~~~~~

The sensible and latent tendencies and the smoke source go into the source of
the level the fire grid refines (:cpp:`erf.fire.anchor_level`, the finest
level by default) at every Runge-Kutta stage of that level's steps. Outside
the refined region, which has no fire grid, the held heat is zero. Coarser
levels receive the heat through average-down (``erf.coupling_type =
"TwoWay"``), which keeps the volume integral, so the heat on level 0 is the
heat injected on the fire's level. With the fire grid on a coarser level than
the finest, the same average-down replaces the heated cells under the finer
level with unheated ones and that heat is lost, and with ``"OneWay"`` the
finer level never sees it; the coupling warns at start-up whenever the level
is not the finest.

A restart that re-makes the level-0 grids (with
:cpp:`erf.regrid_level_0_on_restart`, or on more ranks than the checkpoint has
level-0 boxes) is not supported with a fire running and stops before
reading the checkpoint.

Heat placement around buildings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An engine that knows where buildings stand can report, per atmospheric column,
the open fraction :math:`f` of its footprint and the mean roof height :math:`H`
of the rest, above the ground. With :cpp:`erf.fire.heat_open_fraction` the
coupling places the heat accordingly: a cell's share of the exponential profile
is scaled by :math:`f` where the cell lies below the roof, by one above it, and
linearly across the cell that holds the roof (the cell's faces measured from
the ground), and the column is renormalised over the shares that remain,

.. math::

   w_k = \frac{\bigl(e^{-z_k/\alpha_g} - e^{-z_{k+1}/\alpha_g}\bigr)\, s_k}
              {\sum_j \bigl(e^{-z_j/\alpha_g} - e^{-z_{j+1}/\alpha_g}\bigr)\, s_j}
         \,\bigl(1 - e^{-z_{top}/\alpha_g}\bigr),

so the energy injected per column is what the plain profile injects (exactly,
with neither tendency factor; the factors weigh the lifted heat by the
:math:`\rho` or :math:`1/\Pi` of the cells it moves to). With :math:`B` the
plain profile's share in the cells below the roof and
:math:`L = 1 - e^{-z_{top}/\alpha_g}`, the share placed below the roof falls
from :math:`B` to :math:`f B L / (L - (1 - f) B)`, more than :math:`f B`
because the renormalisation also raises the cells below the roof; the rest is
lifted above it. An engine
that reports no open fraction leaves the plain profile in place, with a
warning at start-up.

Smoke tracer
------------

With :cpp:`erf.fire.smoke_enable` (which needs :cpp:`erf.transport_scalar`)
the smoke tracer is fed at the surface with

.. math::

   F_{smoke} = \epsilon\, \frac{Q}{H_c} \quad [\mathrm{kg\,m^{-2}\,s^{-1}}],

where :math:`\epsilon` = :cpp:`erf.fire.smoke_emission_factor` (kg of smoke
per kg of fuel, default 0.02) and :math:`H_c` is the heat the engine releases
per kilogram of fuel burnt (``FireEngine::smoke_heat_per_kg``). The emission
uses the same held, area-averaged heat flux as the energy injection and is
placed in the lowest atmospheric cell. The tracer is advected and diffused
like any other scalar and is plotted as ``smoke``; its integral is printed
with the other integrated quantities. All passive scalars share one boundary
condition, so on a face that ingests the ``scalar`` boundary plane
(:cpp:`erf.bndry_input_var_names`) the smoke takes the base scalar's plane
values too.

Inputs
------

The coupling's inputs are listed under :ref:`Inputs <sec:InputsFire>`; an
engine's own inputs are its own.

References
----------

- Mandel, J., Beezley, J. D. and Kochanski, A. K. (2011). Coupled atmosphere-wildland fire modeling with WRF 3.3 and SFIRE 2011. Geoscientific Model Development, 4, 591-610.
