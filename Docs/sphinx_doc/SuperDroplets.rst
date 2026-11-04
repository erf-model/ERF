.. _sec:SuperDroplets:

Super-Droplets
==============

The super-droplet method represents the condensed phase with Lagrangian
particles rather than with bulk mixing ratios.  Each particle carries a
multiplicity -- the number of real droplets it stands for -- together with the
mass of each condensate species and of each aerosol it contains.  The method is
selected with ``erf.moisture_model = SuperDroplets``; see
:ref:`Microphysics` for where it sits among the other moisture models, and
:ref:`sec:Plotfile3DReference` for the per-species plotfile names it provides.

All of the inputs below use the prefix ``super_droplets_moisture``, which is the
name of the particle container (``ERF_SuperDropletsMoist.H``).  Inputs that the
container inherits from the generic ERF particle container use the same prefix;
they are listed under :ref:`inputs-particles` for tracer particles and behave
the same way here.

.. _superdroplets-anelastic-compatibility:

Anelastic compatibility
-----------------------

``SuperDroplets`` uses a thermodynamic pressure source selected independently
on each AMR level. On a level with ``erf.anelastic = 1``, the pressure supplied
to the super-droplet microphysics is the hydrostatic base-state pressure at
each cell, rather than a pressure diagnosed from the compressible equation of
state. ERF diagnoses the temperature from dry potential temperature using

.. math::

   \theta = \frac{\rho_d\theta}{\rho_d}, \qquad
   p_{\mathrm{SDM}} = p_{\mathrm{base}}, \qquad
   T = \theta\left(\frac{p_{\mathrm{base}}}{p_{\mathrm{ref}}}\right)^{R_d/c_p}.

Here :math:`\rho_d` is dry-air density, :math:`p_{\mathrm{base}}` is the
cellwise hydrostatic base-state pressure, :math:`p_{\mathrm{ref}}=10^5` Pa is
the fixed reference pressure used to define potential temperature, and
:math:`c_p` is the configured ``erf.c_p``. These pressures are in pascals.
No additional pressure or Exner-function input is required.

During an anelastic SDM microphysics source step, the pressure remains the
base-state pressure while the temperature is refreshed from the updated
potential temperature after each phase-change substep. Latent heating uses
the fixed inverse Exner factor for every process, including sequential
processes within a substep:

.. math::

   \Delta\theta =
   \left(\frac{p_{\mathrm{ref}}}{p_{\mathrm{base}}}\right)^{R_d/c_p}
   \frac{L}{c_p}\,\Delta q.

Here :math:`\Delta q` is the signed, dry-air-referenced mass transfer
and :math:`L` is the corresponding latent heat. This pressure and temperature
are supplied consistently to the scheme's saturation, phase-change,
terminal-velocity, and collision calculations. On compressible levels, SDM
diagnoses pressure and temperature from the local equation of state and uses
the local potential-temperature-to-temperature ratio for latent heating.

For example, with a particles-enabled build:

.. code-block:: text

   erf.anelastic = 1
   erf.moisture_model = SuperDroplets

Choose compatible mesh, solver, particle, and AMR options, and validate the
selected configuration before using results for quantitative cloud-microphysics
analysis.

.. warning::

   Anelastic ``SuperDroplets`` execution has not been validated on CUDA or HIP
   GPUs. Check the selected particle processes, restart behavior, and AMR
   configuration before using GPU results for scientific analysis; do not
   assume reproducibility between CPU and GPU runs.

Species transport and diffusion
-------------------------------

In coupled atmospheric simulations, each additional species specified by
``super_droplets_moisture.species`` has an Eulerian vapor mixing ratio and a
condensate mixing ratio diagnosed from the super-droplets.
ERF advects the vapor with the resolved flow and applies the
configured turbulent diffusivity. Vertical turbulent diffusion includes an
implicit contribution when ``erf.implicit_moisture_diffusion`` is enabled and
the stage factor in ``erf.vert_implicit_fac`` is nonzero. These additional
species have no molecular scalar diffusion, even when molecular diffusion is
enabled for water or passive scalars.

The condensate field is a diagnostic of particle-carried mass, not a second
mass reservoir. It is not independently advected, diffused, or positivity-clipped
by the atmospheric solver. Particle motion and microphysical processes determine
its evolution. In particular, vapor diffusion alone does not represent turbulent
dispersion of condensed particles. These rules apply to the additional species;
the existing water microphysics coupling is unchanged.

For vertical implicit diffusion, the supported species boundary conditions are
zero-gradient extrapolation and prescribed Neumann gradients. If a transported
species has an unsupported vertical boundary condition, ERF issues a warning
and disables the shared moisture implicit-diffusion option. Horizontal diffusion,
explicitly treated stages, and terrain mixed terms can still constrain the
timestep. An implicit vertical solve is not an unconditional positivity or
stability guarantee for the full time integration.

AMR and particle sampling
~~~~~~~~~~~~~~~~~~~~~~~~~

Particles are assigned to the finest available AMR level containing their
positions. Newly injected particles are assigned before phase change or
condensate projection, including when ``include_advection = false``.
This assignment does not require ``split_merge_amr = true``.

Initial populations are seeded only in cells not covered by finer grids, so
refining a region does not add another copy of its prescribed particle mass.
Density scaling and condensate-based attribute initialization are completed
on each level before redistribution. If initial phase-change relaxation is
enabled, it then runs on each particle's owning level and retains the associated
vapor and latent-heat changes in the initial atmospheric state. These
initialization operations are not repeated when restarting from a checkpoint.

A restart restores the condensate fields from the checkpoint, including
fine-level averages in covered coarse cells. Turbulent diffusivities, including
``Kmv``, are not checkpointed. They can be zero in the initial restart plotfile
until the next turbulence update, even when the conserved state has been
restored correctly.

For additional species, cloud-in-cell projection uses only cells on the
particle's level that are not covered by finer grids. The retained weights are
renormalized, so projection at a physical wall or coarse/fine boundary preserves
the particle's represented mass. Near these boundaries the stencil is one-sided;
it does not transfer condensate mass between levels independently of particle
motion. Diagnostic condensate is refreshed after initialization, time-step
redistribution, and regridding without advancing phase change or vapor transport.

With ``split_merge_amr = true``, ERF also adjusts the super-droplet sampling as
particles change levels. When empty fine cells are populated from a neighboring
donor, the donor and its daughters share the original multiplicity. Each retained
daughter represents at least one real droplet; insufficient multiplicity leaves
some cells empty rather than creating physical mass. The allocation uses a fixed
cell ordering for a fixed particle population and grid decomposition. This is
not a guarantee of identical particle distributions across decompositions.

Condensate output and mass budgets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The output fields ``qv_<species>`` and ``qc_<species>`` are mixing ratios per unit
dry-air mass (:math:`\mathrm{kg\,kg^{-1}}`), whereas
``super_droplets_moisture_species_mass_density_<species>`` is particle-derived
condensate mass per physical volume (:math:`\mathrm{kg\,m^{-3}}`). Thus the
condensate density corresponding to ``qc_<species>`` is
:math:`\rho q_{c,s}`, not :math:`q_{c,s}` alone; :math:`\rho` is ERF's dry-air
density (the ``density`` output field).

For a three-dimensional mass budget, sum only uncovered cells of the AMR
hierarchy, using their physical volumes :math:`V_c`:

.. math::

   M_s = \sum_{\mathrm{composite\ cells}} V_c
         \left(\rho q_{v,s} + \rho q_{c,s}\right),
   \qquad
   M_{c,s} = \sum_{\mathrm{active\ particles}} \xi_p m_{p,s}.

Here :math:`\xi_p` is the number of real droplets represented by super-droplet
:math:`p`, and :math:`m_{p,s}` is the mass of species :math:`s` in one real
droplet. The particle sum is an independent measure of condensate mass:
do not add it to the mesh condensate term, which represents the same mass.
Phase change exchanges vapor and condensate. Account separately for injection,
precipitation, and other physical sources, sinks, or boundary fluxes when
assessing conservation; an injection run is not a closed system.

With ``dimensionality = one_d_z``, species condensate density from the reference
column (the lowest horizontal indices in each grid box) is replicated
horizontally before conversion to mixing ratio. Both normal microphysics
updates and diagnostic refreshes use this convention. It is not a horizontal
average: a volume integral over the replicated field must be interpreted using
the represented column area, rather than compared directly with an unreplicated
particle-column sum. The ``three_d``, ``two_d_xz``, and ``two_d_yz`` species
projections do not apply this extra horizontal replication.

Inflow limitations
~~~~~~~~~~~~~~~~~~

There is no per-additional-species prescribed inflow interface. Ordinary inflow
with a specified density uses a zero external value for these species;
``qv`` and ``qc`` inflow settings prescribe water only.
When inflow density is not prescribed, the existing fallback extrapolates
non-ingested conserved fields from the interior. Boundary-plane ingestion
does not supply additional species. Do not interpret water-vapor boundary
settings as non-water inflow concentrations.

Model configuration
-------------------

These inputs select which processes the model includes and what it carries.
They are read by the moisture model itself (``ERF_SuperDropletsMoistInit.cpp``)
rather than by the particle container, but use the same prefix.

+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| Parameter                                       | Definition                                               | Acceptable Values                                | Default                             |
+=================================================+==========================================================+==================================================+=====================================+
| **species**                                     | additional condensate species carried besides water, and | List of species names                            | none                                |
|                                                 | ice when cold processes are on.  Naming a species that   |                                                  |                                     |
|                                                 | is already modelled aborts the run, and at most five     |                                                  |                                     |
|                                                 | extra species may be given                               |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **aerosols**                                    | aerosol materials carried by the super-droplets          | List of species names                            | none                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **include_phase_change**                        | include condensation and evaporation                     | Boolean                                          | true                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **include_advection**                           | advect the super-droplets                                | Boolean                                          | true                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **include_coalescence**                         | include coalescence                                      | Boolean                                          | true                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **include_cold_processes**                      | include ice; adds ice to the species list                | Boolean                                          | true                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **distribution_type**                           | how the initial condensate is distributed over the       | uniform, condensate_density                      | uniform                             |
|                                                 | particles                                                |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **radius_raindrop**                             | radius [m] above which a drop is counted as rain rather  | Real > 0                                         | 4.0e-5                              |
|                                                 | than cloud water                                         |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **rime_mass_ratio**                             | rime mass fraction above which a particle is treated as  | Real                                             | 0.3                                 |
|                                                 | rimed                                                    |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **kinematic_mode**                              | run with a prescribed flow field rather than a coupled   | Boolean                                          | false                               |
|                                                 | one                                                      |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **dimensionality**                              | dimensionality of the super-droplet dynamics             | one_d_z, two_d_xz, two_d_yz, three_d             | three_d                             |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **recycle_particles**                           | return deactivated particles to the domain instead of    | Boolean                                          | false                               |
|                                                 | dropping them                                            |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **num_substeps_phase_change**                   | sub-steps taken per time step for the phase-change       | Integer > 0                                      | 1                                   |
|                                                 | update                                                   |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **initial_phase_change_relaxation**             | relax the particles to an equilibrium size at            | Boolean                                          | false                               |
|                                                 | initialization                                           |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **initial_phase_change_relaxation_time**        | duration [s] of that relaxation                          | Real > 0                                         | 10.0                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **diagnostics_interval**                        | number of steps between writes of the distribution       | Integer > 0                                      | 1                                   |
|                                                 | diagnostics                                              |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+

The species and aerosol names given here determine the per-material keys used by
the initialization and injection blocks below, and the ``qv_<species>`` and
``qc_<species>`` plotfile names.

Solver and time integration
---------------------------

The condensational growth of each particle is an ODE in the droplet mass, solved
per particle per step.

The implicit Newton solve uses the larger of ``newton_solver_stol`` and the
particle arithmetic's machine epsilon as its relative step stopping criterion.
This avoids rejecting a physically resolved update when its remaining Newton
correction is at roundoff. Requesting a tolerance below machine epsilon does
not increase the arithmetic precision of the solve.

+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| Parameter                                       | Definition                                               | Acceptable Values                                | Default                             |
+=================================================+==========================================================+==================================================+=====================================+
| **mass_change_ti_method**                       | time integrator for the mass-change ODE.  An             | backward_euler, crank_nicolson, dirk2, rk3bs,    | backward_euler                      |
|                                                 | unrecognized value aborts the run                        | rk4                                              |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **mass_change_cfl**                             | CFL-like limit on the mass-change sub-step               | Real > 0                                         | 1000.0                              |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **newton_solver_rtol**                          | relative tolerance of the Newton solve for the implicit  | Real > 0                                         | 1.0e-6                              |
|                                                 | integrators                                              |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **newton_solver_atol**                          | absolute tolerance of that solve                         | Real >= 0                                        | 1.0e-99 (0.0 in single precision)   |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **newton_solver_stol**                          | step tolerance of that solve                             | Real > 0                                         | 1.0e-12                             |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **newton_solver_maxits**                        | maximum Newton iterations per particle                   | Integer > 0                                      | 10                                  |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **mass_change_unconverged_log**                 | write a log entry for every particle whose mass-change   | Boolean                                          | false                               |
|                                                 | solve fails to converge                                  |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **mass_change_unconverged_log_filename**        | name of that log file                                    | String                                           | unconverged_superdroplets.log       |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+

Ventilation
-----------

The ventilation factor enhances the mass-change rate of a falling drop.  It
follows Bayley et al. (2025), Eq. 3, and is off by default.

+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| Parameter                                       | Definition                                               | Acceptable Values                                | Default                             |
+=================================================+==========================================================+==================================================+=====================================+
| **mass_change_ventilation**                     | include the ventilation factor in the phase change       | Boolean                                          | false                               |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **ventilation_alpha1**                          | fit coefficient alpha_1                                  | Real                                             | 6.954e7                             |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **ventilation_beta1**                           | fit exponent beta_1                                      | Real                                             | 1.963                               |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **ventilation_alpha2**                          | fit coefficient alpha_2                                  | Real                                             | 1.069e3                             |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **ventilation_beta2**                           | fit exponent beta_2                                      | Real                                             | 0.702                               |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **ventilation_fcap**                            | upper bound applied to the ventilation factor            | Real > 0                                         | 20.0                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+

Coalescence
-----------

+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| Parameter                                       | Definition                                               | Acceptable Values                                | Default                             |
+=================================================+==========================================================+==================================================+=====================================+
| **coalescence_kernel**                          | collision kernel used for coalescence.  An unrecognized  | sedimentation, golovin, Longs, Halls             | sedimentation                       |
|                                                 | value aborts the run.  Note that ``Longs`` and           |                                                  |                                     |
|                                                 | ``Halls`` are matched case-sensitively                   |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **kernel_relative_velocity**                    | which relative velocity enters the kernel                | terminal_velocity, absolute_velocity,            | terminal_velocity                   |
|                                                 |                                                          | radial_velocity                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **include_brownian_coalescence**                | add the Brownian contribution to the kernel              | Boolean                                          | false                               |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **ice_aggregation_efficiency**                  | collection efficiency prefactor for ice-ice aggregation  | Real                                             | 0.1                                 |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **coalescence_bin_size**                        | size, in cells, of the bins within which candidate       | 3 Integers > 0                                   | 1 1 1                               |
|                                                 | collision pairs are drawn                                |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **terminal_velocity_model**                     | terminal-velocity model.  One value sets the water       | RogersYau, AtlasUlbrich, CloudRainShima,         | CloudRainShima IceBohm              |
|                                                 | model; a second, if given, sets the ice model.  Giving   | IceBohm                                          |                                     |
|                                                 | more than two values aborts the run                      |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+

Particle population and motion
------------------------------

+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| Parameter                                       | Definition                                               | Acceptable Values                                | Default                             |
+=================================================+==========================================================+==================================================+=====================================+
| **nucleate_particles**                          | create new super-droplets from the vapour field          | Boolean                                          | false                               |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **density_scaling**                             | scale the initial number density with the local air      | Boolean                                          | false                               |
|                                                 | density                                                  |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **place_randomly_in_cells**                     | place new particles at random positions within each cell | Boolean                                          | true                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **prescribed_advection**                        | move the particles with a prescribed vertical velocity   | Boolean                                          | false                               |
|                                                 | instead of the resolved flow                             |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **prescribed_w**                                | the prescribed updraft [m/s] used in that mode; a value  | Real                                             | 0.0                                 |
|                                                 | of zero selects a sinusoidal pulse instead of a constant |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **split_merge_amr**                             | split and merge particles as they cross AMR level        | Boolean                                          | false                               |
|                                                 | boundaries                                               |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **sigma0**                                      | width parameter of the assumed mass distribution         | Real > 0                                         | 0.62                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+

Deactivation and recycling
--------------------------

A particle whose multiplicity falls below a threshold is deactivated and may be
returned to the domain at a new position.

+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| Parameter                                       | Definition                                               | Acceptable Values                                | Default                             |
+=================================================+==========================================================+==================================================+=====================================+
| **inactive_threshold**                          | multiplicity below which a particle is deactivated       | Real > 0                                         | 0.01                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **write_inactive_plt**                          | include deactivated particles in the particle plotfiles  | Boolean                                          | false                               |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **recycle_xmin**, **recycle_xmax**              | x-bounds of the region into which recycled particles are | Real                                             | domain bounds                       |
|                                                 | returned                                                 |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **recycle_ymin**, **recycle_ymax**              | y-bounds of that region                                  | Real                                             | domain bounds                       |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **recycle_zmin**, **recycle_zmax**              | z-bounds of that region                                  | Real                                             | domain bounds                       |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+

Diagnostics
-----------

+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| Parameter                                       | Definition                                               | Acceptable Values                                | Default                             |
+=================================================+==========================================================+==================================================+=====================================+
| **distribution_grid_size**                      | number of points in the one-dimensional radius grid on   | Integer > 0                                      | 100                                 |
|                                                 | which size distributions are accumulated                 |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **distribution_rmin**                           | smallest radius [m] of that grid.  Read only when built  | Real > 0                                         | 1.0e-6                              |
|                                                 | with ``ERF_USE_ML_UPHYS_DIAGNOSTICS``                    |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **distribution_rmax**                           | largest radius [m] of that grid.  Read only when built   | Real > 0                                         | 5.0e-3                              |
|                                                 | with ``ERF_USE_ML_UPHYS_DIAGNOSTICS``                    |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+

.. _sec:SuperDropletInitialization:

Initialization blocks
---------------------

The initial particle population is described by one or more numbered blocks.
``super_droplets_moisture.num_initializations`` sets how many there are, and
block ``i`` is read under the prefix ``super_droplets_moisture.<i>`` for
``i = 0 ... num_initializations-1``. A block-specific value overrides the same
key under ``super_droplets_moisture``; otherwise that common value is used.

+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| Parameter                                       | Definition                                               | Acceptable Values                                | Default                             |
+=================================================+==========================================================+==================================================+=====================================+
| **num_initializations**                         | how many initialization blocks are read                  | Integer >= 0                                     | 1                                   |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **<i>.initial_distribution_type**               | shape of the seeded region                               | uniform, bubble                                  | uniform                             |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **<i>.initial_particles_per_cell**              | super-droplets created per cell in that region           | Integer > 0                                      | set from the number density         |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **particle_box_lo**, **particle_box_hi**        | corners of the seeded box; read when                     | 3 Reals                                          | domain corners                      |
|                                                 | ``initial_distribution_type`` = ``uniform``              |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **particle_bubble_center**,                     | centre and radii of the seeded bubble; read when         | 3 Reals                                          | none                                |
| **particle_bubble_radius**                      | ``initial_distribution_type`` = ``bubble``               |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **initial_number_density**                      | number density of real droplets represented              | Real > 0                                         | none                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **initial_super_droplet_density**               | number density of super-droplets used to represent them  | Real > 0                                         | none                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **maximum_multiplicity**                        | cap on the multiplicity assigned to a super-droplet      | Real > 0                                         | none                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **multiplicity_type**                           | how multiplicity is distributed across the population    | see ``ERF_SDInitialization.cpp``                 | none                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **ice_apparent_density**                        | apparent density [kg/m^3] used for ice particles         | Real > 0                                         | none                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+

Each material has its own size or mass distribution. For example,
``initial_species_distribution_type_<species> = mass_constant`` uses
``initial_species_mean_mass_<species>`` in kg per real droplet, not per
super-droplet. Other supported distributions are ``mass_exponential``,
``radius_log_normal``, and ``radius_lognormal_autorange``. Their parameters use
the same prefix with ``min_mass``, ``max_mass``, ``min_radius``, ``mean_radius``,
or ``max_radius`` in place of ``mean_mass``, as appropriate. Radii are in metres.
For log-normal radii, specify either ``initial_species_std_radius_<species>``
(the standard deviation of log radius) or
``initial_species_geomstd_radius_<species>`` (the geometric standard deviation),
not both. Aerosols use ``initial_aerosol_`` in place of ``initial_species_``.
Material names must match the configured species and aerosol names.

Injection blocks
----------------

Particles may also be injected during the run.
``super_droplets_moisture.num_injections`` sets how many injection blocks are
read, and block ``i`` uses the prefix ``super_droplets_moisture.injection.<i>``.
Block-specific values override common injection values under
``super_droplets_moisture.injection``. The rate is a volumetric number source
within the injection region, not a total particle count per second for the domain.
The region may lie inside fine coverage or span a coarse/fine boundary;
particles are injected once and then assigned to their finest available level.

+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| Parameter                                       | Definition                                               | Acceptable Values                                | Default                             |
+=================================================+==========================================================+==================================================+=====================================+
| **num_injections**                              | how many injection blocks are read                       | Integer >= 0                                     | 0                                   |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **rate**                                        | real-droplet number source [m^-3 s^-1]                   | Real > 0                                         | 0 (no injection)                    |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **sd_rate**                                     | super-droplet number source [m^-3 s^-1]                  | Real > 0                                         | derived from rate and multiplicity  |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **min_multiplicity**                            | target minimum multiplicity in per-box injection mode    | Real > 0                                         | 100                                 |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **t_start**, **t_stop**                         | simulation times [s] between which injection is active   | Real                                             | whole run                           |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **domain_velocity**                             | translation velocity [m/s] of the injection region       | 3 Reals                                          | 0.0 0.0 0.0                         |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **fractional_tol**                              | tolerance used when converting a fractional injection    | Real > 0                                         | 1.0e-3                              |
|                                                 | count into whole particles                               |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+

Within an injection block, use ``distribution_type = uniform`` for a box
specified by ``particle_box_lo`` and ``particle_box_hi``, or ``bubble`` with
``particle_bubble_center`` and ``particle_bubble_radius``. Material distribution
keys follow the initialization pattern without ``initial_``; for example,
``species_distribution_type_<species>`` and ``species_mean_mass_<species>``.

Specifying ``particles_per_cell`` selects per-cell injection. Otherwise, ERF
accumulates a super-droplet count over the injection region using ``sd_rate``
or ``rate / min_multiplicity``. If both are explicitly specified, the smaller
super-droplet rate is used. Physical-particle counts and integer multiplicities
are accumulated or rounded according to the selected mode, so low rates need
not produce an injection event every time step. Include the actual represented
injected mass in conservation diagnostics rather than assuming a continuously
realized source at every step.
