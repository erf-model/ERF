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

.. note::

   These entries were reconstructed from the source rather than from a
   specification, so the descriptions state what each input controls and where
   it is read.  The physical guidance -- which kernel or which time integrator
   to prefer for a given problem -- is not covered here.

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
``i = 0 ... num_initializations-1`` (``ERF_SuperDropletPCInitializations.cpp``).

+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| Parameter                                       | Definition                                               | Acceptable Values                                | Default                             |
+=================================================+==========================================================+==================================================+=====================================+
| **num_initializations**                         | how many initialization blocks are read                  | Integer > 0                                      | 1                                   |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **<i>.distribution_type**                       | shape of the seeded region                               | box, bubble                                      | box                                 |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **<i>.particles_per_cell**                      | super-droplets created per cell in that region           | Integer > 0                                      | set from the number density         |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **particle_box_lo**, **particle_box_hi**        | corners of the seeded box; read when                     | 3 Reals                                          | domain corners                      |
|                                                 | ``distribution_type`` = ``box``                          |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **particle_bubble_center**,                     | centre and radii of the seeded bubble; read when         | 3 Reals                                          | none                                |
| **particle_bubble_radius**                      | ``distribution_type`` = ``bubble``                       |                                                  |                                     |
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

Within a block, each condensate species and each aerosol is given its own
initial size or mass distribution, using keys named after that material --
``<species>_init_type``, and then the subset of ``<species>_mass_min``,
``<species>_mass_mean``, ``<species>_mass_max``, ``<species>_radius_min``,
``<species>_radius_mean``, ``<species>_radius_max`` and a standard deviation
(``_std``) or geometric standard deviation (``_gstd``) that the chosen type
requires.  Aerosol materials use the same pattern.  The exact key names depend
on the configured species and aerosol names; see ``ERF_SDInitialization.cpp``.

Injection blocks
----------------

Particles may also be injected during the run.
``super_droplets_moisture.num_injections`` sets how many injection blocks are
read, and block ``i`` uses the prefix ``super_droplets_moisture.injection.<i>``.

+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| Parameter                                       | Definition                                               | Acceptable Values                                | Default                             |
+=================================================+==========================================================+==================================================+=====================================+
| **num_injections**                              | how many injection blocks are read                       | Integer >= 0                                     | 0                                   |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **rate**                                        | injection rate of real droplets                          | Real > 0                                         | none                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **sd_rate**                                     | injection rate of super-droplets; an alternative to      | Real > 0                                         | derived from ``rate``               |
|                                                 | deriving it from ``rate``                                |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **min_multiplicity**                            | lower bound on the multiplicity of injected particles    | Real > 0                                         | none                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **t_start**, **t_stop**                         | simulation times [s] between which injection is active   | Real                                             | whole run                           |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **domain_velocity**                             | velocity [m/s] imposed on the injected particles         | 3 Reals                                          | 0.0 0.0 0.0                         |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **fractional_tol**                              | tolerance used when converting a fractional injection    | Real > 0                                         | none                                |
|                                                 | count into whole particles                               |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+

An injection block also accepts the same per-species and per-aerosol
distribution keys as an initialization block.
