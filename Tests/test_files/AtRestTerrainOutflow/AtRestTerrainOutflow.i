# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# A dry, constant-theta atmosphere at rest over a Witch of Agnesi hill, on a
# terrain-fitted mesh.  The initial state is exactly the hydrostatic base state
# with zero velocity and there is no forcing, no diffusion and no damping, so
# the exact solution is that nothing ever happens: w stays zero.
#
# The hill is wide enough (L is a third of the domain) that the terrain still
# has a slope where it meets the lateral boundaries, which is what makes the
# lateral ghost cells interesting: the nodal mesh is extrapolated past the
# domain, so a ghost cell sits at a different height than the cell just inside
# it, and the base state there has to be the same reference atmosphere sampled
# at that different height.
#
# The runner drives this deck twice, with the x boundaries symmetry and with
# them outflow, and compares max|w|.  Under symmetry the mesh extension is an
# exact mirror, the ghost cell and the cell it reflects sit at the same height,
# and the run is well balanced by construction -- it is the control.  The
# outflow run differs from it only in how the boundary is treated.
#
erf.prob_name = "Hydrostatic atmosphere at rest over terrain"

erf.init_type = Isentropic

max_step = 400

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     =   0.0   0.0    0.0
geometry.prob_hi     = 600.0  20.0  600.0
amr.n_cell           =  64     4     64

geometry.is_periodic = 0 1 0

# The runner overrides the two x boundaries
xlo.type = "Outflow"
xhi.type = "Outflow"
zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
# Fixed, and well inside the acoustic CFL, so that the two runs march in step
erf.substepping_type = None
erf.fixed_dt         = 0.005

# DIAGNOSTICS & VERBOSITY
erf.sum_interval = -1
erf.v            = 1
amr.v            = 1

# REFINEMENT / REGRIDDING
amr.max_level = 0

# CHECKPOINT FILES
erf.check_int = -1

# PLOTFILES
erf.plot_file_1 = plt
erf.plot_int_1  = 400
erf.plot_vars_1 = density x_velocity y_velocity z_velocity pressure theta pert_pres pert_dens z_phys

# SOLVER CHOICE -- nothing that could damp or drive the flow
erf.use_gravity  = true
erf.use_coriolis = false
erf.les_type     = "None"

erf.molec_diff_type = "None"

# TERRAIN GRID TYPE
erf.terrain_type      = StaticFittedMesh
erf.terrain_smoothing = 2 # Sullivan TF

# PROBLEM PARAMETERS
prob.custom_terrain_type = "WoA"
prob.hmax  = 100.0
prob.L     = 200.0
prob.T_0   = 300.0
prob.rho_0 = 1.16
