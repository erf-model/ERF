# One step of the k-eqn RANS canopy terms in a horizontally uniform, fully
# forested periodic box without gravity. The forest (forest_uniform.txt) covers
# the whole domain with a 1000 m tree, so f = C_d LAI / h = 0.2 * 250 / 1000
# = 0.05 1/m in every cell, and the canopy cap L = alpha / f = 0.8 m is below the
# geometric length everywhere. The cap needs erf.implicit_tke_dissipation = true.
# With uniform wind and k, every cell above the reach of the bottom wall follows
# the same ODEs:
#   du/dt = -f u |u|
#   dk/dt = beta_p f |u|^3 - beta_d f |u| k - Cmu0^3 k^1.5 / L
# Each RK3 stage (dt_s = dt/3, dt/2, dt, all from the start of the step) takes
# the canopy terms exactly with production and sink frozen at the stage wind,
# and divides the result by the implicit dissipation factor:
#   k_s = [k0 e^-x + (P / c)(1 - e^-x)] / (1 + dt_s Cmu0^3 k0^0.5 / L),
# x = c dt_s, P = beta_p f |u_s|^3, c = beta_d f |u_s|, and the dissipation uses
# k0 = 0.500001 (prob.KE_0 plus erf.tke_min = 1e-6, added at initialisation) as the
# closure is evaluated once per step. One step gives u = 4.8780542453793183 and
# k = 1.0110471279453457. The bottom wall, which the RANS wall distance needs,
# reaches three cells in one step; Tests/CTestList.cmake checks z >= 45 m.
erf.prob_name = "ABL"
max_step = 1
amrex.fpe_trap_invalid = 1
fabarray.mfiter_tile_size = 1024 1024 1024

geometry.prob_extent = 40 40 160
amr.n_cell           =  4  4 16
geometry.is_periodic = 1 1 0
amr.max_level        = 0

zlo.type = "NoSlipWall"
zhi.type = "SlipWall"

erf.use_gravity     = false
erf.molec_diff_type = "None"
erf.les_type        = "None"
erf.rans_type       = "kEqn"
erf.vert_implicit   = false

erf.init_type = "uniform"
prob.rho_0 = 1.2
prob.T_0   = 300.0
prob.U_0   = 5.0
prob.V_0   = 0.0
prob.W_0   = 0.0
prob.KE_0  = 0.5

erf.forest_file       = forest_uniform.txt
erf.forest_canopy_tke = true
erf.implicit_tke_dissipation = true

erf.fixed_dt           = 0.1
erf.fixed_mri_dt_ratio = 4

erf.v         = 1
erf.check_int = -1
erf.plot_file_1 = plt
erf.plot_int_1  = 1
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta KE
