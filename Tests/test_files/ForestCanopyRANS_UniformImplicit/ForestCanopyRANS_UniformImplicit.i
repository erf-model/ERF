# ForestCanopyRANS_Uniform with erf.implicit_tke_dissipation = true, the
# setting the docs recommend with the canopy terms. Each RK3 stage then divides
# the explicit update of k by 1 + dt_s Cmu0^3 k0^0.5 / L:
#   k_s = [k0 e^-x + (P / c)(1 - e^-x)] / (1 + dt_s Cmu0^3 k0^1.5 / (L k0)),
# with the notation of ForestCanopyRANS_Uniform.i. One step gives
# u = 4.8780542453793183 and k = 1.0110471279453457; Tests/CTestList.cmake
# checks z >= 45 m.
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
erf.implicit_tke_dissipation = true

erf.init_type = "uniform"
prob.rho_0 = 1.2
prob.T_0   = 300.0
prob.U_0   = 5.0
prob.V_0   = 0.0
prob.W_0   = 0.0
prob.KE_0  = 0.5

erf.forest_file       = forest_uniform.txt
erf.forest_canopy_tke = true

erf.fixed_dt           = 0.1
erf.fixed_mri_dt_ratio = 4

erf.v         = 1
erf.check_int = -1
erf.plot_file_1 = plt
erf.plot_int_1  = 1
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta KE
