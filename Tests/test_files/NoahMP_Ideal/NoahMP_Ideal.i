# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Regression motivation:
# The Noah-MP land-surface model from an idealized ERF case. It is the only test that
# builds and runs Noah-MP at all (it needs ERF_ENABLE_NOAHMP, which needs NetCDF).
#
# Three things it guards:
#
# 1. The radiation inputs. No erf.radiation_model is set, so nothing writes the
#    downwelling shortwave, downwelling longwave or solar zenith angle Noah-MP reads;
#    they hold the lsm_undefined sentinel NOAHMP::Init filled them with (1e150). Taken
#    raw, Noah-MP integrated on ~2e149 W/m^2 of shortwave until its own energy-budget
#    check ended the run with a Fortran STOP. ERF now replaces an invalid input with zero.
#
# 2. The exit status. That Fortran STOP exits with status 0, so the run died at the
#    first land step with no final plotfile. Launched directly it reported success
#    (Open MPI's launcher does flag it, since the rank never calls MPI_Finalize). ERF
#    now turns an exit() that happens while AMReX is still running into a failure, and
#    the checker requires the last plotfile, so neither depends on the launcher.
#
# 3. The outputs. The checker requires Noah-MP's own fields in physical ranges -- the
#    absorbed shortwave (sav, sag) in particular, which must be zero here.
#
# With no radiation reaching the surface the ground radiates to a 0 K sky and cools
# quickly (about 300 K to 253 K in the one Noah-MP hour this runs). That is the
# correct response to the configuration, not a Noah-MP defect; the start-up warning
# says the same.
#
# The land state comes from wrfinput_d01, which the test generates from
# wrfinput_ideal.cdl with ncgen; namelist.erf names it. NoahmpTable.TBL is copied from
# the Noah-MP submodule by the test harness.

max_step = 2
stop_time = 100.0

amrex.fpe_trap_invalid = 0

geometry.prob_extent = 1000 1000 1000
amr.n_cell           = 4 4 32
geometry.is_periodic = 1 1 0
amr.max_level = 0

# Noah-MP applies its fluxes to the atmosphere through the surface layer, which in
# turn needs a diffusive closure.
zlo.type = "surface_layer"
zhi.type = "SlipWall"
erf.most.z0   = 0.1
# Above the first cell centre (dz = 31.25 m, so 15.6 m); MOSTAverage requires it.
erf.most.zref = 20.0

erf.fixed_dt = 1.0
erf.v = 1

erf.prob_name = "ABL"
erf.init_type = "input_sounding"
erf.sounding_type = Ideal
erf.input_sounding_file = "input_sounding"

erf.use_gravity = true
erf.les_type = "Smagorinsky"
erf.Cs = 0.16
erf.pbl_type = "None"
erf.molec_diff_type = "None"

start_datetime = "2024-08-05 12:00:00"

erf.land_surface_model = "NOAHMP"

erf.check_int = -1
erf.plot_file_1 = plt
erf.plot_int_1 = -1

# Noah-MP's own outputs, which the checker reads.
erf.plot2d_file_1 = plt2d
erf.plot2d_int_1 = 2
erf.plot2d_vars_1 = t_sfc sav sag sensible_heat_flux grdflx fira
