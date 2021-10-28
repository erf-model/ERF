
// This file is automatically created by parse_erf_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh

pp.query("state_interp_order", state_interp_order);
pp.query("lin_limit_state_interp", lin_limit_state_interp);
pp.query("do_reflux", do_reflux);
pp.query("do_avg_down", do_avg_down);
pp.query("load_balance_verbosity", load_balance_verbosity);
pp.query("CI", CI);
pp.query("PrT", PrT);
pp.query("bndry_func_thread_safe", bndry_func_thread_safe);
pp.query("job_name", job_name);

pp.query("use_state_advection", use_state_advection);
pp.query("use_momentum_advection", use_momentum_advection);
pp.query("use_thermal_diffusion", use_thermal_diffusion);
pp.query("alpha_T", alpha_T);
pp.query("use_scalar_diffusion", use_scalar_diffusion);
pp.query("alpha_S", alpha_S);
pp.query("use_momentum_diffusion", use_momentum_diffusion);
pp.query("dynamicViscosity", dynamicViscosity);
pp.query("use_smagorinsky", use_smagorinsky);
pp.query("Cs", Cs);
pp.query("use_pressure", use_pressure);
pp.query("use_gravity", use_gravity);
pp.query("spatial_order", spatial_order);
pp.query("abl_driver_type", abl_driver_type);
pp.queryarr("abl_pressure_grad", abl_pressure_grad);
pp.queryarr("abl_geo_wind", abl_geo_wind);
