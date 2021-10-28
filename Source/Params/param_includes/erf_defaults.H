
// This file is automatically created by parse_erf_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh

int         ERF::state_interp_order = 1;
int         ERF::lin_limit_state_interp = 0;
int         ERF::do_reflux = 1;
int         ERF::do_avg_down = 1;
int         ERF::load_balance_verbosity = 0;
int         ERF::add_ext_src = 0;
int         ERF::add_forcing_src = 0;
amrex::Real ERF::CI = 0.0;
amrex::Real ERF::PrT = 1.0;
amrex::Real ERF::dt_cutoff = 0.0;
amrex::Real ERF::max_dt = 1.e200;
int         ERF::bndry_func_thread_safe = 1;
std::string ERF::job_name = "";

bool        ERF::use_state_advection = true;
bool        ERF::use_momentum_advection = true;
bool        ERF::use_thermal_diffusion = true;
amrex::Real ERF::alpha_T = 1.0;
bool        ERF::use_scalar_diffusion = true;
amrex::Real ERF::alpha_S = 1.0;
bool        ERF::use_momentum_diffusion = true;
amrex::Real ERF::dynamicViscosity = 1.7e-5;
bool        ERF::use_smagorinsky = true;
amrex::Real ERF::Cs = 0.1;
bool        ERF::use_pressure = true;
bool        ERF::use_gravity = true;
int         ERF::spatial_order = 2;
std::string ERF::abl_driver_type = "None"; //"None", "PressureGradient", or "GeostrophicWind"
std::vector<amrex::Real> ERF::abl_pressure_grad = {0.0, 0.0, 0.0}; // {dp/dx, dp/dy, dp/dz}  for "PressureGradient"
std::vector<amrex::Real> ERF::abl_geo_wind = {10.0, 0.0, 0.0}; // {u_geo, v_geo, w_geo} for "GeostrophicWind"
