
// This file is automatically created by parse_erf_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh

int         ERF::levmsk_interior = 0;
int         ERF::levmsk_covered = 1;
int         ERF::levmsk_notcovered = 2;
int         ERF::levmsk_physbnd = 3;
int         ERF::state_interp_order = 1;
int         ERF::lin_limit_state_interp = 0;
int         ERF::do_reflux = 1;
int         ERF::do_avg_down = 1;
int         ERF::load_balance_verbosity = 0;
amrex::Real ERF::difmag = 0.1;
amrex::Real ERF::small_dens = 1.e-200;
amrex::Real ERF::small_temp = 1.e-200;
amrex::Real ERF::small_pres = 1.e-200;
amrex::Real ERF::small_ener = -1.e200;
int         ERF::add_ext_src = 0;
int         ERF::add_forcing_src = 0;
int         ERF::source_term_predictor = 0;
int         ERF::first_order_hydro = 0;
std::string ERF::xl_ext_bc_type = "";
std::string ERF::xr_ext_bc_type = "";
std::string ERF::yl_ext_bc_type = "";
std::string ERF::yr_ext_bc_type = "";
std::string ERF::zl_ext_bc_type = "";
std::string ERF::zr_ext_bc_type = "";
int         ERF::do_les = 0;
int         ERF::use_explicit_filter = 0;
amrex::Real ERF::Cs = 0.0;
amrex::Real ERF::CI = 0.0;
amrex::Real ERF::PrT = 1.0;
amrex::Real ERF::fixed_dt = -1.0;
amrex::Real ERF::initial_dt = -1.0;
amrex::Real ERF::dt_cutoff = 0.0;
amrex::Real ERF::max_dt = 1.e200;
amrex::Real ERF::cfl = 0.8;
amrex::Real ERF::init_shrink = 1.0;
amrex::Real ERF::change_max = 1.1;
int         ERF::bndry_func_thread_safe = 1;
int         ERF::track_grid_losses = 0;
int         ERF::sum_interval = -1;
amrex::Real ERF::sum_per = -1.0e0;
int         ERF::hard_cfl_limit = 1;
std::string ERF::job_name = "";

bool        ERF::use_state_advection = true;
bool        ERF::use_momentum_advection = true;
bool        ERF::use_thermal_diffusion = true;
amrex::Real ERF::alpha_T = 1.0;
bool        ERF::use_scalar_diffusion = true;
amrex::Real ERF::alpha_S = 1.0;
bool        ERF::use_momentum_diffusion = true;
amrex::Real ERF::kinematicViscosity = 1.5e-5;
bool        ERF::use_smagorinsky = true;
bool        ERF::use_pressure = true;
bool        ERF::use_gravity = true;
int         ERF::spatial_order = 2;
