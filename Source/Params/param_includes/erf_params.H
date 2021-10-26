
// This file is automatically created by parse_erf_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh

static int state_interp_order;
static int lin_limit_state_interp;
static int do_reflux;
static int do_avg_down;
static int load_balance_verbosity;
static int add_ext_src;
static int add_forcing_src;
static amrex::Real CI;
static amrex::Real PrT;
static amrex::Real fixed_dt;
static amrex::Real initial_dt;
static amrex::Real dt_cutoff;
static amrex::Real max_dt;
static amrex::Real cfl;
static amrex::Real init_shrink;
static amrex::Real change_max;
static int bndry_func_thread_safe;
static int sum_interval;
static amrex::Real sum_per;
static std::string job_name;

static bool        use_state_advection;
static bool        use_momentum_advection;
static bool        use_thermal_diffusion;
static amrex::Real alpha_T;
static bool        use_scalar_diffusion;
static amrex::Real alpha_S;
static bool        use_momentum_diffusion;
static amrex::Real dynamicViscosity;
static bool        use_smagorinsky;
static amrex::Real Cs;
static bool        use_pressure;
static bool        use_gravity;
static int         spatial_order;
static std::string abl_driver_type;
static std::vector<amrex::Real> abl_pressure_grad;
static std::vector<amrex::Real> abl_geo_wind;
