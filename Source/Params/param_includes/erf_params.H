
// This file is automatically created by parse_erf_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh

static int levmsk_interior;
static int levmsk_covered;
static int levmsk_notcovered;
static int levmsk_physbnd;
static int state_interp_order;
static int lin_limit_state_interp;
static int do_reflux;
static int do_avg_down;
static int load_balance_verbosity;
static amrex::Real difmag;
static amrex::Real small_dens;
static amrex::Real small_temp;
static amrex::Real small_pres;
static amrex::Real small_ener;
static int add_ext_src;
static int add_forcing_src;
static int source_term_predictor;
static int first_order_hydro;
static std::string xl_ext_bc_type;
static std::string xr_ext_bc_type;
static std::string yl_ext_bc_type;
static std::string yr_ext_bc_type;
static std::string zl_ext_bc_type;
static std::string zr_ext_bc_type;
static int do_les;
static int use_explicit_filter;
static amrex::Real Cs;
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
static int track_grid_losses;
static int sum_interval;
static amrex::Real sum_per;
static int hard_cfl_limit;
static std::string job_name;

static bool        use_state_advection;
static bool        use_momentum_advection;
static bool        use_thermal_diffusion;
static amrex::Real alpha_T;
static bool        use_scalar_diffusion;
static amrex::Real alpha_S;
static bool        use_momentum_diffusion;
static amrex::Real kinematicViscosity;
static bool        use_smagorinsky;
static bool        use_pressure;
static bool        use_gravity;
static int         spatial_order;
