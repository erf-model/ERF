#include "ERF.H"
#include "ERF_ShocInterface.H"

SHOCInterface::SHOCInterface (const int& lev,
                              SolverChoice& sc)
{
}

void
ERF::compute_shoc_tendencies (int lev,
                              amrex::MultiFab& cons_in,
                              amrex::MultiFab& xvel_in,
                              amrex::MultiFab& yvel_in,
                              amrex::MultiFab& zvel_in,
                              amrex::MultiFab& source,    amrex::MultiFab& xmom_src,
                              amrex::MultiFab& ymom_src,  amrex::MultiFab& zmom_src,
                              const amrex::Real& dt_advance)
{
#ifdef ERF_USE_NETCDF
    amrex::MultiFab *lat_ptr = lat_m[lev].get();
    amrex::MultiFab *lon_ptr = lon_m[lev].get();
#else
    amrex::MultiFab *lat_ptr = nullptr;
    amrex::MultiFab *lon_ptr = nullptr;
#endif

    amrex::Print() << "Advancing SHOC at level: " << lev << " ...";

    shoc_interface[lev]->set_grids(lev, cons_in.boxArray(), Geom(lev),
                                   &cons_in, z_phys_nd[lev].get(), lat_ptr, lon_ptr);

    shoc_interface[lev]->initialize_impl();
    shoc_interface[lev]->run_impl(dt_advance);
    shoc_interface[lev]->finalize_impl();

    amrex::Print() << "Done advancing SHOC\n";
}

void
SHOCInterface::set_grids (int& level,
                          const amrex::BoxArray& ba,
                          amrex::Geometry& geom,
                          amrex::MultiFab* cons_in,
                          amrex::MultiFab* z_phys,
                          amrex::MultiFab* lat,
                          amrex::MultiFab* lon)
{
    // using namespace ekat::units;

#if 0
    m_grid = grids_manager->get_grid("physics");
    const auto& grid_name = m_grid->name();
#endif

    // Define the different field layouts that will be used for this process

#if 0
    // Layout for horiz_wind field
    FieldLayout vector3d_mid = m_grid->get_3d_vector_layout(true,2);
    add_field<Required>("omega",          scalar3d_mid, Pa/s, grid_name, ps);

    // Layout for 2D (1d horiz X 1d vertical) variable
    // Layout for surf_mom_flux
    FieldLayout vector2d = m_grid->get_2d_vector_layout(2);

    // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces
    FieldLayout scalar2d = m_grid->get_2d_scalar_layout();
    FieldLayout scalar3d_mid = m_grid->get_3d_scalar_layout(true);
    FieldLayout scalar3d_int = m_grid->get_3d_scalar_layout(false);

    // Define fields needed in SHOC.
    // Note: shoc_main is organized by a set of 5 structures, variables below are organized
    //       using the same approach to make it easier to follow.

    constexpr int ps = Spack::n;

    const auto nondim = Units::nondimensional();
    const auto m2 = pow(m,2);
    const auto s2 = pow(s,2);

    // These variables are needed by the interface, but not actually passed to shoc_main.
    add_field<Required>("surf_sens_flux", scalar2d    , W/m2, grid_name);
    add_field<Required>("surf_mom_flux",  vector2d    , N/m2, grid_name);

    add_field<Updated>("surf_evap",       scalar2d    , kg/(m2*s), grid_name);
    add_field<Updated> ("T_mid",          scalar3d_mid, K,         grid_name, ps);
    add_tracer<Updated>("qv", m_grid, kg/kg, ps);

    // If TMS is a process, add surface drag coefficient to required fields
    if (m_params.get<bool>("apply_tms", false)) {
      add_field<Required>("surf_drag_coeff_tms", scalar2d,  kg/(m2*s), grid_name);
  }

    // Input variables
    add_field<Required>("p_mid",          scalar3d_mid, Pa,    grid_name, ps);
    add_field<Required>("p_int",          scalar3d_int, Pa,    grid_name, ps);
    add_field<Required>("pseudo_density", scalar3d_mid, Pa,    grid_name, ps);
    add_field<Required>("phis",           scalar2d    , m2/s2, grid_name, ps);

    // Input/Output variables
    add_field<Updated>("horiz_winds",   vector3d_mid,   m/s,     grid_name, ps);
    add_field<Updated>("sgs_buoy_flux", scalar3d_mid, K*(m/s), grid_name, ps);
    add_field<Updated>("eddy_diff_mom", scalar3d_mid, m2/s,    grid_name, ps);
    add_field<Updated>("cldfrac_liq",   scalar3d_mid, nondim,  grid_name, ps);
    add_tracer<Updated>("tke", m_grid, m2/s2, ps);
    add_tracer<Updated>("qc",  m_grid, kg/kg, ps);

    // Output variables
    add_field<Computed>("pbl_height",    scalar2d    , m,            grid_name);
    add_field<Computed>("inv_qc_relvar", scalar3d_mid, pow(kg/kg,2), grid_name, ps);
    add_field<Computed>("eddy_diff_heat",   scalar3d_mid, m2/s,        grid_name, ps);
    add_field<Computed>("w_variance",       scalar3d_mid, m2/s2,       grid_name, ps);
    add_field<Computed>("cldfrac_liq_prev", scalar3d_mid, nondim,      grid_name, ps);
    add_field<Computed>("ustar",            scalar2d,     m/s,         grid_name, ps);
    add_field<Computed>("obklen",           scalar2d,     m,           grid_name, ps);

    // Extra SHOC output diagnostics
    if (m_params.get<bool>("extra_shoc_diags", false)) {

        // Diagnostic output - mid point grid
        add_field<Computed>("brunt", scalar3d_mid, pow(s,-1), grid_name, ps);
        add_field<Computed>("shoc_mix", scalar3d_mid, m, grid_name, ps);
        add_field<Computed>("isotropy", scalar3d_mid, s, grid_name, ps);
        add_field<Computed>("shoc_cond", scalar3d_mid, kg/kg/s, grid_name, ps);
        add_field<Computed>("shoc_evap", scalar3d_mid, kg/kg/s, grid_name, ps);

        // Diagnostic output - interface grid
        add_field<Computed>("wthl_sec", scalar3d_int, K*(m/s), grid_name, ps);
        add_field<Computed>("thl_sec", scalar3d_int, pow(K,2), grid_name, ps);
        add_field<Computed>("wqw_sec", scalar3d_int, (kg/kg)*(m/s), grid_name, ps);
        add_field<Computed>("qw_sec", scalar3d_int, pow(kg/kg,2), grid_name, ps);
        add_field<Computed>("uw_sec", scalar3d_int, pow(m/s,2), grid_name, ps);
        add_field<Computed>("vw_sec", scalar3d_int, pow(m/s,2), grid_name, ps);
        add_field<Computed>("w3", scalar3d_int, pow(m/s,3), grid_name, ps);

    } // Extra SHOC output diagnostics

    // Tracer group
    add_group<Updated>("turbulence_advected_tracers", grid_name, ps, MonolithicAlloc::Required);

    // Boundary flux fields for energy and mass conservation checks
    if (has_column_conservation_check()) {
      add_field<Computed>("vapor_flux", scalar2d, kg/(m2*s), grid_name);
      add_field<Computed>("water_flux", scalar2d, m/s,     grid_name);
      add_field<Computed>("ice_flux",   scalar2d, m/s,     grid_name);
      add_field<Computed>("heat_flux",  scalar2d, W/m2,    grid_name);
    }
#endif

    // Set data members that may change
    m_lev            = level;
    m_geom           = geom;
    // m_cons_in        = cons_in;
    // m_z_phys         = z_phys;
    // m_lat            = lat;
    // m_lon            = lon;

    // Reset vector of offsets for columnar data
    m_num_layers = geom.Domain().length(2);

    m_num_cols = 0;
    m_col_offsets.clear();
    m_col_offsets.resize(int(ba.size()));

#if 0
    for (amrex::MFIter mfi(*m_cons_in); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        int nx = vbx.length(0);
        int ny = vbx.length(1);
        m_col_offsets[mfi.index()] = m_ncol;
        m_ncol += nx * ny;
    }
#endif

    // Allocate the buffer arrays
    alloc_buffers();

    // Fill the KOKKOS Views from AMReX MFs
    mf_to_kokkos_buffers();
}

void
SHOCInterface::alloc_buffers ()
{}

void
SHOCInterface::dealloc_buffers ()
{}

void
SHOCInterface::mf_to_kokkos_buffers ()
{}

void
SHOCInterface::kokkos_buffers_to_mf ()
{}

void
SHOCInterface::run_impl (const amrex::Real dt)
{
  EKAT_REQUIRE_MSG (dt<=300,
      "Error! SHOC is intended to run with a timestep no longer than 5 minutes.\n"
      "       Please, reduce timestep (perhaps increasing subcycling iterations).\n");

#if 0
  const auto nlev_packs  = ekat::npack<Spack>(m_num_layers);
  const auto scan_policy    = ekat::ExeSpaceUtils<KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(m_num_cols, nlev_packs);
  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nlev_packs);

  // Preprocessing of SHOC inputs. Kernel contains a parallel_scan,
  // so a special TeamPolicy is required.
  Kokkos::parallel_for("shoc_preprocess",
                       scan_policy,
                       shoc_preprocess);
  Kokkos::fence();

  auto wtracer_sfc = shoc_preprocess.wtracer_sfc;
  Kokkos::deep_copy(wtracer_sfc, 0);

  if (m_params.get<bool>("apply_tms", false)) {
    apply_turbulent_mountain_stress();
  }

  if (m_params.get<bool>("check_flux_state_consistency", false)) {
    check_flux_state_consistency(dt);
  }

  // For now set the host timestep to the shoc timestep. This forces
  // number of SHOC timesteps (nadv) to be 1.
  // TODO: input parameter?
  hdtime = dt;
  m_nadv = std::max(static_cast<int>(round(hdtime/dt)),1);

  // Reset internal WSM variables.
  workspace_mgr.reset_internals();

  // Run shoc main
  SHF::shoc_main(m_num_cols, m_num_layers, m_num_layers+1, m_npbl, m_nadv, m_num_tracers, dt,
                 workspace_mgr,runtime_options,input,input_output,output,history_output
#ifdef SCREAM_SHOC_SMALL_KERNELS
                 , temporaries
#endif
                 );

  // Postprocessing of SHOC outputs
  Kokkos::parallel_for("shoc_postprocess",
                       default_policy,
                       shoc_postprocess);
  Kokkos::fence();

  // Extra SHOC output diagnostics
  if (m_params.get<bool>("extra_shoc_diags", false)) {

    const auto& shoc_mix = get_field_out("shoc_mix").get_view<Spack**>();
    Kokkos::deep_copy(shoc_mix,history_output.shoc_mix);

    const auto& brunt = get_field_out("brunt").get_view<Spack**>();
    Kokkos::deep_copy(brunt,history_output.brunt);

    const auto& w3 = get_field_out("w3").get_view<Spack**>();
    Kokkos::deep_copy(w3,history_output.w3);

    const auto& isotropy = get_field_out("isotropy").get_view<Spack**>();
    Kokkos::deep_copy(isotropy,history_output.isotropy);

    const auto& wthl_sec = get_field_out("wthl_sec").get_view<Spack**>();
    Kokkos::deep_copy(wthl_sec,history_output.wthl_sec);

    const auto& wqw_sec = get_field_out("wqw_sec").get_view<Spack**>();
    Kokkos::deep_copy(wqw_sec,history_output.wqw_sec);

    const auto& uw_sec = get_field_out("uw_sec").get_view<Spack**>();
    Kokkos::deep_copy(uw_sec,history_output.uw_sec);

    const auto& vw_sec = get_field_out("vw_sec").get_view<Spack**>();
    Kokkos::deep_copy(vw_sec,history_output.vw_sec);

    const auto& qw_sec = get_field_out("qw_sec").get_view<Spack**>();
    Kokkos::deep_copy(qw_sec,history_output.qw_sec);

    const auto& thl_sec = get_field_out("thl_sec").get_view<Spack**>();
    Kokkos::deep_copy(thl_sec,history_output.thl_sec);

  } // Extra SHOC output diagnostics
#endif
}

#if 0
int3d_k
get_subcolumn_mask (const int ncol,
                    const int nlay,
                    const int ngpt,
                    real2d_k& cldf,
                    const int overlap_option,
                    int1d_k& seeds)
{
    // Routine will return subcolumn mask with values of 0 indicating no cloud, 1 indicating cloud
    int3d_k subcolumn_mask("subcolumn_mask", ncol, nlay, ngpt);

    // Subcolumn generators are a means for producing a variable x(i,j,k), where
    //
    //     c(i,j,k) = 1 for x(i,j,k) >  1 - cldf(i,j)
    //     c(i,j,k) = 0 for x(i,j,k) <= 1 - cldf(i,j)
    //
    // I am going to call this "cldx" to be just slightly less ambiguous
    real3d_k cldx("cldx", ncol, nlay, ngpt);

    // Apply overlap assumption to set cldx
    if (overlap_option == 0) {  // Dummy mask, always cloudy
        Kokkos::deep_copy(cldx, 1);
    } else {  // Default case, maximum-random overlap
        // Maximum-random overlap:
        // Uses essentially the algorithm described in eq (14) in Raisanen et al. 2004,
        // https://rmets.onlinelibrary.wiley.com/doi/epdf/10.1256/qj.03.99. Also the same
        // algorithm used in RRTMG implementation of maximum-random overlap (see
        // https://github.com/AER-RC/RRTMG_SW/blob/master/src/mcica_subcol_gen_sw.f90)
        //
        // First, fill cldx with random numbers.
        Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, nlay, ngpt}),
                             KOKKOS_LAMBDA (int icol, int ilay, int igpt)
        {
            conv::Random rand(seeds(icol) + ilay*ngpt + igpt);
            cldx(icol,ilay,igpt) = rand.genFP<RealT>();
        });

        // Step down columns and apply algorithm from eq (14)
        Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, ngpt}),
                             KOKKOS_LAMBDA (int icol, int igpt)
        {
            for (int ilay = 1; ilay < nlay; ilay++) {
                // Check cldx in level above and see if it satisfies conditions to create a cloudy subcolumn
                if (cldx(icol,ilay-1,igpt) > 1.0 - cldf(icol,ilay-1)) {
                    // Cloudy subcolumn above, use same random number here so that clouds in these two adjacent
                    // layers are maximimally overlapped
                    cldx(icol,ilay,igpt) = cldx(icol,ilay-1,igpt);
                } else {
                    // Cloud-less above, use new random number so that clouds are distributed
                    // randomly in this layer. Need to scale new random number to range
                    // [0, 1.0 - cldf(ilay-1)] because we have artificially changed the distribution
                    // of random numbers in this layer with the above branch of the conditional,
                    // which would otherwise inflate cloud fraction in this layer.
                    cldx(icol,ilay,igpt) = cldx(icol,ilay  ,igpt) * (1.0 - cldf(icol,ilay-1));
                }
            }
        });
    }

    // Use cldx array to create subcolumn mask
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, nlay, ngpt}),
                             KOKKOS_LAMBDA (int icol, int ilay, int igpt)
    {
        if (cldx(icol,ilay,igpt) > 1.0 - cldf(icol,ilay)) {
            subcolumn_mask(icol,ilay,igpt) = 1;
        } else {
            subcolumn_mask(icol,ilay,igpt) = 0;
        }
    });
    return subcolumn_mask;
}

#endif

void
SHOCInterface::initialize_impl ()
{
#if 0
  // Gather runtime options
  runtime_options.lambda_low    = m_params.get<double>("lambda_low");
  runtime_options.lambda_high   = m_params.get<double>("lambda_high");
  runtime_options.lambda_slope  = m_params.get<double>("lambda_slope");
  runtime_options.lambda_thresh = m_params.get<double>("lambda_thresh");
  runtime_options.thl2tune      = m_params.get<double>("thl2tune");
  runtime_options.qw2tune       = m_params.get<double>("qw2tune");
  runtime_options.qwthl2tune    = m_params.get<double>("qwthl2tune");
  runtime_options.w2tune        = m_params.get<double>("w2tune");
  runtime_options.length_fac    = m_params.get<double>("length_fac");
  runtime_options.c_diag_3rd_mom = m_params.get<double>("c_diag_3rd_mom");
  runtime_options.Ckh           = m_params.get<double>("coeff_kh");
  runtime_options.Ckm           = m_params.get<double>("coeff_km");
  runtime_options.shoc_1p5tke   = m_params.get<bool>("shoc_1p5tke");
  runtime_options.extra_diags   = m_params.get<bool>("extra_shoc_diags");

  // Initialize all of the structures that are passed to shoc_main in run_impl.
  // Note: Some variables in the structures are not stored in the field manager.  For these
  //       variables a local view is constructed.

  const auto& T_mid               = get_field_out("T_mid").get_view<Spack**>();
  const auto& p_mid               = get_field_in("p_mid").get_view<const Spack**>();
  const auto& p_int               = get_field_in("p_int").get_view<const Spack**>();
  const auto& pseudo_density      = get_field_in("pseudo_density").get_view<const Spack**>();
  const auto& omega               = get_field_in("omega").get_view<const Spack**>();
  const auto& surf_sens_flux      = get_field_in("surf_sens_flux").get_view<const Real*>();
  const auto& surf_evap           = get_field_in("surf_evap").get_view<const Real*>();
  const auto& surf_mom_flux       = get_field_in("surf_mom_flux").get_view<const Real**>();
  const auto& qtracers            = get_group_out("turbulence_advected_tracers").m_monolithic_field->get_strided_view<Spack***>();
  const auto& qc                  = get_field_out("qc").get_view<Spack**>();
  const auto& qv                  = get_field_out("qv").get_view<Spack**>();
  const auto& tke                 = get_field_out("tke").get_view<Spack**>();
  const auto& cldfrac_liq         = get_field_out("cldfrac_liq").get_view<Spack**>();
  const auto& cldfrac_liq_prev    = get_field_out("cldfrac_liq_prev").get_view<Spack**>();
  const auto& sgs_buoy_flux       = get_field_out("sgs_buoy_flux").get_view<Spack**>();
  const auto& tk                  = get_field_out("eddy_diff_mom").get_view<Spack**>();
  const auto& inv_qc_relvar       = get_field_out("inv_qc_relvar").get_view<Spack**>();
  const auto& phis                = get_field_in("phis").get_view<const Real*>();

  // Alias local variables from temporary buffer
  auto z_mid       = m_buffer.z_mid;
  auto z_int       = m_buffer.z_int;
  auto wpthlp_sfc  = m_buffer.wpthlp_sfc;
  auto wprtp_sfc   = m_buffer.wprtp_sfc;
  auto upwp_sfc    = m_buffer.upwp_sfc;
  auto vpwp_sfc    = m_buffer.vpwp_sfc;
  auto rrho        = m_buffer.rrho;
  auto rrho_i      = m_buffer.rrho_i;
  auto thv         = m_buffer.thv;
  auto dz          = m_buffer.dz;
  auto zt_grid     = m_buffer.zt_grid;
  auto zi_grid     = m_buffer.zi_grid;
  auto wtracer_sfc = m_buffer.wtracer_sfc;
  auto wm_zt       = m_buffer.wm_zt;
  auto inv_exner   = m_buffer.inv_exner;
  auto thlm        = m_buffer.thlm;
  auto qw          = m_buffer.qw;
  auto dse         = m_buffer.dse;
  auto tke_copy    = m_buffer.tke_copy;
  auto qc_copy     = m_buffer.qc_copy;
  auto shoc_ql2    = m_buffer.shoc_ql2;

  // For now, set z_int(i,nlevs) = z_surf = 0
  const Real z_surf = 0.0;

  // Some SHOC variables should be initialized uniformly if an Initial run
  if (run_type==RunType::Initial){
    Kokkos::deep_copy(sgs_buoy_flux,0.0);
    Kokkos::deep_copy(tk,0.0);
    Kokkos::deep_copy(tke,0.0004);
    Kokkos::deep_copy(tke_copy,0.0004);
    Kokkos::deep_copy(cldfrac_liq,0.0);
  }

  shoc_preprocess.set_variables(m_num_cols,m_num_layers,z_surf,
                                T_mid,p_mid,p_int,pseudo_density,omega,phis,surf_sens_flux,surf_evap,
                                surf_mom_flux,qtracers,qv,qc,qc_copy,tke,tke_copy,z_mid,z_int,
                                dse,rrho,rrho_i,thv,dz,zt_grid,zi_grid,wpthlp_sfc,wprtp_sfc,upwp_sfc,vpwp_sfc,
                                wtracer_sfc,wm_zt,inv_exner,thlm,qw, cldfrac_liq, cldfrac_liq_prev);

  // Input Variables:
  input.zt_grid     = shoc_preprocess.zt_grid;
  input.zi_grid     = shoc_preprocess.zi_grid;
  input.pres        = p_mid;
  input.presi       = p_int;
  input.pdel        = pseudo_density;
  input.thv         = shoc_preprocess.thv;
  input.w_field     = shoc_preprocess.wm_zt;
  input.wthl_sfc    = shoc_preprocess.wpthlp_sfc;
  input.wqw_sfc     = shoc_preprocess.wprtp_sfc;
  input.uw_sfc      = shoc_preprocess.upwp_sfc;
  input.vw_sfc      = shoc_preprocess.vpwp_sfc;
  input.wtracer_sfc = shoc_preprocess.wtracer_sfc;
  input.inv_exner   = shoc_preprocess.inv_exner;
  input.phis        = phis;

  // Input/Output Variables
  input_output.host_dse     = shoc_preprocess.shoc_s;
  input_output.tke          = shoc_preprocess.tke_copy;
  input_output.thetal       = shoc_preprocess.thlm;
  input_output.qw           = shoc_preprocess.qw;
  input_output.horiz_wind   = get_field_out("horiz_winds").get_view<Spack***>();
  input_output.wthv_sec     = sgs_buoy_flux;
  input_output.qtracers     = shoc_preprocess.qtracers;
  input_output.tk           = tk;
  input_output.shoc_cldfrac = cldfrac_liq;
  input_output.shoc_ql      = qc_copy;

  // Output Variables
  output.pblh     = get_field_out("pbl_height").get_view<Real*>();
  output.shoc_ql2 = shoc_ql2;
  output.tkh      = get_field_out("eddy_diff_heat").get_view<Spack**>();
  output.ustar    = get_field_out("ustar").get_view<Real*>();
  output.obklen   = get_field_out("obklen").get_view<Real*>();

  // Output (diagnostic)
  history_output.shoc_mix  = m_buffer.shoc_mix;
  history_output.isotropy  = m_buffer.isotropy;
  if (m_params.get<bool>("extra_shoc_diags", false)) {
    history_output.shoc_cond = get_field_out("shoc_cond").get_view<Spack**>();
    history_output.shoc_evap = get_field_out("shoc_evap").get_view<Spack**>();
  } else {
    history_output.shoc_cond = m_buffer.unused;
    history_output.shoc_evap = m_buffer.unused;
  }
  history_output.w_sec     = get_field_out("w_variance").get_view<Spack**>();
  history_output.thl_sec   = m_buffer.thl_sec;
  history_output.qw_sec    = m_buffer.qw_sec;
  history_output.qwthl_sec = m_buffer.qwthl_sec;
  history_output.wthl_sec  = m_buffer.wthl_sec;
  history_output.wqw_sec   = m_buffer.wqw_sec;
  history_output.wtke_sec  = m_buffer.wtke_sec;
  history_output.uw_sec    = m_buffer.uw_sec;
  history_output.vw_sec    = m_buffer.vw_sec;
  history_output.w3        = m_buffer.w3;
  history_output.wqls_sec  = m_buffer.wqls_sec;
  history_output.brunt     = m_buffer.brunt;

#ifdef SCREAM_SHOC_SMALL_KERNELS
  temporaries.se_b = m_buffer.se_b;
  temporaries.ke_b = m_buffer.ke_b;
  temporaries.wv_b = m_buffer.wv_b;
  temporaries.wl_b = m_buffer.wl_b;
  temporaries.se_a = m_buffer.se_a;
  temporaries.ke_a = m_buffer.ke_a;
  temporaries.wv_a = m_buffer.wv_a;
  temporaries.wl_a = m_buffer.wl_a;
  temporaries.kbfs = m_buffer.kbfs;
  temporaries.ustar2 = m_buffer.ustar2;
  temporaries.wstar = m_buffer.wstar;

  temporaries.rho_zt = m_buffer.rho_zt;
  temporaries.shoc_qv = m_buffer.shoc_qv;
  temporaries.tabs = m_buffer.tabs;
  temporaries.dz_zt = m_buffer.dz_zt;
  temporaries.dz_zi = m_buffer.dz_zi;
#endif

  shoc_postprocess.set_variables(m_num_cols,m_num_layers,
                                 rrho,qv,qw,qc,qc_copy,tke,tke_copy,qtracers,shoc_ql2,
                                 cldfrac_liq,inv_qc_relvar,
                                 T_mid, dse, z_mid, phis);

  if (has_column_conservation_check()) {
    const auto& vapor_flux = get_field_out("vapor_flux").get_view<Real*>();
    const auto& water_flux = get_field_out("water_flux").get_view<Real*>();
    const auto& ice_flux   = get_field_out("ice_flux").get_view<Real*>();
    const auto& heat_flux  = get_field_out("heat_flux").get_view<Real*>();
    shoc_postprocess.set_mass_and_energy_fluxes (surf_evap, surf_sens_flux,
                                                 vapor_flux, water_flux,
                                                 ice_flux, heat_flux);
  }

  // Set field property checks for the fields in this process
  using Interval = FieldWithinIntervalCheck;
  using LowerBound = FieldLowerBoundCheck;
  add_postcondition_check<Interval>(get_field_out("T_mid"),m_grid,100.0,500.0,false);
  add_postcondition_check<Interval>(get_field_out("qc"),m_grid,0.0,0.1,false);
  add_postcondition_check<Interval>(get_field_out("horiz_winds"),m_grid,-400.0,400.0,false);
  add_postcondition_check<LowerBound>(get_field_out("pbl_height"),m_grid,0);
  add_postcondition_check<Interval>(get_field_out("cldfrac_liq"),m_grid,0.0,1.0,false);
  add_postcondition_check<LowerBound>(get_field_out("tke"),m_grid,0);
  // For qv, ensure it doesn't get negative, by allowing repair of any neg value.
  // TODO: use a repairable lb that clips only "small" negative values
  add_postcondition_check<Interval>(get_field_out("qv"),m_grid,0,0.2,true);

  // Setup WSM for internal local variables
  const auto nlev_packs  = ekat::npack<Spack>(m_num_layers);
  const auto nlevi_packs = ekat::npack<Spack>(m_num_layers+1);
  const int n_wind_slots = ekat::npack<Spack>(2)*Spack::n;
  const int n_trac_slots = ekat::npack<Spack>(m_num_tracers+3)*Spack::n;
  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nlev_packs);
  workspace_mgr.setup(m_buffer.wsm_data, nlevi_packs, 14+(n_wind_slots+n_trac_slots), default_policy);

  // Calculate pref_mid, and use that to calculate
  // maximum number of levels in pbl from surface
  const auto pref_mid = m_buffer.pref_mid;
  const auto s_pref_mid = ekat::scalarize(pref_mid);
  const auto hyam = m_grid->get_geometry_data("hyam").get_view<const Real*>();
  const auto hybm = m_grid->get_geometry_data("hybm").get_view<const Real*>();
  const auto ps0 = C::P0;
  const auto psref = ps0;
  Kokkos::parallel_for(Kokkos::RangePolicy<>(0, m_num_layers), KOKKOS_LAMBDA (const int lev) {
    s_pref_mid(lev) = ps0*hyam(lev) + psref*hybm(lev);
  });
  Kokkos::fence();

  const int ntop_shoc = 0;
  const int nbot_shoc = m_num_layers;
  m_npbl = SHF::shoc_init(nbot_shoc,ntop_shoc,pref_mid);

  // Compute cell length for input dx and dy.
  const auto ncols = m_num_cols;
  view_1d cell_length("cell_length", ncols);
  if (m_grid->has_geometry_data("dx_short")) {
    // In this case IOP is running with a planar geometry
    auto dx = m_grid->get_geometry_data("dx_short").get_view<const Real,Host>()();
    Kokkos::deep_copy(cell_length, dx*1000); // convert km -> m
  } else {
    const auto area = m_grid->get_geometry_data("area").get_view<const Real*>();
    const auto lat  = m_grid->get_geometry_data("lat").get_view<const Real*>();
    Kokkos::parallel_for(ncols, KOKKOS_LAMBDA (const int icol) {
      // For now, we are considering dy=dx. Here, we
      // will need to compute dx/dy instead of cell_length
      // if we have dy!=dx.
      cell_length(icol) = PF::calculate_dx_from_area(area(icol),lat(icol));;
    });
  }
  input.dx = cell_length;
  input.dy = cell_length;

    // Initialize Kokkos SHOC pool allocator
    // const size_t nvar = 300;
    // const size_t nbnd = std::max(k_dist_sw_k->get_nband(),k_dist_sw_k->get_nband());
    // const size_t ncol = gas_concs_k.ncol;
    // const size_t nlay = gas_concs_k.nlay;
    // auto my_size_ref = static_cast<unsigned long>(nvar * ncol * nlay * nbnd);

    // pool_t::init(my_size_ref);

    // We are now initialized!
    initialized = true;
#endif
}

void
SHOCInterface::finalize_impl ()
{
    // Do nothing (per SHOCMacrophysics::finalize_impl())

    // Fill the AMReX MFs from Kokkos Views
    kokkos_buffers_to_mf();

    // Deallocate the buffer arrays
    dealloc_buffers();
}

