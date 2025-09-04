#include "ERF.H"
#include "ERF_ShocInterface.H"

using namespace amrex;

SHOCInterface::SHOCInterface (const int& lev,
                              SolverChoice& sc)
{
    // Construct parser object for following reads
    ParmParse pp("erf.shoc");

    // Parse runtime inputs at start up
    pp.get("lambda_low"      , runtime_options.lambda_low    );
    pp.get("lambda_high"     , runtime_options.lambda_high   );
    pp.get("lambda_slope"    , runtime_options.lambda_slope  );
    pp.get("lambda_thresh"   , runtime_options.lambda_thresh );
    pp.get("thl2tune"        , runtime_options.thl2tune      );
    pp.get("qw2tune"         , runtime_options.qw2tune       );
    pp.get("qwthl2tune"      , runtime_options.qwthl2tune    );
    pp.get("w2tune"          , runtime_options.w2tune        );
    pp.get("length_fac"      , runtime_options.length_fac    );
    pp.get("c_diag_3rd_mom"  , runtime_options.c_diag_3rd_mom);
    pp.get("coeff_kh"        , runtime_options.Ckh           );
    pp.get("coeff_km"        , runtime_options.Ckm           );
    pp.get("shoc_1p5tke"     , runtime_options.shoc_1p5tke   );
    pp.get("extra_shoc_diags", runtime_options.extra_diags   );
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
    // Set data members that may change
    m_lev            = level;
    m_geom           = geom;
    m_cons_in        = cons_in;
    m_z_phys         = z_phys;

    // Ensure the boxes span klo -> khi
    int klo = geom.Domain().smallEnd(2);
    int khi = geom.Domain().bigEnd(2);

    // Reset vector of offsets for columnar data
    m_num_layers = geom.Domain().length(2);

    m_num_cols = 0;
    m_col_offsets.clear();
    m_col_offsets.resize(int(ba.size()));
    for (amrex::MFIter mfi(*m_cons_in); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE((klo == vbx.smallEnd(2)) &&
                                         (khi == vbx.bigEnd(2)),
                                         "Vertical decomposition with shoc is not allowed.");
        int nx = vbx.length(0);
        int ny = vbx.length(1);
        m_col_offsets[mfi.index()] = m_num_cols;
        m_num_cols += nx * ny;
    }

    // Allocate the buffer arrays
    alloc_buffers();

    // Fill the KOKKOS Views from AMReX MFs
    mf_to_kokkos_buffers();
}


void
SHOCInterface::alloc_buffers ()
{
    // SHOCPreprocess data structures
    //=======================================================
    phis             = view_1d("Phis"           , m_num_cols);
    wpthlp_sfc       = view_1d("Wpthlp_sfc"     , m_num_cols);
    wprtp_sfc        = view_1d("Wprtp_sfc"      , m_num_cols);
    upwp_sfc         = view_1d("Upwp_sfc"       , m_num_cols);
    vpwp_sfc         = view_1d("Vpwp_sfc"       , m_num_cols);
    surf_sens_flux   = view_1d("Sfc sens flux"  , m_num_cols);
    surf_evap        = view_1d("Sfc latent flux", m_num_cols);
    surf_mom_flux    = sview_2d("Sfc mom flux"  , m_num_cols, m_num_vel_comp);
    wtracer_sfc      = view_2d("Wtracer_sfc"    , m_num_cols, m_num_tracers );

    T_mid            = view_2d("T_mid"             , m_num_cols, m_num_layers  );
    p_mid            = view_2d("P_mid"             , m_num_cols, m_num_layers  );
    p_int            = view_2d("P_int"             , m_num_cols, m_num_layers+1);
    dens             = view_2d("Density"           , m_num_cols, m_num_layers  );
    omega            = view_2d("Omega"             , m_num_cols, m_num_layers  );
    qv               = view_2d("Qv"                , m_num_cols, m_num_layers  );
    qc               = view_2d("Qc"                , m_num_cols, m_num_layers  );
    qc_copy          = view_2d("Qc_copy"           , m_num_cols, m_num_layers  );
    z_mid            = view_2d("Z_mid"             , m_num_cols, m_num_layers  );
    z_int            = view_2d("Z_int"             , m_num_cols, m_num_layers+1);
    shoc_s           = view_2d("Shoc_s"            , m_num_cols, m_num_layers  );
    tke              = view_2d("Tke"               , m_num_cols, m_num_layers  );
    tke_copy         = view_2d("Tke_copy"          , m_num_cols, m_num_layers  );
    rrho             = view_2d("Rrho"              , m_num_cols, m_num_layers  );
    rrho_i           = view_2d("Rrho_i"            , m_num_cols, m_num_layers  );
    thv              = view_2d("Thv"               , m_num_cols, m_num_layers  );
    dz               = view_2d("dz"                , m_num_cols, m_num_layers  );
    dse              = view_2d("dse"               , m_num_cols, m_num_layers  );
    zt_grid          = view_2d("Zt_grid"           , m_num_cols, m_num_layers  );
    zi_grid          = view_2d("Zi_grid"           , m_num_cols, m_num_layers  );
    wm_zt            = view_2d("wm_zt"             , m_num_cols, m_num_layers  );
    inv_exner        = view_2d("Inv_exner"         , m_num_cols, m_num_layers  );
    thlm             = view_2d("Thlm"              , m_num_cols, m_num_layers  );
    qw               = view_2d("Qw"                , m_num_cols, m_num_layers  );
    cloud_frac       = view_2d("Cld_frac"          , m_num_cols, m_num_layers  );
    cldfrac_liq      = view_2d("Cld_frac_liq"      , m_num_cols, m_num_layers  );
    cldfrac_liq_prev = view_2d("cld_frac_liq_prev" , m_num_cols, m_num_layers  );

    // NOTE: Use layoutright format
    Kokkos::LayoutStride layout(m_num_cols   , m_num_cols*m_num_layers,   // stride for dim0
                                m_num_layers , m_num_layers,              // stride for dim1
                                m_num_tracers, 1                          // stride for dim2
                                );
    qtracers         = view_3d_strided("Qtracers"  , layout);


    // SHOCPostprocess data structures
    //=======================================================
    shoc_ql2      = view_2d("qc^2"         , m_num_cols, m_num_layers);
    inv_qc_relvar = view_2d("inv_qc_relvar", m_num_cols, m_num_layers);

    // SHOC InputOutput data structures
    //=======================================================
    sgs_buoy_flux = view_2d("sgs_buoy_flux", m_num_cols, m_num_layers);
    tk            = view_2d("eddy_diff_mom", m_num_cols, m_num_layers);

    horiz_wind    = view_3d("horiz_wind"   , m_num_cols, m_num_layers, m_num_vel_comp);

    // SHOC Output data structures
    //=======================================================
    pblh   = view_1d("pbl_height", m_num_cols);
    ustar  = view_1d("ustar"     , m_num_cols);
    obklen = view_1d("obklen"    , m_num_cols);

    tkh    = view_2d("eddy_diff_heat", m_num_cols, m_num_layers);

    // SHOC HistoryOutput data structures
    //=======================================================
    shoc_mix  = view_2d("shoc_mix" , m_num_cols, m_num_layers);
    w_sec     = view_2d("w_sec"    , m_num_cols, m_num_layers);
    thl_sec   = view_2d("thl_sec"  , m_num_cols, m_num_layers);
    qw_sec    = view_2d("qw_sec"   , m_num_cols, m_num_layers);
    qwthl_sec = view_2d("qwthl_sec", m_num_cols, m_num_layers);
    wthl_sec  = view_2d("wthl_sec" , m_num_cols, m_num_layers);
    wqw_sec   = view_2d("wqw_sec"  , m_num_cols, m_num_layers);
    wtke_sec  = view_2d("wtke_sec" , m_num_cols, m_num_layers);
    uw_sec    = view_2d("uw_sec"   , m_num_cols, m_num_layers);
    vw_sec    = view_2d("vw_sec"   , m_num_cols, m_num_layers);
    w3        = view_2d("w3"       , m_num_cols, m_num_layers);
    wqls_sec  = view_2d("wqls_sec" , m_num_cols, m_num_layers);
    brunt     = view_2d("brunt"    , m_num_cols, m_num_layers);
    isotropy  = view_2d("isotropy" , m_num_cols, m_num_layers);
    shoc_cond = view_2d("shoc_cond", m_num_cols, m_num_layers);
    shoc_evap = view_2d("shoc_evap", m_num_cols, m_num_layers);

    // SHOC Miscellaneous data structures
    //=======================================================
}


void
SHOCInterface::dealloc_buffers ()
{
    // SHOCPreprocess data structures
    //=======================================================
    phis             = view_1d();
    wpthlp_sfc       = view_1d();
    wprtp_sfc        = view_1d();
    upwp_sfc         = view_1d();
    vpwp_sfc         = view_1d();
    surf_sens_flux   = view_1d();
    surf_evap        = view_1d();
    surf_mom_flux    = sview_2d();
    wtracer_sfc      = view_2d();

    T_mid            = view_2d();
    p_mid            = view_2d();
    p_int            = view_2d();
    dens             = view_2d();
    omega            = view_2d();
    qv               = view_2d();
    qc               = view_2d();
    qc_copy          = view_2d();
    z_mid            = view_2d();
    z_int            = view_2d();
    shoc_s           = view_2d();
    tke              = view_2d();
    tke_copy         = view_2d();
    rrho             = view_2d();
    rrho_i           = view_2d();
    thv              = view_2d();
    dz               = view_2d();
    dse              = view_2d();
    zt_grid          = view_2d();
    zi_grid          = view_2d();
    wm_zt            = view_2d();
    inv_exner        = view_2d();
    thlm             = view_2d();
    qw               = view_2d();
    cloud_frac       = view_2d();
    cldfrac_liq      = view_2d();
    cldfrac_liq_prev = view_2d();

    qtracers         = view_3d_strided();

    // SHOCPostprocess data structures
    //=======================================================
    shoc_ql2      = view_2d();
    inv_qc_relvar = view_2d();

    // SHOC InputOutput data structures
    //=======================================================
    sgs_buoy_flux = view_2d();
    tk            = view_2d();

    horiz_wind    = view_3d();

    // SHOC Output data structures
    //=======================================================
    pblh   = view_1d();
    ustar  = view_1d();
    obklen = view_1d();

    tkh    = view_2d();

    // SHOC HistoryOutput data structures
    //=======================================================
    shoc_mix  = view_2d();
    w_sec     = view_2d();
    thl_sec   = view_2d();
    qw_sec    = view_2d();
    qwthl_sec = view_2d();
    wthl_sec  = view_2d();
    wqw_sec   = view_2d();
    wtke_sec  = view_2d();
    uw_sec    = view_2d();
    vw_sec    = view_2d();
    w3        = view_2d();
    wqls_sec  = view_2d();
    brunt     = view_2d();
    isotropy  = view_2d();
    shoc_cond = view_2d();
    shoc_evap = view_2d();

    // SHOC Miscellaneous data structures
    //=======================================================
}


void
SHOCInterface::mf_to_kokkos_buffers ()
{
    // Expose for device (shoc preprocess)
    //=======================================================
    auto T_mid_d            = T_mid;
    auto p_mid_d            = p_mid;
    auto p_int_d            = p_int;
    auto dens_d             = dens;
    auto omega_d            = omega;
    auto phis_d             = phis;
    auto surf_sens_flux_d   = surf_sens_flux;
    auto surf_evap_d        = surf_evap;
    auto surf_mom_flux_d    = surf_mom_flux;
    auto qtracers_d         = qtracers;
    auto qv_d               = qv;
    auto qc_d               = qc;
    auto qc_copy_d          = qc_copy;
    auto z_mid_d            = z_mid;
    auto z_int_d            = z_int;
    auto shoc_s_d           = shoc_s;
    auto tke_d              = tke;
    auto tke_copy_d         = tke_copy;
    auto rrho_d             = rrho;
    auto rrho_i_d           = rrho_i;
    auto thv_d              = thv;
    auto dz_d               = dz;
    auto dse_d              = dse;
    auto zt_grid_d          = zt_grid;
    auto zi_grid_d          = zi_grid;
    auto wpthlp_sfc_d       = wpthlp_sfc;
    auto wprtp_sfc_d        = wprtp_sfc;
    auto upwp_sfc_d         = upwp_sfc;
    auto vpwp_sfc_d         = vpwp_sfc;
    auto wtracer_sfc_d      = wtracer_sfc;
    auto wm_zt_d            = wm_zt;
    auto inv_exner_d        = inv_exner;
    auto thlm_d             = thlm;
    auto qw_d               = qw;
    auto cloud_frac_d       = cloud_frac;
    auto cldfrac_liq_d      = cldfrac_liq;
    auto cldfrac_liq_prev_d = cldfrac_liq_prev;

    int  ncol  = m_num_cols;
    int  nlay  = m_num_layers;
    Real dz    = m_geom.CellSize(2);
    bool moist = (m_cons_in->nComp() > RhoQ1_comp);
    bool ice   = (m_cons_in->nComp() > RhoQ3_comp);
    auto ProbLoArr = m_geom.ProbLoArray();
    for (MFIter mfi(*m_cons_in); mfi.isValid(); ++mfi) {
        const auto& vbx  = mfi.validbox();
        const int nx     = vbx.length(0);
        const int imin   = vbx.smallEnd(0);
        const int jmin   = vbx.smallEnd(1);
        const int offset = m_col_offsets[mfi.index()];
        const Array4<const Real>& cons_arr = m_cons_in->const_array(mfi);
        const Array4<const Real>& z_arr    = (m_z_phys) ? m_z_phys->const_array(mfi) :
                                                          Array4<const Real>{};
        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // map [i,j,k] 0-based to [icol, ilay] 0-based
            const int icol   = (j-jmin)*nx + (i-imin) + offset;
            const int ilay   = k;

            // EOS input (at CC)
            Real r  = cons_arr(i,j,k,Rho_comp);
            Real rt = cons_arr(i,j,k,RhoTheta_comp);
            Real qv = (moist) ? cons_arr(i,j,k,RhoQ1_comp)/r : 0.0;
            Real qc = (moist) ? cons_arr(i,j,k,RhoQ2_comp)/r : 0.0;
            Real qi = (ice)   ? cons_arr(i,j,k,RhoQ3_comp)/r : 0.0;

            // EOS avg to z-face
            Real r_lo  = cons_arr(i,j,k-1,Rho_comp);
            Real rt_lo = cons_arr(i,j,k-1,RhoTheta_comp);
            Real qv_lo = (moist) ? cons_arr(i,j,k-1,RhoQ1_comp)/r_lo : 0.0;
            Real r_avg  = 0.5 * (r  + r_lo);
            Real rt_avg = 0.5 * (rt + rt_lo);
            Real qv_avg = 0.5 * (qv + qv_lo);

            // Fill view1d
            phis_d(icol)             = 0.0;
            surf_sens_flux_d(icol)   = 0.0;
            surf_evap_d(icol)        = 0.0;
            surf_mom_flux_d(icol,0)  = 0.0;
            surf_mom_flux_d(icol,1)  = 0.0;

            // Fill view2d at CC
            dens_d(icol,ilay)     = r;
            T_mid_d(icol,ilay)    = getTgivenRandRTh(r, rt, qv);
            p_mid_d(icol,ilay)    = getPgivenRTh(rt, qv);
            qv_d(icol,ilay)       = qv;
            qc_d(icol,ilay)       = qc;
            qc_copy_d(icol,ilay)  = qc;
            tke_d(icol,ilay)      = std::max(cons_arr(i,j,k,RhoKE_comp)/r, 0.0);
            tke_copy_d(icol,ilay) = std::max(cons_arr(i,j,k,RhoKE_comp)/r, 0.0);

            dz_d(icol,ilay)      = (z_arr) ? 0.25 * ( (z_arr(i  ,j  ,k+1) - z_arr(i  ,j  ,k))
                                                    + (z_arr(i+1,j  ,k+1) - z_arr(i+1,j  ,k))
                                                    + (z_arr(i  ,j+1,k+1) - z_arr(i  ,j+1,k))
                                                    + (z_arr(i+1,j+1,k+1) - z_arr(i+1,j+1,k)) ) : dz;
            z_mid_d(icol,ilay)   = (z_arr) ? 0.125 * ( z_arr(i  ,j  ,k+1) + z_arr(i  ,j  ,k)
                                                     + z_arr(i+1,j  ,k+1) + z_arr(i+1,j  ,k)
                                                     + z_arr(i  ,j+1,k+1) + z_arr(i  ,j+1,k)
                                                     + z_arr(i+1,j+1,k+1) + z_arr(i+1,j+1,k) ) :
                                             ProbLoArr[2] + (k + 0.5) * dz;

            // Fill view2d at w-face
            p_int_d(icol,ilay) = getPgivenRTh(rt_avg, qv_avg);
            z_int_d(icol,ilay) = (z_arr) ? 0.25 * ( z_arr(i  ,j  ,k) + z_arr(i+1,j  ,k)
                                                  + z_arr(i  ,j+1,k) + z_arr(i+1,j+1,k) ) :
                                           ProbLoArr[2] + (k) * dz;
            if (ilay==(nlay-1)) {
                Real r_hi  = cons_arr(i,j,k+1,Rho_comp);
                Real rt_hi = cons_arr(i,j,k+1,RhoTheta_comp);
                Real qv_hi = (moist) ? std::max(cons_arr(i,j,k+1,RhoQ1_comp)/r_hi,0.0) : 0.0;
                r_avg  = 0.5 * (r  + r_hi);
                rt_avg = 0.5 * (rt + rt_hi);
                qv_avg = 0.5 * (qv + qv_hi);
                p_int_d(icol,ilay+1) = getPgivenRTh(rt_avg, qv_avg);
                z_int_d(icol,ilay+1) = (z_arr) ? 0.25 * ( z_arr(i  ,j  ,k+1) + z_arr(i+1,j  ,k+1)
                                                        + z_arr(i  ,j+1,k+1) + z_arr(i+1,j+1,k+1) ) :
                                                 ProbLoArr[2] + (k+1) * dz;
            }

            // Questionable guesses
            zt_grid_d(icol,ilay)          = z_mid_d(icol,ilay); // Our thermo grid is the height grid?
            zi_grid_d(icol,ilay)          = z_int_d(icol,ilay); // Heights of interface grid?
            cloud_frac_d(icol,ilay)       = ((qc+qi)>0.0) ? 1. : 0.;
            cldfrac_liq_d(icol,ilay)      = (qc>0.0)      ? 1. : 0.;
            cldfrac_liq_prev_d(icol,ilay) = (qc>0.0)      ? 1. : 0.;

            // Don't know
            wpthlp_sfc_d(icol)     = 0.0;
            wprtp_sfc_d(icol)      = 0.0;
            upwp_sfc_d(icol)       = 0.0;
            vpwp_sfc_d(icol)       = 0.0;
            wtracer_sfc_d(icol,0)  = 0.0;

            omega_d(icol,ilay)     = 0.0;
            shoc_s_d(icol,ilay)    = 0.0;
            rrho_d(icol,ilay)      = 0.0;
            rrho_i_d(icol,ilay)    = 0.0;
            thv_d(icol,ilay)       = 0.0;
            wm_zt_d(icol,ilay)     = 0.0;
            inv_exner_d(icol,ilay) = 0.0;
            thlm_d(icol,ilay)      = 0.0;
            qw_d(icol,ilay)        = 0.0;
            dse(icol,ilay)         = 0.0;


            // Fill view3d
            qtracers_d(icol,ilay,0) = 0.0;
        });
    }
}


void
SHOCInterface::kokkos_buffers_to_mf ()
{}


void
SHOCInterface::initialize_impl ()
{
    // For now, set z_int(i,nlevs) = z_surf = 0
    const Real z_surf = 0.0;

    // Set preprocess variables
    shoc_preprocess.set_variables(m_num_cols, m_num_layers, z_surf,
                                  T_mid, p_mid, p_int, dens, omega, phis, surf_sens_flux, surf_evap,
                                  surf_mom_flux, qtracers, qv, qc, qc_copy, tke, tke_copy, z_mid, z_int,
                                  dse, rrho, rrho_i, thv, dz, zt_grid, zi_grid, wpthlp_sfc, wprtp_sfc, upwp_sfc,
                                  vpwp_sfc, wtracer_sfc, wm_zt, inv_exner, thlm, qw, cldfrac_liq, cldfrac_liq_prev);

    // Input Variables:
    input.zt_grid     = zt_grid;
    input.zi_grid     = zi_grid;
    input.pres        = p_mid;
    input.presi       = p_int;
    input.pdel        = dens;
    input.thv         = thv;
    input.w_field     = wm_zt;
    input.wthl_sfc    = wpthlp_sfc;
    input.wqw_sfc     = wprtp_sfc;
    input.uw_sfc      = upwp_sfc;
    input.vw_sfc      = vpwp_sfc;
    input.wtracer_sfc = wtracer_sfc;
    input.inv_exner   = inv_exner;
    input.phis        = phis;

    // Input/Output Variables
    input_output.host_dse     = shoc_s;
    input_output.tke          = tke_copy;
    input_output.thetal       = thlm;
    input_output.qw           = qw;
    input_output.horiz_wind   = horiz_wind;
    input_output.wthv_sec     = sgs_buoy_flux;
    input_output.qtracers     = qtracers;
    input_output.tk           = tk;
    input_output.shoc_cldfrac = cldfrac_liq;
    input_output.shoc_ql      = qc_copy;

    // Output Variables
    output.pblh     = pblh;
    output.shoc_ql2 = shoc_ql2;
    output.tkh      = tkh;
    output.ustar    = ustar;
    output.obklen   = obklen;

    // Output (diagnostic)
    history_output.shoc_mix  = shoc_mix;
    history_output.isotropy  = isotropy;
    history_output.shoc_cond = shoc_cond;
    history_output.shoc_evap = shoc_evap;
    history_output.w_sec     = w_sec;
    history_output.thl_sec   = thl_sec;
    history_output.qw_sec    = qw_sec;
    history_output.qwthl_sec = qwthl_sec;
    history_output.wthl_sec  = wthl_sec;
    history_output.wqw_sec   = wqw_sec;
    history_output.wtke_sec  = wtke_sec;
    history_output.uw_sec    = uw_sec;
    history_output.vw_sec    = vw_sec;
    history_output.w3        = w3;
    history_output.wqls_sec  = wqls_sec;
    history_output.brunt     = brunt;

    // Set postprocess variables
    shoc_postprocess.set_variables(m_num_cols, m_num_layers,
                                   rrho, qv, qw, qc, qc_copy, tke,
                                   tke_copy, qtracers, shoc_ql2,
                                   cldfrac_liq, inv_qc_relvar,
                                   T_mid, dse, z_mid, phis);

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


void
SHOCInterface::finalize_impl ()
{
    // Do nothing (per SHOCMacrophysics::finalize_impl())

    // Fill the AMReX MFs from Kokkos Views
    kokkos_buffers_to_mf();

    // Deallocate the buffer arrays
    dealloc_buffers();
}

