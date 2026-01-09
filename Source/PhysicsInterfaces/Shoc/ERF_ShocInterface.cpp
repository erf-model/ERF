#include "ERF.H"
#include "ERF_ShocInterface.H"

using namespace amrex;

SHOCInterface::SHOCInterface (const int& lev,
                              SolverChoice& sc)
{
    //
    // Defaults set from E3SM/components/eamxx/cime_config/namelist_defaults_eamxx.xml
    //
    // Turn off SGS variability in SHOC, effectively reducing it to a 1.5 TKE closure?
    bool def_shoc_1p5tke = false;
    // Minimum value of stability correction
    Real def_lambda_low = 0.001;
    // Maximum value of stability correction
    Real def_lambda_high = 0.04;
    // Slope of change from lambda_low to lambda_high
    Real def_lambda_slope = 2.65;
    // stability threshold for which to apply more stability correction
    Real def_lambda_thresh = 0.02;
    // Temperature variance tuning factor
    Real def_thl2tune = 1.0;
    // Moisture variance tuning factor
    Real def_qw2tune = 1.0;
    // Temperature moisture covariance
    Real def_qwthl2tune = 1.0;
    // Vertical velocity variance
    Real def_w2tune = 1.0;
    // Length scale factor
    Real def_length_fac = 0.5;
    // Third moment vertical velocity damping factor
    Real def_c_diag_3rd_mom = 7.0;
    // Eddy diffusivity coefficient for heat
    Real def_coeff_kh = 0.1;
    // Eddy diffusivity coefficient for momentum
    Real def_coeff_km = 0.1;

    runtime_options.lambda_low    = def_lambda_low;
    runtime_options.lambda_high   = def_lambda_high;
    runtime_options.lambda_slope  = def_lambda_slope;
    runtime_options.lambda_thresh = def_lambda_thresh;

    runtime_options.thl2tune   = def_thl2tune;
    runtime_options.qwthl2tune = def_qwthl2tune;
    runtime_options.qw2tune    = def_qw2tune;
    runtime_options.w2tune     = def_w2tune;

    runtime_options.length_fac     = def_length_fac;
    runtime_options.c_diag_3rd_mom = def_c_diag_3rd_mom;
    runtime_options.Ckh            = def_coeff_kh;
    runtime_options.Ckm            = def_coeff_km;
    runtime_options.shoc_1p5tke    = def_shoc_1p5tke;
    runtime_options.extra_diags    = extra_shoc_diags;

    // Construct parser object for following reads
    ParmParse pp("erf.shoc");

    // Parse runtime inputs at start up
    pp.query("lambda_low"      , runtime_options.lambda_low    );
    pp.query("lambda_high"     , runtime_options.lambda_high   );
    pp.query("lambda_slope"    , runtime_options.lambda_slope  );
    pp.query("lambda_thresh"   , runtime_options.lambda_thresh );
    pp.query("thl2tune"        , runtime_options.thl2tune      );
    pp.query("qw2tune"         , runtime_options.qw2tune       );
    pp.query("qwthl2tune"      , runtime_options.qwthl2tune    );
    pp.query("w2tune"          , runtime_options.w2tune        );
    pp.query("length_fac"      , runtime_options.length_fac    );
    pp.query("c_diag_3rd_mom"  , runtime_options.c_diag_3rd_mom);
    pp.query("coeff_kh"        , runtime_options.Ckh           );
    pp.query("coeff_km"        , runtime_options.Ckm           );
    pp.query("shoc_1p5tke"     , runtime_options.shoc_1p5tke   );
    pp.query("extra_shoc_diags", runtime_options.extra_diags   );

    // Set to default but allow us to change it through the inputs file
    pp.query("apply_tms", apply_tms);
    pp.query("check_flux_state", check_flux_state);
    pp.query("extra_shoc_diags", extra_shoc_diags);
    pp.query("column_conservation_check", column_conservation_check);
}


void
ERF::compute_shoc_tendencies (int lev,
                              MultiFab* cons,
                              MultiFab* xvel,
                              MultiFab* yvel,
                              MultiFab* zvel,
                              Real* w_subsid,
                              MultiFab* tau13,
                              MultiFab* tau23,
                              MultiFab* hfx3,
                              MultiFab* qfx3,
                              MultiFab* eddyDiffs,
                              MultiFab* z_phys_nd,
                              const Real& dt_advance)
{
    Print() << "Advancing SHOC at level: " << lev << " ...";

    shoc_interface[lev]->set_grids(lev, cons->boxArray(), Geom(lev),
                                   cons , xvel , yvel, zvel, w_subsid,
                                   tau13, tau23, hfx3, qfx3,
                                   eddyDiffs, z_phys_nd);

    auto t0 = amrex::second();
    shoc_interface[lev]->initialize_impl();
    auto t1 = amrex::second();
    shoc_interface[lev]->run_impl(dt_advance);
    auto t2 = amrex::second();
    shoc_interface[lev]->finalize_impl(dt_advance);
    auto t3 = amrex::second();

    Print() << "Time in SHOC: initialize_impl() " << t1 - t0 << std::endl;
    Print() << "Time in SHOC: run_impl()        " << t2 - t1 << std::endl;
    Print() << "Time in SHOC: finalize_impl()   " << t3 - t2 << std::endl;
    Print() << "Done advancing SHOC\n";
}


void
SHOCInterface::set_grids (int& level,
                          const BoxArray& ba,
                          Geometry& geom,
                          MultiFab* cons,
                          MultiFab* xvel,
                          MultiFab* yvel,
                          MultiFab* zvel,
                          Real* w_subsid,
                          MultiFab* tau13,
                          MultiFab* tau23,
                          MultiFab* hfx3,
                          MultiFab* qfx3,
                          MultiFab* eddyDiffs,
                          MultiFab* z_phys)
{
    // Set data members that may change
    m_lev     = level;
    m_geom    = geom;
    m_cons    = cons;
    m_xvel    = xvel;
    m_yvel    = yvel;
    m_zvel    = zvel;
    m_w_subsid = w_subsid;
    m_tau13   = tau13;
    m_tau23   = tau23;
    m_hfx3    = hfx3;
    m_qfx3    = qfx3;
    m_mu      = eddyDiffs;
    m_z_phys  = z_phys;

    // Ensure the boxes span klo -> khi
    int klo = geom.Domain().smallEnd(2);
    int khi = geom.Domain().bigEnd(2);

    // Reset vector of offsets for columnar data
    m_num_layers = geom.Domain().length(2);

    int num_cols = 0;
    m_col_offsets.clear();
    m_col_offsets.resize(int(ba.size()));
    for (MFIter mfi(*m_cons); mfi.isValid(); ++mfi) {
        // NOTE: Get lateral ghost cells for CC <--> FC
        const auto& gbx = mfi.tilebox(IntVect(0,0,0),IntVect(1,1,0));
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE((klo == gbx.smallEnd(2)) &&
                                         (khi == gbx.bigEnd(2)),
                                         "Vertical decomposition with shoc is not allowed.");
        int nx = gbx.length(0);
        int ny = gbx.length(1);
        m_col_offsets[mfi.index()] = num_cols;
        num_cols += nx * ny;
    }

    // Resize the Kokkos variables that persist in memory
    if (num_cols != m_num_cols) {
        sgs_buoy_flux = view_2d();
        tk            = view_2d();
        sgs_buoy_flux = view_2d("sgs_buoy_flux", num_cols, m_num_layers);
        tk            = view_2d("eddy_diff_mom", num_cols, m_num_layers);
    }
    m_num_cols = num_cols;

    // Allocate the tendency MultiFabs
    c_tend.define(m_cons->boxArray(), m_cons->DistributionMap(), m_cons->nComp(), 0);
    u_tend.define(m_xvel->boxArray(), m_xvel->DistributionMap(), m_xvel->nComp(), 0);
    v_tend.define(m_yvel->boxArray(), m_yvel->DistributionMap(), m_yvel->nComp(), 0);

    // Allocate the buffer arrays in ERF
    alloc_buffers();

    // Allocate the m_buffer struct
    init_buffers();

    // Fill the KOKKOS Views from AMReX MFs
    mf_to_kokkos_buffers();
}


void
SHOCInterface::alloc_buffers ()
{
    // Interface data structures
    //=======================================================
    omega            = view_2d("Omega"             , m_num_cols, m_num_layers  );
    surf_sens_flux   = view_1d("Sfc sens flux"     , m_num_cols);
    surf_mom_flux    = sview_2d("Sfc mom flux"     , m_num_cols, m_num_vel_comp);
    surf_evap        = view_1d("Sfc evap"          , m_num_cols);
    T_mid            = view_2d("T_mid"             , m_num_cols, m_num_layers  );
    qv               = view_2d("Qv"                , m_num_cols, m_num_layers  );
    surf_drag_coeff_tms = view_1d("surf_drag_coeff", m_num_cols);

    // Input data structures
    //=======================================================
    p_mid            = view_2d("P_mid"             , m_num_cols, m_num_layers  );
    p_int            = view_2d("P_int"             , m_num_cols, m_num_layers+1);
    pseudo_dens      = view_2d("Pseudo density"    , m_num_cols, m_num_layers  );
    phis             = view_1d("Phis"              , m_num_cols);

    // Input/Output data structures
    //=======================================================
    horiz_wind       = view_3d("horiz_wind"        , m_num_cols, m_num_vel_comp, m_num_layers);
    cldfrac_liq      = view_2d("Cld_frac_liq"      , m_num_cols, m_num_layers  );
    tke              = view_2d("Tke"               , m_num_cols, m_num_layers  );
    qc               = view_2d("Qc"                , m_num_cols, m_num_layers  );

    // Output data structures
    //=======================================================
    pblh             = view_1d("pbl_height"        , m_num_cols);
    inv_qc_relvar    = view_2d("inv_qc_relvar"     , m_num_cols, m_num_layers);
    tkh              = view_2d("eddy_diff_heat"    , m_num_cols, m_num_layers);
    w_sec            = view_2d("w_sec"             , m_num_cols, m_num_layers);
    cldfrac_liq_prev = view_2d("cld_frac_liq_prev" , m_num_cols, m_num_layers  );
    ustar            = view_1d("ustar"             , m_num_cols);
    obklen           = view_1d("obklen"            , m_num_cols);

    // Extra diagnostic data structures
    //=======================================================
    if (extra_shoc_diags) {
        brunt            = view_2d("brunt"    , m_num_cols, m_num_layers);
        shoc_mix         = view_2d("shoc_mix" , m_num_cols, m_num_layers);
        isotropy         = view_2d("isotropy" , m_num_cols, m_num_layers);
        shoc_cond        = view_2d("shoc_cond", m_num_cols, m_num_layers);
        shoc_evap        = view_2d("shoc_evap", m_num_cols, m_num_layers);

        wthl_sec         = view_2d("wthl_sec" , m_num_cols, m_num_layers+1);
        thl_sec          = view_2d("thl_sec"  , m_num_cols, m_num_layers+1);
        wqw_sec          = view_2d("wqw_sec"  , m_num_cols, m_num_layers+1);
        qw_sec           = view_2d("qw_sec"   , m_num_cols, m_num_layers+1);
        uw_sec           = view_2d("uw_sec"   , m_num_cols, m_num_layers+1);
        vw_sec           = view_2d("vw_sec"   , m_num_cols, m_num_layers+1);
        w3               = view_2d("w3"       , m_num_cols, m_num_layers+1);
    }

    // Tracer data structures
    //=======================================================
    // NOTE: Use layoutright format
    Kokkos::LayoutStride layout(m_num_cols   , m_num_cols*m_num_layers,   // stride for dim0
                                m_num_layers , m_num_layers,              // stride for dim1
                                m_num_tracers, 1                          // stride for dim2
                                );
    qtracers             = view_3d_strided("Qtracers"  , layout);

    // Boundary flux data structures
    //=======================================================
    if (column_conservation_check) {
        vapor_flux       = view_1d("vapor_flux", m_num_cols);
        water_flux       = view_1d("water_flux", m_num_cols);
        ice_flux         = view_1d("ice_flux"  , m_num_cols);
        heat_flux        = view_1d("heat_flux" , m_num_cols);
    }
}


void
SHOCInterface::dealloc_buffers ()
{
    // Contiguous memory buffer view
    //=======================================================
    tot_buff_view    = view_1d();

    // Interface data structures
    //=======================================================
    omega            = view_2d();
    surf_sens_flux   = view_1d();
    surf_mom_flux    = sview_2d();
    surf_evap        = view_1d();
    T_mid            = view_2d();
    qv               = view_2d();
    surf_drag_coeff_tms = view_1d();

    // Input data structures
    //=======================================================
    p_mid            = view_2d();
    p_int            = view_2d();
    pseudo_dens      = view_2d();
    phis             = view_1d();

    // Input/Output data structures
    //=======================================================
    horiz_wind       = view_3d();
    cldfrac_liq      = view_2d();
    tke              = view_2d();
    qc               = view_2d();

    // Output data structures
    //=======================================================
    pblh             = view_1d();
    inv_qc_relvar    = view_2d();
    tkh              = view_2d();
    w_sec            = view_2d();
    cldfrac_liq_prev = view_2d();
    ustar            = view_1d();
    obklen           = view_1d();

    // Extra diagnostic data structures
    //=======================================================
    if (extra_shoc_diags) {
        brunt            = view_2d();
        shoc_mix         = view_2d();
        isotropy         = view_2d();
        shoc_cond        = view_2d();
        shoc_evap        = view_2d();

        wthl_sec         = view_2d();
        thl_sec          = view_2d();
        wqw_sec          = view_2d();
        qw_sec           = view_2d();
        uw_sec           = view_2d();
        vw_sec           = view_2d();
        w3               = view_2d();
    }

    // Tracer data structures
    //=======================================================
    qtracers             = view_3d_strided();

    // Boundary flux data structures
    //=======================================================
    if (column_conservation_check) {
        vapor_flux       = view_1d();
        water_flux       = view_1d();
        ice_flux         = view_1d();
        heat_flux        = view_1d();
    }
}


void
SHOCInterface::mf_to_kokkos_buffers ()
{
    // FillBoundary for internal ghost cells for u/v averaging
    m_tau13->FillBoundary(m_geom.periodicity());
    m_tau23->FillBoundary(m_geom.periodicity());
    m_hfx3->FillBoundary(m_geom.periodicity());
    m_qfx3->FillBoundary(m_geom.periodicity());

    //
    // Expose for device capture
    //

    // Interface data structures
    //=======================================================
    auto omega_d = omega;
    auto surf_sens_flux_d = surf_sens_flux;
    auto surf_mom_flux_d = surf_mom_flux;
    auto surf_evap_d = surf_evap;
    auto T_mid_d = T_mid;
    auto qv_d = qv;
    auto surf_drag_coeff_tms_d = surf_drag_coeff_tms;

    // Input data structures
    //=======================================================
    auto p_mid_d = p_mid;
    auto p_int_d = p_int;
    auto pseudo_dens_d = pseudo_dens;
    auto phis_d = phis;

    // Input/Output data structures
    //=======================================================
    auto horiz_wind_d = horiz_wind;
    auto cldfrac_liq_d = cldfrac_liq;
    auto tke_d = tke;
    auto qc_d = qc;

    // Enforce the correct grid heights and density
    //=======================================================
    auto dz_d = m_buffer.dz;

    // Subsidence pointer to device vector data
    Real* w_sub = m_w_subsid;

    int  nlay  = m_num_layers;
    Real dz    = m_geom.CellSize(2);
    bool moist = (m_cons->nComp() > RhoQ1_comp);
    auto ProbLoArr = m_geom.ProbLoArray();

    auto domain    = m_geom.Domain();
    int ilo        = domain.smallEnd(0);
    int ihi        = domain.bigEnd(0);
    int jlo        = domain.smallEnd(1);
    int jhi        = domain.bigEnd(1);

    for (MFIter mfi(*m_cons); mfi.isValid(); ++mfi) {
        // NOTE: Grown box to get ghost cells in views
        const auto& gbx  = mfi.tilebox(IntVect(0,0,0),IntVect(1,1,0));
        const int nx     = gbx.length(0);
        const int imin   = gbx.smallEnd(0);
        const int jmin   = gbx.smallEnd(1);
        const int kmax   = gbx.bigEnd(2);
        const int offset = m_col_offsets[mfi.index()];

        const Array4<const Real>& cons_arr = m_cons->const_array(mfi);

        const Array4<const Real>& u_arr    = m_xvel->const_array(mfi);
        const Array4<const Real>& v_arr    = m_yvel->const_array(mfi);
        const Array4<const Real>& w_arr    = m_zvel->const_array(mfi);

        const Array4<const Real>& t13_arr   = m_tau13->const_array(mfi);
        const Array4<const Real>& t23_arr   = m_tau23->const_array(mfi);
        const Array4<const Real>& hfx3_arr  = m_hfx3->const_array(mfi);
        const Array4<const Real>& qfx3_arr  = m_qfx3->const_array(mfi);

        const Array4<const Real>& z_arr    = (m_z_phys) ? m_z_phys->const_array(mfi) :
                                                          Array4<const Real>{};
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // NOTE: k gets permuted with ilay
            // map [i,j,k] 0-based to [icol, ilay] 0-based
            const int icol   = (j-jmin)*nx + (i-imin) + offset;
            const int ilay   = kmax - k;
            const int ilayi  = kmax + 1 - k;

            // EOS input (at CC)
            Real r  = cons_arr(i,j,k,Rho_comp);
            Real rt = cons_arr(i,j,k,RhoTheta_comp);
            Real qv = (moist) ? cons_arr(i,j,k,RhoQ1_comp)/r : 0.0;
            Real qc = (moist) ? cons_arr(i,j,k,RhoQ2_comp)/r : 0.0;

            // EOS avg to z-face
            Real r_lo  = cons_arr(i,j,k-1,Rho_comp);
            Real rt_lo = cons_arr(i,j,k-1,RhoTheta_comp);
            Real qv_lo = (moist) ? cons_arr(i,j,k-1,RhoQ1_comp)/r_lo : 0.0;
            Real rt_avg = 0.5 * (rt + rt_lo);
            Real qv_avg = 0.5 * (qv + qv_lo);

            // Delta z
            Real delz = (z_arr) ? 0.25 * ( (z_arr(i  ,j  ,k+1) - z_arr(i  ,j  ,k))
                                         + (z_arr(i+1,j  ,k+1) - z_arr(i+1,j  ,k))
                                         + (z_arr(i  ,j+1,k+1) - z_arr(i  ,j+1,k))
                                         + (z_arr(i+1,j+1,k+1) - z_arr(i+1,j+1,k)) ) : dz;

            // W at cc (cannot be 0?; inspection of shoc code...)
            Real w_cc = 0.5 * (w_arr(i,j,k) + w_arr(i,j,k+1));
            w_cc += (w_sub) ? w_sub[k] : 0.0;
            Real w_limited = std::copysign(std::max(std::fabs(w_cc),1.0e-6),w_cc);

            // Input/Output data structures
            //=======================================================
            horiz_wind_d(icol,0,ilay)  = 0.5 * (u_arr(i,j,k) + u_arr(i+1,j  ,k));
            horiz_wind_d(icol,1,ilay)  = 0.5 * (v_arr(i,j,k) + v_arr(i  ,j+1,k));
            cldfrac_liq_d(icol,ilay)   = (qc>0.0) ? 1. : 0.;
            tke_d(icol,ilay)           = std::max(cons_arr(i,j,k,RhoKE_comp)/r, 0.0);
            qc_d(icol,ilay)            = qc;

            // Interface data structures
            //=======================================================
            // eamxx_common_physics_functions_impl.hpp: calculate_vertical_velocity
            omega_d(icol,ilay)           = -w_limited * r * CONST_GRAV;
            if (k==0) {
                int ii  = std::min(std::max(i,ilo),ihi);
                int jj  = std::min(std::max(j,jlo),jhi);

                surf_mom_flux_d(icol,0)  = 0.5 * (t13_arr(ii,jj,k) + t13_arr(ii+1,jj  ,k));
                surf_mom_flux_d(icol,1)  = 0.5 * (t23_arr(ii,jj,k) + t23_arr(ii  ,jj+1,k));
                // No unit conversion to W/m^2 (ERF_ShocInterface.H L224)
                surf_sens_flux_d(icol)   = hfx3_arr(ii,jj,k);
                surf_evap_d(icol)        = (moist) ? qfx3_arr(ii,jj,k) : 0.0;
                // Back out the drag coeff
                Real wsp = sqrt( horiz_wind_d(icol,0,ilay)[0]*horiz_wind_d(icol,0,ilay)[0]
                               + horiz_wind_d(icol,1,ilay)[0]*horiz_wind_d(icol,1,ilay)[0] );
                surf_drag_coeff_tms_d(icol) = surf_mom_flux_d(icol,0) /
                                              (-r * wsp * horiz_wind_d(icol,0,ilay)[0]);
            }
            T_mid_d(icol,ilay)          = getTgivenRandRTh(r, rt, qv);
            qv_d(icol,ilay)             = qv;

            // Input data structures
            //=======================================================
            p_mid_d(icol,ilay)       = getPgivenRTh(rt, qv);
            p_int_d(icol,ilayi)      = getPgivenRTh(rt_avg, qv_avg);
            // eamxx_common_physics_functions_impl.hpp: calculate_density
            pseudo_dens_d(icol,ilay) = r * CONST_GRAV * delz;
            // Enforce the grid spacing
            dz_d(icol,ilay) = delz;
            // Surface geopotential
            if (k==0) {
                Real z = (z_arr) ? 0.125 * ( (z_arr(i  ,j  ,k+1) + z_arr(i  ,j  ,k))
                                           + (z_arr(i+1,j  ,k+1) + z_arr(i+1,j  ,k))
                                           + (z_arr(i  ,j+1,k+1) + z_arr(i  ,j+1,k))
                                           + (z_arr(i+1,j+1,k+1) + z_arr(i+1,j+1,k)) ) :
                                   ProbLoArr[2];
                phis_d(icol) = CONST_GRAV * z;
            }

            if (ilay==0) {
                Real r_hi  = cons_arr(i,j,k+1,Rho_comp);
                Real rt_hi = cons_arr(i,j,k+1,RhoTheta_comp);
                Real qv_hi = (moist) ? std::max(cons_arr(i,j,k+1,RhoQ1_comp)/r_hi,0.0) : 0.0;
                rt_avg = 0.5 * (rt + rt_hi);
                qv_avg = 0.5 * (qv + qv_hi);
                p_int_d(icol,0) = getPgivenRTh(rt_avg, qv_avg);
            }
        });
    }
}


void
SHOCInterface::kokkos_buffers_to_mf (const Real dt)
{
    //
    // Expose for device capture
    //

    // Buffer data structures
    //=======================================================
    auto thlm_d = m_buffer.thlm;

    // Interface data structures
    //=======================================================
    auto T_mid_d = T_mid;
    auto qv_d = qv;

    // Input/Output data structures
    //=======================================================
    auto horiz_wind_d = horiz_wind;
    auto tke_d = tke;
    auto qc_d  = qc;

    bool moist = (m_cons->nComp() > RhoQ1_comp);
    for (MFIter mfi(*m_cons); mfi.isValid(); ++mfi) {
        // NOTE: No ghost cells when going back to MFs
        const auto& vbx_cc = mfi.validbox();
        const auto& vbx_x  = convert(vbx_cc,IntVect(1,0,0));
        const auto& vbx_y  = convert(vbx_cc,IntVect(0,1,0));

        // NOTE: Grown box only for mapping
        const auto& gbx  = mfi.tilebox(IntVect(0,0,0),IntVect(1,1,0));
        const int nx     = gbx.length(0);
        const int imin   = gbx.smallEnd(0);
        const int jmin   = gbx.smallEnd(1);
        const int kmax   = gbx.bigEnd(2);
        const int offset = m_col_offsets[mfi.index()];

        const Array4<const Real>& cons_arr = m_cons->const_array(mfi);
        const Array4<const Real>& u_arr    = m_xvel->const_array(mfi);
        const Array4<const Real>& v_arr    = m_yvel->const_array(mfi);

        const Array4<Real>& c_tend_arr     = c_tend.array(mfi);
        const Array4<Real>& u_tend_arr     = u_tend.array(mfi);
        const Array4<Real>& v_tend_arr     = v_tend.array(mfi);

        ParallelFor(vbx_cc, vbx_x, vbx_y,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // NOTE: k gets permuted with ilay
            // map [i,j,k] 0-based to [icol, ilay] 0-based
            const int icol   = (j-jmin)*nx + (i-imin) + offset;
            const int ilay   = kmax - k;

            // Density at CC
            Real r  = cons_arr(i,j,k,Rho_comp);

            // Theta at CC (eamxx_common_physics_functions_impl.hpp L123)
            Real Th = thlm_d(icol,ilay)[0] / ( 1.0 - (1.0 / T_mid_d(icol,ilay)[0])
                                             * (C::LatVap/C::Cpair) * qc_d(icol,ilay)[0] );

            // Populate the tendencies
            c_tend_arr(i,j,k,RhoTheta_comp)  = ( Th                  - cons_arr(i,j,k,RhoTheta_comp)/r ) / dt;
            c_tend_arr(i,j,k,RhoKE_comp)     = ( tke_d(icol,ilay)[0] - cons_arr(i,j,k,RhoKE_comp   )/r ) / dt;
            if (moist) {
                c_tend_arr(i,j,k,RhoQ1_comp) = ( qv_d(icol,ilay)[0] - cons_arr(i,j,k,RhoQ1_comp)/r ) / dt;
                c_tend_arr(i,j,k,RhoQ2_comp) = ( qc_d(icol,ilay)[0] - cons_arr(i,j,k,RhoQ2_comp)/r ) / dt;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // NOTE: k gets permuted with ilay
            // map [i,j,k] 0-based to [icol, ilay] 0-based
            const int icol   = (j-jmin)*nx + (i-imin) + offset;
            const int ilay   = kmax - k;

            int icolim = (j-jmin)*nx + (i-1-imin) + offset;
            Real uvel  = 0.5 * (horiz_wind_d(icol,0,ilay)[0] + horiz_wind_d(icolim,0,ilay)[0]);
            u_tend_arr(i,j,k) = ( uvel - u_arr(i,j,k) ) / dt;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // NOTE: k gets permuted with ilay
            // map [i,j,k] 0-based to [icol, ilay] 0-based
            const int icol   = (j-jmin)*nx + (i-imin) + offset;
            const int ilay   = kmax - k;

            int icoljm = (j-1-jmin)*nx + (i-imin) + offset;
            Real vvel  = 0.5 * (horiz_wind_d(icol,1,ilay)[0] + horiz_wind_d(icoljm,1,ilay)[0]);
            v_tend_arr(i,j,k) = ( vvel - v_arr(i,j,k) ) / dt;
        });
    }
}


void
SHOCInterface::set_eddy_diffs ()
{
    //
    // Expose for device capture
    //

    // Input/Output data structures
    //=======================================================
    auto tk_d = tk;

    // NOTE: Loop on grown box to fill ghost cells but limit
    //       to valid box where views are defined.
    for (MFIter mfi(*m_mu); mfi.isValid(); ++mfi) {
        const auto& gbx_cc = mfi.growntilebox();
        const auto& vbx_cc = mfi.validbox();

        // NOTE: Grown box only for mapping
        const auto& gbx  = mfi.tilebox(IntVect(0,0,0),IntVect(1,1,0));
        const int nx     = gbx.length(0);
        const int imin   = gbx.smallEnd(0);
        const int jmin   = gbx.smallEnd(1);
        const int kmax   = gbx.bigEnd(2);
        const int offset = m_col_offsets[mfi.index()];

        // Limiting to validbox
        const int iminv  = vbx_cc.smallEnd(0);
        const int imaxv  = vbx_cc.bigEnd(0);
        const int jminv  = vbx_cc.smallEnd(1);
        const int jmaxv  = vbx_cc.bigEnd(1);
        const int kminv  = vbx_cc.smallEnd(2);
        const int kmaxv  = vbx_cc.bigEnd(2);

        const Array4<Real>& mu_arr = m_mu->array(mfi);

        ParallelFor(gbx_cc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Limiting
            int ii = std::min(std::max(i,iminv),imaxv);
            int jj = std::min(std::max(j,jminv),jmaxv);
            int kk = std::min(std::max(k,kminv),kmaxv);

            // NOTE: k gets permuted with ilay
            // map [i,j,k] 0-based to [icol, ilay] 0-based
            const int icol   = (jj-jmin)*nx + (ii-imin) + offset;
            const int ilay   = kmax - kk;

            // NOTE: Set mom_v for tau_33, all other vertical comps are 0
            mu_arr(i,j,k,EddyDiff::Mom_v)   = tk_d(icol,ilay)[0];
            mu_arr(i,j,k,EddyDiff::Theta_v) = 0.0;
            mu_arr(i,j,k,EddyDiff::KE_v)    = 0.0;
            mu_arr(i,j,k,EddyDiff::Q_v)     = 0.0;
        });
    }

    // Correct the internal ghost cells that have foextrap
    m_mu->FillBoundary(m_geom.periodicity());
}


void
SHOCInterface::set_diff_stresses ()
{
    for (MFIter mfi(*m_hfx3); mfi.isValid(); ++mfi) {
        const auto& vbx_cc  = mfi.validbox();
        const auto& vbx_xz  = convert(vbx_cc,IntVect(1,0,1));
        const auto& vbx_yz  = convert(vbx_cc,IntVect(0,1,1));

        const Array4<Real>& hfx_arr = m_hfx3->array(mfi);
        const Array4<Real>& qfx_arr = m_qfx3->array(mfi);

        const Array4<Real>& t13_arr = m_tau13->array(mfi);
        const Array4<Real>& t23_arr = m_tau23->array(mfi);

        ParallelFor(vbx_cc, vbx_xz, vbx_yz,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            hfx_arr(i,j,k) = 0.;
            qfx_arr(i,j,k) = 0.;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            t13_arr(i,j,k) = 0.;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            t23_arr(i,j,k) = 0.;
        });
    }
}


void
SHOCInterface::add_fast_tend (Vector<MultiFab>& S_rhs)
{
    for (MFIter mfi(*m_cons); mfi.isValid(); ++mfi) {
        const auto& vbx_cc  = mfi.validbox();
        const auto& vbx_x   = convert(vbx_cc,IntVect(1,0,0));
        const auto& vbx_y   = convert(vbx_cc,IntVect(0,1,0));

        const Array4<const Real>& c_arr = m_cons->const_array(mfi);

        const Array4<Real>& cc_rhs_arr = S_rhs[IntVars::cons].array(mfi);
        const Array4<Real>& ru_rhs_arr = S_rhs[IntVars::xmom].array(mfi);
        const Array4<Real>& rv_rhs_arr = S_rhs[IntVars::ymom].array(mfi);

        const Array4<const Real>& c_tend_arr = c_tend.const_array(mfi);
        const Array4<const Real>& u_tend_arr = u_tend.const_array(mfi);
        const Array4<const Real>& v_tend_arr = v_tend.const_array(mfi);

        ParallelFor(vbx_cc, vbx_x, vbx_y,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cc_rhs_arr(i,j,k,RhoTheta_comp) += c_arr(i,j,k,Rho_comp) * c_tend_arr(i,j,k,RhoTheta_comp);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            ru_rhs_arr(i,j,k) += c_arr(i,j,k,Rho_comp) * u_tend_arr(i,j,k);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            rv_rhs_arr(i,j,k) += c_arr(i,j,k,Rho_comp) * v_tend_arr(i,j,k);
        });
    }
}


void
SHOCInterface::add_slow_tend (const MFIter& mfi,
                              const Box& tbx,
                              const Array4<Real>& cc_rhs_arr)
{
    bool moist = (m_cons->nComp() > RhoQ1_comp);

    const Array4<const Real>& c_arr = m_cons->const_array(mfi);

    const Array4<const Real>& c_tend_arr = c_tend.const_array(mfi);

    ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        cc_rhs_arr(i,j,k,RhoKE_comp) += c_arr(i,j,k,Rho_comp) * c_tend_arr(i,j,k,RhoKE_comp);
        if (moist) {
            cc_rhs_arr(i,j,k,RhoQ1_comp) += c_arr(i,j,k,Rho_comp) * c_tend_arr(i,j,k,RhoQ1_comp);
            cc_rhs_arr(i,j,k,RhoQ2_comp) += c_arr(i,j,k,Rho_comp) * c_tend_arr(i,j,k,RhoQ2_comp);
        }
    });
}


size_t
SHOCInterface::requested_buffer_size_in_bytes() const
{
    using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;

    const int nlev_packs       = ekat::npack<Spack>(m_num_layers);
    const int nlevi_packs      = ekat::npack<Spack>(m_num_layers+1);
    const int num_tracer_packs = ekat::npack<Spack>(m_num_tracers);

    // Number of Reals needed by local views in the interface
    const size_t interface_request = Buffer::num_1d_scalar_ncol*m_num_cols*sizeof(Real) +
                                     Buffer::num_1d_scalar_nlev*nlev_packs*sizeof(Spack) +
                                     Buffer::num_2d_vector_mid*m_num_cols*nlev_packs*sizeof(Spack) +
                                     Buffer::num_2d_vector_int*m_num_cols*nlevi_packs*sizeof(Spack) +
                                     Buffer::num_2d_vector_tr*m_num_cols*num_tracer_packs*sizeof(Spack);

    // Number of Reals needed by the WorkspaceManager passed to shoc_main
    const auto policy        = TPF::get_default_team_policy(m_num_cols, nlev_packs);
    const int n_wind_slots   = ekat::npack<Spack>(m_num_vel_comp)*Spack::n;
    const int n_trac_slots   = ekat::npack<Spack>(m_num_tracers) *Spack::n;
    const size_t wsm_request = WSM::get_total_bytes_needed(nlevi_packs, 14+(n_wind_slots+n_trac_slots), policy);

    return ( (interface_request + wsm_request)/sizeof(Real) );
}


void
SHOCInterface::init_buffers()
{

    // Buffer of contiguous memory
    auto buffer_size = requested_buffer_size_in_bytes();
    tot_buff_view = view_1d("contiguous shoc_buffer",buffer_size);
    Real* mem = reinterpret_cast<Real*>(tot_buff_view.data());

    // 1d scalar views
    using scalar_view_t = decltype(m_buffer.wpthlp_sfc);
    scalar_view_t* _1d_scalar_view_ptrs[Buffer::num_1d_scalar_ncol] =
        {&m_buffer.wpthlp_sfc, &m_buffer.wprtp_sfc, &m_buffer.upwp_sfc, &m_buffer.vpwp_sfc
#ifdef SCREAM_SHOC_SMALL_KERNELS
         , &m_buffer.se_b, &m_buffer.ke_b, &m_buffer.wv_b, &m_buffer.wl_b
         , &m_buffer.se_a, &m_buffer.ke_a, &m_buffer.wv_a, &m_buffer.wl_a
         , &m_buffer.kbfs, &m_buffer.ustar2, &m_buffer.wstar
#endif
        };
    for (int i = 0; i < Buffer::num_1d_scalar_ncol; ++i) {
        *_1d_scalar_view_ptrs[i] = scalar_view_t(mem, m_num_cols);
        mem += _1d_scalar_view_ptrs[i]->size();
    }

    Spack* s_mem = reinterpret_cast<Spack*>(mem);

    // 2d packed views
    const int nlev_packs       = ekat::npack<Spack>(m_num_layers);
    const int nlevi_packs      = ekat::npack<Spack>(m_num_layers+1);
    const int num_tracer_packs = ekat::npack<Spack>(m_num_tracers);

    m_buffer.pref_mid = decltype(m_buffer.pref_mid)(s_mem, nlev_packs);
    s_mem += m_buffer.pref_mid.size();

    using spack_2d_view_t = decltype(m_buffer.z_mid);
    spack_2d_view_t* _2d_spack_mid_view_ptrs[Buffer::num_2d_vector_mid] =
        {&m_buffer.z_mid, &m_buffer.rrho, &m_buffer.thv, &m_buffer.dz, &m_buffer.zt_grid, &m_buffer.wm_zt, &m_buffer.unused,
         &m_buffer.inv_exner, &m_buffer.thlm, &m_buffer.qw, &m_buffer.dse, &m_buffer.tke_copy, &m_buffer.qc_copy,
         &m_buffer.shoc_ql2, &m_buffer.shoc_mix, &m_buffer.isotropy, &m_buffer.w_sec, &m_buffer.wqls_sec, &m_buffer.brunt
#ifdef SCREAM_SHOC_SMALL_KERNELS
         , &m_buffer.rho_zt, &m_buffer.shoc_qv, &m_buffer.tabs, &m_buffer.dz_zt
#endif
        };

    spack_2d_view_t* _2d_spack_int_view_ptrs[Buffer::num_2d_vector_int] =
        {&m_buffer.z_int, &m_buffer.rrho_i, &m_buffer.zi_grid, &m_buffer.thl_sec, &m_buffer.qw_sec,
         &m_buffer.qwthl_sec, &m_buffer.wthl_sec, &m_buffer.wqw_sec, &m_buffer.wtke_sec, &m_buffer.uw_sec,
         &m_buffer.vw_sec, &m_buffer.w3
#ifdef SCREAM_SHOC_SMALL_KERNELS
         , &m_buffer.dz_zi
#endif
        };

    for (int i = 0; i < Buffer::num_2d_vector_mid; ++i) {
        *_2d_spack_mid_view_ptrs[i] = spack_2d_view_t(s_mem, m_num_cols, nlev_packs);
        s_mem += _2d_spack_mid_view_ptrs[i]->size();
    }

    for (int i = 0; i < Buffer::num_2d_vector_int; ++i) {
        *_2d_spack_int_view_ptrs[i] = spack_2d_view_t(s_mem, m_num_cols, nlevi_packs);
        s_mem += _2d_spack_int_view_ptrs[i]->size();
    }
    m_buffer.wtracer_sfc = decltype(m_buffer.wtracer_sfc)(s_mem, m_num_cols, num_tracer_packs);
    s_mem += m_buffer.wtracer_sfc.size();

    // WSM data
    m_buffer.wsm_data = s_mem;
}


void
SHOCInterface::initialize_impl ()
{
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

    // Set preprocess variables
    shoc_preprocess.set_variables(m_num_cols, m_num_layers, z_surf,
                                  T_mid, p_mid, p_int, pseudo_dens, omega, phis, surf_sens_flux, surf_evap,
                                  surf_mom_flux, qtracers, qv, qc, qc_copy, tke, tke_copy, z_mid, z_int,
                                  dse, rrho, rrho_i, thv, dz, zt_grid, zi_grid, wpthlp_sfc, wprtp_sfc, upwp_sfc,
                                  vpwp_sfc, wtracer_sfc, wm_zt, inv_exner, thlm, qw, cldfrac_liq, cldfrac_liq_prev);

    // Input Variables:
    input.zt_grid     = shoc_preprocess.zt_grid;
    input.zi_grid     = shoc_preprocess.zi_grid;
    input.pres        = p_mid;
    input.presi       = p_int;
    input.pdel        = pseudo_dens;
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
    input_output.horiz_wind   = horiz_wind;
    input_output.wthv_sec     = sgs_buoy_flux;
    input_output.qtracers     = shoc_preprocess.qtracers;
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
    history_output.shoc_mix  = m_buffer.shoc_mix;
    history_output.isotropy  = m_buffer.isotropy;
    if (extra_shoc_diags) {
        history_output.shoc_cond = shoc_cond;
        history_output.shoc_evap = shoc_evap;
    } else {
        history_output.shoc_cond = m_buffer.unused;
        history_output.shoc_evap = m_buffer.unused;
    }
    history_output.w_sec     = w_sec;
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

    // Set postprocess variables
    shoc_postprocess.set_variables(m_num_cols, m_num_layers,
                                   rrho, qv, qw, qc, qc_copy, tke,
                                   tke_copy, qtracers, shoc_ql2,
                                   cldfrac_liq, inv_qc_relvar,
                                   T_mid, dse, z_mid, phis);

    // Setup WSM for internal local variables
    using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;
    const auto nlev_packs  = ekat::npack<Spack>(m_num_layers);
    const auto nlevi_packs = ekat::npack<Spack>(m_num_layers+1);
    const int n_wind_slots = ekat::npack<Spack>(m_num_vel_comp)*Spack::n;
    const int n_trac_slots = ekat::npack<Spack>(m_num_tracers)*Spack::n;
    const auto default_policy = TPF::get_default_team_policy(m_num_cols, nlev_packs);
    workspace_mgr.setup(m_buffer.wsm_data, nlevi_packs, 14+(n_wind_slots+n_trac_slots), default_policy);

    // NOTE: Vertical indices were permuted, so top and bottom are correct
    // Maximum number of levels in pbl from surface
    const int ntop_shoc = 0;
    const int nbot_shoc = m_num_layers;
    auto p_mid_d = p_mid;
    view_1d pref_mid("pref_mid", m_num_layers);
    Spack* s_mem = reinterpret_cast<Spack*>(pref_mid.data());
    SHF::view_1d<Spack> pref_mid_um(s_mem, m_num_layers);
    const auto policy     =  TPF::get_default_team_policy(m_num_cols, nlev_packs);
    Kokkos::parallel_for("pref_mid",
                         policy,
                         KOKKOS_LAMBDA (const KT::MemberType& team)
    {
        const auto i = team.league_rank();
        if (i==0) {
            const auto& pmid_i = ekat::subview(p_mid_d, i);
            Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_packs), [&](const int& k)
            {
                pref_mid_um(k) = pmid_i(k);
            });
        }
    });
    Kokkos::fence();
    m_npbl = SHF::shoc_init(nbot_shoc,ntop_shoc,pref_mid_um);

    // Cell length for input dx and dy
    view_1d cell_length_x("cell_length_x", m_num_cols);
    view_1d cell_length_y("cell_length_y", m_num_cols);
    Kokkos::deep_copy(cell_length_x, m_geom.CellSize(0));
    Kokkos::deep_copy(cell_length_y, m_geom.CellSize(1));
    input.dx = cell_length_x;
    input.dy = cell_length_y;
}


void
SHOCInterface::run_impl (const Real dt)
{
    using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;

    EKAT_REQUIRE_MSG (dt<=300,
                      "Error! SHOC is intended to run with a timestep no longer than 5 minutes.\n"
                      "       Please, reduce timestep (perhaps increasing subcycling iterations).\n");

    const auto nlev_packs     = ekat::npack<Spack>(m_num_layers);
    const auto scan_policy    = TPF::get_thread_range_parallel_scan_team_policy(m_num_cols, nlev_packs);
    const auto default_policy = TPF::get_default_team_policy(m_num_cols, nlev_packs);

    // Preprocessing of SHOC inputs. Kernel contains a parallel_scan,
    // so a special TeamPolicy is required.
    Kokkos::parallel_for("shoc_preprocess",
                         scan_policy,
                         shoc_preprocess);
    Kokkos::fence();

    auto wtracer_sfc = shoc_preprocess.wtracer_sfc;
    Kokkos::deep_copy(wtracer_sfc, 0);

    if (apply_tms) {
        apply_turbulent_mountain_stress();
    }

    if (check_flux_state) {
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
    if (runtime_options.extra_diags) {
        Kokkos::deep_copy(shoc_mix,history_output.shoc_mix);
        Kokkos::deep_copy(brunt,history_output.brunt);
        Kokkos::deep_copy(w3,history_output.w3);
        Kokkos::deep_copy(isotropy,history_output.isotropy);
        Kokkos::deep_copy(wthl_sec,history_output.wthl_sec);
        Kokkos::deep_copy(wqw_sec,history_output.wqw_sec);
        Kokkos::deep_copy(uw_sec,history_output.uw_sec);
        Kokkos::deep_copy(vw_sec,history_output.vw_sec);
        Kokkos::deep_copy(qw_sec,history_output.qw_sec);
        Kokkos::deep_copy(thl_sec,history_output.thl_sec);
    }
}


void
SHOCInterface::finalize_impl (const Real dt)
{
    // Do nothing (per SHOCMacrophysics::finalize_impl())

    // Fill the AMReX MFs from Kokkos Views
    kokkos_buffers_to_mf(dt);

    // Deallocate the buffer arrays
    dealloc_buffers();
}


void
SHOCInterface::apply_turbulent_mountain_stress()
{
    auto rrho_i   = m_buffer.rrho_i;
    auto upwp_sfc = m_buffer.upwp_sfc;
    auto vpwp_sfc = m_buffer.vpwp_sfc;
    auto surf_drag_coeff_tms_d = surf_drag_coeff_tms;
    auto horiz_wind_d = horiz_wind;

    const int nlev_v  = (m_num_layers-1)/Spack::n;
    const int nlev_p  = (m_num_layers-1)%Spack::n;
    const int nlevi_v = m_num_layers/Spack::n;
    const int nlevi_p = m_num_layers%Spack::n;

    Kokkos::parallel_for("apply_tms", KT::RangePolicy(0, m_num_cols), KOKKOS_LAMBDA (const int i)
    {
        upwp_sfc(i) -= surf_drag_coeff_tms_d(i)*horiz_wind_d(i,0,nlev_v)[nlev_p]/rrho_i(i,nlevi_v)[nlevi_p];
        vpwp_sfc(i) -= surf_drag_coeff_tms_d(i)*horiz_wind_d(i,1,nlev_v)[nlev_p]/rrho_i(i,nlevi_v)[nlevi_p];
    });
}


void
SHOCInterface::check_flux_state_consistency(const double dt)
{
    using PC  = scream::physics::Constants<Real>;
    using RU  = ekat::ReductionUtils<KT::ExeSpace>;
    using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;

    const Real gravit = PC::gravit;
    const Real qmin   = 1e-12; // minimum permitted constituent concentration (kg/kg)

    const auto& pseudo_density = pseudo_dens;
    auto qv_d = qv;
    auto surf_evap_d = surf_evap;

    const auto nlevs           = m_num_layers;
    const auto nlev_packs      = ekat::npack<Spack>(nlevs);
    const auto last_pack_idx   = (nlevs-1)/Spack::n;
    const auto last_pack_entry = (nlevs-1)%Spack::n;
    const auto policy          = TPF::get_default_team_policy(m_num_cols, nlev_packs);
    Kokkos::parallel_for("check_flux_state_consistency",
                         policy,
                         KOKKOS_LAMBDA (const KT::MemberType& team)
    {
        const auto i = team.league_rank();

        const auto& pseudo_density_i = ekat::subview(pseudo_density, i);
        const auto& qv_i             = ekat::subview(qv_d, i);

        // reciprocal of pseudo_density at the bottom layer
        const auto rpdel = 1.0/pseudo_density_i(last_pack_idx)[last_pack_entry];

        // Check if the negative surface latent heat flux can exhaust
        // the moisture in the lowest model level. If so, apply fixer.
        const auto condition = surf_evap_d(i) - (qmin - qv_i(last_pack_idx)[last_pack_entry])/(dt*gravit*rpdel);
        if (condition < 0) {
            const auto cc = std::fabs(surf_evap_d(i)*dt*gravit);

            auto tracer_mass = [&](const int k)
            {
                return qv_i(k)*pseudo_density_i(k);
            };
            Real mm = RU::view_reduction(team, 0, nlevs, tracer_mass);

            EKAT_KERNEL_ASSERT_MSG(mm >= cc, "Error! Total mass of column vapor should be greater than mass of surf_evap.\n");

            Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_packs), [&](const int& k)
            {
                const auto adjust = cc*qv_i(k)*pseudo_density_i(k)/mm;
                qv_i(k) = (qv_i(k)*pseudo_density_i(k) - adjust)/pseudo_density_i(k);
            });

            surf_evap_d(i) = 0;
        }
    });
}
