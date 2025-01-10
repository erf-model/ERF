/*
 * RTE-RRTMGP radiation model interface to ERF
 * The original code is developed by RobertPincus, and the code is open source available at:
 *                        https://github.com/earth-system-radiation/rte-rrtmgp
 * Please reference to the following paper,
 *                        https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019MS001621
 * NOTE: we use the C++ version of RTE-RRTMGP, which is reimplemented the original Fortran
 * code using C++ YAKL for CUDA, HiP and SYCL application by E3SM ECP team, the C++ version
 * of the rte-rrtmgp code is located at:
 *                       https://github.com/E3SM-Project/rte-rrtmgp
 * The RTE-RRTMGP uses BSD-3-Clause Open Source License, if you want to make changes,
 * and modifications to the code, please refer to BSD-3-Clause Open Source License.
 */

#include <ERF_Radiation.H>

using namespace amrex;
using yakl::intrinsics::size;
using yakl::fortran::parallel_for;
using yakl::fortran::SimpleBounds;


Radiation::Radiation (SolverChoice& sc)
{
    // Initialize YAKL
    if (!yakl::isInitialized()) { yakl::init(); }

    // Check if we have a valid moisture model
    if (sc.moisture_type != MoistureType::None) { m_moist = true; }

    // Check if we have a moisture model with ice
    if (sc.moisture_type == MoistureType::SAM)  { m_ice = true; }

    // Construct parser object for following reads
    ParmParse pp("erf");

    // Radiation timestep, as a number of atm steps
    pp.query("rad_freq_in_steps", m_rad_freq_in_steps);

    // Do MCICA subcolumn sampling
    pp.query("do_subcol_sampling", m_do_subcol_sampling);

    // Determine orbital year. If orbital_year is negative, use current year
    // from timestamp for orbital year; if positive, use provided orbital year
    // for duration of simulation.
    m_fixed_orbital_year = pp.query("orbital_year",m_orbital_year);

    // Get orbital parameters from yaml file
    pp.query("orbital_eccentricity", m_orbital_eccen);
    pp.query("orbital_obliquity"   , m_orbital_obliq);
    pp.query("orbital_mvelp"       , m_orbital_mvelp);

    // Value for prescribing an invariant solar constant (i.e. total solar irradiance at
    // TOA).  Used for idealized experiments such as RCE. Disabled when value is less than 0.
    pp.query("fixed_total_solar_irradiance", m_fixed_total_solar_irradiance);

    // Determine whether or not we are using a fixed solar zenith angle (positive value)
    pp.query("Fixed Solar Zenith Angle", m_fixed_solar_zenith_angle);

    // Get prescribed surface values of greenhouse gases
    pp.query("o3vmr" , m_o3vmr );
    pp.query("co2vmr", m_co2vmr);
    pp.query("n2ovmr", m_n2ovmr);
    pp.query("ch4vmr", m_ch4vmr);
    pp.query("f11vmr", m_f11vmr);
    pp.query("f12vmr", m_f12vmr);
    pp.query("n2vmr" , m_n2vmr );
    pp.query("covmr" , m_covmr );

    // Required aerosol optical properties from SPA
    pp.query("do_aerosol_rad", m_do_aerosol_rad);

    // Whether we do extra clean/clear sky calculations
    pp.query("extra_clnclrsky_diag", m_extra_clnclrsky_diag);
    pp.query("extra_clnsky_diag"   , m_extra_clnsky_diag);

    // Parse path and file names
    pp.query("rrtmgp_file_path"      , rrtmgp_file_path);
    pp.query("rrtmgp_coeffs_sw"      , rrtmgp_coeffs_sw  );
    pp.query("rrtmgp_coeffs_lw"      , rrtmgp_coeffs_lw  );
    pp.query("rrtmgp_cloud_optics_sw", rrtmgp_cloud_optics_sw);
    pp.query("rrtmgp_cloud_optics_lw", rrtmgp_cloud_optics_lw);

    // Append file names to path
    rrtmgp_coeffs_file_sw       = rrtmgp_file_path + "/" + rrtmgp_coeffs_sw;
    rrtmgp_coeffs_file_lw       = rrtmgp_file_path + "/" + rrtmgp_coeffs_lw;
    rrtmgp_cloud_optics_file_sw = rrtmgp_file_path + "/" + rrtmgp_cloud_optics_sw;
    rrtmgp_cloud_optics_file_lw = rrtmgp_file_path + "/" + rrtmgp_cloud_optics_lw;
}

void
Radiation::set_grids (int& level,
                      int& step,
                      const amrex::Real& time,
                      const amrex::Real& dt,
                      const amrex::BoxArray& ba,
                      const amrex::Geometry& geom,
                      const amrex::MultiFab* cons_in,
                      amrex::MultiFab* lsm_fluxes,
                      amrex::MultiFab* lsm_zenith,
                      amrex::MultiFab* qheating_rates,
                      amrex::MultiFab* z_phys,
                      amrex::MultiFab* lat,
                      amrex::MultiFab* lon)

{
    // Reset data members for AMR
    m_lev            = level;
    m_step           = step;
    m_dt             = dt;
    m_geom           = geom;
    m_cons_in        = cons_in;
    m_lsm_fluxes     = lsm_fluxes;
    m_qheating_rates = qheating_rates;
    m_z_phys         = z_phys;
    m_lat            = lat;
    m_lon            = lon;

    // Update the day and month
    time_t timestamp = time_t(int(time));
    struct tm *timeinfo = localtime(&timestamp);
    if (m_fixed_orbital_year) {
        m_orbital_mon  = timeinfo->tm_mon + 1;
        m_orbital_day  = timeinfo->tm_mday;
    } else {
        m_orbital_year = timeinfo->tm_year + 1900;
        m_orbital_mon  = timeinfo->tm_mon  + 1;
        m_orbital_day  = timeinfo->tm_mday;
    }

    // Only allocate and proceed if we are going to update radiation
    m_update_rad = false;
    if (m_rad_freq_in_steps > 0) { m_update_rad = ( (m_step == 0) || (m_step % m_rad_freq_in_steps == 0) ); }

    if (update_rad) {
        // Reset vector of offsets for columnar data
        m_nlay = geom.Domain().length(2);

        m_ncol = 0;
        m_col_offsets.clear();
        m_col_offsets.resize(int(ba.size()));
        for (MFIter mfi(cons_in, TileNoZ()); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();
            int nx = vbx.length(0);
            int ny = vbx.length(1);
            m_col_offsets[mfi.index()] = ncol;
            m_ncol += nx * ny;
        }

        // Allocate the buffer arrays
        alloc_buffers();

        // Fill the YAKL Arrays from AMReX MFs
        mf_to_yakl_buffers();
    }
}

void
Radiation::alloc_buffers ()
{
    // 1d size (m_ngas)
    m_gas_mol_weights = real1d("m_gas_mol_weights", m_ngas);
    realHost1d m_gas_mol_weights_h("m_gas_mol_weights_h", m_ngas);
    parallel_for(m_ngas, YAKL_LAMBDA (int igas)
    {
        m_gas_mol_weights_h(igas)   = m_mol_weight_gas[igas];
        gas_names_yakl_offset(igas) = m_gas_names[igas];
    });
    m_gas_mol_weights.deep_copy_to(m_gas_mol_weights_h);


    // 1d size (ncol)
    cosine_zenith    = real1d("cosine_zenith"   , m_ncol);
    mu0              = real1d("mu0"             , m_ncol);
    sfc_alb_dir_vis  = real1d("sfc_alb_dir_vis" , m_ncol);
    sfc_alb_dir_nir  = real1d("sfc_alb_dir_nir" , m_ncol);
    sfc_alb_dif_vis  = real1d("sfc_alb_dif_vis" , m_ncol);
    sfc_alb_dif_nir  = real1d("sfc_alb_dif_nir" , m_ncol);
    sfc_flux_dir_vis = real1d("sfc_flux_dir_vis", m_ncol);
    sfc_flux_dir_nir = real1d("sfc_flux_dir_nir", m_ncol);
    sfc_flux_dif_vis = real1d("sfc_flux_dif_vis", m_ncol);
    sfc_flux_dif_nir = real1d("sfc_flux_dif_nir", m_ncol);

    // 2d size (ncol, nlay)
    d_dz          = real2d("d_dz"         , m_ncol, m_nlay);
    r_lay         = real2d("r_lay"        , m_ncol, m_nlay);
    p_lay         = real2d("p_lay"        , m_ncol, m_nlay);
    t_lay         = real2d("t_lay"        , m_ncol, m_nlay);
    z_del         = real2d("z_del"        , m_ncol, m_nlay);
    p_del         = real2d("p_del"        , m_ncol, m_nlay);
    qc            = real2d("qc"           , m_ncol, m_nlay);
    nc            = real2d("nc"           , m_ncol, m_nlay);
    qi            = real2d("qi"           , m_ncol, m_nlay);
    cldfrac_tot   = real2d("cldfrac_tot"  , m_ncol, m_nlay);
    eff_radius_qc = real2d("eff_radius_qc", m_ncol, m_nlay);
    eff_radius_qi = real2d("eff_radius_qi", m_ncol, m_nlay);
    tmp2d         = real2d("tmp2d"        , m_ncol, m_nlay);
    lwp           = real2d("lwp"          , m_ncol, m_nlay);
    iwp           = real2d("iwp"          , m_ncol, m_nlay);
    sw_heating    = real2d("sw_heating"   , m_ncol, m_nlay);
    lw_heating    = real2d("lw_heating"   , m_ncol, m_nlay);

    // 2d size (ncol, nlay+1)
    d_tint                   = real2d("d_tint"                  , m_ncol, m_nlay+1);
    p_lev                    = real2d("p_lev"                   , m_ncol, m_nlay+1);
    t_lev                    = real2d("t_lev"                   , m_ncol, m_nlay+1);
    sw_flux_up               = real2d("sw_flux_up"              , m_ncol, m_nlay+1);
    sw_flux_dn               = real2d("sw_flux_dn"              , m_ncol, m_nlay+1);
    sw_flux_dn_dir           = real2d("sw_flux_dn_dir"          , m_ncol, m_nlay+1);
    lw_flux_up               = real2d("sw_flux_up"              , m_ncol, m_nlay+1);
    lw_flux_dn               = real2d("sw_flux_dn"              , m_ncol, m_nlay+1);
    sw_clnclrsky_flux_up     = real2d("sw_clnclrsky_flux_up"    , m_ncol, m_nlay+1);
    sw_clnclrsky_flux_dn     = real2d("sw_clnclrsky_flux_dn"    , m_ncol, m_nlay+1);
    sw_clnclrsky_flux_dn_dir = real2d("sw_clnclrsky_flux_dn_dir", m_ncol, m_nlay+1);
    sw_clrsky_flux_up        = real2d("sw_clrsky_flux_up"       , m_ncol, m_nlay+1);
    sw_clrsky_flux_dn        = real2d("sw_clrsky_flux_dn"       , m_ncol, m_nlay+1);
    sw_clrsky_flux_dn_dir    = real2d("sw_clrsky_flux_dn_dir"   , m_ncol, m_nlay+1);
    sw_clnsky_flux_up        = real2d("sw_clnsky_flux_up"       , m_ncol, m_nlay+1);
    sw_clnsky_flux_dn        = real2d("sw_clnsky_flux_dn"       , m_ncol, m_nlay+1);
    sw_clnsky_flux_dn_dir    = real2d("sw_clnsky_flux_dn_dir"   , m_ncol, m_nlay+1);
    lw_clnclrsky_flux_up     = real2d("lw_clnclrsky_flux_up"    , m_ncol, m_nlay+1);
    lw_clnclrsky_flux_dn     = real2d("lw_clnclrsky_flux_dn"    , m_ncol, m_nlay+1);
    lw_clrsky_flux_up        = real2d("lw_clrsky_flux_up"       , m_ncol, m_nlay+1);
    lw_clrsky_flux_dn        = real2d("lw_clrsky_flux_dn"       , m_ncol, m_nlay+1);
    lw_clnsky_flux_up        = real2d("lw_clnsky_flux_up"       , m_ncol, m_nlay+1);
    lw_clnsky_flux_dn        = real2d("lw_clnsky_flux_dn"       , m_ncol, m_nlay+1);

    // 3d size (ncol, nlay+1, nswbands)
    sw_bnd_flux_up  = real3d("sw_bnd_flux_up" , m_ncol, m_nlay+1, m_nswbands);
    sw_bnd_flux_dn  = real3d("sw_bnd_flux_dn" , m_ncol, m_nlay+1, m_nswbands);
    sw_bnd_flux_dir = real3d("sw_bnd_flux_dir", m_ncol, m_nlay+1, m_nswbands);
    sw_bnd_flux_dif = real3d("sw_bnd_flux_dif", m_ncol, m_nlay+1, m_nswbands);

    // 2d size (ncol, nswbands)
    sfc_alb_dir = real2d("sfc_alb_dir", m_ncol, m_nswbands);
    sfc_alb_dif = real2d("sfc_alb_dif", m_ncol, m_nswbands);

    // 3d size (ncol, nlay, n[sw,lw]bands)
    aero_tau_sw = real3d("aero_tau_sw", m_ncol, m_nlay, m_nswbands);
    aero_ssa_sw = real3d("aero_ssa_sw", m_ncol, m_nlay, m_nswbands);
    aero_g_sw   = real3d("aero_g_sw"  , m_ncol, m_nlay, m_nswbands);
    aero_tau_lw = real3d("aero_tau_lw", m_ncol, m_nlay, m_nlwbands);

    // 3d size (ncol, nlay, n[sw,lw]bnds)
    cld_tau_sw_bnd = real3d("cld_tau_sw_bnd", m_ncol, m_nlay, m_nswbands);
    cld_tau_lw_bnd = real3d("cld_tau_lw_bnd", m_ncol, m_nlay, m_nlwbands);

    // 3d size (ncol, nlay, n[sw,lw]gpts)
    cld_tau_sw_gpt = real3d("cld_tau_sw_gpt", m_ncol, m_nlay, m_nswgpts);
    cld_tau_lw_gpt = real3d("cld_tau_lw_gpt", m_ncol, m_nlay, m_nlwgpts);
}

void
Radiation::dealloc_buffers ()
{
    // 1d size (m_ngas)
    m_gas_mol_weights.deallocate();

    // 1d size (ncol)
    cosine_zenith.deallocate();
    mu0.deallocate();
    sfc_alb_dir_vis.deallocate();
    sfc_alb_dir_nir.deallocate();
    sfc_alb_dif_vis.deallocate();
    sfc_alb_dif_nir.deallocate();
    sfc_flux_dir_vis.deallocate();
    sfc_flux_dir_nir.deallocate();
    sfc_flux_dif_vis.deallocate();
    sfc_flux_dif_nir.deallocate();

    // 2d size (ncol, nlay)
    d_dz.deallocate();
    r_lay.deallocate();
    p_lay.deallocate();
    t_lay.deallocate();
    z_del.deallocate();
    p_del.deallocate();
    qc.deallocate();
    nc.deallocate();
    qi.deallocate();
    cldfrac_tot.deallocate();
    eff_radius_qc.deallocate();
    eff_radius_qi.deallocate();
    tmp2d.deallocate();
    lwp.deallocate();
    iwp.deallocate();

    sw_heating.deallocate();
    lw_heating.deallocate();

    // 2d size (ncol, nlay+1)
    d_tint.deallocate();
    p_lev.deallocate();
    t_lev.deallocate();

    sw_flux_up.deallocate();
    sw_flux_dn.deallocate();
    sw_flux_dn_dir.deallocate();
    lw_flux_up.deallocate();
    lw_flux_dn.deallocate();
    sw_clnclrsky_flux_up.deallocate();
    sw_clnclrsky_flux_dn.deallocate();
    sw_clnclrsky_flux_dn_dir.deallocate();
    sw_clrsky_flux_up.deallocate();
    sw_clrsky_flux_dn.deallocate();
    sw_clrsky_flux_dn_dir.deallocate();
    sw_clnsky_flux_up.deallocate();
    sw_clnsky_flux_dn.deallocate();
    sw_clnsky_flux_dn_dir.deallocate();
    lw_clnclrsky_flux_up.deallocate();
    lw_clnclrsky_flux_dn.deallocate();
    lw_clrsky_flux_up.deallocate();
    lw_clrsky_flux_dn.deallocate();
    lw_clnsky_flux_up.deallocate();
    lw_clnsky_flux_dn.deallocate();

    // 3d size (ncol, nlay+1, nswbands)
    sw_bnd_flux_up.deallocate();
    sw_bnd_flux_dn.deallocate();
    sw_bnd_flux_dir.deallocate();
    sw_bnd_flux_dif.deallocate();

    // 2d size (ncol, nswbands)
    sfc_alb_dir.deallocate();
    sfc_alb_dif.deallocate();

    // 3d size (ncol, nlay, n[sw,lw]bands)
    aero_tau_sw.deallocate();
    aero_ssa_sw.deallocate();
    aero_g_sw.deallocate();
    aero_tau_lw.deallocate();

    // 3d size (ncol, nlay, n[sw,lw]bnds)
    cld_tau_sw_bnd.deallocate();
    cld_tau_lw_bnd.deallocate();

    // 3d size (ncol, nlay, n[sw,lw]gpts)
    cld_tau_sw_gpt.deallocate();
    cld_tau_lw_gpt.deallocate();
}


void
Radiation::mf_to_yakl_buffers ()
{
    bool moist = m_moist;
    bool ice   = m_ice;
    int  ncol  = m_ncol;
    int  nlay  = m_nlay;
    Real dz    = geom.CellSize(2);
    for (MFIter mfi(*m_cons_in); mfi.isValid(); ++mfi) {
        const auto& vbx  = mfi.validbox();
        const int nx     = vbx.length(0);
        const int imin   = vbx.smallEnd(0);
        const int jmin   = vbx.smallEnd(1);
        const int offset = m_col_offsets[mfi.Index()];
        const Array4<      Real>& cons_arr = m_cons_in->const_array(mfi);
        const Array4<const Real>& z_arr    = (m_z_phys) ? m_z_phys->const_array(mfi) :
                                                          Array4<const Real>{};
        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // map [i,j,k] 0-based to [icol, ilay] 1-based
            const int icol = (j-jmin)*nx + (i-imin) + 1 + offset;
            const int ilay = k+1;

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
            Real qc_lo = (moist) ? cons_arr(i,j,k-1,RhoQ2_comp)/r_lo : 0.0;
            Real qi_lo = (ice)   ? cons_arr(i,j,k-1,RhoQ3_comp)/r_lo : 0.0;
            Real r_avg  = 0.5 * (r  + r_lo);
            Real rt_avg = 0.5 * (rt + rt_lo);
            Real qv_avg = 0.5 * (qv + qv_lo);
            Real qc_avg = 0.5 * (qc + qc_lo);
            Real qi_avg = 0.5 * (qi + qi_lo);

            // Buffers at CC
            r_lay(icol,ilay) = r;
            p_lay(icol,ilay) = getPgivenRTh(rt, qv);
            t_lay(icol,ilay) = getTgivenRandRTh(r, rt, qv);
            z_del(icol,ilay) = (z_arr) ? 0.25 * ( (z_arr(i  ,j  ,k+1) - z_arr(i  ,j  ,k))
                                                + (z_arr(i+1,j  ,k+1) - z_arr(i+1,j  ,k))
                                                + (z_arr(i  ,j+1,k+1) - z_arr(i  ,j+1,k))
                                                + (z_arr(i+1,j  ,k+1) - z_arr(i+1,j  ,k)) ) : dz;
            qc(icol,ilay)    = qc;
            qi(icol,ilay)    = qi;
            cld_frac_tot(icol,ilay) = ((qc+qi)>1.0e-5) ? 1. : 0.;

            // HACK HACK HACK
            lwp(icol,ilay) = 0.0;
            iwp(icol,ilay) = 0.0;

            // Buffers on z-faces (nlay+1)
            p_lev(icol,ilay) = getPgivenRTh(rt_avg, qv_avg);
            t_lev(icol,ilay) = getTgivenRandRTh(r_avg, rt_avg, qv_avg);
            if (ilay==nlay) {
                Real r_hi  = cons_arr(i,j,k+1,Rho_comp);
                Real rt_hi = cons_arr(i,j,k+1,RhoTheta_comp);
                Real qv_hi = (moist) ? cons_arr(i,j,k+1,RhoQ1_comp)/r_hi : 0.0;
                r_avg  = 0.5 * (r  + r_hi);
                rt_avg = 0.5 * (rt + rt_hi);
                qv_avg = 0.5 * (qv + qv_hi);
                p_lev(icol,ilay+1) = getPgivenRTh(rt_avg, qv_avg);
                t_lev(icol,ilay+1) = getTgivenRandRTh(r_avg, rt_avg, qv_avg);
            }
        });
    }

    // Separate YAKL kernel for derived quantities
    parallel_for(SimpleBounds<2>(ncol, nlay), YAKL_LAMBDA (int icol, int ilay)
    {
        p_del(icol,ilay) = p_lev(icol,ilay+1) - p_lev(icol,ilay);
        nc(icol,ilay)    = 0.0;
        rel(icol,ilay)   = 0.0;
        rei(icol,ilay)   = 0.0;
    });

    // HACK HACK HACK
    // No LSM, so follow EAMXX dummy atmos and set constants
    yakl::memset(mu0, 0.86);
    yakl::memset(sfc_alb_dir_vis, 0.06);
    yakl::memset(sfc_alb_dir_nir, 0.06);
    yakl::memset(sfc_alb_dif_vis, 0.06);
    yakl::memset(sfc_alb_dif_nir, 0.06);

    // HACK HACK HACK
    yakl::memset(aero_tau_sw, 0.0);
    yakl::memset(aero_ssa_sw, 0.0);
    yakl::memset(aero_g_sw  , 0.0);
    yakl::memset(aero_tau_lw, 0.0);
}


void
Radiation::yakl_buffers_to_mf ()
{
    for (MFIter mfi(*m_cons_in); mfi.isValid(); ++mfi) {
        const auto& vbx      = mfi.validbox();
        const auto& sbx      = makeSlab(vbx,2,vbx.smallEnd(2));
        const int nx         = vbx.length(0);
        const int imin       = vbx.smallEnd(0);
        const int jmin       = vbx.smallEnd(1);
        const int offset     = m_col_offsets[mfi.Index()];
        const Array4<Real>& q_arr   = m_qheating_rates->const_array(mfi);
        const Array4<Real>& lsm_arr = (m_lsm_fluxes) m_lsm_fluxes->array(mfi) :
                                                     Array4<Real>{};
        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // map [i,j,k] 0-based to [icol, ilay] 1-based
            const int icol = (j-jmin)*nx + (i-imin) + 1 + offset;
            const int ilay = k+1;

            // Heating rate for SW and LW (this is for Temperature)
            q_arr(i,j,k,0) = sw_heating(icol,ilay);
            q_arr(i,j,k,1) = lw_heating(icol,ilay);
        });
        ParallelFor(sbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // map [i,j,k] 0-based to [icol, ilay] 1-based
            const int icol = (j-jmin)*nx + (i-imin) + 1 + offset;

            // SW fluxes for LSM
            lsm_arr(i,j,k,0) = sfc_flux_dir_vis(icol);
            lsm_arr(i,j,k,1) = sfc_flux_dir_nir(icol);
            lsm_arr(i,j,k,2) = sfc_flux_dif_vis(icol);
            lsm_arr(i,j,k,3) = sfc_flux_dif_nir(icol);

            // New SW flux for LSM
            lsm_arr(i,j,k,4) = sfc_flux_dir_vis(icol) + sfc_flux_dir_nir(icol)
                             + sfc_flux_dif_vis(icol) + sfc_flux_dif_nir(icol);

            // LW fluxe for LSM (at bottom surface)
            lsm_arr(i,j,k,5) = lw_flux_dn(icol,1);
        });
    }
}

void
Radiation::initialize_impl ()
{
    // Call API to initialize
    m_gas_concs.init(gas_names_yakl_offset, m_ncol, m_nlay);
    rrtmgp::rrtmgp_initialize(m_gas_concs,
                              rrtmgp_coeffs_file_sw      , rrtmgp_coeffs_file_lw      ,
                              rrtmgp_cloud_optics_file_sw, rrtmgp_cloud_optics_file_lw);
}


void
Radiation::run_impl ()
{
    // Local copies
    const auto ncol     = m_ncol;
    const auto nlay     = m_nlay;
    const auto nlwbands = m_nlwbands;
    const auto nswbands = m_nswbands;
    const auto nlwgpts  = m_nlwgpts;
    const auto do_aerosol_rad = m_do_aerosol_rad;

    // Compute orbital parameters; these are used both for computing
    // the solar zenith angle and also for computing total solar
    // irradiance scaling (tsi_scaling).
    Real obliqr, lambm0, mvelpp;
    auto orbital_year = m_orbital_year;
    auto eccen        = m_orbital_eccen;
    auto obliq        = m_orbital_obliq;
    auto mvelp        = m_orbital_mvelp;
    if (eccen >= 0 && obliq >= 0 && mvelp >= 0) {
      // fixed orbital parameters forced with orbital_year == ORB_UNDEF_INT
      orbital_year = ORB_UNDEF_INT;
    }
    orbital_params(&orbital_year, &eccen, &obliq,
                   &mvelp, &obliqr, &lambm0, &mvelpp);

    // Use the orbital parameters to calculate the solar declination and eccentricity factor
    Real delta, eccf;
    // TODO: Generalize this.
    auto calday = (m_orbital_mon-1.0)*365.0/12.0 + m_orbital_day;  // Want day + fraction; calday 1 == Jan 1 0Z
    orbital_decl(calday, eccen, mvelpp,
                 lambm0, obliqr, &delta, &eccf);

    // Overwrite eccf if using a fixed solar constant.
    auto fixed_total_solar_irradiance = m_fixed_total_solar_irradiance;
    if (fixed_total_solar_irradiance >= 0){
       eccf = fixed_total_solar_irradiance/1360.9;
    }

    // Precompute volume mixing ratio (VMR) for all gases
    //
    // H2O is obtained from qv.
    // All other comps are set to constants for now
    for (int igas(0); igas < m_ngas; ++igas) {
      auto name = m_gas_names[igas];
      auto gas_mol_weight = m_mol_weight_gas[igas];
      if (name == "H2O") {
          parallel_for(SimpleBounds<2>(ncol, nlay), YAKL_LAMBDA (int icol, int ilay)
          {
              //h2o_vmr(icol,ilay) = qv(icol,ilay) / (1.0 - qv(icol,ilay)) * mwdair/gas_mol_weight;
              tmp2d(icol,ilay) = qv(icol,ilay) * mwdair/ gas_mol_weight;
          });
      } else if (name == "CO2") {
          yakl::memset(tmp2d, m_co2vmr);
      } else if (name == "O3")  {
          yakl::memset(tmp2d, m_o3vmr );
      } else if (name == "N2O") {
          yakl::memset(tmp2d, m_n2ovmr);
      } else if (name == "CO")  {
          yakl::memset(tmp2d, m_covmr );
      } else if (name == "CH4") {
          yakl::memset(tmp2d, m_ch4vmr);
      } else if (name == "O2") {
          yakl::memset(tmp2d, m_o2vmr );
      } else if (name == "N2") {
          yakl::memset(tmp2d, m_n2vmr );
      } else {
          Abort("Radiation: Unknown gas component.");
      }

      // Populate GasConcs object
      m_gas_concs.set_vmr(name, tmp2d);
    }


    // TODO: No LSM so leaving comment for code
    // Calculate T_int from longwave flux up from the surface, assuming
    // blackbody emission with emissivity of 1.

    /*
    // NOTE: mu0 is HACKED to a constant
    //====================================

    // Determine the cosine zenith angle.
    // This must be done on HOST and copied to device.
    realHost1d h_mu0("h_mu0", ncol);
    if (m_fixed_solar_zenith_angle > 0) {
        yakl::memset(h_mu0, m_fixed_solar_zenith_angle);
    } else {
        parallel_for(ncol, YAKL_LAMBDA (int icol)
        {
            // Convert lat/lon to radians
            double lat = h_lat(icol)*PC::Pi/180.0;
            double lon = h_lon(icol)*PC::Pi/180.0;
            h_mu0(icol) = orbital_cos_zenith(calday, lat, lon, delta, m_rad_freq_in_steps * dt);
        });
    }
    mu0.deep_copy_to(h_mu0);
    */

    // Compute layer cloud mass (per unit area)
    rrtmgp::mixing_ratio_to_cloud_mass(qc, cldfrac_tot, p_del, lwp);
    rrtmgp::mixing_ratio_to_cloud_mass(qi, cldfrac_tot, p_del, iwp);

    // Convert to g/m2 (needed by RRTMGP)
    parallel_for(SimpleBounds<2>(ncol, nlay), YAKL_LAMBDA (int icol, int ilay)
    {
        lwp(icol,ilay) *= 1.e3;
        iwp(icol,ilay) *= 1.e3;
    });

    // Compute band-by-band surface_albedos. This is needed since
    // the AD passes broadband albedos, but rrtmgp require band-by-band.
    rrtmgp::compute_band_by_band_surface_albedos(ncol, nswbands,
                                                 sfc_alb_dir_vis, sfc_alb_dir_nir,
                                                 sfc_alb_dif_vis, sfc_alb_dif_nir,
                                                 sfc_alb_dir, sfc_alb_dif);

    // Run RRTMGP driver
    rrtmgp::rrtmgp_main(ncol, m_nlay,
                        p_lay, t_lay, p_lev, t_lev,
                        m_gas_concs,
                        sfc_alb_dir, sfc_alb_dif, mu0,
                        lwp, iwp, rel, rei, cldfrac_tot,
                        aero_tau_sw, aero_ssa_sw, aero_g_sw, aero_tau_lw,
                        cld_tau_sw_bnd, cld_tau_lw_bnd,
                        cld_tau_sw_gpt, cld_tau_lw_gpt,
                        sw_flux_up, sw_flux_dn, sw_flux_dn_dir, lw_flux_up, lw_flux_dn,
                        sw_clnclrsky_flux_up, sw_clnclrsky_flux_dn, sw_clnclrsky_flux_dn_dir,
                        sw_clrsky_flux_up, sw_clrsky_flux_dn, sw_clrsky_flux_dn_dir,
                        sw_clnsky_flux_up, sw_clnsky_flux_dn, sw_clnsky_flux_dn_dir,
                        lw_clnclrsky_flux_up, lw_clnclrsky_flux_dn,
                        lw_clrsky_flux_up, lw_clrsky_flux_dn,
                        lw_clnsky_flux_up, lw_clnsky_flux_dn,
                        sw_bnd_flux_up, sw_bnd_flux_dn, sw_bnd_flux_dir, lw_bnd_flux_up, lw_bnd_flux_dn,
                        eccf, m_atm_logger,
                        m_extra_clnclrsky_diag, m_extra_clnsky_diag);


    // Update heating tendency
    rrtmgp::compute_heating_rate(sw_flux_up, sw_flux_dn, r_lay, d_dz, sw_heating);
    rrtmgp::compute_heating_rate(lw_flux_up, lw_flux_dn, r_lay, d_dz, lw_heating);


    // Compute surface fluxes
    const int kbot = nlay + 1; // Should this be 1 for our layout?
    parallel_for(SimpleBounds<3>(ncol, nlay+1, nswbands), YAKL_LAMBDA (int icol, int ilay, int ibnd)
    {
        sw_bnd_flux_dif(icol,ilay,ibnd) = sw_bnd_flux_dn(icol,ilay,ibnd) - sw_bnd_flux_dir(icol,ilay,ibnd);
    });
    rrtmgp::compute_broadband_surface_fluxes(ncol, kbot, nswbands,
                                             sw_bnd_flux_dir , sw_bnd_flux_dif ,
                                             sfc_flux_dir_vis, sfc_flux_dir_nir,
                                             sfc_flux_dif_vis, sfc_flux_dif_nir);


    // TODO: Verify these are not needed, we don't have such output variables
    //=======================================================================

    // Compute diagnostic total cloud area (vertically-projected cloud cover)
    real1d cldlow ("cldlow", ncol);
    real1d cldmed ("cldmed", ncol);
    real1d cldhgh ("cldhgh", ncol);
    real1d cldtot ("cldtot", ncol);
    // NOTE: limits for low, mid, and high clouds are mostly taken from EAM F90 source, with the
    // exception that I removed the restriction on low clouds to be above (numerically lower pressures)
    // 1200 hPa, and on high clouds to be below (numerically high pressures) 50 hPa. This probably
    // does not matter in practice, as clouds probably should not be produced above 50 hPa and we
    // should not be encountering surface pressure above 1200 hPa, but in the event that things go off
    // the rails we might want to look at these still.
    rrtmgp::compute_cloud_area(ncol, nlay, nlwgpts, 700e2, std::numeric_limits<Real>::max(), p_lay, cld_tau_lw_gpt, cldlow);
    rrtmgp::compute_cloud_area(ncol, nlay, nlwgpts, 400e2,                            700e2, p_lay, cld_tau_lw_gpt, cldmed);
    rrtmgp::compute_cloud_area(ncol, nlay, nlwgpts,     0,                            400e2, p_lay, cld_tau_lw_gpt, cldhgh);
    rrtmgp::compute_cloud_area(ncol, nlay, nlwgpts,     0, std::numeric_limits<Real>::max(), p_lay, cld_tau_lw_gpt, cldtot);


    // Compute cloud-top diagnostics following AeroCOM recommendation
    auto idx_067 = rrtmgp::get_wavelength_index_sw(0.67e-6); // Get visible 0.67 micron band for COSP
    auto idx_105 = rrtmgp::get_wavelength_index_lw(10.5e-6); // Get IR 10.5 micron band for COSP
    // Compute cloud-top diagnostics following AeroCom recommendation
    real1d cdnc_at_cldtop  ("cdnc_at_cldtop" , ncol);
    real1d T_mid_at_cldtop ("T_mid_at_cldtop", ncol);
    real1d p_mid_at_cldtop ("p_mid_at_cldtop", ncol);
    real1d cldfrac_ice_at_cldtop ("cldfrac_ice_at_cldtop", ncol);
    real1d cldfrac_liq_at_cldtop ("cldfrac_liq_at_cldtop", ncol);
    real1d cldfrac_tot_at_cldtop ("cldfrac_tot_at_cldtop", ncol);
    real1d eff_radius_qc_at_cldtop ("eff_radius_qc_at_cldtop", ncol);
    real1d eff_radius_qi_at_cldtop ("eff_radius_qi_at_cldtop", ncol);
    rrtmgp::compute_aerocom_cloudtop(ncol, nlay, t_lay, p_lay, p_del, z_del, qc, qi, rel, rei, cldfrac_tot,
                                     nc, T_mid_at_cldtop, p_mid_at_cldtop, cldfrac_ice_at_cldtop,
                                     cldfrac_liq_at_cldtop, cldfrac_tot_at_cldtop, cdnc_at_cldtop,
                                     eff_radius_qc_at_cldtop, eff_radius_qi_at_cldtop);
}


void
Radiation::finalize_impl ()
{
    // Finish rrtmgp
    m_gas_concs.reset();
    rrtmgp::rrtmgp_finalize();

    // Fill the AMReX MFs from YAKL Arrays
    yakl_buffers_to_mf();

    // Deallocate the buffer arrays
    dealloc_buffers();
}
