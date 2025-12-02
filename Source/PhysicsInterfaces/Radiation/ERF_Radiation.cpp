/*
 * RTE-RRTMGP radiation model interface to ERF
 * The original code is developed by RobertPincus, and the code is open source available at:
 *                        https://github.com/earth-system-radiation/rte-rrtmgp
 * Please reference to the following paper,
 *                        https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019MS001621
 * NOTE: we use the C++ version of RTE-RRTMGP, which is reimplemented the original Fortran
 * code using C++ KOKKOS for CUDA, HiP and SYCL application by E3SM ECP team, the C++ version
 * of the rte-rrtmgp code is located at:
 *                       https://github.com/E3SM-Project/rte-rrtmgp
 * The RTE-RRTMGP uses BSD-3-Clause Open Source License, if you want to make changes,
 * and modifications to the code, please refer to BSD-3-Clause Open Source License.
 */

#include "ERF_NCInterface.H"
#include "ERF_Radiation.H"

using namespace amrex;

Radiation::Radiation (const int& lev,
                      SolverChoice& sc)
{
    // Note that Kokkos is now initialized in main.cpp

    // Check if we have a valid moisture model
    if (sc.moisture_type != MoistureType::None) { m_moist = true; }

    // Check if we have a moisture model with ice
    if (sc.moisture_type == MoistureType::SAM)  { m_ice = true; }

    // Check if we have a land surface model enabled
    if (sc.lsm_type != LandSurfaceType::None) { m_lsm = true; }

    // Construct parser object for following reads
    ParmParse pp("erf");

    // Must specify a surface temp (LSM can overwrite)
    pp.get("rad_t_sfc",m_rad_t_sfc);

    // Radiation timestep, as a number of atm steps
    pp.query("rad_freq_in_steps", m_rad_freq_in_steps);

    // Flag to write fluxes to plt file
    pp.query("rad_write_fluxes", m_rad_write_fluxes);

    // Do MCICA subcolumn sampling
    pp.query("rad_do_subcol_sampling", m_do_subcol_sampling);

    // Determine orbital year. If orbital_year is negative, use current year
    // from timestamp for orbital year; if positive, use provided orbital year
    // for duration of simulation.
    m_fixed_orbital_year = pp.query("rad_orbital_year", m_orbital_year);

    // Get orbital parameters from inputs file
    pp.query("rad_orbital_eccentricity", m_orbital_eccen);
    pp.query("rad_orbital_obliquity"   , m_orbital_obliq);
    pp.query("rad_orbital_mvelp"       , m_orbital_mvelp);

    // Get a constant lat/lon for idealized simulations
    pp.query("rad_cons_lat", m_lat_cons);
    pp.query("rad_cons_lon", m_lon_cons);

    // Value for prescribing an invariant solar constant (i.e. total solar irradiance at
    // TOA).  Used for idealized experiments such as RCE. Disabled when value is less than 0.
    pp.query("fixed_total_solar_irradiance", m_fixed_total_solar_irradiance);

    // Determine whether or not we are using a fixed solar zenith angle (positive value)
    pp.query("fixed_solar_zenith_angle", m_fixed_solar_zenith_angle);

    // Get prescribed surface values of greenhouse gases
    pp.query("co2vmr", m_co2vmr);
    pp.queryarr("o3vmr" , m_o3vmr );
    pp.query("n2ovmr", m_n2ovmr);
    pp.query("covmr" , m_covmr );
    pp.query("ch4vmr", m_ch4vmr);
    pp.query("o2vmr" , m_o2vmr );
    pp.query("n2vmr" , m_n2vmr );

    // Required aerosol optical properties from SPA
    pp.query("rad_do_aerosol", m_do_aerosol_rad);

    // Whether we do extra clean/clear sky calculations
    pp.query("rad_extra_clnclrsky_diag", m_extra_clnclrsky_diag);
    pp.query("rad_extra_clnsky_diag"   , m_extra_clnsky_diag);

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

    // Get dimensions from lookup data
    if (ParallelDescriptor::IOProcessor()) {
        auto ncf_sw = ncutils::NCFile::open(rrtmgp_coeffs_file_sw, NC_CLOBBER | NC_NETCDF4);
        m_nswbands = ncf_sw.dim("bnd").len();
        m_nswgpts  = ncf_sw.dim("gpt").len();
        ncf_sw.close();

        auto ncf_lw = ncutils::NCFile::open(rrtmgp_coeffs_file_lw, NC_CLOBBER | NC_NETCDF4);
        m_nlwbands = ncf_lw.dim("bnd").len();
        m_nlwgpts  = ncf_lw.dim("gpt").len();
        ncf_lw.close();
    }
    int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank
    ParallelDescriptor::Bcast(&m_nswbands, 1, ioproc);
    ParallelDescriptor::Bcast(&m_nlwbands, 1, ioproc);
    ParallelDescriptor::Bcast(&m_nswgpts,  1, ioproc);
    ParallelDescriptor::Bcast(&m_nlwgpts,  1, ioproc);

    // Output for user
    if (lev == 0) {
        Print() << "Radiation interface constructed:\n";
        Print() << "========================================================\n";
        Print() << "Coeff SW file: " << rrtmgp_coeffs_file_sw << "\n";
        Print() << "Coeff LW file: " << rrtmgp_coeffs_file_lw << "\n";
        Print() << "Cloud SW file: " << rrtmgp_cloud_optics_file_sw << "\n";
        Print() << "Cloud LW file: " << rrtmgp_cloud_optics_file_lw << "\n";
        Print() << "Number of short/longwave bands: "
                << m_nswbands << " " << m_nlwbands << "\n";
        Print() << "Number of short/longwave gauss points: "
                << m_nswgpts  << " " << m_nlwgpts  << "\n";
        Print() << "========================================================\n";
    }
}

void
Radiation::set_grids (int& level,
                      int& step,
                      Real& time,
                      const Real& dt,
                      const BoxArray& ba,
                      Geometry& geom,
                      MultiFab* cons_in,
                      iMultiFab* lmask,
                      MultiFab*  t_surf,
                      MultiFab* lsm_fluxes,
                      MultiFab* lsm_zenith,
                      Vector<MultiFab*>& lsm_input_ptrs,
                      MultiFab* qheating_rates,
                      MultiFab* rad_fluxes,
                      MultiFab* z_phys,
                      MultiFab* lat,
                      MultiFab* lon)

{
    // Set data members that may change
    m_lev            = level;
    m_step           = step;
    m_time           = time;
    m_dt             = dt;
    m_geom           = geom;
    m_cons_in        = cons_in;
    m_lsm_fluxes     = lsm_fluxes;
    m_lsm_zenith     = lsm_zenith;
    m_qheating_rates = qheating_rates;
    m_rad_fluxes     = rad_fluxes;
    m_z_phys         = z_phys;
    m_lat            = lat;
    m_lon            = lon;

    // Update the day and month
    time_t timestamp = time_t(time);
    struct tm *timeinfo = gmtime(&timestamp);
    if (m_fixed_orbital_year) {
        m_orbital_mon  = timeinfo->tm_mon + 1;
        m_orbital_day  = timeinfo->tm_mday;
        m_orbital_sec  = timeinfo->tm_hour*3600 + timeinfo->tm_min*60 + timeinfo->tm_sec;
    } else {
        m_orbital_year = timeinfo->tm_year + 1900;
        m_orbital_mon  = timeinfo->tm_mon  + 1;
        m_orbital_day  = timeinfo->tm_mday;
        m_orbital_sec  = timeinfo->tm_hour*3600 + timeinfo->tm_min*60 + timeinfo->tm_sec;
    }

    // Only allocate and proceed if we are going to update radiation
    m_update_rad = false;
    if (m_rad_freq_in_steps > 0) { m_update_rad = ( (m_step == 0) || (m_step % m_rad_freq_in_steps == 0) ); }

    if (m_update_rad) {
        // Call to Init() has set the dimensions: ncol & nlay

        // Allocate the buffer arrays
        alloc_buffers();

        // Fill the KOKKOS Views from AMReX MFs
        mf_to_kokkos_buffers(lmask, t_surf, lsm_input_ptrs);

        // Initialize datalog MF on first step
        if (m_first_step) {
            m_first_step = false;
            datalog_mf.define(cons_in->boxArray(), cons_in->DistributionMap(), 25, 0);
            datalog_mf.setVal(0.0);
        }
    }
}

void
Radiation::alloc_buffers ()
{
    // 1d size (m_ngas)
    const Real* mol_weight_gas_p = m_mol_weight_gas.data();
    const std::string* gas_names_p = m_gas_names.data();
    m_gas_mol_weights = real1d_k("m_gas_mol_weights", m_ngas);
    realHost1d_k m_gas_mol_weights_h("m_gas_mol_weights_h", m_ngas);
    gas_names_offset.clear(); gas_names_offset.resize(m_ngas);
    std::string* gas_names_offset_p = gas_names_offset.data();
    Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::Serial>(0, m_ngas),
                         [&] (int igas)
    {
        m_gas_mol_weights_h(igas) = mol_weight_gas_p[igas];
        gas_names_offset_p[igas]  = gas_names_p[igas];
    });
    Kokkos::deep_copy(m_gas_mol_weights, m_gas_mol_weights_h);

    // 1d size (1 or nlay)
    m_o3_size = m_o3vmr.size();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(((m_o3_size==1) || (m_o3_size==m_nlay)),
                                     "O3 VMR array must be length 1 or nlay");
    Real* o3vmr_p = m_o3vmr.data();
    o3_lay = real1d_k("o3_lay", m_o3_size);
    realHost1d_k o3_lay_h("o3_lay_h", m_o3_size);
    Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::Serial>(0, m_o3_size),
                         [&] (int io3)
    {
        o3_lay_h(io3) = o3vmr_p[io3];
    });
    Kokkos::deep_copy(o3_lay, o3_lay_h);

    // 1d size (ncol)
    mu0              = real1d_k("mu0"             , m_ncol);
    sfc_alb_dir_vis  = real1d_k("sfc_alb_dir_vis" , m_ncol);
    sfc_alb_dir_nir  = real1d_k("sfc_alb_dir_nir" , m_ncol);
    sfc_alb_dif_vis  = real1d_k("sfc_alb_dif_vis" , m_ncol);
    sfc_alb_dif_nir  = real1d_k("sfc_alb_dif_nir" , m_ncol);
    sfc_flux_dir_vis = real1d_k("sfc_flux_dir_vis", m_ncol);
    sfc_flux_dir_nir = real1d_k("sfc_flux_dir_nir", m_ncol);
    sfc_flux_dif_vis = real1d_k("sfc_flux_dif_vis", m_ncol);
    sfc_flux_dif_nir = real1d_k("sfc_flux_dif_nir", m_ncol);
    lat              = real1d_k("lat"             , m_ncol);
    lon              = real1d_k("lon"             , m_ncol);
    sfc_emis         = real1d_k("sfc_emis"        , m_ncol);
    t_sfc            = real1d_k("t_sfc"           , m_ncol);
    lw_src           = real1d_k("lw_src"          , m_ncol);

    // 2d size (ncol, nlay)
    r_lay         = real2d_k("r_lay"        , m_ncol, m_nlay);
    p_lay         = real2d_k("p_lay"        , m_ncol, m_nlay);
    t_lay         = real2d_k("t_lay"        , m_ncol, m_nlay);
    z_del         = real2d_k("z_del"        , m_ncol, m_nlay);
    qv_lay        = real2d_k("qv"           , m_ncol, m_nlay);
    qc_lay        = real2d_k("qc"           , m_ncol, m_nlay);
    qi_lay        = real2d_k("qi"           , m_ncol, m_nlay);
    cldfrac_tot   = real2d_k("cldfrac_tot"  , m_ncol, m_nlay);
    eff_radius_qc = real2d_k("eff_radius_qc", m_ncol, m_nlay);
    eff_radius_qi = real2d_k("eff_radius_qi", m_ncol, m_nlay);
    lwp           = real2d_k("lwp"          , m_ncol, m_nlay);
    iwp           = real2d_k("iwp"          , m_ncol, m_nlay);
    sw_heating    = real2d_k("sw_heating"   , m_ncol, m_nlay);
    lw_heating    = real2d_k("lw_heating"   , m_ncol, m_nlay);
    sw_clrsky_heating = real2d_k("sw_clrsky_heating", m_ncol, m_nlay);
    lw_clrsky_heating = real2d_k("lw_clrsky_heating", m_ncol, m_nlay);

    // 2d size (ncol, nlay+1)
    d_tint                   = real2d_k("d_tint"                  , m_ncol, m_nlay+1);
    p_lev                    = real2d_k("p_lev"                   , m_ncol, m_nlay+1);
    t_lev                    = real2d_k("t_lev"                   , m_ncol, m_nlay+1);

    sw_flux_up               = real2d_k("sw_flux_up"              , m_ncol, m_nlay+1);
    sw_flux_dn               = real2d_k("sw_flux_dn"              , m_ncol, m_nlay+1);
    sw_flux_dn_dir           = real2d_k("sw_flux_dn_dir"          , m_ncol, m_nlay+1);

    lw_flux_up               = real2d_k("lw_flux_up"              , m_ncol, m_nlay+1);
    lw_flux_dn               = real2d_k("lw_flux_dn"              , m_ncol, m_nlay+1);

    sw_clnclrsky_flux_up     = real2d_k("sw_clnclrsky_flux_up"    , m_ncol, m_nlay+1);
    sw_clnclrsky_flux_dn     = real2d_k("sw_clnclrsky_flux_dn"    , m_ncol, m_nlay+1);
    sw_clnclrsky_flux_dn_dir = real2d_k("sw_clnclrsky_flux_dn_dir", m_ncol, m_nlay+1);
    sw_clrsky_flux_up        = real2d_k("sw_clrsky_flux_up"       , m_ncol, m_nlay+1);
    sw_clrsky_flux_dn        = real2d_k("sw_clrsky_flux_dn"       , m_ncol, m_nlay+1);
    sw_clrsky_flux_dn_dir    = real2d_k("sw_clrsky_flux_dn_dir"   , m_ncol, m_nlay+1);
    sw_clnsky_flux_up        = real2d_k("sw_clnsky_flux_up"       , m_ncol, m_nlay+1);
    sw_clnsky_flux_dn        = real2d_k("sw_clnsky_flux_dn"       , m_ncol, m_nlay+1);
    sw_clnsky_flux_dn_dir    = real2d_k("sw_clnsky_flux_dn_dir"   , m_ncol, m_nlay+1);

    lw_clnclrsky_flux_up     = real2d_k("lw_clnclrsky_flux_up"    , m_ncol, m_nlay+1);
    lw_clnclrsky_flux_dn     = real2d_k("lw_clnclrsky_flux_dn"    , m_ncol, m_nlay+1);
    lw_clrsky_flux_up        = real2d_k("lw_clrsky_flux_up"       , m_ncol, m_nlay+1);
    lw_clrsky_flux_dn        = real2d_k("lw_clrsky_flux_dn"       , m_ncol, m_nlay+1);
    lw_clnsky_flux_up        = real2d_k("lw_clnsky_flux_up"       , m_ncol, m_nlay+1);
    lw_clnsky_flux_dn        = real2d_k("lw_clnsky_flux_dn"       , m_ncol, m_nlay+1);

    // 3d size (ncol, nlay+1, nswbands)
    sw_bnd_flux_up  = real3d_k("sw_bnd_flux_up" , m_ncol, m_nlay+1, m_nswbands);
    sw_bnd_flux_dn  = real3d_k("sw_bnd_flux_dn" , m_ncol, m_nlay+1, m_nswbands);
    sw_bnd_flux_dir = real3d_k("sw_bnd_flux_dir", m_ncol, m_nlay+1, m_nswbands);
    sw_bnd_flux_dif = real3d_k("sw_bnd_flux_dif", m_ncol, m_nlay+1, m_nswbands);

    // 3d size (ncol, nlay+1, nlwbands)
    lw_bnd_flux_up = real3d_k("lw_bnd_flux_up" , m_ncol, m_nlay+1, m_nlwbands);
    lw_bnd_flux_dn = real3d_k("lw_bnd_flux_dn" , m_ncol, m_nlay+1, m_nlwbands);

    // 2d size (ncol, nswbands)
    sfc_alb_dir = real2d_k("sfc_alb_dir", m_ncol, m_nswbands);
    sfc_alb_dif = real2d_k("sfc_alb_dif", m_ncol, m_nswbands);

    // 2d size (ncol, nlwbands)
    //emis_sfc    = real2d_k("emis_sfc", m_ncol, m_nlwbands);

    /*
    // 3d size (ncol, nlay, n[sw,lw]bands)
    aero_tau_sw = real3d_k("aero_tau_sw", m_ncol, m_nlay, m_nswbands);
    aero_ssa_sw = real3d_k("aero_ssa_sw", m_ncol, m_nlay, m_nswbands);
    aero_g_sw   = real3d_k("aero_g_sw"  , m_ncol, m_nlay, m_nswbands);
    aero_tau_lw = real3d_k("aero_tau_lw", m_ncol, m_nlay, m_nlwbands);

    // 3d size (ncol, nlay, n[sw,lw]bnds)
    cld_tau_sw_bnd = real3d_k("cld_tau_sw_bnd", m_ncol, m_nlay, m_nswbands);
    cld_tau_lw_bnd = real3d_k("cld_tau_lw_bnd", m_ncol, m_nlay, m_nlwbands);

    // 3d size (ncol, nlay, n[sw,lw]gpts)
    cld_tau_sw_gpt = real3d_k("cld_tau_sw_gpt", m_ncol, m_nlay, m_nswgpts);
    cld_tau_lw_gpt = real3d_k("cld_tau_lw_gpt", m_ncol, m_nlay, m_nlwgpts);
    */
}

void
Radiation::dealloc_buffers ()
{
    // 1d size (m_ngas)
    m_gas_mol_weights = real1d_k();

    // 1d size (1 or nlay)
    o3_lay = real1d_k();

    // 1d size (ncol)
    mu0              = real1d_k();
    sfc_alb_dir_vis  = real1d_k();
    sfc_alb_dir_nir  = real1d_k();
    sfc_alb_dif_vis  = real1d_k();
    sfc_alb_dif_nir  = real1d_k();
    sfc_flux_dir_vis = real1d_k();
    sfc_flux_dir_nir = real1d_k();
    sfc_flux_dif_vis = real1d_k();
    sfc_flux_dif_nir = real1d_k();
    lat              = real1d_k();
    lon              = real1d_k();
    sfc_emis         = real1d_k();
    t_sfc            = real1d_k();
    lw_src           = real1d_k();

    // 2d size (ncol, nlay)
    r_lay             = real2d_k();
    p_lay             = real2d_k();
    t_lay             = real2d_k();
    z_del             = real2d_k();
    qv_lay            = real2d_k();
    qc_lay            = real2d_k();
    qi_lay            = real2d_k();
    cldfrac_tot       = real2d_k();
    eff_radius_qc     = real2d_k();
    eff_radius_qi     = real2d_k();
    lwp               = real2d_k();
    iwp               = real2d_k();
    sw_heating        = real2d_k();
    lw_heating        = real2d_k();
    sw_clrsky_heating = real2d_k();
    lw_clrsky_heating = real2d_k();
    sw_heating        = real2d_k();
    lw_heating        = real2d_k();
    sw_clrsky_heating = real2d_k();
    lw_clrsky_heating = real2d_k();

    // 2d size (ncol, nlay+1)
    d_tint                   = real2d_k();
    p_lev                    = real2d_k();
    t_lev                    = real2d_k();
    sw_flux_up               = real2d_k();
    sw_flux_dn               = real2d_k();
    sw_flux_dn_dir           = real2d_k();
    lw_flux_up               = real2d_k();
    lw_flux_dn               = real2d_k();
    sw_clnclrsky_flux_up     = real2d_k();
    sw_clnclrsky_flux_dn     = real2d_k();
    sw_clnclrsky_flux_dn_dir = real2d_k();
    sw_clrsky_flux_up        = real2d_k();
    sw_clrsky_flux_dn        = real2d_k();
    sw_clrsky_flux_dn_dir    = real2d_k();
    sw_clnsky_flux_up        = real2d_k();
    sw_clnsky_flux_dn        = real2d_k();
    sw_clnsky_flux_dn_dir    = real2d_k();
    lw_clnclrsky_flux_up     = real2d_k();
    lw_clnclrsky_flux_dn     = real2d_k();
    lw_clrsky_flux_up        = real2d_k();
    lw_clrsky_flux_dn        = real2d_k();
    lw_clnsky_flux_up        = real2d_k();
    lw_clnsky_flux_dn        = real2d_k();

    // 3d size (ncol, nlay+1, nswbands)
    sw_bnd_flux_up  = real3d_k();
    sw_bnd_flux_dn  = real3d_k();
    sw_bnd_flux_dir = real3d_k();
    sw_bnd_flux_dif = real3d_k();

    // 3d size (ncol, nlay+1, nlwbands)
    lw_bnd_flux_up = real3d_k();
    lw_bnd_flux_dn = real3d_k();

    // 2d size (ncol, nswbands)
    sfc_alb_dir = real2d_k();
    sfc_alb_dif = real2d_k();

    // 2d size (ncol, nlwbands)
    //emis_sfc = real2d_k();

    /*
    // 3d size (ncol, nlay, n[sw,lw]bands)
    aero_tau_sw = real3d_k();
    aero_ssa_sw = real3d_k();
    aero_g_sw   = real3d_k();
    aero_tau_lw = real3d_k();

    // 3d size (ncol, nlay, n[sw,lw]bnds)
    cld_tau_sw_bnd = real3d_k();
    cld_tau_lw_bnd = real3d_k();

    // 3d size (ncol, nlay, n[sw,lw]gpts)
    cld_tau_sw_gpt = real3d_k();
    cld_tau_lw_gpt = real3d_k();
    */
}


void
Radiation::mf_to_kokkos_buffers (iMultiFab* lmask,
                                 MultiFab*  t_surf,
                                 Vector<MultiFab*>& lsm_input_ptrs)
{
    Table2D<Real,Order::C> r_lay_tab(r_lay.data(), {0,0}, {static_cast<int>(r_lay.extent(0)),static_cast<int>(r_lay.extent(1))});
    Table2D<Real,Order::C> p_lay_tab(p_lay.data(), {0,0}, {static_cast<int>(p_lay.extent(0)),static_cast<int>(p_lay.extent(1))});
    Table2D<Real,Order::C> t_lay_tab(t_lay.data(), {0,0}, {static_cast<int>(t_lay.extent(0)),static_cast<int>(t_lay.extent(1))});
    Table2D<Real,Order::C> z_del_tab(z_del.data(), {0,0}, {static_cast<int>(z_del.extent(0)),static_cast<int>(z_del.extent(1))});
    Table2D<Real,Order::C> qv_lay_tab(qv_lay.data(), {0,0}, {static_cast<int>(qv_lay.extent(0)),static_cast<int>(qv_lay.extent(1))});
    Table2D<Real,Order::C> qc_lay_tab(qc_lay.data(), {0,0}, {static_cast<int>(qc_lay.extent(0)),static_cast<int>(qc_lay.extent(1))});
    Table2D<Real,Order::C> qi_lay_tab(qi_lay.data(), {0,0}, {static_cast<int>(qi_lay.extent(0)),static_cast<int>(qi_lay.extent(1))});
    Table2D<Real,Order::C> cldfrac_tot_tab(cldfrac_tot.data(), {0,0}, {static_cast<int>(cldfrac_tot.extent(0)),static_cast<int>(cldfrac_tot.extent(1))});

    Table2D<Real,Order::C> lwp_tab(lwp.data(), {0,0}, {static_cast<int>(lwp.extent(0)),static_cast<int>(lwp.extent(1))});
    Table2D<Real,Order::C> iwp_tab(iwp.data(), {0,0}, {static_cast<int>(iwp.extent(0)),static_cast<int>(iwp.extent(1))});
    Table2D<Real,Order::C> eff_radius_qc_tab(eff_radius_qc.data(), {0,0}, {static_cast<int>(eff_radius_qc.extent(0)),static_cast<int>(eff_radius_qc.extent(1))});
    Table2D<Real,Order::C> eff_radius_qi_tab(eff_radius_qi.data(), {0,0}, {static_cast<int>(eff_radius_qi.extent(0)),static_cast<int>(eff_radius_qi.extent(1))});

    Table2D<Real,Order::C> p_lev_tab(p_lev.data(), {0,0}, {static_cast<int>(p_lev.extent(0)),static_cast<int>(p_lev.extent(1))});
    Table2D<Real,Order::C> t_lev_tab(t_lev.data(), {0,0}, {static_cast<int>(t_lev.extent(0)),static_cast<int>(t_lev.extent(1))});

    Table1D<Real> lat_tab(lat.data(), {0}, {static_cast<int>(lat.extent(0))});
    Table1D<Real> lon_tab(lon.data(), {0}, {static_cast<int>(lon.extent(0))});
    Table1D<Real> t_sfc_tab(t_sfc.data(), {0}, {static_cast<int>(t_sfc.extent(0))});

    bool moist = m_moist;
    bool ice   = m_ice;
    const bool has_lsm = m_lsm;
    const bool has_lat = m_lat;
    const bool has_lon = m_lon;
    const bool has_surflayer = (t_surf);
    int  ncol  = m_ncol;
    int  nlay  = m_nlay;
    Real dz    = m_geom.CellSize(2);
    Real cons_lat = m_lat_cons;
    Real cons_lon = m_lon_cons;
    Real rad_t_sfc = m_rad_t_sfc;

    for (MFIter mfi(*m_cons_in); mfi.isValid(); ++mfi) {
        const auto& vbx  = mfi.validbox();
        const int nx     = vbx.length(0);
        const int imin   = vbx.smallEnd(0);
        const int jmin   = vbx.smallEnd(1);
        const int offset = m_col_offsets[mfi.index()];
        const Array4<const Real>& cons_arr = m_cons_in->const_array(mfi);
        const Array4<const Real>& z_arr    = (m_z_phys) ? m_z_phys->const_array(mfi) :
                                                          Array4<const Real>{};
        const Array4<const Real>& lat_arr  = (m_lat)    ? m_lat->const_array(mfi) :
                                                          Array4<const Real>{};
        const Array4<const Real>& lon_arr  = (m_lon)    ? m_lon->const_array(mfi) :
                                                          Array4<const Real>{};
        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // map [i,j,k] 0-based to [icol, ilay] 0-based
            const int icol   = (j-jmin)*nx + (i-imin) + offset;
            const int ilay   = k;

            // EOS input (at CC)
            Real r  = cons_arr(i,j,k,Rho_comp);
            Real rt = cons_arr(i,j,k,RhoTheta_comp);
            Real qv = (moist) ? std::max(cons_arr(i,j,k,RhoQ1_comp)/r,0.0) : 0.0;
            Real qc = (moist) ? std::max(cons_arr(i,j,k,RhoQ2_comp)/r,0.0) : 0.0;
            Real qi = (ice)   ? std::max(cons_arr(i,j,k,RhoQ3_comp)/r,0.0) : 0.0;

            // EOS avg to z-face
            Real r_lo   = cons_arr(i,j,k-1,Rho_comp);
            Real rt_lo  = cons_arr(i,j,k-1,RhoTheta_comp);
            Real qv_lo  = (moist) ? cons_arr(i,j,k-1,RhoQ1_comp)/r_lo : 0.0;
            Real r_avg  = 0.5 * (r  + r_lo);
            Real rt_avg = 0.5 * (rt + rt_lo);
            Real qv_avg = 0.5 * (qv + qv_lo);

            // Views at CC
            r_lay_tab(icol,ilay) = r;

            p_lay_tab(icol,ilay) = getPgivenRTh(rt, qv);
            t_lay_tab(icol,ilay) = getTgivenRandRTh(r, rt, qv);
            z_del_tab(icol,ilay) = (z_arr) ? 0.25 * ( (z_arr(i  ,j  ,k+1) - z_arr(i  ,j  ,k))
                                                    + (z_arr(i+1,j  ,k+1) - z_arr(i+1,j  ,k))
                                                    + (z_arr(i  ,j+1,k+1) - z_arr(i  ,j+1,k))
                                                    + (z_arr(i+1,j+1,k+1) - z_arr(i+1,j+1,k)) ) : dz;
            qv_lay_tab(icol,ilay) = qv;
            qc_lay_tab(icol,ilay) = qc;
            qi_lay_tab(icol,ilay) = qi;
            cldfrac_tot_tab(icol,ilay) = ((qc+qi)>0.0) ? 1. : 0.;

            // NOTE: These are populated in 'mixing_ratio_to_cloud_mass'
            lwp_tab(icol,ilay) = 0.0;
            iwp_tab(icol,ilay) = 0.0;

            // NOTE: These would be populated from P3 (we use the constants in p3_main_impl.hpp)
            eff_radius_qc_tab(icol,ilay) = (qc>0.0) ? 10.0e-6 : 0.0;
            eff_radius_qi_tab(icol,ilay) = (qi>0.0) ? 25.0e-6 : 0.0;

            // Buffers on z-faces (nlay+1)
            p_lev_tab(icol,ilay) = getPgivenRTh(rt_avg, qv_avg);
            t_lev_tab(icol,ilay) = getTgivenRandRTh(r_avg, rt_avg, qv_avg);
            if (ilay==(nlay-1)) {
                Real r_hi  = cons_arr(i,j,k+1,Rho_comp);
                Real rt_hi = cons_arr(i,j,k+1,RhoTheta_comp);
                Real qv_hi = (moist) ? std::max(cons_arr(i,j,k+1,RhoQ1_comp)/r_hi,0.0) : 0.0;
                r_avg  = 0.5 * (r  + r_hi);
                rt_avg = 0.5 * (rt + rt_hi);
                qv_avg = 0.5 * (qv + qv_hi);
                p_lev_tab(icol,ilay+1) = getPgivenRTh(rt_avg, qv_avg);
                t_lev_tab(icol,ilay+1) = getTgivenRandRTh(r_avg, rt_avg, qv_avg);
            }

            // 1D data structures
            if (k==0) {
                lat_tab(icol) = (has_lat) ? lat_arr(i,j,0) : cons_lat;
                lon_tab(icol) = (has_lon) ? lon_arr(i,j,0) : cons_lon;
            }

        });
    } // mfi

    // Populate vars LSM would provide
    if (!has_lsm && !has_surflayer) {
        // Parsed surface temp
        Kokkos::deep_copy(t_sfc, rad_t_sfc);

        // EAMXX dummy atmos constants
        Kokkos::deep_copy(sfc_alb_dir_vis, 0.06);
        Kokkos::deep_copy(sfc_alb_dir_nir, 0.06);
        Kokkos::deep_copy(sfc_alb_dif_vis, 0.06);
        Kokkos::deep_copy(sfc_alb_dif_nir, 0.06);

        // AML NOTE: These are not used in current EAMXX, I've left
        //           the code to plug into these if we need it.
        //
        // Current EAMXX constants
        Kokkos::deep_copy(sfc_emis, 0.98);
        Kokkos::deep_copy(lw_src  , 0.0 );
    } else {
        Vector<real1d_k> rrtmgp_in_vars = {t_sfc, sfc_emis,
                                           sfc_alb_dir_vis, sfc_alb_dir_nir,
                                           sfc_alb_dif_vis, sfc_alb_dif_nir};
        Vector<Real> rrtmgp_default_vals = {rad_t_sfc, 0.98,
                                            0.06, 0.06,
                                            0.06, 0.06};
        for (int ivar(0); ivar<lsm_input_ptrs.size(); ivar++) {
            auto rrtmgp_default_val = rrtmgp_default_vals[ivar];
            auto rrtmgp_to_fill = rrtmgp_in_vars[ivar];
            for (MFIter mfi(*m_cons_in); mfi.isValid(); ++mfi) {
                const auto& vbx  = mfi.validbox();
                const auto& sbx  = makeSlab(vbx,2,vbx.smallEnd(2));
                const int nx     = vbx.length(0);
                const int imin   = vbx.smallEnd(0);
                const int jmin   = vbx.smallEnd(1);
                const int offset = m_col_offsets[mfi.index()];
                const Array4<const int>& lmask_arr  = (lmask)   ? lmask->const_array(mfi) :
                                                                  Array4<const int> {};
                const Array4<const Real>& tsurf_arr  = (t_surf) ? t_surf->const_array(mfi) :
                                                                  Array4<const Real> {};
                const Array4<const Real>& lsm_in_arr = (lsm_input_ptrs[ivar]) ? lsm_input_ptrs[ivar]->const_array(mfi) :
                                                                                Array4<const Real> {};
                ParallelFor(sbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    // map [i,j,k] 0-based to [icol, ilay] 0-based
                    const int icol   = (j-jmin)*nx + (i-imin) + offset;

                    // Check if over land
                    bool is_land = (lmask_arr) ? lmask_arr(i,j,k) : 1;

                    // Check if valid LSM data
                    bool valid_lsm_data{false};
                    if (lsm_in_arr) { valid_lsm_data = (lsm_in_arr(i,j,k) > 0.); }

                    // Have LSM and are over land
                    if (is_land && valid_lsm_data) {
                        rrtmgp_to_fill(icol) = lsm_in_arr(i,j,k);
                    }
                    // We have a SurfLayer (enforce consistency with temperature)
                    else if (tsurf_arr && (ivar==0)) {
                        rrtmgp_to_fill(icol) = tsurf_arr(i,j,k);
                    }
                    // Use the default value
                    else {
                        rrtmgp_to_fill(icol) = rrtmgp_default_val;
                    }
                });
            } //mfi
        } // ivar
        Kokkos::deep_copy(lw_src, 0.0 );
    } // have lsm

    // Enforce consistency between t_sfc and t_lev at bottom surface
    Kokkos::parallel_for(Kokkos::RangePolicy(0, ncol),
                         KOKKOS_LAMBDA (int icol)
    {
        t_lev_tab(icol,0) = t_sfc_tab(icol);
    });
}


void
Radiation::kokkos_buffers_to_mf (Vector<MultiFab*>& lsm_output_ptrs)
{
    // Heating rate, fluxes, zenith, lsm ptrs
    Vector<real2d_k> rrtmgp_out_vars = {sw_flux_dn, lw_flux_dn};

    Table2D<Real,Order::C> p_lay_tab(p_lay.data(), {0,0}, {static_cast<int>(p_lay.extent(0)),static_cast<int>(p_lay.extent(1))});
    Table2D<Real,Order::C> sw_heating_tab(sw_heating.data(), {0,0}, {static_cast<int>(sw_heating.extent(0)),static_cast<int>(sw_heating.extent(1))});
    Table2D<Real,Order::C> lw_heating_tab(lw_heating.data(), {0,0}, {static_cast<int>(lw_heating.extent(0)),static_cast<int>(lw_heating.extent(1))});
    Table2D<Real,Order::C> sw_flux_up_tab(sw_flux_up.data(), {0,0}, {static_cast<int>(sw_flux_up.extent(0)),static_cast<int>(sw_flux_up.extent(1))});
    Table2D<Real,Order::C> sw_flux_dn_tab(sw_flux_dn.data(), {0,0}, {static_cast<int>(sw_flux_dn.extent(0)),static_cast<int>(sw_flux_dn.extent(1))});
    Table2D<Real,Order::C> lw_flux_up_tab(lw_flux_up.data(), {0,0}, {static_cast<int>(lw_flux_up.extent(0)),static_cast<int>(lw_flux_up.extent(1))});
    Table2D<Real,Order::C> lw_flux_dn_tab(lw_flux_dn.data(), {0,0}, {static_cast<int>(lw_flux_dn.extent(0)),static_cast<int>(lw_flux_dn.extent(1))});

    Table1D<Real> sfc_flux_dir_vis_tab(sfc_flux_dir_vis.data(), {0}, {static_cast<int>(sfc_flux_dir_vis.extent(0))});
    Table1D<Real> sfc_flux_dir_nir_tab(sfc_flux_dir_nir.data(), {0}, {static_cast<int>(sfc_flux_dir_nir.extent(0))});
    Table1D<Real> sfc_flux_dif_vis_tab(sfc_flux_dif_vis.data(), {0}, {static_cast<int>(sfc_flux_dif_vis.extent(0))});
    Table1D<Real> sfc_flux_dif_nir_tab(sfc_flux_dif_nir.data(), {0}, {static_cast<int>(sfc_flux_dif_nir.extent(0))});
    Table1D<Real>              mu0_tab(mu0.data(),              {0}, {static_cast<int>(mu0.extent(0))});

    // Expose for device

    for (MFIter mfi(*m_cons_in); mfi.isValid(); ++mfi) {
        const auto& vbx      = mfi.validbox();
        const auto& sbx      = makeSlab(vbx,2,vbx.smallEnd(2));
        const int nx         = vbx.length(0);
        const int imin       = vbx.smallEnd(0);
        const int jmin       = vbx.smallEnd(1);
        const int offset     = m_col_offsets[mfi.index()];
        const Array4<Real>& q_arr = m_qheating_rates->array(mfi);
        const Array4<Real>& f_arr = m_rad_fluxes->array(mfi);
        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // map [i,j,k] 0-based to [icol, ilay] 0-based
            const int icol = (j-jmin)*nx + (i-imin) + offset;
            const int ilay = k;

            // Temperature heating rate for SW and LW
            q_arr(i,j,k,0) = sw_heating_tab(icol,ilay);
            q_arr(i,j,k,1) = lw_heating_tab(icol,ilay);

            // Convert the dT/dz to dTheta/dz
            Real iexner = 1./getExnergivenP(Real(p_lay_tab(icol,ilay)), R_d/Cp_d);
            q_arr(i,j,k,0) *= iexner;
            q_arr(i,j,k,1) *= iexner;

            // Populate the fluxes
            f_arr(i,j,k,0) = sw_flux_up_tab(icol,ilay);
            f_arr(i,j,k,1) = sw_flux_dn_tab(icol,ilay);
            f_arr(i,j,k,2) = lw_flux_up_tab(icol,ilay);
            f_arr(i,j,k,3) = lw_flux_dn_tab(icol,ilay);
        });
        if (m_lsm_fluxes) {
            const Array4<Real>& lsm_arr =  m_lsm_fluxes->array(mfi);
            ParallelFor(sbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // map [i,j,k] 0-based to [icol, ilay] 0-based
                const int icol = (j-jmin)*nx + (i-imin) + offset;

                // SW fluxes for LSM
                lsm_arr(i,j,k,0) = sfc_flux_dir_vis_tab(icol);
                lsm_arr(i,j,k,1) = sfc_flux_dir_nir_tab(icol);
                lsm_arr(i,j,k,2) = sfc_flux_dif_vis_tab(icol);
                lsm_arr(i,j,k,3) = sfc_flux_dif_nir_tab(icol);

                // Net SW flux for LSM
                lsm_arr(i,j,k,4) = sfc_flux_dir_vis_tab(icol) + sfc_flux_dir_nir_tab(icol)
                                 + sfc_flux_dif_vis_tab(icol) + sfc_flux_dif_nir_tab(icol);

                // LW flux for LSM (at bottom surface)
                lsm_arr(i,j,k,5) = lw_flux_dn_tab(icol,0);
            });
        }
        if (m_lsm_zenith) {
            const Array4<Real>& lsm_zenith_arr =  m_lsm_zenith->array(mfi);
            ParallelFor(sbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // map [i,j,k] 0-based to [icol, ilay] 0-based
                const int icol = (j-jmin)*nx + (i-imin) + offset;

                // export cosine zenith angle for LSM
                lsm_zenith_arr(i,j,k) = mu0_tab(icol);
            });
        }
        for (int ivar(0); ivar<lsm_output_ptrs.size(); ivar++) {
            if (lsm_output_ptrs[ivar]) {
                const Array4<Real>& lsm_out_arr = lsm_output_ptrs[ivar]->array(mfi);
                if (ivar==0) {
                    ParallelFor(sbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        // map [i,j,k] 0-based to [icol, ilay] 0-based
                        const int icol   = (j-jmin)*nx + (i-imin) + offset;

                        // export the desired variable at surface
                        lsm_out_arr(i,j,k) = mu0_tab(icol);
                    });
                } else {
                    auto rrtmgp_for_fill = rrtmgp_out_vars[ivar-1];
                    ParallelFor(sbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        // map [i,j,k] 0-based to [icol, ilay] 0-based
                        const int icol   = (j-jmin)*nx + (i-imin) + offset;

                        // export the desired variable at surface
                        lsm_out_arr(i,j,k) = rrtmgp_for_fill(icol,0);
                    });
                } // ivar
            } // valid ptr
        } // ivar
    }// mfi
}

void
Radiation::write_rrtmgp_fluxes ()
{
    Table2D<Real,Order::C> sw_flux_up_tab(sw_flux_up.data(), {0,0}, {static_cast<int>(sw_flux_up.extent(0)),static_cast<int>(sw_flux_up.extent(1))});
    Table2D<Real,Order::C> sw_flux_dn_tab(sw_flux_dn.data(), {0,0}, {static_cast<int>(sw_flux_dn.extent(0)),static_cast<int>(sw_flux_dn.extent(1))});
    Table2D<Real,Order::C> sw_flux_dn_dir_tab(sw_flux_dn_dir.data(), {0,0}, {static_cast<int>(sw_flux_dn_dir.extent(0)),static_cast<int>(sw_flux_dn_dir.extent(1))});
    Table2D<Real,Order::C> lw_flux_up_tab(lw_flux_up.data(), {0,0}, {static_cast<int>(lw_flux_up.extent(0)),static_cast<int>(lw_flux_up.extent(1))});
    Table2D<Real,Order::C> lw_flux_dn_tab(lw_flux_dn.data(), {0,0}, {static_cast<int>(lw_flux_dn.extent(0)),static_cast<int>(lw_flux_dn.extent(1))});

    int n_fluxes = 5;
    MultiFab mf_flux(m_cons_in->boxArray(), m_cons_in->DistributionMap(), n_fluxes, 0);

   for (MFIter mfi(mf_flux); mfi.isValid(); ++mfi) {
        const auto& vbx      = mfi.validbox();
        const int nx         = vbx.length(0);
        const int imin       = vbx.smallEnd(0);
        const int jmin       = vbx.smallEnd(1);
        const int offset     = m_col_offsets[mfi.index()];
        const Array4<Real>& dst_arr = mf_flux.array(mfi);
        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // map [i,j,k] 0-based to [icol, ilay] 0-based
            const int icol = (j-jmin)*nx + (i-imin) + offset;
            const int ilay = k;

            // SW and LW fluxes
            dst_arr(i,j,k,0) = sw_flux_up_tab(icol,ilay);
            dst_arr(i,j,k,1) = sw_flux_dn_tab(icol,ilay);
            dst_arr(i,j,k,2) = sw_flux_dn_dir_tab(icol,ilay);
            dst_arr(i,j,k,3) = lw_flux_up_tab(icol,ilay);
            dst_arr(i,j,k,4) = lw_flux_dn_tab(icol,ilay);
        });
   }


   std::string plotfilename = amrex::Concatenate("plt_rad", m_step, 5);
   Vector<std::string> flux_names = {"sw_flux_up", "sw_flux_dn", "sw_flux_dir",
                                     "lw_flux_up", "lw_flux_dn"};
   WriteSingleLevelPlotfile(plotfilename, mf_flux, flux_names, m_geom, m_time, m_step);
}

void Radiation::populateDatalogMF ()
{
    Table2D<Real,Order::C> sw_flux_up_tab(sw_flux_up.data(), {0,0}, {static_cast<int>(sw_flux_up.extent(0)),static_cast<int>(sw_flux_up.extent(1))});
    Table2D<Real,Order::C> sw_flux_dn_tab(sw_flux_dn.data(), {0,0}, {static_cast<int>(sw_flux_dn.extent(0)),static_cast<int>(sw_flux_dn.extent(1))});
    Table2D<Real,Order::C> sw_flux_dn_dir_tab(sw_flux_dn_dir.data(), {0,0}, {static_cast<int>(sw_flux_dn_dir.extent(0)),static_cast<int>(sw_flux_dn_dir.extent(1))});
    Table2D<Real,Order::C> lw_flux_up_tab(lw_flux_up.data(), {0,0}, {static_cast<int>(lw_flux_up.extent(0)),static_cast<int>(lw_flux_up.extent(1))});
    Table2D<Real,Order::C> lw_flux_dn_tab(lw_flux_dn.data(), {0,0}, {static_cast<int>(lw_flux_dn.extent(0)),static_cast<int>(lw_flux_dn.extent(1))});

    Table2D<Real,Order::C> sw_clrsky_flux_up_tab(sw_clrsky_flux_up.data(), {0,0},
                                                 {static_cast<int>(sw_clrsky_flux_up.extent(0)),static_cast<int>(sw_clrsky_flux_up.extent(1))});
    Table2D<Real,Order::C> sw_clrsky_flux_dn_tab(sw_clrsky_flux_dn.data(), {0,0},
                                                 {static_cast<int>(sw_clrsky_flux_dn.extent(0)),static_cast<int>(sw_clrsky_flux_dn.extent(1))});
    Table2D<Real,Order::C> sw_clrsky_flux_dn_dir_tab(sw_clrsky_flux_dn_dir.data(), {0,0},
                                                 {static_cast<int>(sw_clrsky_flux_dn_dir.extent(0)),static_cast<int>(sw_clrsky_flux_dn_dir.extent(1))});
    Table2D<Real,Order::C> lw_clrsky_flux_up_tab(lw_clrsky_flux_up.data(), {0,0},
                                                 {static_cast<int>(lw_clrsky_flux_up.extent(0)),static_cast<int>(lw_clrsky_flux_up.extent(1))});
    Table2D<Real,Order::C> lw_clrsky_flux_dn_tab(lw_clrsky_flux_dn.data(), {0,0},
                                                 {static_cast<int>(lw_clrsky_flux_dn.extent(0)),static_cast<int>(lw_clrsky_flux_dn.extent(1))});
    Table2D<Real,Order::C> sw_clrsky_heating_tab(sw_clrsky_heating.data(), {0,0},
                                                 {static_cast<int>(sw_clrsky_heating.extent(0)),static_cast<int>(sw_clrsky_heating.extent(1))});
    Table2D<Real,Order::C> lw_clrsky_heating_tab(lw_clrsky_heating.data(), {0,0},
                                                 {static_cast<int>(lw_clrsky_heating.extent(0)),static_cast<int>(lw_clrsky_heating.extent(1))});
    Table2D<Real,Order::C> sw_clnsky_flux_up_tab(sw_clnsky_flux_up.data(), {0,0},
                                                 {static_cast<int>(sw_clnsky_flux_up.extent(0)),static_cast<int>(sw_clnsky_flux_up.extent(1))});
    Table2D<Real,Order::C> sw_clnsky_flux_dn_tab(sw_clnsky_flux_dn.data(), {0,0},
                                                 {static_cast<int>(sw_clnsky_flux_dn.extent(0)),static_cast<int>(sw_clnsky_flux_dn.extent(1))});
    Table2D<Real,Order::C> sw_clnsky_flux_dn_dir_tab(sw_clnsky_flux_dn_dir.data(), {0,0},
                                                 {static_cast<int>(sw_clnsky_flux_dn_dir.extent(0)),static_cast<int>(sw_clnsky_flux_dn_dir.extent(1))});
    Table2D<Real,Order::C> lw_clnsky_flux_up_tab(lw_clnsky_flux_up.data(), {0,0},
                                                 {static_cast<int>(lw_clnsky_flux_up.extent(0)),static_cast<int>(lw_clnsky_flux_up.extent(1))});
    Table2D<Real,Order::C> lw_clnsky_flux_dn_tab(lw_clnsky_flux_dn.data(), {0,0},
                                                 {static_cast<int>(lw_clnsky_flux_dn.extent(0)),static_cast<int>(lw_clnsky_flux_dn.extent(1))});
    Table2D<Real,Order::C> sw_clnclrsky_flux_up_tab(sw_clnclrsky_flux_up.data(), {0,0},
                                                 {static_cast<int>(sw_clnclrsky_flux_up.extent(0)),static_cast<int>(sw_clnclrsky_flux_up.extent(1))});
    Table2D<Real,Order::C> sw_clnclrsky_flux_dn_tab(sw_clnclrsky_flux_dn.data(), {0,0},
                                                 {static_cast<int>(sw_clnclrsky_flux_dn.extent(0)),static_cast<int>(sw_clnclrsky_flux_dn.extent(1))});
    Table2D<Real,Order::C> sw_clnclrsky_flux_dn_dir_tab(sw_clnclrsky_flux_dn_dir.data(), {0,0},
                                                 {static_cast<int>(sw_clnclrsky_flux_dn_dir.extent(0)),static_cast<int>(sw_clnclrsky_flux_dn_dir.extent(1))});
    Table2D<Real,Order::C> lw_clnclrsky_flux_up_tab(lw_clnclrsky_flux_up.data(), {0,0},
                                                 {static_cast<int>(lw_clnclrsky_flux_up.extent(0)),static_cast<int>(lw_clnclrsky_flux_up.extent(1))});
    Table2D<Real,Order::C> lw_clnclrsky_flux_dn_tab(lw_clnclrsky_flux_dn.data(), {0,0},
                                                 {static_cast<int>(lw_clnclrsky_flux_dn.extent(0)),static_cast<int>(lw_clnclrsky_flux_dn.extent(1))});

    Table1D<Real>              mu0_tab(mu0.data(),              {0}, {static_cast<int>(mu0.extent(0))});

    auto extra_clnsky_diag = m_extra_clnsky_diag;
    auto extra_clnclrsky_diag = m_extra_clnclrsky_diag;

    for (MFIter mfi(datalog_mf); mfi.isValid(); ++mfi) {
        const auto& vbx      = mfi.validbox();
        const int nx         = vbx.length(0);
        const int imin       = vbx.smallEnd(0);
        const int jmin       = vbx.smallEnd(1);
        const int offset     = m_col_offsets[mfi.index()];
        const Array4<Real>& dst_arr = datalog_mf.array(mfi);
        const Array4<Real>&   q_arr = m_qheating_rates->array(mfi);
        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // map [i,j,k] 0-based to [icol, ilay] 0-based
            const int icol = (j-jmin)*nx + (i-imin) + offset;
            const int ilay = k;

            dst_arr(i,j,k,0) = q_arr(i, j, k, 0);
            dst_arr(i,j,k,1) = q_arr(i, j, k, 1);

            // SW and LW fluxes
            dst_arr(i,j,k,2) = sw_flux_up_tab(icol,ilay);
            dst_arr(i,j,k,3) = sw_flux_dn_tab(icol,ilay);
            dst_arr(i,j,k,4) = sw_flux_dn_dir_tab(icol,ilay);
            dst_arr(i,j,k,5) = lw_flux_up_tab(icol,ilay);
            dst_arr(i,j,k,6) = lw_flux_dn_tab(icol,ilay);

            // Cosine zenith angle
            dst_arr(i,j,k,7) = mu0_tab(icol);

            // Clear sky heating rates and fluxes:
            dst_arr(i,j,k,8) = sw_clrsky_heating_tab(icol, ilay);
            dst_arr(i,j,k,9) = lw_clrsky_heating_tab(icol, ilay);

            dst_arr(i,j,k,10) = sw_clrsky_flux_up_tab(icol,ilay);
            dst_arr(i,j,k,11) = sw_clrsky_flux_dn_tab(icol,ilay);
            dst_arr(i,j,k,12) = sw_clrsky_flux_dn_dir_tab(icol,ilay);
            dst_arr(i,j,k,13) = lw_clrsky_flux_up_tab(icol,ilay);
            dst_arr(i,j,k,14) = lw_clrsky_flux_dn_tab(icol,ilay);

            // Clean sky fluxes:
            if (extra_clnsky_diag) {
                dst_arr(i,j,k,15) = sw_clnsky_flux_up_tab(icol,ilay);
                dst_arr(i,j,k,16) = sw_clnsky_flux_dn_tab(icol,ilay);
                dst_arr(i,j,k,17) = sw_clnsky_flux_dn_dir_tab(icol,ilay);
                dst_arr(i,j,k,18) = lw_clnsky_flux_up_tab(icol,ilay);
                dst_arr(i,j,k,19) = lw_clnsky_flux_dn_tab(icol,ilay);
            }

            // Clean-clear sky fluxes:
            if (extra_clnclrsky_diag) {
                dst_arr(i,j,k,20) = sw_clnclrsky_flux_up_tab(icol,ilay);
                dst_arr(i,j,k,21) = sw_clnclrsky_flux_dn_tab(icol,ilay);
                dst_arr(i,j,k,22) = sw_clnclrsky_flux_dn_dir_tab(icol,ilay);
                dst_arr(i,j,k,23) = lw_clnclrsky_flux_up_tab(icol,ilay);
                dst_arr(i,j,k,24) = lw_clnclrsky_flux_dn_tab(icol,ilay);
            }
        });
   }
}

void Radiation::WriteDataLog (const Real &time)
{
    constexpr int datwidth = 14;
    constexpr int datprecision = 9;
    constexpr int timeprecision = 13;

    Gpu::HostVector<Real> h_avg_radqrsw, h_avg_radqrlw, h_avg_sw_up, h_avg_sw_dn, h_avg_sw_dn_dir, h_avg_lw_up, h_avg_lw_dn, h_avg_zenith;
    // Clear sky
    Gpu::HostVector<Real> h_avg_radqrcsw, h_avg_radqrclw, h_avg_sw_clr_up, h_avg_sw_clr_dn, h_avg_sw_clr_dn_dir, h_avg_lw_clr_up, h_avg_lw_clr_dn;
    // Clean sky
    Gpu::HostVector<Real> h_avg_sw_cln_up, h_avg_sw_cln_dn, h_avg_sw_cln_dn_dir, h_avg_lw_cln_up, h_avg_lw_cln_dn;
    // Clean clear sky
    Gpu::HostVector<Real> h_avg_sw_clnclr_up, h_avg_sw_clnclr_dn, h_avg_sw_clnclr_dn_dir, h_avg_lw_clnclr_up, h_avg_lw_clnclr_dn;


    auto domain = m_geom.Domain();
    h_avg_radqrsw    = sumToLine(datalog_mf, 0, 1, domain, 2);
    h_avg_radqrlw    = sumToLine(datalog_mf, 1, 1, domain, 2);
    h_avg_sw_up      = sumToLine(datalog_mf, 2, 1, domain, 2);
    h_avg_sw_dn      = sumToLine(datalog_mf, 3, 1, domain, 2);
    h_avg_sw_dn_dir  = sumToLine(datalog_mf, 4, 1, domain, 2);
    h_avg_lw_up      = sumToLine(datalog_mf, 5, 1, domain, 2);
    h_avg_lw_dn      = sumToLine(datalog_mf, 6, 1, domain, 2);
    h_avg_zenith     = sumToLine(datalog_mf, 7, 1, domain, 2);

    h_avg_radqrcsw       = sumToLine(datalog_mf, 8, 1, domain, 2);
    h_avg_radqrclw       = sumToLine(datalog_mf, 9, 1, domain, 2);
    h_avg_sw_clr_up      = sumToLine(datalog_mf, 10, 1, domain, 2);
    h_avg_sw_clr_dn      = sumToLine(datalog_mf, 11, 1, domain, 2);
    h_avg_sw_clr_dn_dir  = sumToLine(datalog_mf, 12, 1, domain, 2);
    h_avg_lw_clr_up      = sumToLine(datalog_mf, 13, 1, domain, 2);
    h_avg_lw_clr_dn      = sumToLine(datalog_mf, 14, 1, domain, 2);

    if (m_extra_clnsky_diag) {
        h_avg_sw_cln_up      = sumToLine(datalog_mf, 15, 1, domain, 2);
        h_avg_sw_cln_dn      = sumToLine(datalog_mf, 16, 1, domain, 2);
        h_avg_sw_cln_dn_dir  = sumToLine(datalog_mf, 17, 1, domain, 2);
        h_avg_lw_cln_up      = sumToLine(datalog_mf, 18, 1, domain, 2);
        h_avg_lw_cln_dn      = sumToLine(datalog_mf, 19, 1, domain, 2);
    }

    if (m_extra_clnclrsky_diag) {
        h_avg_sw_clnclr_up      = sumToLine(datalog_mf, 20, 1, domain, 2);
        h_avg_sw_clnclr_dn      = sumToLine(datalog_mf, 21, 1, domain, 2);
        h_avg_sw_clnclr_dn_dir  = sumToLine(datalog_mf, 22, 1, domain, 2);
        h_avg_lw_clnclr_up      = sumToLine(datalog_mf, 23, 1, domain, 2);
        h_avg_lw_clnclr_dn      = sumToLine(datalog_mf, 24, 1, domain, 2);
    }

    Real area_z = static_cast<Real>(domain.length(0)*domain.length(1));
    int nz = domain.length(2);
    for (int k = 0; k < nz; k++) {
        h_avg_radqrsw[k] /= area_z;
        h_avg_radqrlw[k] /= area_z;
        h_avg_sw_up[k] /= area_z;
        h_avg_sw_dn[k] /= area_z;
        h_avg_sw_dn_dir[k] /= area_z;
        h_avg_lw_up[k] /= area_z;
        h_avg_lw_dn[k] /= area_z;
        h_avg_zenith[k] /= area_z;

        h_avg_radqrcsw[k] /= area_z;
        h_avg_radqrclw[k] /= area_z;
        h_avg_sw_clr_up[k] /= area_z;
        h_avg_sw_clr_dn[k] /= area_z;
        h_avg_sw_clr_dn_dir[k] /= area_z;
        h_avg_lw_clr_up[k] /= area_z;
        h_avg_lw_clr_dn[k] /= area_z;
    }

    if (m_extra_clnsky_diag) {
        for (int k = 0; k < nz; k++) {
            h_avg_sw_cln_up[k] /= area_z;
            h_avg_sw_cln_dn[k] /= area_z;
            h_avg_sw_cln_dn_dir[k] /= area_z;
            h_avg_lw_cln_up[k] /= area_z;
            h_avg_lw_cln_dn[k] /= area_z;
        }
    }

    if (m_extra_clnclrsky_diag) {
        for (int k = 0; k < nz; k++) {
            h_avg_sw_clnclr_up[k] /= area_z;
            h_avg_sw_clnclr_dn[k] /= area_z;
            h_avg_sw_clnclr_dn_dir[k] /= area_z;
            h_avg_lw_clnclr_up[k] /= area_z;
            h_avg_lw_clnclr_dn[k] /= area_z;
        }
    }

    if (ParallelDescriptor::IOProcessor()) {
        std::ostream& log = *datalog;
        if (log.good()) {

            for (int k = 0; k < nz; k++)
            {
                Real z = k * m_geom.CellSize(2);
                log << std::setw(datwidth) << std::setprecision(timeprecision) << time << " "
                    << std::setw(datwidth) << std::setprecision(datprecision) << z << " "
                    << h_avg_radqrsw[k]   << " " << h_avg_radqrlw[k]       << " " << h_avg_sw_up[k] << " "
                    << h_avg_sw_dn[k]     << " " << h_avg_sw_dn_dir[k]     << " " << h_avg_lw_up[k] << " "
                    << h_avg_lw_dn[k]     << " " << h_avg_zenith[k]        << " "
                    << h_avg_radqrcsw[k]  << " " << h_avg_radqrclw[k]      << " " << h_avg_sw_clr_up[k] << " "
                    << h_avg_sw_clr_dn[k] << " " << h_avg_sw_clr_dn_dir[k] << " " << h_avg_lw_clr_up[k] << " "
                    << h_avg_lw_clr_dn[k] << " ";
                    if (m_extra_clnsky_diag) {
                        log << h_avg_sw_cln_up[k] << " " << h_avg_sw_cln_dn[k] << " " << h_avg_sw_cln_dn_dir[k] << " "
                            << h_avg_lw_cln_up[k] << " " << h_avg_lw_cln_dn[k] << " ";
                    } else {
                        log << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " ";
                    }

                    if (m_extra_clnclrsky_diag) {
                        log << h_avg_sw_clnclr_up[k] << " " << h_avg_sw_clnclr_dn[k] << " " << h_avg_sw_clnclr_dn_dir[k] << " "
                            << h_avg_lw_clnclr_up[k] << " " << h_avg_lw_clnclr_dn[k] << std::endl;
                    } else {
                        log << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;
                    }
            }
            // Write top face values
            Real z = nz * m_geom.CellSize(2);
            log << std::setw(datwidth) << std::setprecision(timeprecision) << time << " "
                << std::setw(datwidth) << std::setprecision(datprecision) << z << " "
                << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " "
                << 0.0 << " " << 0.0 << " "
                << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " "
                << 0.0 << " "
                << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " "
                << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0
                << std::endl;
        }
    }
}

void
Radiation::initialize_impl ()
{
    // Call API to initialize
    m_gas_concs.init(gas_names_offset, m_ncol, m_nlay);
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
    const auto nswbands = m_nswbands;

    // Compute orbital parameters; these are used both for computing
    // the solar zenith angle and also for computing total solar
    // irradiance scaling (tsi_scaling).
    double obliqr, lambm0, mvelpp;
    int  orbital_year = m_orbital_year;
    double eccen      = m_orbital_eccen;
    double obliq      = m_orbital_obliq;
    double mvelp      = m_orbital_mvelp;
    if (eccen >= 0 && obliq >= 0 && mvelp >= 0) {
      // fixed orbital parameters forced with orbital_year == ORB_UNDEF_INT
      orbital_year = ORB_UNDEF_INT;
    }
    orbital_params(orbital_year, eccen, obliq,
                   mvelp, obliqr, lambm0, mvelpp);

    // Use the orbital parameters to calculate the solar declination and eccentricity factor
    double delta, eccf;
    // Want day + fraction; calday 1 == Jan 1 0Z
    static constexpr double dpy[] = {0.0  ,  31.0,  59.0,  90.0, 120.0, 151.0,
                                     181.0, 212.0, 243.0, 273.0, 304.0, 334.0};
    bool leap = (m_orbital_year % 4 == 0 && (!(m_orbital_year % 100 == 0) || (m_orbital_year % 400 == 0))) ? true : false;
    double calday = dpy[m_orbital_mon-1] + (m_orbital_day-1.0) + m_orbital_sec/86400.0;
    // add extra day if leap year
    if (leap) { calday += 1.0; }
    orbital_decl(calday, eccen, mvelpp, lambm0, obliqr, delta, eccf);

    // Overwrite eccf if using a fixed solar constant.
    auto fixed_total_solar_irradiance = m_fixed_total_solar_irradiance;
    if (fixed_total_solar_irradiance >= 0){
       eccf = fixed_total_solar_irradiance/1360.9;
    }

    // Precompute volume mixing ratio (VMR) for all gases
    //
    // H2O is obtained from qv.
    // O3 may be a constant or a 1D vector
    // All other comps are set to constants for now
    for (int igas(0); igas < m_ngas; ++igas) {
        auto tmp2d = Kokkos::View<RealT**,layout_t,KokkosDefaultMem>("tmp2d", ncol, nlay);
        auto name = m_gas_names[igas];
        auto gas_mol_weight = m_mol_weight_gas[igas];
        if (name == "H2O") {
            auto qv_lay_d = qv_lay;
            Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, nlay}),
                                 KOKKOS_LAMBDA (int icol, int ilay)
            {
                tmp2d(icol,ilay) = qv_lay_d(icol,ilay) * mwdair/gas_mol_weight;
            });
        } else if (name == "CO2") {
            Kokkos::deep_copy(tmp2d, m_co2vmr);
        } else if (name == "O3")  {
            if (m_o3_size==1) {
                Kokkos::deep_copy(tmp2d, m_o3vmr[0] );
            } else {
                auto o3_lay_d = o3_lay;
                Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, nlay}),
                                     KOKKOS_LAMBDA (int icol, int ilay)
                {
                    tmp2d(icol,ilay) = o3_lay_d(ilay);
                });
            }
        } else if (name == "N2O") {
            Kokkos::deep_copy(tmp2d, m_n2ovmr);
        } else if (name == "CO")  {
            Kokkos::deep_copy(tmp2d, m_covmr );
        } else if (name == "CH4") {
            Kokkos::deep_copy(tmp2d, m_ch4vmr);
        } else if (name == "O2") {
            Kokkos::deep_copy(tmp2d, m_o2vmr );
        } else if (name == "N2") {
            Kokkos::deep_copy(tmp2d, m_n2vmr );
        } else {
            Abort("Radiation: Unknown gas component.");
        }

        // Populate GasConcs object
        m_gas_concs.set_vmr(name, tmp2d);
        Kokkos::fence();
    }

    // Populate mu0 1D array
    // This must be done on HOST and copied to device.
    auto h_mu0 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), mu0);
    if (m_fixed_solar_zenith_angle > 0) {
        Kokkos::deep_copy(h_mu0, m_fixed_solar_zenith_angle);
    } else {
        auto h_lat = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), lat);
        auto h_lon = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), lon);
        double dt  = double(m_dt);
        auto rad_freq_in_steps = m_rad_freq_in_steps;
        Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::Serial>(0, ncol),
                             [&] (int icol)
        {
            // Convert lat/lon to radians
            double lat_col = h_lat(icol)*PI/180.0;
            double lon_col = h_lon(icol)*PI/180.0;
            double lcalday = calday;
            double ldelta  = delta;
            h_mu0(icol)    = Real(orbital_cos_zenith(lcalday, lat_col, lon_col, ldelta, rad_freq_in_steps * dt));
        });
    }
    Kokkos::deep_copy(mu0, h_mu0);

    // Compute layer cloud mass per unit area (populates lwp/iwp)
    rrtmgp::mixing_ratio_to_cloud_mass(qc_lay, cldfrac_tot, r_lay, z_del, lwp);
    rrtmgp::mixing_ratio_to_cloud_mass(qi_lay, cldfrac_tot, r_lay, z_del, iwp);

    // Convert to g/m2 (needed by RRTMGP)
    Table2D<Real,Order::C> lwp_tab(lwp.data(), {0,0}, {static_cast<int>(lwp.extent(0)),static_cast<int>(lwp.extent(1))});
    Table2D<Real,Order::C> iwp_tab(iwp.data(), {0,0}, {static_cast<int>(iwp.extent(0)),static_cast<int>(iwp.extent(1))});
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, nlay}),
                         KOKKOS_LAMBDA (int icol, int ilay)
    {
        lwp_tab(icol,ilay) *= 1.e3;
        iwp_tab(icol,ilay) *= 1.e3;
    });

    // Expand surface_albedos along nswbands.
    // This is needed since rrtmgp require band-by-band.
    rrtmgp::compute_band_by_band_surface_albedos(ncol, nswbands,
                                                 sfc_alb_dir_vis, sfc_alb_dir_nir,
                                                 sfc_alb_dif_vis, sfc_alb_dif_nir,
                                                 sfc_alb_dir    , sfc_alb_dif);
    // Run RRTMGP driver
    rrtmgp::rrtmgp_main(ncol, m_nlay,
                        p_lay, t_lay,
                        p_lev, t_lev,
                        m_gas_concs,
                        sfc_alb_dir, sfc_alb_dif, mu0,
                        t_sfc, sfc_emis, lw_src,
                        lwp, iwp, eff_radius_qc, eff_radius_qi, cldfrac_tot,
                        aero_tau_sw, aero_ssa_sw, aero_g_sw, aero_tau_lw,
                        cld_tau_sw_bnd, cld_tau_lw_bnd,
                        cld_tau_sw_gpt, cld_tau_lw_gpt,
                        sw_flux_up, sw_flux_dn, sw_flux_dn_dir,
                        lw_flux_up, lw_flux_dn,
                        sw_clnclrsky_flux_up, sw_clnclrsky_flux_dn, sw_clnclrsky_flux_dn_dir,
                        sw_clrsky_flux_up, sw_clrsky_flux_dn, sw_clrsky_flux_dn_dir,
                        sw_clnsky_flux_up, sw_clnsky_flux_dn, sw_clnsky_flux_dn_dir,
                        lw_clnclrsky_flux_up, lw_clnclrsky_flux_dn,
                        lw_clrsky_flux_up, lw_clrsky_flux_dn,
                        lw_clnsky_flux_up, lw_clnsky_flux_dn,
                        sw_bnd_flux_up, sw_bnd_flux_dn, sw_bnd_flux_dir,
                        lw_bnd_flux_up, lw_bnd_flux_dn,
                        eccf, m_extra_clnclrsky_diag, m_extra_clnsky_diag);

#if 0
    // UNIT TEST
    //================================================================================
    Kokkos::deep_copy(mu0, 0.86);
    Kokkos::deep_copy(sfc_alb_dir_vis, 0.06);
    Kokkos::deep_copy(sfc_alb_dir_nir, 0.06);
    Kokkos::deep_copy(sfc_alb_dif_vis, 0.06);
    Kokkos::deep_copy(sfc_alb_dif_nir, 0.06);

    Kokkos::deep_copy(aero_tau_sw, 0.0);
    Kokkos::deep_copy(aero_ssa_sw, 0.0);
    Kokkos::deep_copy(aero_g_sw  , 0.0);
    Kokkos::deep_copy(aero_tau_lw, 0.0);

    // Generate some fake liquid and ice water data. We pick values to be midway between
    // the min and max of the valid lookup table values for effective radii
    real rel_val = 0.5 * (rrtmgp::cloud_optics_sw_k->get_min_radius_liq()
                        + rrtmgp::cloud_optics_sw_k->get_max_radius_liq());
    real rei_val = 0.5 * (rrtmgp::cloud_optics_sw_k->get_min_radius_ice()
                        + rrtmgp::cloud_optics_sw_k->get_max_radius_ice());

    // Restrict clouds to troposphere (> 100 hPa = 100*100 Pa) and not very close to the ground (< 900 hPa), and
    // put them in 2/3 of the columns since that's roughly the total cloudiness of earth.
    // Set sane values for liquid and ice water path.
    // NOTE: these "sane" values are in g/m2!
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, nlay}),
                         KOKKOS_LAMBDA (int icol, int ilay)
    {
        cldfrac_tot(icol,ilay) = (p_lay(icol,ilay) > 100. * 100.) &&
                                 (p_lay(icol,ilay) < 900. * 100.) &&
                                 (icol%3 != 0);
        // Ice and liquid will overlap in a few layers
        lwp(icol,ilay) = (cldfrac_tot(icol,ilay) && t_lay(icol,ilay) > 263.) ? 10. : 0.;
        iwp(icol,ilay) = (cldfrac_tot(icol,ilay) && t_lay(icol,ilay) < 273.) ? 10. : 0.;
        eff_radius_qc(icol,ilay) = (lwp(icol,ilay) > 0.) ? rel_val : 0.;
        eff_radius_qi(icol,ilay) = (iwp(icol,ilay) > 0.) ? rei_val : 0.;
    });

    rrtmgp::compute_band_by_band_surface_albedos(ncol, nswbands,
                                                 sfc_alb_dir_vis, sfc_alb_dir_nir,
                                                 sfc_alb_dif_vis, sfc_alb_dif_nir,
                                                 sfc_alb_dir    , sfc_alb_dif);

    rrtmgp::rrtmgp_main(ncol, m_nlay,
                        p_lay, t_lay,
                        p_lev, t_lev,
                        m_gas_concs,
                        sfc_alb_dir, sfc_alb_dif, mu0,
                        t_sfc, sfc_emis, lw_src,
                        lwp, iwp, eff_radius_qc, eff_radius_qi, cldfrac_tot,
                        aero_tau_sw, aero_ssa_sw, aero_g_sw, aero_tau_lw,
                        cld_tau_sw_bnd, cld_tau_lw_bnd,
                        cld_tau_sw_gpt, cld_tau_lw_gpt,
                        sw_flux_up, sw_flux_dn, sw_flux_dn_dir,
                        lw_flux_up, lw_flux_dn,
                        sw_clnclrsky_flux_up, sw_clnclrsky_flux_dn, sw_clnclrsky_flux_dn_dir,
                        sw_clrsky_flux_up, sw_clrsky_flux_dn, sw_clrsky_flux_dn_dir,
                        sw_clnsky_flux_up, sw_clnsky_flux_dn, sw_clnsky_flux_dn_dir,
                        lw_clnclrsky_flux_up, lw_clnclrsky_flux_dn,
                        lw_clrsky_flux_up, lw_clrsky_flux_dn,
                        lw_clnsky_flux_up, lw_clnsky_flux_dn,
                        sw_bnd_flux_up, sw_bnd_flux_dn, sw_bnd_flux_dir,
                        lw_bnd_flux_up, lw_bnd_flux_dn,
                        1.0, false, false);
    //================================================================================
#endif

    // Update heating tendency
    rrtmgp::compute_heating_rate(sw_flux_up, sw_flux_dn, r_lay, z_del, sw_heating);
    rrtmgp::compute_heating_rate(lw_flux_up, lw_flux_dn, r_lay, z_del, lw_heating);

    /*
    // AML DEBUG
    Kokkos::parallel_for(nlay+1, KOKKOS_LAMBDA (int ilay)
    {
        printf("Fluxes: %i %e %e %e %e %e\n",ilay,
                                             sw_flux_up(5,ilay), sw_flux_dn(5,ilay), sw_flux_dn_dir(5,ilay),
                                             lw_flux_up(5,ilay), lw_flux_dn(5,ilay));
    });
    Kokkos::parallel_for(nlay, KOKKOS_LAMBDA (int ilay)
    {
        printf("Heating Rate: %i %e %e\n",ilay,sw_heating(5,ilay),lw_heating(5,ilay));
    });
    */

    // Compute surface fluxes
    const int kbot = 0;
    Table3D<Real,Order::C> sw_bnd_flux_dif_tab(sw_bnd_flux_dif.data(), {0,0,0},
                           {static_cast<int>(sw_bnd_flux_dif.extent(0)),static_cast<int>(sw_bnd_flux_dif.extent(1)),static_cast<int>(sw_bnd_flux_dif.extent(2))});
    Table3D<Real,Order::C> sw_bnd_flux_dn_tab(sw_bnd_flux_dn.data(), {0,0,0},
                           {static_cast<int>(sw_bnd_flux_dn.extent(0)),static_cast<int>(sw_bnd_flux_dn.extent(1)),static_cast<int>(sw_bnd_flux_dn.extent(2))});
    Table3D<Real,Order::C> sw_bnd_flux_dir_tab(sw_bnd_flux_dir.data(), {0,0,0},
                           {static_cast<int>(sw_bnd_flux_dir.extent(0)),static_cast<int>(sw_bnd_flux_dir.extent(1)),static_cast<int>(sw_bnd_flux_dir.extent(2))});
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, nlay+1, nswbands}),
                         KOKKOS_LAMBDA (int icol, int ilay, int ibnd)
    {
        sw_bnd_flux_dif_tab(icol,ilay,ibnd) = sw_bnd_flux_dn_tab(icol,ilay,ibnd) - sw_bnd_flux_dir_tab(icol,ilay,ibnd);
    });
    rrtmgp::compute_broadband_surface_fluxes(ncol, kbot, nswbands,
                                             sw_bnd_flux_dir , sw_bnd_flux_dif ,
                                             sfc_flux_dir_vis, sfc_flux_dir_nir,
                                             sfc_flux_dif_vis, sfc_flux_dif_nir);
}


void
Radiation::finalize_impl (Vector<MultiFab*>& lsm_output_ptrs)
{
    // Finish rrtmgp
    m_gas_concs.reset();
    rrtmgp::rrtmgp_finalize();

    // Fill the AMReX MFs from Kokkos Views
    kokkos_buffers_to_mf(lsm_output_ptrs);

    // Write fluxes if requested
    if (m_rad_write_fluxes) { write_rrtmgp_fluxes(); }

    // Fill output data for datalog before deallocating
    if (datalog_int > 0) {
        rrtmgp::compute_heating_rate(sw_clrsky_flux_up, sw_clrsky_flux_dn, r_lay, z_del, sw_clrsky_heating);
        rrtmgp::compute_heating_rate(lw_clrsky_flux_up, lw_clrsky_flux_dn, r_lay, z_del, lw_clrsky_heating);

        populateDatalogMF();
    }

    // Deallocate the buffer arrays
    dealloc_buffers();
}
