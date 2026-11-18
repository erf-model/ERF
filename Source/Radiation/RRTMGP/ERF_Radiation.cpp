/*
 * RTE-RRTMGP radiation model interface to ERF
 * The original code is developed by RobertPincus, and the code is open source available at:
 *                        https://github.com/earth-system-radiation/rte-rrtmgp
 * Please reference to the following paper,
 *                        https://agupubs.onlinelibrary.wiley.com/doi/Real(10.1029)/2019MS001621
 * NOTE: we use the C++ version of RTE-RRTMGP, which is reimplemented the original Fortran
 * code using C++ KOKKOS for CUDA, HiP and SYCL application by E3SM ECP team, the C++ version
 * of the rte-rrtmgp code is located at:
 *                       https://github.com/E3SM-Project/rte-rrtmgp
 * The RTE-RRTMGP uses BSD-3-Clause Open Source License, if you want to make changes,
 * and modifications to the code, please refer to BSD-3-Clause Open Source License.
 */

#include <cmath>
#include <filesystem>
#include <limits>
#include "ERF_Constants.H"
#include <sstream>

#include "ERF_NCInterface.H"
#include "ERF_Radiation.H"
#include "ERF_RRTMGP_SurfaceTemperature.H"
#include "ERF_TerrainMetrics.H"

using namespace amrex;

namespace {
/**
 * Verify that the directory named by erf.rrtmgp_file_path exists and holds the four
 * netCDF lookup tables RRTMGP needs. Without this check, a missing directory or file
 * surfaces as a bare "NetCDF: No such file or directory" from deep inside the netCDF
 * layer, with no indication of which file was wanted or which input controls it.
 *
 * @param[in] path   value of erf.rrtmgp_file_path
 * @param[in] files  (input parameter name, file name) pairs for the required data files
 */
void
check_rrtmgp_data_files (const std::string& path,
                         const std::vector<std::pair<std::string,std::string>>& files)
{
    namespace fs = std::filesystem;

    // Report the absolute path as well, since a relative erf.rrtmgp_file_path is
    // resolved against the run directory and not against the inputs file.
    std::error_code ec;
    const fs::path dir(path);
    const fs::path abs_dir = fs::absolute(dir, ec);
    const std::string resolved = ec ? path : abs_dir.string();

    std::ostringstream problem;

    if (!fs::exists(dir, ec) || !fs::is_directory(dir, ec)) {
        problem << "The RRTMGP data directory does not exist:\n"
                << "    erf.rrtmgp_file_path = " << path << "\n"
                << "    (resolved to " << resolved << ")\n";
    } else {
        std::ostringstream missing;
        int nmissing = 0;
        for (const auto& [param, file] : files) {
            if (!fs::is_regular_file(dir / file, ec)) {
                missing << "    " << file << "  (set by erf." << param << ")\n";
                ++nmissing;
            }
        }
        if (nmissing == 0) { return; }
        problem << "The RRTMGP data directory\n"
                << "    erf.rrtmgp_file_path = " << path << "\n"
                << "    (resolved to " << resolved << ")\n"
                << "exists, but " << nmissing << " of the required data file(s) were not found there:\n"
                << missing.str();
    }

    Abort("RRTMGP radiation was requested, but its lookup data could not be found.\n"
          + problem.str()
          + "Set erf.rrtmgp_file_path in the inputs file to the directory that holds the\n"
            "RRTMGP netCDF lookup tables. The k-distribution files are shipped with the\n"
            "RRTMGP submodule in Submodules/RRTMGP/rrtmgp/data and the cloud optics files\n"
            "in Submodules/RRTMGP/extensions/cloud_optics (run\n"
            "'git submodule update --init --recursive' if those directories are empty);\n"
            "the full data package may also be downloaded from\n"
            "https://doi.org/10.22002/ppv8a-4q131 .\n");
}
} // namespace

Radiation::Radiation (const int& lev,
                      SolverChoice& sc)
{
    // Note that Kokkos is now initialized in main.cpp

    // Stored now, not only in set_grids, so it is valid from construction on
    m_lev = lev;

    m_rdOcp = sc.rdOcp;

    // Check if we have a valid moisture model
    if (sc.moisture_type != MoistureType::None) { m_moist = true; }

    // Cloud-ice support follows the configured moisture-component mapping.
    m_qi_comp = sc.moisture_indices.qi;
    m_ice = (m_qi_comp >= 0);

    // Check if we have a land surface model enabled
    if (sc.lsm_type != LandSurfaceType::None) { m_lsm = true; }

    // Construct parser object for following reads
    ParmParse pp("erf");

    // Must specify a surface temp (LSM can overwrite)
    pp.get("rad_t_sfc", m_rad_t_sfc);

    // SolverChoice reads the sun, site and update interval once for RRTMGP
    // and the optional topographic surface correction.
    const auto& solar = sc.rrtmgpSolarChoice;
    m_rad_freq_in_steps = solar.radiation_interval;
    m_lat_cons = solar.latitude;
    m_lon_cons = solar.longitude;
    m_fixed_solar_zenith_angle = solar.fixed_mu;

    // Get nvar if specified
    pp.queryAdd("rad_nvar", m_rad_nvar);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_rad_nvar >= 0,
        "erf.rad_nvar must be greater than 0. "
        "It controls the amount of memory allocated for temporaries with RRTMGP; "
        "a value of 0 would allocate no memory.");

    // Number of columns per RRTMGP chunk (controls peak GPU memory)
    pp.queryAdd("rad_ncol_chunk", m_ncol_chunk_requested);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_ncol_chunk_requested > 0,
        "erf.rad_ncol_chunk must be a positive integer (default 1024). "
        "It controls the number of columns processed per RRTMGP kernel launch; "
        "a value of 0 or negative would produce an infinite loop.");
    m_ncol_chunk = m_ncol_chunk_requested;

    // Flag to write fluxes to plt file
    pp.queryAdd("rad_write_fluxes", m_rad_write_fluxes);

    // Do MCICA subcolumn sampling
    pp.queryAdd("rad_do_subcol_sampling", m_do_subcol_sampling);

    // Determine orbital year. If orbital_year is negative, use current year
    // from timestamp for orbital year; if non-negative, use provided orbital year
    // for duration of simulation.  Note that this is keyed off the value itself,
    // not off whether the input was present, so that a negative value defers to
    // the timestamp as documented.
    pp.queryAdd("rad_orbital_year", m_orbital_year);
    m_fixed_orbital_year = (m_orbital_year >= 0);

    // Get orbital parameters from inputs file
    pp.queryAdd("rad_orbital_eccentricity", m_orbital_eccen);
    pp.queryAdd("rad_orbital_obliquity"   , m_orbital_obliq);
    pp.queryAdd("rad_orbital_mvelp"       , m_orbital_mvelp);

    // Value for prescribing an invariant solar constant (i.e. total solar irradiance at
    // TOA).  Used for idealized experiments such as RCE. Disabled when value is less than zero
    pp.queryAdd("fixed_total_solar_irradiance", m_fixed_total_solar_irradiance);

    // The solar-input range checks live in RRTMGPSolarChoice::read. Surface
    // temperature is still read and validated by this interface.
    if (!std::isfinite(m_rad_t_sfc) || m_rad_t_sfc <= Real(0.0)) {
        amrex::Abort("erf.rad_t_sfc = " + std::to_string(m_rad_t_sfc) +
                     " must be a positive temperature [K].");
    }

    // Get prescribed surface values of greenhouse gases
    pp.queryAdd("co2vmr", m_co2vmr);
    // The O3 default is a single near-surface value applied to every layer, which
    // is far too low in the stratosphere; warn so deep runs do not use it silently
    if (!pp.queryAdd("o3vmr", m_o3vmr) && lev == 0) {
        // amrex::Warning prints on every rank, so print once from the I/O rank instead
        amrex::Print() << "WARNING: erf.o3vmr is not set: using " << m_o3vmr[0]
                       << " at every layer. This near-surface value underestimates ozone"
                          " (and shortwave heating) aloft; set erf.o3vmr to a profile for"
                          " domains that reach the stratosphere.\n";
    }
    pp.queryAdd("n2ovmr", m_n2ovmr);
    pp.queryAdd("covmr" , m_covmr );
    pp.queryAdd("ch4vmr", m_ch4vmr);
    pp.queryAdd("o2vmr" , m_o2vmr );
    pp.queryAdd("n2vmr" , m_n2vmr );

    // Aerosol forcing hook (not implemented). The aerosol arrays that used to be
    // passed through rrtmgp_main were never populated with real data, so enabling
    // this flag only ever multiplied radiation by zero aerosol optics. The hook is
    // kept so a future SPA/prescribed-aerosol scheme can wire in without touching
    // the ParmParse surface.
    pp.queryAdd("rad_do_aerosol", m_do_aerosol_rad);
    if (m_do_aerosol_rad) {
        amrex::Abort("erf.rad_do_aerosol = true is not supported: aerosol forcing is "
                     "currently not implemented in the ERF RRTMGP interface. The hook "
                     "is retained for a future aerosol coupling; set rad_do_aerosol = "
                     "false (or remove it) to continue.");
    }

    // Whether we do extra clean/clear sky calculations
    pp.queryAdd("rad_extra_clnclrsky_diag", m_extra_clnclrsky_diag);
    pp.queryAdd("rad_extra_clnsky_diag"   , m_extra_clnsky_diag);

    // Parse path and file names
    pp.queryAdd("rrtmgp_file_path"      , rrtmgp_file_path);
    pp.queryAdd("rrtmgp_coeffs_sw"      , rrtmgp_coeffs_sw  );
    pp.queryAdd("rrtmgp_coeffs_lw"      , rrtmgp_coeffs_lw  );
    pp.queryAdd("rrtmgp_cloud_optics_sw", rrtmgp_cloud_optics_sw);
    pp.queryAdd("rrtmgp_cloud_optics_lw", rrtmgp_cloud_optics_lw);

    // Fail early, and with a message that names the missing files and the inputs that
    // control them, rather than letting netCDF report a bare "No such file or directory"
    // when the coefficients are opened just below (or, for the cloud optics files, much
    // later inside rrtmgp_initialize). Every rank does this so that all of them abort.
    check_rrtmgp_data_files(rrtmgp_file_path,
                            {{"rrtmgp_coeffs_sw"      , rrtmgp_coeffs_sw      },
                             {"rrtmgp_coeffs_lw"      , rrtmgp_coeffs_lw      },
                             {"rrtmgp_cloud_optics_sw", rrtmgp_cloud_optics_sw},
                             {"rrtmgp_cloud_optics_lw", rrtmgp_cloud_optics_lw}});

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
                      double& time,
                      const double& dt,
                      const BoxArray& /*ba*/,
                      Geometry& geom,
                      MultiFab* cons_in,
                      iMultiFab* lmask,
                      MultiFab*  t_surf,
                      const Vector<const MultiFab*>& lsm_input_ptrs,
                      MultiFab* qheating_rates,
                      MultiFab* rad_fluxes,
                      MultiFab* z_phys,
                      MultiFab* z_sfc,
                      MultiFab* lat,
                      MultiFab* lon,
                      const bool force_update)

{
    // Set data members that may change
    m_lev            = level;
    m_step           = step;
    m_time           = time;
    m_dt             = dt;
    m_geom           = geom;
    m_cons_in        = cons_in;
    m_qheating_rates = qheating_rates;
    m_rad_fluxes     = rad_fluxes;
    m_z_phys         = z_phys;
    m_z_sfc          = z_sfc;
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
    m_update_rad = force_update;
    if (m_rad_freq_in_steps > 0) {
        m_update_rad = m_update_rad || (m_step == 0) || (m_step % m_rad_freq_in_steps == 0);
    }

    if (m_update_rad) {
        // Call to Init() has set the dimensions: ncol & nlay

        // Allocate the buffer arrays
        alloc_buffers();

        // Fill the KOKKOS Views from AMReX MFs
        mf_to_kokkos_buffers(lmask, t_surf, lsm_input_ptrs);

        // (Re)define the datalog MF whenever the grids change; this must always
        // match the layout of cons_in since populateDatalogMF() iterates over it
        // while indexing m_col_offsets and m_qheating_rates.
        if (datalog_int > 0) {
            bool needs_define = ( (datalog_mf.boxArray()        != cons_in->boxArray()) ||
                                  (datalog_mf.DistributionMap() != cons_in->DistributionMap()) ||
                                  (datalog_mf.nComp()           != nDatalogComp()) ||
                                  (datalog_bdy_mf.nComp()       != nDatalogBdyComp()) );
            if (needs_define) {
                datalog_mf.define(cons_in->boxArray(), cons_in->DistributionMap(), nDatalogComp(), 0);
                datalog_mf.setVal(0.0);

                // Every box spans the full column (asserted in Init()), so collapsing each one to a
                // single cell in z keeps the same box ordering and lets populateDatalogMF() share one MFIter
                BoxList bl_bdy = cons_in->boxArray().boxList();
                for (auto& b : bl_bdy) { b.setRange(2, 0, 1); }
                datalog_bdy_mf.define(BoxArray(std::move(bl_bdy)), cons_in->DistributionMap(),
                                      nDatalogBdyComp(), 0);
                datalog_bdy_mf.setVal(0.0);
            }
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
                                     "erf.o3vmr must have 1 or nlay values");
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
    if (datalog_int > 0) {
        sw_clrsky_heating = real2d_k("sw_clrsky_heating", m_ncol, m_nlay);
        lw_clrsky_heating = real2d_k("lw_clrsky_heating", m_ncol, m_nlay);
    }

    // 2d size (ncol, nlay+1)
    d_tint                   = real2d_k("d_tint"               , m_ncol, m_nlay+1);
    p_lev                    = real2d_k("p_lev"                , m_ncol, m_nlay+1);
    t_lev                    = real2d_k("t_lev"                , m_ncol, m_nlay+1);

    sw_flux_up               = real2d_k("sw_flux_up"           , m_ncol, m_nlay+1);
    sw_flux_dn               = real2d_k("sw_flux_dn"           , m_ncol, m_nlay+1);
    sw_flux_dn_dir           = real2d_k("sw_flux_dn_dir"       , m_ncol, m_nlay+1);

    lw_flux_up               = real2d_k("lw_flux_up"           , m_ncol, m_nlay+1);
    lw_flux_dn               = real2d_k("lw_flux_dn"           , m_ncol, m_nlay+1);

    // Clear-sky flux arrays are always needed
    if (datalog_int > 0) {
        sw_clrsky_flux_up        = real2d_k("sw_clrsky_flux_up"    , m_ncol, m_nlay+1);
        sw_clrsky_flux_dn        = real2d_k("sw_clrsky_flux_dn"    , m_ncol, m_nlay+1);
        sw_clrsky_flux_dn_dir    = real2d_k("sw_clrsky_flux_dn_dir", m_ncol, m_nlay+1);
        lw_clrsky_flux_up        = real2d_k("lw_clrsky_flux_up"    , m_ncol, m_nlay+1);
        lw_clrsky_flux_dn        = real2d_k("lw_clrsky_flux_dn"    , m_ncol, m_nlay+1);
    } else {
        // Use m_ncol_chunk_requested to prevent pool shrinkage after regrid
        sw_clrsky_flux_up        = real2d_k("sw_clrsky_flux_up"    , m_ncol_chunk_requested, m_nlay+1);
        sw_clrsky_flux_dn        = real2d_k("sw_clrsky_flux_dn"    , m_ncol_chunk_requested, m_nlay+1);
        sw_clrsky_flux_dn_dir    = real2d_k("sw_clrsky_flux_dn_dir", m_ncol_chunk_requested, m_nlay+1);
        lw_clrsky_flux_up        = real2d_k("lw_clrsky_flux_up"    , m_ncol_chunk_requested, m_nlay+1);
        lw_clrsky_flux_dn        = real2d_k("lw_clrsky_flux_dn"    , m_ncol_chunk_requested, m_nlay+1);
    }

    // Clean-clear-sky diagnostic fluxes (only when enabled)
    if (m_extra_clnclrsky_diag) {
        sw_clnclrsky_flux_up     = real2d_k("sw_clnclrsky_flux_up"    , m_ncol, m_nlay+1);
        sw_clnclrsky_flux_dn     = real2d_k("sw_clnclrsky_flux_dn"    , m_ncol, m_nlay+1);
        sw_clnclrsky_flux_dn_dir = real2d_k("sw_clnclrsky_flux_dn_dir", m_ncol, m_nlay+1);
        lw_clnclrsky_flux_up     = real2d_k("lw_clnclrsky_flux_up"    , m_ncol, m_nlay+1);
        lw_clnclrsky_flux_dn     = real2d_k("lw_clnclrsky_flux_dn"    , m_ncol, m_nlay+1);
    } else {
        sw_clnclrsky_flux_up     = real2d_k("sw_clnclrsky_flux_up"    , 1, 1);
        sw_clnclrsky_flux_dn     = real2d_k("sw_clnclrsky_flux_dn"    , 1, 1);
        sw_clnclrsky_flux_dn_dir = real2d_k("sw_clnclrsky_flux_dn_dir", 1, 1);
        lw_clnclrsky_flux_up     = real2d_k("lw_clnclrsky_flux_up"    , 1, 1);
        lw_clnclrsky_flux_dn     = real2d_k("lw_clnclrsky_flux_dn"    , 1, 1);
    }

    // Clean-sky diagnostic fluxes (only when enabled)
    if (m_extra_clnsky_diag) {
        sw_clnsky_flux_up        = real2d_k("sw_clnsky_flux_up"       , m_ncol, m_nlay+1);
        sw_clnsky_flux_dn        = real2d_k("sw_clnsky_flux_dn"       , m_ncol, m_nlay+1);
        sw_clnsky_flux_dn_dir    = real2d_k("sw_clnsky_flux_dn_dir"   , m_ncol, m_nlay+1);
        lw_clnsky_flux_up        = real2d_k("lw_clnsky_flux_up"       , m_ncol, m_nlay+1);
        lw_clnsky_flux_dn        = real2d_k("lw_clnsky_flux_dn"       , m_ncol, m_nlay+1);
    } else {
        sw_clnsky_flux_up        = real2d_k("sw_clnsky_flux_up"       , 1, 1);
        sw_clnsky_flux_dn        = real2d_k("sw_clnsky_flux_dn"       , 1, 1);
        sw_clnsky_flux_dn_dir    = real2d_k("sw_clnsky_flux_dn_dir"   , 1, 1);
        lw_clnsky_flux_up        = real2d_k("lw_clnsky_flux_up"       , 1, 1);
        lw_clnsky_flux_dn        = real2d_k("lw_clnsky_flux_dn"       , 1, 1);
    }

    // 3d size (ncol_chunk, nlay+1, nswbands)
    // Use m_ncol_chunk_requested to prevent pool shrinkage after regrid
    sw_bnd_flux_up  = real3d_k("sw_bnd_flux_up" , m_ncol_chunk_requested, m_nlay+1, m_nswbands);
    sw_bnd_flux_dn  = real3d_k("sw_bnd_flux_dn" , m_ncol_chunk_requested, m_nlay+1, m_nswbands);
    sw_bnd_flux_dir = real3d_k("sw_bnd_flux_dir", m_ncol_chunk_requested, m_nlay+1, m_nswbands);
    sw_bnd_flux_dif = real3d_k("sw_bnd_flux_dif", m_ncol_chunk_requested, m_nlay+1, m_nswbands);

    // 3d size (ncol_chunk, nlay+1, nlwbands)
    // Use m_ncol_chunk_requested to prevent pool shrinkage after regrid
    lw_bnd_flux_up = real3d_k("lw_bnd_flux_up" , m_ncol_chunk_requested, m_nlay+1, m_nlwbands);
    lw_bnd_flux_dn = real3d_k("lw_bnd_flux_dn" , m_ncol_chunk_requested, m_nlay+1, m_nlwbands);

    // 2d size (ncol, nswbands)
    sfc_alb_dir = real2d_k("sfc_alb_dir", m_ncol, m_nswbands);
    sfc_alb_dif = real2d_k("sfc_alb_dif", m_ncol, m_nswbands);

    // Aerosol optical properties — allocated only when aerosol coupling is on.
    // The flag gates allocation so today (coupling not implemented, abort fires
    // in the constructor) these stay as empty Views and cost nothing. When a
    // future aerosol scheme populates them, hook up the plumbing into
    // rrtmgp_main as well.
    if (m_do_aerosol_rad) {
        aero_tau_sw = real3d_k("aero_tau_sw", m_ncol, m_nlay, m_nswbands);
        aero_ssa_sw = real3d_k("aero_ssa_sw", m_ncol, m_nlay, m_nswbands);
        aero_g_sw   = real3d_k("aero_g_sw",   m_ncol, m_nlay, m_nswbands);
        aero_tau_lw = real3d_k("aero_tau_lw", m_ncol, m_nlay, m_nlwbands);
    }
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

    // Aerosol scaffolding (no-op unless m_do_aerosol_rad enabled allocation above)
    aero_tau_sw = real3d_k();
    aero_ssa_sw = real3d_k();
    aero_g_sw   = real3d_k();
    aero_tau_lw = real3d_k();
}


void
Radiation::mf_to_kokkos_buffers (iMultiFab* lmask,
                                 MultiFab*  t_surf,
                                 const Vector<const MultiFab*>& lsm_input_ptrs)
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
    const int qi_comp = m_qi_comp;
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
    Real rdOcp = m_rdOcp;

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
            Real qv = (moist) ? std::max(cons_arr(i,j,k,RhoQ1_comp)/r,Real(0.)) : Real(0.);
            Real qc = (moist) ? std::max(cons_arr(i,j,k,RhoQ2_comp)/r,Real(0.)) : Real(0.);
            Real qi = (ice && qi_comp >= 0) ? std::max(cons_arr(i,j,k,qi_comp)/r,Real(0.)) : Real(0.);

            // EOS avg to z-face
            Real r_lo   = cons_arr(i,j,k-1,Rho_comp);
            Real rt_lo  = cons_arr(i,j,k-1,RhoTheta_comp);
            Real qv_lo  = (moist) ? cons_arr(i,j,k-1,RhoQ1_comp)/r_lo : Real(0.);
            Real dz_k   = (z_arr) ? Real(0.125) * ( (z_arr(i  ,j  ,k+1) - z_arr(i  ,j  ,k))
                                                  + (z_arr(i+1,j  ,k+1) - z_arr(i+1,j  ,k))
                                                  + (z_arr(i  ,j+1,k+1) - z_arr(i  ,j+1,k))
                                                  + (z_arr(i+1,j+1,k+1) - z_arr(i+1,j+1,k)) ) : Real(0.5)*dz; // Dist from w-face to CC at k
            Real dz_km1 = (z_arr) ? Real(0.125) * ( (z_arr(i  ,j  ,k  ) - z_arr(i  ,j  ,k-1))
                                                  + (z_arr(i+1,j  ,k  ) - z_arr(i+1,j  ,k-1))
                                                  + (z_arr(i  ,j+1,k  ) - z_arr(i  ,j+1,k-1))
                                                  + (z_arr(i+1,j+1,k  ) - z_arr(i+1,j+1,k-1)) ) : Real(0.5)*dz; // Dist from w-face to CC at k-1
            // NOTE: Linear interpolation to the w-face weights each CC value by the
            //       distance from the face to the *opposite* CC (inverse distance)
            Real r_avg  = (dz_km1*r  + dz_k*r_lo ) / (dz_k + dz_km1);
            Real rt_avg = (dz_km1*rt + dz_k*rt_lo) / (dz_k + dz_km1);
            Real qv_avg = (dz_km1*qv + dz_k*qv_lo) / (dz_k + dz_km1);

            // Views at CC
            r_lay_tab(icol,ilay) = r;

            p_lay_tab(icol,ilay) = getPgivenRTh(rt, qv);
            t_lay_tab(icol,ilay) = getTgivenRandRTh(r, rt, qv);
            z_del_tab(icol,ilay) = (z_arr) ? Real(0.25) * ( (z_arr(i  ,j  ,k+1) - z_arr(i  ,j  ,k))
                                                          + (z_arr(i+1,j  ,k+1) - z_arr(i+1,j  ,k))
                                                          + (z_arr(i  ,j+1,k+1) - z_arr(i  ,j+1,k))
                                                          + (z_arr(i+1,j+1,k+1) - z_arr(i+1,j+1,k)) ) : dz;
            qv_lay_tab(icol,ilay) = qv;
            qc_lay_tab(icol,ilay) = qc;
            qi_lay_tab(icol,ilay) = qi;
            // NOTE: Binary (all-or-nothing) cloud fraction, consistent with saturation adjustment and suited to
            //       resolved clouds; with it erf.rad_do_subcol_sampling has no effect (see #5212)
            cldfrac_tot_tab(icol,ilay) = ((qc+qi)>Real(0.)) ? Real(1.) : Real(0.);

            // NOTE: These are populated in 'mixing_ratio_to_cloud_mass'
            lwp_tab(icol,ilay) = Real(0.);
            iwp_tab(icol,ilay) = Real(0.);

            // NOTE: These would be populated from P3 (we use the constants in p3_main_impl.hpp)
            // NOTE: These are in units of micron!
            eff_radius_qc_tab(icol,ilay) = (qc>Real(0.)) ? Real(10.0) : Real(0.);
            eff_radius_qi_tab(icol,ilay) = (qi>Real(0.)) ? Real(25.0) : Real(0.);

            // Buffers on z-faces (nlay+1)
            p_lev_tab(icol,ilay) = getPgivenRTh(rt_avg, qv_avg);
            t_lev_tab(icol,ilay) = getTgivenRandRTh(r_avg, rt_avg, qv_avg);
            if (ilay==(nlay-1)) {
                Real r_hi  = cons_arr(i,j,k+1,Rho_comp);
                Real rt_hi = cons_arr(i,j,k+1,RhoTheta_comp);
                Real qv_hi = (moist) ? std::max(cons_arr(i,j,k+1,RhoQ1_comp)/r_hi,Real(0.)) : Real(0.);
                Real dz_kp1 = (z_arr) ? Real(0.125) * ( (z_arr(i  ,j  ,k+2) - z_arr(i  ,j  ,k+1))
                                                      + (z_arr(i+1,j  ,k+2) - z_arr(i+1,j  ,k+1))
                                                      + (z_arr(i  ,j+1,k+2) - z_arr(i  ,j+1,k+1))
                                                      + (z_arr(i+1,j+1,k+2) - z_arr(i+1,j+1,k+1)) ) : Real(0.5)*dz; // Dist from w-face to CC at k+1
                r_avg  = (dz_kp1*r  + dz_k*r_hi ) / (dz_k + dz_kp1);
                rt_avg = (dz_kp1*rt + dz_k*rt_hi) / (dz_k + dz_kp1);
                qv_avg = (dz_kp1*qv + dz_k*qv_hi) / (dz_k + dz_kp1);
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
        Kokkos::deep_copy(sfc_alb_dir_vis, Real(0.06));
        Kokkos::deep_copy(sfc_alb_dir_nir, Real(0.06));
        Kokkos::deep_copy(sfc_alb_dif_vis, Real(0.06));
        Kokkos::deep_copy(sfc_alb_dif_nir, Real(0.06));

        // AML NOTE: These are not used in current EAMXX, I've left
        //           the code to plug into these if we need it.
        //
        // Current EAMXX constants
        Kokkos::deep_copy(sfc_emis, Real(0.98));
        Kokkos::deep_copy(lw_src  , zero );
    } else {
        Vector<real1d_k> rrtmgp_in_vars = {t_sfc, sfc_emis,
                                           sfc_alb_dir_vis, sfc_alb_dir_nir,
                                           sfc_alb_dif_vis, sfc_alb_dif_nir};
        Vector<Real> rrtmgp_default_vals = {rad_t_sfc, Real(0.98),
                                            Real(0.06), Real(0.06),
                                            Real(0.06), Real(0.06)};
        for (int ivar(0); ivar<lsm_input_ptrs.size(); ivar++) {
            auto rrtmgp_default_val = rrtmgp_default_vals[ivar];
            auto rrtmgp_to_fill_k   = rrtmgp_in_vars[ivar];
            amrex::Table1D<amrex::Real> rrtmgp_to_fill(rrtmgp_to_fill_k.data(),
                                                       0, rrtmgp_to_fill_k.extent(0));
            for (MFIter mfi(*m_cons_in); mfi.isValid(); ++mfi) {
                const auto& vbx  = mfi.validbox();
                const auto& sbx  = makeSlab(vbx,2,vbx.smallEnd(2));
                const int nx     = vbx.length(0);
                const int imin   = vbx.smallEnd(0);
                const int jmin   = vbx.smallEnd(1);
                const int offset = m_col_offsets[mfi.index()];
                const int k_surface = vbx.smallEnd(2);
                const Array4<const Real>& cons_arr = m_cons_in->const_array(mfi);
                const Array4<const Real>& z_arr = (m_z_phys) ? m_z_phys->const_array(mfi) :
                                                               Array4<const Real>{};
                const Array4<const Real>& zsfc_arr = (m_z_phys) ? m_z_sfc->const_array(mfi) :
                                                                  Array4<const Real>{};
                const Array4<const int>& lmask_arr   = (lmask)   ? lmask->const_array(mfi) :
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

                    // Surface temperature has a distinct contract: LSM and
                    // RRTMGP use absolute temperature, while SurfaceLayer
                    // supplies potential temperature.  The other LSM inputs
                    // retain their existing validity/default handling.
                    if (ivar == 0) {
                        const bool has_lsm_t_sfc = static_cast<bool>(lsm_in_arr);
                        const bool valid_lsm_t_sfc =
                            has_lsm_t_sfc && amrex::Math::isfinite(lsm_in_arr(i,j,k)) &&
                            (lsm_in_arr(i,j,k) > Real(0.)) &&
                            (lsm_in_arr(i,j,k) < lsm_undefined);
                        // Match TwoStream: convert SurfaceLayer theta with
                        // the physical surface pressure diagnosed from the
                        // lowest atmospheric cell.
                        rrtmgp::resolve_surface_temperature(
                            is_land,
                            has_lsm_t_sfc,
                            valid_lsm_t_sfc,
                            has_lsm_t_sfc ? lsm_in_arr(i,j,k) : Real(0.),
                            static_cast<bool>(tsurf_arr),
                            tsurf_arr ? tsurf_arr(i,j,k) : Real(0.),
                            erf_surface_temperature::pressure_at_surface(
                                cons_arr(i,j,k_surface,Rho_comp),
                                cons_arr(i,j,k_surface,RhoTheta_comp),
                                moist ? std::max(cons_arr(i,j,k_surface,RhoQ1_comp) /
                                                 cons_arr(i,j,k_surface,Rho_comp), Real(0.)) : Real(0.),
                                z_arr ? Compute_Zrel_AtCellCenter(i,j,k_surface,zsfc_arr,z_arr) : Real(0.5)*dz),
                            rdOcp,
                            rrtmgp_default_val,
                            rrtmgp_to_fill(icol));
                    } else {
                        // Have LSM and are over land.
                        const bool valid_lsm_data =
                            (lsm_in_arr && (lsm_in_arr(i,j,k) < lsm_undefined));
                        if (is_land && valid_lsm_data) {
                            rrtmgp_to_fill(icol) = lsm_in_arr(i,j,k);
                        } else {
                            // Use the default value.
                            rrtmgp_to_fill(icol) = rrtmgp_default_val;
                        }
                    }
                });
            } //mfi
        } // ivar
        Kokkos::deep_copy(lw_src, zero );
    } // have lsm

    // Enforce consistency between t_sfc and t_lev at bottom surface
    Kokkos::parallel_for(Kokkos::RangePolicy(0, ncol),
                         KOKKOS_LAMBDA (int icol)
    {
        t_lev_tab(icol,0) = t_sfc_tab(icol);
    });
}


void
Radiation::kokkos_buffers_to_mf (const Vector<MultiFab*>& lsm_output_ptrs)
{
    // The heating rates and fluxes read below were written by Kokkos kernels in run_impl().
    // Kokkos launches on its own default instance while the ParallelFors here launch on
    // amrex::Gpu::gpuStream(), and nothing orders the consumer against the producer: MFIter
    // synchronizes the streams it used, but only after its own kernels, which is too late.
    // This is the same cross-stream race that made radqrclw vary run to run in the datalog
    // path (see section 7 of RRTMGP_Memory_Reduction.md); that fence guarded one call site,
    // this one guards the hand-off that feeds qheating_rates and rad_fluxes back into the
    // solution.  No-op on a CPU build.
    Kokkos::fence();

    // Heating rate, fluxes, zenith, lsm ptrs

    Table2D<Real,Order::C> p_lay_tab(p_lay.data(), {0,0}, {static_cast<int>(p_lay.extent(0)),static_cast<int>(p_lay.extent(1))});
    Table2D<Real,Order::C> sw_heating_tab(sw_heating.data(), {0,0}, {static_cast<int>(sw_heating.extent(0)),static_cast<int>(sw_heating.extent(1))});
    Table2D<Real,Order::C> lw_heating_tab(lw_heating.data(), {0,0}, {static_cast<int>(lw_heating.extent(0)),static_cast<int>(lw_heating.extent(1))});
    Table2D<Real,Order::C> sw_flux_up_tab(sw_flux_up.data(), {0,0}, {static_cast<int>(sw_flux_up.extent(0)),static_cast<int>(sw_flux_up.extent(1))});
    Table2D<Real,Order::C> sw_flux_dn_tab(sw_flux_dn.data(), {0,0}, {static_cast<int>(sw_flux_dn.extent(0)),static_cast<int>(sw_flux_dn.extent(1))});
    Table2D<Real,Order::C> lw_flux_up_tab(lw_flux_up.data(), {0,0}, {static_cast<int>(lw_flux_up.extent(0)),static_cast<int>(lw_flux_up.extent(1))});
    Table2D<Real,Order::C> lw_flux_dn_tab(lw_flux_dn.data(), {0,0}, {static_cast<int>(lw_flux_dn.extent(0)),static_cast<int>(lw_flux_dn.extent(1))});

    TableData<Real,1> sfc_flux_sw_dn; sfc_flux_sw_dn.resize({0}, {static_cast<int>(sw_flux_dn.extent(0))});
    TableData<Real,1> sfc_flux_lw_dn; sfc_flux_lw_dn.resize({0}, {static_cast<int>(lw_flux_dn.extent(0))});
    Table1D<Real> sfc_flux_sw_dn_tab = sfc_flux_sw_dn.table();
    Table1D<Real> sfc_flux_lw_dn_tab = sfc_flux_lw_dn.table();
    Table1D<Real> sfc_flux_sw_dir_vis_tab(sfc_flux_dir_vis.data(), {0}, {static_cast<int>(sfc_flux_dir_vis.extent(0))});
    Table1D<Real> sfc_flux_sw_dir_nir_tab(sfc_flux_dir_nir.data(), {0}, {static_cast<int>(sfc_flux_dir_nir.extent(0))});
    Table1D<Real> sfc_flux_sw_dif_vis_tab(sfc_flux_dif_vis.data(), {0}, {static_cast<int>(sfc_flux_dif_vis.extent(0))});
    Table1D<Real> sfc_flux_sw_dif_nir_tab(sfc_flux_dif_nir.data(), {0}, {static_cast<int>(sfc_flux_dif_nir.extent(0))});
    Table1D<Real>              mu0_tab(mu0.data(),              {0}, {static_cast<int>(mu0.extent(0))});
    Vector<Table1D<Real>> rrtmgp_out_vars = {mu0_tab                , sfc_flux_sw_dn_tab     ,
                                             sfc_flux_sw_dir_vis_tab, sfc_flux_sw_dir_nir_tab,
                                             sfc_flux_sw_dif_vis_tab, sfc_flux_sw_dif_nir_tab,
                                             sfc_flux_lw_dn_tab     };

    const Real rdOcp = m_rdOcp;
    for (MFIter mfi(*m_cons_in); mfi.isValid(); ++mfi) {
        const auto& vbx      = mfi.validbox();
        const auto& sbx      = makeSlab(vbx,2,vbx.smallEnd(2));
        const int nx         = vbx.length(0);
        const int imin       = vbx.smallEnd(0);
        const int jmin       = vbx.smallEnd(1);
        const int ktop       = vbx.bigEnd(2);
        const int offset     = m_col_offsets[mfi.index()];
        const Array4<Real>& q_arr = m_qheating_rates->array(mfi);
        const Array4<Real>& f_arr = m_rad_fluxes->array(mfi);
        ParallelFor(vbx, [=]
                    AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // map [i,j,k] 0-based to [icol, ilay] 0-based
            const int icol = (j-jmin)*nx + (i-imin) + offset;
            const int ilay = k;

            // Temperature heating rate for SW and LW
            q_arr(i,j,k,0) = sw_heating_tab(icol,ilay);
            q_arr(i,j,k,1) = lw_heating_tab(icol,ilay);

            // Convert the dT/dz to dTheta/dz
            Real iexner = rrtmgp::inverse_exner(Real(p_lay_tab(icol,ilay)), rdOcp);
            q_arr(i,j,k,0) *= iexner;
            q_arr(i,j,k,1) *= iexner;

            // Populate the fluxes: level ilay is the lower interface of
            // layer k, and the top-of-atmosphere level nlay goes into the
            // z-ghost cell above the top layer (rad_fluxes has one).
            f_arr(i,j,k,0) = sw_flux_up_tab(icol,ilay);
            f_arr(i,j,k,1) = sw_flux_dn_tab(icol,ilay);
            f_arr(i,j,k,2) = lw_flux_up_tab(icol,ilay);
            f_arr(i,j,k,3) = lw_flux_dn_tab(icol,ilay);
            if (k == ktop) {
                f_arr(i,j,k+1,0) = sw_flux_up_tab(icol,ilay+1);
                f_arr(i,j,k+1,1) = sw_flux_dn_tab(icol,ilay+1);
                f_arr(i,j,k+1,2) = lw_flux_up_tab(icol,ilay+1);
                f_arr(i,j,k+1,3) = lw_flux_dn_tab(icol,ilay+1);
            }

            if (k==0) {
                sfc_flux_sw_dn_tab(icol) = sw_flux_dn_tab(icol,ilay);
                sfc_flux_lw_dn_tab(icol) = lw_flux_dn_tab(icol,ilay);
            }
        });
        for (int ivar(0); ivar<lsm_output_ptrs.size(); ivar++) {
            if (lsm_output_ptrs[ivar]) {
                auto rrtmgp_for_fill = rrtmgp_out_vars[ivar];
                const Array4<Real>& lsm_out_arr = lsm_output_ptrs[ivar]->array(mfi);
                ParallelFor(sbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    // map [i,j,k] 0-based to [icol, ilay] 0-based
                    const int icol   = (j-jmin)*nx + (i-imin) + offset;

                    // export the desired variable at surface
                    lsm_out_arr(i,j,k) = rrtmgp_for_fill(icol);
                });
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


   std::string plotfilename = amrex::Concatenate("plt_rad_lev" + std::to_string(m_lev) + "_", m_step, 5);
   Vector<std::string> flux_names = {"sw_flux_up", "sw_flux_dn", "sw_flux_dir",
                                     "lw_flux_up", "lw_flux_dn"};
   WriteSingleLevelPlotfile(plotfilename, mf_flux, flux_names, m_geom, static_cast<Real>(m_time), m_step);
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

    // Layer pressure, to convert the clear-sky heating rates to potential temperature tendencies
    Table2D<Real,Order::C> p_lay_tab(p_lay.data(), {0,0}, {static_cast<int>(p_lay.extent(0)),static_cast<int>(p_lay.extent(1))});
    const Real rdOcp = m_rdOcp;

    auto extra_clnsky_diag    = m_extra_clnsky_diag;
    auto extra_clnclrsky_diag = m_extra_clnclrsky_diag;

    Real dz = m_geom.CellSize(2);
    const int ktop = m_geom.Domain().bigEnd(2);

    // Local copies of the layout constants for the kernel. The per-cell flux arrays are sized
    // by max_datalog_flux, so the enabled flux count must fit in them.
    constexpr int max_flux = max_datalog_flux;
    constexpr int nset     = n_flux_per_set;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nDatalogFluxComp() <= max_flux,
                                     "RRTMGP datalog flux count exceeds max_datalog_flux");

    // Every write goes through a guard, so a count that does not match the allocation never
    // indexes past the end; it raises this flag instead, which the host checks after the loop
    // since device-side asserts are compiled out of release GPU builds
    const int ncomp = datalog_mf.nComp();
    const int nflux = nDatalogFluxComp();
    Gpu::Buffer<int> d_err({0});
    int* err = d_err.data();

    for (MFIter mfi(datalog_mf); mfi.isValid(); ++mfi) {
        const auto& vbx      = mfi.validbox();
        const int nx         = vbx.length(0);
        const int imin       = vbx.smallEnd(0);
        const int jmin       = vbx.smallEnd(1);
        const int offset     = m_col_offsets[mfi.index()];
        const Array4<Real>& dst_arr = datalog_mf.array(mfi);
        const Array4<Real>& bdy_arr = datalog_bdy_mf.array(mfi);
        const Array4<Real>&   q_arr = m_qheating_rates->array(mfi);
        // Heights above ground need both the node heights and the surface heights
        const bool has_z = (m_z_phys != nullptr) && (m_z_sfc != nullptr);
        const Array4<const Real>& z_arr    = (has_z) ? m_z_phys->const_array(mfi) :
                                                       Array4<const Real>{};
        const Array4<const Real>& zsfc_arr = (has_z) ? m_z_sfc->const_array(mfi) :
                                                       Array4<const Real>{};

        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // map [i,j,k] 0-based to [icol, ilay] 0-based
            const int icol = (j-jmin)*nx + (i-imin) + offset;
            const int ilay = k;

            // Fluxes on the lower (flo) and upper (fhi) faces of this layer, in datalog order:
            // all-sky, clear-sky, then the optional clean-sky and clean-clear-sky sets
            Real flo[max_flux];
            Real fhi[max_flux];
            int nf = 0;
            auto add_flux = [&] (const Table2D<Real,Order::C>& tab)
            {
                if (nf < max_flux) {
                    flo[nf] = tab(icol,ilay);
                    fhi[nf] = tab(icol,ilay+1);
                }
                ++nf;
            };

            add_flux(sw_flux_up_tab);
            add_flux(sw_flux_dn_tab);
            add_flux(sw_flux_dn_dir_tab);
            add_flux(lw_flux_up_tab);
            add_flux(lw_flux_dn_tab);

            add_flux(sw_clrsky_flux_up_tab);
            add_flux(sw_clrsky_flux_dn_tab);
            add_flux(sw_clrsky_flux_dn_dir_tab);
            add_flux(lw_clrsky_flux_up_tab);
            add_flux(lw_clrsky_flux_dn_tab);

            if (extra_clnsky_diag) {
                add_flux(sw_clnsky_flux_up_tab);
                add_flux(sw_clnsky_flux_dn_tab);
                add_flux(sw_clnsky_flux_dn_dir_tab);
                add_flux(lw_clnsky_flux_up_tab);
                add_flux(lw_clnsky_flux_dn_tab);
            }

            if (extra_clnclrsky_diag) {
                add_flux(sw_clnclrsky_flux_up_tab);
                add_flux(sw_clnclrsky_flux_dn_tab);
                add_flux(sw_clnclrsky_flux_dn_dir_tab);
                add_flux(lw_clnclrsky_flux_up_tab);
                add_flux(lw_clnclrsky_flux_dn_tab);
            }

            const int nf_set = amrex::min(nf, nflux);

            int idx = 0;
            auto put = [&] (Real val)
            {
                if (idx < ncomp) { dst_arr(i,j,k,idx) = val; }
                ++idx;
            };

            // Height above ground level at cell center
            put( (z_arr) ? Compute_Zrel_AtCellCenter(i,j,k,zsfc_arr,z_arr) : (k + myhalf)*dz );

            // Heating rates
            put( q_arr(i, j, k, 0) );
            put( q_arr(i, j, k, 1) );

            // All-sky fluxes, averaged from the faces to the cell center
            for (int n = 0; n < nset && n < nf_set; ++n) { put( myhalf * (flo[n] + fhi[n]) ); }

            // Cosine zenith angle
            put( mu0_tab(icol) );

            // Clear sky heating rates, converted from dT/dt to dTheta/dt to match the
            // all-sky rates, which kokkos_buffers_to_mf() has already converted
            const Real iexner = rrtmgp::inverse_exner(Real(p_lay_tab(icol,ilay)), rdOcp);
            put( sw_clrsky_heating_tab(icol, ilay) * iexner );
            put( lw_clrsky_heating_tab(icol, ilay) * iexner );

            // Clear, clean and clean-clear sky fluxes, averaged to the cell center
            for (int n = nset; n < nf_set; ++n) { put( myhalf * (flo[n] + fhi[n]) ); }

            // Fluxes on the surface and top-of-atmosphere faces, and the top-face height
            if (k == 0) {
                for (int n = 0; n < nf_set; ++n) { bdy_arr(i,j,0,bdySfcComp(n)) = flo[n]; }
            }
            if (k == ktop) {
                for (int n = 0; n < nf_set; ++n) { bdy_arr(i,j,0,bdyToaComp(n,nflux)) = fhi[n]; }
                bdy_arr(i,j,0,bdyZtopComp(nflux)) = (z_arr) ?
                    fourth*( z_arr(i,j,k+1) + z_arr(i+1,j,k+1) + z_arr(i,j+1,k+1) + z_arr(i+1,j+1,k+1) )
                  - fourth*( zsfc_arr(i,j,0) + zsfc_arr(i+1,j,0) + zsfc_arr(i,j+1,0) + zsfc_arr(i+1,j+1,0) )
                  : (k + 1)*dz;
            }

            // The write counts must match the components allocated in set_grids
            if (idx != ncomp || nf != nflux) { *err = 1; }
        });
    }

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(*(d_err.copyToHost()) == 0,
                                     "RRTMGP datalog write count does not match the datalog allocation");
}

void Radiation::WriteDataLog (const double &time)
{
    // Columns after time, in file order. The datalog_mf components follow this order with the
    // disabled sets left out, and the flux columns also index the fluxes in datalog_bdy_mf.
    enum class Kind { Height, Heating, Mu0, Flux };
    struct Column { const char* name; Kind kind; bool enabled; };
    const bool cln    = m_extra_clnsky_diag;
    const bool clnclr = m_extra_clnclrsky_diag;
    const Vector<Column> cols = {
        {"z_agl"                 , Kind::Height , true  },
        {"sw_heating"            , Kind::Heating, true  },
        {"lw_heating"            , Kind::Heating, true  },
        {"sw_flux_up"            , Kind::Flux   , true  },
        {"sw_flux_dn"            , Kind::Flux   , true  },
        {"sw_flux_dn_dir"        , Kind::Flux   , true  },
        {"lw_flux_up"            , Kind::Flux   , true  },
        {"lw_flux_dn"            , Kind::Flux   , true  },
        {"mu0"                   , Kind::Mu0    , true  },
        {"sw_clrsky_heating"     , Kind::Heating, true  },
        {"lw_clrsky_heating"     , Kind::Heating, true  },
        {"sw_clrsky_flux_up"     , Kind::Flux   , true  },
        {"sw_clrsky_flux_dn"     , Kind::Flux   , true  },
        {"sw_clrsky_flux_dn_dir" , Kind::Flux   , true  },
        {"lw_clrsky_flux_up"     , Kind::Flux   , true  },
        {"lw_clrsky_flux_dn"     , Kind::Flux   , true  },
        {"sw_clnsky_flux_up"     , Kind::Flux   , cln   },
        {"sw_clnsky_flux_dn"     , Kind::Flux   , cln   },
        {"sw_clnsky_flux_dn_dir" , Kind::Flux   , cln   },
        {"lw_clnsky_flux_up"     , Kind::Flux   , cln   },
        {"lw_clnsky_flux_dn"     , Kind::Flux   , cln   },
        {"sw_clnclrsky_flux_up"    , Kind::Flux , clnclr},
        {"sw_clnclrsky_flux_dn"    , Kind::Flux , clnclr},
        {"sw_clnclrsky_flux_dn_dir", Kind::Flux , clnclr},
        {"lw_clnclrsky_flux_up"    , Kind::Flux , clnclr},
        {"lw_clnclrsky_flux_dn"    , Kind::Flux , clnclr}
    };
    const int ncols = static_cast<int>(cols.size());

    // Map each enabled column to its datalog_mf component and, for fluxes, its flux index
    Vector<int> cc_comp(ncols, -1);
    Vector<int> flux_comp(ncols, -1);
    int icc   = 0;
    int iflux = 0;
    for (int c = 0; c < ncols; ++c) {
        if (cols[c].enabled) {
            cc_comp[c] = icc++;
            if (cols[c].kind == Kind::Flux) { flux_comp[c] = iflux++; }
        }
    }
    const int nflux = nDatalogFluxComp();
    AMREX_ALWAYS_ASSERT(icc == datalog_mf.nComp() && iflux == nflux &&
                        datalog_bdy_mf.nComp() == nDatalogBdyComp());

    auto domain = m_geom.Domain();
    const Real area_z = static_cast<Real>(domain.length(0)*domain.length(1));
    const int nz = domain.length(2);

    // Planar averages of the cell-centred profiles
    Vector<Gpu::HostVector<Real>> cc_avg(icc);
    for (int n = 0; n < icc; ++n) {
        cc_avg[n] = sumToLine(datalog_mf, n, 1, domain, 2);
        for (auto& val : cc_avg[n]) { val /= area_z; }
    }

    // Planar averages of the surface and top-of-atmosphere values
    Vector<Real> bdy_avg(datalog_bdy_mf.nComp());
    for (int n = 0; n < datalog_bdy_mf.nComp(); ++n) { bdy_avg[n] = datalog_bdy_mf.sum(n, true); }
    ParallelDescriptor::ReduceRealSum(bdy_avg.data(), static_cast<int>(bdy_avg.size()));
    for (auto& val : bdy_avg) { val /= area_z; }

    if (ParallelDescriptor::IOProcessor()) {
        std::ostream& log = *datalog;
        if (log.good()) {
            // The log is opened in append mode, so only a new file gets the header
            log.seekp(0, std::ios::end);
            if (log.tellp() == 0) {
                log << "# Planar averages; each output is a block of nz+2 rows: the surface face,\n"
                    << "# the nz cell centers, then the top-of-atmosphere face. z_agl is the height\n"
                    << "# above ground. Heating rates exist only at cell centers and are nan on the\n"
                    << "# face rows. Disabled clean-sky and clean-clear-sky columns are 0.\n"
                    << "# time";
                for (const auto& col : cols) { log << " " << col.name; }
                log << std::endl;
            }

            const Real nan = std::numeric_limits<Real>::quiet_NaN();

            // Row -1 is the surface face, rows 0..nz-1 the cell centers, row nz the top face
            for (int row = -1; row <= nz; ++row)
            {
                const bool is_face = (row < 0) || (row == nz);
                const bool is_top  = (row == nz);

                log << std::setw(datalog_width) << std::setprecision(datalog_time_precision) << time
                    << std::setprecision(datalog_precision);
                for (int c = 0; c < ncols; ++c) {
                    Real val = zero;
                    if (cols[c].enabled) {
                        if (!is_face) {
                            val = cc_avg[cc_comp[c]][row];
                        } else {
                            switch (cols[c].kind) {
                                case Kind::Height:
                                    val = (is_top) ? bdy_avg[bdyZtopComp(nflux)] : zero;
                                    break;
                                case Kind::Heating:
                                    val = nan;
                                    break;
                                case Kind::Mu0:
                                    val = cc_avg[cc_comp[c]][0];
                                    break;
                                case Kind::Flux:
                                    val = bdy_avg[(is_top) ? bdyToaComp(flux_comp[c], nflux)
                                                           : bdySfcComp(flux_comp[c])];
                                    break;
                            }
                        }
                    }
                    log << " ";
                    if (c == 0) { log << std::setw(datalog_width); }
                    log << val;
                }
                log << std::endl;
            } // for row
        } // log is good
    } // IOProc
}

void
Radiation::initialize_impl ()
{
    // Initialize gas concentrations for this step
    m_gas_concs.init(gas_names_offset, m_ncol, m_nlay);

    // Load k-distribution and cloud optics data only once.
    // These are static lookup tables that never change.
    // Size the memory pool for the requested chunk size (not the effective one, and
    // not min with the current m_ncol) so that the pool remains valid even if m_ncol
    // grows after regridding/load balancing. The pool is created once and never
    // resized, whereas the effective chunk size is recomputed at every Init() and is
    // bounded above by the request.
    if (!rrtmgp::initialized) {
        gas_concs_t gas_concs_pool;
        gas_concs_pool.init(gas_names_offset, m_ncol_chunk_requested, m_nlay);
        rrtmgp::rrtmgp_initialize(gas_concs_pool,
                                  rrtmgp_coeffs_file_sw      , rrtmgp_coeffs_file_lw      ,
                                  rrtmgp_cloud_optics_file_sw, rrtmgp_cloud_optics_file_lw,
                                  m_rad_nvar);
        gas_concs_pool.reset();
    }
}


void
Radiation::run_impl ()
{
    // A rank that owns no boxes on this level has no columns and therefore no
    // radiation work to do. Bail out before the chunk loop; there are no MPI
    // collectives in this routine, so returning early cannot deadlock.
    if (m_ncol == 0) { return; }

    // Local copies
    const auto ncol     = m_ncol;
    const auto nlay     = m_nlay;
    const auto nswbands = m_nswbands;

    // Compute orbital parameters; these are used both for computing
    // the solar zenith angle and also for computing total solar
    // irradiance scaling (tsi_scaling).
    double obliqr, lambm0, mvelpp;
    // Any of eccen/obliq/mvelp that the user set (i.e. is non-negative) is used
    // as is by orbital_params; the rest are computed from the orbital year.
    int  orbital_year = m_orbital_year;
    double eccen      = m_orbital_eccen;
    double obliq      = m_orbital_obliq;
    double mvelp      = m_orbital_mvelp;
    orbital_params(orbital_year, eccen, obliq,
                   mvelp, obliqr, lambm0, mvelpp);

    // Use the orbital parameters to calculate the solar declination and eccentricity factor
    double delta, eccf;
    // Day of the year plus fraction, calday 1 == Jan 1 0Z (leap-aware)
    double calday = orbital_calday(m_orbital_year, m_orbital_mon, m_orbital_day, m_orbital_sec);
    orbital_decl(calday, eccen, mvelpp, lambm0, obliqr, delta, eccf);
    m_calday = calday;
    m_declin = delta;

    // Overwrite eccf if using a fixed solar constant.
    auto fixed_total_solar_irradiance = m_fixed_total_solar_irradiance;
    if (fixed_total_solar_irradiance >= 0){
       eccf = fixed_total_solar_irradiance/Real(1360.9);
    }

    // Precompute volume mixing ratio (VMR) for all gases
    //
    // H2O is obtained from qv.
    // O3 may be a constant or a 1D vector
    // All other comps are set to constants for now
    Vector<real2d_k> vmr_full_vec(m_ngas);
    for (int igas(0); igas < m_ngas; ++igas) {
        auto name = m_gas_names[igas];
        vmr_full_vec[igas] = real2d_k("vmr_full_" + name, ncol, nlay);
        auto tmp2d = vmr_full_vec[igas];
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
    }

    // Determine the cosine zenith angle.
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
                             [&,PI_d=PI] (int icol)
        {
            // Convert lat/lon to radians
            double lat_col = h_lat(icol)*PI_d/Real(180.0);
            double lon_col = h_lon(icol)*PI_d/Real(180.0);
            double lcalday = calday;
            double ldelta  = delta;
            double dt_avg  = static_cast<double>(rad_freq_in_steps) * dt;
            h_mu0(icol)    = Real(orbital_cos_zenith(lcalday, lat_col, lon_col, ldelta, dt_avg));
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
        lwp_tab(icol,ilay) *= Real(1.e3);
        iwp_tab(icol,ilay) *= Real(1.e3);
    });

    // -----------------------------------------------------------------------
    // Process radiation in column chunks to limit peak GPU memory.
    // Radiation columns are independent (no horizontal coupling), so
    // chunking produces bit-identical results.
    // -----------------------------------------------------------------------
    const int ncol_chunk = std::min(m_ncol_chunk, ncol);
    const int kbot = 0;

    for (int col_s = 0; col_s < ncol; col_s += ncol_chunk) {
        const int ncol_c = std::min(ncol_chunk, ncol - col_s);

        // --- Chunk subviews: 1D (ncol) ---
        real1d_k mu0_c              (mu0.data()              + col_s, ncol_c);
        real1d_k sfc_alb_dir_vis_c  (sfc_alb_dir_vis.data()  + col_s, ncol_c);
        real1d_k sfc_alb_dir_nir_c  (sfc_alb_dir_nir.data()  + col_s, ncol_c);
        real1d_k sfc_alb_dif_vis_c  (sfc_alb_dif_vis.data()  + col_s, ncol_c);
        real1d_k sfc_alb_dif_nir_c  (sfc_alb_dif_nir.data()  + col_s, ncol_c);
        real1d_k sfc_flux_dir_vis_c (sfc_flux_dir_vis.data() + col_s, ncol_c);
        real1d_k sfc_flux_dir_nir_c (sfc_flux_dir_nir.data() + col_s, ncol_c);
        real1d_k sfc_flux_dif_vis_c (sfc_flux_dif_vis.data() + col_s, ncol_c);
        real1d_k sfc_flux_dif_nir_c (sfc_flux_dif_nir.data() + col_s, ncol_c);
        real1d_k t_sfc_c            (t_sfc.data()            + col_s, ncol_c);
        real1d_k sfc_emis_c         (sfc_emis.data()         + col_s, ncol_c);
        real1d_k lw_src_c           (lw_src.data()           + col_s, ncol_c);

        // --- Chunk subviews: 2D (ncol, nlay) via LayoutRight pointer offset ---
        const int stride2_nlay   = nlay;
        const int stride2_nlayp1 = nlay + 1;
        real2d_k p_lay_c        (p_lay.data()         + col_s*stride2_nlay, ncol_c, nlay);
        real2d_k t_lay_c        (t_lay.data()         + col_s*stride2_nlay, ncol_c, nlay);
        real2d_k r_lay_c        (r_lay.data()         + col_s*stride2_nlay, ncol_c, nlay);
        real2d_k z_del_c        (z_del.data()         + col_s*stride2_nlay, ncol_c, nlay);
        real2d_k lwp_c          (lwp.data()           + col_s*stride2_nlay, ncol_c, nlay);
        real2d_k iwp_c          (iwp.data()           + col_s*stride2_nlay, ncol_c, nlay);
        real2d_k eff_radius_qc_c(eff_radius_qc.data() + col_s*stride2_nlay, ncol_c, nlay);
        real2d_k eff_radius_qi_c(eff_radius_qi.data() + col_s*stride2_nlay, ncol_c, nlay);
        real2d_k cldfrac_tot_c  (cldfrac_tot.data()   + col_s*stride2_nlay, ncol_c, nlay);
        real2d_k sw_heating_c   (sw_heating.data()    + col_s*stride2_nlay, ncol_c, nlay);
        real2d_k lw_heating_c   (lw_heating.data()    + col_s*stride2_nlay, ncol_c, nlay);

        // --- Chunk subviews: 2D (ncol, nlay+1) ---
        real2d_k p_lev_c                   (p_lev.data()                 + col_s*stride2_nlayp1, ncol_c, nlay+1);
        real2d_k t_lev_c                   (t_lev.data()                 + col_s*stride2_nlayp1, ncol_c, nlay+1);
        real2d_k sw_flux_up_c              (sw_flux_up.data()            + col_s*stride2_nlayp1, ncol_c, nlay+1);
        real2d_k sw_flux_dn_c              (sw_flux_dn.data()            + col_s*stride2_nlayp1, ncol_c, nlay+1);
        real2d_k sw_flux_dn_dir_c          (sw_flux_dn_dir.data()        + col_s*stride2_nlayp1, ncol_c, nlay+1);
        real2d_k lw_flux_up_c              (lw_flux_up.data()            + col_s*stride2_nlayp1, ncol_c, nlay+1);
        real2d_k lw_flux_dn_c              (lw_flux_dn.data()            + col_s*stride2_nlayp1, ncol_c, nlay+1);
        // Clear-sky flux subviews (always active)
        // NOTE: once on m_ncol_chunk if not writing a datalog
        real2d_k sw_clrsky_flux_up_c, sw_clrsky_flux_dn_c, sw_clrsky_flux_dn_dir_c;
        real2d_k lw_clrsky_flux_up_c, lw_clrsky_flux_dn_c;
        if (datalog_int > 0) {
            sw_clrsky_flux_up_c     = real2d_k(sw_clrsky_flux_up.data()     + col_s*stride2_nlayp1, ncol_c, nlay+1);
            sw_clrsky_flux_dn_c     = real2d_k(sw_clrsky_flux_dn.data()     + col_s*stride2_nlayp1, ncol_c, nlay+1);
            sw_clrsky_flux_dn_dir_c = real2d_k(sw_clrsky_flux_dn_dir.data() + col_s*stride2_nlayp1, ncol_c, nlay+1);
            lw_clrsky_flux_up_c     = real2d_k(lw_clrsky_flux_up.data()     + col_s*stride2_nlayp1, ncol_c, nlay+1);
            lw_clrsky_flux_dn_c     = real2d_k(lw_clrsky_flux_dn.data()     + col_s*stride2_nlayp1, ncol_c, nlay+1);
        } else {
            sw_clrsky_flux_up_c     = real2d_k(sw_clrsky_flux_up.data()     , ncol_c, nlay+1);
            sw_clrsky_flux_dn_c     = real2d_k(sw_clrsky_flux_dn.data()     , ncol_c, nlay+1);
            sw_clrsky_flux_dn_dir_c = real2d_k(sw_clrsky_flux_dn_dir.data() , ncol_c, nlay+1);
            lw_clrsky_flux_up_c     = real2d_k(lw_clrsky_flux_up.data()     , ncol_c, nlay+1);
            lw_clrsky_flux_dn_c     = real2d_k(lw_clrsky_flux_dn.data()     , ncol_c, nlay+1);
        }

        // Diagnostic flux subviews (placeholder when disabled)
        real2d_k sw_clnclrsky_flux_up_c, sw_clnclrsky_flux_dn_c, sw_clnclrsky_flux_dn_dir_c;
        real2d_k lw_clnclrsky_flux_up_c, lw_clnclrsky_flux_dn_c;
        if (m_extra_clnclrsky_diag) {
            sw_clnclrsky_flux_up_c     = real2d_k(sw_clnclrsky_flux_up.data()     + col_s*stride2_nlayp1, ncol_c, nlay+1);
            sw_clnclrsky_flux_dn_c     = real2d_k(sw_clnclrsky_flux_dn.data()     + col_s*stride2_nlayp1, ncol_c, nlay+1);
            sw_clnclrsky_flux_dn_dir_c = real2d_k(sw_clnclrsky_flux_dn_dir.data() + col_s*stride2_nlayp1, ncol_c, nlay+1);
            lw_clnclrsky_flux_up_c     = real2d_k(lw_clnclrsky_flux_up.data()     + col_s*stride2_nlayp1, ncol_c, nlay+1);
            lw_clnclrsky_flux_dn_c     = real2d_k(lw_clnclrsky_flux_dn.data()     + col_s*stride2_nlayp1, ncol_c, nlay+1);
        } else {
            sw_clnclrsky_flux_up_c     = real2d_k("sw_clnclrsky_flux_up_c"    , 1, 1);
            sw_clnclrsky_flux_dn_c     = real2d_k("sw_clnclrsky_flux_dn_c"    , 1, 1);
            sw_clnclrsky_flux_dn_dir_c = real2d_k("sw_clnclrsky_flux_dn_dir_c", 1, 1);
            lw_clnclrsky_flux_up_c     = real2d_k("lw_clnclrsky_flux_up_c"    , 1, 1);
            lw_clnclrsky_flux_dn_c     = real2d_k("lw_clnclrsky_flux_dn_c"    , 1, 1);
        }

        real2d_k sw_clnsky_flux_up_c, sw_clnsky_flux_dn_c, sw_clnsky_flux_dn_dir_c;
        real2d_k lw_clnsky_flux_up_c, lw_clnsky_flux_dn_c;
        if (m_extra_clnsky_diag) {
            sw_clnsky_flux_up_c     = real2d_k(sw_clnsky_flux_up.data()     + col_s*stride2_nlayp1, ncol_c, nlay+1);
            sw_clnsky_flux_dn_c     = real2d_k(sw_clnsky_flux_dn.data()     + col_s*stride2_nlayp1, ncol_c, nlay+1);
            sw_clnsky_flux_dn_dir_c = real2d_k(sw_clnsky_flux_dn_dir.data() + col_s*stride2_nlayp1, ncol_c, nlay+1);
            lw_clnsky_flux_up_c     = real2d_k(lw_clnsky_flux_up.data()     + col_s*stride2_nlayp1, ncol_c, nlay+1);
            lw_clnsky_flux_dn_c     = real2d_k(lw_clnsky_flux_dn.data()     + col_s*stride2_nlayp1, ncol_c, nlay+1);
        } else {
            sw_clnsky_flux_up_c     = real2d_k("sw_clnsky_flux_up_c"    , 1, 1);
            sw_clnsky_flux_dn_c     = real2d_k("sw_clnsky_flux_dn_c"    , 1, 1);
            sw_clnsky_flux_dn_dir_c = real2d_k("sw_clnsky_flux_dn_dir_c", 1, 1);
            lw_clnsky_flux_up_c     = real2d_k("lw_clnsky_flux_up_c"    , 1, 1);
            lw_clnsky_flux_dn_c     = real2d_k("lw_clnsky_flux_dn_c"    , 1, 1);
        }

        // --- Chunk subviews: 2D (ncol, nswbands) ---
        real2d_k sfc_alb_dir_c(sfc_alb_dir.data() + col_s*nswbands, ncol_c, nswbands);
        real2d_k sfc_alb_dif_c(sfc_alb_dif.data() + col_s*nswbands, ncol_c, nswbands);

        // --- Chunk subviews: 3D (ncol, nlay+1, nbands) ---
        // NOTE: Allocate these once on m_ncol_chunk and use what we need in the chunk loop
        real3d_k sw_bnd_flux_up_c (sw_bnd_flux_up.data() , ncol_c, nlay+1, nswbands);
        real3d_k sw_bnd_flux_dn_c (sw_bnd_flux_dn.data() , ncol_c, nlay+1, nswbands);
        real3d_k sw_bnd_flux_dir_c(sw_bnd_flux_dir.data(), ncol_c, nlay+1, nswbands);
        real3d_k sw_bnd_flux_dif_c(sw_bnd_flux_dif.data(), ncol_c, nlay+1, nswbands);
        real3d_k lw_bnd_flux_up_c (lw_bnd_flux_up.data() , ncol_c, nlay+1, m_nlwbands);
        real3d_k lw_bnd_flux_dn_c (lw_bnd_flux_dn.data() , ncol_c, nlay+1, m_nlwbands);

        // --- Create chunk gas concentrations by subsetting from pre-fetched VMR ---
        gas_concs_t gas_concs_c;
        gas_concs_c.init(gas_names_offset, ncol_c, nlay);
        for (int igas = 0; igas < m_ngas; ++igas) {
            real2d_k vmr_c("vmr_c", ncol_c, nlay);
            auto vmr_full = vmr_full_vec[igas];
            auto cs = col_s;
            Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol_c, nlay}),
                                 KOKKOS_LAMBDA (int i, int j) {
                vmr_c(i, j) = vmr_full(cs + i, j);
            });
            gas_concs_c.set_vmr(m_gas_names[igas], vmr_c);
        }

        // Expand surface albedos along nswbands for this chunk
        rrtmgp::compute_band_by_band_surface_albedos(ncol_c, nswbands,
                                                     sfc_alb_dir_vis_c, sfc_alb_dir_nir_c,
                                                     sfc_alb_dif_vis_c, sfc_alb_dif_nir_c,
                                                     sfc_alb_dir_c    , sfc_alb_dif_c);

        // Run RRTMGP driver for this column chunk
        rrtmgp::rrtmgp_main(ncol_c, m_nlay,
                            p_lay_c, t_lay_c,
                            p_lev_c, t_lev_c,
                            gas_concs_c,
                            sfc_alb_dir_c, sfc_alb_dif_c, mu0_c,
                            t_sfc_c, sfc_emis_c, lw_src_c,
                            lwp_c, iwp_c, eff_radius_qc_c, eff_radius_qi_c, cldfrac_tot_c,
                            sw_flux_up_c, sw_flux_dn_c, sw_flux_dn_dir_c,
                            lw_flux_up_c, lw_flux_dn_c,
                            sw_clnclrsky_flux_up_c, sw_clnclrsky_flux_dn_c, sw_clnclrsky_flux_dn_dir_c,
                            sw_clrsky_flux_up_c, sw_clrsky_flux_dn_c, sw_clrsky_flux_dn_dir_c,
                            sw_clnsky_flux_up_c, sw_clnsky_flux_dn_c, sw_clnsky_flux_dn_dir_c,
                            lw_clnclrsky_flux_up_c, lw_clnclrsky_flux_dn_c,
                            lw_clrsky_flux_up_c, lw_clrsky_flux_dn_c,
                            lw_clnsky_flux_up_c, lw_clnsky_flux_dn_c,
                            sw_bnd_flux_up_c, sw_bnd_flux_dn_c, sw_bnd_flux_dir_c,
                            lw_bnd_flux_up_c, lw_bnd_flux_dn_c,
                            eccf, m_do_subcol_sampling,
                            m_extra_clnclrsky_diag, m_extra_clnsky_diag);

        // Compute heating rates for this chunk
        rrtmgp::compute_heating_rate(sw_flux_up_c, sw_flux_dn_c, r_lay_c, z_del_c, sw_heating_c);
        rrtmgp::compute_heating_rate(lw_flux_up_c, lw_flux_dn_c, r_lay_c, z_del_c, lw_heating_c);

        // Compute diffuse band fluxes and broadband surface fluxes for this chunk
        Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol_c, nlay+1, nswbands}),
                             KOKKOS_LAMBDA (int icol, int ilay, int ibnd)
        {
            sw_bnd_flux_dif_c(icol,ilay,ibnd) = sw_bnd_flux_dn_c(icol,ilay,ibnd) - sw_bnd_flux_dir_c(icol,ilay,ibnd);
        });
        rrtmgp::compute_broadband_surface_fluxes(ncol_c, kbot, nswbands,
                                                 sw_bnd_flux_dir_c , sw_bnd_flux_dif_c ,
                                                 sfc_flux_dir_vis_c, sfc_flux_dir_nir_c,
                                                 sfc_flux_dif_vis_c, sfc_flux_dif_nir_c);

        gas_concs_c.reset();
    } // end column chunk loop
}


void
Radiation::finalize_impl (const Vector<MultiFab*>& lsm_output_ptrs)
{
    // Reset gas concentrations (k-dist data persists across steps)
    m_gas_concs.reset();

    // Fill the AMReX MFs from Kokkos Views
    kokkos_buffers_to_mf(lsm_output_ptrs);

    // Write fluxes if requested
    if (m_rad_write_fluxes) { write_rrtmgp_fluxes(); }

    // Fill output data for datalog before deallocating
    if (datalog_int > 0) {
        rrtmgp::compute_heating_rate(sw_clrsky_flux_up, sw_clrsky_flux_dn, r_lay, z_del, sw_clrsky_heating);
        rrtmgp::compute_heating_rate(lw_clrsky_flux_up, lw_clrsky_flux_dn, r_lay, z_del, lw_clrsky_heating);
        Kokkos::fence();
        populateDatalogMF();
    }

    // Deallocate the buffer arrays
    dealloc_buffers();
}
