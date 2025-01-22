#ifndef ERF_RADIATION_H
#define ERF_RADIATION_H

/*
 * RTE-RRTMGP radiation model interface to ERF
 * The original code is developed by RobertPincus, and the code is open source available at:
 *                        https://github.com/earth-system-radiation/rte-rrtmgp
 * Please reference to the following paper,
 *                        https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019MS001621
 * NOTE: we use the C++ version of RTE-RRTMGP, which is the implementation of the original Fortran
 * code using C++ YAKL for CUDA, HiP and SYCL application by E3SM ECP team, the C++ version
 * of the rte-rrtmgp code is located at:
 *                       https://github.com/E3SM-Project/rte-rrtmgp
 * The RTE-RRTMGP uses BSD-3-Clause Open Source License, if you want to make changes,
 * and modifications to the code, please refer to BSD-3-Clause Open Source License.
 */

#include <ctime>
#include <string>
#include <vector>
#include <memory>

#include <mo_gas_concentrations.h>

#include <YAKL.h>

#include <AMReX_ParmParse.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_TableData.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

#include <ERF_RRTMGP_Interface.H>
#include <ERF_RRTMGP_Utils.H>
#include <ERF_OrbCosZenith.H>
#include <ERF_Constants.H>
#include <ERF_IndexDefines.H>
#include <ERF_DataStruct.H>
#include <ERF_EOS.H>

class Radiation {
public:

    // Constructor
    Radiation (const int& lev,
               SolverChoice& sc);

    // Destructor
    ~Radiation () = default;

    // Set the grid info for columnar data in YAKL
    void
    set_grids (int& level,
               int& step,
               amrex::Real& time,
               const amrex::Real& dt,
               const amrex::BoxArray& ba,
               amrex::Geometry& geom,
               amrex::MultiFab* cons_in,
               amrex::MultiFab* lsm_fluxes,
               amrex::MultiFab* lsm_zenith,
               amrex::MultiFab* qheating_rates,
               amrex::MultiFab* z_phys,
               amrex::MultiFab* lat,
               amrex::MultiFab* lon);

    // Initialize the temporary variables
    void
    alloc_buffers ();

    // Clear the temporary variables
    void
    dealloc_buffers ();

    // Fill YAKL Arrays from AMReX MultiFabs
    void
    mf_to_yakl_buffers();

    // Fill AMReX MultiFabs from YAKL Arrays
    void
    yakl_buffers_to_mf ();

    // Write the rrtmgp fluxes
    void
    write_rrtmgp_fluxes ();

    // Initialize the implementation
    void
    initialize_impl ();

    // Run the implementation
    void
    run_impl ();

    // Finalize the implementation
    void
    finalize_impl ();

    // Wrapper for implementation steps
    void
    rad_run_impl ()
    {
        if (m_update_rad) {
            amrex::Print() << "Advancing radiation at level: " << m_lev << " ...";
            this->initialize_impl();
            this->run_impl();
            this->finalize_impl();
            amrex::Print() << "DONE\n";
        }
    }

private:

    // Process interface vars from ERF/AMReX
    //===================================================================================

    // Grid level
    int m_lev;

    // Step number
    int m_step;

    // Current time
    amrex::Real m_time;

    // Timestep at given level
    amrex::Real m_dt;

    // Geometry at given level
    amrex::Geometry m_geom;

    // Boxarray at given level
    amrex::BoxArray m_ba;

    // Are we updating radiation?
    bool m_update_rad = false;

    // Are we writing out the fluxes?
    bool m_rad_write_fluxes = false;

    // Do we have moisture and cold comps?
    bool m_moist = false;
    bool m_ice   = false;

    // Pointer to the CC conserved vars
    amrex::MultiFab* m_cons_in = nullptr;

    // Pointer to the radiation source terms
    amrex::MultiFab* m_qheating_rates = nullptr;

    // Pointer to the terrain heights
    amrex::MultiFab* m_z_phys = nullptr;

    // Pointer to latitude and longitude
    amrex::MultiFab* m_lat = nullptr;
    amrex::MultiFab* m_lon = nullptr;

    // Pointer to output data for LSM
    amrex::MultiFab* m_lsm_fluxes = nullptr;
    amrex::MultiFab* m_lsm_zenith = nullptr;

    // Path, data file, and coefficient file for K-distribution
    std::string rrtmgp_file_path       = ".";
    std::string rrtmgp_coeffs_sw       = "rrtmgp-data-sw-g224-2018-12-04.nc";
    std::string rrtmgp_coeffs_lw       = "rrtmgp-data-lw-g224-2018-12-04.nc";
    std::string rrtmgp_cloud_optics_sw = "rrtmgp-cloud-optics-coeffs-sw.nc";
    std::string rrtmgp_cloud_optics_lw = "rrtmgp-cloud-optics-coeffs-lw.nc";
    std::string rrtmgp_coeffs_file_sw;
    std::string rrtmgp_coeffs_file_lw;
    std::string rrtmgp_cloud_optics_file_sw;
    std::string rrtmgp_cloud_optics_file_lw;

    // Active gases
    int m_ngas = 8;
    const std::vector<std::string> m_gas_names = {"H2O", "CO2", "O3", "N2O",
                                                  "CO" , "CH4", "O2", "N2" };
    const std::vector<amrex::Real> m_mol_weight_gas = {18.01528, 44.00950, 47.9982, 44.0128,
                                                       28.01010, 16.04246, 31.9980, 28.0134}; // g/mol

    // Prescribed greenhouse gas surface concentrations in moles / moles air
    amrex::Real m_co2vmr = 388.717e-6;
    amrex::Vector<amrex::Real> m_o3vmr;
    amrex::Real m_n2ovmr = 323.141e-9;
    amrex::Real m_covmr  = 1.0e-7;
    amrex::Real m_ch4vmr = 1807.851e-9;
    amrex::Real m_o2vmr  = 0.209448;
    amrex::Real m_n2vmr  = 0.7906;
    //amrex::Real m_f11vmr = 768.7644e-12;
    //amrex::Real m_f12vmr = 531.2820e-12;

    int m_o3_size;
    real1d    m_gas_mol_weights;
    string1dv gas_names_yakl_offset;
    GasConcs  m_gas_concs;

    // Process interface vars modeled after EAMXX
    //===================================================================================

    // Keep track of number of columns and levels
    int m_ncol;
    int m_nlay;

    // Offsets for MultiFab <-> YAKL transfer
    amrex::Vector<int> m_col_offsets;

    // Whether we use aerosol forcing in radiation
    bool m_do_aerosol_rad = true;

    // Whether we do extra aerosol forcing calls
    bool m_extra_clnsky_diag    = false;
    bool m_extra_clnclrsky_diag = false;

    // The orbital year, used for zenith angle calculations:
    // If > 0, use constant orbital year for duration of simulation
    // If < 0, use year from timestamp for orbital parameters
    int m_orbital_year = -9999;
    int m_orbital_mon  = -9999;
    int m_orbital_day  = -9999;
    int m_orbital_sec  = -9999;

    // Orbital parameters, used for zenith angle calculations.
    // If >= 0, bypass computation based on orbital year and use fixed parameters
    // If <  0, compute based on orbital year, specified above
    bool m_fixed_orbital_year   = false;
    amrex::Real m_orbital_eccen = -9999.;  // Eccentricity
    amrex::Real m_orbital_obliq = -9999.;  // Obliquity
    amrex::Real m_orbital_mvelp = -9999.;  // Vernal Equinox Mean Longitude of Perihelion

    // Value for prescribing an invariant solar constant (i.e. total solar irradiance
    // at TOA).  Used for idealized experiments such as RCE. This is only used when a
    // positive value is supplied.
    amrex::Real m_fixed_total_solar_irradiance = -9999.;

    // Fixed solar zenith angle to use for shortwave calculations
    // This is only used if a positive value is supplied
    amrex::Real m_fixed_solar_zenith_angle = -9999.;

    // Need to hard-code some dimension sizes for now.
    // TODO: find a better way of configuring this
    int m_nswbands = 14;
    int m_nlwbands = 16;
    int m_nswgpts  = 112;
    int m_nlwgpts  = 128;

    // Rad frequency in number of steps
    int m_rad_freq_in_steps = 1;

    // Whether or not to do subcolumn sampling of cloud state for MCICA
    bool m_do_subcol_sampling = true;

    // 1d size (1 or nlay)
    real1d o3_lay;

    // 1d size (ncol)
    real1d cosine_zenith;
    real1d mu0;
    real1d sfc_alb_dir_vis;
    real1d sfc_alb_dir_nir;
    real1d sfc_alb_dif_vis;
    real1d sfc_alb_dif_nir;
    real1d sfc_flux_dir_vis;
    real1d sfc_flux_dir_nir;
    real1d sfc_flux_dif_vis;
    real1d sfc_flux_dif_nir;
    real1d lat;
    real1d lon;

    // 2d size (ncol, nlay)
    real2d r_lay;
    real2d p_lay;
    real2d t_lay;
    real2d z_del;
    real2d p_del;
    real2d qv_lay;
    real2d qc_lay;
    real2d qi_lay;
    real2d cldfrac_tot;
    real2d eff_radius_qc;
    real2d eff_radius_qi;
    real2d tmp2d;
    real2d lwp;
    real2d iwp;
    real2d sw_heating;
    real2d lw_heating;

    // 2d size (ncol, nlay+1)
    real2d d_tint;
    real2d p_lev;
    real2d t_lev;
    real2d sw_flux_up;
    real2d sw_flux_dn;
    real2d sw_flux_dn_dir;
    real2d lw_flux_up;
    real2d lw_flux_dn;
    real2d sw_clnclrsky_flux_up;
    real2d sw_clnclrsky_flux_dn;
    real2d sw_clnclrsky_flux_dn_dir;
    real2d sw_clrsky_flux_up;
    real2d sw_clrsky_flux_dn;
    real2d sw_clrsky_flux_dn_dir;
    real2d sw_clnsky_flux_up;
    real2d sw_clnsky_flux_dn;
    real2d sw_clnsky_flux_dn_dir;
    real2d lw_clnclrsky_flux_up;
    real2d lw_clnclrsky_flux_dn;
    real2d lw_clrsky_flux_up;
    real2d lw_clrsky_flux_dn;
    real2d lw_clnsky_flux_up;
    real2d lw_clnsky_flux_dn;

    // 3d size (ncol, nlay+1, nswbands)
    real3d sw_bnd_flux_up;
    real3d sw_bnd_flux_dn;
    real3d sw_bnd_flux_dir;
    real3d sw_bnd_flux_dif;

    // 3d size (ncol, nlay+1, nlwbands)
    real3d lw_bnd_flux_up;
    real3d lw_bnd_flux_dn;

    // 2d size (ncol, nswbands)
    real2d sfc_alb_dir;
    real2d sfc_alb_dif;

    // 3d size (ncol, nlay, n[sw,lw]bands)
    real3d aero_tau_sw;
    real3d aero_ssa_sw;
    real3d aero_g_sw;
    real3d aero_tau_lw;

    // 3d size (ncol, nlay, n[sw,lw]bnds)
    real3d cld_tau_sw_bnd;
    real3d cld_tau_lw_bnd;

    // 3d size (ncol, nlay, n[sw,lw]gpts)
    real3d cld_tau_sw_gpt;
    real3d cld_tau_lw_gpt;
};
#endif
