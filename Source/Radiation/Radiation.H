/*
 * RTE-RRTMGP radiation model interface to ERF
 * The orginal code is developed by RobertPincus, and the code is open source available at:
 *                        https://github.com/earth-system-radiation/rte-rrtmgp
 * Please reference to the following paper,
 *                        https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019MS001621
 * NOTE: we use the C++ version of RTE-RRTMGP, which is the implemention of the original Fortran
 * code using C++ YAKL for CUDA, HiP and SYCL application by E3SM ECP team, the C++ version
 * of the rte-rrtmgp code is located at:
 *                       https://github.com/E3SM-Project/rte-rrtmgp
 * The RTE-RRTMGP uses BSD-3-Clause Open Source License, if you want to make changes,
 * and modifications to the code, please refer to BSD-3-Clause Open Source License.
 */
#ifndef ERF_RADIATION_H
#define ERF_RADIATION_H

#include <string>
#include <vector>
#include <memory>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_TableData.H>
#include <AMReX_MultiFabUtil.H>

#include "ERF_Config.H"
#include "ERF_Constants.H"
#include "Rad_constants.H"
#include "Rrtmgp.H"
#include "Optics.H"
#include "Aero_rad_props.H"
#include "Parameterizations.H"
#include "Albedo.H"

// Radiation code interface class
class Radiation {
  public:
    Radiation () {
        // First, make sure yakl has been initialized
        if (!yakl::isInitialized()) yakl::init();
    }

    ~Radiation () = default;

    // init
    void initialize (const amrex::MultiFab& cons_in,
                     const amrex::MultiFab& qmoist,
                     const amrex::BoxArray& grids,
                     const amrex::Geometry& geom,
                     const amrex::Real& dt_advance,
                     const bool& do_sw_rad,
                     const bool& do_lw_rad,
                     const bool& do_aero_rad,
                     const bool& do_snow_opt,
                     const bool& is_cmip6_volcano);

    // run radiation model
    void run ();

    // call back
    void on_complete ();

    void radiation_driver_lw (int ncol, int nlev,
                              const real3d& gas_vmr,
                              const real2d& pmid, const real2d& pint, const real2d& tmid, const real2d& tint,
                              const real3d& cld_tau_gpt, const real3d& aer_tau_bnd,
                              FluxesByband& fluxes_clrsky, FluxesByband& fluxes_allsky,
                              const real2d& qrl, const real2d& qrlc);

    void radiation_driver_sw (int ncol,
                              const real3d& gas_vmr, const real2d& pmid, const real2d& pint, const real2d& tmid,
                              const real2d& albedo_dir, const real2d& albedo_dif, const real1d& coszrs,
                              const real3d& cld_tau_gpt, const real3d& cld_ssa_gpt, const real3d& cld_asm_gpt,
                              const real3d& aer_tau_bnd, const real3d& aer_ssa_bnd, const real3d& aer_asm_bnd,
                              FluxesByband& fluxes_clrsky, FluxesByband& fluxes_allsky,
                              const real2d& qrs, const real2d& qrsc);

    void set_daynight_indices (const real1d& coszrs,
                               const int1d& day_indices,
                               const int1d& night_indices);

    void get_gas_vmr (const std::vector<std::string>& gas_names,
                      const real3d& gas_vmr);

    void calculate_heating_rate (const real2d& flux_up,
                                 const real2d& flux_dn,
                                 const real2d& pint,
                                 const real2d& heating_rate);
  private:
    // geometry
    amrex::Geometry m_geom;

    // valid boxes on which to evolve the solution
    amrex::BoxArray m_box;

    // number of vertical levels
    int nlev, zlo, zhi;

    // number of columns in horizontal plane
    int ncol;

    int nlwgpts, nswgpts;
    int nlwbands, nswbands;

    // radiation options
    bool do_short_wave_rad;
    bool do_long_wave_rad;
    bool do_snow_optics;

    // Flag to indicate whether to do aerosol optical calculations. This
    // zeroes out the aerosol optical properties if False
    bool do_aerosol_rad = true;

    // rrtmgp
    Rrtmgp radiation;

    // optics radation properties
    Optics optics;

    // aerosol optics properties
    AerRadProps aer_rad;

    // input/output of radiation model
    //
    // Net flux calculated in this routine; used to check energy conservation in
    // the physics package driver?
    real1d net_flux;

    // This should be module data or something specific to aerosol where it is used?
    bool is_cmip6_volc;    // true if cmip6 style volcanic file is read otherwise false

    real dt;               // time step(s) - needed for aerosol optics call

    // Surface and top fluxes
    real1d fsns;  // Surface solar absorbed flux
    real1d fsnt;  // Net column abs solar flux at model top
    real1d flns;  // Srf longwave cooling (up-down) flux
    real1d flnt;  // Net outgoing lw flux at model top
    real1d fsds;  // Surface solar down flux

    // radiation data
    const std::vector<std::string> active_gases = {
                                                   "H2O", "CO2", "O3", "N2O",
                                                   "CO", "CH4", "O2", "N2" };

    bool spectralflux  = false;  // calculate fluxes (up and down) per band.

    // Flag to indicate whether or not to use the radiation timestep for solar zenith
    // angle calculations. If true, use the radiation timestep for all solar zenith
    // angle (cosz) calculations.
    // TODO: How does this differ if value is .false.?
    bool use_rad_dt_cosz  = false;

    // Value for prescribing an invariant solar constant (i.e. total solar
    // irradiance at TOA). Used for idealized experiments such as RCE.
    // Disabled when value is less than 0.
    real fixed_total_solar_irradiance = -1.;

    // The RRTMGP warnings are printed when the state variables need to be limted,
    // such as when the temperature drops too low. This is not normally an issue,
    // but in aquaplanet and RCE configurations these situations occur much more
    // frequently, so this flag was added to be able to disable those messages.
    bool rrtmgp_enable_temperature_warnings = true;

    // Output diagnostic brightness temperatures at the top of the
    // atmosphere for 7 TOVS/HIRS channels (2,4,6,8,10,11,12) and 4 TOVS/MSU
    // channels (1,2,3,4).
    // TODO: where are these options set?
    bool dohirs = false;
    int ihirsfq = 1;       // frequency (timesteps) of brightness temperature calcs

    // time step to use for the shr_orb_cosz calculation, if use_rad_dt_cosz set to true
    // TODO: where is this set, and what is shr_orb_cosz? Alternative solar zenith
    // angle calculation? What is the other behavior?
    real dt_avg = 0.0;

    // k-distribution coefficients files to read from. These are set via namelist
    // variables.
    std::string rrtmgp_coefficients_file_sw = "/ccs/home/yuanx/erf.git_develop/Source/Radiation/data/rrtmgp_coefficients_sw_20181204.nc";
    std::string rrtmgp_coefficients_file_lw = "/ccs/home/yuanx/erf.git_develop/Source/Radiation/data/rrtmgp_coefficients_lw_20181204.nc";

    // Band midpoints; these need to be module variables because of how cam_history works;
    // add_hist_coord sets up pointers to these, so they need to persist.
    real1d sw_band_midpoints;
    real1d lw_band_midpoints;

    // rad constituents mixing ratios
    int ngas, naer;
    std::vector<std::string> gasnames;
    std::vector<std::string> aernames;

    int1d rrtmg_to_rrtmgp;

    // Pointers to heating rates on physics buffer
    real2d qrs; // shortwave radiative heating rate
    real2d qrl; // longwave  radiative heating rate

    // Pointers to fields on the physics buffer
    real2d zi;
    real2d clear_rh;

    // Clear-sky heating rates are not on the physics buffer, and we have no
    // reason to put them there, so declare these are regular arrays here
    real2d qrsc;
    real2d qrlc;

    real2d qt, qi, qc, qn;
    real2d tmid, pmid, pdel;
    real2d pint, tint;
    real2d albedo_dir, albedo_dif;
};
#endif // ERF_RADIATION_H
