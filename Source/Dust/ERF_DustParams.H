/**
 * @file ERF_DustParams.H
 * @brief Input parameter container for the ERF-Dust 2D surface emission layer.
 *
 * Parameters are read from the "erf.dust" ParmParse prefix at runtime.
 * Default values correspond to bare quartz soil at an arid site.
 *
 * References:
 *   Bagnold, R. A. (1941). The Physics of Blown Sand and Desert Dunes.
 *     Methuen, London.
 *   Marticorena, B., & Bergametti, G. (1995). Modeling the atmospheric dust
 *     cycle: 1. Design of a soil-derived dust emission scheme.
 *     J. Geophys. Res., 100, 16415-16430.
 *     https://doi.org/10.1029/95JD00690
 *   Shao, Y., & Lu, H. (2000). A simple expression for wind erosion threshold
 *     friction velocity.
 *     J. Geophys. Res., 105, 22437-22443.
 *     https://doi.org/10.1029/2000JD900304
 */

#ifndef ERF_DUST_PARAMS_H
#define ERF_DUST_PARAMS_H

#include <AMReX_ParmParse.H>
#include <string>

/**
 * @struct DustParams
 * @brief Container for ERF-Dust module parameters.
 *
 * All parameters are read from the "erf.dust" ParmParse prefix at runtime.
 */
struct DustParams {

    bool        enable              = false;    ///< Enable dust module
    bool        dust_debug          = false;    ///< Write per-timestep debug fields to plotfile
    int         n_size_bins         = 3;        ///< Number of particle size bins
                                                 ///<   Bin 0: coarse  (2.5 um < d <= 10 um)
                                                 ///<   Bin 1: PM10    (1.0 um < d <= 2.5 um)
                                                 ///<   Bin 2: PM2.5   (d <= 1.0 um)
    amrex::Real bin_diameter_um[3]  = {7.0, 3.5, 0.7}; ///< Representative aerodynamic diameter per bin [um]
    amrex::Real particle_density    = 2650.0;  ///< Bulk particle density [kg/m^3]; quartz default
    amrex::Real z0_dust             = 0.01;    ///< Aerodynamic roughness length over bare tailings [m]
    amrex::Real silt_fraction       = 0.10;    ///< Default surface silt mass fraction [dimensionless]
                                                 ///<   Used when no surface map file is provided.
                                                 ///<   Typical range for mine tailings: 0.05-0.30.
                                                 ///<   Reference: Marticorena & Bergametti (1995),
                                                 ///<   J. Geophys. Res., 100, 16415-16430.
                                                 ///<   https://doi.org/10.1029/95JD00690
    amrex::Real threshold_A_coeff  = 0.0123;  ///< Bagnold threshold coefficient A [dimensionless]
                                                 ///<   u*_t = A * sqrt(rho_p * g * d / rho_a)
                                                 ///<   Reference: Bagnold (1941), The Physics of Blown
                                                 ///<   Sand and Desert Dunes, Methuen, London.
    amrex::Real crust_index         = 0.0;     ///< Default surface crust strength index [dimensionless]
                                                 ///<   Range [0,1]: 0 = uncrusted, 1 = fully crusted.
                                                 ///<   Updated by PHREEQC reader in Phase 4.
    amrex::Real phreeqc_update_interval_s = 86400.0; ///< Interval at which PHREEQC output files
                                                 ///<   are re-read [s]. Default = 1 day.
                                                 ///<   Geochemical processes (crust formation,
                                                 ///<   efflorescence) evolve on day-to-week timescales,
                                                 ///<   much slower than the atmospheric timestep.
    std::string surface_map_file    = "";      ///< Path to NetCDF surface property map; empty = uniform
    std::string phreeqc_output_file = "";      ///< Path to PHREEQC output NetCDF; empty = disabled

    /**
     * Constructor that reads all parameters from ParmParse.
     */
    DustParams() {
        amrex::ParmParse pp("erf.dust");

        pp.query("enable", enable);
        pp.query("dust_debug", dust_debug);
        pp.query("n_size_bins", n_size_bins);
        pp.queryarr("bin_diameter_um", bin_diameter_um, n_size_bins);
        pp.query("particle_density", particle_density);
        pp.query("z0_dust", z0_dust);
        pp.query("silt_fraction", silt_fraction);
        pp.query("threshold_A_coeff", threshold_A_coeff);
        pp.query("crust_index", crust_index);
        pp.query("phreeqc_update_interval_s", phreeqc_update_interval_s);
        pp.query("surface_map_file", surface_map_file);
        pp.query("phreeqc_output_file", phreeqc_output_file);
    }
};

#endif // ERF_DUST_PARAMS_H
