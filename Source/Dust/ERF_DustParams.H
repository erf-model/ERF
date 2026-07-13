/**
 * @file ERF_DustParams.H
 * @brief Input parameter container for the ERF-Dust 2D surface emission layer.
 *
 * Parameters are read from the "erf.dust" ParmParse prefix at runtime.
 * Default values correspond to bare quartz soil at an arid site.
 *
 * Phase 1 creates this struct and registers it with CMake and GNUmake.
 * Phase 2 adds grid_ratio and uses this struct in ERF_DustGrid and ERF_DustLayer.
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
#include <AMReX_Vector.H>
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
    int         n_size_bins         = 3;        ///< Number of particle size bins; each bin maps to
                                                 ///<   one RhoQ*_comp scalar component in ERF.
                                                 ///<   Bin 0: coarse  (2.5 um < d <= 10 um)
                                                 ///<   Bin 1: PM10    (1.0 um < d <= 2.5 um)
                                                 ///<   Bin 2: PM2.5   (d <= 1.0 um)
    int         grid_ratio          = 1;        ///< Refinement factor of dust grid relative to ERF
                                                 ///<   level-0 horizontal grid [dimensionless].
                                                 ///<   grid_ratio = 1 places one dust cell per
                                                 ///<   atmospheric cell. Values > 1 subdivide each
                                                 ///<   atmospheric cell into grid_ratio^2 dust cells.
                                                 ///<   All atmospheric box x,y sizes must be divisible
                                                 ///<   by grid_ratio. Set erf.dust.grid_ratio to match
                                                 ///<   the desired surface resolution.
    amrex::Vector<amrex::Real> bin_diameter_um  = {7.0, 3.5, 0.7}; ///< Representative aerodynamic diameter per bin [um]
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
    amrex::Real alpha_crust     = 0.5;  ///< Crust strength reduction coefficient for u*_t [-].
                                         ///<   u*_t_new = u*_t_base * (1 - alpha_crust * crust_index).
                                         ///<   Reference: Marticorena & Bergametti (1995),
                                         ///<     https://doi.org/10.1029/95JD00690
    amrex::Real alpha_efflor    = 0.3;  ///< Efflorescence reduction coefficient for u*_t [-].
    amrex::Real ustar_t_base    = -1.0; ///< Base u*_t before chemistry [m/s].
                                         ///<   Negative means compute from Bagnold at init.
    amrex::Real test_ustar      = 0.0;   ///< Placeholder u* [m/s] used in Phase 6 before
                                         ///<   Phase 9 extracts u* from the MRF surface layer.
                                         ///<   Set > 0 to test emission physics.
                                         ///<   0 = no emission (all cells below threshold).
    amrex::Real rho_air         = 1.225; ///< Air density used in emission flux [kg/m^3].
                                         ///<   Standard atmosphere default.
    std::string phreeqc_crust_var    = "crust_index";    ///< Variable/column name for crust index.
    std::string phreeqc_silt_var     = "silt_fraction";  ///< Variable/column name for silt fraction.
    std::string phreeqc_efflor_var   = "efflorescence";  ///< Variable/column name for efflorescence.
    std::string phreeqc_supp_var     = "suppression_mod";///< Variable/column name for suppression modifier.
    std::string phreeqc_metal_var    = "metal_as_bin0";  ///< Variable/column name for As mass fraction bin 0.
    std::string surface_map_file    = "";      ///< Path to NetCDF surface property map; empty = uniform
    std::string phreeqc_output_file = "";      ///< Path to PHREEQC output NetCDF; empty = disabled
    std::string soil_type_file      = "";      ///< Path to soil type raster (.asc or .nc).
                                                 ///<   Code 0 = undefined; 1-16 = STATSGO natural soils;
                                                 ///<   100 = mine tailings; 101 = Li brine pond;
                                                 ///<   102 = REE tailings; 103 = Cu tailings;
                                                 ///<   104 = unpaved haul road.
                                                 ///<   STATSGO: https://www.nrcs.usda.gov/resources/
                                                 ///<     data-and-reports/ssurgo/statsgo2-data
                                                 ///<   AP-42 road emission: https://www.epa.gov/
                                                 ///<     air-emissions-factors-and-quantification/
                                                 ///<     ap-42-compilation-air-emission-factors
    std::string silt_fraction_file  = "";      ///< Path to silt fraction raster. Empty = uniform silt_fraction.
    std::string crust_index_file    = "";      ///< Path to crust index raster. Empty = uniform crust_index.
    std::string moisture_flag_file  = "";      ///< Path to moisture inhibition flag raster. Empty = uniform 0.
    std::string suppression_file    = "";      ///< Path to suppression coverage raster. Empty = uniform 0.
                                                 ///<   Supported formats: ESRI ASCII (.asc), NetCDF (.nc).
                                                 ///<   NetCDF requires ERF_ENABLE_NETCDF=ON.

    /**
     * Constructor that reads all parameters from ParmParse.
     */
    DustParams() {
        amrex::ParmParse pp("erf.dust");

        pp.query("enable", enable);
        pp.query("dust_debug", dust_debug);
        pp.query("n_size_bins", n_size_bins);
        pp.queryarr("bin_diameter_um", bin_diameter_um);
        pp.query("particle_density", particle_density);
        pp.query("z0_dust", z0_dust);
        pp.query("silt_fraction", silt_fraction);
        pp.query("threshold_A_coeff", threshold_A_coeff);
        pp.query("crust_index", crust_index);
        pp.query("grid_ratio", grid_ratio);
        pp.query("phreeqc_update_interval_s", phreeqc_update_interval_s);
        pp.query("alpha_crust", alpha_crust);
        pp.query("alpha_efflor", alpha_efflor);
        pp.query("ustar_t_base", ustar_t_base);
        pp.query("test_ustar", test_ustar);
        pp.query("rho_air", rho_air);
        pp.query("phreeqc_crust_var", phreeqc_crust_var);
        pp.query("phreeqc_silt_var", phreeqc_silt_var);
        pp.query("phreeqc_efflor_var", phreeqc_efflor_var);
        pp.query("phreeqc_supp_var", phreeqc_supp_var);
        pp.query("phreeqc_metal_var", phreeqc_metal_var);
        pp.query("surface_map_file", surface_map_file);
        pp.query("phreeqc_output_file", phreeqc_output_file);
        pp.query("soil_type_file",     soil_type_file);
        pp.query("silt_fraction_file", silt_fraction_file);
        pp.query("crust_index_file",   crust_index_file);
        pp.query("moisture_flag_file", moisture_flag_file);
        pp.query("suppression_file",   suppression_file);
    }
};

#endif // ERF_DUST_PARAMS_H