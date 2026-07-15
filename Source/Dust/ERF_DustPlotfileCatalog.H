/**
 * @file ERF_DustPlotfileCatalog.H
 * @brief Variable name list for dust 2D plotfile output.
 *
 * Pattern follows ERF_FirePlotfileCatalog.H in Source/Fire/:
 *   https://github.com/hgopalan/ERF/blob/ERF-Fire/Source/Fire/ERF_FirePlotfileCatalog.H
 */
#ifdef ERF_USE_DUST
#include <AMReX_Vector.H>
#include <string>

inline amrex::Vector<std::string> dust_plotfile_var_names()
{
    return {
        "dust_emission_flux",  ///< Total emission flux [kg/m^2/s]. Marticorena & Bergametti (1995).
        "dust_ustar_in",       ///< Friction velocity from SurfaceLayer [m/s].
        "dust_ustar_t",        ///< Threshold friction velocity [m/s]. Shao & Lu (2000).
        "dust_deposition_rate",///< Accumulated dry deposition [kg/m^2]. Zhang et al. (2001).
        "dust_conc_sfc",       ///< Near-surface dust concentration [kg/m^3].
        "dust_surf_moist",     ///< Surface moisture flux [kg/m^2/s]. Null when no moisture scheme.
        "dust_suppression",    ///< Suppression agent coverage [0-1].
        "dust_retreat_flag",   ///< Re-treatment flag [0/1].
        // Phase 17: EPA NAAQS fields
        "dust_pm25_ug_m3",      ///< Instantaneous PM2.5 [µg/m^3].
        "dust_pm10_ug_m3",      ///< Instantaneous PM10 [µg/m^3].
        "dust_pm25_24h_ug_m3",  ///< 24-hour running average PM2.5 [µg/m^3].
        "dust_pm10_24h_ug_m3",  ///< 24-hour running average PM10 [µg/m^3].
        "dust_pm25_exceed",     ///< Exceedance flag: 1 where PM2.5 24h > 35 µg/m^3.
        "dust_pm10_exceed",     ///< Exceedance flag: 1 where PM10 24h > 150 µg/m^3.
        // Phase 18: MSHA worker exposure tracking
        "dust_msha_dose_mg_m3_h",  ///< Cumulative 8h TWA dose [mg/m^3/h].
        "dust_msha_twa_mg_m3",     ///< 8-hour time-weighted average [mg/m^3]. 30 CFR 56.5001.
        "dust_msha_exceed",        ///< PEL exceedance flag [0/1].
        "dust_msha_shift_twa",     ///< TWA at last shift end [mg/m^3].
        // Phase 20: Multi-mine site index
        "dust_site_id"             ///< Mine site index [0=unassigned, 1..N].
    };
}

inline int dust_plotfile_ncomp() { return 20; }
#endif

