#ifndef ERF_IBSEB_PARAMS_H
#define ERF_IBSEB_PARAMS_H

#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX.H>
#include <string>

/**
 * \file ERF_IBSEBParams.H
 * \brief Inputs of the immersed-boundary surface energy balance (``erf.ibseb.*``).
 *
 * The balance runs on the faces of buildings represented by immersed forcing
 * (``erf.buildings_type = ImmersedForcing``). Every input lives under the
 * ``erf.ibseb`` prefix and everything is off unless ``erf.ibseb.enable`` is set,
 * in which case ERF::init_ibseb() builds one IBFaceSet per level.
 *
 * The parameters are parsed once, in ERF::ReadParameters(), and copied into
 * each IBFaceSet, so a face set never reads ParmParse itself. New inputs of
 * later phases go here, with their defaults, so that the whole input surface
 * of the balance is in one place. The phases and their inputs are listed in
 * Source/ImmersedBoundarySEB/IBSEB_DEVELOPMENT.md.
 */
struct IBSEBParams
{
    bool        enable        = false;                 ///< Master switch
    int         n_slab_layers = 4;                     ///< Conduction layers per face
    amrex::Real T_skin_init   = 300.0;                 ///< Initial skin temperature [K]
    amrex::Real T_interior    = 293.0;                 ///< Interior (deep) temperature [K]
    std::string csv_file      = "ibseb_buildings.csv"; ///< Per-building report
    int         csv_int       = 100;                   ///< Steps between report rows; <= 0 disables
    bool        debug         = false;                 ///< [IBSEB DEBUG] prints: set description at build, per-step summary, per-building rows
    std::string dump_faces_file;                       ///< Per-rank face dump at report time (empty = off)

    // ---- Radiation provider (phase 2) --------------------------------------
    // Where the downwelling shortwave (and, in phase 3, longwave) comes from.
    // "prescribed" is built in: clear-sky formulas with a fixed or computed
    // sun. "two_stream" is reserved for the ERF-Radiation column fluxes.
    std::string radiation     = "prescribed";
    // "fixed": the sun stays at sun_zenith_deg / sun_azimuth_deg with the
    // direct-normal and diffuse irradiances given below (analytic tests).
    // "solar": sun from the site and time; irradiances from the clear-sky
    // formulas.
    std::string sun_mode      = "fixed";
    amrex::Real sun_zenith_deg    = 45.0;   ///< fixed: zenith angle [deg]
    amrex::Real sun_azimuth_deg   = 180.0;  ///< fixed: azimuth clockwise from north [deg]; 180 = sun in the south
    amrex::Real sw_direct_normal  = 800.0;  ///< fixed: direct-normal irradiance [W/m2]
    amrex::Real sw_diffuse        = 100.0;  ///< fixed: diffuse irradiance on a horizontal surface [W/m2]
    amrex::Real latitude_deg      = 40.0;   ///< solar: site latitude, north positive
    amrex::Real longitude_deg     = -105.0; ///< solar: site longitude, east positive
    amrex::Real day_of_year       = 172.0;  ///< solar: day of year
    amrex::Real utc_offset_hours  = -7.0;   ///< solar: time zone offset from UTC [h]
    amrex::Real time_zero_utc_s   = 12.0 * 3600.0; ///< solar: seconds after 00:00 UTC at simulation time zero
    amrex::Real solar_constant    = 1361.0; ///< solar: [W/m2]
    amrex::Real sw_transmission   = 0.7;    ///< solar: bulk clear-sky transmission (Bird)
    amrex::Real sw_diffuse_coeff  = 0.5;    ///< solar: share of the attenuated beam that arrives as diffuse light
    // Optical properties, uniform until the material library of phase 5.
    amrex::Real albedo            = 0.3;    ///< Face shortwave albedo
    amrex::Real albedo_ground     = 0.2;    ///< Ground albedo for the reflected diffuse term

    /**
     * Read the ``erf.ibseb.*`` inputs and validate them.
     *
     * Only queried, never required: a deck without any ``erf.ibseb`` key gets
     * the defaults above with the balance off. Validation is skipped when the
     * balance is off so an unrelated deck cannot abort on an unused key.
     */
    void init_params ()
    {
        amrex::ParmParse pp("erf.ibseb");
        pp.query("enable",        enable);
        pp.query("n_slab_layers", n_slab_layers);
        pp.query("T_skin_init",   T_skin_init);
        pp.query("T_interior",    T_interior);
        pp.query("csv_file",      csv_file);
        pp.query("csv_int",       csv_int);
        pp.query("debug",         debug);
        pp.query("dump_faces_file", dump_faces_file);
        pp.query("radiation",        radiation);
        pp.query("sun_mode",         sun_mode);
        pp.query("sun_zenith_deg",   sun_zenith_deg);
        pp.query("sun_azimuth_deg",  sun_azimuth_deg);
        pp.query("sw_direct_normal", sw_direct_normal);
        pp.query("sw_diffuse",       sw_diffuse);
        pp.query("latitude_deg",     latitude_deg);
        pp.query("longitude_deg",    longitude_deg);
        pp.query("day_of_year",      day_of_year);
        pp.query("utc_offset_hours", utc_offset_hours);
        pp.query("time_zero_utc_s",  time_zero_utc_s);
        pp.query("solar_constant",   solar_constant);
        pp.query("sw_transmission",  sw_transmission);
        pp.query("sw_diffuse_coeff", sw_diffuse_coeff);
        pp.query("albedo",           albedo);
        pp.query("albedo_ground",    albedo_ground);
        if (enable) {
            if (radiation != "prescribed") {
                amrex::Abort("erf.ibseb.radiation: only \"prescribed\" is available in this build");
            }
            if (sun_mode != "fixed" && sun_mode != "solar") {
                amrex::Abort("erf.ibseb.sun_mode must be \"fixed\" or \"solar\"");
            }
            if (sun_zenith_deg < 0.0 || sun_zenith_deg > 180.0) { amrex::Abort("erf.ibseb.sun_zenith_deg must be in [0, 180]"); }
            if (sw_direct_normal < 0.0 || sw_diffuse < 0.0) { amrex::Abort("erf.ibseb.sw_direct_normal and sw_diffuse must be >= 0"); }
            if (sw_transmission <= 0.0 || sw_transmission > 1.0) { amrex::Abort("erf.ibseb.sw_transmission must be in (0, 1]"); }
            if (albedo < 0.0 || albedo > 1.0 || albedo_ground < 0.0 || albedo_ground > 1.0) {
                amrex::Abort("erf.ibseb.albedo and albedo_ground must be in [0, 1]");
            }
            if (n_slab_layers < 1) { amrex::Abort("erf.ibseb.n_slab_layers must be >= 1"); }
            if (T_skin_init <= 0.0 || T_interior <= 0.0) {
                amrex::Abort("erf.ibseb.T_skin_init and T_interior must be positive [K]");
            }
        }
    }
};

#endif
