#ifndef ERF_FIRE_COUPLING_PARAMS_H
#define ERF_FIRE_COUPLING_PARAMS_H

#include <AMReX_REAL.H>
#include <string>

/**
 * @file ERF_FireCouplingParams.H
 * @brief The inputs of the atmosphere side of the fire coupling
 *
 * The erf.fire.* keys that the shared coupling layer (FireCoupler) reads:
 * whether the fire runs and which engine, where its grid sits, when and how
 * its fluxes enter the atmosphere, and the smoke tracer. Every fire engine
 * shares them; each engine reads its own keys besides. read() parses them
 * with their range checks, aborting with a message that names the key.
 */
struct FireCouplingParams
{
    /// erf.fire.enable: run a fire
    bool        enable                = false;
    /// erf.fire.engine: the fire engine behind the coupling (default: the only
    /// engine compiled in; required when there are several)
    std::string engine;
    /// erf.fire.grid_ratio: fire cells per atmospheric cell in x and y
    int         grid_ratio            = 5;
    /// erf.fire.anchor_level: the level the fire grid refines, -1 = the finest
    int         anchor_level          = -1;
    /// erf.fire.coupling_type: passive, lagged or synchronous
    std::string coupling_type         = "lagged";
    /// erf.fire.heat_flux_alfg: e-folding height of the injected heat [m]
    amrex::Real heat_flux_alfg        = 45.0;
    /// erf.fire.fire_atm_feedback: multiplier on the injected fluxes
    amrex::Real fire_atm_feedback     = 1.0;
    /// erf.fire.inject_latent: inject the engine's latent flux as vapour
    bool        inject_latent         = true;
    /// erf.fire.source_mode: overwrite or add to the theta and vapour sources
    std::string source_mode           = "overwrite";
    /// erf.fire.source_mode was given in the deck (the default overwrite must be chosen explicitly
    /// when another theta source would be discarded)
    bool        source_mode_given     = false;
    /// erf.fire.heat_open_fraction: place the heat by the engine's open fraction
    bool        heat_open_fraction    = false;
    /// erf.fire.heat_tendency_density: density factor of the theta tendency
    bool        heat_tendency_density = false;
    /// erf.fire.heat_tendency_exner: Exner factor of the theta tendency
    bool        heat_tendency_exner   = true;
    /// erf.fire.heat_profile: exponential (WRF-SFIRE, e-folding heat_flux_alfg)
    /// or bottom_cell (all of the flux into the lowest cell)
    std::string heat_profile          = "exponential";
    /// erf.fire.smoke_enable: inject the smoke tracer
    bool        smoke_enable          = false;
    /// erf.fire.smoke_emission_factor: kg smoke per kg fuel burned [-]
    amrex::Real smoke_emission_factor = 0.02;
    /// erf.fire.fire_debug: debug output of the coupling
    bool        fire_debug            = false;

    bool is_passive()     const { return coupling_type == "passive"; }
    bool is_lagged()      const { return coupling_type == "lagged"; }
    bool is_synchronous() const { return coupling_type == "synchronous"; }
    /// Heat (and vapour) reach the atmosphere: lagged and synchronous coupling.
    bool injects_flux()   const { return is_lagged() || is_synchronous(); }
    bool heat_in_bottom_cell() const { return heat_profile == "bottom_cell"; }

    /// Parse and check the keys above from erf.fire.*. With erf.fire.enable,
    /// erf.fire.engine must name an engine compiled into this executable
    /// (fire_engine_names(), ERF_FireEngineRegistry.H).
    static FireCouplingParams read ();
};

#endif
