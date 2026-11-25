#ifndef ERF_FLOOD_COUPLING_PARAMS_H
#define ERF_FLOOD_COUPLING_PARAMS_H

#include <string>

/**
 * @file ERF_FloodCouplingParams.H
 * @brief The inputs of the atmosphere side of the flood coupling
 *
 * The erf.flood.* keys that the shared coupling layer (FloodCoupler) reads:
 * whether a flood engine runs, which one, and where its grid sits. Every flood
 * engine shares them; each engine reads its own keys besides (flood.<engine>.*
 * is recommended). read() parses them with their range checks, aborting with a
 * message that names the key.
 *
 * The flood plotfile keys (erf.flood.plot_file, plot_int, plot_per) are read
 * by ERF with its other output keys (ERF::ReadParameters).
 */
struct FloodCouplingParams
{
    /// erf.flood.enable: run a flood engine
    bool        enable       = false;
    /// erf.flood.engine: the flood engine behind the coupling (default: the
    /// only engine compiled in; required when there are several)
    std::string engine;
    /// erf.flood.grid_ratio: flood cells per atmospheric cell in x and y
    int         grid_ratio   = 1;
    /// erf.flood.anchor_level: the level the flood grid refines, -1 = the finest
    int         anchor_level = -1;

    /// Parse and check the keys above from erf.flood.*. With erf.flood.enable,
    /// erf.flood.engine must name an engine compiled into this executable
    /// (flood_engine_names(), ERF_FloodEngineRegistry.H).
    static FloodCouplingParams read ();
};

#endif
