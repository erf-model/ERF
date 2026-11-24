#ifndef ERF_SURFACE_COUPLING_NAMES_H
#define ERF_SURFACE_COUPLING_NAMES_H

#include <string>

/**
 * @file ERF_SurfaceCouplingNames.H
 * @brief The names a surface coupling uses in its inputs and messages
 *
 * The fire and the flood couplings build on the same surface-coupling layer
 * (Source/SurfaceCoupling): the 2D grid that refines one atmospheric level,
 * the maps to and from its columns, the engine registry, the engine choice
 * and the start-up checks. Each coupling hands the layer its names, so every
 * input and message keeps the coupling's own words: erf.fire.grid_ratio,
 * "[FIRE] ... the fire grid ...", erf.flood.anchor_level, and so on.
 */
struct SurfaceCouplingNames
{
    const char* lower;  ///< "fire": input prefix erf.<lower>, "the <lower> grid", "<lower> engine"
    const char* title;  ///< "Fire": Source/<title>Coupling, ERF_<title>EngineTemplate.H
    const char* upper;  ///< "FIRE": the message tag [<upper>], <upper>_HOME, ERF_<upper>_ENGINE_DIRS

    /// "[FIRE]"
    std::string tag () const { return std::string("[") + upper + "]"; }
    /// "erf.fire." followed by key
    std::string key (const std::string& k) const { return std::string("erf.") + lower + "." + k; }
    /// "fire grid"
    std::string grid () const { return std::string(lower) + " grid"; }
};

#endif
