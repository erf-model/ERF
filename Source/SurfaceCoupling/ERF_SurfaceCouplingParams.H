#ifndef ERF_SURFACE_COUPLING_PARAMS_H
#define ERF_SURFACE_COUPLING_PARAMS_H

#include <ERF_SurfaceCouplingNames.H>
#include <AMReX_ParmParse.H>

#include <string>
#include <vector>

/**
 * @file ERF_SurfaceCouplingParams.H
 * @brief The inputs every surface coupling reads: enable, engine, grid_ratio,
 *        anchor_level
 */

/// "a", "b", "c" (each name in double quotes, comma separated)
std::string surface_quoted_list (const std::vector<std::string>& v);

/**
 * @brief Read erf.<coupling>.enable, grid_ratio, anchor_level and engine
 *
 * @param[in]  pp            ParmParse with the coupling's prefix (erf.fire, erf.flood)
 * @param[in]  names         The coupling's names
 * @param[in]  engines       The engines compiled into this executable, sorted
 * @param[in,out] enable     erf.<coupling>.enable (default on entry)
 * @param[in,out] grid_ratio erf.<coupling>.grid_ratio (default on entry); aborts below 1
 * @param[in,out] anchor_level erf.<coupling>.anchor_level (default on entry); aborts below -1
 * @param[out] engine        The engine: the only one compiled in, or the one the
 *                           deck names. Aborts when the deck enables the coupling
 *                           or names an engine and no engine, or no single one,
 *                           fits.
 */
void read_surface_coupling_common (amrex::ParmParse& pp,
                                   const SurfaceCouplingNames& names,
                                   const std::vector<std::string>& engines,
                                   bool& enable,
                                   int& grid_ratio,
                                   int& anchor_level,
                                   std::string& engine);

#endif
