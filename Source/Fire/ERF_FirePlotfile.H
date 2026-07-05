#ifndef ERF_FIRE_PLOTFILE_H
#define ERF_FIRE_PLOTFILE_H

#include <string>

// Forward declarations
class FireLayer;

namespace amrex {
    class Real;
}

/**
 * @file ERF_FirePlotfile.H
 * @brief Fire plotfile writer for output on native fire grid
 *
 * Writes fire state variables (phi, ros, wind, slopes, fuel moisture) to
 * AMReX plotfile format on the native fire grid (C x finer than atmospheric).
 */

/**
 * @brief Write a fire plotfile at the current simulation time.
 *
 * Gathers fire state from FireLayer, assembles into a single MultiFab in
 * canonical order, and writes to AMReX plotfile format on the fire grid.
 * Includes a JSON sidecar with metadata.
 *
 * @param[in] plotfile_prefix  Base name for plotfile (e.g., "plt_fire_")
 * @param[in] fire_layer       Reference to active FireLayer
 * @param[in] time             Current simulation time [s]
 * @param[in] step             Current timestep number
 *
 * All output is on the native fire grid geometry.
 * Variable order is defined by fire_plotfile_var_names().
 */
void WriteFirePlotfile(const std::string& plotfile_prefix,
                       const FireLayer&   fire_layer,
                       amrex::Real        time,
                       int                step);

#endif
