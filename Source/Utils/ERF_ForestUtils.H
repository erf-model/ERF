#ifndef ERF_FOREST_UTILS_H_
#define ERF_FOREST_UTILS_H_

#include <cmath>
#include <sstream>
#include <string>

#include <AMReX_REAL.H>

namespace erf_forest_utils {

inline std::string
validate_cartesian_coordinates (bool has_x, bool has_y,
                                bool has_lon, bool has_lat,
                                const std::string& source)
{
    if (has_x != has_y) {
        return source + " must provide both projected Cartesian x and y coordinate variables";
    }
    if (!has_x && (has_lon || has_lat)) {
        return source + " provides lon/lat coordinates, but projected Cartesian x/y are required; "
               "no lon/lat-to-Cartesian projection is applied";
    }
    if (!has_x) {
        return source + " must provide projected Cartesian x and y coordinate variables";
    }
    return {};
}

/**
 * Return an actionable diagnostic for an invalid type-2 canopy LAI maximum.
 * The caller supplies the source description so file-backed and input-backed
 * errors can identify where the bad value came from without aborting here.
 */
inline std::string
validate_laimax (const amrex::Real laimax, const std::string& source)
{
    if (std::isfinite(static_cast<double>(laimax)) &&
        laimax >= amrex::Real(0.0) && laimax < amrex::Real(1.0)) {
        return {};
    }

    std::ostringstream message;
    message << source << " must satisfy 0 <= laimax < 1 (got " << laimax << ")";
    return message.str();
}

} // namespace erf_forest_utils

#endif
