//
// parameterization provides the tracer particle properties
// that are not provided by the microphysics model
//
#ifndef ERF_PARAMETERIZATIONS_H
#define ERF_PARAMETERIZATIONS_H_

#include "rrtmgp_const.h"
#include "Rad_constants.H"
#include "ERF_Constants.H"

 AMREX_FORCE_INLINE
 real liquid_effective_radius(real temperature, real landfrac, real icefrac = 0., real snowh = 0.) {
    const real rliqland = 8.0;   // liquid drop size if over land
    const real rliqocean = 14.0; // liquid drop size if over ocean
    const real rliqice = 14.0;   // liquid drop size if over sea ice

    real reffective = 0.;
    // jrm Reworked effective radius algorithm
    // Start with temperature-dependent value appropriate for continental air
    reffective = rliqland+(rliqocean-rliqland)*std::min(1.0, std::max(0.0, (tmelt-temperature)*0.05));

    // Modify for snow depth over land
    reffective = reffective + (rliqocean-reffective)*std::min(1.0, std::max(0.0, snowh*10.));

    // Ramp between polluted value over land to clean value over ocean.
    reffective = reffective + (rliqocean-reffective) * std::min(1.0, std::max(0.0, 1.0 - landfrac));

    // Ramp between the resultant value and a sea ice value in the presence of ice.
    reffective = reffective + (rliqice-reffective) * std::min(1.0, std::max(0.0, icefrac));
    return reffective;
 }

 AMREX_FORCE_INLINE
 real ice_effective_radius(real temperature) {
    int index = int(temperature - (RadConstants::icesize_table_min_temp - 1.));
    index = std::min(std::max(index, 1), 94);
    real fraction = temperature-int(temperature);
    return RadConstants::retab[index]*(1.-fraction)+RadConstants::retab[index+1]*fraction;
 }
#endif
