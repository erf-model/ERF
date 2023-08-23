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

  // Computes surface albedos over ocean and the surface 
  // Two spectral surface albedos for direct (dir) and diffuse (dif)
  // incident radiation are calculated. The spectral intervals are:
  //   s (shortwave)  = 0.2-0.7 micro-meters
  //   l (longwave)   = 0.7-5.0 micro-meters
  //
  // Uses knowledge of surface type to specify albedo, as follows:
  //
  // Ocean           Uses solar zenith angle to compute albedo for direct
  //                 radiation; diffuse radiation values constant; albedo
  //                 independent of spectral interval and other physical
  //                 factors such as ocean surface wind speed.
  //
  // For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
  // Approximation for Solar Radiation in the NCAR Community Climate Model,
  // Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
 AMREX_FORCE_INLINE 
 void albedo(bool ocean, real1d& coszrs, real1d& tsurface, real1d& asdir, 
              real1d& aldir, real1d& asdif, real1d& aldif)
  {
    using yakl::intrinsics::size;
    const real adif = 0.06;
    int ncol = size(tsurface, 1);
    if (ocean) {
      // Ice-free ocean albedos function of solar zenith angle only, and
      // independent of spectral interval:
      yakl::c::parallel_for(ncol, YAKL_LAMBDA (int icol) {
        if (coszrs(icol) <= 0) {
          aldir(icol) = 0; 
          asdir(icol) = 0; 
          aldif(icol) = 0; 
          asdif(icol) = 0;
        } else if (tsurface(icol) > 271.) {
          aldir(icol) = ( .026 / (std::pow(coszrs(icol), 1.7) + .065)) + 
             (.15*(coszrs(icol) - 0.10) * (coszrs(icol) - 0.50) * (coszrs(icol) - 1.00) );
          asdir(icol) = aldir(icol);
          aldif(icol) = adif;
          asdif(icol) = adif;
        } else {
         // albedo of sea-ice/snow surface
           aldir(icol) = 0.45;
           asdir(icol) = 0.75;
           aldif(icol) = 0.45;
           asdif(icol) = 0.75;
       }
      });
    } else  {// land
       yakl::c::parallel_for(ncol, YAKL_LAMBDA (int icol) {
         if (coszrs(icol) <= 0.) {
           aldir(icol) = 0.; 
           asdir(icol) = 0.; 
           aldif(icol) = 0.; 
           asdif(icol) = 0.;
         } else {
           // Albedos for land type I (Briegleb)
           asdir(icol) = 1.4 * 0.06 / ( 1. + 0.8 * coszrs(icol));
           asdif(icol) = 1.2 * 0.06;
           aldir(icol) = 1.4 * 0.24 / ( 1. + 0.8 * coszrs(icol));
           aldif(icol) = 1.2 * 0.24;
        }
     });
   }
 }
#endif
