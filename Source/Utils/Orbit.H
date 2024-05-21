#ifndef _ORBIT_H_
#define _ORBIT_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <Rrtmgp.H>
#include <ERF_Constants.H>

void
zenith (int& calday,
        amrex::MultiFab* clat,
        amrex::MultiFab* clon,
        real1d& coszrs,
        int& ncol,
        const amrex::Real& eccen,
        const amrex::Real& mvelpp,
        const amrex::Real& lambm0,
        const amrex::Real& obliqr,
        amrex::Real uniform_angle=-1.0);


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
shr_orb_cosz (const amrex::Real& jday,
              const amrex::Real& lat,
              const amrex::Real& lon,
              const amrex::Real& declin,
              amrex::Real uniform_angle)
{
    amrex::Real cos_sol_zen_ang;

    // If uniform angle is specified we use that.
    // NOTE: uniform angle is in degrees
    if (uniform_angle>0.0) {
        cos_sol_zen_ang = std::cos(uniform_angle * PI/180.0);
    }
    // Default method
    else {
        cos_sol_zen_ang = std::sin(lat)*std::sin(declin) - std::cos(lat)*std::cos(declin) *
                          std::cos((jday-std::floor(jday))*2.0*PI + lon);
    }

    return cos_sol_zen_ang;
}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
shr_orb_decl (const amrex::Real& calday,
              const amrex::Real& eccen ,
              const amrex::Real& mvelpp,
              const amrex::Real& lambm0,
              const amrex::Real& obliqr,
              amrex::Real& delta,
              amrex::Real& eccf)
{

    static constexpr amrex::Real day_p_yr = 365.0;
    static constexpr amrex::Real ve = 80.5;

    // Compute eccentricity factor and solar declination using
    // day value where a round day (such as 213.0) refers to 0z at
    // Greenwich longitude.

    // Use formulas from Berger, Andre 1978: Long-Term Variations of Daily
    // Insolation and Quaternary Climatic Changes. J. of the Atmo. Sci.
    // 35:2362-2367.

    // To get the earths true longitude (position in orbit; lambda in Berger
    // 1978) which is necessary to find the eccentricity factor and declination,
    // must first calculate the mean longitude (lambda m in Berger 1978) at
    // the present day.  This is done by adding to lambm0 (the mean longitude
    // at the vernal equinox, set as March 21 at noon, when lambda=0; in radians)
    // an increment (delta lambda m in Berger 1978) that is the number of
    // days past or before (a negative increment) the vernal equinox divided by
    // the days in a model year times the 2*pi radians in a complete orbit.
    amrex::Real lambm = lambm0 + (calday - ve)*2.0*PI/day_p_yr;
    amrex::Real lmm   = lambm  - mvelpp;

    // The earths true longitude, in radians, is then found from
    // the formula in Berger 1978:
    amrex::Real sinl  = std::sin(lmm);
    amrex::Real lamb  = lambm  + eccen*(2.0*sinl + eccen*(1.25*sin(2.0*lmm)
                               + eccen*((13.0/12.0)*sin(3.0*lmm) - 0.25*sinl)));


    // Using the obliquity, eccentricity, moving vernal equinox longitude of
    // perihelion (plus), and earths true longitude, the declination (delta)
    // and the normalized earth/sun distance (rho in Berger 1978; actually inverse
    // rho will be used), and thus the eccentricity factor (eccf), can be
    // calculated from formulas given in Berger 1978.
    amrex::Real invrho = (1.0 + eccen*std::cos(lamb - mvelpp)) / (1.0 - eccen*eccen);

    // Set solar declination and eccentricity factor
    delta  = std::asin(std::sin(obliqr)*std::sin(lamb));
    eccf   = invrho*invrho;
}
#endif
