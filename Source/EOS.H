#ifndef _EOS_H_
#define _EOS_H_

#include <Constants.H>
#include <AMReX.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE 
amrex::Real getTgivenRTh(amrex::Real density, amrex::Real theta) 
{
    amrex::Real p_loc = p_0 * std::pow(density * R_d * theta / p_0, Gamma);
    return p_loc / (R_d * density); 
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE 
amrex::Real getPgivenRTh(amrex::Real density, amrex::Real theta) 
{
    // diagnostic relation for the full pressure (dry air + water vapor) (arw_v4, Eq. 2.16)
    //   density: density of _dry_ air
    //   theta: _moist_ potential temperature
    return p_0 * std::pow(density * R_d * theta / p_0, Gamma);
}
#endif
