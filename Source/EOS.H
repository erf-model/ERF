#ifndef _EOS_H_
#define _EOS_H_

#include <Constants.H>
#include <AMReX.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE 
amrex::Real getTgivenRTh(amrex::Real density, amrex::Real theta) 
{
    amrex::Real p_loc = std::pow(density * R_d * theta / p_0, Gamma) ;
    return p_loc / (R_d * density); 
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE 
amrex::Real getPgivenRTh(amrex::Real density, amrex::Real theta) 
{
    return std::pow(density * R_d * theta / p_0, Gamma) ;
}
#endif
