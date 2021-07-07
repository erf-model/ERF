#ifndef _EOS_H_
#define _EOS_H_

#include <Constants.H>
#include <AMReX.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE 
amrex::Real getTgivenRandRTh(amrex::Real rho, amrex::Real rhotheta) 
{
    amrex::Real p_loc = p_0 * std::pow(R_d * rhotheta / p_0, Gamma);
    return p_loc / (R_d * rho); 
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE 
amrex::Real getPgivenRTh(amrex::Real rhotheta) 
{
    // diagnostic relation for the full pressure
    // see https://erf.readthedocs.io/en/latest/Euler_Equations.html
    return p_0 * std::pow(R_d * rhotheta / p_0, Gamma);
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getRhoThetagivenP(amrex::Real p)
{
    // diagnostic relation for the full pressure
    // see https://erf.readthedocs.io/en/latest/Euler_Equations.html
    return std::pow(p*std::pow(p_0, Gamma-1), 1.0 / Gamma)/ R_d;
}
#endif
