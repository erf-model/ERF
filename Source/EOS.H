#ifndef ERF_EOS_H_
#define ERF_EOS_H_
#include <ERF_Constants.H>
#include <AMReX.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getTgivenRandRTh(const amrex::Real rho, const amrex::Real rhotheta)
{
    amrex::Real p_loc = p_0 * std::pow(R_d * rhotheta * ip_0, Gamma);
    return p_loc / (R_d * rho);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getThgivenRandT(const amrex::Real rho, const amrex::Real T, const amrex::Real rdOcp)
{
    amrex::Real p_loc = rho * R_d * T;
    return T * std::pow((p_0/p_loc),rdOcp);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getPgivenRTh(const amrex::Real rhotheta, const amrex::Real qv = 0.)
{
#if defined(ERF_USE_MOISTURE)
    amrex::Real R_t        =  R_d + qv*R_v;
    amrex::Real Cp_t       = Cp_d + qv*Cp_v;
    amrex::Real Gamma_t    = Cp_t/(Cp_t-R_t);
    amrex::Real rhotheta_t = rhotheta*(1.0+qv);
    return p_0 * std::pow(R_t * rhotheta_t * ip_0, Gamma_t);
#elif defined(ERF_USE_WARM_NO_PRECIP)
    amrex::Real rhotheta_t = rhotheta*(1.0+(R_v/R_d)*qv);
    return p_0 * std::pow(R_d * rhotheta_t * ip_0, Gamma);
#else
    amrex::ignore_unused(qv);
    return p_0 * std::pow(R_d * rhotheta * ip_0, Gamma);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getRhogivenThetaPress (const amrex::Real theta, const amrex::Real p, const amrex::Real rdOcp)
{
    return std::pow(p_0, rdOcp) * std::pow(p, iGamma) / (R_d * theta);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getdPdRgivenConstantTheta(const amrex::Real rho, const amrex::Real theta)
{
    return Gamma * p_0 * std::pow( (R_d * theta * ip_0), Gamma) * std::pow(rho, Gamma-1.0) ;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getExnergivenP(const amrex::Real P, const amrex::Real rdOcp)
{
    // Exner function pi in terms of P
    return std::pow(P * ip_0, rdOcp);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getExnergivenRTh(const amrex::Real rhotheta, const amrex::Real rdOcp)
{
    // Exner function pi in terms of (rho theta)
    return std::pow(R_d * rhotheta * ip_0, Gamma * rdOcp);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getRhoThetagivenP(const amrex::Real p)
{
    // diagnostic relation for the full pressure
    // see https://erf.readthedocs.io/en/latest/theory/NavierStokesEquations.html
    return std::pow(p*std::pow(p_0, Gamma-1), iGamma) * iR_d;
}

#endif

