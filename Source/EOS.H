#ifndef ERF_EOS_H_
#define ERF_EOS_H_
#include <ERF_Constants.H>
#include <AMReX.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <cmath>

/**
 * Function to return potential temperature given pressure and temperature
 *
 * @params[in] pressure
 * @params[in] temperature
 * @params[in] rd0cp ratio of R_d to c_p
*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getThgivenPandT(const amrex::Real T, const amrex::Real P, const amrex::Real rdOcp)
{
    return T*std::pow(p_0/P, rdOcp);
}

/**
 * Function to return temperature given density and potential temperature
 *
 * @params[in] rho density
 * @params[in] rhotheta density times potential temperature theta
*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getTgivenRandRTh(const amrex::Real rho, const amrex::Real rhotheta, const amrex::Real qv=0.0)
{
    // rho and rhotheta are dry values. We should be using moist value of theta when using moisture
    // theta_m = theta * (1 + R_v/R_d*qv)
    amrex::Real p_loc = p_0 * std::pow(R_d * rhotheta * (1.0 + R_v/R_d*qv) * ip_0, Gamma);
    // p = rho_d * R_d * T_v (not T)
    // T_v = T * (1 + R_v/R_d*qv)
    return p_loc / (R_d * rho * (1.0 + R_v/R_d*qv) );
}

/**
 * Function to return potential temperature given density and temperature
 *
 * @params[in] rho density
 * @params[in] T temperature
 * @params[in] rd0cp ratio of R_d to c_p
*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getThgivenRandT(const amrex::Real rho, const amrex::Real T, const amrex::Real rdOcp, const amrex::Real qv=0.0)
{
    // p = rho_d * R_d * T_moist
    amrex::Real p_loc = rho * R_d * T * (1.0 + R_v/R_d*qv);
    // theta_d = T * (p0/p)^(R_d/C_p)
    return T * std::pow((p_0/p_loc),rdOcp);
}

/**
 * Function to return pressure given density times theta
 *
 * @params[in] rhotheta density times potential temperature
 * @params[in] qv water vapor
*/
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

/**
 * Function to return density given theta and pressure
 *
 * @params[in] theta potential temperature
 * @params[in] p pressure
 * @params[in] rd0cp ratio of R_d to c_p
*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getRhogivenThetaPress (const amrex::Real theta, const amrex::Real p, const amrex::Real rdOcp, const amrex::Real qv=0.0)
{
    // We should be using moist value of theta when using moisture
    // theta_m = theta * (1 + R_v/R_d*qv)
    return std::pow(p_0, rdOcp) * std::pow(p, iGamma) / (R_d * theta * (1.0 + R_v/R_d*qv) );
}

/**
 * Function to return dP/drho at constant theta
 *
 * @params[in] theta potential temperature
 * @params[in] p pressure
 * @params[in] rd0cp ratio of R_d to c_p
*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getdPdRgivenConstantTheta(const amrex::Real rho, const amrex::Real theta, const amrex::Real qv=0.0)
{
    // We should be using moist value of theta when using moisture
    // theta_m = theta * (1 + R_v/R_d*qv)
    return Gamma * p_0 * std::pow( (R_d * theta * (1.0 + R_v/R_d*qv) * ip_0), Gamma) * std::pow(rho, Gamma-1.0) ;
}

/**
 * Function to return the Exner function pi given pressure
 * @params[in] p pressure
 * @params[in] rd0cp ratio of R_d to c_p
*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getExnergivenP(const amrex::Real P, const amrex::Real rdOcp)
{
    // Exner function pi in terms of P
    return std::pow(P * ip_0, rdOcp);
}

/**
 * Function to return the Exner function pi given densith times potential temperature
 *
 * @params[in] rhotheta density times potential temperature
 * @params[in] rd0cp ratio of R_d to c_p
*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getExnergivenRTh(const amrex::Real rhotheta, const amrex::Real rdOcp, const amrex::Real qv=0.0 )
{
    // Exner function pi in terms of (rho theta)
    // We should be using moist value of theta when using moisture
    // theta_m = theta * (1 + R_v/R_d*qv)
    return std::pow(R_d * rhotheta *  (1.0 + R_v/R_d*qv) * ip_0, Gamma * rdOcp);
}

/**
 * Function to return the density given pressure
 *
 * @params[in] p pressure
*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getRhoThetagivenP(const amrex::Real p, const amrex::Real qv=0.0)
{
    // diagnostic relation for the full pressure
    // see https://erf.readthedocs.io/en/latest/theory/NavierStokesEquations.html
    // For cases with mositure theta = theta_m / (1 + R_v/R_d*qv)
    return std::pow(p*std::pow(p_0, Gamma-1), iGamma) * iR_d / (1.0 + R_v/R_d*qv) ;
}

#endif

