#ifndef ERF_EOS_H_
#define ERF_EOS_H_
#include <ERF_Constants.H>
#include <AMReX.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <cmath>

/**
 * Function to return temperature given density and potential temperatue
 *
 * @params[in] rho density
 * @params[in] rhotheta density times potential temperature theta
*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getTgivenRandRTh(const amrex::Real rho, const amrex::Real rhotheta, const amrex::Real qv=0.0)
{
	// rho and rhotheta are dry values. We should be using moist values with moisture
    amrex::Real p_loc = p_0 * std::pow(R_d * rhotheta * (1.0 + R_v/R_d*qv) * ip_0, Gamma);
    return p_loc / (R_d * rho * (1.0 + R_v/R_d*qv) );
}

/**
 * Function to return potential temperature given density and temperatue
 *
 * @params[in] rho density
 * @params[in] T temperature
 * @params[in] rd0cp ratio of R_d to c_p
*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getThgivenRandT(const amrex::Real rho, const amrex::Real T, const amrex::Real rdOcp, const amrex::Real qv=0.0)
{
    amrex::Real p_loc = rho * R_d * T * (1.0 + R_v/R_d*qv);
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
amrex::Real getdPdRgivenConstantTheta(const amrex::Real rho, const amrex::Real theta)
{
    return Gamma * p_0 * std::pow( (R_d * theta * ip_0), Gamma) * std::pow(rho, Gamma-1.0) ;
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
amrex::Real getExnergivenRTh(const amrex::Real rhotheta, const amrex::Real rdOcp)
{
    // Exner function pi in terms of (rho theta)
    return std::pow(R_d * rhotheta * ip_0, Gamma * rdOcp);
}

/**
 * Function to return the density given pressure
 *
 * @params[in] p pressure
*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getRhoThetagivenP(const amrex::Real p)
{
    // diagnostic relation for the full pressure
    // see https://erf.readthedocs.io/en/latest/theory/NavierStokesEquations.html
    return std::pow(p*std::pow(p_0, Gamma-1), iGamma) * iR_d;
}

#endif

