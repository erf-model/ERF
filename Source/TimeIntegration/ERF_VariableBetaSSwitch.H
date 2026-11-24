#ifndef ERF_VARIABLE_BETA_S_SWITCH_H_
#define ERF_VARIABLE_BETA_S_SWITCH_H_

#include <AMReX_Algorithm.H>
#include <AMReX_REAL.H>

/**
 * Whether the constant beta_s is below the large-Courant stability limit min(sigma^2, 1) of the
 * implicit acoustic substepping somewhere on a level whose largest terrain slope over the w-faces
 * is sigma2_max = max(h_xi^2 + h_eta^2).  ERF::CheckVariableBetaS then turns on
 * erf.use_variable_beta_s (unless the inputs set it).
 *
 * @param[in] beta_s      constant beta_s (erf.beta_s)
 * @param[in] sigma2_max  largest h_xi^2 + h_eta^2 over the w-faces of the level
 */
inline bool
VariableBetaSSwitchNeeded (amrex::Real beta_s, amrex::Real sigma2_max)
{
    return beta_s < amrex::min(sigma2_max, amrex::Real(1.0));
}

#endif
