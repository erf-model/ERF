#ifndef ERF_FIXED_DT_CFL_LIMIT_H_
#define ERF_FIXED_DT_CFL_LIMIT_H_

#include <AMReX_Algorithm.H>
#include <AMReX_REAL.H>

/**
 * Courant number a fixed step may not exceed (erf-model/ERF#5238): min(erf.cfl, 0.8).  An
 * erf.cfl below 0.8 is honored; a larger one is capped at 0.8, a margin below the stability
 * limit.
 *
 * @param[in] cfl  erf.cfl
 */
inline double
FixedDtCFLLimitValue (double cfl)
{
    return amrex::min(cfl, 0.8);
}

/**
 * Factor by which to scale erf.fixed_dt (and erf.fixed_fast_dt) on a level so that the fixed
 * step respects the CFL limit (erf-model/ERF#5238).  The fixed step exceeds the limit when its
 * Courant number fixed_dt * cfl_rate is above FixedDtCFLLimitValue(erf.cfl) = min(erf.cfl, 0.8),
 * and is then scaled to exactly that Courant number, so the scaled step is continuous in the
 * Courant number.  Returns 1 when the fixed step is within the limit, or when there is no fixed
 * step, no CFL limit (cfl <= 0) or no constraining rate on the level.
 *
 * @param[in] fixed_dt  fixed time step of the level
 * @param[in] cfl_rate  inverse time step at a Courant number of one on the level (the rate the
 *                      adaptive dt of estTimeStep is based on)
 * @param[in] cfl       erf.cfl
 */
inline double
FixedDtCFLScale (double fixed_dt, double cfl_rate, double cfl)
{
    if (fixed_dt <= 0.0 || cfl_rate <= 0.0 || cfl <= 0.0) { return 1.0; }
    const double limit   = FixedDtCFLLimitValue(cfl);
    const double courant = fixed_dt * cfl_rate;
    // Relative tolerance, so a step already set to the limit is not scaled again by round-off
    if (courant <= limit * (1.0 + 1.e-10)) { return 1.0; }
    return limit / courant;
}

#endif
