#ifndef ERF_FIXED_DT_CFL_LIMIT_H_
#define ERF_FIXED_DT_CFL_LIMIT_H_

#include <AMReX_Algorithm.H>
#include <AMReX_REAL.H>

/**
 * Factor by which to scale erf.fixed_dt (and erf.fixed_fast_dt) on a level so that the fixed
 * step respects the CFL limit (erf-model/ERF#5238).  The fixed step exceeds the limit when its
 * Courant number fixed_dt * cfl_rate is above erf.cfl; it is then scaled to the Courant number
 * min(erf.cfl, 0.8).  Returns 1 when the fixed step is within the limit, or when there is no
 * fixed step, no CFL limit (cfl <= 0) or no constraining rate on the level.
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
    const double courant = fixed_dt * cfl_rate;
    // Relative tolerance, so a step already set to the limit is not scaled again by round-off
    if (courant <= cfl * (1.0 + 1.e-10)) { return 1.0; }
    return amrex::min(cfl, 0.8) / courant;
}

#endif
