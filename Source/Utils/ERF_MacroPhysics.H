#ifndef ERF_MACROPHYS_H
#define ERF_MACROPHYS_H

#include "ERF_MicrophysicsUtils.H"

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
static void
NewtonIterSat (const amrex::Real& fac_cond,
               amrex::Real& tabs,
               const amrex::Real& pres,
               amrex::Real& qv)
{
    // Solution tolerance
    amrex::Real tol = amrex::Real(1.0e-8);

    // Saturation moisture fractions
    amrex::Real qsat, dqsat;

    // Newton iteration vars
    int  niter;
    amrex::Real fff, dfff;
    amrex::Real dtabs;

    // Initial guess for temperature
    amrex::Real tabs_old = tabs;

    niter = 0;
    dtabs = 1;

    //==================================================
    // Newton iteration to qv=qsat (cloud phase only)
    //==================================================
    do {
        // Saturation moisture fractions
        erf_qsatw(tabs, pres, qsat);
        erf_dtqsatw(tabs, pres, dqsat);

        // Function for root finding:
        // 0 = -T_new + T_old + L_v/C_p * (qv - qsat)
        fff   = -tabs + tabs_old + fac_cond*(qv - qsat);

        // Derivative of function (T_new iterated on)
        dfff  = -one - fac_cond*dqsat;

        // Update the temperature
        dtabs = -fff/dfff;
        tabs += dtabs;

        // Update iteration
        niter = niter+1;
    } while(std::abs(dtabs) > tol && niter < 20);

    // Update from last iteration (dq = dq/dt * dt)
    qv  = qsat + dqsat*dtabs;
}
#endif
