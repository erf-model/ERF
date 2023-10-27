#ifndef Initialization_1D_H
#define Initialization_1D_H
#include <AMReX_REAL.H>

/**
 * Interpolates 1D array beta at the location alpha_interp in the array alpha
 * requiring 1D arrays alpha and beta to be the same size alpha_size.
 * This routine assumes the alpha data are monotonic from lowest to highest.
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real interpolate_1d (const amrex::Real* alpha, const amrex::Real* beta,
                            const amrex::Real alpha_interp, const int alpha_size)
{
    using namespace amrex;

    amrex::Real beta_interp;

    // If the interpolation point already exists in the array,
    // just return it
    for (int i = 0; i < alpha_size; ++i) {
        if (alpha[i] == alpha_interp) {
            beta_interp = beta[i];
            return beta_interp;
        }
    }

    // We didn't find an exact match for alpha_interp in alpha,
    // so we need linear interpolation/extrapolation
    amrex::Real alpha_min = alpha[0];
    amrex::Real alpha_max = alpha[alpha_size-1];
    if (alpha_interp >= alpha_min && alpha_interp <= alpha_max) //interpolate
    {
        for (int i = 0; i < alpha_size; ++i)
        {
            if (alpha_interp >= alpha[i] && alpha_interp <= alpha[i + 1])
            {
                //y = y0 + (y1-y0)*(x-x0)/(x1-x0);
                amrex::Real y0 = beta[i];
                amrex::Real y1 = beta[i + 1];
                amrex::Real x = alpha_interp;
                amrex::Real x0 = alpha[i];
                amrex::Real x1 = alpha[i + 1];
                beta_interp = y0 + (y1 - y0)*(x - x0) / (x1 - x0);
            }
        }
    }
    else //extrapolate
    {
        // y = y0 + ((x - x0) / (x1 - x0)) * (y1 - y0)
        int i = 0;
        if (alpha_interp > alpha_max)
        {
            i = alpha_size-2;
        }
        else if (alpha_interp >= alpha_min)
        {
            amrex::Abort("interpolate_1d: we shouldn't be here!");
        }

        amrex::Real y0 = beta[i];
        amrex::Real y1 = beta[i + 1];
        amrex::Real x = alpha_interp;
        amrex::Real x0 = alpha[i];
        amrex::Real x1 = alpha[i + 1];
        beta_interp = y0 + ((x - x0) / (x1 - x0)) * (y1 - y0);
    }

    return beta_interp;
}
#endif
