#ifndef ERF_Math_H
#define ERF_Math_H
#include <AMReX_REAL.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real interpolate_1d(const amrex::Real* alpha, const amrex::Real* beta,
                           const amrex::Real alpha_interp, const int alpha_size)
{
    using namespace amrex;

    /*
    Interpolates 1D array beta at the location alpha_interp in the array alpha
    requiring 1D arrays alpha and beta to be the same size alpha_size.
    */

    amrex::Real beta_interp;

    //in case the interpolation point alraedy exists in the array
    //just return it
    for (int i = 0; i < alpha_size; ++i) {
        if (alpha[i] == alpha_interp) {
            beta_interp = beta[i];
            return beta_interp;
        }
    }

    // we didn't find an exact match for alpha_interp in alpha,
    // so we need linear interpolation/extrapolation
    const amrex::Real* alpha_begin = &alpha[0];
    const amrex::Real* alpha_end   = &alpha[alpha_size];
    amrex::Real max = *std::max_element(alpha_begin, alpha_end);
    amrex::Real min = *std::min_element(alpha_begin, alpha_end);
    if (alpha_interp >= min && alpha_interp <= max) //interpolate
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
        //y = y0 + ((x - x0) / (x1 - x0)) * (y1 - y0)
        int i;
        if (alpha_interp >= *alpha_end)
        {
            i = static_cast<int>(alpha_end - alpha_begin - 1);
        }
        else
        {
            i = 0;
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
