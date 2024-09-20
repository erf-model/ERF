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
    amrex::Real beta_interp = -1.0e34;

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
                return beta_interp;
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
        return beta_interp;
    }

    return beta_interp;
}

/*
 * Interpolate between values from a vector of reals (e.g., zlevels_stag) to a
 * new, expanded vector with refine_fac-1 new uniformly spaced points between
 * the original values.
 *
 * Example:
 *   Vector<Real> tmp_orig = {0,1,2};
 *   Vector<Real> tmp(4);
 *   expand_and_interpolate_1d(tmp, tmp_orig, 2, false); // tmp == {0, 0.5, 1, 1.5, 2}
 *   expand_and_interpolate_1d(tmp, tmp_orig, 2, true); // tmp == {0.25, 0.75, 1.25, 1.75}
 */
AMREX_FORCE_INLINE
void
expand_and_interpolate_1d (amrex::Vector<amrex::Real>& znew,
                           const amrex::Vector<amrex::Real>& zorig,
                           int refine_fac,
                           bool destag=false)
{
    const int Nz = zorig.size();
    const int Nz_refined = refine_fac*(Nz-1) + 1;
    if (destag) {
        AMREX_ALWAYS_ASSERT(znew.size() == Nz_refined-1);
    } else {
        AMREX_ALWAYS_ASSERT(znew.size() == Nz_refined);
    }
    amrex::Vector<amrex::Real> znew_stag(Nz_refined);
    for (int k=0; k < Nz; ++k) {
        znew_stag[refine_fac*k] = zorig[k];
    }
    if (refine_fac > 1) {
        for (int k=0; k < Nz-1; ++k) {
            amrex::Real z1 = zorig[k];
            amrex::Real dz = zorig[k+1] - z1;
            for (int koff=1; koff < refine_fac; ++koff) {
                znew_stag[refine_fac*k + koff] = amrex::Real(koff)/refine_fac * dz + z1;
            }
        }
    }
    if (destag) {
        for (int k=0; k < Nz_refined-1; ++k) {
            znew[k] = 0.5*(znew_stag[k] + znew_stag[k+1]);
        }
    } else {
        znew = znew_stag;
    }
}
#endif
