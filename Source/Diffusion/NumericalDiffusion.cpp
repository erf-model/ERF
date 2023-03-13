#include <NumericalDiffusion.H>

using namespace amrex;

void
NumericalDiffusion (const Box& bx,
                    const int  n_start,
                    const int  n_end,
                    const Real dt,
                    const Array4<const Real>& data,
                    const Array4<      Real>& rhs,
                    const SolverChoice& solverChoice)
{
    BL_PROFILE_VAR("NumericalDiffusion()",NumericalDiffusion);

    // Set indices for components
    const int ncomp = n_end - n_start + 1;

    // Capture diffusion coeff
    Real coeff6 = solverChoice.NumDiffCoeff / (2.0 * dt);

    // Compute 5th order derivative and augment RHS
    amrex::ParallelFor(bx,ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int m) noexcept
    {
        int n = n_start + m;
        Real xflux_lo = 10. * (data(i  ,j,k,n) - data(i-1,j,k,n))
                       - 5. * (data(i+1,j,k,n) - data(i-2,j,k,n))
                            + (data(i+2,j,k,n) - data(i-3,j,k,n));
        if ( (xflux_lo * (data(i,j,k,n) - data(i-1,j,k,n)) ) > 0.) xflux_lo = 0.;
        Real xflux_hi = 10. * (data(i+1,j,k,n) - data(i  ,j,k,n))
                       - 5. * (data(i+2,j,k,n) - data(i-1,j,k,n))
                            + (data(i+3,j,k,n) - data(i-2,j,k,n));
        if ( (xflux_hi * (data(i+1,j,k,n) - data(i,j,k,n)) ) > 0.) xflux_hi = 0.;
        Real yflux_lo = 10. * (data(i,j  ,k,n) - data(i,j-1,k,n))
                       - 5. * (data(i,j+1,k,n) - data(i,j-2,k,n))
                            + (data(i,j+2,k,n) - data(i,j-3,k,n));
        if ( (yflux_lo * (data(i,j,k,n) - data(i,j-1,k,n)) ) > 0.) yflux_lo = 0.;
        Real yflux_hi = 10. * (data(i,j+1,k,n) - data(i,j  ,k,n))
                       - 5. * (data(i,j+2,k,n) - data(i,j-1,k,n))
                            + (data(i,j+3,k,n) - data(i,j-2,k,n));
        if ( (yflux_hi * (data(i,j+1,k,n) - data(i,j,k,n)) ) > 0.) yflux_hi = 0.;
        rhs(i,j,k,n) += coeff6 * (xflux_hi - xflux_lo
                                + yflux_hi - yflux_lo);
    });
}
