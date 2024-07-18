#include <NumericalDiffusion.H>

using namespace amrex;

/**
 * Function to compute 6th order numerical diffusion RHS.
 *
 * @param[in]  bx box to loop over
 * @param[in]  start_comp staring component index
 * @param[in]  num_comp number of total components
 * @param[in]  num_diff_coeff
 * @param[in]  data variables used to compute RHS
 * @param[out] rhs store the right hand side
 * @param[in]  mf_x map factor at x-face
 * @param[in]  mf_y map factor at y-face
 * @param[in]  avg_mf_x_y flag to average map factor x in y-dir
 * @param[in]  avg_mf_y_x flag to average map factor y in x-dir
 */
void
NumericalDiffusion (const Box& bx,
                    const int  start_comp,
                    const int  num_comp,
                    const Real dt,
                    const Real num_diff_coeff,
                    const Array4<const Real>& data,
                    const Array4<      Real>& rhs,
                    const Array4<const Real>& mf_x,
                    const Array4<const Real>& mf_y,
                    const bool avg_mf_x_y,
                    const bool avg_mf_y_x)
{
    BL_PROFILE_VAR("NumericalDiffusion()",NumericalDiffusion);

    // Average map factors to correct locations
    Box planebx(bx); planebx.setSmall(2,0); planebx.setBig(2,0);
    FArrayBox mf_x_bar; mf_x_bar.resize(planebx,1,The_Async_Arena());
    FArrayBox mf_y_bar; mf_y_bar.resize(planebx,1,The_Async_Arena());
    const Array4<Real>& mfx_arr  = mf_x_bar.array();
    const Array4<Real>& mfy_arr  = mf_y_bar.array();
    ParallelFor(planebx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (avg_mf_x_y) {
            mfx_arr(i,j,k) = 0.5 * ( mf_x(i,j-1,k) + mf_x(i,j,k) );
        } else {
            mfx_arr(i,j,k) = mf_x(i,j,k);
        }
        if (avg_mf_y_x) {
            mfy_arr(i,j,k) = 0.5 * ( mf_y(i-1,j,k) + mf_y(i,j,k) );
        } else {
            mfy_arr(i,j,k) = mf_y(i,j,k);
        }
    });

    // Capture diffusion coeff
    Real coeff6 = num_diff_coeff / (2.0 * dt);

    // Compute 5th order derivative and augment RHS
    ParallelFor(bx, num_comp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int m) noexcept
    {
        int n = start_comp + m;
        Real xflux_lo = 10. * (data(i  ,j,k,n) - data(i-1,j,k,n))
                       - 5. * (data(i+1,j,k,n) - data(i-2,j,k,n))
                            + (data(i+2,j,k,n) - data(i-3,j,k,n));
        if ( (xflux_lo * (data(i,j,k,n) - data(i-1,j,k,n)) ) < 0.) xflux_lo = 0.;
        Real xflux_hi = 10. * (data(i+1,j,k,n) - data(i  ,j,k,n))
                       - 5. * (data(i+2,j,k,n) - data(i-1,j,k,n))
                            + (data(i+3,j,k,n) - data(i-2,j,k,n));
        if ( (xflux_hi * (data(i+1,j,k,n) - data(i,j,k,n)) ) < 0.) xflux_hi = 0.;
        Real yflux_lo = 10. * (data(i,j  ,k,n) - data(i,j-1,k,n))
                       - 5. * (data(i,j+1,k,n) - data(i,j-2,k,n))
                            + (data(i,j+2,k,n) - data(i,j-3,k,n));
        if ( (yflux_lo * (data(i,j,k,n) - data(i,j-1,k,n)) ) < 0.) yflux_lo = 0.;
        Real yflux_hi = 10. * (data(i,j+1,k,n) - data(i,j  ,k,n))
                       - 5. * (data(i,j+2,k,n) - data(i,j-1,k,n))
                            + (data(i,j+3,k,n) - data(i,j-2,k,n));
        if ( (yflux_hi * (data(i,j+1,k,n) - data(i,j,k,n)) ) < 0.) yflux_hi = 0.;
        rhs(i,j,k,n) += coeff6 * ( (xflux_hi - xflux_lo) * mfx_arr(i,j,0)
                                 + (yflux_hi - yflux_lo) * mfy_arr(i,j,0) );
    });
}

/**
 * Function to compute 6th order numerical diffusion RHS.
 *
 * @param[in]  bx box to loop over
 * @param[in]  start_comp staring component index
 * @param[in]  num_comp number of total components
 * @param[in]  num_diff_coeff
 * @param[in]  data variables used to compute RHS
 * @param[out] rhs store the right hand side
 */
void
NumericalDiffusionVert (const Box& bx,
                        const int  start_comp,
                        const int  num_comp,
                        const Real dt,
                        const Real num_diff_coeff,
                        const Array4<const Real>& data,
                        const Array4<      Real>& rhs)
{
    BL_PROFILE_VAR("NumericalDiffusionVert()",NumericalDiffusionVert);

    // Capture diffusion coeff
    Real coeff6 = num_diff_coeff / (2.0 * dt);

    // Compute 5th order derivative and augment RHS
    ParallelFor(bx, num_comp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int m) noexcept
    {
        int n = start_comp + m;
        Real zflux_lo = 10. * (data(i,j,k  ,n) - data(i,j,k-1,n))
                       - 5. * (data(i,j,k+1,n) - data(i,j,k-2,n))
                            + (data(i,j,k+2,n) - data(i,j,k-3,n));
        if ( (zflux_lo * (data(i,j,k,n) - data(i,j,k-1,n)) ) < 0.) zflux_lo = 0.;
        Real zflux_hi = 10. * (data(i,j,k+1,n) - data(i,j,k  ,n))
                       - 5. * (data(i,j,k+2,n) - data(i,j,k-1,n))
                            + (data(i,j,k+3,n) - data(i,j,k-2,n));
        if ( (zflux_hi * (data(i,j,k+1,n) - data(i,j,k,n)) ) < 0.) zflux_hi = 0.;
        rhs(i,j,k,n) += coeff6 * ( (zflux_hi - zflux_lo) );
    });
}
