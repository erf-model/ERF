#include <ERF_NumericalDiffusion.H>

using namespace amrex;

/**
 * Function to compute 6th order numerical diffusion RHS.
 *
 * @param[in]  bx box to loop over
 * @param[in]  start_comp staring component index
 * @param[in]  num_comp number of total components
 * @param[in]  num_diff_coeff
 * @param[in]  prim_data primitive variables
 * @param[in]  cell_data cell center variables
 * @param[out] rhs store the right hand side
 * @param[in]  mf map factor
 */
void
NumericalDiffusion_Scal (const Box& bx,
                         const int  start_comp,
                         const int  num_comp,
                         const Real dt,
                         const Real num_diff_coeff,
                         const Array4<const Real>& prim_data,
                         const Array4<const Real>& cell_data,
                         const Array4<      Real>& rhs,
                         const Array4<const Real>& mfx_arr,
                         const Array4<const Real>& mfy_arr)
{
    BL_PROFILE_VAR("NumericalDiffusion_Scal()",NumericalDiffusion_Scal);

    // Capture diffusion coeff
    Real coeff6 = num_diff_coeff / (2.0 * dt);

    // Compute 5th order derivative and augment RHS
    ParallelFor(bx, num_comp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int m) noexcept
    {
        int n   = start_comp + m;     // conserved index
        int nm1 = (n==0) ? 0 : n - 1; // prim index
        Real rho_x_lo = (n==0) ? 1.0 : 0.5 * ( cell_data(i-1,j,k,Rho_comp) + cell_data(i  ,j,k,Rho_comp) );
        Real xflux_lo = rho_x_lo * calc_fifth_order_deriv(prim_data(i+2,j,k,nm1), prim_data(i+1,j,k,nm1),
                                                          prim_data(i  ,j,k,nm1), prim_data(i-1,j,k,nm1),
                                                          prim_data(i-2,j,k,nm1), prim_data(i-3,j,k,nm1));
        if ( (xflux_lo * (prim_data(i,j,k,nm1) - prim_data(i-1,j,k,nm1)) ) < 0.) xflux_lo = 0.;


        Real rho_x_hi = (n==0) ? 1.0 : 0.5 * ( cell_data(i+1,j,k,Rho_comp) + cell_data(i  ,j,k,Rho_comp) );
        Real xflux_hi = rho_x_hi * calc_fifth_order_deriv(prim_data(i+3,j,k,nm1), prim_data(i+2,j,k,nm1),
                                                          prim_data(i+1,j,k,nm1), prim_data(i  ,j,k,nm1),
                                                          prim_data(i-1,j,k,nm1), prim_data(i-2,j,k,nm1));
        if ( (xflux_hi * (prim_data(i+1,j,k,nm1) - prim_data(i,j,k,nm1)) ) < 0.) xflux_hi = 0.;


        Real rho_y_lo = (n==0) ? 1.0 : 0.5 * ( cell_data(i,j-1,k,Rho_comp) + cell_data(i,j  ,k,Rho_comp) );
        Real yflux_lo = rho_y_lo * calc_fifth_order_deriv(prim_data(i,j+2,k,nm1), prim_data(i,j+1,k,nm1),
                                                          prim_data(i,j  ,k,nm1), prim_data(i,j-1,k,nm1),
                                                          prim_data(i,j-2,k,nm1), prim_data(i,j-3,k,nm1));
        if ( (yflux_lo * (prim_data(i,j,k,nm1) - prim_data(i,j-1,k,nm1)) ) < 0.) yflux_lo = 0.;


        Real rho_y_hi = (n==0) ? 1.0 : 0.5 * ( cell_data(i,j+1,k,Rho_comp) + cell_data(i,j  ,k,Rho_comp) );
        Real yflux_hi = rho_y_hi * calc_fifth_order_deriv(prim_data(i,j+3,k,nm1), prim_data(i,j+2,k,nm1),
                                                          prim_data(i,j+1,k,nm1), prim_data(i,j  ,k,nm1),
                                                          prim_data(i,j-1,k,nm1), prim_data(i,j-2,k,nm1));
        if ( (yflux_hi * (prim_data(i,j+1,k,nm1) - prim_data(i,j,k,nm1)) ) < 0.) yflux_hi = 0.;


        rhs(i,j,k,n) += coeff6 * ( mfx_arr(i,j,0) * (xflux_hi - xflux_lo)
                                  +mfy_arr(i,j,0) * (yflux_hi - yflux_lo) );
    });
}

/**
 * Function to compute 6th order numerical diffusion RHS.
 *
 * @param[in]  bx box to loop over
 * @param[in]  start_comp staring component index
 * @param[in]  num_comp number of total components
 * @param[in]  num_diff_coeff
 * @param[in]  prim_data primitive variables
 * @param[in]  cell_data cell center variables
 * @param[out] rhs store the right hand side
 * @param[in]  mf map factor
 */
void
NumericalDiffusion_Xmom (const Box& bx,
                         const Real dt,
                         const Real num_diff_coeff,
                         const Array4<const Real>& prim_data,
                         const Array4<const Real>& cell_data,
                         const Array4<      Real>& rhs,
                         const Array4<const Real>& mfx_arr,
                         const Array4<const Real>& mfy_arr)
{
    BL_PROFILE_VAR("NumericalDiffusion_Xmom()",NumericalDiffusion_Xmom);

    // Capture diffusion coeff
    Real coeff6 = num_diff_coeff / (2.0 * dt);

    // Compute 5th order derivative and augment RHS
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real rho_x_lo = cell_data(i-1,j,k,Rho_comp);
        Real xflux_lo = rho_x_lo * calc_fifth_order_deriv(prim_data(i+2,j,k), prim_data(i+1,j,k),
                                                          prim_data(i  ,j,k), prim_data(i-1,j,k),
                                                          prim_data(i-2,j,k), prim_data(i-3,j,k));
        if ( (xflux_lo * (prim_data(i,j,k) - prim_data(i-1,j,k)) ) < 0.) xflux_lo = 0.;


        Real rho_x_hi = cell_data(i  ,j,k,Rho_comp);
        Real xflux_hi = rho_x_hi * calc_fifth_order_deriv(prim_data(i+3,j,k), prim_data(i+2,j,k),
                                                          prim_data(i+1,j,k), prim_data(i  ,j,k),
                                                          prim_data(i-1,j,k), prim_data(i-2,j,k));
        if ( (xflux_hi * (prim_data(i+1,j,k) - prim_data(i,j,k)) ) < 0.) xflux_hi = 0.;


        Real rho_y_lo = 0.25 * ( cell_data(i  ,j  ,k,Rho_comp) + cell_data(i-1,j  ,k,Rho_comp)
                               + cell_data(i  ,j-1,k,Rho_comp) + cell_data(i-1,j-1,k,Rho_comp) );
        Real yflux_lo = rho_y_lo * calc_fifth_order_deriv(prim_data(i,j+2,k), prim_data(i,j+1,k),
                                                          prim_data(i,j  ,k), prim_data(i,j-1,k),
                                                          prim_data(i,j-2,k), prim_data(i,j-3,k));
        if ( (yflux_lo * (prim_data(i,j,k) - prim_data(i,j-1,k)) ) < 0.) yflux_lo = 0.;


        Real rho_y_hi = 0.25 * ( cell_data(i  ,j  ,k,Rho_comp) + cell_data(i-1,j  ,k,Rho_comp)
                               + cell_data(i  ,j+1,k,Rho_comp) + cell_data(i-1,j+1,k,Rho_comp) );
        Real yflux_hi = rho_y_hi * calc_fifth_order_deriv(prim_data(i,j+3,k), prim_data(i,j+2,k),
                                                          prim_data(i,j+1,k), prim_data(i,j  ,k),
                                                          prim_data(i,j-1,k), prim_data(i,j-2,k));
        if ( (yflux_hi * (prim_data(i,j+1,k) - prim_data(i,j,k)) ) < 0.) yflux_hi = 0.;


        rhs(i,j,k,0) += coeff6 * ( mfx_arr(i,j,0) * (xflux_hi - xflux_lo)
                                  +mfy_arr(i,j,0) * (yflux_hi - yflux_lo) );
    });
}


/**
 * Function to compute 6th order numerical diffusion RHS.
 *
 * @param[in]  bx box to loop over
 * @param[in]  start_comp staring component index
 * @param[in]  num_comp number of total components
 * @param[in]  num_diff_coeff
 * @param[in]  prim_data primitive variables
 * @param[in]  cell_data cell center variables
 * @param[out] rhs store the right hand side
 * @param[in]  mf map factor
 */
void
NumericalDiffusion_Ymom (const Box& bx,
                         const Real dt,
                         const Real num_diff_coeff,
                         const Array4<const Real>& prim_data,
                         const Array4<const Real>& cell_data,
                         const Array4<      Real>& rhs,
                         const Array4<const Real>& mfx_arr,
                         const Array4<const Real>& mfy_arr)
{
    BL_PROFILE_VAR("NumericalDiffusion_Ymom()",NumericalDiffusion_Ymom);

    // Capture diffusion coeff
    Real coeff6 = num_diff_coeff / (2.0 * dt);

    // Compute 5th order derivative and augment RHS
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real rho_x_lo = 0.25 * ( cell_data(i  ,j  ,k,Rho_comp) + cell_data(i  ,j-1,k,Rho_comp)
                               + cell_data(i-1,j  ,k,Rho_comp) + cell_data(i-1,j-1,k,Rho_comp) );
        Real xflux_lo = rho_x_lo * calc_fifth_order_deriv(prim_data(i+2,j,k), prim_data(i+1,j,k),
                                                          prim_data(i  ,j,k), prim_data(i-1,j,k),
                                                          prim_data(i-2,j,k), prim_data(i-3,j,k));
        if ( (xflux_lo * (prim_data(i,j,k) - prim_data(i-1,j,k)) ) < 0.) xflux_lo = 0.;


        Real rho_x_hi = 0.25 * ( cell_data(i  ,j  ,k,Rho_comp) + cell_data(i  ,j-1,k,Rho_comp)
                               + cell_data(i+1,j  ,k,Rho_comp) + cell_data(i+1,j-1,k,Rho_comp) );
        Real xflux_hi = rho_x_hi * calc_fifth_order_deriv(prim_data(i+3,j,k), prim_data(i+2,j,k),
                                                          prim_data(i+1,j,k), prim_data(i  ,j,k),
                                                          prim_data(i-1,j,k), prim_data(i-2,j,k));
        if ( (xflux_hi * (prim_data(i+1,j,k) - prim_data(i,j,k)) ) < 0.) xflux_hi = 0.;


        Real rho_y_lo = cell_data(i,j-1,k,Rho_comp);
        Real yflux_lo = rho_y_lo * calc_fifth_order_deriv(prim_data(i,j+2,k), prim_data(i,j+1,k),
                                                          prim_data(i,j  ,k), prim_data(i,j-1,k),
                                                          prim_data(i,j-2,k), prim_data(i,j-3,k));
        if ( (yflux_lo * (prim_data(i,j,k) - prim_data(i,j-1,k)) ) < 0.) yflux_lo = 0.;


        Real rho_y_hi = cell_data(i,j  ,k,Rho_comp);
        Real yflux_hi = rho_y_hi * calc_fifth_order_deriv(prim_data(i,j+3,k), prim_data(i,j+2,k),
                                                          prim_data(i,j+1,k), prim_data(i,j  ,k),
                                                          prim_data(i,j-1,k), prim_data(i,j-2,k));
        if ( (yflux_hi * (prim_data(i,j ,k) - prim_data(i,j,k)) ) < 0.) yflux_hi = 0.;


        rhs(i,j,k,0) += coeff6 * ( mfx_arr(i,j,0) * (xflux_hi - xflux_lo)
                                  +mfy_arr(i,j,0) * (yflux_hi - yflux_lo) );
    });
}
