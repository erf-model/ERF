#include <AMReX.H>
#include <ERF_Diffusion.H>
#include <ERF_IndexDefines.H>

using namespace amrex;

/**
 * Function for computing the momentum RHS for diffusion operator without terrain.
 *
 * @param[in]  bxx nodal x box for x-mom
 * @param[in]  bxy nodal y box for y-mom
 * @param[in]  bxz nodal z box for z-mom
 * @param[out] rho_u_rhs RHS for x-mom
 * @param[out] rho_v_rhs RHS for y-mom
 * @param[out] rho_w_rhs RHS for z-mom
 * @param[in]  tau11 11 stress
 * @param[in]  tau22 22 stress
 * @param[in]  tau33 33 stress
 * @param[in]  tau12 12 stress
 * @param[in]  tau13 13 stress
 * @param[in]  tau23 23 stress
 * @param[in]  dxInv inverse cell size array
 * @param[in]  mf_m map factor at cell center
 */
void
DiffusionSrcForMom_N (const Box& bxx, const Box& bxy , const Box& bxz,
                      const Array4<Real>& rho_u_rhs  ,
                      const Array4<Real>& rho_v_rhs  ,
                      const Array4<Real>& rho_w_rhs  ,
                      const Array4<const Real>& tau11,
                      const Array4<const Real>& tau22,
                      const Array4<const Real>& tau33,
                      const Array4<const Real>& tau12,
                      const Array4<const Real>& tau13,
                      const Array4<const Real>& tau23,
                      const GpuArray<Real, AMREX_SPACEDIM>& dxInv,
                      const Array4<const Real>& mf_mx,
                      const Array4<const Real>& mf_ux,
                      const Array4<const Real>& mf_vx,
                      const Array4<const Real>& mf_my,
                      const Array4<const Real>& mf_uy,
                      const Array4<const Real>& mf_vy)
{
    BL_PROFILE_VAR("DiffusionSrcForMom_N()",DiffusionSrcForMom_N);

    auto dxinv = dxInv[0], dyinv = dxInv[1], dzinv = dxInv[2];

    ParallelFor(bxx, bxy, bxz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        // Inv Jacobian
        Real mfsq = mf_ux(i,j,0) * mf_uy(i,j,0);

        // Area corrections
        Real Imfy_hi = 1. / mf_my(i  ,j,0);
        Real Imfy_lo = 1. / mf_my(i-1,j,0);
        Real Imfx_hi = 1. / (0.5 * (mf_vx(i,j+1,0) + mf_vx(i-1,j+1,0)));
        Real Imfx_lo = 1. / (0.5 * (mf_vx(i,j  ,0) + mf_vx(i-1,j  ,0)));
        rho_u_rhs(i,j,k) -= ( (tau11(i  , j  , k  )*Imfy_hi - tau11(i-1, j  ,k  )*Imfy_lo) * dxinv * mfsq   // Contribution to x-mom eqn from diffusive flux in x-dir
                            + (tau12(i  , j+1, k  )*Imfx_hi - tau12(i  , j  ,k  )*Imfx_lo) * dyinv * mfsq   // Contribution to x-mom eqn from diffusive flux in y-dir
                            + (tau13(i  , j  , k+1)         - tau13(i  , j  ,k  )        ) * dzinv );       // Contribution to x-mom eqn from diffusive flux in z-dir;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        // Inv Jacobian
        Real mfsq = mf_vx(i,j,0) * mf_vy(i,j,0);

        // Area corrections
        Real Imfy_hi = 1. / (0.5 * (mf_uy(i+1,j,0) + mf_uy(i+1,j-1,0)));
        Real Imfy_lo = 1. / (0.5 * (mf_uy(i  ,j,0) + mf_uy(i  ,j-1,0)));
        Real Imfx_hi = 1. / mf_mx(i  ,j,0);
        Real Imfx_lo = 1. / mf_mx(i-1,j,0);
        rho_v_rhs(i,j,k) -= ( (tau12(i+1, j  , k  )*Imfy_hi - tau12(i  , j  , k  )*Imfy_lo) * dxinv * mfsq  // Contribution to y-mom eqn from diffusive flux in x-dir
                            + (tau22(i  , j  , k  )*Imfx_hi - tau22(i  , j-1, k  )*Imfx_lo) * dyinv * mfsq  // Contribution to y-mom eqn from diffusive flux in y-dir
                            + (tau23(i  , j  , k+1)         - tau23(i  , j  , k  )        ) * dzinv );      // Contribution to y-mom eqn from diffusive flux in z-dir;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        // Inv Jacobian
        Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);

        // Area corrections
        Real Imfy_hi = 1. / mf_uy(i+1,j  ,0);
        Real Imfy_lo = 1. / mf_uy(i  ,j  ,0);
        Real Imfx_hi = 1. / mf_vx(i  ,j+1,0);
        Real Imfx_lo = 1. / mf_vx(i  ,j  ,0);
        rho_w_rhs(i,j,k) -= ( (tau13(i+1, j  , k  )*Imfy_hi - tau13(i  , j  , k  )*Imfy_lo) * dxinv * mfsq  // Contribution to z-mom eqn from diffusive flux in x-dir
                            + (tau23(i  , j+1, k  )*Imfx_hi - tau23(i  , j  , k  )*Imfx_lo) * dyinv * mfsq  // Contribution to z-mom eqn from diffusive flux in y-dir
                            + (tau33(i  , j  , k  )         - tau33(i  , j  , k-1)        ) * dzinv );      // Contribution to z-mom eqn from diffusive flux in z-dir;
    });
}
