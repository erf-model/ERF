#include <AMReX.H>
#include <Diffusion.H>
#include <IndexDefines.H>

using namespace amrex;

/**
 * Function for computing the momentum RHS for diffusion operator with terrain.
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
 * @param[in]  tau21 21 stress
 * @param[in]  tau23 23 stress
 * @param[in]  tau31 31 stress
 * @param[in]  tau32 32 stress
 * @param[in]  cons conserved cell center quantities
 * @param[in]  detJ Jacobian determinant
 * @param[in]  differChoice container with diffusion parameters
 * @param[in]  dxInv inverse cell size array
 * @param[in]  mf_m map factor at cell center
 */
void
DiffusionSrcForMom_T (const Box& bxx, const Box& bxy , const Box& bxz,
                      const Array4<Real>& rho_u_rhs  ,
                      const Array4<Real>& rho_v_rhs  ,
                      const Array4<Real>& rho_w_rhs  ,
                      const Array4<const Real>& tau11, const Array4<const Real>& tau22, const Array4<const Real>& tau33,
                      const Array4<const Real>& tau12, const Array4<const Real>& tau13,
                      const Array4<const Real>& tau21, const Array4<const Real>& tau23,
                      const Array4<const Real>& tau31, const Array4<const Real>& tau32,
                      const Array4<const Real>& cons , const Array4<const Real>& detJ ,
                      const DiffChoice& diffChoice,
                      const GpuArray<Real, AMREX_SPACEDIM>& dxInv,
                      const Array4<const Real>& mf_m,
                      const Array4<const Real>& /*mf_u*/,
                      const Array4<const Real>& /*mf_v*/)
{
    BL_PROFILE_VAR("DiffusionSrcForMom_T()",DiffusionSrcForMom_T);

    auto dxinv = dxInv[0], dyinv = dxInv[1], dzinv = dxInv[2];

    ParallelFor(bxx, bxy, bxz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real mf   = mf_m(i,j,0);

        Real diffContrib  = ( (tau11(i  , j  , k  ) - tau11(i-1, j  ,k  )) * dxinv * mf   // Contribution to x-mom eqn from diffusive flux in x-dir
                            + (tau12(i  , j+1, k  ) - tau12(i  , j  ,k  )) * dyinv * mf   // Contribution to x-mom eqn from diffusive flux in y-dir
                            + (tau13(i  , j  , k+1) - tau13(i  , j  ,k  )) * dzinv );     // Contribution to x-mom eqn from diffusive flux in z-dir;
        diffContrib      /= 0.5*(detJ(i,j,k) + detJ(i-1,j,k));
        rho_u_rhs(i,j,k) -= diffContrib;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real mf   = mf_m(i,j,0);

        Real diffContrib  = ( (tau21(i+1, j  , k  ) - tau21(i  , j  , k  )) * dxinv * mf   // Contribution to y-mom eqn from diffusive flux in x-dir
                            + (tau22(i  , j  , k  ) - tau22(i  , j-1, k  )) * dyinv * mf   // Contribution to y-mom eqn from diffusive flux in y-dir
                            + (tau23(i  , j  , k+1) - tau23(i  , j  , k  )) * dzinv );     // Contribution to y-mom eqn from diffusive flux in z-dir;
        diffContrib      /= 0.5*(detJ(i,j,k) + detJ(i,j-1,k));
        rho_v_rhs(i,j,k) -= diffContrib;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real mf   = mf_m(i,j,0);

        Real diffContrib  = ( (tau31(i+1, j  , k  ) - tau31(i  , j  , k  )) * dxinv * mf   // Contribution to z-mom eqn from diffusive flux in x-dir
                            + (tau32(i  , j+1, k  ) - tau32(i  , j  , k  )) * dyinv * mf   // Contribution to z-mom eqn from diffusive flux in y-dir
                            + (tau33(i  , j  , k  ) - tau33(i  , j  , k-1)) * dzinv );     // Contribution to z-mom eqn from diffusive flux in z-dir;
        diffContrib      /= 0.5*(detJ(i,j,k) + detJ(i,j,k-1));
        rho_w_rhs(i,j,k) -= diffContrib;
    });
}
