#include <AMReX.H>
#include <Diffusion.H>
#include <IndexDefines.H>

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
 * @param[in]  cons conserved cell center quantities
 * @param[in]  diffChoice container with diffusion parameters
 * @param[in]  dxInv inverse cell size array
 * @param[in]  mf_m map factor at cell center
 */
void
DiffusionSrcForMom_N (const Box& bxx, const Box& bxy , const Box& bxz,
                      const Array4<Real>& rho_u_rhs  ,
                      const Array4<Real>& rho_v_rhs  ,
                      const Array4<Real>& rho_w_rhs  ,
                      const Array4<const Real>& tau11, const Array4<const Real>& tau22,
                      const Array4<const Real>& tau33, const Array4<const Real>& tau12,
                      const Array4<const Real>& tau13, const Array4<const Real>& tau23,
                      const Array4<const Real>& cons , const DiffChoice& diffChoice,
                      const GpuArray<Real, AMREX_SPACEDIM>& dxInv,
                      const Array4<const Real>& mf_m,
                      const Array4<const Real>& /*mf_u*/,
                      const Array4<const Real>& /*mf_v*/)
{
    BL_PROFILE_VAR("DiffusionSrcForMom_N()",DiffusionSrcForMom_N);

    auto dxinv = dxInv[0], dyinv = dxInv[1], dzinv = dxInv[2];

    if (diffChoice.molec_diff_type == MolecDiffType::ConstantAlpha)
    {
        auto rho0_trans = diffChoice.rho0_trans;
        ParallelFor(bxx, bxy, bxz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real mf   = mf_m(i,j,0);

            Real diffContrib  = ( (tau11(i  , j  , k  ) - tau11(i-1, j  ,k  )) * dxinv * mf   // Contribution to x-mom eqn from diffusive flux in x-dir
                                + (tau12(i  , j+1, k  ) - tau12(i  , j  ,k  )) * dyinv * mf   // Contribution to x-mom eqn from diffusive flux in y-dir
                                + (tau13(i  , j  , k+1) - tau13(i  , j  ,k  )) * dzinv );     // Contribution to x-mom eqn from diffusive flux in z-dir;
            diffContrib      *= 0.5 * (cons(i,j,k,Rho_comp) + cons(i-1,j,k,Rho_comp))  / rho0_trans;
            rho_u_rhs(i,j,k) -= diffContrib;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real mf   = mf_m(i,j,0);

            Real diffContrib  = ( (tau12(i+1, j  , k  ) - tau12(i  , j  , k  )) * dxinv * mf   // Contribution to y-mom eqn from diffusive flux in x-dir
                                + (tau22(i  , j  , k  ) - tau22(i  , j-1, k  )) * dyinv * mf   // Contribution to y-mom eqn from diffusive flux in y-dir
                                + (tau23(i  , j  , k+1) - tau23(i  , j  , k  )) * dzinv );     // Contribution to y-mom eqn from diffusive flux in z-dir;
            diffContrib      *= 0.5 * (cons(i,j,k,Rho_comp) + cons(i,j-1,k,Rho_comp))  / rho0_trans;
            rho_v_rhs(i,j,k) -= diffContrib;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real mf   = mf_m(i,j,0);

            Real diffContrib  = ( (tau13(i+1, j  , k  ) - tau13(i  , j  , k  )) * dxinv * mf   // Contribution to z-mom eqn from diffusive flux in x-dir
                                + (tau23(i  , j+1, k  ) - tau23(i  , j  , k  )) * dyinv * mf   // Contribution to z-mom eqn from diffusive flux in y-dir
                                + (tau33(i  , j  , k  ) - tau33(i  , j  , k-1)) * dzinv );     // Contribution to z-mom eqn from diffusive flux in z-dir;
            diffContrib      *= 0.5 * (cons(i,j,k,Rho_comp) + cons(i,j,k-1,Rho_comp))  / rho0_trans;
            rho_w_rhs(i,j,k) -= diffContrib;
        });

    } else {
        ParallelFor(bxx, bxy, bxz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real mf   = mf_m(i,j,0);

            rho_u_rhs(i,j,k) -= ( (tau11(i  , j  , k  ) - tau11(i-1, j  ,k  )) * dxinv * mf   // Contribution to x-mom eqn from diffusive flux in x-dir
                                + (tau12(i  , j+1, k  ) - tau12(i  , j  ,k  )) * dyinv * mf   // Contribution to x-mom eqn from diffusive flux in y-dir
                                + (tau13(i  , j  , k+1) - tau13(i  , j  ,k  )) * dzinv );     // Contribution to x-mom eqn from diffusive flux in z-dir;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real mf   = mf_m(i,j,0);

            rho_v_rhs(i,j,k) -= ( (tau12(i+1, j  , k  ) - tau12(i  , j  , k  )) * dxinv * mf   // Contribution to y-mom eqn from diffusive flux in x-dir
                                + (tau22(i  , j  , k  ) - tau22(i  , j-1, k  )) * dyinv * mf   // Contribution to y-mom eqn from diffusive flux in y-dir
                                + (tau23(i  , j  , k+1) - tau23(i  , j  , k  )) * dzinv );     // Contribution to y-mom eqn from diffusive flux in z-dir;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real mf   = mf_m(i,j,0);

            rho_w_rhs(i,j,k) -= ( (tau13(i+1, j  , k  ) - tau13(i  , j  , k  )) * dxinv * mf   // Contribution to z-mom eqn from diffusive flux in x-dir
                                + (tau23(i  , j+1, k  ) - tau23(i  , j  , k  )) * dyinv * mf   // Contribution to z-mom eqn from diffusive flux in y-dir
                                + (tau33(i  , j  , k  ) - tau33(i  , j  , k-1)) * dzinv );     // Contribution to z-mom eqn from diffusive flux in z-dir;
        });
    }
}
