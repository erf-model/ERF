#include <AMReX.H>
#include <DiffusionSrcForMom_N.H>
#include <IndexDefines.H>

using namespace amrex;

void
DiffusionSrcForMom_N (const Box& bxx, const Box& bxy , const Box& bxz,
                      const Array4<Real>& rho_u_rhs  ,
                      const Array4<Real>& rho_v_rhs  ,
                      const Array4<Real>& rho_w_rhs  ,
                      const Array4<const Real>& tau11, const Array4<const Real>& tau22,
                      const Array4<const Real>& tau33, const Array4<const Real>& tau12,
                      const Array4<const Real>& tau13, const Array4<const Real>& tau23,
                      const Array4<const Real>& cons , const SolverChoice& solverChoice,
                      const GpuArray<Real, AMREX_SPACEDIM>& dxInv,
                      const Array4<const Real>& mf_m,
                      const Array4<const Real>& mf_u,
                      const Array4<const Real>& mf_v)
{
    BL_PROFILE_VAR("DiffusionSrcForMom_N()",DiffusionSrcForMom_N);

    auto dxinv = dxInv[0], dyinv = dxInv[1], dzinv = dxInv[2];

    if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha)
    {
        amrex::ParallelFor(bxx, bxy, bxz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real mf   = mf_m(i,j,0);

            Real diffContrib  = ( (tau11(i  , j  , k  ) - tau11(i-1, j  ,k  )) * dxinv * mf   // Contribution to x-mom eqn from diffusive flux in x-dir
                                + (tau12(i  , j+1, k  ) - tau12(i  , j  ,k  )) * dyinv * mf   // Contribution to x-mom eqn from diffusive flux in y-dir
                                + (tau13(i  , j  , k+1) - tau13(i  , j  ,k  )) * dzinv );     // Contribution to x-mom eqn from diffusive flux in z-dir;
            diffContrib      *= 0.5 * (cons(i,j,k,Rho_comp) + cons(i-1,j,k,Rho_comp))  / solverChoice.rho0_trans;
            rho_u_rhs(i,j,k) += diffContrib;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real mf   = mf_m(i,j,0);

            Real diffContrib  = ( (tau12(i+1, j  , k  ) - tau12(i  , j  , k  )) * dxinv * mf   // Contribution to y-mom eqn from diffusive flux in x-dir
                                + (tau22(i  , j  , k  ) - tau22(i  , j-1, k  )) * dyinv * mf   // Contribution to y-mom eqn from diffusive flux in y-dir
                                + (tau23(i  , j  , k+1) - tau23(i  , j  , k  )) * dzinv );     // Contribution to y-mom eqn from diffusive flux in z-dir;
            diffContrib      *= 0.5 * (cons(i,j,k,Rho_comp) + cons(i,j-1,k,Rho_comp))  / solverChoice.rho0_trans;
            rho_v_rhs(i,j,k) += diffContrib;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real mf   = mf_m(i,j,0);

            Real diffContrib  = ( (tau13(i+1, j  , k  ) - tau13(i  , j  , k  )) * dxinv * mf   // Contribution to z-mom eqn from diffusive flux in x-dir
                                + (tau23(i  , j+1, k  ) - tau23(i  , j  , k  )) * dyinv * mf   // Contribution to z-mom eqn from diffusive flux in y-dir
                                + (tau33(i  , j  , k  ) - tau33(i  , j  , k-1)) * dzinv );     // Contribution to z-mom eqn from diffusive flux in z-dir;
            diffContrib      *= 0.5 * (cons(i,j,k,Rho_comp) + cons(i,j,k-1,Rho_comp))  / solverChoice.rho0_trans;
            rho_w_rhs(i,j,k) += diffContrib;
        });

    } else {
        amrex::ParallelFor(bxx, bxy, bxz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real mf   = mf_m(i,j,0);
            
            rho_u_rhs(i,j,k) += ( (tau11(i  , j  , k  ) - tau11(i-1, j  ,k  )) * dxinv * mf   // Contribution to x-mom eqn from diffusive flux in x-dir
                                + (tau12(i  , j+1, k  ) - tau12(i  , j  ,k  )) * dyinv * mf   // Contribution to x-mom eqn from diffusive flux in y-dir
                                + (tau13(i  , j  , k+1) - tau13(i  , j  ,k  )) * dzinv );     // Contribution to x-mom eqn from diffusive flux in z-dir;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real mf   = mf_m(i,j,0);

            rho_v_rhs(i,j,k) += ( (tau12(i+1, j  , k  ) - tau12(i  , j  , k  )) * dxinv * mf   // Contribution to y-mom eqn from diffusive flux in x-dir
                                + (tau22(i  , j  , k  ) - tau22(i  , j-1, k  )) * dyinv * mf   // Contribution to y-mom eqn from diffusive flux in y-dir
                                + (tau23(i  , j  , k+1) - tau23(i  , j  , k  )) * dzinv );     // Contribution to y-mom eqn from diffusive flux in z-dir;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real mf   = mf_m(i,j,0);

            rho_w_rhs(i,j,k) += ( (tau13(i+1, j  , k  ) - tau13(i  , j  , k  )) * dxinv * mf   // Contribution to z-mom eqn from diffusive flux in x-dir
                                + (tau23(i  , j+1, k  ) - tau23(i  , j  , k  )) * dyinv * mf   // Contribution to z-mom eqn from diffusive flux in y-dir
                                + (tau33(i  , j  , k  ) - tau33(i  , j  , k-1)) * dzinv );     // Contribution to z-mom eqn from diffusive flux in z-dir;
        });
    }
}

