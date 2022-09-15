#include <AMReX.H>
#include <DiffusionSrcForMom_N.H>
#include <IndexDefines.H>

using namespace amrex;

#ifdef ERF_DIFF_OPTIM
void
DiffusionSrcForMom_N (const Box& bxx, const Box& bxy , const Box& bxz,
                      const Array4<Real>& rho_u_rhs  ,
                      const Array4<Real>& rho_v_rhs  ,
                      const Array4<Real>& rho_w_rhs  ,
                      const Array4<const Real>& tau11, const Array4<const Real>& tau22,
                      const Array4<const Real>& tau33, const Array4<const Real>& tau12,
                      const Array4<const Real>& tau13, const Array4<const Real>& tau23,
                      const Array4<const Real>& cons , const SolverChoice& solverChoice,
                      const GpuArray<Real, AMREX_SPACEDIM>& dxInv)
{
    BL_PROFILE_VAR("DiffusionSrcForMom_N()",DiffusionSrcForMom_N);

    auto dxinv = dxInv[0], dyinv = dxInv[1], dzinv = dxInv[2];

    if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha)
    {
        amrex::ParallelFor(bxx, bxy, bxz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real diffContrib  = ( (tau11(i  , j  , k  ) - tau11(i-1, j  ,k  )) * dxinv    // Contribution to x-mom eqn from diffusive flux in x-dir
                                + (tau12(i  , j+1, k  ) - tau12(i  , j  ,k  )) * dyinv    // Contribution to x-mom eqn from diffusive flux in y-dir
                                + (tau13(i  , j  , k+1) - tau13(i  , j  ,k  )) * dzinv ); // Contribution to x-mom eqn from diffusive flux in z-dir;
            diffContrib      *= 0.5 * (cons(i,j,k,Rho_comp) + cons(i-1,j,k,Rho_comp))  / solverChoice.rho0_trans;
            rho_u_rhs(i,j,k) += diffContrib;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real diffContrib  = ( (tau12(i+1, j  , k  ) - tau12(i  , j  , k  )) * dxinv    // Contribution to y-mom eqn from diffusive flux in x-dir
                                + (tau22(i  , j  , k  ) - tau22(i  , j-1, k  )) * dyinv    // Contribution to y-mom eqn from diffusive flux in y-dir
                                + (tau23(i  , j  , k+1) - tau23(i  , j  , k  )) * dzinv ); // Contribution to y-mom eqn from diffusive flux in z-dir;
            diffContrib      *= 0.5 * (cons(i,j,k,Rho_comp) + cons(i,j-1,k,Rho_comp))  / solverChoice.rho0_trans;
            rho_v_rhs(i,j,k) += diffContrib;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real diffContrib  = ( (tau13(i+1, j  , k  ) - tau13(i  , j  , k  )) * dxinv    // Contribution to z-mom eqn from diffusive flux in x-dir
                                + (tau23(i  , j+1, k  ) - tau23(i  , j  , k  )) * dyinv    // Contribution to z-mom eqn from diffusive flux in y-dir
                                + (tau33(i  , j  , k  ) - tau33(i  , j  , k-1)) * dzinv ); // Contribution to z-mom eqn from diffusive flux in z-dir;
            diffContrib      *= 0.5 * (cons(i,j,k,Rho_comp) + cons(i,j,k-1,Rho_comp))  / solverChoice.rho0_trans;
            rho_w_rhs(i,j,k) += diffContrib;
        });

    } else {
        amrex::ParallelFor(bxx, bxy, bxz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            rho_u_rhs(i,j,k) += ( (tau11(i  , j  , k  ) - tau11(i-1, j  ,k  )) * dxinv    // Contribution to x-mom eqn from diffusive flux in x-dir
                                + (tau12(i  , j+1, k  ) - tau12(i  , j  ,k  )) * dyinv    // Contribution to x-mom eqn from diffusive flux in y-dir
                                + (tau13(i  , j  , k+1) - tau13(i  , j  ,k  )) * dzinv ); // Contribution to x-mom eqn from diffusive flux in z-dir;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            rho_v_rhs(i,j,k) += ( (tau12(i+1, j  , k  ) - tau12(i  , j  , k  )) * dxinv    // Contribution to y-mom eqn from diffusive flux in x-dir
                                + (tau22(i  , j  , k  ) - tau22(i  , j-1, k  )) * dyinv    // Contribution to y-mom eqn from diffusive flux in y-dir
                                + (tau23(i  , j  , k+1) - tau23(i  , j  , k  )) * dzinv ); // Contribution to y-mom eqn from diffusive flux in z-dir;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            rho_w_rhs(i,j,k) += ( (tau13(i+1, j  , k  ) - tau13(i  , j  , k  )) * dxinv    // Contribution to z-mom eqn from diffusive flux in x-dir
                                + (tau23(i  , j+1, k  ) - tau23(i  , j  , k  )) * dyinv    // Contribution to z-mom eqn from diffusive flux in y-dir
                                + (tau33(i  , j  , k  ) - tau33(i  , j  , k-1)) * dzinv ); // Contribution to z-mom eqn from diffusive flux in z-dir;
        });
    }
}

#else
void
DiffusionSrcForMom_N (const Box& bxx, const Box& bxy, const Box& bxz, const Box& domain,
                      const Array4<      Real>& rho_u_rhs, const Array4<      Real>& rho_v_rhs,
                      const Array4<      Real>& rho_w_rhs,
                      const Array4<const Real>& u        , const Array4<const Real>& v,
                      const Array4<const Real>& w        , const Array4<const Real>& K_turb,
                      const Array4<const Real>& cell_data, const Array4<const Real>& er_arr,
                      const SolverChoice& solverChoice   , const BCRec* bc_ptr,
                      const GpuArray<Real, AMREX_SPACEDIM>& dxInv)
{
    BL_PROFILE_VAR("DiffusionSrcForMom_N()",DiffusionSrcForMom_N);

    if ( (solverChoice.molec_diff_type != MolecDiffType::None) ||
         (solverChoice.les_type        !=       LESType::None) ||
         (solverChoice.pbl_type        !=       PBLType::None) )
    {

    const int l_use_terrain   = solverChoice.use_terrain;

    int domhi_z = domain.bigEnd(2);

    AMREX_ALWAYS_ASSERT(!l_use_terrain);

    // *********************************************************************
    // Define diffusive updates in the RHS of {x, y, z}-momentum equations
    // *********************************************************************
    amrex::ParallelFor(bxx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        rho_u_rhs(i, j, k) += DiffusionSrcForXMom(i, j, k, u, v, w, cell_data,
                                                  dxInv, K_turb, solverChoice,
                                                  domain, bc_ptr, er_arr);
    });
    amrex::ParallelFor(bxy,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        rho_v_rhs(i, j, k) += DiffusionSrcForYMom(i, j, k, u, v, w, cell_data,
                                                  dxInv, K_turb, solverChoice,
                                                  domain, bc_ptr, er_arr);
    });
    amrex::ParallelFor(bxz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        int k_diff = (k == domhi_z+1) ? domhi_z : k;
        Real temp = DiffusionSrcForZMom(i, j, k_diff, u, v, w, cell_data,
                                        dxInv, K_turb, solverChoice,
                                        domain, bc_ptr, er_arr);
        rho_w_rhs(i, j, k) += temp;
    });
    }
}
#endif
