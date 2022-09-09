#include <SpatialStencils.H>
#include <StressTerm.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
DiffusionSrcForXMom (int i, int j, int k,
                     const amrex::Array4<const amrex::Real>& u,
                     const amrex::Array4<const amrex::Real>& v,
                     const amrex::Array4<const amrex::Real>& w,
                     const amrex::Array4<const amrex::Real>& cons,
                     const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                     const amrex::Array4<const amrex::Real>& K_turb,
                     const SolverChoice &solverChoice,
                     const amrex::Box& domain, const amrex::BCRec* bc_ptr,
                     const amrex::Array4<const amrex::Real>& er_arr)
{
    amrex::Real diffContrib = 0.0;

    if ( (solverChoice.molec_diff_type != MolecDiffType::None) ||
         (solverChoice.les_type        !=       LESType::None) ||
         (solverChoice.pbl_type        !=       PBLType::None) )
    {
        auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

        amrex::Real tau11Next, tau11Prev, tau12Next, tau12Prev, tau13Next, tau13Prev;
        tau11Next = ComputeStressTerm(i+1, j, k, u, v, w, MomentumEqn::x,
                                      DiffusionDir::x, cellSizeInv, K_turb, solverChoice, domain, bc_ptr, er_arr(i  ,j,k) );
        tau11Prev = ComputeStressTerm(i  , j, k, u, v, w, MomentumEqn::x,
                                      DiffusionDir::x, cellSizeInv, K_turb, solverChoice, domain, bc_ptr, er_arr(i-1,j,k) );
        tau12Next = ComputeStressTerm(i, j+1, k, u, v, w, MomentumEqn::x,
                                      DiffusionDir::y, cellSizeInv, K_turb, solverChoice, domain, bc_ptr, 0.);
        tau12Prev = ComputeStressTerm(i, j  , k, u, v, w, MomentumEqn::x,
                                      DiffusionDir::y, cellSizeInv, K_turb, solverChoice, domain, bc_ptr, 0.);
        tau13Next = ComputeStressTerm(i, j, k+1, u, v, w, MomentumEqn::x,
                                      DiffusionDir::z, cellSizeInv, K_turb, solverChoice, domain, bc_ptr, 0.);
        tau13Prev = ComputeStressTerm(i, j, k  , u, v, w, MomentumEqn::x,
                                      DiffusionDir::z, cellSizeInv, K_turb, solverChoice, domain, bc_ptr, 0.);

        diffContrib = (tau11Next - tau11Prev) * dxInv  // Contribution to x-mom eqn from diffusive flux in x-dir
                    + (tau12Next - tau12Prev) * dyInv  // Contribution to x-mom eqn from diffusive flux in y-dir
                    + (tau13Next - tau13Prev) * dzInv; // Contribution to x-mom eqn from diffusive flux in z-dir
        if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha)
        {
            diffContrib *= 0.5 * (cons(i,j,k,Rho_comp) + cons(i-1,j,k,Rho_comp))  / solverChoice.rho0_trans;
        }
    }
    return diffContrib;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
DiffusionSrcForYMom (int i, int j, int k,
                     const amrex::Array4<const amrex::Real>& u,
                     const amrex::Array4<const amrex::Real>& v,
                     const amrex::Array4<const amrex::Real>& w,
                     const amrex::Array4<const amrex::Real>& cons,
                     const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                     const amrex::Array4<const amrex::Real>& K_turb,
                     const SolverChoice &solverChoice,
                     const amrex::Box& domain, const amrex::BCRec* bc_ptr,
                     const amrex::Array4<const amrex::Real>& er_arr)
{
    amrex::Real diffContrib = 0.0;

    if ( (solverChoice.molec_diff_type != MolecDiffType::None) ||
         (solverChoice.les_type        !=       LESType::None) ||
         (solverChoice.pbl_type        !=       PBLType::None) )
    {
        auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

        amrex::Real tau21Next, tau21Prev, tau22Next, tau22Prev, tau23Next, tau23Prev;
        tau21Next = ComputeStressTerm(i+1, j, k, u, v, w, MomentumEqn::y,
                                      DiffusionDir::x, cellSizeInv, K_turb, solverChoice, domain, bc_ptr, 0.);
        tau21Prev = ComputeStressTerm(i  , j, k, u, v, w, MomentumEqn::y,
                                      DiffusionDir::x, cellSizeInv, K_turb, solverChoice, domain, bc_ptr, 0.);
        tau22Next = ComputeStressTerm(i, j+1, k, u, v, w, MomentumEqn::y,
                                      DiffusionDir::y, cellSizeInv, K_turb, solverChoice, domain, bc_ptr, er_arr(i,j  ,k) );
        tau22Prev = ComputeStressTerm(i, j  , k, u, v, w, MomentumEqn::y,
                                      DiffusionDir::y, cellSizeInv, K_turb, solverChoice, domain, bc_ptr, er_arr(i,j-1,k) );
        tau23Next = ComputeStressTerm(i, j, k+1, u, v, w, MomentumEqn::y,
                                      DiffusionDir::z, cellSizeInv, K_turb, solverChoice, domain, bc_ptr, 0.);
        tau23Prev = ComputeStressTerm(i, j, k  , u, v, w, MomentumEqn::y,
                                      DiffusionDir::z, cellSizeInv, K_turb, solverChoice, domain, bc_ptr, 0.);

        diffContrib = (tau21Next - tau21Prev) * dxInv  // Contribution to y-mom eqn from diffusive flux in x-dir
                    + (tau22Next - tau22Prev) * dyInv  // Contribution to y-mom eqn from diffusive flux in y-dir
                    + (tau23Next - tau23Prev) * dzInv; // Contribution to y-mom eqn from diffusive flux in z-dir
        if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha)
        {
            diffContrib *= 0.5 * (cons(i,j,k,Rho_comp) + cons(i,j-1,k,Rho_comp))  / solverChoice.rho0_trans;
        }
    }
    return diffContrib;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
DiffusionSrcForZMom (int i, int j, int k,
                     const amrex::Array4<const amrex::Real>& u,
                     const amrex::Array4<const amrex::Real>& v,
                     const amrex::Array4<const amrex::Real>& w,
                     const amrex::Array4<const amrex::Real>& cons,
                     const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                     const amrex::Array4<const amrex::Real>& K_turb,
                     const SolverChoice &solverChoice,
                     const amrex::Box& domain, const amrex::BCRec* bc_ptr,
                     const amrex::Array4<const amrex::Real>& er_arr)
{
    amrex::Real diffContrib = 0.0;

    if ( (solverChoice.molec_diff_type != MolecDiffType::None) ||
         (solverChoice.les_type        !=       LESType::None) ||
         (solverChoice.pbl_type        !=       PBLType::None) )
    {
        auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

        amrex::Real tau31Next, tau31Prev, tau32Next, tau32Prev, tau33Next, tau33Prev;

        tau31Next = ComputeStressTerm(i+1, j, k, u, v, w, MomentumEqn::z,
                                      DiffusionDir::x, cellSizeInv, K_turb, solverChoice, domain, bc_ptr, 0.);
        tau31Prev = ComputeStressTerm(i  , j, k, u, v, w, MomentumEqn::z,
                                      DiffusionDir::x, cellSizeInv, K_turb, solverChoice, domain, bc_ptr, 0.);
        tau32Next = ComputeStressTerm(i, j+1, k, u, v, w, MomentumEqn::z,
                                      DiffusionDir::y, cellSizeInv, K_turb, solverChoice, domain, bc_ptr, 0.);
        tau32Prev = ComputeStressTerm(i, j  , k, u, v, w, MomentumEqn::z,
                                      DiffusionDir::y, cellSizeInv, K_turb, solverChoice, domain, bc_ptr, 0.);
        tau33Next = ComputeStressTerm(i, j, k+1, u, v, w, MomentumEqn::z,
                                      DiffusionDir::z, cellSizeInv, K_turb, solverChoice, domain, bc_ptr, er_arr(i,j,k  ) );
        tau33Prev = ComputeStressTerm(i, j, k  , u, v, w, MomentumEqn::z,
                                      DiffusionDir::z, cellSizeInv, K_turb, solverChoice, domain, bc_ptr, er_arr(i,j,k-1) );

        diffContrib = (tau31Next - tau31Prev) * dxInv  // Contribution to z-mom eqn from diffusive flux in x-dir
                    + (tau32Next - tau32Prev) * dyInv  // Contribution to z-mom eqn from diffusive flux in y-dir
                    + (tau33Next - tau33Prev) * dzInv; // Contribution to z-mom eqn from diffusive flux in z-dir
        if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha)
        {
            diffContrib *= 0.5 * (cons(i,j,k,Rho_comp) + cons(i,j,k-1,Rho_comp))  / solverChoice.rho0_trans;
        }
    }
    return diffContrib;
}
