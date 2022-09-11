#include <AMReX.H>
#include <DiffusionSrcForMom_N.H>
#include <DiffusionSrcForMom_T.H>
#include <IndexDefines.H>

using namespace amrex;

void
DiffusionSrcForMom_T (const Box& bxx, const Box& bxy, const Box& bxz, const Box& domain,
                      const Array4<      Real>& rho_u_rhs, const Array4<      Real>& rho_v_rhs,
                      const Array4<      Real>& rho_w_rhs,
                      const Array4<const Real>& u        , const Array4<const Real>& v,
                      const Array4<const Real>& w        , const Array4<const Real>& K_turb,
                      const Array4<const Real>& cell_data, const Array4<const Real>& er_arr,
                      const SolverChoice& solverChoice   , const BCRec* bc_ptr,
                      const Array4<const Real>& z_nd     , const Array4<const Real>& detJ,
                      const GpuArray<Real, AMREX_SPACEDIM>& dxInv)
{
    BL_PROFILE_VAR("DiffusionSrcForMom_T()",DiffusionSrcForMom_T);

    if ( (solverChoice.molec_diff_type != MolecDiffType::None) ||
         (solverChoice.les_type        !=       LESType::None) ||
         (solverChoice.pbl_type        !=       PBLType::None) )
    {

    const int l_use_terrain   = solverChoice.use_terrain;

    AMREX_ALWAYS_ASSERT(l_use_terrain);

    int domhi_z = domain.bigEnd(2);

    // *********************************************************************
    // Define diffusive updates in the RHS of {x, y, z}-momentum equations
    // *********************************************************************
    amrex::ParallelFor(bxx, bxy, bxz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        amrex::Real diff_update;
        if (k < domhi_z) {
            diff_update = DiffusionSrcForXMomWithTerrain(i, j, k, u, v, w, cell_data,
                                                         dxInv, K_turb, solverChoice,
                                                         z_nd, detJ,
                                                         domain, bc_ptr, er_arr);
        } else {
            const GpuArray<Real, AMREX_SPACEDIM>& dxInv_terr = {dxInv[0], dxInv[1], dxInv[2]/detJ(i,j,k)};
            diff_update = DiffusionSrcForXMom(i, j, k, u, v, w, cell_data,
                                              dxInv_terr, K_turb, solverChoice,
                                              domain, bc_ptr, er_arr);
        }
        rho_u_rhs(i, j, k) += diff_update;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        amrex::Real diff_update;
        if (k < domhi_z) {
            diff_update = DiffusionSrcForYMomWithTerrain(i, j, k, u, v, w, cell_data,
                                                         dxInv, K_turb, solverChoice,
                                                         z_nd, detJ,
                                                         domain, bc_ptr, er_arr);
        } else {
            const GpuArray<Real, AMREX_SPACEDIM>& dxInv_terr = {dxInv[0], dxInv[1], dxInv[2]/detJ(i,j,k)};
            diff_update = DiffusionSrcForYMom(i, j, k, u, v, w, cell_data,
                                              dxInv_terr, K_turb, solverChoice,
                                              domain, bc_ptr, er_arr);
        }
        rho_v_rhs(i, j, k) += diff_update;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        amrex::Real diff_update;
        if (k < domhi_z) {
            diff_update = DiffusionSrcForZMomWithTerrain(i, j, k, u, v, w, cell_data,
                                                         dxInv, K_turb, solverChoice,
                                                         z_nd, detJ,
                                                         domain, bc_ptr, er_arr);
        } else {
            int k_diff = domhi_z;
            const GpuArray<Real, AMREX_SPACEDIM>& dxInv_terr = {dxInv[0], dxInv[1], dxInv[2]/detJ(i,j,k_diff)};
            diff_update = DiffusionSrcForZMom(i, j, k_diff, u, v, w, cell_data,
                                              dxInv_terr, K_turb, solverChoice,
                                              domain, bc_ptr, er_arr);
        }
        rho_w_rhs(i, j, k) += diff_update;
    });
    }
}
