#include <AMReX.H>
#include <DiffusionSrcForMom_N.H>
#include <DiffusionSrcForMom_T.H>
#include <IndexDefines.H>

using namespace amrex;

void
DiffusionSrcForMom_T (int level, const Box& bx, const Box& valid_bx, const Box& domain,
                      const Array4<const  int>& mlo_x    , const Array4<const int>& mlo_y,
                      const Array4<const  int>& mhi_x    , const Array4<const int>& mhi_y,
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

    const int l_use_terrain   = solverChoice.use_terrain;

    AMREX_ALWAYS_ASSERT(l_use_terrain);

    const int domhi_z = domain.bigEnd(2);

    int vlo_x = valid_bx.smallEnd(0); int vhi_x = valid_bx.bigEnd(0);
    int vlo_y = valid_bx.smallEnd(1); int vhi_y = valid_bx.bigEnd(1);
    int vlo_z = valid_bx.smallEnd(2); int vhi_z = valid_bx.bigEnd(2);

    Box bxx = surroundingNodes(bx,0);
    Box bxy = surroundingNodes(bx,1);
    Box bxz = surroundingNodes(bx,2);

    // ******************************************************************
    // This assumes that refined regions are always rectangular
    // ******************************************************************
    bool left_edge_dirichlet = ( level > 0 && mlo_x(vlo_x,vlo_y,vlo_z) );
    bool rght_edge_dirichlet = ( level > 0 && mhi_x(vhi_x,vhi_y,vhi_z) );
    bool  bot_edge_dirichlet = ( level > 0 && mlo_y(vlo_x,vlo_y,vlo_z) );
    bool  top_edge_dirichlet = ( level > 0 && mhi_y(vhi_x,vhi_y,vhi_z) );

    if (left_edge_dirichlet) bxx.growLo(0,-1);
    if (rght_edge_dirichlet) bxx.growHi(0,-1);
    if ( bot_edge_dirichlet) bxy.growLo(1,-1);
    if ( top_edge_dirichlet) bxy.growHi(1,-1);

    // We don't compute diffusion src for w at k = 0
    bxz.setSmall(2,1);

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
                                                         z_nd, detJ, domain, bc_ptr);
        } else {
            diff_update = DiffusionSrcForXMom(i, j, k, u, v, w, cell_data,
                                              dxInv, K_turb, solverChoice,
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
                                                         z_nd, detJ, domain, bc_ptr);
        } else {
            diff_update = DiffusionSrcForYMom(i, j, k, u, v, w, cell_data,
                                              dxInv, K_turb, solverChoice,
                                              domain, bc_ptr, er_arr);
        }
        rho_v_rhs(i, j, k) += diff_update;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        int k_diff = (k == domhi_z+1) ? domhi_z : k;
        amrex::Real diff_update;
        if (k < domhi_z) {
            diff_update = DiffusionSrcForZMomWithTerrain(i, j, k_diff, u, v, w, cell_data,
                                                         dxInv, K_turb, solverChoice,
                                                         z_nd, detJ, domain, bc_ptr);
        } else {
            diff_update = DiffusionSrcForZMom(i, j, k_diff, u, v, w, cell_data,
                                              dxInv, K_turb, solverChoice,
                                              domain, bc_ptr, er_arr);
        }
        rho_w_rhs(i, j, k) += diff_update;
    });
}
