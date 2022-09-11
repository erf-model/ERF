#include <AMReX.H>
#include <DiffusionSrcForMom_N.H>
#include <IndexDefines.H>

using namespace amrex;

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
