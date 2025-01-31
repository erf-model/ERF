#include "ERF.H"
#include "ERF_Utils.H"

using namespace amrex;

/**
 * Project the single-level velocity field to enforce incompressibility
 * Note that the level may or may not be level 0.
 */
void ERF::compute_divergence (int lev, MultiFab& rhs, Array<MultiFab const*,AMREX_SPACEDIM> rho0_u_const, Geometry const& geom_at_lev)
{
    BL_PROFILE("ERF::compute_divergence()");

    auto dxInv = geom_at_lev.InvCellSizeArray();

    // ****************************************************************************
    // Compute divergence which will form RHS
    // Note that we replace "rho0w" with the contravariant momentum, Omega
    // ****************************************************************************
    if (solverChoice.terrain_type == TerrainType::EB)
    {
        bool already_on_centroids = true;
        EB_computeDivergence(rhs, rho0_u_const, geom_at_lev, already_on_centroids);
    }
    else if (SolverChoice::mesh_type == MeshType::ConstantDz)
    {
        computeDivergence(rhs, rho0_u_const, geom_at_lev);
    }
    else
    {
        for ( MFIter mfi(rhs,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.tilebox();
            const Array4<Real const>& rho0u_arr = rho0_u_const[0]->const_array(mfi);
            const Array4<Real const>& rho0v_arr = rho0_u_const[1]->const_array(mfi);
            const Array4<Real const>& rho0w_arr = rho0_u_const[2]->const_array(mfi);
            const Array4<Real      >&  rhs_arr = rhs.array(mfi);

            if (SolverChoice::mesh_type == MeshType::StretchedDz) {
                Real* stretched_dz_d_ptr = stretched_dz_d[lev].data();
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    Real dz = stretched_dz_d_ptr[k];
                    rhs_arr(i,j,k) =  (rho0u_arr(i+1,j,k) - rho0u_arr(i,j,k)) * dxInv[0]
                                     +(rho0v_arr(i,j+1,k) - rho0v_arr(i,j,k)) * dxInv[1]
                                     +(rho0w_arr(i,j,k+1) - rho0w_arr(i,j,k)) / dz;
                });
            } else {

                //
                // Note we compute the divergence using "rho0w" == Omega
                //
                const Array4<Real const>& ax_arr = ax[lev]->const_array(mfi);
                const Array4<Real const>& ay_arr = ay[lev]->const_array(mfi);
                const Array4<Real const>& dJ_arr = detJ_cc[lev]->const_array(mfi);
                //
                // az == 1 for terrain-fitted coordinates
                //
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    rhs_arr(i,j,k) =   ((ax_arr(i+1,j,k)*rho0u_arr(i+1,j,k) - ax_arr(i,j,k)*rho0u_arr(i,j,k)) * dxInv[0]
                                       +(ay_arr(i,j+1,k)*rho0v_arr(i,j+1,k) - ay_arr(i,j,k)*rho0v_arr(i,j,k)) * dxInv[1]
                                       +(                rho0w_arr(i,j,k+1) -               rho0w_arr(i,j,k)) * dxInv[2]) / dJ_arr(i,j,k);
                });
            }
        } // mfi
    }
}
