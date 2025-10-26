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
            const Array4<Real      >&   rhs_arr = rhs.array(mfi);

            const Array4<Real const>&      mf_mx = mapfac[lev][MapFacType::m_x]->const_array(mfi);
            const Array4<Real const>&      mf_my = mapfac[lev][MapFacType::m_y]->const_array(mfi);
            const Array4<Real const>&      mf_vx = mapfac[lev][MapFacType::v_x]->const_array(mfi);
            const Array4<Real const>&      mf_uy = mapfac[lev][MapFacType::u_y]->const_array(mfi);

            if (SolverChoice::mesh_type == MeshType::StretchedDz)
            {
                Real* stretched_dz_d_ptr = stretched_dz_d[lev].data();
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real inv_dz = 1.0/stretched_dz_d_ptr[k];
                    Real mfsq   = mf_mx(i,j,0) * mf_my(i,j,0);
                    rhs_arr(i,j,k) = (  (rho0u_arr(i+1,j  ,k  )/mf_uy(i+1,j,0) - rho0u_arr(i,j,k)/mf_uy(i,j,0)) * dxInv[0]
                                       +(rho0v_arr(i  ,j+1,k  )/mf_vx(i,j+1,0) - rho0v_arr(i,j,k)/mf_vx(i,j,0)) * dxInv[1]
                                       +(rho0w_arr(i  ,j  ,k+1)/mfsq           - rho0w_arr(i,j,k)/mfsq        ) * inv_dz  ) * mfsq;
                });
            }
            else
            {
                //
                // Note we compute the divergence using "rho0w" == Omega
                //
                const Array4<Real const>& ax_arr = ax[lev]->const_array(mfi);
                const Array4<Real const>& ay_arr = ay[lev]->const_array(mfi);
                const Array4<Real const>& dJ_arr = detJ_cc[lev]->const_array(mfi);

                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);
                    rhs_arr(i,j,k) = ( ( ax_arr(i+1,j,k)*rho0u_arr(i+1,j,k)/mf_uy(i+1,j,0)
                                        -ax_arr(i  ,j,k)*rho0u_arr(i  ,j,k)/mf_uy(i  ,j,0)  ) * dxInv[0]
                                     + ( ay_arr(i,j+1,k)*rho0v_arr(i,j+1,k)/mf_vx(i,j+1,0)
                                        -ay_arr(i,j  ,k)*rho0v_arr(i,j  ,k)/mf_vx(i,j  ,0)  ) * dxInv[1]
                                      +(                 rho0w_arr(i,j,k+1)/mfsq
                                        -                rho0w_arr(i,j,k  )/mfsq            ) * dxInv[2] ) * mfsq / dJ_arr(i,j,k);
                });
            }
        } // mfi
    }
}
