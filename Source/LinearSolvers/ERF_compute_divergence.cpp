#include "ERF.H"
#include "ERF_Utils.H"

using namespace amrex;

/**
 * Project the single-level velocity field to enforce incompressibility
 * Note that the level may or may not be level 0.
 */
void ERF::compute_divergence (int lev, MultiFab& rhs, Vector<MultiFab>& mom_mf, Geometry const& geom_at_lev)
{
    BL_PROFILE("ERF::compute_divergence()");

    bool l_use_terrain = (solverChoice.terrain_type != TerrainType::None);

    auto dxInv = geom[lev].InvCellSizeArray();

    Array<MultiFab const*, AMREX_SPACEDIM> rho0_u_const;
    rho0_u_const[0] = &mom_mf[IntVars::xmom];
    rho0_u_const[1] = &mom_mf[IntVars::ymom];
    rho0_u_const[2] = &mom_mf[IntVars::zmom];

    // ****************************************************************************
    // Compute divergence which will form RHS
    // Note that we replace "rho0w" with the contravariant momentum, Omega
    // ****************************************************************************
#ifdef ERF_USE_EB
    bool already_on_centroids = true;
    EB_computeDivergence(rhs, rho0_u_const, geom_at_lev, already_on_centroids);
#else
    if (l_use_terrain && SolverChoice::terrain_is_flat)
    {
        for ( MFIter mfi(rhs,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.tilebox();
            const Array4<Real const>& rho0u_arr = mom_mf[IntVars::xmom].const_array(mfi);
            const Array4<Real const>& rho0v_arr = mom_mf[IntVars::ymom].const_array(mfi);
            const Array4<Real const>& rho0w_arr = mom_mf[IntVars::zmom].const_array(mfi);
            const Array4<Real      >&  rhs_arr = rhs.array(mfi);

            Real* stretched_dz_d_ptr = stretched_dz_d[lev].data();
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Real dz = stretched_dz_d_ptr[k];
                rhs_arr(i,j,k) =  (rho0u_arr(i+1,j,k) - rho0u_arr(i,j,k)) * dxInv[0]
                                 +(rho0v_arr(i,j+1,k) - rho0v_arr(i,j,k)) * dxInv[1]
                                 +(rho0w_arr(i,j,k+1) - rho0w_arr(i,j,k)) / dz;
            });
        } // mfi
    }
    else if (l_use_terrain) // terrain is not flat
    {
        //
        // Note we compute the divergence using "rho0w" == Omega
        //
        for ( MFIter mfi(rhs,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.tilebox();
            const Array4<Real >& rho0u_arr = mom_mf[IntVars::xmom].array(mfi);
            const Array4<Real >& rho0v_arr = mom_mf[IntVars::ymom].array(mfi);
            const Array4<Real >& rho0w_arr = mom_mf[IntVars::zmom].array(mfi);
            const Array4<Real      >&  rhs_arr = rhs.array(mfi);

            const Array4<Real const>& ax_arr = ax[lev]->const_array(mfi);
            const Array4<Real const>& ay_arr = ay[lev]->const_array(mfi);
            const Array4<Real const>& az_arr = az[lev]->const_array(mfi);
            const Array4<Real const>& dJ_arr = detJ_cc[lev]->const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                rhs_arr(i,j,k) =   ((ax_arr(i+1,j,k)*rho0u_arr(i+1,j,k) - ax_arr(i,j,k)*rho0u_arr(i,j,k)) * dxInv[0]
                                   +(ay_arr(i,j+1,k)*rho0v_arr(i,j+1,k) - ay_arr(i,j,k)*rho0v_arr(i,j,k)) * dxInv[1]
                                   +(az_arr(i,j,k+1)*rho0w_arr(i,j,k+1) - az_arr(i,j,k)*rho0w_arr(i,j,k)) * dxInv[2]) / dJ_arr(i,j,k);
            });
        } // mfi

    }
    else // no terrain
    {
        computeDivergence(rhs, rho0_u_const, geom_at_lev);
    }
#endif
}
