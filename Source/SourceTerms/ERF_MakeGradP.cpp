#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

#include "ERF_SrcHeaders.H"
#include "ERF_DataStruct.H"
#include "ERF_Utils.H"

using namespace amrex;

/**
 * Function for computing the pressure gradient
 *
 * @param[in]  level     level of resolution
 * @param[in]  geom      geometry container at this level
 * @param[in]  S_data    current solution
 * @param[in]  p0        base ststa pressure
 * @param[in]  z_phys_nd z on nodes
 * @param[in]  z_phys_cc z on cell centers
 * @param[out] gradp     pressure gradient
 */

void make_gradp_pert (int level,
                      const SolverChoice& solverChoice,
                      const Geometry& geom,
                      MultiFab& S_data,
                      MultiFab& p0,
                      std::unique_ptr<MultiFab>& z_phys_nd,
                      std::unique_ptr<MultiFab>& z_phys_cc,
                      Vector<MultiFab>& gradp)
{
    const bool l_use_moisture  = (solverChoice.moisture_type != MoistureType::None);
    //
    // Note that we only recompute gradp if compressible;
    //      if anelastic then we have computed gradp in the projection
    //      and we can reuse it, no need to recompute it
    //
    if (solverChoice.anelastic[level] == 0)
    {
        MultiFab p(S_data.boxArray(), S_data.DistributionMap(), 1, 1);

        // *****************************************************************************
        // Compute pressure or perturbataional pressure
        // *****************************************************************************
        for ( MFIter mfi(S_data); mfi.isValid(); ++mfi)
        {
            Box gbx = mfi.tilebox(); gbx.grow(IntVect(1,1,1));
            if (gbx.smallEnd(2) < 0) gbx.setSmall(2,0);
            const Array4<const Real>& cell_data = S_data.array(mfi);
            const Array4<const Real>& p0_arr = p0.const_array(mfi);
            const Array4<      Real>& pptemp_arr = p.array(mfi);
            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real qv_for_p = (l_use_moisture) ? cell_data(i,j,k,RhoQ1_comp)/cell_data(i,j,k,Rho_comp) : 0.0;
                pptemp_arr(i,j,k) = getPgivenRTh(cell_data(i,j,k,RhoTheta_comp),qv_for_p) - p0_arr(i,j,k);
            });
        }

        compute_gradp(p,geom,z_phys_nd,z_phys_cc,gradp,solverChoice);
    }
}

void
compute_gradp (const MultiFab& p,
               const Geometry& geom,
               std::unique_ptr<MultiFab>& z_phys_nd,
               std::unique_ptr<MultiFab>& z_phys_cc,
               Vector<MultiFab>& gradp,
               const SolverChoice& solverChoice)
{
    const bool l_use_terrain_fitted_coords = (solverChoice.mesh_type != MeshType::ConstantDz);

    const Box domain = geom.Domain();
    const int klo = domain.smallEnd(2);
    const int khi = domain.bigEnd(2);

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();

    // *****************************************************************************
    // Take gradient of relevant quantity (p0, pres, or pert_pres = pres - p0)
    // *****************************************************************************
    for ( MFIter mfi(p); mfi.isValid(); ++mfi)
    {
        Box tbx = mfi.nodaltilebox(0);
        Box tby = mfi.nodaltilebox(1);
        Box tbz = mfi.nodaltilebox(2);

        // We don't compute gpz on the bottom or top domain boundary
        if (tbz.smallEnd(2) == klo) {
            tbz.growLo(2,-1);
        }
        if (tbz.bigEnd(2) == khi+1) {
            tbz.growHi(2,-1);
        }

        // Terrain metrics
        const Array4<const Real>& z_nd_arr = z_phys_nd->const_array(mfi);
        const Array4<const Real>& z_cc_arr = z_phys_cc->const_array(mfi);

        const Array4<const Real>& p_arr = p.const_array(mfi);

        const Array4<      Real>& gpx_arr = gradp[GpVars::gpx].array(mfi);
        const Array4<      Real>& gpy_arr = gradp[GpVars::gpy].array(mfi);
        const Array4<      Real>& gpz_arr = gradp[GpVars::gpz].array(mfi);

        ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            //Note : mx/my == 1, so no map factor needed here
            Real gpx = dxInv[0] * (p_arr(i,j,k) - p_arr(i-1,j,k));

            if (l_use_terrain_fitted_coords) {
                Real met_h_xi = (z_cc_arr(i,j,k) - z_cc_arr(i-1,j,k)) * dxInv[0];

                Real dz_phys_hi, dz_phys_lo;
                Real gpz_lo, gpz_hi;
                if (k==klo) {
                    dz_phys_hi = z_cc_arr(i  ,j,k+1) -   z_cc_arr(i  ,j,k  );
                    dz_phys_lo = z_cc_arr(i-1,j,k+1) -   z_cc_arr(i-1,j,k  );
                    gpz_hi  = (p_arr(i  ,j,k+1) - p_arr(i  ,j,k  )) / dz_phys_hi;
                    gpz_lo  = (p_arr(i-1,j,k+1) - p_arr(i-1,j,k  )) / dz_phys_lo;
                } else if (k==khi) {
                    dz_phys_hi = z_cc_arr(i  ,j,k  ) -   z_cc_arr(i  ,j,k-1);
                    dz_phys_lo = z_cc_arr(i-1,j,k  ) -   z_cc_arr(i-1,j,k-1);
                    gpz_hi  = (p_arr(i  ,j,k  ) - p_arr(i  ,j,k-1)) / dz_phys_hi;
                    gpz_lo  = (p_arr(i-1,j,k  ) - p_arr(i-1,j,k-1)) / dz_phys_lo;
                } else {
                    dz_phys_hi = z_cc_arr(i  ,j,k+1) -   z_cc_arr(i  ,j,k-1);
                    dz_phys_lo = z_cc_arr(i-1,j,k+1) -   z_cc_arr(i-1,j,k-1);
                    gpz_hi  = (p_arr(i  ,j,k+1) - p_arr(i  ,j,k-1)) / dz_phys_hi;
                    gpz_lo  = (p_arr(i-1,j,k+1) - p_arr(i-1,j,k-1)) / dz_phys_lo;
                }
                Real gpx_metric = met_h_xi * 0.5 * (gpz_hi + gpz_lo);
                gpx -= gpx_metric;
            }
            gpx_arr(i,j,k) = gpx;
        });

        ParallelFor(tby, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            //Note : mx/my == 1, so no map factor needed here
            Real gpy = dxInv[1] * (p_arr(i,j,k) - p_arr(i,j-1,k));

            if (l_use_terrain_fitted_coords) {
                Real met_h_eta = (z_cc_arr(i,j,k) - z_cc_arr(i,j-1,k)) * dxInv[1];

                Real dz_phys_hi, dz_phys_lo;
                Real gpz_lo, gpz_hi;
                if (k==klo) {
                    dz_phys_hi = z_cc_arr(i,j  ,k+1) -   z_cc_arr(i,j  ,k  );
                    dz_phys_lo = z_cc_arr(i,j-1,k+1) -   z_cc_arr(i,j-1,k  );
                    gpz_hi  = (p_arr(i,j  ,k+1) - p_arr(i,j  ,k  )) / dz_phys_hi;
                    gpz_lo  = (p_arr(i,j-1,k+1) - p_arr(i,j-1,k  )) / dz_phys_lo;
                } else if (k==khi) {
                    dz_phys_hi = z_cc_arr(i,j  ,k  ) -   z_cc_arr(i,j  ,k-1);
                    dz_phys_lo = z_cc_arr(i,j-1,k  ) -   z_cc_arr(i,j-1,k-1);
                    gpz_hi  = (p_arr(i,j  ,k  ) - p_arr(i,j  ,k-1)) / dz_phys_hi;
                    gpz_lo  = (p_arr(i,j-1,k  ) - p_arr(i,j-1,k-1)) / dz_phys_lo;
                } else {
                    dz_phys_hi = z_cc_arr(i,j  ,k+1) -   z_cc_arr(i,j  ,k-1);
                    dz_phys_lo = z_cc_arr(i,j-1,k+1) -   z_cc_arr(i,j-1,k-1);
                    gpz_hi  = (p_arr(i,j  ,k+1) - p_arr(i,j  ,k-1)) / dz_phys_hi;
                    gpz_lo  = (p_arr(i,j-1,k+1) - p_arr(i,j-1,k-1)) / dz_phys_lo;
                }
                Real gpy_metric = met_h_eta * 0.5 * (gpz_hi + gpz_lo);
                gpy -= gpy_metric;
            } // l_use_terrain_fitted_coords

            gpy_arr(i,j,k) = gpy;
        });

        ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real met_h_zeta = (l_use_terrain_fitted_coords) ? Compute_h_zeta_AtKface(i, j, k, dxInv, z_nd_arr) : 1;
            gpz_arr(i,j,k) = dxInv[2] * ( p_arr(i,j,k)-p_arr(i,j,k-1) )  / met_h_zeta;
        });
    } // mfi
}
