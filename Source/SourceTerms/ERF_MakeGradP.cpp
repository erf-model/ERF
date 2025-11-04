#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include "AMReX_BCRec.H"

#include "ERF.H"
#include "ERF_SrcHeaders.H"
#include "ERF_DataStruct.H"
#include "ERF_Utils.H"
#include "ERF_EB.H"
#include <ERF_EBSlopes.H>

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
 * @param[in]  d_bcrec_ptr Boundary Condition Record
 * @param[in]  ebfact    EB factory container at this level
 * @param[out] gradp     pressure gradient
 */

void make_gradp_pert (int level,
                      const SolverChoice& solverChoice,
                      const Geometry& geom,
                      Vector<MultiFab>& S_data,
                      const MultiFab& p0,
                      const MultiFab& z_phys_nd,
                      const MultiFab& z_phys_cc,
                      Vector<std::unique_ptr<MultiFab>>& mapfac,
                      const eb_& ebfact,
                      Vector<MultiFab>& gradp)
{
    const bool l_use_moisture  = (solverChoice.moisture_type != MoistureType::None);
    const bool l_eb_terrain    = (solverChoice.terrain_type == TerrainType::EB);
    //
    // Note that we only recompute gradp if compressible;
    //      if anelastic then we have computed gradp in the projection
    //      and we can reuse it, no need to recompute it
    //
    if (solverChoice.anelastic[level] == 0)
    {
        const int ngrow = (l_eb_terrain) ? 3 : 1;
        MultiFab p(S_data[Vars::cons].boxArray(), S_data[Vars::cons].DistributionMap(), 1, ngrow);

        // *****************************************************************************
        // Compute pressure or perturbataional pressure
        // *****************************************************************************
        for ( MFIter mfi(S_data[Vars::cons]); mfi.isValid(); ++mfi)
        {
            Box gbx = mfi.tilebox();
            gbx.grow(IntVect(ngrow,ngrow,ngrow));

            if (gbx.smallEnd(2) < 0) gbx.setSmall(2,0);
            const Array4<const Real>& cell_data = S_data[Vars::cons].array(mfi);
            const Array4<const Real>& p0_arr = p0.const_array(mfi);
            const Array4<      Real>& pptemp_arr = p.array(mfi);
            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real qv_for_p = (l_use_moisture) ? cell_data(i,j,k,RhoQ1_comp)/cell_data(i,j,k,Rho_comp) : 0.0;
                pptemp_arr(i,j,k) = getPgivenRTh(cell_data(i,j,k,RhoTheta_comp),qv_for_p) - p0_arr(i,j,k);
            });
        }

        if (solverChoice.gradp_type == 0) {
            compute_gradp(p,geom,z_phys_nd,z_phys_cc,mapfac,ebfact,gradp,solverChoice);
        } else if (solverChoice.gradp_type == 1) {
            AMREX_ASSERT_WITH_MESSAGE(solverChoice.terrain_type != TerrainType::EB,
                "gradp_type==1 not implemented for EB");
            compute_gradp_interpz(p,geom,z_phys_nd,z_phys_cc,mapfac,gradp,solverChoice);
        }
    }
}

void
compute_gradp (const MultiFab& p,
               const Geometry& geom,
               const MultiFab& z_phys_nd,
               const MultiFab& z_phys_cc,
               Vector<std::unique_ptr<MultiFab>>& mapfac,
               const eb_& ebfact,
               Vector<MultiFab>& gradp,
               const SolverChoice& solverChoice)
{
    const bool l_use_terrain_fitted_coords = (solverChoice.mesh_type != MeshType::ConstantDz);

    const Box domain = geom.Domain();
    const int domain_klo = domain.smallEnd(2);
    const int domain_khi = domain.bigEnd(2);

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
        if (tbz.smallEnd(2) == domain_klo) {
            tbz.growLo(2,-1);
        }
        if (tbz.bigEnd(2) == domain_khi+1) {
            tbz.growHi(2,-1);
        }

        // Terrain metrics
        const Array4<const Real>& z_nd_arr = z_phys_nd.const_array(mfi);
        const Array4<const Real>& z_cc_arr = z_phys_cc.const_array(mfi);

        const Array4<const Real>& p_arr = p.const_array(mfi);

        const Array4<      Real>& gpx_arr = gradp[GpVars::gpx].array(mfi);
        const Array4<      Real>& gpy_arr = gradp[GpVars::gpy].array(mfi);
        const Array4<      Real>& gpz_arr = gradp[GpVars::gpz].array(mfi);

        const Array4<const Real>& mf_ux_arr = mapfac[MapFacType::u_x]->const_array(mfi);
        const Array4<const Real>& mf_vy_arr = mapfac[MapFacType::v_y]->const_array(mfi);

        if (solverChoice.terrain_type != TerrainType::EB) {

            ParallelFor(tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                //Note : mx/my == 1, so no map factor needed here
                Real gpx = dxInv[0] * (p_arr(i,j,k) - p_arr(i-1,j,k));

                if (l_use_terrain_fitted_coords) {
                    Real met_h_xi = (z_cc_arr(i,j,k) - z_cc_arr(i-1,j,k)) * dxInv[0];

                    Real dz_phys_hi, dz_phys_lo;
                    Real gpz_lo, gpz_hi;
                    if (k==domain_klo) {
                        dz_phys_hi = z_cc_arr(i  ,j,k+1) -   z_cc_arr(i  ,j,k  );
                        dz_phys_lo = z_cc_arr(i-1,j,k+1) -   z_cc_arr(i-1,j,k  );
                        gpz_hi  = (p_arr(i  ,j,k+1) - p_arr(i  ,j,k  )) / dz_phys_hi;
                        gpz_lo  = (p_arr(i-1,j,k+1) - p_arr(i-1,j,k  )) / dz_phys_lo;
                    } else if (k==domain_khi) {
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

                // NOTE that the gradp array now carries the map factor!
                gpx_arr(i,j,k) *= mf_ux_arr(i,j,0);
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                //Note : mx/my == 1, so no map factor needed here
                Real gpy = dxInv[1] * (p_arr(i,j,k) - p_arr(i,j-1,k));

                if (l_use_terrain_fitted_coords) {
                    Real met_h_eta = (z_cc_arr(i,j,k) - z_cc_arr(i,j-1,k)) * dxInv[1];

                    Real dz_phys_hi, dz_phys_lo;
                    Real gpz_lo, gpz_hi;
                    if (k==domain_klo) {
                        dz_phys_hi = z_cc_arr(i,j  ,k+1) -   z_cc_arr(i,j  ,k  );
                        dz_phys_lo = z_cc_arr(i,j-1,k+1) -   z_cc_arr(i,j-1,k  );
                        gpz_hi  = (p_arr(i,j  ,k+1) - p_arr(i,j  ,k  )) / dz_phys_hi;
                        gpz_lo  = (p_arr(i,j-1,k+1) - p_arr(i,j-1,k  )) / dz_phys_lo;
                    } else if (k==domain_khi) {
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
                }
                gpy_arr(i,j,k) = gpy;

                // NOTE that the gradp array now carries the map factor!
                gpy_arr(i,j,k) *= mf_vy_arr(i,j,0);
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                Real met_h_zeta = (l_use_terrain_fitted_coords) ? Compute_h_zeta_AtKface(i, j, k, dxInv, z_nd_arr) : 1;
                gpz_arr(i,j,k) = dxInv[2] * ( p_arr(i,j,k)-p_arr(i,j,k-1) )  / met_h_zeta;
            });

        } else {

            // Pressure gradients are fitted at the centroids of cut cells, if EB and Compressible.
            // Least-Squares Fitting: Compute slope using 3x3x3 stencil

            const bool l_fitting = false;

            const Real* dx_arr = geom.CellSize();
            const Real dx = dx_arr[0];
            const Real dy = dx_arr[1];
            const Real dz = dx_arr[2];

            // EB factory
            Array4<const EBCellFlag> cellflg = (ebfact.get_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();

            // EB u-factory
            Array4<const EBCellFlag> u_cellflg = (ebfact.get_u_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
            Array4<const Real      > u_volfrac = (ebfact.get_u_const_factory())->getVolFrac().const_array(mfi);
            Array4<const Real      > u_volcent = (ebfact.get_u_const_factory())->getCentroid().const_array(mfi);

            // EB v-factory
            Array4<const EBCellFlag> v_cellflg = (ebfact.get_v_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
            Array4<const Real      > v_volfrac = (ebfact.get_v_const_factory())->getVolFrac().const_array(mfi);
            Array4<const Real      > v_volcent = (ebfact.get_v_const_factory())->getCentroid().const_array(mfi);

            // EB w-factory
            Array4<const EBCellFlag> w_cellflg = (ebfact.get_w_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
            Array4<const Real      > w_volfrac = (ebfact.get_w_const_factory())->getVolFrac().const_array(mfi);
            Array4<const Real      > w_volcent = (ebfact.get_w_const_factory())->getCentroid().const_array(mfi);

            if (l_fitting) {

                ParallelFor(tbx, tby, tbz,
                [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    if (u_volfrac(i,j,k) > 0.0) {

                        if (u_cellflg(i,j,k).isSingleValued()) {

                            GpuArray<Real,AMREX_SPACEDIM> slopes;
                            slopes = erf_calc_slopes_eb_staggered(Vars::xvel, Vars::cons, dx, dy, dz, i, j, k, p_arr, u_volcent, u_cellflg);

                            gpx_arr(i,j,k) = slopes[0];

                        } else {
                            gpx_arr(i,j,k) = dxInv[0] * (p_arr(i,j,k) - p_arr(i-1,j,k));
                        }

                    } else {
                        gpx_arr(i,j,k) = 0.0;
                    }
                },
                [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    if (v_volfrac(i,j,k) > 0.0) {

                        if (v_cellflg(i,j,k).isSingleValued()) {

                            GpuArray<Real,AMREX_SPACEDIM> slopes;
                            slopes = erf_calc_slopes_eb_staggered(Vars::yvel, Vars::cons, dx, dy, dz, i, j, k, p_arr, v_volcent, v_cellflg);

                            gpy_arr(i,j,k) = slopes[1];

                        } else {
                            gpy_arr(i,j,k) = dxInv[1] * (p_arr(i,j,k) - p_arr(i,j-1,k));
                        }
                    } else {
                        gpy_arr(i,j,k) = 0.0;
                    }
                },
                [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    if (w_volfrac(i,j,k) > 0.0) {

                        if (w_cellflg(i,j,k).isSingleValued()) {

                            GpuArray<Real,AMREX_SPACEDIM> slopes;
                            slopes = erf_calc_slopes_eb_staggered(Vars::zvel, Vars::cons, dx, dy, dz, i, j, k, p_arr, w_volcent, w_cellflg);

                            gpz_arr(i,j,k) = slopes[2];

                        } else {
                            gpz_arr(i,j,k) = dxInv[2] * (p_arr(i,j,k) - p_arr(i,j,k-1));
                        }
                    } else {
                        gpz_arr(i,j,k) = 0.0;
                    }
                });

            } else {

                // Simple calculation: assuming pressures at cell centers

                ParallelFor(tbx, tby, tbz,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if (u_volfrac(i,j,k) > 0.0) {
                        if (cellflg(i,j,k).isCovered()) {
                            gpx_arr(i,j,k) = dxInv[0] * (p_arr(i-3,j,k) - 3.*p_arr(i-2,j,k) + 2.*p_arr(i-1,j,k));
                        } else if (cellflg(i-1,j,k).isCovered()) {
                            gpx_arr(i,j,k) = dxInv[0] * (3.*p_arr(i+1,j,k) - p_arr(i+2,j,k) - 2.*p_arr(i,j,k));
                        } else {
                            gpx_arr(i,j,k) = dxInv[0] * (p_arr(i,j,k) - p_arr(i-1,j,k));
                        }
                    } else {
                        gpx_arr(i,j,k) = 0.0;
                    }
                },
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if (v_volfrac(i,j,k) > 0.0) {
                        if (cellflg(i,j,k).isCovered()) {
                            gpy_arr(i,j,k) = dxInv[1] * (p_arr(i,j-3,k) - 3.*p_arr(i,j-2,k) + 2.*p_arr(i,j-1,k));
                        } else if (cellflg(i,j-1,k).isCovered()) {
                            gpy_arr(i,j,k) = dxInv[1] * (3.*p_arr(i,j+1,k) - p_arr(i,j+2,k) - 2.*p_arr(i,j,k));
                        } else {
                            gpy_arr(i,j,k) = dxInv[1] * (p_arr(i,j,k) - p_arr(i,j-1,k));
                        }
                    } else {
                        gpy_arr(i,j,k) = 0.0;
                    }
                },
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if (w_volfrac(i,j,k) > 0.0) {
                        if (cellflg(i,j,k).isCovered()) {
                            gpz_arr(i,j,k) = dxInv[2] * ( p_arr(i,j,k-3) - 3.*p_arr(i,j,k-2) + 2.*p_arr(i,j,k-1) );
                        } else if (cellflg(i,j,k-1).isCovered()) {
                            gpz_arr(i,j,k) = dxInv[2] * ( 3.*p_arr(i,j,k+1) - p_arr(i,j,k+2) - 2.*p_arr(i,j,k) );
                        } else {
                            gpz_arr(i,j,k) = dxInv[2] * ( p_arr(i,j,k)-p_arr(i,j,k-1) );
                        }
                    } else {
                        gpz_arr(i,j,k) = 0.0;
                    }
                });

            } // l_fitting

        } // TerrainType::EB

    } // mfi
}

void
compute_gradp_interpz (const MultiFab& p,
                       const Geometry& geom,
                       const MultiFab& z_phys_nd,
                       const MultiFab& z_phys_cc,
                       Vector<std::unique_ptr<MultiFab>>& mapfac,
                       Vector<MultiFab>& gradp,
                       const SolverChoice& solverChoice)
{
    const bool l_use_terrain_fitted_coords = (solverChoice.mesh_type != MeshType::ConstantDz);

    const Box domain = geom.Domain();
    const int domain_klo = domain.smallEnd(2);
    const int domain_khi = domain.bigEnd(2);

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
        if (tbz.smallEnd(2) == domain_klo) {
            tbz.growLo(2,-1);
        }
        if (tbz.bigEnd(2) == domain_khi+1) {
            tbz.growHi(2,-1);
        }

        // Terrain metrics
        const Array4<const Real>& z_nd_arr = z_phys_nd.const_array(mfi);
        const Array4<const Real>& z_cc_arr = z_phys_cc.const_array(mfi);

        const Array4<const Real>& p_arr = p.const_array(mfi);

        const Array4<      Real>& gpx_arr = gradp[GpVars::gpx].array(mfi);
        const Array4<      Real>& gpy_arr = gradp[GpVars::gpy].array(mfi);
        const Array4<      Real>& gpz_arr = gradp[GpVars::gpz].array(mfi);

        const Array4<const Real>& mf_ux_arr = mapfac[MapFacType::u_x]->const_array(mfi);
        const Array4<const Real>& mf_vy_arr = mapfac[MapFacType::v_y]->const_array(mfi);

        ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if (l_use_terrain_fitted_coords) {
                Real p_lo = p_arr(i-1,j,k);
                Real p_hi = p_arr(i,j,k);
                Real dz_int = 0.5 * (z_cc_arr(i,j,k) - z_cc_arr(i-1,j,k));
                if (dz_int > 0) {
                    // Klemp 2011, Eqn. 16: s = 1/2
                    if (k==domain_klo) {
                        p_hi = quad_interp_1d(z_cc_arr(i,j,k) - dz_int,
                                              z_cc_arr(i,j,k  ), p_arr(i,j,k  ),
                                              z_cc_arr(i,j,k+1), p_arr(i,j,k+1),
                                              z_cc_arr(i,j,k+2), p_arr(i,j,k+2));
                    } else {
                        p_hi -= dz_int * ( (   p_arr(i  ,j,k  ) -    p_arr(i  ,j,k-1))
                                         / (z_cc_arr(i  ,j,k  ) - z_cc_arr(i  ,j,k-1)) );
                    }
                    if (k==domain_khi) {
                        p_lo = quad_interp_1d(z_cc_arr(i-1,j,k) + dz_int,
                                              z_cc_arr(i-1,j,k-2), p_arr(i-1,j,k-2),
                                              z_cc_arr(i-1,j,k-1), p_arr(i-1,j,k-1),
                                              z_cc_arr(i-1,j,k  ), p_arr(i-1,j,k  ));
                    } else {
                        p_lo += dz_int * ( (   p_arr(i-1,j,k+1) -    p_arr(i-1,j,k  ))
                                         / (z_cc_arr(i-1,j,k+1) - z_cc_arr(i-1,j,k  )) );
                    }
                } else if (dz_int < 0) {
                    // Klemp 2011, Eqn. 16: s = -1/2
                    if (k==domain_khi) {
                        p_hi = quad_interp_1d(z_cc_arr(i,j,k) - dz_int,
                                              z_cc_arr(i,j,k-2), p_arr(i,j,k-2),
                                              z_cc_arr(i,j,k-1), p_arr(i,j,k-1),
                                              z_cc_arr(i,j,k  ), p_arr(i,j,k  ));
                    } else {
                        p_hi -= dz_int * ( (   p_arr(i  ,j,k+1) -    p_arr(i  ,j,k  ))
                                         / (z_cc_arr(i  ,j,k+1) - z_cc_arr(i  ,j,k  )) );
                    }
                    if (k==domain_klo) {
                        p_lo = quad_interp_1d(z_cc_arr(i-1,j,k) + dz_int,
                                              z_cc_arr(i-1,j,k  ), p_arr(i-1,j,k  ),
                                              z_cc_arr(i-1,j,k+1), p_arr(i-1,j,k+1),
                                              z_cc_arr(i-1,j,k+2), p_arr(i-1,j,k+2));
                    } else {
                        p_lo += dz_int * ( (   p_arr(i-1,j,k  ) -    p_arr(i-1,j,k-1))
                                         / (z_cc_arr(i-1,j,k  ) - z_cc_arr(i-1,j,k-1)) );
                    }
                }
                gpx_arr(i,j,k) = dxInv[0] * (p_hi - p_lo);
            } else {
                gpx_arr(i,j,k) = dxInv[0] * (p_arr(i,j,k) - p_arr(i-1,j,k));
            }

            // NOTE that the gradp array now carries the map factor!
            gpx_arr(i,j,k) *= mf_ux_arr(i,j,0);
        },
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if (l_use_terrain_fitted_coords) {
                Real p_lo = p_arr(i,j-1,k);
                Real p_hi = p_arr(i,j,k);
                Real dz_int = 0.5 * (z_cc_arr(i,j,k) - z_cc_arr(i,j-1,k));
                if (dz_int > 0) {
                    // Klemp 2011, Eqn. 16: s = 1/2
                    if (k==domain_klo) {
                        p_hi = quad_interp_1d(z_cc_arr(i,j,k) - dz_int,
                                              z_cc_arr(i,j,k  ), p_arr(i,j,k  ),
                                              z_cc_arr(i,j,k+1), p_arr(i,j,k+1),
                                              z_cc_arr(i,j,k+2), p_arr(i,j,k+2));
                    } else {
                        p_hi -= dz_int * ( (   p_arr(i,j  ,k  ) -    p_arr(i,j  ,k-1))
                                         / (z_cc_arr(i,j  ,k  ) - z_cc_arr(i,j  ,k-1)) );
                    }
                    if (k==domain_khi) {
                        p_lo = quad_interp_1d(z_cc_arr(i,j-1,k) + dz_int,
                                              z_cc_arr(i,j-1,k-2), p_arr(i,j-1,k-2),
                                              z_cc_arr(i,j-1,k-1), p_arr(i,j-1,k-1),
                                              z_cc_arr(i,j-1,k  ), p_arr(i,j-1,k  ));
                    } else {
                        p_lo += dz_int * ( (   p_arr(i,j-1,k+1) -    p_arr(i,j-1,k  ))
                                         / (z_cc_arr(i,j-1,k+1) - z_cc_arr(i,j-1,k  )) );
                    }
                } else if (dz_int < 0) {
                    // Klemp 2011, Eqn. 16: s = -1/2
                    if (k==domain_khi) {
                        p_hi = quad_interp_1d(z_cc_arr(i,j,k) - dz_int,
                                              z_cc_arr(i,j,k-2), p_arr(i,j,k-2),
                                              z_cc_arr(i,j,k-1), p_arr(i,j,k-1),
                                              z_cc_arr(i,j,k  ), p_arr(i,j,k  ));
                    } else {
                        p_hi -= dz_int * ( (   p_arr(i,j  ,k+1) -    p_arr(i,j  ,k  ))
                                         / (z_cc_arr(i,j  ,k+1) - z_cc_arr(i,j  ,k  )) );
                    }
                    if (k==domain_klo) {
                        p_lo = quad_interp_1d(z_cc_arr(i,j-1,k) + dz_int,
                                              z_cc_arr(i,j-1,k  ), p_arr(i,j-1,k  ),
                                              z_cc_arr(i,j-1,k+1), p_arr(i,j-1,k+1),
                                              z_cc_arr(i,j-1,k+2), p_arr(i,j-1,k+2));
                    } else {
                        p_lo += dz_int * ( (   p_arr(i,j-1,k  ) -    p_arr(i,j-1,k-1))
                                         / (z_cc_arr(i,j-1,k  ) - z_cc_arr(i,j-1,k-1)) );
                    }
                }
                gpx_arr(i,j,k) = dxInv[1] * (p_hi - p_lo);
            } else {
                gpy_arr(i,j,k) = dxInv[1] * (p_arr(i,j,k) - p_arr(i,j-1,k));
            }

            // NOTE that the gradp array now carries the map factor!
            gpy_arr(i,j,k) *= mf_vy_arr(i,j,0);
        },
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Note: identical to gradp_type == 0
            Real met_h_zeta = (l_use_terrain_fitted_coords) ? Compute_h_zeta_AtKface(i, j, k, dxInv, z_nd_arr) : 1;
            gpz_arr(i,j,k) = dxInv[2] * ( p_arr(i,j,k)-p_arr(i,j,k-1) )  / met_h_zeta;
        });
    } // mfi
}
