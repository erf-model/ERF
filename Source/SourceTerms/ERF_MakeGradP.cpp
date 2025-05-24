#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_EB_Slopes_K.H>
#include "AMReX_BCRec.H"

#include "ERF.H"
#include "ERF_SrcHeaders.H"
#include "ERF_DataStruct.H"
#include "ERF_Utils.H"
#include "ERF_EB.H"

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
                      MultiFab& p0,
                      MultiFab& z_phys_nd,
                      MultiFab& z_phys_cc,
                      BCRec const* d_bcrec_ptr,
                      const eb_& ebfact,
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
        MultiFab p(S_data[Vars::cons].boxArray(), S_data[Vars::cons].DistributionMap(), 1, 1);

        // *****************************************************************************
        // Compute pressure or perturbataional pressure
        // *****************************************************************************
        for ( MFIter mfi(S_data[Vars::cons]); mfi.isValid(); ++mfi)
        {
            Box gbx = mfi.tilebox(); gbx.grow(IntVect(1,1,1));
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

        compute_gradp(p,geom,z_phys_nd,z_phys_cc,d_bcrec_ptr,ebfact,gradp,solverChoice);
    }
}

void
compute_gradp (const MultiFab& p,
               const Geometry& geom,
               MultiFab& z_phys_nd,
               MultiFab& z_phys_cc,
               BCRec const* d_bcrec_ptr,
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

        if (solverChoice.terrain_type != TerrainType::EB) {

            ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
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
            });

            ParallelFor(tby, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
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
            });

            ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                Real met_h_zeta = (l_use_terrain_fitted_coords) ? Compute_h_zeta_AtKface(i, j, k, dxInv, z_nd_arr) : 1;
                gpz_arr(i,j,k) = dxInv[2] * ( p_arr(i,j,k)-p_arr(i,j,k-1) )  / met_h_zeta;
            });

        } else {

            // Pressure gradients are fitted at the centroids of cut cells, if EB and Compressible.
            // Least-Squares Fitting: Compute slope using 3x3x3 stencil

            const bool l_fitting = false;

            const int domain_ilo = domain.smallEnd(0);
            const int domain_ihi = domain.bigEnd(0);
            const int domain_jlo = domain.smallEnd(1);
            const int domain_jhi = domain.bigEnd(1);

            int n = 0;

            // Need to check d_bcrec_ptr[n]

            bool extdir_ilo = (d_bcrec_ptr[n].lo(0)==BCType::ext_dir || d_bcrec_ptr[n].lo(0)==BCType::hoextrap);
            bool extdir_ihi = (d_bcrec_ptr[n].hi(0)==BCType::ext_dir || d_bcrec_ptr[n].hi(0)==BCType::hoextrap);
            bool extdir_jlo = (d_bcrec_ptr[n].lo(1)==BCType::ext_dir || d_bcrec_ptr[n].lo(1)==BCType::hoextrap);
            bool extdir_jhi = (d_bcrec_ptr[n].hi(1)==BCType::ext_dir || d_bcrec_ptr[n].hi(1)==BCType::hoextrap);
            bool extdir_klo = (d_bcrec_ptr[n].lo(2)==BCType::ext_dir || d_bcrec_ptr[n].lo(2)==BCType::hoextrap);
            bool extdir_khi = (d_bcrec_ptr[n].hi(2)==BCType::ext_dir || d_bcrec_ptr[n].hi(2)==BCType::hoextrap);

            // Pressure in staggered grids
            MultiFab u_p((ebfact.get_u_const_factory())->getVolFrac().boxArray(), (ebfact.get_u_const_factory())->getVolFrac().DistributionMap(), 1, 1);
            MultiFab v_p((ebfact.get_v_const_factory())->getVolFrac().boxArray(), (ebfact.get_v_const_factory())->getVolFrac().DistributionMap(), 1, 1);
            MultiFab w_p((ebfact.get_w_const_factory())->getVolFrac().boxArray(), (ebfact.get_w_const_factory())->getVolFrac().DistributionMap(), 1, 1);

            const Array4<Real>& u_p_arr = u_p.array(mfi);
            const Array4<Real>& v_p_arr = v_p.array(mfi);
            const Array4<Real>& w_p_arr = w_p.array(mfi);

            // EB u-factory
            Array4<const EBCellFlag> u_cflag = (ebfact.get_u_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
            Array4<const Real      > u_vfrac = (ebfact.get_u_const_factory())->getVolFrac().const_array(mfi);
            Array4<const Real      > u_vcent = (ebfact.get_u_const_factory())->getCentroid().const_array(mfi);
            Array4<const Real      > u_fcx   = (ebfact.get_u_const_factory())->getFaceCent()[0]->const_array(mfi);
            Array4<const Real      > u_fcy   = (ebfact.get_u_const_factory())->getFaceCent()[1]->const_array(mfi);
            Array4<const Real      > u_fcz   = (ebfact.get_u_const_factory())->getFaceCent()[2]->const_array(mfi);

            // EB v-factory
            Array4<const EBCellFlag> v_cflag = (ebfact.get_v_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
            Array4<const Real      > v_vfrac = (ebfact.get_v_const_factory())->getVolFrac().const_array(mfi);
            Array4<const Real      > v_vcent = (ebfact.get_v_const_factory())->getCentroid().const_array(mfi);
            Array4<const Real      > v_fcx   = (ebfact.get_v_const_factory())->getFaceCent()[0]->const_array(mfi);
            Array4<const Real      > v_fcy   = (ebfact.get_v_const_factory())->getFaceCent()[1]->const_array(mfi);
            Array4<const Real      > v_fcz   = (ebfact.get_v_const_factory())->getFaceCent()[2]->const_array(mfi);

            // EB w-factory
            Array4<const EBCellFlag> w_cflag = (ebfact.get_w_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
            Array4<const Real      > w_vfrac = (ebfact.get_w_const_factory())->getVolFrac().const_array(mfi);
            Array4<const Real      > w_vcent = (ebfact.get_w_const_factory())->getCentroid().const_array(mfi);
            Array4<const Real      > w_fcx   = (ebfact.get_w_const_factory())->getFaceCent()[0]->const_array(mfi);
            Array4<const Real      > w_fcy   = (ebfact.get_w_const_factory())->getFaceCent()[1]->const_array(mfi);
            Array4<const Real      > w_fcz   = (ebfact.get_w_const_factory())->getFaceCent()[2]->const_array(mfi);

            if (l_fitting) {

                // STEP 1: Compute pressure in the staggered grids

                const Box tbx_g1 = amrex::grow(tbx,1);
                const Box tby_g1 = amrex::grow(tby,1);
                const Box tbz_g1 = amrex::grow(tbz,1);

                int tbx_g1_ilo = tbx_g1.smallEnd(0);
                int tby_g1_jlo = tby_g1.smallEnd(1);
                int tbz_g1_klo = tbz_g1.smallEnd(2);
                int tbx_g1_ihi = tbx_g1.bigEnd(0);
                int tby_g1_jhi = tby_g1.bigEnd(1);
                int tbz_g1_khi = tbz_g1.bigEnd(2);

                // Print()<<"SK: tbx    = "<<tbx   <<", tby    = "<<tby   <<", tbz    = "<<tbz<<std::endl;
                // Print()<<"SK: tbx_g1 = "<<tbx_g1<<", tby_g1 = "<<tby_g1<<", tbz_g1 = "<<tbz_g1<<std::endl;

                ParallelFor(tbx_g1, tby_g1, tbz_g1,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (u_vfrac(i,j,k) > 0.0) {
                        // Print()<<"SK: tbx_g1/ i,j,k = "<<i<<" "<<j<<" "<<k<<" "<<std::endl;
                        if (i == tbx_g1_ilo) {
                            u_p_arr(i,j,k) = p_arr(i,j,k);
                        } else if (i == tbx_g1_ihi) {
                            u_p_arr(i,j,k) = p_arr(i-1,j,k);
                        } else {
                            u_p_arr(i,j,k) = 0.5 * ( p_arr(i-1,j,k) + p_arr(i,j,k) );
                        }
                    } else {
                        u_p_arr(i,j,k) = 0.0;
                    }
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (u_vfrac(i,j,k) > 0.0) {
                        // Print()<<"SK: tby_g1/ i,j,k = "<<i<<" "<<j<<" "<<k<<" "<<std::endl;
                        if (j == tby_g1_jlo) {
                            v_p_arr(i,j,k) = p_arr(i,j,k);
                        } else if (j == tby_g1_jhi) {
                            v_p_arr(i,j,k) = p_arr(i,j-1,k);
                        } else {
                            v_p_arr(i,j,k) = 0.5 * ( p_arr(i,j-1,k) + p_arr(i,j,k) );
                        }

                    } else {
                        v_p_arr(i,j,k) = 0.0;
                    }
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (w_vfrac(i,j,k) > 0.0) {
                        if (k == tbz_g1_klo) {
                            w_p_arr(i,j,k) = p_arr(i,j,k);
                        } else if (k == tbz_g1_khi) {
                            w_p_arr(i,j,k) = p_arr(i,j,k-1);
                        } else {
                            w_p_arr(i,j,k) = 0.5 * ( p_arr(i,j,k-1) + p_arr(i,j,k) );
                        }
                    } else {
                        w_p_arr(i,j,k) = 0.0;
                    }
                });

                // STEP 2: Compute Least-Squares slopes

                ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if (u_vfrac(i,j,k) > 0.0) {
                        GpuArray<Real,AMREX_SPACEDIM> slopes_eb;

                        slopes_eb = amrex_calc_slopes_extdir_eb(
                            i, j, k, n, u_p_arr, u_vcent, u_vfrac,
                            AMREX_D_DECL(u_fcx,u_fcy,u_fcz),u_cflag,
                            AMREX_D_DECL(extdir_ilo, extdir_jlo, extdir_klo),
                            AMREX_D_DECL(extdir_ihi, extdir_jhi, extdir_khi),
                            AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                            AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                            2);

                        gpx_arr(i,j,k) = slopes_eb[0];
                    } else {
                        gpx_arr(i,j,k) = 0.0;
                    }
                });

                ParallelFor(tby, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    // Need to interpolate pressure from CC to FC grid here.

                    if (v_vfrac(i,j,k) > 0.0) {
                        GpuArray<Real,AMREX_SPACEDIM> slopes_eb;

                        slopes_eb = amrex_calc_slopes_extdir_eb(
                            i, j, k, n, v_p_arr, v_vcent, v_vfrac,
                            AMREX_D_DECL(v_fcx,v_fcy,v_fcz),v_cflag,
                            AMREX_D_DECL(extdir_ilo, extdir_jlo, extdir_klo),
                            AMREX_D_DECL(extdir_ihi, extdir_jhi, extdir_khi),
                            AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                            AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                            2);

                        gpy_arr(i,j,k) = slopes_eb[1];
                    } else {
                        gpy_arr(i,j,k) = 0.0;
                    }
                });

                ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    // Need to interpolate pressure from CC to FC grid here.

                    if (w_vfrac(i,j,k) > 0.0) {
                        GpuArray<Real,AMREX_SPACEDIM> slopes_eb;

                        slopes_eb = amrex_calc_slopes_extdir_eb(
                            i, j, k, n, w_p_arr, w_vcent, w_vfrac,
                            AMREX_D_DECL(w_fcx,w_fcy,w_fcz),w_cflag,
                            AMREX_D_DECL(extdir_ilo, extdir_jlo, extdir_klo),
                            AMREX_D_DECL(extdir_ihi, extdir_jhi, extdir_khi),
                            AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                            AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                            2);

                        gpz_arr(i,j,k) = slopes_eb[2];
                    } else {
                        gpz_arr(i,j,k) = 0.0;
                    }
                });

            } else {

                // Simple calculation: assuming pressures at cell centers

                ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if (u_vfrac(i,j,k) > 0.0) {
                        gpx_arr(i,j,k) = dxInv[0] * (p_arr(i,j,k) - p_arr(i-1,j,k));
                    } else {
                        gpx_arr(i,j,k) = 0.0;
                    }
                });

                ParallelFor(tby, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if (v_vfrac(i,j,k) > 0.0) {
                        gpy_arr(i,j,k) = dxInv[1] * (p_arr(i,j,k) - p_arr(i,j-1,k));
                    } else {
                        gpy_arr(i,j,k) = 0.0;
                    }
                });

                ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if (w_vfrac(i,j,k) > 0.0) {
                        gpz_arr(i,j,k) = dxInv[2] * ( p_arr(i,j,k)-p_arr(i,j,k-1) )  ;
                    } else {
                        gpz_arr(i,j,k) = 0.0;
                    }
                });

            } // l_fitting

        } // TerrainType::EB

    } // mfi
}
