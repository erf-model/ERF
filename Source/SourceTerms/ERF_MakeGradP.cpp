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
 * @param[in]  pp_inc    pressure perturbation if anelastic
 * @param[in]  z_phys_nd z on nodes
 * @param[in]  z_phys_cc z on cell centers
 * @param[out] gradp     pressure gradient
 */

void make_gradp (int level,
                 const SolverChoice& solverChoice,
                 const Geometry& geom,
                 MultiFab& S_data,
                 MultiFab& p0,
                 const MultiFab& pp_inc,
                 std::unique_ptr<MultiFab>& z_phys_nd,
                 std::unique_ptr<MultiFab>& z_phys_cc,
                 BCRec const* d_bcrec_ptr,
                 const eb_& ebfact,
                 Vector<MultiFab>& gradp)
{
    const bool l_use_terrain_fitted_coords = (solverChoice.mesh_type != MeshType::ConstantDz);
    const bool l_use_moisture  = (solverChoice.moisture_type != MoistureType::None);
    const bool l_anelastic = solverChoice.anelastic[level];

    const Box domain = geom.Domain();
    const int domain_klo = domain.smallEnd(2);
    const int domain_khi = domain.bigEnd(2);

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();

    for ( MFIter mfi(S_data); mfi.isValid(); ++mfi)
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

        const Array4<const Real>& cell_data = S_data.array(mfi);

        const Array4<      Real>& gpx_arr = gradp[GpVars::gpx].array(mfi);
        const Array4<      Real>& gpy_arr = gradp[GpVars::gpy].array(mfi);
        const Array4<      Real>& gpz_arr = gradp[GpVars::gpz].array(mfi);

        // Base state
        const Array4<const Real>& p0_arr = p0.const_array(mfi);

        // Terrain metrics
        const Array4<const Real>& z_nd = z_phys_nd->const_array(mfi);
        const Array4<const Real>& z_cc = z_phys_cc->const_array(mfi);

        // *****************************************************************************
        // Perturbational pressure field
        // *****************************************************************************
        FArrayBox pprime;
        if (!l_anelastic) {
            Box gbx = mfi.tilebox(); gbx.grow(IntVect(1,1,1));
            if (gbx.smallEnd(2) < 0) gbx.setSmall(2,0);
            pprime.resize(gbx,1,The_Async_Arena());
            const Array4<Real>& pptemp_arr = pprime.array();
            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
#ifdef AMREX_USE_GPU
                if (cell_data(i,j,k,RhoTheta_comp) <= 0.) AMREX_DEVICE_PRINTF("BAD THETA AT %d %d %d %e %e \n",
                    i,j,k,cell_data(i,j,k,RhoTheta_comp),cell_data(i,j,k+1,RhoTheta_comp));
#else
                if (cell_data(i,j,k,RhoTheta_comp) <= 0.) {
                    printf("BAD THETA AT %d %d %d %e %e \n",
                    i,j,k,cell_data(i,j,k,RhoTheta_comp),cell_data(i,j,k+1,RhoTheta_comp));
                    Abort("Bad theta in ERF_slow_rhs_pre");
                }
#endif
                Real qv_for_p = (l_use_moisture) ? cell_data(i,j,k,RhoQ1_comp)/cell_data(i,j,k,Rho_comp) : 0.0;
                pptemp_arr(i,j,k) = getPgivenRTh(cell_data(i,j,k,RhoTheta_comp),qv_for_p) - p0_arr(i,j,k);
            });
        }

        const Array4<const Real>& pp_arr = (l_anelastic) ? pp_inc.const_array(mfi) : pprime.const_array();

        // Only if EB and Compressible, pressure gradients are fitted at the centroids of cut cells.

        if (solverChoice.terrain_type != TerrainType::EB || l_anelastic) {

            ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                //Note : mx/my == 1, so no map factor needed here
                Real gpx = dxInv[0] * (pp_arr(i,j,k) - pp_arr(i-1,j,k));

                if (l_use_terrain_fitted_coords) {
                    Real met_h_xi = (z_cc(i,j,k) - z_cc(i-1,j,k)) * dxInv[0];

                    Real dz_phys_hi, dz_phys_lo;
                    Real gpz_lo, gpz_hi;
                    if (k==domain_klo) {
                        dz_phys_hi = z_cc(i  ,j,k+1) -   z_cc(i  ,j,k  );
                        dz_phys_lo = z_cc(i-1,j,k+1) -   z_cc(i-1,j,k  );
                        gpz_hi  = (pp_arr(i  ,j,k+1) - pp_arr(i  ,j,k  )) / dz_phys_hi;
                        gpz_lo  = (pp_arr(i-1,j,k+1) - pp_arr(i-1,j,k  )) / dz_phys_lo;
                    } else if (k==domain_khi) {
                        dz_phys_hi = z_cc(i  ,j,k  ) -   z_cc(i  ,j,k-1);
                        dz_phys_lo = z_cc(i-1,j,k  ) -   z_cc(i-1,j,k-1);
                        gpz_hi  = (pp_arr(i  ,j,k  ) - pp_arr(i  ,j,k-1)) / dz_phys_hi;
                        gpz_lo  = (pp_arr(i-1,j,k  ) - pp_arr(i-1,j,k-1)) / dz_phys_lo;
                    } else {
                        dz_phys_hi = z_cc(i  ,j,k+1) -   z_cc(i  ,j,k-1);
                        dz_phys_lo = z_cc(i-1,j,k+1) -   z_cc(i-1,j,k-1);
                        gpz_hi  = (pp_arr(i  ,j,k+1) - pp_arr(i  ,j,k-1)) / dz_phys_hi;
                        gpz_lo  = (pp_arr(i-1,j,k+1) - pp_arr(i-1,j,k-1)) / dz_phys_lo;
                    }
                    Real gpx_metric = met_h_xi * 0.5 * (gpz_hi + gpz_lo);
                    gpx -= gpx_metric;
                }
                gpx_arr(i,j,k) = gpx;
            });

            ParallelFor(tby, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                //Note : mx/my == 1, so no map factor needed here
                Real gpy = dxInv[1] * (pp_arr(i,j,k) - pp_arr(i,j-1,k));

                if (l_use_terrain_fitted_coords) {
                    Real met_h_eta = (z_cc(i,j,k) - z_cc(i,j-1,k)) * dxInv[1];

                    Real dz_phys_hi, dz_phys_lo;
                    Real gpz_lo, gpz_hi;
                    if (k==domain_klo) {
                        dz_phys_hi = z_cc(i,j  ,k+1) -   z_cc(i,j  ,k  );
                        dz_phys_lo = z_cc(i,j-1,k+1) -   z_cc(i,j-1,k  );
                        gpz_hi  = (pp_arr(i,j  ,k+1) - pp_arr(i,j  ,k  )) / dz_phys_hi;
                        gpz_lo  = (pp_arr(i,j-1,k+1) - pp_arr(i,j-1,k  )) / dz_phys_lo;
                    } else if (k==domain_khi) {
                        dz_phys_hi = z_cc(i,j  ,k  ) -   z_cc(i,j  ,k-1);
                        dz_phys_lo = z_cc(i,j-1,k  ) -   z_cc(i,j-1,k-1);
                        gpz_hi  = (pp_arr(i,j  ,k  ) - pp_arr(i,j  ,k-1)) / dz_phys_hi;
                        gpz_lo  = (pp_arr(i,j-1,k  ) - pp_arr(i,j-1,k-1)) / dz_phys_lo;
                    } else {
                        dz_phys_hi = z_cc(i,j  ,k+1) -   z_cc(i,j  ,k-1);
                        dz_phys_lo = z_cc(i,j-1,k+1) -   z_cc(i,j-1,k-1);
                        gpz_hi  = (pp_arr(i,j  ,k+1) - pp_arr(i,j  ,k-1)) / dz_phys_hi;
                        gpz_lo  = (pp_arr(i,j-1,k+1) - pp_arr(i,j-1,k-1)) / dz_phys_lo;
                    }
                    Real gpy_metric = met_h_eta * 0.5 * (gpz_hi + gpz_lo);
                    gpy -= gpy_metric;
                } // l_use_terrain_fitted_coords

                gpy_arr(i,j,k) = gpy;
            });

            ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                Real met_h_zeta = (l_use_terrain_fitted_coords) ? Compute_h_zeta_AtKface(i, j, k, dxInv, z_nd) : 1;
                gpz_arr(i,j,k) = dxInv[2] * ( pp_arr(i,j,k)-pp_arr(i,j,k-1) )  / met_h_zeta;
            });

        } else {

            // Least-Squares fitting of pressure gradient for Compressible-EB case
            // Compute slope using 3x3x3 stencil

            // BCRec const* d_bcrec_ptr = domain_bcs_type_d.data();

            // for (int n = 0; n < AMREX_SPACEDIM+NBCVAR_max; n++) {
            //     Print()<<"SK: MakeGradP/ n = "<< n << ", d_bcrec_ptr[n].lo(0) = "<< d_bcrec_ptr[n].lo(0) << std::endl;
            //     Print()<<"SK: MakeGradP/ n = "<< n << ", d_bcrec_ptr[n].hi(0) = "<< d_bcrec_ptr[n].hi(0) << std::endl;
            //     Print()<<"SK: MakeGradP/ n = "<< n << ", d_bcrec_ptr[n].lo(1) = "<< d_bcrec_ptr[n].lo(1) << std::endl;
            //     Print()<<"SK: MakeGradP/ n = "<< n << ", d_bcrec_ptr[n].hi(1) = "<< d_bcrec_ptr[n].hi(1) << std::endl;
            //     Print()<<"SK: MakeGradP/ n = "<< n << ", d_bcrec_ptr[n].lo(2) = "<< d_bcrec_ptr[n].lo(2) << std::endl;
            //     Print()<<"SK: MakeGradP/ n = "<< n << ", d_bcrec_ptr[n].hi(2) = "<< d_bcrec_ptr[n].hi(2) << std::endl;
            // }
            

            // bool extdir_ilo = (d_bcrec_ptr[n].lo(0)==BCType::ext_dir || d_bcrec_ptr[n].lo(0)==BCType::hoextrap);
            // bool extdir_ihi = (d_bcrec_ptr[n].hi(0)==BCType::ext_dir || d_bcrec_ptr[n].hi(0)==BCType::hoextrap);
            // bool extdir_jlo = (d_bcrec_ptr[n].lo(1)==BCType::ext_dir || d_bcrec_ptr[n].lo(1)==BCType::hoextrap);
            // bool extdir_jhi = (d_bcrec_ptr[n].hi(1)==BCType::ext_dir || d_bcrec_ptr[n].hi(1)==BCType::hoextrap);
            // bool extdir_klo = (d_bcrec_ptr[n].lo(2)==BCType::ext_dir || d_bcrec_ptr[n].lo(2)==BCType::hoextrap);
            // bool extdir_khi = (d_bcrec_ptr[n].hi(2)==BCType::ext_dir || d_bcrec_ptr[n].hi(2)==BCType::hoextrap);

            // const int domain_ilo = domain.smallEnd(0);
            // const int domain_ihi = domain.bigEnd(0);
            // const int domain_jlo = domain.smallEnd(1);
            // const int domain_jhi = domain.bigEnd(1);

            // int n = 0;

            // // EB u-factory
            // Array4<const EBCellFlag> u_cflag = (ebfact.get_u_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
            // Array4<const Real      > u_vfrac = (ebfact.get_u_const_factory())->getVolFrac().const_array(mfi);
            // Array4<const Real      > u_vcent = (ebfact.get_u_const_factory())->getCentroid().const_array(mfi);
            // Array4<const Real      > u_fcx   = (ebfact.get_u_const_factory())->getFaceCent()[0]->const_array(mfi);
            // Array4<const Real      > u_fcy   = (ebfact.get_u_const_factory())->getFaceCent()[1]->const_array(mfi);
            // Array4<const Real      > u_fcz   = (ebfact.get_u_const_factory())->getFaceCent()[2]->const_array(mfi);

            // // EB v-factory
            // Array4<const EBCellFlag> v_cflag = (ebfact.get_v_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
            // Array4<const Real      > v_vfrac = (ebfact.get_v_const_factory())->getVolFrac().const_array(mfi);
            // Array4<const Real      > v_vcent = (ebfact.get_v_const_factory())->getCentroid().const_array(mfi);
            // Array4<const Real      > v_fcx   = (ebfact.get_v_const_factory())->getFaceCent()[0]->const_array(mfi);
            // Array4<const Real      > v_fcy   = (ebfact.get_v_const_factory())->getFaceCent()[1]->const_array(mfi);
            // Array4<const Real      > v_fcz   = (ebfact.get_v_const_factory())->getFaceCent()[2]->const_array(mfi);

            // // EB w-factory
            // Array4<const EBCellFlag> w_cflag = (ebfact.get_w_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
            // Array4<const Real      > w_vfrac = (ebfact.get_w_const_factory())->getVolFrac().const_array(mfi);
            // Array4<const Real      > w_vcent = (ebfact.get_w_const_factory())->getCentroid().const_array(mfi);
            // Array4<const Real      > w_fcx   = (ebfact.get_w_const_factory())->getFaceCent()[0]->const_array(mfi);
            // Array4<const Real      > w_fcy   = (ebfact.get_w_const_factory())->getFaceCent()[1]->const_array(mfi);
            // Array4<const Real      > w_fcz   = (ebfact.get_w_const_factory())->getFaceCent()[2]->const_array(mfi);

            // GpuArray<Real,AMREX_SPACEDIM> slopes_eb;

            // ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            // {

            //     // Interpolate pressure from CC to FC grid.

            //     if (vfrac(i,j,k) > 0.0)
            //     {
                    

            //         slopes_eb = amrex_calc_slopes_extdir_eb(
            //             i,j,k,n,pp_arr,u_vcent,u_vfrac,
            //             AMREX_D_DECL(u_fcx,u_fcy,u_fcz),u_cflag,
            //             AMREX_D_DECL(extdir_ilo, extdir_jlo, extdir_klo),
            //             AMREX_D_DECL(extdir_ihi, extdir_jhi, extdir_khi),
            //             AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
            //             AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
            //             2);

            //     }
            // });
        }

    }
}
