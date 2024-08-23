#ifndef _PBL_HEIGHT_H_
#define _PBL_HEIGHT_H_

#include <AMReX_MultiFabUtil.H>
#include <TileNoZ.H>
#include <Thetav.H>

struct MYNNPBLH {
    /*
     * Diagnose the PBL height
     *
     * Approach follows WRF, which uses a hybrid of the theta-increase
     * method for CBLs and a TKE threshold method for SBLs. The TKE method is
     * focused on PBL heights below ~500 m; above 1 km, the tanh blending makes
     * the TKE-based contribution negligible.
     *
     * See Nielsen-Gammon et al. 2008, JAS
     *     Banta 2008, Acta Geophys.
     */
    AMREX_GPU_HOST
    AMREX_FORCE_INLINE
    void compute_pblh(const amrex::Geometry& geom,
                      const amrex::MultiFab* z_phys_cc,
                      amrex::MultiFab* pblh,
                      const amrex::MultiFab& cons,
                      const amrex::iMultiFab* lmask,
                      const int moist_flag) const
    {
#if 0
        // NOTE: Cannot use ReduceToPlane because it clips the box to the
        // validbox only, i.e., lateral ghost cells aren't updated
        int dir = 2; // z
        auto const& cons_arrs = cons.const_arrays();
        auto thetav_min = amrex::ReduceToPlane<amrex::ReduceOpMin,amrex::Real>(dir, bxlow, cons,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) -> amrex::Real
            {
                return Thetav(i,j,k,cons_arrs[box_no],moist_flag);
            });
#endif

        // Create 2D multifabs like pblh
        auto const& ba = pblh->boxArray();
        auto const& dm = pblh->DistributionMap();
        auto const& ng = pblh->nGrowVect();

        amrex::MultiFab min_thetav(ba,dm,1,ng);
        min_thetav.setVal(1.E34);

        amrex::MultiFab pblh_tke(ba,dm,1,ng);
        pblh_tke.setVal(0);

        pblh->setVal(0);

        // Now, loop over columns...
        for (amrex::MFIter mfi(cons,TileNoZ()); mfi.isValid(); ++mfi)
        {
            amrex::Box gtbx = mfi.tilebox();
            gtbx.grow(ng); // grow in x,y only
            gtbx.setBig(2,gtbx.bigEnd(2)-1); // only loop to nz-1

            auto min_thv_arr  = min_thetav.array(mfi);
            auto pblh_arr     = pblh->array(mfi);
            auto pblh_tke_arr = pblh_tke.array(mfi);

            const auto cons_arr  = cons.const_array(mfi);
            const auto lmask_arr = (lmask) ? lmask->const_array(mfi) : amrex::Array4<int> {};

            // -----------------------------------------------------
            // WITH terrain/grid stretching
            // -----------------------------------------------------
            if (z_phys_cc)
            {
                const auto zphys_arr = z_phys_cc->const_array(mfi);

                // Need to sort out ghost cell differences (z_phys_cc has ng=1)
                const amrex::Box& domain = geom.Domain();
                amrex::IntVect ng_zphys = z_phys_cc->nGrowVect();

                int imin = domain.smallEnd(0); if (geom.isPeriodic(0)) imin -= ng_zphys[0];
                int jmin = domain.smallEnd(1); if (geom.isPeriodic(1)) jmin -= ng_zphys[1];

                int imax = domain.bigEnd(0); if (geom.isPeriodic(0)) imax += ng_zphys[0];
                int jmax = domain.bigEnd(1); if (geom.isPeriodic(1)) jmax += ng_zphys[1];

                // Find minimum thetav in the surface layer (this updates
                // ghost cells, too)
                ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    int ii = amrex::max(amrex::min(i,imax),imin);
                    int jj = amrex::max(amrex::min(j,jmax),jmin);

                    if (zphys_arr(ii,jj,k) < thetamin_height) {
                        amrex::Real thv = Thetav(i,j,k,cons_arr,moist_flag);
                        if (min_thv_arr(i,j,0) > thv) min_thv_arr(i,j,0) = thv;
                    }
                });

                // This depends on TileNoZ and k increasing monotonically
                ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    int ii = amrex::max(amrex::min(i,imax),imin);
                    int jj = amrex::max(amrex::min(j,jmax),jmin);

                    if (pblh_arr(i,j,0) == 0)
                    {
                        //
                        // Find PBL height based on thetav increase (best for CBLs)
                        //
                        amrex::Real thv  = Thetav(i,j,k  ,cons_arr,moist_flag);
                        amrex::Real thv1 = Thetav(i,j,k+1,cons_arr,moist_flag);

                        int is_land = (lmask_arr) ? lmask_arr(i,j,0) : 1;
                        if (is_land && (thv1 >= min_thv_arr(i,j,0) + theta_incr_land)
                                    && (thv  <  min_thv_arr(i,j,0) + theta_incr_land))
                        {
                            // Interpolate to get lowest height where theta = min_theta + theta_incr
                            pblh_arr(i,j,0) = zphys_arr(ii,jj,k)
                                            + (zphys_arr(ii,jj,k+1)-zphys_arr(ii,jj,k))/(thv1-thv)
                                              * (min_thv_arr(i,j,0) + theta_incr_land - thv);
                        }
                        else if (!is_land && (thv1 >= min_thv_arr(i,j,0) + theta_incr_water)
                                          && (thv  <  min_thv_arr(i,j,0) + theta_incr_water))
                        {
                            // Interpolate to get lowest height where theta = min_theta + theta_incr
                            pblh_arr(i,j,0) = zphys_arr(ii,jj,k)
                                            + (zphys_arr(ii,jj,k+1)-zphys_arr(ii,jj,k))/(thv1-thv)
                                              * (min_thv_arr(i,j,0) + theta_incr_water - thv);
                        }
                    }
                    if (pblh_tke_arr(i,j,0) == 0)
                    {
                        //
                        // Find PBL height based on TKE (for SBLs only)
                        //
                        amrex::Real tke    = 0.5 * cons_arr(i,j,k  ,RhoQKE_comp) / cons_arr(i,j,k  ,Rho_comp);
                        amrex::Real tke1   = 0.5 * cons_arr(i,j,k+1,RhoQKE_comp) / cons_arr(i,j,k+1,Rho_comp);
                        amrex::Real maxtke = 0.5 * cons_arr(i,j,0  ,RhoQKE_comp) / cons_arr(i,j,0  ,Rho_comp);
                        // - threshold is 5% of max TKE (Kosovic & Curry 2000, JAS)
                        amrex::Real TKEeps = 0.05 * maxtke;
                        TKEeps = amrex::max(TKEeps, 0.02); // min val from WRF

                        if ((tke1 <= TKEeps) && (tke > TKEeps))
                        {
                            // Interpolate to get lowest height where TKE -> 0
                            pblh_tke_arr(i,j,0) = zphys_arr(ii,jj,k)
                                                + (zphys_arr(ii,jj,k+1)-zphys_arr(ii,jj,k))/(tke1-tke)
                                                  * (TKEeps - tke);
                        }
                    }
                });
            }
            else
            // -----------------------------------------------------
            // NO terrain
            // -----------------------------------------------------
            {
                const amrex::Real dz_no_terrain = geom.CellSize(2);

                // Find minimum thetav in the surface layer (this updates
                // ghost cells, too)
                // - box size is known a priori
                int kmax = static_cast<int>(thetamin_height / dz_no_terrain);
                AMREX_ASSERT(kmax > 0);
                amrex::Box gtbxlow = gtbx;
                gtbxlow.setBig(2,kmax);

                ParallelFor(gtbxlow, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    amrex::Real thv = Thetav(i,j,k,cons_arr,moist_flag);
                    if (min_thv_arr(i,j,0) > thv) min_thv_arr(i,j,0) = thv;
                });

                // This depends on TileNoZ and k increasing monotonically
                ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if (pblh_arr(i,j,0) == 0)
                    {
                        //
                        // Find PBL height based on thetav increase (best for CBLs)
                        //
                        amrex::Real thv  = Thetav(i,j,k  ,cons_arr,moist_flag);
                        amrex::Real thv1 = Thetav(i,j,k+1,cons_arr,moist_flag);

                        int is_land = (lmask_arr) ? lmask_arr(i,j,0) : 1;
                        if (is_land && (thv1 >= min_thv_arr(i,j,0) + theta_incr_land)
                                    && (thv  <  min_thv_arr(i,j,0) + theta_incr_land))
                        {
                            // Interpolate to get lowest height where theta = min_theta + theta_incr
                            pblh_arr(i,j,0) = (k+0.5)*dz_no_terrain
                                            + dz_no_terrain/(thv1-thv)
                                              * (min_thv_arr(i,j,0) + theta_incr_land - thv);
                        }
                        else if (!is_land && (thv1 >= min_thv_arr(i,j,0) + theta_incr_water)
                                          && (thv  <  min_thv_arr(i,j,0) + theta_incr_water))
                        {
                            // Interpolate to get lowest height where theta = min_theta + theta_incr
                            pblh_arr(i,j,0) = (k+0.5)*dz_no_terrain
                                            + dz_no_terrain/(thv1-thv)
                                              * (min_thv_arr(i,j,0) + theta_incr_water - thv);
                        }
                    }
                    if (pblh_tke_arr(i,j,0) == 0)
                    {
                        //
                        // Find PBL height based on TKE (for SBLs only)
                        //
                        amrex::Real tke    = 0.5 * cons_arr(i,j,k  ,RhoQKE_comp) / cons_arr(i,j,k  ,Rho_comp);
                        amrex::Real tke1   = 0.5 * cons_arr(i,j,k+1,RhoQKE_comp) / cons_arr(i,j,k+1,Rho_comp);
                        amrex::Real maxtke = 0.5 * cons_arr(i,j,0  ,RhoQKE_comp) / cons_arr(i,j,0  ,Rho_comp);
                        // - threshold is 5% of max TKE (Kosovic & Curry 2000, JAS)
                        amrex::Real TKEeps = 0.05 * maxtke;
                        TKEeps = amrex::max(TKEeps, 0.02); // min val from WRF

                        if ((tke1 <= TKEeps) && (tke > TKEeps))
                        {
                            // Interpolate to get lowest height where TKE -> 0
                            pblh_tke_arr(i,j,0) = (k+0.5)*dz_no_terrain
                                                + dz_no_terrain/(tke1-tke) * (TKEeps - tke);
                        }
                    }
                });
            }
        }// MFIter

        //
        // Calculate hybrid PBL height
        //
        for (amrex::MFIter mfi(*pblh); mfi.isValid(); ++mfi)
        {
            const auto cons_arr = cons.const_array(mfi);
            auto pblh_tke_arr = pblh_tke.array(mfi);
            auto pblh_arr     = pblh->array(mfi);

            amrex::Box gtbx = mfi.growntilebox();
            ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
            {
                //
                // Clip PBLH_TKE to more realistic values
                //
                // Note from WRF MYNN-EDMF: TKE-based PBLH can be very large in cells
                //   with convective precipiation (> 8km!), so an artificial limit is
                //   imposed to not let PBLH_TKE exceed the theta_v-based PBL height
                //   +/- 350 m. This has no impact on 98-99% of the domain, but is the
                //   simplest patch that adequately addresses these extremely large
                //   PBLHs.
                amrex::Real zi = pblh_arr(i,j,0);
                pblh_tke_arr(i,j,0) = amrex::max(
                    amrex::min(pblh_tke_arr(i,j,0), zi+350.),
                    amrex::max(zi-350., 10.));

                //
                // Finally, blend between the two PBLH estimates
                //
                amrex::Real maxqke = cons_arr(i,j,0  ,RhoQKE_comp) / cons_arr(i,j,0  ,Rho_comp);
                if (maxqke > 0.05) {
                    amrex::Real wt = 0.5*std::tanh((zi - sbl_lim)/sbl_damp) + 0.5;
                    pblh_arr(i,j,0) = (1.-wt)*pblh_tke_arr(i,j,0) + wt*pblh_arr(i,j,0);
                }
            });
        }//MFIter
    }

    static constexpr amrex::Real thetamin_height  = 200.0; // [m] height below which min thetav is determined
    static constexpr amrex::Real theta_incr_land  = 1.25;  // [K] theta increase that determines the height of the capping inversion over land
    static constexpr amrex::Real theta_incr_water = 1.0;   // [K] theta increase that determines the height of the capping inversion over water
    static constexpr amrex::Real sbl_lim          = 200.0; // [m] upper limit of the stable BL height
    static constexpr amrex::Real sbl_damp         = 400.0; // [m] transition length for blending
};
#endif
