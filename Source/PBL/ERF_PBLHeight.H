#ifndef ERF_PBL_HEIGHT_H_
#define ERF_PBL_HEIGHT_H_

#include <AMReX_MultiFabUtil.H>
#include <ERF_TileNoZ.H>
#include <ERF_MoistUtils.H>

/**
 * @brief Diagnostic utility for the planetary boundary layer height.
 */
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
                      const MoistureComponentIndices& moisture_indices) const
    {
#if 0
        // NOTE: Cannot use ReduceToPlane because it clips the box to the
        // validbox only, i.e., lateral ghost cells aren't updated
        int dir = 2; // z
        auto const& cons_arrs = cons.const_arrays();
        auto thetav_min = amrex::ReduceToPlane<amrex::ReduceOpMin,amrex::Real>(dir, bxlow, cons,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) -> amrex::Real
            {
                return GetThetav(i,j,k,cons_arrs[box_no],moisture_indices);
            });
#endif

        // Create 2D multifabs like pblh
        auto const& ba = pblh->boxArray();
        auto const& dm = pblh->DistributionMap();
        auto const& ng = pblh->nGrowVect();

        amrex::MultiFab min_thetav(ba,dm,1,ng);
        min_thetav.setVal(bogus_large_value);

        amrex::MultiFab pblh_tke(ba,dm,1,ng);
        pblh_tke.setVal(zero);

        pblh->setVal(zero);

        // Now, loop over columns...
        for (amrex::MFIter mfi(cons,TileNoZ()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& domain = geom.Domain();

            // NOTE: this MFIter is over cons, which carries more ghost cells than the
            //       2D fields we write here (pblh, pblh_tke, min_thetav all have the
            //       ghost cells of pblh).  Grow by *their* ghost vector, not by the
            //       default, or we write past the end of their allocations.
            amrex::Box gtbx = mfi.growntilebox(ng);
            gtbx.setSmall(2,domain.smallEnd(2)); // don't loop over ghost cells
            gtbx.setBig(2,domain.bigEnd(2));     // in z

            // Columns are parallelized over; the k loops below are serialized within
            // each thread.  Everything in this routine is a per-column scan -- the
            // minimum theta_v, and the lowest k at which each PBLH criterion is met --
            // so parallelizing over k would both race on the 2D output and lose the
            // "lowest k wins" semantics.
            amrex::Box gtbx2d = gtbx; gtbx2d.setRange(2,0);
            const int klo = gtbx.smallEnd(2);
            const int khi = gtbx.bigEnd(2);

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
                int imin = lbound(zphys_arr).x;
                int jmin = lbound(zphys_arr).y;
                int imax = ubound(zphys_arr).x;
                int jmax = ubound(zphys_arr).y;

                // Find minimum thetav in the surface layer (this updates
                // ghost cells, too)
                ParallelFor(gtbx2d, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
                {
                    int ii = amrex::max(amrex::min(i,imax),imin);
                    int jj = amrex::max(amrex::min(j,jmax),jmin);

                    amrex::Real min_thv = min_thv_arr(i,j,0);
                    for (int k(klo); k <= khi; ++k) {
                        if (zphys_arr(ii,jj,k) < thetamin_height) {
                            amrex::Real thv = GetThetav(i, j, k, cons_arr, moisture_indices);
                            min_thv = amrex::min(min_thv, thv);
                        }
                    }
                    min_thv_arr(i,j,0) = min_thv;
                });

                // This depends on TileNoZ and k increasing monotonically
                ParallelFor(gtbx2d, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
                {
                    int ii = amrex::max(amrex::min(i,imax),imin);
                    int jj = amrex::max(amrex::min(j,jmax),jmin);

                    const amrex::Real min_thv = min_thv_arr(i,j,0);
                    const int is_land = (lmask_arr) ? lmask_arr(i,j,0) : 1;

                    amrex::Real zi     = pblh_arr(i,j,0);
                    amrex::Real zi_tke = pblh_tke_arr(i,j,0);

                    for (int k(klo); k <= khi; ++k)
                    {
                        if (zi == 0)
                        {
                            //
                            // Find PBL height based on thetav increase (best for CBLs)
                            //
                            amrex::Real thv  = GetThetav(i, j, k  , cons_arr, moisture_indices);
                            amrex::Real thv1 = GetThetav(i, j, k+1, cons_arr, moisture_indices);

                            if (is_land && (thv1 >= min_thv + theta_incr_land)
                                        && (thv  <  min_thv + theta_incr_land))
                            {
                                // Interpolate to get lowest height where theta = min_theta + theta_incr
                                zi = zphys_arr(ii,jj,k)
                                   + (zphys_arr(ii,jj,k+1)-zphys_arr(ii,jj,k))/(thv1-thv)
                                     * (min_thv + theta_incr_land - thv);
                            }
                            else if (!is_land && (thv1 >= min_thv + theta_incr_water)
                                              && (thv  <  min_thv + theta_incr_water))
                            {
                                // Interpolate to get lowest height where theta = min_theta + theta_incr
                                zi = zphys_arr(ii,jj,k)
                                   + (zphys_arr(ii,jj,k+1)-zphys_arr(ii,jj,k))/(thv1-thv)
                                     * (min_thv + theta_incr_water - thv);
                            }
                        }
                        if (zi_tke == 0)
                        {
                            //
                            // Find PBL height based on TKE (for SBLs only)
                            //
                            amrex::Real tke    = cons_arr(i,j,k  ,RhoKE_comp) / cons_arr(i,j,k  ,Rho_comp);
                            amrex::Real tke1   = cons_arr(i,j,k+1,RhoKE_comp) / cons_arr(i,j,k+1,Rho_comp);
                            amrex::Real maxtke = cons_arr(i,j,0  ,RhoKE_comp) / cons_arr(i,j,0  ,Rho_comp);
                            // - threshold is 5% of max TKE (Kosovic & Curry 2000, JAS)
                            amrex::Real TKEeps = amrex::Real(0.05) * maxtke;
                            TKEeps = amrex::max(TKEeps, amrex::Real(0.02)); // min val from WRF

                            if ((tke1 <= TKEeps) && (tke > TKEeps))
                            {
                                // Interpolate to get lowest height where TKE -> 0
                                zi_tke = zphys_arr(ii,jj,k)
                                       + (zphys_arr(ii,jj,k+1)-zphys_arr(ii,jj,k))/(tke1-tke)
                                         * (TKEeps - tke);
                            }
                        }
                        // Both criteria are satisfied; the remaining levels are no-ops
                        if ((zi != 0) && (zi_tke != 0)) { break; }
                    }

                    pblh_arr(i,j,0)     = zi;
                    pblh_tke_arr(i,j,0) = zi_tke;
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
                const int khi_low = amrex::min(kmax, khi);

                ParallelFor(gtbx2d, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
                {
                    amrex::Real min_thv = min_thv_arr(i,j,0);
                    for (int k(klo); k <= khi_low; ++k) {
                        amrex::Real thv = GetThetav(i, j, k, cons_arr, moisture_indices);
                        min_thv = amrex::min(min_thv, thv);
                    }
                    min_thv_arr(i,j,0) = min_thv;
                });

                // This depends on TileNoZ and k increasing monotonically
                ParallelFor(gtbx2d, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
                {
                    const amrex::Real min_thv = min_thv_arr(i,j,0);
                    const int is_land = (lmask_arr) ? lmask_arr(i,j,0) : 1;

                    amrex::Real zi     = pblh_arr(i,j,0);
                    amrex::Real zi_tke = pblh_tke_arr(i,j,0);

                    for (int k(klo); k <= khi; ++k)
                    {
                        if (zi == 0)
                        {
                            //
                            // Find PBL height based on thetav increase (best for CBLs)
                            //
                            amrex::Real thv  = GetThetav(i, j, k  , cons_arr, moisture_indices);
                            amrex::Real thv1 = GetThetav(i, j, k+1, cons_arr, moisture_indices);

                            if (is_land && (thv1 >= min_thv + theta_incr_land)
                                        && (thv  <  min_thv + theta_incr_land))
                            {
                                // Interpolate to get lowest height where theta = min_theta + theta_incr
                                zi = (k+myhalf)*dz_no_terrain
                                   + dz_no_terrain/(thv1-thv)
                                     * (min_thv + theta_incr_land - thv);
                            }
                            else if (!is_land && (thv1 >= min_thv + theta_incr_water)
                                              && (thv  <  min_thv + theta_incr_water))
                            {
                                // Interpolate to get lowest height where theta = min_theta + theta_incr
                                zi = (k+myhalf)*dz_no_terrain
                                   + dz_no_terrain/(thv1-thv)
                                     * (min_thv + theta_incr_water - thv);
                            }
                        }
                        if (zi_tke == 0)
                        {
                            //
                            // Find PBL height based on TKE (for SBLs only)
                            //
                            amrex::Real tke    = cons_arr(i,j,k  ,RhoKE_comp) / cons_arr(i,j,k  ,Rho_comp);
                            amrex::Real tke1   = cons_arr(i,j,k+1,RhoKE_comp) / cons_arr(i,j,k+1,Rho_comp);
                            amrex::Real maxtke = cons_arr(i,j,0  ,RhoKE_comp) / cons_arr(i,j,0  ,Rho_comp);
                            // - threshold is 5% of max TKE (Kosovic & Curry 2000, JAS)
                            amrex::Real TKEeps = amrex::Real(0.05) * maxtke;
                            TKEeps = amrex::max(TKEeps, amrex::Real(0.02)); // min val from WRF

                            if ((tke1 <= TKEeps) && (tke > TKEeps))
                            {
                                // Interpolate to get lowest height where TKE -> 0
                                zi_tke = (k+myhalf)*dz_no_terrain
                                       + dz_no_terrain/(tke1-tke) * (TKEeps - tke);
                            }
                        }
                        // Both criteria are satisfied; the remaining levels are no-ops
                        if ((zi != 0) && (zi_tke != 0)) { break; }
                    }

                    pblh_arr(i,j,0)     = zi;
                    pblh_tke_arr(i,j,0) = zi_tke;
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
                //   with convective precipitation (> 8km!), so an artificial limit is
                //   imposed to not let PBLH_TKE exceed the theta_v-based PBL height
                //   +/- 350 m. This has no impact on 98-99% of the domain, but is the
                //   simplest patch that adequately addresses these extremely large
                //   PBLHs.
                amrex::Real zi = pblh_arr(i,j,0);
                pblh_tke_arr(i,j,0) = amrex::max(
                    amrex::min(pblh_tke_arr(i,j,0), zi+amrex::Real(350.)),
                    amrex::max(zi-amrex::Real(350.), amrex::Real(10.)));

                //
                // Finally, blend between the two PBLH estimates
                //
                amrex::Real maxqke = two * cons_arr(i,j,0  ,RhoKE_comp) / cons_arr(i,j,0  ,Rho_comp);
                if (maxqke > amrex::Real(0.05)) {
                    amrex::Real wt = myhalf*std::tanh((zi - sbl_lim)/sbl_damp) + myhalf;
                    pblh_arr(i,j,0) = (one-wt)*pblh_tke_arr(i,j,0) + wt*pblh_arr(i,j,0);
                }
            });
        }//MFIter
    }

    /**
     * @var MYNNPBLH::thetamin_height
     * @brief Height below which minimum theta-v is determined [m].
     * @var MYNNPBLH::theta_incr_land
     * @brief Theta increase determining the capping inversion height over land [K].
     * @var MYNNPBLH::theta_incr_water
     * @brief Theta increase determining the capping inversion height over water [K].
     * @var MYNNPBLH::sbl_lim
     * @brief Upper limit of the stable boundary layer height [m].
     * @var MYNNPBLH::sbl_damp
     * @brief Transition length for blending [m].
     */
    static constexpr amrex::Real thetamin_height  = amrex::Real(200.0); // [m] height below which min thetav is determined
    static constexpr amrex::Real theta_incr_land  = amrex::Real(1.25);  // [K] theta increase that determines the height of the capping inversion over land
    static constexpr amrex::Real theta_incr_water = one;   // [K] theta increase that determines the height of the capping inversion over water
    static constexpr amrex::Real sbl_lim          = amrex::Real(200.0); // [m] upper limit of the stable BL height
    static constexpr amrex::Real sbl_damp         = amrex::Real(400.0); // [m] transition length for blending
};
#endif
