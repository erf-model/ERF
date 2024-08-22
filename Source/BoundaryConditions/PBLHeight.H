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
     * method for CBLs and a TKE threshold method for SBLs.
     *
     * See Nielsen-Gammon et al. 2008, JAS
     *     Banta 2008, Acta Geophys.
     */
    AMREX_GPU_HOST
    AMREX_FORCE_INLINE
    void compute_pblh(const amrex::Box& domain,
                      const amrex::Real dz_no_terrain,
                      const amrex::MultiFab* z_phys_cc,
                      amrex::MultiFab* pblh,
                      const amrex::MultiFab* u_star,
                      const amrex::MultiFab* t_star,
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
        const amrex::BoxArray ba             = pblh->boxArray();
        const amrex::DistributionMapping& dm = pblh->DistributionMap();
        //const amrex::IntVect ng              = pblh->nGrowVect();
        // - should only need 1 layer; also, z_phys_cc is created with ng=1 and
        //   looping over a larger box can cause the z_phys_cc array to take
        //   the wrong stride
        const amrex::IntVect ng(1,1,0);

        amrex::MultiFab min_thetav(ba,dm,1,ng);
        min_thetav.setVal(1.E34);

        amrex::MultiFab pblh_tke(ba,dm,1,ng);
        pblh_tke.setVal(0);

        pblh->setVal(0);

        // Now, loop over columns
        for (amrex::MFIter mfi(cons,TileNoZ()); mfi.isValid(); ++mfi)
        {
            amrex::Box gtbx = mfi.tilebox();
            gtbx.grow(ng); // grow in x,y only
            gtbx.setBig(2,gtbx.bigEnd(2)-1); // only loop to nz-1

            auto min_thv_arr = min_thetav.array(mfi);
            auto pblh_arr    = pblh->array(mfi);

            const auto cons_arr  = cons.const_array(mfi);
            const auto lmask_arr = (lmask) ? lmask->const_array(mfi) : amrex::Array4<int> {};

            // Find thetav-increase based PBL height (best for CBLs)
            // - WITH terrain/grid stretching
            // -----------------------------------------------------
            if (z_phys_cc)
            {
                const auto zphys_arr = z_phys_cc->const_array(mfi);

                // Find minimum thetav in the surface layer (this updates
                // ghost cells, too)
                ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if (zphys_arr(i,j,k) < thetamin_height) {
                        amrex::Real thv = Thetav(i,j,k,cons_arr,moist_flag);
                        if (min_thv_arr(i,j,0) > thv) min_thv_arr(i,j,0) = thv;
                    }
                });

                // This depends on TileNoZ and k increasing monotonically
                ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if (pblh_arr(i,j,0) == 0)
                    {
                        amrex::Real thv  = Thetav(i,j,k  ,cons_arr,moist_flag);
                        amrex::Real thv1 = Thetav(i,j,k+1,cons_arr,moist_flag);
                        int is_land = (lmask_arr) ? lmask_arr(i,j,0) : 1;
                        if (is_land && (thv1 >= min_thv_arr(i,j,0) + theta_incr_land)
                                    && (thv  <  min_thv_arr(i,j,0) + theta_incr_land))
                        {
                            pblh_arr(i,j,0) = zphys_arr(i,j,k)
                                            + (zphys_arr(i,j,k+1)-zphys_arr(i,j,k))/(thv1-thv)*(min_thv_arr(i,j,0) + theta_incr_land - thv);
                        }
                        else if (!is_land && (thv1 >= min_thv_arr(i,j,0) + theta_incr_water)
                                          && (thv  <  min_thv_arr(i,j,0) + theta_incr_water))
                        {
                            pblh_arr(i,j,0) = zphys_arr(i,j,k)
                                            + (zphys_arr(i,j,k+1)-zphys_arr(i,j,k))/(thv1-thv)*(min_thv_arr(i,j,0) + theta_incr_water - thv);
                        }
                    }
                });
            }
            // Find thetav-increase based PBL height (best for CBLs)
            // - NO terrain
            // -----------------------------------------------------
            else
            {
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
                        amrex::Real thv  = Thetav(i,j,k  ,cons_arr,moist_flag);
                        amrex::Real thv1 = Thetav(i,j,k+1,cons_arr,moist_flag);
                        int is_land = (lmask_arr) ? lmask_arr(i,j,0) : 1;
                        if (is_land && (thv1 >= min_thv_arr(i,j,0) + theta_incr_land)
                                    && (thv  <  min_thv_arr(i,j,0) + theta_incr_land))
                        {
                            pblh_arr(i,j,0) = (k+0.5)*dz_no_terrain + dz_no_terrain/(thv1-thv)*(min_thv_arr(i,j,0) + theta_incr_land - thv);
                        }
                        else if (!is_land && (thv1 >= min_thv_arr(i,j,0) + theta_incr_water)
                                          && (thv  <  min_thv_arr(i,j,0) + theta_incr_water))
                        {
                            pblh_arr(i,j,0) = (k+0.5)*dz_no_terrain + dz_no_terrain/(thv1-thv)*(min_thv_arr(i,j,0) + theta_incr_water - thv);
                        }
                    }
                });
            }

        }

        for (amrex::MFIter mfi(*pblh); mfi.isValid(); ++mfi)
        {
            const auto pblh_arr    = pblh->const_array(mfi);
            const auto min_thv_arr = min_thetav.const_array(mfi);
            const auto u_star_arr  = u_star->const_array(mfi);
            const auto t_star_arr  = t_star->const_array(mfi);
            const auto lmask_arr   = (lmask) ? lmask->const_array(mfi) : amrex::Array4<int> {};

            amrex::Box gtbx = mfi.growntilebox();
            ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                amrex::Real hfx = -u_star_arr(i,j,k) * t_star_arr(i,j,k);
                amrex::Print() << "Estimated PBLH" << amrex::IntVect(i,j,k) << " = " << pblh_arr(i,j,k)
                    << " with thetav_min, hfx = " << min_thv_arr(i,j,k) << " " << hfx << std::endl;
            });
        }//MFIter
    }

    amrex::Real thetamin_height  = 200.0; // [m]
    amrex::Real theta_incr_land  = 1.25;  // [K]
    amrex::Real theta_incr_water = 1.0;   // [K]
};
#endif
