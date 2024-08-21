#ifndef _PBL_HEIGHT_H_
#define _PBL_HEIGHT_H_

#include <AMReX_MultiFabUtil.H>
#include <Utils/Thetav.H>

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
        amrex::Print() << "moist_flag=" << moist_flag
            << " have_terrain=" << bool(z_phys_cc)
            << " nGrowVect()=" << cons.nGrowVect()
            << std::endl;
        // Find thetav_min (a BaseFab<Real>) in the lowest 200 m AGL
        int dir = 2; // z
        auto const& cons_arrs = cons.const_arrays();
        amrex::Box gtbx = domain; gtbx.grow(cons.nGrowVect());
        amrex::BaseFab<amrex::Real> thetav_min;
        if (z_phys_cc) {
            auto const& zphys_arrs = z_phys_cc->const_arrays();
            thetav_min = amrex::ReduceToPlane<amrex::ReduceOpMin,amrex::Real>(dir, gtbx, cons,
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) -> amrex::Real
                {
                    if (zphys_arrs[box_no](i,j,k) <= thetamin_height) {
                        return Thetav(i,j,k,cons_arrs[box_no],moist_flag);
                    } else {
                        return 9e9;
                    }
                });
        } else {
            // box size is known a priori
            int kmax = static_cast<int>(thetamin_height / dz_no_terrain);
            AMREX_ASSERT(kmax > 0);
            amrex::Box gtbxlow = gtbx;
            gtbxlow.setBig(2,kmax);
            thetav_min = amrex::ReduceToPlane<amrex::ReduceOpMin,amrex::Real>(dir, gtbxlow, cons,
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) -> amrex::Real
                {
                    return Thetav(i,j,k,cons_arrs[box_no],moist_flag);
                });
        }

        for (amrex::MFIter mfi(*pblh); mfi.isValid(); ++mfi)
        {
            amrex::Box gtbx = mfi.growntilebox();

            auto pblh_arr = pblh->array(mfi);

            const auto u_star_arr = u_star->const_array(mfi);
            const auto t_star_arr = t_star->const_array(mfi);
            const auto thv_min_arr = thetav_min.const_array();
            const auto lmask_arr  = (lmask) ? lmask->const_array(mfi) : amrex::Array4<int> {};
            const auto zphys_arr  = (z_phys_cc) ? z_phys_cc->const_array(mfi) : amrex::Array4<const amrex::Real>{};

            ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                amrex::Real hfx = -u_star_arr(i,j,k) * t_star_arr(i,j,k);
                amrex::Print() << "Computing PBLH at " << amrex::IntVect(i,j,k)
                    << " with thetav_min, hfx = " << thv_min_arr(i,j,k) << " " << hfx << std::endl;
                pblh_arr(i,j,k) = 0.0;
            });
        }//MFIter
    }

    amrex::Real thetamin_height = 200.0; // [m]
    amrex::Real theta_incr      = 1.5;   // [K]
};
#endif
