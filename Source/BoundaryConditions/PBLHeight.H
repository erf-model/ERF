#ifndef _PBL_HEIGHT_H_
#define _PBL_HEIGHT_H_

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
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real compute_pblh(int i, int j, int k,
                             const amrex::Real dz,
                             const amrex::Array4<const amrex::Real>& zphys_arr,
                             const amrex::Array4<const amrex::Real>& cons_arr,
                             const amrex::Array4<const amrex::Real>& u_star_arr,
                             const amrex::Array4<const amrex::Real>& t_star_arr,
                             const amrex::Array4<const int        >& lmask_arr,
                             const int moist_flag) const
    {
        amrex::Real thv = Thetav(i,j,k, cons_arr, moist_flag);
        amrex::Real hfx = -u_star_arr(i,j,k) * t_star_arr(i,j,k);
        amrex::Print() << "Computing PBLH at " << amrex::IntVect(i,j,k)
            << " with thetav, hfx = " << thv << " " << hfx << std::endl;
        return 0;
    }

    amrex::Real theta_incr = 1.5; // [K]
};
#endif
