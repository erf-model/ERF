#ifndef ERF_PBL_HEIGHT_YSU_MRF_H_
#define ERF_PBL_HEIGHT_YSU_MRF_H_

#include <AMReX_MultiFabUtil.H>
#include <ERF_TileNoZ.H>
#include <ERF_Constants.H>
#include <ERF_TerrainMetrics.H>
#include <ERF_IndexDefines.H>

/*
 * Scheme-native PBL-height diagnostics for the YSU and MRF PBL schemes.
 *
 * Both schemes already compute a bulk-Richardson PBL height internally inside
 * their diffusivity routines (Source/PBL/ERF_ComputeDiffusivityYSU.cpp,
 * ERF_ComputeDiffusivityMRF.cpp) but discard it into a scratch FArrayBox after
 * the K-profile is built. These two estimator structs reproduce that exact
 * bulk-Ri logic and write it into the persistent `pblh` MultiFab so it can be
 * output as the `pblh` plot variable for runs using erf.pbl_type = YSU / MRF.
 *
 * Design notes
 * ------------
 *  - Unlike MYNNPBLH (which needs only cons + heights), the YSU/MRF bulk-Ri
 *    method needs the MOST surface-layer state (t_surf, the 10 m averages
 *    t10av/ws10av, z0, u_star, t_star, Obukhov length) plus the horizontal
 *    velocities for the wind-shear denominator. These are passed in explicitly
 *    by SurfaceLayer::update_pblh rather than discovered via a template.
 *  - This is a DIAGNOSTIC: it must never Abort. The diffusivity YSU routine
 *    aborts on unstable columns (it only builds a stable K-profile); here we
 *    instead always run the bulk-Ri search and fall back to the first cell
 *    height if the critical Rib is never reached.
 *  - Heights use z_phys_cc (cell-centered, same field MYNNPBLH uses). For the
 *    no-terrain case (z_phys_cc == nullptr) we use (k+1/2)*dz. Over the (near
 *    sea-level) coastal/ocean domains this matches the AGL heights the
 *    diffusivity routines use via Compute_Zrel to within the terrain height.
 */

namespace ERFPBLHUtil {

    // Cell-center height: prefer z_phys_cc, else uniform grid.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real zcc (int i, int j, int k,
                     const amrex::Array4<const amrex::Real>& zphys_arr,
                     bool have_terrain, amrex::Real dz, amrex::Real problo_z,
                     int imin, int jmin, int imax, int jmax) noexcept
    {
        if (have_terrain) {
            int ii = amrex::max(amrex::min(i,imax),imin);
            int jj = amrex::max(amrex::min(j,jmax),jmin);
            return zphys_arr(ii,jj,k);
        } else {
            return problo_z + (k + myhalf)*dz;
        }
    }

    // Cell-center wind speed squared from face-centered velocities.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real ws2_cc (int i, int j, int k,
                        const amrex::Array4<const amrex::Real>& u,
                        const amrex::Array4<const amrex::Real>& v) noexcept
    {
        return fourth*( (u(i,j,k)+u(i+1,j  ,k))*(u(i,j,k)+u(i+1,j  ,k))
                      + (v(i,j,k)+v(i  ,j+1,k))*(v(i,j,k)+v(i  ,j+1,k)) );
    }

} // namespace ERFPBLHUtil

/*
 * YSU PBL height (bulk-Richardson, stable-conditions form from
 * ERF_ComputeDiffusivityYSU.cpp). Reports the height where the bulk Ri reaches
 * the critical value (land vs. water Rib_cr); falls back to the first cell when
 * the column never becomes critical.
 */
struct YSUPBLH {

    amrex::Real m_land_Ribcr;
    amrex::Real m_coriolis_freq;
    bool        m_force_over_water;
    amrex::Real m_most_zref;

    YSUPBLH (amrex::Real land_Ribcr, amrex::Real coriolis_freq,
             bool force_over_water, amrex::Real most_zref)
        : m_land_Ribcr(land_Ribcr), m_coriolis_freq(coriolis_freq),
          m_force_over_water(force_over_water), m_most_zref(most_zref) {}

    AMREX_GPU_HOST
    AMREX_FORCE_INLINE
    void compute_pblh (const amrex::Geometry& geom,
                       const amrex::MultiFab* z_phys_cc,
                       amrex::MultiFab* pblh,
                       const amrex::MultiFab& cons,
                       const amrex::MultiFab& xvel,
                       const amrex::MultiFab& yvel,
                       const amrex::MultiFab& t_surf,
                       const amrex::MultiFab& t10av,
                       const amrex::MultiFab& ws10av,
                       const amrex::MultiFab& z0,
                       const amrex::iMultiFab* lmask) const
    {
        using namespace amrex;
        const Real f0           = m_coriolis_freq;
        const bool fow          = m_force_over_water;
        const Real land_Ribcr   = m_land_Ribcr;
        const Real most_zref    = m_most_zref;
        const Box& domain       = geom.Domain();
        const int  klo          = domain.smallEnd(2);
        const int  khi          = domain.bigEnd(2);
        const Real dz           = geom.CellSize(2);
        const Real problo_z     = geom.ProbLo(2);
        const bool have_terrain = (z_phys_cc != nullptr);

        pblh->setVal(zero);

        for (MFIter mfi(cons); mfi.isValid(); ++mfi)
        {
            Box vbx = mfi.validbox();
            Box xybx(IntVect(vbx.smallEnd(0),vbx.smallEnd(1),0),
                     IntVect(vbx.bigEnd(0)  ,vbx.bigEnd(1)  ,0));

            const auto cons_arr = cons.const_array(mfi);
            const auto uvel     = xvel.const_array(mfi);
            const auto vvel     = yvel.const_array(mfi);
            auto       pblh_arr = pblh->array(mfi);

            const auto t_surf_arr = t_surf.const_array(mfi);
            const auto t10av_arr  = t10av.const_array(mfi);
            const auto ws10av_arr = ws10av.const_array(mfi);
            const auto z0_arr     = z0.const_array(mfi);
            const auto lmask_arr  = (lmask) ? lmask->const_array(mfi) : Array4<int>{};

            Array4<const Real> zphys_arr{};
            int imin=0,jmin=0,imax=0,jmax=0;
            if (have_terrain) {
                zphys_arr = z_phys_cc->const_array(mfi);
                imin = lbound(zphys_arr).x; jmin = lbound(zphys_arr).y;
                imax = ubound(zphys_arr).x; jmax = ubound(zphys_arr).y;
            }

            ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                // Reconstruct the surface bulk Richardson number from the MOST
                // surface-layer state (in WRF this is supplied to YSU by the
                // surface-layer model).
                const Real t_surf  = t_surf_arr(i,j,0);
                const Real t_layer = t10av_arr(i,j,0);
                const Real ws_layer= ws10av_arr(i,j,0);
                const Real ws2_ref = amrex::max(ws_layer*ws_layer, Real(1.0e-8));
                const Real Rib_layer = CONST_GRAV * most_zref / ws2_ref
                                     * (t_layer - t_surf) / t_layer;

                // Critical bulk Ri: land constant, water Rossby-number based (H10).
                Real Rib_cr;
                bool over_land = (lmask_arr) ? lmask_arr(i,j,0) : 1;
                if (over_land && !fow) {
                    Rib_cr = land_Ribcr;
                } else {
                    const Real z0v = amrex::max(z0_arr(i,j,0), Real(1.0e-6));
                    const Real Rossby = ws_layer/(amrex::max(f0,Real(1.0e-12))*z0v);
                    Rib_cr = amrex::min(Real(0.16)*std::pow(Real(1.0e-7)*Rossby,-Real(0.18)),
                                        Real(0.3));
                }

                // March up the column until the bulk Ri exceeds the critical
                // value, interpolating to the crossing height.
                bool above_critical = false;
                int  kpbl = klo;
                Real Rib_up = Rib_layer, Rib_dn = Rib_layer;
                const Real base_theta = cons_arr(i,j,klo,RhoTheta_comp)
                                      / cons_arr(i,j,klo,Rho_comp);
                while (!above_critical && (kpbl+1) <= khi) {
                    kpbl += 1;
                    const Real zval = ERFPBLHUtil::zcc(i,j,kpbl,zphys_arr,have_terrain,
                                                       dz,problo_z,imin,jmin,imax,jmax);
                    const Real ws2  = amrex::max(ERFPBLHUtil::ws2_cc(i,j,kpbl,uvel,vvel),
                                                 Real(1.0e-8));
                    const Real theta= cons_arr(i,j,kpbl,RhoTheta_comp)
                                    / cons_arr(i,j,kpbl,Rho_comp);
                    Rib_dn = Rib_up;
                    Rib_up = (theta-base_theta)/base_theta * CONST_GRAV * zval / ws2;
                    above_critical = (Rib_up >= Rib_cr);
                }

                Real interp_fact;
                if (Rib_dn >= Rib_cr)      { interp_fact = zero; }
                else if (Rib_up <= Rib_cr) { interp_fact = one;  }
                else { interp_fact = (Rib_cr - Rib_dn)/(Rib_up - Rib_dn); }

                const Real zval_up = ERFPBLHUtil::zcc(i,j,kpbl  ,zphys_arr,have_terrain,
                                                      dz,problo_z,imin,jmin,imax,jmax);
                const Real zval_dn = ERFPBLHUtil::zcc(i,j,kpbl-1,zphys_arr,have_terrain,
                                                      dz,problo_z,imin,jmin,imax,jmax);
                Real h = zval_dn + interp_fact*(zval_up - zval_dn);

                // Diagnostic fallback: never let PBLH collapse below the first
                // cell (the diffusivity routine aborts here; we degrade instead).
                const Real zval_0 = ERFPBLHUtil::zcc(i,j,klo,zphys_arr,have_terrain,
                                                     dz,problo_z,imin,jmin,imax,jmax);
                pblh_arr(i,j,0) = amrex::max(h, zval_0);
            });
        }
    }
};

/*
 * MRF PBL height (Hong & Pan 1996), predictor + thermal-excess corrector,
 * ported from ERF_ComputeDiffusivityMRF.cpp. Reports the corrected height.
 */
struct MRFPBLH {

    amrex::Real m_Ribcr;
    amrex::Real m_const_b;
    amrex::Real m_sf;

    MRFPBLH (amrex::Real Ribcr, amrex::Real const_b, amrex::Real sf)
        : m_Ribcr(Ribcr), m_const_b(const_b), m_sf(sf) {}

    AMREX_GPU_HOST
    AMREX_FORCE_INLINE
    void compute_pblh (const amrex::Geometry& geom,
                       const amrex::MultiFab* z_phys_cc,
                       amrex::MultiFab* pblh,
                       const amrex::MultiFab& cons,
                       const amrex::MultiFab& xvel,
                       const amrex::MultiFab& yvel,
                       const amrex::MultiFab& u_star,
                       const amrex::MultiFab& t_star,
                       const amrex::MultiFab& olen,
                       const amrex::MultiFab& t10av) const
    {
        using namespace amrex;
        const Real Ribcr       = m_Ribcr;
        const Real const_b     = m_const_b;
        const Real sf          = m_sf;
        const Box& domain      = geom.Domain();
        const int  klo         = domain.smallEnd(2);
        const int  khi         = domain.bigEnd(2);
        const Real dz          = geom.CellSize(2);
        const Real problo_z    = geom.ProbLo(2);
        const bool have_terrain= (z_phys_cc != nullptr);

        pblh->setVal(zero);

        for (MFIter mfi(cons); mfi.isValid(); ++mfi)
        {
            Box vbx = mfi.validbox();
            Box xybx(IntVect(vbx.smallEnd(0),vbx.smallEnd(1),0),
                     IntVect(vbx.bigEnd(0)  ,vbx.bigEnd(1)  ,0));

            const auto cons_arr = cons.const_array(mfi);
            const auto uvel     = xvel.const_array(mfi);
            const auto vvel     = yvel.const_array(mfi);
            auto       pblh_arr = pblh->array(mfi);

            const auto u_star_arr = u_star.const_array(mfi);
            const auto t_star_arr = t_star.const_array(mfi);
            const auto l_obuk_arr = olen.const_array(mfi);
            const auto t10av_arr  = t10av.const_array(mfi);

            Array4<const Real> zphys_arr{};
            int imin=0,jmin=0,imax=0,jmax=0;
            if (have_terrain) {
                zphys_arr = z_phys_cc->const_array(mfi);
                imin = lbound(zphys_arr).x; jmin = lbound(zphys_arr).y;
                imax = ubound(zphys_arr).x; jmax = ubound(zphys_arr).y;
            }

            ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                const Real t_layer = t10av_arr(i,j,0);

                // ---- Step 1: predictor PBL height (no thermal excess) ----
                Real zval_p = problo_z + (klo+myhalf)*dz;
                int  kpbl_p = klo;
                bool crit_p = false;
                while (!crit_p && (kpbl_p+1) <= khi) {
                    kpbl_p += 1;
                    zval_p = ERFPBLHUtil::zcc(i,j,kpbl_p,zphys_arr,have_terrain,
                                              dz,problo_z,imin,jmin,imax,jmax);
                    const Real theta = cons_arr(i,j,kpbl_p,RhoTheta_comp)
                                     / cons_arr(i,j,kpbl_p,Rho_comp);
                    const Real ws2   = amrex::max(ERFPBLHUtil::ws2_cc(i,j,kpbl_p,uvel,vvel),
                                                  Real(1.0e-8));
                    const Real Rib   = CONST_GRAV * zval_p * (theta - t_layer)/(ws2*t_layer);
                    crit_p = (Rib >= Ribcr);
                }
                const Real zval0_cell = ERFPBLHUtil::zcc(i,j,klo,zphys_arr,have_terrain,
                                                         dz,problo_z,imin,jmin,imax,jmax);
                const Real pblh_pred = (crit_p) ? zval_p : zval0_cell;

                // ---- Step 2: thermal-excess corrector ----
                const Real lobuk = l_obuk_arr(i,j,0);
                const Real phiM  = (lobuk > 0)
                                 ? (one + Real(5.0)*sf*pblh_pred/lobuk)
                                 : std::pow(one - Real(8.0)*sf*pblh_pred/lobuk, -one/three);
                const Real wstar = u_star_arr(i,j,0)/phiM;
                const Real t_excess = -const_b * u_star_arr(i,j,0) * t_star_arr(i,j,0)
                                    / amrex::max(wstar, Real(1.0e-8));
                const Real t_surf = t_layer
                                  + amrex::max(amrex::min(t_excess, Real(3.0)), Real(0.0));

                int  kpbl = klo;
                Real zval0, zval, Rib0, Rib;
                {
                    zval = ERFPBLHUtil::zcc(i,j,kpbl,zphys_arr,have_terrain,
                                            dz,problo_z,imin,jmin,imax,jmax);
                    const Real theta = cons_arr(i,j,kpbl,RhoTheta_comp)
                                     / cons_arr(i,j,kpbl,Rho_comp);
                    const Real ws2   = amrex::max(ERFPBLHUtil::ws2_cc(i,j,kpbl,uvel,vvel),
                                                  Real(1.0e-8));
                    Rib = CONST_GRAV * zval * (theta - t_surf)/(ws2*t_layer);
                }
                bool above_critical = false;
                while (!above_critical && (kpbl+1) <= khi) {
                    zval0 = zval; Rib0 = Rib; kpbl += 1;
                    zval = ERFPBLHUtil::zcc(i,j,kpbl,zphys_arr,have_terrain,
                                            dz,problo_z,imin,jmin,imax,jmax);
                    const Real theta = cons_arr(i,j,kpbl,RhoTheta_comp)
                                     / cons_arr(i,j,kpbl,Rho_comp);
                    const Real ws2   = amrex::max(ERFPBLHUtil::ws2_cc(i,j,kpbl,uvel,vvel),
                                                  Real(1.0e-8));
                    Rib = CONST_GRAV * zval * (theta - t_surf)/(ws2*t_layer);
                    above_critical = (Rib >= Ribcr);
                }

                if (above_critical) {
                    pblh_arr(i,j,0) = zval0 + (zval-zval0)/(Rib-Rib0)*(Ribcr-Rib0);
                } else {
                    pblh_arr(i,j,0) = zval0_cell;
                }
            });
        }
    }
};

#endif
