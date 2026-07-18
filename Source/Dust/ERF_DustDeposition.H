/**
 * @file ERF_DustDeposition.H
 * @brief Dry deposition lower boundary condition for the dust scalar.
 *
 * Replaces the zero lower-boundary flux in ERF_DustSettling.H with the
 * dry deposition velocity v_d from Seinfeld & Pandis (2006) Chapter 19.
 *
 * v_d = v_s + 1 / (r_a + r_s + r_a * r_s * v_s)
 *   r_a = 1 / (u* * kappa)    aerodynamic resistance [s/m]
 *   r_s = 1 / (u* * E_0)      surface resistance [s/m]
 *   kappa = 0.4 (von Karman)
 *   E_0 = 3e-3 (bare mine surface, Zhang et al. 2001)
 *
 * Structural reference: ERF_FireAtmCoupling.H in Source/Fire/:
 *   https://github.com/hgopalan/ERF/blob/ecb9b7656bbb98fcc5d4f781e486825ee9478a73/
 *   Source/Fire/ERF_FireAtmCoupling.H
 *
 * References:
 *   Seinfeld, J.H., & Pandis, S.N. (2006). Atmospheric Chemistry and Physics,
 *     2nd ed. Wiley, Ch. 19. ISBN 978-0-471-72018-8.
 *   Zhang, L. et al. (2001). Atmos. Environ., 35, 549-560.
 *     https://doi.org/10.1016/S1352-2310(00)00326-5
 *   Slinn, W.G.N. (1982). Atmos. Environ., 16, 1785-1794.
 *     https://doi.org/10.1016/0004-6981(82)90271-2
 */

#ifndef ERF_DUST_DEPOSITION_H
#define ERF_DUST_DEPOSITION_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>
#include <AMReX_Print.H>
#include <cmath>

// Dust settling constants and functions
#include <ERF_DustSettling.H>

// Density component index (needed for S_old array access)
#ifndef Rho_comp
#include <ERF_IndexDefines.H>
#endif

namespace DustDepositionConst {
    constexpr amrex::Real KAPPA   = 0.4;      ///< von Karman constant [-]
    constexpr amrex::Real E0_BARE = 3.0e-3;   ///< Collection efficiency for bare surface [-]
    constexpr amrex::Real R_MAX   = 1.0e6;    ///< Maximum resistance [s/m]
    constexpr amrex::Real VD_MAX  = 0.1;      ///< Maximum deposition velocity [m/s]
    constexpr amrex::Real VD_MIN  = 0.0;      ///< Minimum deposition velocity [m/s]
}

/**
 * @brief Aerodynamic resistance r_a [s/m].
 *
 * r_a = 1 / (u* * kappa)
 *
 * @param u_star Friction velocity [m/s].
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_r_a(amrex::Real u_star) noexcept {
    return amrex::min(1.0 / (amrex::max(u_star, 1.0e-4) * DustDepositionConst::KAPPA),
                      DustDepositionConst::R_MAX);
}

/**
 * @brief Surface resistance r_s [s/m].
 *
 * r_s = 1 / (u* * E_0)
 *
 * @param u_star Friction velocity [m/s].
 * @param E_0    Collection efficiency [-].
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_r_s(amrex::Real u_star, amrex::Real E_0) noexcept {
    return amrex::min(1.0 / (amrex::max(u_star, 1.0e-4) * amrex::max(E_0, 1.0e-6)),
                      DustDepositionConst::R_MAX);
}

/**
 * @brief Dry deposition velocity v_d [m/s].
 *
 * v_d = v_s + 1 / (r_a + r_s + r_a * r_s * v_s)
 *
 * @param v_s    Stokes settling velocity [m/s] from ERF_DustSettling.H.
 * @param u_star Friction velocity [m/s] from dust_ustar_in (Phase 9).
 * @param E_0    Collection efficiency. Default DustDepositionConst::E0_BARE.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_deposition_velocity(amrex::Real v_s, amrex::Real u_star,
                                         amrex::Real E_0) noexcept {
    amrex::Real ra    = compute_r_a(u_star);
    amrex::Real rs    = compute_r_s(u_star, E_0);
    amrex::Real denom = ra + rs + ra * rs * v_s;
    amrex::Real vd    = (denom > 1.0e-10) ? v_s + 1.0/denom : v_s;
    return amrex::min(amrex::max(vd, DustDepositionConst::VD_MIN),
                      DustDepositionConst::VD_MAX);
}

/**
 * @brief Apply dry deposition tendency at k=klo and write deposition flux.
 *
 * At k=klo only:
 *   d(dust_comp)/dt += -v_d * RhoDust(k=klo) / dz_k0
 *   dep_flux_out(i,j,0) = v_d * RhoDust(k=klo)   [kg/m^2/s]
 *
 * dep_flux_out is later mapped to the 2D dust grid by
 * apply_deposition_to_dust_grid. All k > klo: no change from this function.
 *
 * u_star is taken from ustar_atm (a 2D MultiFab at k=0 on the atm grid,
 * filled from SurfaceLayer::get_u_star(0) by fill_dust_ustar_from_surface_layer
 * in Phase 9). When grid_ratio=1, ustar_atm is on the same 2D BoxArray as
 * the dust grid and can be used directly.
 */
inline void apply_dust_deposition_bc(
    amrex::MultiFab&       cc_source,
    amrex::MultiFab&       dep_flux_out,
    const amrex::MultiFab& S_old,
    const amrex::MultiFab& ustar_atm,
    const amrex::MultiFab& z_phys_cc,
    const amrex::Geometry& geom_atm,
    amrex::Real d_m, amrex::Real rho_p, amrex::Real E_0,
    int dust_comp, bool dust_debug)
{
    using namespace amrex;
    const Box& domain = geom_atm.Domain();
    const int klo = domain.smallEnd(2);
    const int khi = domain.bigEnd(2);
    const Real dz_avg = geom_atm.CellSize(2);

    dep_flux_out.setVal(0.0);

    for (MFIter mfi(cc_source, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        if (bx.smallEnd(2) > klo) continue;

        auto src = cc_source.array(mfi);
        auto dep = dep_flux_out.array(mfi);
        auto s   = S_old.const_array(mfi);
        auto us  = ustar_atm.const_array(mfi);
        auto z   = z_phys_cc.const_array(mfi);

        Box bx_sfc = bx; bx_sfc.setSmall(2,klo); bx_sfc.setBig(2,klo);

        const int   comp  = dust_comp;
        const Real  dp    = d_m;
        const Real  rhop  = rho_p;
        const Real  e0    = E_0;

        ParallelFor(bx_sfc, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            Real rho_a = s(i,j,k,Rho_comp);
            if (rho_a <= 0.0) rho_a = 1.2;

            Real dz = (k < khi) ? (z(i,j,k+1)-z(i,j,k)) : dz_avg;
            if (dz <= 1.0e-10) dz = dz_avg;

            Real u_star = amrex::max(us(i,j,0), 1.0e-4);
            Real vs = compute_stokes_settling(dp, rhop, rho_a,
                                              DustSettlingConst::MU_AIR_STD);
            Real vd = compute_deposition_velocity(vs, u_star, e0);

            Real rho_dust = amrex::max(src(i,j,k,comp), 0.0);
            dep(i,j,0)    = vd * rho_dust;
            src(i,j,k,comp) -= vd * rho_dust / dz;
        });
    }

    if (dust_debug) {
        Real us_max  = ustar_atm.max(0);
        Real vs_ref  = compute_stokes_settling(d_m, rho_p, 1.225,
                                               DustSettlingConst::MU_AIR_STD);
        Real vd_ref  = compute_deposition_velocity(vs_ref, us_max, E_0);
        Real dep_max = dep_flux_out.max(0);
        amrex::Print() << "[DUST DEPOSITION] d=" << d_m*1e6
                       << " um  u*_max=" << us_max
                       << " m/s  v_d=" << vd_ref
                       << " m/s  dep_flux_max=" << dep_max << " kg/m^2/s\n";
    }
}

/**
 * @brief Accumulate deposition flux from atm grid onto 2D dust grid.
 *
 * Adds dep_flux_atm(i/C, j/C, 0) * dt to dust_deposition_rate(i,j,0)
 * each timestep. dust_deposition_rate is initialised to zero and
 * accumulates without reset (total deposited mass [kg/m^2] over the run).
 */
inline void apply_deposition_to_dust_grid(
    amrex::MultiFab&       dust_deposition_rate,
    const amrex::MultiFab& dep_flux_atm,
    int C, amrex::Real dt)
{
    for (amrex::MFIter mfi(dust_deposition_rate, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto dep_d = dust_deposition_rate.array(mfi);
        auto dep_a = dep_flux_atm.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            dep_d(i,j,k) += dep_a(i/C, j/C, 0) * dt;
        });
    }
}

#endif // ERF_DUST_DEPOSITION_H
