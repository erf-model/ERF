/**
 * @file ERF_DustSettling.H
 * @brief Stokes settling velocity and vertical settling tendency for dust scalars.
 *
 * Computes per-bin Stokes settling velocity and applies the resulting
 * downward mass flux divergence as a tendency in cc_source.
 *
 * Settling velocity (Stokes law, valid for d < 100 um):
 *   v_s = (rho_p - rho_a) * g * d^2 / (18 * mu_a)
 *
 * Cunningham slip correction (Stokes-Cunningham):
 *   C_c = 1 + (2 * lambda / d) * (A1 + A2 * exp(-A3 * d / lambda))
 *   v_s_corrected = v_s * C_c
 * where lambda = 0.066 um (mean free path at STP),
 * A1=1.257, A2=0.400, A3=0.55 (Allen & Raabe 1985).
 * For d > 1 um the correction is < 15% and may be omitted for coarse bins.
 *
 * Vertical tendency (first-order upwind flux divergence):
 *   d(RhoDust)/dt += -(v_s(k)*RhoDust(k) - v_s(k-1)*RhoDust(k-1)) / dz
 * Downward flux at k=0 is held at zero until Phase 12 dry deposition.
 *
 * Pattern reference: ERF_FireAtmCoupling.H in Source/Fire/:
 *   https://github.com/hgopalan/ERF/blob/7856affedd8bb3c878b62c308f39c36bf2228791/
 *   Source/Fire/ERF_FireAtmCoupling.H
 *
 * References:
 *   Stokes, G.G. (1851). Trans. Cambridge Phil. Soc., 9, 8-106.
 *   Allen, M.D., & Raabe, O.G. (1985). Aerosol Sci. Tech., 4, 269-286.
 *     https://doi.org/10.1080/02786828508959055
 *   Seinfeld, J.H., & Pandis, S.N. (2006). Atmospheric Chemistry and Physics,
 *     2nd ed. Wiley. Chapter 9. ISBN 978-0-471-72018-8.
 *   Ginoux, P. et al. (2001). J. Geophys. Res., 106, 20255.
 *     https://doi.org/10.1029/2000JD901323
 */

#ifndef ERF_DUST_SETTLING_H
#define ERF_DUST_SETTLING_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>
#include <AMReX_Print.H>
#include <cmath>

// Density component index (needed for S_old array access)
#ifndef Rho_comp
#include <ERF_IndexDefines.H>
#endif

namespace DustSettlingConst {
    constexpr amrex::Real G_ACCEL    = 9.81;     ///< Gravitational acceleration [m/s^2]
    constexpr amrex::Real MU_AIR_STD = 1.81e-5;  ///< Dynamic viscosity of air at 293 K [Pa s]
    constexpr amrex::Real LAMBDA_MFP = 6.6e-8;   ///< Mean free path at STP [m] (0.066 um)
    constexpr amrex::Real CC_A1      = 1.257;    ///< Cunningham coefficient A1
    constexpr amrex::Real CC_A2      = 0.400;    ///< Cunningham coefficient A2
    constexpr amrex::Real CC_A3      = 0.550;    ///< Cunningham coefficient A3
    constexpr amrex::Real V_S_MAX    = 1.0;      ///< Settling velocity cap [m/s]
    constexpr amrex::Real V_S_MIN    = 0.0;      ///< Settling velocity floor [m/s]
}

/**
 * @brief Cunningham slip correction factor for sub-micron particles.
 *
 * C_c = 1 + (2*lambda/d) * (A1 + A2*exp(-A3*d/lambda))
 * For d >> lambda, C_c approaches 1.
 * Reference: Allen & Raabe (1985), https://doi.org/10.1080/02786828508959055
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_cunningham_cc(amrex::Real d_m) noexcept
{
    if (d_m <= 0.0) return 1.0;
    amrex::Real ratio = 2.0 * DustSettlingConst::LAMBDA_MFP / d_m;
    return 1.0 + ratio * (DustSettlingConst::CC_A1
                        + DustSettlingConst::CC_A2
                          * std::exp(-DustSettlingConst::CC_A3 * d_m
                                     / DustSettlingConst::LAMBDA_MFP));
}

/**
 * @brief Stokes settling velocity with Cunningham correction [m/s].
 *
 * v_s = (rho_p - rho_a) * g * d^2 * C_c / (18 * mu_a)
 * Clamped to [V_S_MIN, V_S_MAX].
 *
 * @param d_m    Particle diameter [m].
 * @param rho_p  Particle density [kg/m^3].
 * @param rho_a  Local air density [kg/m^3].
 * @param mu_a   Dynamic viscosity [Pa s]. Uses MU_AIR_STD if <= 0.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_stokes_settling(amrex::Real d_m, amrex::Real rho_p,
                                     amrex::Real rho_a, amrex::Real mu_a) noexcept
{
    if (d_m <= 0.0 || rho_p <= 0.0) return 0.0;
    amrex::Real mu = (mu_a > 0.0) ? mu_a : DustSettlingConst::MU_AIR_STD;
    amrex::Real rho_diff = rho_p - amrex::max(rho_a, 0.0);
    if (rho_diff <= 0.0) return 0.0;
    amrex::Real Cc  = compute_cunningham_cc(d_m);
    amrex::Real vs  = rho_diff * DustSettlingConst::G_ACCEL * d_m * d_m * Cc / (18.0 * mu);
    return amrex::min(amrex::max(vs, DustSettlingConst::V_S_MIN),
                      DustSettlingConst::V_S_MAX);
}

/**
 * @brief Apply Stokes settling tendency to cc_source for one dust bin.
 *
 * First-order upwind flux divergence:
 *   flux_down(k)   = v_s * RhoDust(k)        (downward flux at top of cell k)
 *   flux_down(k-1) = v_s * RhoDust(k-1)      (downward flux at bottom of cell k)
 *   tend(k) = -(flux_down(k) - flux_down(k-1)) / dz
 *
 * At k = klo: flux_down(k-1) = 0 (no deposition; Phase 12 adds this).
 * At k = khi: uses cell values without extrapolation.
 *
 * v_s is computed per-cell from local air density S_old(Rho_comp).
 * Pattern: follows apply_fire_tendency_to_cc_source in ERF_FireAtmCoupling.H.
 *
 * @param[in,out] cc_source     Cell-centred source MultiFab on atm grid.
 * @param[in]     S_old         Conserved state at time n (for air density).
 * @param[in]     z_phys_cc     Cell-centre terrain height [m].
 * @param[in]     geom_atm      Atmospheric geometry.
 * @param[in]     d_m           Particle diameter for this bin [m].
 * @param[in]     rho_p         Particle density [kg/m^3].
 * @param[in]     dust_comp     cc_source component index for this bin.
 * @param[in]     dust_debug    Print diagnostics when true.
 */
inline void apply_dust_settling_to_cc_source(
    amrex::MultiFab&       cc_source,
    const amrex::MultiFab& S_old,
    const amrex::MultiFab& z_phys_cc,
    const amrex::Geometry& geom_atm,
    amrex::Real            d_m,
    amrex::Real            rho_p,
    int                    dust_comp,
    bool                   dust_debug = false)
{
    using namespace amrex;
    const Box& domain = geom_atm.Domain();
    int klo = domain.smallEnd(2);
    int khi = domain.bigEnd(2);
    const auto& dx = geom_atm.CellSize();
    Real dz_avg = dx[2];

    for (MFIter mfi(cc_source, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto src_arr = cc_source.array(mfi);
        auto s_arr   = S_old.const_array(mfi);
        auto z_arr   = z_phys_cc.const_array(mfi);
        const int comp  = dust_comp;
        const Real dp   = d_m;
        const Real rhop = rho_p;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real rho_a = s_arr(i, j, k, Rho_comp);
            if (rho_a <= 0.0) rho_a = 1.2;

            Real dz = (k < khi) ? (z_arr(i,j,k+1) - z_arr(i,j,k)) : dz_avg;
            if (dz <= 1.0e-10) dz = dz_avg;

            // Settling velocity at this level using local air density.
            Real vs_k = compute_stokes_settling(dp, rhop, rho_a,
                                                DustSettlingConst::MU_AIR_STD);

            // Downward flux at top face of cell k.
            Real flux_k = vs_k * src_arr(i, j, k, comp);

            // Downward flux at bottom face of cell k (from cell k-1).
            Real flux_km1 = 0.0;
            if (k > klo) {
                Real rho_a_km1 = s_arr(i, j, k-1, Rho_comp);
                if (rho_a_km1 <= 0.0) rho_a_km1 = 1.2;
                Real vs_km1 = compute_stokes_settling(dp, rhop, rho_a_km1,
                                                      DustSettlingConst::MU_AIR_STD);
                flux_km1 = vs_km1 * src_arr(i, j, k-1, comp);
            }
            // At k=klo: flux_km1 = 0 (Phase 12 will add deposition here).

            // First-order upwind: tend = -(flux_k - flux_km1) / dz
            src_arr(i, j, k, comp) -= (flux_k - flux_km1) / dz;
        });
    }

    if (dust_debug) {
        Real vs_std = compute_stokes_settling(d_m, rho_p, 1.225,
                                              DustSettlingConst::MU_AIR_STD);
        Real tend_max = cc_source.max(dust_comp);
        amrex::Print() << "[DUST SETTLING] bin d=" << d_m*1e6
                       << " um  v_s=" << vs_std << " m/s"
                       << "  tend_max=" << tend_max << " [kg/m^3/s]\n";
    }
}

#endif // ERF_DUST_SETTLING_H
