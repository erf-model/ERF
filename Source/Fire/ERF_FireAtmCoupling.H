// One-way fire-to-atmosphere coupling: distribute fire heat flux through
// the atmospheric column and inject into cc_source.
//
// This mirrors the WRF-Fire approach in module_fr_sfire_atm.F
// (subroutine fire_tendency), which:
//   1. Builds a vertical flux profile with exponential decay:
//        hfx(z) = (Q_sfc / Cp) * exp(-z / alfg)
//      where z is height above terrain, alfg is the e-folding scale [m].
//   2. Converts the flux profile to a theta tendency via flux divergence:
//        d(RhoTheta)/dt = -rho * d(hfx)/dz
//   3. Does the same for latent heat (grnqfx -> rqvfrten):
//        qfx(z) = (Q_lat / Lv) * exp(-z / alfg)
//        d(RhoQ1)/dt   = -rho * d(qfx)/dz
//
// The e-folding scale alfg controls how deep into the atmosphere the fire
// heating penetrates. With alfg = 45 m, the flux has decayed to 37% at 45 m,
// 5% at 135 m, and 1% at 225 m above the surface.
//
// Reference:
//   Mandel, J. et al. (2011). A wildland fire model with data assimilation.
//   Mathematics and Computers in Simulation, 79(3), 584-606.
//   WRF-SFIRE source: phys/module_fr_sfire_atm.F, subroutine fire_tendency.
//
// Coupling time lag:
//   ERF uses one-step explicit lag: fire heat flux computed at time n is
//   injected into cc_source at the beginning of step n+1, before advance_dycore().
//   This is consistent with WRF-Fire's explicit coupling approach.
//
// Latent heat:
//   WRF-Fire computes grnqfx (ground latent heat flux) alongside grnhfx and
//   passes both through fire_tendency. ERF-Fire follows the same approach:
//   Q_latent = Lv * bmst * (dmass/dt)   where bmst = M_f / (1 + M_f)
//   This is stored as fire_latent_flux [W/m2] and distributed with the same
//   exponential profile.
//
// Assumptions:
//   - Single AMR level (no refinement in z between fire and atmospheric grids).
//   - fire_heat_flux and fire_latent_flux are on the fire grid (C*nx x C*ny x 1).
//   - The coarsening from fire grid to atmospheric grid uses average_down().
//   - z_phys_cc is the cell-centre terrain-following height MultiFab [m].
//   - The atmospheric BoxArray is not split in z (max_grid_size_z = nz enforced).
//   - Latent heat injection into RhoQ1_comp is skipped when moisture is inactive.
//   - fire_atm_feedback multiplier (equivalent to WRF-Fire's fire_atm_feedback)
//     scales both sensible and latent flux. Set to 0 for one-way coupling only
//     (fire responds to wind but atmosphere does not respond to fire heat).

#ifndef ERF_FIRE_ATM_COUPLING_H
#define ERF_FIRE_ATM_COUPLING_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_Box.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFabUtil.H>
#include <cmath>
#include <ERF_Constants.H>
#include <ERF_IndexDefines.H>

// Default e-folding height for ground fire heat distribution [m].
// WRF-Fire default: alfg ~ 45-75 m. At z = alfg: flux decays to 1/e ~ 37%.
constexpr amrex::Real FIRE_ATM_ALFG_DEFAULT = 45.0;

// Default e-folding height for canopy fire contribution [m].
// Only used when canopy fire heat is non-zero.
constexpr amrex::Real FIRE_ATM_ALFC_DEFAULT = 60.0;

// Coarsen fire_heat_flux from the fire grid (C*nx x C*ny x 1) to the
// atmospheric level-0 grid (nx x ny x 1) using area-weighted averaging.
// Uses amrex::average_down() with refinement ratio (C, C, 1).
inline void coarsen_fire_flux_to_atm(
    amrex::MultiFab& Q_atm_out,
    const amrex::MultiFab& Q_fire,
    const amrex::Geometry& geom_atm,
    const amrex::Geometry& geom_fire,
    int C)
{
    amrex::average_down(Q_fire, Q_atm_out, geom_fire, geom_atm, 0, Q_fire.nComp(),
                        amrex::IntVect(C, C, 1));
}

// Compute fire latent heat flux [W/m2] from sensible heat flux and fuel moisture.
inline void compute_fire_latent_flux(
    amrex::MultiFab& Q_lat_out,
    const amrex::MultiFab& Q_sfc,
    amrex::Real M_f,
    amrex::Real h_fuel_Jkg)
{
    if (h_fuel_Jkg <= 0.0) {
        Q_lat_out.setVal(0.0);
        return;
    }

    constexpr amrex::Real L_V = 2.501e6;
    amrex::Real bmst = M_f / (1.0 + M_f);
    constexpr amrex::Real f_water = 0.56;
    amrex::Real ratio = L_V / h_fuel_Jkg;

    for (amrex::MFIter mfi(Q_sfc, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const& q_sfc_a = Q_sfc.array(mfi);
        auto const& q_lat_a = Q_lat_out.array(mfi);
        auto const& bx = mfi.tilebox();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(const amrex::IntVect& iv) {
            amrex::Real factor = ratio * (bmst + (1.0 - bmst) * f_water);
            q_lat_a(iv) = factor * q_sfc_a(iv);
        });
    }
}

// Distribute fire surface heat flux through the atmospheric column and add
// the resulting tendency to cc_source.
//
// Optional fire_debug argument prints per-call tendency diagnostics to stdout:
//   [FIRE COUPLING] RhoTheta tendency sum = ...  max = ... [kg K m^-3 s^-1]
//   expected_surface_max = rho * Q_max / Cp * (1 - exp(-dz/alfg)) / dz
//
// With Q_max=14330 W/m2, rho~1.2, Cp=1004.64, alfg=45 m, dz=4 m:
//   expected_surface_max ~ 0.38 kg K m^-3 s^-1  (~0.32 K/s / rho)
inline void apply_fire_tendency_to_cc_source(
    amrex::MultiFab& cc_source,
    const amrex::MultiFab& Q_atm,
    const amrex::MultiFab* Q_lat_atm,
    const amrex::MultiFab& z_phys_cc,
    const amrex::MultiFab& S_old,
    const amrex::Geometry& geom_atm,
    amrex::Real alfg,
    amrex::Real feedback,
    bool has_moisture,
    bool fire_debug = false)
{
    using namespace amrex;

    if (alfg <= 0.0) { return; }
    if (feedback <= 0.0) { return; }

    constexpr Real Cp_d = 1004.64;
    constexpr Real L_V  = 2.501e6;

    const Box& domain = geom_atm.Domain();
    int klo = domain.smallEnd(2);
    int khi = domain.bigEnd(2);
    const auto& dx = geom_atm.CellSize();
    Real dz_avg = dx[2];

    for (MFIter mfi(cc_source, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const& cc_src_a = cc_source.array(mfi);
        auto const& q_atm_a  = Q_atm.array(mfi);
        auto const& z_cc_a   = z_phys_cc.array(mfi);
        auto const& s_old_a  = S_old.array(mfi);

        const Box& bx = mfi.tilebox();

        bool has_latent = (Q_lat_atm != nullptr) && has_moisture;
        Array4<Real const> q_lat_arr;
        if (has_latent) {
            q_lat_arr = Q_lat_atm->const_array(mfi);
        }

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(const IntVect& iv) {
            int i = iv[0];
            int j = iv[1];
            int k = iv[2];

            Real q_sfc_sensible = feedback * q_atm_a(i, j, 0);
            Real q_sfc_latent = 0.0;
            if (has_latent) {
                q_sfc_latent = feedback * q_lat_arr(i, j, 0);
            }

            Real z_sfc  = z_cc_a(i, j, klo);
            Real z_k    = z_cc_a(i, j, k) - z_sfc;
            Real exp_k  = std::exp(-z_k / alfg);

            Real z_kp1  = (k < khi) ? (z_cc_a(i, j, k + 1) - z_sfc) : (z_k + dz_avg);
            Real exp_kp1 = std::exp(-z_kp1 / alfg);

            Real dz = (k < khi) ? (z_cc_a(i, j, k + 1) - z_cc_a(i, j, k)) : dz_avg;
            if (dz <= 1.0e-10) { dz = dz_avg; }

            Real rho = s_old_a(i, j, k, Rho_comp);
            if (rho <= 0.0) { rho = 1.2; }

            Real hfx_k    = (q_sfc_sensible / Cp_d) * exp_k;
            Real hfx_kp1  = (q_sfc_sensible / Cp_d) * exp_kp1;
            Real theta_tend = -rho * (hfx_kp1 - hfx_k) / dz;
            cc_src_a(i, j, k, RhoTheta_comp) += theta_tend;

            if (has_latent) {
                Real qfx_k   = (q_sfc_latent / L_V) * exp_k;
                Real qfx_kp1 = (q_sfc_latent / L_V) * exp_kp1;
                Real q_tend  = -rho * (qfx_kp1 - qfx_k) / dz;
                cc_src_a(i, j, k, RhoQ1_comp) += q_tend;
            }
        });
    }

    // Tendency diagnostics: runs on host after GPU kernel completes.
    // Prints once per apply call (3x per timestep for RK3).
    // Expected rt_tend_max ~ rho * Q_max / Cp * (1 - exp(-dz/alfg)) / dz
    if (fire_debug) {
        Real Q_max        = Q_atm.max(0);
        Real rt_tend_sum  = cc_source.sum(RhoTheta_comp);
        Real rt_tend_max  = cc_source.max(RhoTheta_comp);
        Real expected_max = Q_max / (Cp_d * 1.2)
                            * (1.0 - std::exp(-dz_avg / alfg)) / dz_avg;
        amrex::Print() << "[FIRE COUPLING] Q_atm_max=" << Q_max
                       << " W/m2  alfg=" << alfg << " m\n"
                       << "[FIRE COUPLING] RhoTheta tendency"
                       << "  sum="      << rt_tend_sum
                       << "  max="      << rt_tend_max
                       << "  expected_surface_max=" << expected_max
                       << " [kg K m^-3 s^-1]\n";
    }
}

#endif  // ERF_FIRE_ATM_COUPLING_H