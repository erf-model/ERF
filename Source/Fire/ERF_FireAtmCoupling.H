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
//   3. Does the same for latent heat (grnqfx → rqvfrten):
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
//   This is stored as fire_latent_flux [W/m²] and distributed with the same
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
// WRF-Fire default: alfg ≈ 45-75 m. At z = alfg: flux decays to 1/e ≈ 37%.
constexpr amrex::Real FIRE_ATM_ALFG_DEFAULT = 45.0;

// Default e-folding height for canopy fire contribution [m].
// Only used when canopy fire heat is non-zero.
constexpr amrex::Real FIRE_ATM_ALFC_DEFAULT = 60.0;

// Coarsen fire_heat_flux from the fire grid (C*nx × C*ny × 1) to the
// atmospheric level-0 grid (nx × ny × 1) using area-weighted averaging.
// Uses amrex::average_down() with refinement ratio (C, C, 1).
// This is a local operation: each MPI rank owns the same x-y columns.
//
// Parameters:
//   Q_atm_out  — output: atmospheric-resolution heat flux [W/m²], nx×ny×1
//   Q_fire     — input:  fire-resolution heat flux [W/m²], C*nx×C*ny×1
//   geom_atm   — atmospheric geometry
//   geom_fire  — fire grid geometry
//   C          — refinement factor
inline void coarsen_fire_flux_to_atm(
    amrex::MultiFab& Q_atm_out,
    const amrex::MultiFab& Q_fire,
    const amrex::Geometry& geom_atm,
    const amrex::Geometry& geom_fire,
    int C)
{
    // Use AMReX's average_down with refinement ratio (C, C, 1) in z.
    amrex::average_down(Q_fire, Q_atm_out, geom_fire, geom_atm, 0, Q_fire.nComp(),
                        amrex::IntVect(C, C, 1));
}

// Compute fire latent heat flux [W/m²] from sensible heat flux and fuel moisture.
// Follows WRF-Fire: grnqfx = (bmst + (1-bmst)*0.56) * (dmass/dt) * xlv
// where bmst = M_f / (1 + M_f) is the ratio of water to wet fuel mass.
//
// In ERF-Fire, we approximate from the surface sensible heat flux Q_sfc:
//   Q_latent = Q_sfc * Lv/h_fuel * bmst / (1 - bmst)
//
// This is equivalent to WRF-Fire when h_fuel is the net combustion heat.
// The 0.56 factor (fraction of combustion product water from burning)
// is included per WRF-Fire.
//
// Parameters:
//   Q_lat_out  — output: latent heat flux [W/m²], nx×ny×1
//   Q_sfc      — input:  sensible heat flux [W/m²], nx×ny×1
//   M_f        — mean fuel moisture [fraction], scalar (domain average)
//   h_fuel_Jkg — heat content of dry fuel [J/kg]
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

    // Latent heat of vaporization [J/kg] at standard conditions
    constexpr amrex::Real L_V = 2.501e6;

    // Fuel moisture ratio: bmst = M_f / (1 + M_f)
    amrex::Real bmst = M_f / (1.0 + M_f);

    // Fraction of combustion product water from burning (WRF-Fire)
    constexpr amrex::Real f_water = 0.56;

    // Effective latent heat ratio from fuel heat content
    amrex::Real ratio = L_V / h_fuel_Jkg;

    // Per-cell computation of latent flux from sensible flux
    for (amrex::MFIter mfi(Q_sfc, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const& q_sfc_a = Q_sfc.array(mfi);
        auto const& q_lat_a = Q_lat_out.array(mfi);
        auto const& bx = mfi.tilebox();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(const amrex::IntVect& iv) {
            // Q_latent = Q_sensible * (Lv / h_fuel) * (bmst + (1-bmst)*0.56)
            amrex::Real factor = ratio * (bmst + (1.0 - bmst) * f_water);
            q_lat_a(iv) = factor * q_sfc_a(iv);
        });
    }
}

// Distribute fire surface heat flux through the atmospheric column and add
// the resulting tendency to cc_source.
//
// For each atmospheric column (i,j) and each k-level:
//   z_k   = z_phys_cc(i,j,k) - z_phys_cc(i,j,klo)   [height above surface, m]
//
//   Sensible flux profile (following WRF-Fire fire_tendency):
//     hfx(i,k,j) = (Q_atm(i,j) / Cp_d) * exp(-z_k / alfg)
//
//   Latent flux profile (same scale height):
//     qfx(i,k,j) = (Q_lat_atm(i,j) / Lv) * exp(-z_k / alfg)
//
//   Tendency from vertical divergence (WRF-Fire Eq., fire_tendency):
//     d(RhoTheta)/dt = -rho * (hfx(k+1) - hfx(k)) / dz_k
//     d(RhoQ1)/dt    = -rho * (qfx(k+1) - qfx(k)) / dz_k
//
//   These are added (accumulated) into cc_source[RhoTheta_comp] and
//   cc_source[RhoQ1_comp] respectively.
//
// The flux at k=0 (bottom face, k=klo) is the full surface flux Q_atm/Cp.
// The flux at k=khi+1 (top face) is exp(-z_top/alfg) ≈ 0 for z_top >> alfg.
//
// Parameters:
//   cc_source     — in/out: 3D source array for conserved variables
//   Q_atm         — input:  coarsened sensible heat flux [W/m²], nx×ny×1
//   Q_lat_atm     — input:  coarsened latent heat flux [W/m²], nx×ny×1 (may be null)
//   z_phys_cc     — input:  cell-centre terrain-following height [m], nx×ny×nz
//   S_old         — input:  atmospheric conserved state at time n (for density)
//   geom_atm      — atmospheric geometry
//   alfg          — e-folding height scale [m]
//   feedback      — multiplier applied to both fluxes (0 = no feedback, 1 = full)
//   has_moisture  — whether to inject into RhoQ1_comp (false if no moisture model)
inline void apply_fire_tendency_to_cc_source(
    amrex::MultiFab& cc_source,
    const amrex::MultiFab& Q_atm,
    const amrex::MultiFab* Q_lat_atm,
    const amrex::MultiFab& z_phys_cc,
    const amrex::MultiFab& S_old,
    const amrex::Geometry& geom_atm,
    amrex::Real alfg,
    amrex::Real feedback,
    bool has_moisture)
{
    using namespace amrex;

    if (alfg <= 0.0) { return; }
    if (feedback <= 0.0) { return; }

    // Physical constants
    constexpr Real Cp_d = 1004.64;  // specific heat at constant pressure [J/(kg·K)]
    constexpr Real L_V = 2.501e6;   // latent heat of vaporization [J/kg]

    const Box& domain = geom_atm.Domain();
    int klo = domain.smallEnd(2);
    int khi = domain.bigEnd(2);
    const auto& dx = geom_atm.CellSize();
    Real dz_avg = dx[2];  // average cell height

    // Iterate over the 3D atmospheric grid
    for (MFIter mfi(cc_source, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const& cc_src_a = cc_source.array(mfi);
        auto const& q_atm_a  = Q_atm.array(mfi);
        auto const& z_cc_a   = z_phys_cc.array(mfi);
        auto const& s_old_a  = S_old.array(mfi);

        const Box& bx = mfi.tilebox();

        // FIX: store Array4 by value to avoid taking address of temporary
        bool has_latent = (Q_lat_atm != nullptr) && has_moisture;
        Array4<Real const> q_lat_arr;
        if (has_latent) {
            q_lat_arr = Q_lat_atm->const_array(mfi);
        }

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(const IntVect& iv) {
            int i = iv[0];
            int j = iv[1];
            int k = iv[2];

            // Surface values (assume 2D fire flux is constant in z, at k=klo)
            Real q_sfc_sensible = feedback * q_atm_a(i, j, 0);
            Real q_sfc_latent = 0.0;
            if (has_latent) {
                q_sfc_latent = feedback * q_lat_arr(i, j, 0);
            }

            // Surface height (lowest cell centre)
            Real z_sfc = z_cc_a(i, j, klo);

            // Height of this cell centre above surface
            Real z_k = z_cc_a(i, j, k) - z_sfc;

            // Exponential decay factor at this cell
            Real exp_k = std::exp(-z_k / alfg);

            // Height of next cell centre above surface (or use dz_avg if at top)
            Real z_kp1 = (k < khi) ? (z_cc_a(i, j, k + 1) - z_sfc)
                                    : (z_k + dz_avg);
            Real exp_kp1 = std::exp(-z_kp1 / alfg);

            // Vertical grid spacing
            Real dz = (k < khi) ? (z_cc_a(i, j, k + 1) - z_cc_a(i, j, k))
                                : dz_avg;

            // Prevent division by zero
            if (dz <= 1.0e-10) { dz = dz_avg; }

            // Density at this level [kg/m³]
            // FIX: use Rho_comp macro from ERF_IndexDefines.H
            Real rho = s_old_a(i, j, k, Rho_comp);
            if (rho <= 0.0) { rho = 1.2; }  // fallback

            // Sensible heat flux divergence tendency
            // hfx(k) = (Q_sfc / Cp) * exp(-z_k / alfg)
            // d(hfx)/dz = (Q_sfc / Cp) * (-1/alfg) * exp(-z / alfg)
            // d(RhoTheta)/dt = -rho * d(hfx)/dz = rho * (Q_sfc / Cp / alfg) * exp(-z/alfg)
            // Using flux values at cell faces: hfx(k+1/2) - hfx(k-1/2)
            Real hfx_k   = (q_sfc_sensible / Cp_d) * exp_k;
            Real hfx_kp1 = (q_sfc_sensible / Cp_d) * exp_kp1;
            Real theta_tend = -rho * (hfx_kp1 - hfx_k) / dz;

            // FIX: use RhoTheta_comp macro from ERF_IndexDefines.H
            cc_src_a(i, j, k, RhoTheta_comp) += theta_tend;

            // Latent heat flux divergence tendency (same structure)
            // FIX: use RhoQ1_comp macro from ERF_IndexDefines.H
            if (has_latent) {
                Real qfx_k   = (q_sfc_latent / L_V) * exp_k;
                Real qfx_kp1 = (q_sfc_latent / L_V) * exp_kp1;
                Real q_tend  = -rho * (qfx_kp1 - qfx_k) / dz;

                cc_src_a(i, j, k, RhoQ1_comp) += q_tend;
            }
        });
    }
}

#endif  // ERF_FIRE_ATM_COUPLING_H