// One-way fire-to-atmosphere coupling: distribute fire heat flux through
// the atmospheric column and inject into cc_source.
//
// This mirrors the WRF-Fire approach in module_fr_sfire_atm.F
// (subroutine fire_tendency), which:
//   1. Builds a vertical flux profile with exponential decay:
//        hfx(z) = (Q_sfc / Cp) * exp(-z / alfg)
//      where z is the height above the lowest cell centre of the column (half
//      a cell above the terrain: the profile's datum is that centre, so the
//      first cell takes the share 1 - e^{-dz/a}, 0.105 for dz = 5 m and
//      a = 45 m, and nothing is lost below it) and alfg is the e-folding
//      scale [m].
//   2. Converts the flux profile to a theta tendency via flux divergence:
//        d(RhoTheta)/dt = -d(hfx)/dz
//      with dz the thickness of the cell (face to face), so the column
//      receives the flux whatever the vertical spacing; erf.fire.
//      heat_tendency_density multiplies by rho and erf.fire.heat_tendency_exner
//      divides by the Exner function (see below).
//   3. Does the same for latent heat (grnqfx -> rqvfrten):
//        qfx(z) = (Q_lat / Lv) * exp(-z / alfg)
//        d(RhoQ1)/dt   = -d(qfx)/dz
//
// The e-folding scale alfg controls how deep into the atmosphere the fire
// heating penetrates. With alfg = 45 m, the flux has decayed to 37% at 45 m,
// 5% at 135 m, and 1% at 225 m above the surface.
//
// Reference:
//   Mandel, J., Beezley, J. D., Coen, J. L. and Kim, M. (2009). Data assimilation
//   for wildland fires. IEEE Control Systems Magazine, 29(3), 47-65; and
//   Mandel, J., Beezley, J. D. and Kochanski, A. K. (2011). Coupled
//   atmosphere-wildland fire modeling with WRF 3.3 and SFIRE 2011. Geoscientific
//   Model Development, 4, 591-610.
//   WRF-SFIRE source: phys/module_fr_sfire_atm.F, subroutine fire_tendency.
//
// Coupling time lag:
//   ERF uses one-step explicit lag: the fire heat flux computed at the end of
//   step n is held in Q_atm_prev and injected into cc_source during step n+1,
//   from the slow right-hand side of every RK stage (ERF_TI_slow_rhs_pre.H),
//   after make_sources has rebuilt cc_source for that stage.
//   This is consistent with WRF-Fire's explicit coupling approach.
//
// RK-stage safety:
//   apply_fire_tendency_to_cc_source is invoked once per RK stage of the
//   slow-RHS integrator. Because the lagged Q_atm_prev is fixed for the whole
//   timestep, the fire tendency it produces must also be identical at every
//   RK stage. With erf.fire.source_mode = overwrite (the default) we OVERWRITE
//   cc_source[RhoTheta_comp] (and RhoQ1_comp when latent heat is on): the
//   values are cleared to zero at the top of the routine and then set via '='
//   inside the kernel, which discards every other source make_sources put
//   there (radiative heating, Rayleigh damping, nudging). source_mode = add
//   adds to the rebuilt source instead. Accumulating with '+=' across RK
//   stages without the per-stage rebuild would triple the injected heat.
//
// Latent heat:
//   WRF-Fire computes grnqfx (ground latent heat flux) alongside grnhfx and
//   passes both through fire_tendency. Here the fire engine computes the
//   latent flux (from the fuel moisture and the water of combustion of the fuel
//   it burns) and hands it over in FireFluxes::latent [W/m2]; the coupling
//   distributes it with the same exponential profile.
//
// Assumptions:
//   - The tendency goes into the source of the level the fire grid refines
//     (erf.fire.anchor_level, the finest level by default); coarser levels see
//     the result through average-down with erf.coupling_type = TwoWay (the
//     default), not with OneWay; finer levels not at all.
//   - The engine's sensible and latent fluxes (FireFluxes) are on the fire grid
//     (C*nx x C*ny x 1 over the fire region).
//   - The coarsening from fire grid to atmospheric grid is the area average of
//     average_down(), with the region's lower corner added (FireGrid::atm_lo).
//   - z_phys_cc is the cell-centre terrain-following height MultiFab [m].
//   - The atmospheric BoxArray is not split in z (max_grid_size_z = nz enforced).
//   - Latent heat injection into RhoQ1_comp is skipped when moisture is inactive.
//   - fire_atm_feedback multiplier (equivalent to WRF-Fire's fire_atm_feedback)
//     scales both sensible and latent flux. Set to 0 for one-way coupling only
//     (fire responds to wind but atmosphere does not respond to fire heat).
//   - The RhoTheta / RhoQ1 slots of cc_source are OWNED by the fire coupling
//     for the duration of a timestep. If any other physics writes into those
//     same slots per stage, this routine must be moved out of the per-stage
//     slow-RHS path and called once before mri_integrator.advance() instead.

#ifndef ERF_FIRE_ATM_COUPLING_H
#define ERF_FIRE_ATM_COUPLING_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_Box.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFabUtil.H>
#include <cmath>
#include <limits>
#include <ERF_Constants.H>
#include <ERF_IndexDefines.H>
#include <ERF_EOS.H>
#include <ERF_FireGrid.H>

// Default e-folding height for ground fire heat distribution [m].
// WRF-Fire default: alfg ~ 45-75 m. At z = alfg: flux decays to 1/e ~ 37%.
// The crown-fire release is added to the same surface flux and injected with
// this one profile; there is no separate canopy e-folding height.
constexpr amrex::Real FIRE_ATM_ALFG_DEFAULT = amrex::Real(45.0);

// Coarsen a fire-grid field (C*nx x C*ny x 1) onto the atmospheric columns of
// the level the fire grid refines (nx x ny x 1) by area averaging.
//
// Q_atm_out lives on the k = 0 slabs of that level's boxes, in their order and
// with their owners; Q_fire on the fire grid, whose index space starts at the
// region's corner. On level 0 the corner is the origin and this is
// amrex::average_down() with ratio (C, C, 1), as it always was. Elsewhere the
// same average is taken with the corner added: fire cells (i - lo_x) C + r,
// (j - lo_y) C + q for r, q in [0, C) average into atmospheric cell (i, j).
// average_down() itself would place the result at the fire grid's own indices.
inline void coarsen_fire_flux_to_atm(
    amrex::MultiFab& Q_atm_out,
    const amrex::MultiFab& Q_fire,
    const amrex::Geometry& geom_atm,
    const FireGrid& fg)
{
    const int C = fg.C;
    if (fg.atm_lo == amrex::IntVect::TheZeroVector()) {
        amrex::average_down(Q_fire, Q_atm_out, fg.geom, geom_atm, 0, Q_fire.nComp(),
                            amrex::IntVect(C, C, 1));
        return;
    }

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        Q_fire.boxArray().size() == Q_atm_out.boxArray().size() &&
        Q_fire.DistributionMap() == Q_atm_out.DistributionMap(),
        "[FIRE] coarsen_fire_flux_to_atm: the fire grid and the atmospheric slabs must pair box for box");
    AMREX_ALWAYS_ASSERT(Q_atm_out.nComp() >= Q_fire.nComp());

    const int lo_x  = fg.atm_lo[0];
    const int lo_y  = fg.atm_lo[1];
    const int ncomp = Q_fire.nComp();
    const amrex::Real volfrac = amrex::Real(1.0) / amrex::Real(C * C);
    for (amrex::MFIter mfi(Q_atm_out); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto crse = Q_atm_out.array(mfi);
        auto fine = Q_fire.const_array(mfi);
        amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int ii = (i - lo_x) * C;
            const int jj = (j - lo_y) * C;
            amrex::Real c = amrex::Real(0.0);
            for (int q = 0; q < C; ++q) {
                for (int r = 0; r < C; ++r) {
                    c += fine(ii + r, jj + q, 0, n);
                }
            }
            crse(i, j, k, n) = volfrac * c;
        });
    }
}

// Height of the lower (upper = 0) or upper (upper = 1) face of cell (i, j, k):
// the mean of its four nodes where the nodal heights exist, half a cell below
// (above) the centre on a uniform grid without them.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real fire_cell_face (amrex::Array4<amrex::Real const> const& z_nd,
                            amrex::Array4<amrex::Real const> const& z_cc,
                            bool has_nd, amrex::Real dz_avg,
                            int i, int j, int k, int upper) noexcept
{
    if (has_nd) {
        const int kf = k + upper;
        return amrex::Real(0.25) * (z_nd(i, j, kf) + z_nd(i+1, j, kf)
                                  + z_nd(i, j+1, kf) + z_nd(i+1, j+1, kf));
    }
    return z_cc(i, j, k) + (upper ? amrex::Real(0.5) : amrex::Real(-0.5)) * dz_avg;
}

// Thickness of a cell from its faces. Only a degenerate (zero or inverted)
// cell is floored, at 1e-6 of the mean spacing, so that the division stays
// finite; a thin first cell on a stretched grid keeps its own thickness (a
// floor of 0.1 dz_avg put 59 % of the heat into a column with a 1 m first
// cell stretched by 1.1 over 64 cells, at the default alfg of 45 m).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real fire_cell_thickness (amrex::Real z_lo, amrex::Real z_hi, amrex::Real dz_avg) noexcept
{
    return amrex::max(z_hi - z_lo, amrex::Real(1.0e-6) * dz_avg);
}

// Distribute fire surface heat flux through the atmospheric column and write
// the resulting tendency into cc_source.
//
// Behaviour: OVERWRITES cc_source[RhoTheta_comp] and (when latent heat is on)
// cc_source[RhoQ1_comp]. See file-level "RK-stage safety" note above.
//
// Optional fire_debug argument prints per-call tendency diagnostics to stdout:
//   [FIRE COUPLING] RhoTheta tendency sum = ...  max = ... [kg K m^-3 s^-1]
//   expected_surface_max = feedback * Q_max / Cp * (1 - exp(-dz/alfg)) / dz
//   (the default heat_tendency_density = false form, d(rho theta)/dt =
//   -(1/Cp) dQ/dz, for the lowest cell on a uniform grid; with the density
//   factor the tendency carries an extra rho, with the Exner factor 1/Pi)
//
// With Q_max=14330 W/m2, Cp=Cp_d (1004.5 J/(kg K)), alfg=45 m, dz=5 m:
//   expected_surface_max ~ 0.30 kg K m^-3 s^-1
//
// z_phys_nd, the nodal heights of the level, gives each cell's faces: its
// thickness divides the tendency, and the open-fraction placement measures the
// roof height from the ground. ERF always passes it; with nullptr (callers
// without nodal heights, such as the unit tests) the faces are half a uniform
// cell from the centres, which is right only on a uniform grid.
inline void apply_fire_tendency_to_cc_source(
    amrex::MultiFab& cc_source,
    const amrex::MultiFab& Q_atm,
    const amrex::MultiFab* Q_lat_atm,
    const amrex::MultiFab& z_phys_cc,
    const amrex::MultiFab* z_phys_nd,
    const amrex::MultiFab& S_old,
    const amrex::Geometry& geom_atm,
    amrex::Real alfg,
    amrex::Real feedback,
    bool has_moisture,
    bool fire_debug = false,
    bool additive = false,
    const amrex::MultiFab* open_frac_atm = nullptr,
    const amrex::MultiFab* roof_h_atm = nullptr,
    bool use_open_weights = false,
    bool density_factor = false,
    bool exner_factor = false,
    bool bottom_cell = false)
{
    using namespace amrex;

    if (alfg <= Real(0.0)) { return; }
    if (feedback <= Real(0.0)) { return; }

    // Source conversion must use the same dry-air heat capacity as ERF's EOS.
    constexpr Real Cp_fire = Cp_d;
    // and the latent heat ERF uses (an engine converting a vapour mass flux to
    // its latent flux multiplies by this same L_v, so the water round-trips exactly)
    constexpr Real L_V  = L_v;
    // Exner function of the cell from its rho theta (erf.fire.heat_tendency_exner):
    // heating q [W/m3] gives d(theta)/dt = q / (rho cp Pi), so the rho theta
    // tendency carries 1/Pi; without it the enthalpy injected is Pi times the
    // flux, 4.5 % short at 850 hPa. WRF-SFIRE omits the factor; it is on by
    // default here.
    const Real rdOcp = R_d / Cp_d;

    const Box& domain = geom_atm.Domain();
    int klo = domain.smallEnd(2);
    int khi = domain.bigEnd(2);
    const auto& dx = geom_atm.CellSize();
    Real dz_avg = dx[2];

    const bool has_latent_global = (Q_lat_atm != nullptr) && has_moisture;
    // The open fraction and roof height are handed in whenever the engine
    // reports them; the weights only use them when asked.
    const bool have_open = (open_frac_atm != nullptr) && (roof_h_atm != nullptr);
    const bool use_open  = have_open && use_open_weights;

    // Debug only: the fire tendency alone, times the cell thickness (in additive mode
    // cc_source also holds the other sources), for the energy check below.
    MultiFab tend_dbg;
    if (fire_debug) {
        tend_dbg.define(cc_source.boxArray(), cc_source.DistributionMap(), 1, 0);
        tend_dbg.setVal(0.0);
    }

    // Legacy (overwrite) mode: zero the components we own before writing them,
    // so any state left over from a previous stage does not leak through. In
    // additive mode make_sources has just rebuilt the source for this stage,
    // so the fire tendency is accumulated once per stage on top of it.
    if (!additive) {
        cc_source.setVal(0.0, RhoTheta_comp, 1, cc_source.nGrowVect());
        if (has_latent_global) {
            cc_source.setVal(0.0, RhoQ1_comp, 1, cc_source.nGrowVect());
        }
    }

    // Column normalisation of the open-fraction weights: one k = klo slab per
    // box, filled once per column (the boxes hold whole columns, which the
    // fire prerequisites assert), before the tiled loop that reads it.
    // The engine's open fraction and roof height are read with the level's
    // MFIter below, so they must hold one k = klo slab per box of the level,
    // in its order and on its DistributionMapping (only a debug build checks
    // that inside AMReX).
    if (have_open) {
        const BoxArray& cba = cc_source.boxArray();
        bool same_layout = open_frac_atm->DistributionMap() == cc_source.DistributionMap() &&
                           roof_h_atm->DistributionMap()    == cc_source.DistributionMap() &&
                           open_frac_atm->boxArray().size() == cba.size() &&
                           roof_h_atm->boxArray().size()    == cba.size();
        for (int ib = 0; same_layout && ib < cba.size(); ++ib) {
            const Box slab = makeSlab(cba[ib], 2, klo);
            same_layout = (open_frac_atm->boxArray()[ib] == slab) && (roof_h_atm->boxArray()[ib] == slab);
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(same_layout,
            "[FIRE] the engine's heat_open_fraction_atm and roof_height_atm must hold one k = 0 slab "
            "per box of the fire grid's level, in its box order and on its DistributionMapping");
    }
    MultiFab norm_mf;
    if (use_open) {
        norm_mf.define(open_frac_atm->boxArray(), cc_source.DistributionMap(), 1, 0);
        for (MFIter mfi(norm_mf); mfi.isValid(); ++mfi) {
            auto const& norm_a = norm_mf.array(mfi);
            auto const& open_arr = open_frac_atm->const_array(mfi);
            auto const& roof_arr = roof_h_atm->const_array(mfi);
            auto const& z_cc_a   = z_phys_cc.const_array(mfi);
            const bool has_nd = (z_phys_nd != nullptr);
            Array4<Real const> z_nd_a;
            if (has_nd) { z_nd_a = z_phys_nd->const_array(mfi); }
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int /*kslab*/) noexcept {
                const Real f_open = open_arr(i, j, 0);
                const Real h_roof = roof_arr(i, j, 0);   // above the ground
                const Real z_sfc  = z_cc_a(i, j, klo);
                const Real ground = fire_cell_face(z_nd_a, z_cc_a, has_nd, dz_avg, i, j, klo, 0);
                Real norm = Real(0.0);
                for (int kk = klo; kk <= khi; ++kk) {
                    const Real zz    = z_cc_a(i, j, kk) - z_sfc;
                    const Real zzp1  = (kk < khi) ? (z_cc_a(i, j, kk + 1) - z_sfc) : (zz + dz_avg);
                    Real s = Real(1.0);
                    if (h_roof > Real(0.0)) {
                        const Real zlo2 = fire_cell_face(z_nd_a, z_cc_a, has_nd, dz_avg, i, j, kk, 0) - ground;
                        const Real zhi2 = fire_cell_face(z_nd_a, z_cc_a, has_nd, dz_avg, i, j, kk, 1) - ground;
                        const Real dzz  = fire_cell_thickness(zlo2, zhi2, dz_avg);
                        if (zhi2 <= h_roof)     { s = f_open; }
                        else if (zlo2 < h_roof) { s = f_open + (Real(1.0) - f_open) * (zhi2 - h_roof) / dzz; }
                    }
                    norm += (std::exp(-zz / alfg) - std::exp(-zzp1 / alfg)) * s;
                }
                norm_a(i, j, klo) = norm;
            });
        }
    }

    for (MFIter mfi(cc_source, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const& cc_src_a = cc_source.array(mfi);
        auto const& q_atm_a  = Q_atm.array(mfi);
        auto const& z_cc_a   = z_phys_cc.array(mfi);
        auto const& s_old_a  = S_old.array(mfi);
        // Faces of cell k: the mean of its four lower (upper) nodes where the
        // nodal heights exist, half a cell below (above) the centre otherwise
        const bool has_nd = (z_phys_nd != nullptr);
        Array4<Real const> z_nd_a;
        if (has_nd) { z_nd_a = z_phys_nd->const_array(mfi); }

        const Box& bx = mfi.tilebox();

        const bool has_latent = has_latent_global;
        Array4<Real const> q_lat_arr;
        if (has_latent) {
            q_lat_arr = Q_lat_atm->const_array(mfi);
        }
        Array4<Real const> open_arr, roof_arr;
        if (use_open) {
            open_arr = open_frac_atm->const_array(mfi);
            roof_arr = roof_h_atm->const_array(mfi);
        }
        Array4<Real> dbg_arr;
        if (fire_debug) { dbg_arr = tend_dbg.array(mfi); }

        // Column normalisation of the open-fraction weights, computed once per
        // column above the tile loop.
        Array4<Real const> norm_arr;
        if (use_open) { norm_arr = norm_mf.const_array(mfi); }

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(const IntVect& iv) {
            int i = iv[0];
            int j = iv[1];
            int k = iv[2];

            Real q_sfc_sensible = feedback * q_atm_a(i, j, 0);
            Real q_sfc_latent = Real(0.0);
            if (has_latent) {
                q_sfc_latent = feedback * q_lat_arr(i, j, 0);
            }

            Real z_sfc  = z_cc_a(i, j, klo);
            Real z_k    = z_cc_a(i, j, k) - z_sfc;
            Real exp_k  = std::exp(-z_k / alfg);

            Real z_kp1  = (k < khi) ? (z_cc_a(i, j, k + 1) - z_sfc) : (z_k + dz_avg);
            Real exp_kp1 = std::exp(-z_kp1 / alfg);

            // The cell's own thickness, face to face: the column then receives
            // sum_k w_k of the flux, 1 - exp(-z_top/alfg), on any vertical
            // spacing (the centre-to-centre spacing used before put a share
            // dz_cell / dz_cc into each cell, 2 / (1 + r) of the flux on a grid
            // stretched by r)
            const Real z_lo = fire_cell_face(z_nd_a, z_cc_a, has_nd, dz_avg, i, j, k, 0);
            const Real z_hi = fire_cell_face(z_nd_a, z_cc_a, has_nd, dz_avg, i, j, k, 1);
            const Real dz   = fire_cell_thickness(z_lo, z_hi, dz_avg);

            Real rho = s_old_a(i, j, k, Rho_comp);
            if (rho <= Real(0.0)) { rho = Real(1.2); }
            // Heating gives d(rho theta)/dt = -(1/Cp) dQ/dz, with no density
            // factor. The legacy form carries an extra rho, so its injected
            // energy is rho times the flux; density_factor = false drops it.
            // The Exner factor divides the tendency by Pi(p) so that the
            // enthalpy, not rho theta, integrates to the flux.
            // Preserve the explicit legacy density option, but apply Exner
            // only to heat. Q_lat / L_V is a water mass flux: converting it
            // to d(rho_d q_v)/dt must not multiply that mass by 1/Pi.
            const Real mass_fac = density_factor ? rho : Real(1.0);
            Real heat_fac = mass_fac;
            if (exner_factor) {
                const Real rhotheta = amrex::max(s_old_a(i, j, k, RhoTheta_comp), Real(1.0e-10));
                // the moist pressure where there is vapour (about 1 % in Pi at 15 g/kg)
                const Real qv = has_moisture ? amrex::max(s_old_a(i, j, k, RhoQ1_comp), Real(0.0)) / rho : Real(0.0);
                const Real p = getPgivenRTh(rhotheta, qv);
                const Real exner = std::pow(p / p_0, rdOcp);
                heat_fac /= amrex::max(exner, Real(1.0e-6));
            }

            // Layer weight: the fraction of the column's exponential profile
            // that falls in this cell. Without an open fraction this is the flux
            // difference form used since the coupling was written. With the
            // open-fraction option, cells below the sampled roof height carry
            // only their open share and the column is renormalised over the
            // cells that remain, so the injected energy per column is unchanged
            // and none of it lands inside a building.
            Real w_k = exp_k - exp_kp1;
            // erf.fire.heat_profile = bottom_cell: the whole flux enters the
            // lowest cell of the column (the LSFIRE coupling)
            if (bottom_cell) { w_k = (k == klo) ? Real(1.0) : Real(0.0); }
            if (use_open) {
                const Real f_open = open_arr(i, j, 0);
                const Real h_roof = roof_arr(i, j, 0);
                // Open share of this cell: fully open above the roof, f_open
                // below it, linear across the cell that holds the roof.
                // The roof height is above the ground, so are the cell's faces.
                Real share = Real(1.0);
                if (h_roof > Real(0.0)) {
                    const Real ground = fire_cell_face(z_nd_a, z_cc_a, has_nd, dz_avg, i, j, klo, 0);
                    const Real zlo_g  = z_lo - ground;
                    const Real zhi_g  = z_hi - ground;
                    if (zhi_g <= h_roof)      { share = f_open; }
                    else if (zlo_g < h_roof)  { share = f_open + (Real(1.0) - f_open) * (zhi_g - h_roof) / dz; }
                }
                // Column normalisation over the open shares of the same
                // profile, precomputed per column above.
                const Real norm = norm_arr(i, j, klo);
                // The legacy profile integrates to (1 - exp(-z_top/alfg)); keep
                // that total so the two modes inject the same energy per column.
                const Real z_top   = (z_cc_a(i, j, khi) - z_sfc) + dz_avg;
                const Real legacy  = Real(1.0) - std::exp(-z_top / alfg);
                w_k = (norm > Real(1.0e-30)) ? w_k * share * legacy / amrex::max(norm, Real(1.0e-30)) : Real(0.0);
            }

            // Without the open-fraction option the flux-difference expression
            // is kept (exp_k - exp_kp1 times the flux, over the cell thickness).
            Real theta_tend, q_tend = Real(0.0);
            if (use_open || bottom_cell) {
                theta_tend = heat_fac * (q_sfc_sensible / Cp_fire) * w_k / dz;
                if (has_latent) { q_tend = mass_fac * (q_sfc_latent / L_V) * w_k / dz; }
            } else {
                Real hfx_k    = (q_sfc_sensible / Cp_fire) * exp_k;
                Real hfx_kp1  = (q_sfc_sensible / Cp_fire) * exp_kp1;
                theta_tend = -heat_fac * (hfx_kp1 - hfx_k) / dz;
                if (has_latent) {
                    Real qfx_k   = (q_sfc_latent / L_V) * exp_k;
                    Real qfx_kp1 = (q_sfc_latent / L_V) * exp_kp1;
                    q_tend  = -mass_fac * (qfx_kp1 - qfx_k) / dz;
                }
            }
            if (additive) {
                cc_src_a(i, j, k, RhoTheta_comp) += theta_tend;
            } else {
                // OVERWRITE (see file-level RK-stage safety note).
                cc_src_a(i, j, k, RhoTheta_comp) = theta_tend;
            }
            if (has_latent) {
                if (additive) {
                    cc_src_a(i, j, k, RhoQ1_comp) += q_tend;
                } else {
                    cc_src_a(i, j, k, RhoQ1_comp) = q_tend;
                }
            }
            // the cell's true thickness (not the floored one), so a degenerate
            // cell's lost energy shows in the ratio
            if (dbg_arr) { dbg_arr(i, j, k) = theta_tend * (z_hi - z_lo); }
        });
    }

    // Tendency diagnostics: runs on host after GPU kernel completes.
    // Prints once per apply call (3x per timestep for RK3).
    // Expected rt_tend_max ~ feedback * Q_max / Cp * (1 - exp(-dz/alfg)) / dz
    // (times rho with the density factor) on a uniform grid
    if (fire_debug) {
        Real Q_max        = Q_atm.max(0);
        Real rt_tend_sum  = cc_source.sum(RhoTheta_comp);
        Real rt_tend_max  = cc_source.max(RhoTheta_comp);
        Real expected_max = feedback * Q_max / Cp_fire
                            * (Real(1.0) - std::exp(-dz_avg / alfg)) / dz_avg;
        amrex::Print() << "[FIRE COUPLING] Q_atm_max=" << Q_max
                       << " W/m2  alfg=" << alfg << " m\n"
                       << "[FIRE COUPLING] RhoTheta tendency"
                       << "  sum="      << rt_tend_sum
                       << "  max="      << rt_tend_max
                       << "  expected_surface_max=" << expected_max
                       << " [kg K m^-3 s^-1]\n";

        // Energy check: Cp times the sum over cells of (rho theta tendency) dz
        // per unit area, against the flux handed in. The profile carries
        // 1 - exp(-z_top/alfg) of the flux, essentially all of it for any
        // realistic domain height, so with neither the density nor the Exner
        // factor the ratio is that (the density factor multiplies it by rho,
        // the Exner factor by 1/Pi, both profile-weighted means, essentially
        // their values near the surface). The largest fire tendency inside a fully blocked cell below the roof must be
        // zero with the open-fraction option.
        const Real area = dx[0] * dx[1];
        const Real e_out = Cp_fire * tend_dbg.sum(0) * area;
        const Real e_in  = feedback * Q_atm.sum(0) * area;
        // Columns that straddle a footprint (0 < open fraction < 1) are where
        // the placement matters: fully covered columns carry no flux once their
        // fire cells are masked. For those partial columns report the largest
        // share of the column's heating that lands below the roof H (above the
        // ground). With B the plain profile's share below the roof (measured
        // from the first cell centre) and L = 1 - exp(-z_top/alfg), the
        // open-fraction weights give f B L / (L - (1 - f) B), more than f B
        // because the column is renormalised; and the largest potential temperature of the
        // state in cells below a roof, which shows heat that earlier steps left
        // inside buildings.
        Real below_roof_max = -Real(1.0), theta_blocked_max = -Real(1.0);
        if (have_open) {
            MultiFab col(open_frac_atm->boxArray(), open_frac_atm->DistributionMap(), 1, 0);
            MultiFab thb(tend_dbg.boxArray(), tend_dbg.DistributionMap(), 1, 0);
            col.setVal(0.0); thb.setVal(0.0);
            for (MFIter mfi(col, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx2 = mfi.tilebox();
                auto const& c   = col.array(mfi);
                auto const& t   = tend_dbg.const_array(mfi);
                auto const& z   = z_phys_cc.const_array(mfi);
                auto const& ofrac  = open_frac_atm->const_array(mfi);
                auto const& hr  = roof_h_atm->const_array(mfi);
                const bool has_nd = (z_phys_nd != nullptr);
                Array4<Real const> znd;
                if (has_nd) { znd = z_phys_nd->const_array(mfi); }
                amrex::ParallelFor(bx2, [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/) noexcept {
                    const Real f = ofrac(i, j, 0);
                    if (f <= Real(1.0e-6) || f >= Real(1.0) - Real(1.0e-6)) { return; }
                    const Real ground = fire_cell_face(znd, z, has_nd, dz_avg, i, j, klo, 0);
                    Real below = Real(0.0), total = Real(0.0);
                    for (int kk = klo; kk <= khi; ++kk) {
                        const Real zk = z(i, j, kk) - ground;
                        total += t(i, j, kk);
                        if (zk < hr(i, j, 0)) { below += t(i, j, kk); }
                    }
                    c(i, j, 0) = (total > Real(1.0e-30)) ? below / amrex::max(total, Real(1.0e-30)) : Real(0.0);
                });
            }
            for (MFIter mfi(thb, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                auto const& b   = thb.array(mfi);
                auto const& z   = z_phys_cc.const_array(mfi);
                auto const& so  = S_old.const_array(mfi);
                auto const& ofrac  = open_frac_atm->const_array(mfi);
                auto const& hr  = roof_h_atm->const_array(mfi);
                const bool has_nd = (z_phys_nd != nullptr);
                Array4<Real const> znd;
                if (has_nd) { znd = z_phys_nd->const_array(mfi); }
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    const Real zk = z(i, j, k) - fire_cell_face(znd, z, has_nd, dz_avg, i, j, klo, 0);
                    if (ofrac(i, j, 0) < Real(1.0) - Real(1.0e-6) && zk < hr(i, j, 0)) {
                        b(i, j, k) = so(i, j, k, RhoTheta_comp) / so(i, j, k, Rho_comp);
                    }
                });
            }
            below_roof_max    = col.max(0);
            theta_blocked_max = thb.max(0);
        }
        // The ratio divides only when e_in > 0, but the optimiser evaluates the division
        // speculatively for both arms of the selection, and 0/0 on the first lagged step
        // raised the invalid flag that aborts the run under amrex.fpe_trap_invalid; the
        // floored denominator (the smallest normal Real, so it holds in single
        // precision too) keeps the discarded value finite and the printed one unchanged.
        amrex::Print() << "[FIRE COUPLING] energy_in=" << e_in << " W  energy_out=" << e_out
                       << " W  ratio=" << ((e_in > 0.0) ? e_out / amrex::max(e_in, std::numeric_limits<Real>::min()) : 0.0)
                       << "  below_roof_share_max=" << below_roof_max
                       << "  theta_blocked_max=" << theta_blocked_max << "\n";
    }
}

#endif  // ERF_FIRE_ATM_COUPLING_H