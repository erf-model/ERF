#include <ERF_FireWindExtract.H>


using namespace amrex;

void compute_terrain_curvature(
    MultiFab& curvature,
    const MultiFab& fire_slopes,
    const Geometry& geom_fire)
{
    // Curvature = d²z/dx² + d²z/dy²
    // For now, implement as a stub that computes from slope field
    // using finite differences

    Real dx = geom_fire.CellSize(0);
    Real dy = geom_fire.CellSize(1);

    for (MFIter mfi(curvature, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> curv = curvature.array(mfi);
        Array4<const Real> slopes = fire_slopes.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            int i = iv[0];
            int j = iv[1];
            int k = 0;

            // dz/dx at (i±1, j)
            Real dzdx_ip = slopes(i+1, j, k, 0);  // dz/dx component
            Real dzdx_im = slopes(i-1, j, k, 0);
            Real d2zdx2 = (dzdx_ip - dzdx_im) / (2.0 * dx);

            // dz/dy at (i, j±1)
            Real dzdy_jp = slopes(i, j+1, k, 1);  // dz/dy component
            Real dzdy_jm = slopes(i, j-1, k, 1);
            Real d2zdy2 = (dzdy_jp - dzdy_jm) / (2.0 * dy);

            curv(i, j, k) = d2zdx2 + d2zdy2;
        });
    }
}

void apply_farsite_terrain_wind(
    MultiFab& fire_wind_eff,
    const MultiFab& fire_slopes,
    const MultiFab& fire_curvature,
    Real k_ridge, Real k_shelter,
    Real k_valley, Real k_deflect)
{
    // Apply FARSITE terrain wind corrections per Finney (1998)
    // Ridge speed-up, lee sheltering, valley channeling, directional deflection
    
    for (MFIter mfi(fire_wind_eff, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> wind = fire_wind_eff.array(mfi);
        Array4<const Real> slopes = fire_slopes.array(mfi);
        Array4<const Real> curv = fire_curvature.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            int i = iv[0];
            int j = iv[1];
            int k = 0;

            // Read terrain properties
            Real sx = slopes(i, j, k, 0);  // dz/dx
            Real sy = slopes(i, j, k, 1);  // dz/dy
            Real curv_val = curv(i, j, k);
            Real ux = wind(i, j, k, 0);
            Real uy = wind(i, j, k, 1);

            // Compute magnitudes
            Real slope_mag = std::sqrt(sx*sx + sy*sy);
            Real wind_mag = std::sqrt(ux*ux + uy*uy);

            // Compute wind-upslope alignment: cos(angle between wind and slope gradient)
            // cos_upslope > 0 means wind pointing upslope (against gravity, climbing)
            Real cos_upslope = 0.0;
            if (slope_mag * wind_mag > 1.0e-10) {
                cos_upslope = (sx*ux + sy*uy) / (slope_mag * wind_mag);
            }

            // Classify terrain position and compute speed factor.
            // Slope-based classification uses wind-slope alignment (cos_upslope)
            // to identify windward (upslope) and lee (downslope) faces.
            // Curvature is used only for valley detection (concave terrain).
            // NOTE: uniform planar slopes have zero curvature, so curvature
            // must NOT be required for ridge/shelter classification.
            Real factor = 1.0;

            // Windward slope: wind climbing the slope (upslope face / ridge)
            if (slope_mag > 0.05 && cos_upslope > 0.0) {
                factor = k_ridge;
            }
            // Shelter: lee-side slope (wind blowing downslope)
            else if (slope_mag > 0.05 && cos_upslope < -0.5) {
                factor = k_shelter;
            }
            // Valley: concave terrain (negative curvature)
            else if (curv_val < -0.01) {
                factor = k_valley + (1.0_rt - k_valley) * (1.0_rt - amrex::min(slope_mag, 1.0_rt));
            }
            // Flat/other: no correction
            else {
                factor = 1.0;
            }

            // Clamp factor to reasonable range
            factor = amrex::max(0.1_rt, amrex::min(factor, 3.0_rt));

            // Apply speed scaling: scale wind magnitude
            ux *= factor;
            uy *= factor;
            Real wind_mag_new = wind_mag * factor;

            // Apply valley wind deflection toward slope aspect
            // Only when in valley (curv < -0.01) and on sloped terrain (slope_mag > 0.05)
            if (curv_val < -0.01 && slope_mag > 0.05) {
                // Compute z-component of slope × wind cross product
                Real sin_cross = sx * uy - sy * ux;
                
                // Compute deflection angle
                Real denom = amrex::max(slope_mag * wind_mag_new * factor, 1.0e-10);
                Real sin_cross_norm = amrex::max(-1.0_rt, amrex::min(1.0_rt, sin_cross / denom));
                Real deflect_angle = k_deflect * std::asin(sin_cross_norm);

                // Apply 2D rotation matrix to deflect wind direction
                Real cos_angle = std::cos(deflect_angle);
                Real sin_angle = std::sin(deflect_angle);
                Real ux_new = cos_angle * ux - sin_angle * uy;
                Real uy_new = sin_angle * ux + cos_angle * uy;
                ux = ux_new;
                uy = uy_new;
            }

            // Write back
            wind(i, j, k, 0) = ux;
            wind(i, j, k, 1) = uy;
        });
    }
}

/**
 * @brief Wind at a target height within one atmospheric column
 *
 * Brackets the target height by bisection on the column's terrain-following
 * cell-centre heights, averages the face velocities to cell centres at the two
 * bracketing levels, and interpolates linearly between them. A target below the
 * lowest cell centre or above the highest is clamped to that level rather than
 * extrapolated.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void column_wind_at_height(const Array4<const Real>& xvel,
                           const Array4<const Real>& yvel,
                           const Array4<const Real>& z_cc,
                           int ia, int ja, int nz,
                           Real z_target,
                           Real& u_out, Real& v_out) noexcept
{
    // Bracket by bisection; z_cc is monotonically increasing in k.
    int k_lo;
    if (z_target <= z_cc(ia, ja, 0)) {
        // Below the lowest cell centre. Clamping to level 0 keeps the fire on
        // the near-surface wind; the previous linear scan fell through to the
        // top of the domain here, which is reachable whenever the requested
        // reference height is below half the first cell thickness.
        k_lo = 0;
    } else if (z_target >= z_cc(ia, ja, nz - 1)) {
        k_lo = nz - 2;
    } else {
        int lo = 0;
        int hi = nz - 1;
        while (hi - lo > 1) {
            const int mid = (lo + hi) / 2;
            if (z_cc(ia, ja, mid) <= z_target) { lo = mid; } else { hi = mid; }
        }
        k_lo = lo;
    }
    k_lo = amrex::max(0, amrex::min(k_lo, nz - 2));
    const int k_hi = k_lo + 1;

    const Real z_lo = z_cc(ia, ja, k_lo);
    const Real z_hi = z_cc(ia, ja, k_hi);
    Real alpha = 0.0;
    if (z_hi > z_lo) {
        alpha = (z_target - z_lo) / (z_hi - z_lo);
        alpha = amrex::max(Real(0.0), amrex::min(Real(1.0), alpha));
    }

    const Real u_lo = 0.5 * (xvel(ia, ja, k_lo) + xvel(ia + 1, ja, k_lo));
    const Real v_lo = 0.5 * (yvel(ia, ja, k_lo) + yvel(ia, ja + 1, k_lo));
    const Real u_hi = 0.5 * (xvel(ia, ja, k_hi) + xvel(ia + 1, ja, k_hi));
    const Real v_hi = 0.5 * (yvel(ia, ja, k_hi) + yvel(ia, ja + 1, k_hi));

    u_out = u_lo + alpha * (u_hi - u_lo);
    v_out = v_lo + alpha * (v_hi - v_lo);
}

void fill_fire_wind_from_interpolation(
    MultiFab&       fire_wind_ref,
    MultiFab&       fire_wind_extract_z_mf,
    const MultiFab& xvel_mf,
    const MultiFab& yvel_mf,
    const MultiFab& z_phys_cc_mf,
    const MultiFab& fire_surface_z_mf,
    const MultiFab& fire_col_ground_mf,
    const FireGrid&        fg,
    Real            z_ref,
    int             nz,
    const MultiFab* fuel_model_mf,
    const Real*     d_fcwh,
    int             nfuelcats,
    int             wind_interp)
{
    const int C = fg.C;
    const bool bilinear = (wind_interp == 1);

    for (MFIter mfi(fire_wind_ref, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> fire_wind = fire_wind_ref.array(mfi);
        Array4<Real> extract_z = fire_wind_extract_z_mf.array(mfi);
        Array4<const Real> xvel = xvel_mf.array(mfi);
        Array4<const Real> yvel = yvel_mf.array(mfi);
        Array4<const Real> z_phys_cc = z_phys_cc_mf.array(mfi);
        Array4<const Real> z_surf_arr = fire_surface_z_mf.const_array(mfi);
        Array4<const Real> z_col_arr  = fire_col_ground_mf.const_array(mfi);

        // Phase 13A: Per-fuel wind height support
        Array4<const Real> fuel_model;
        if (fuel_model_mf != nullptr) {
            fuel_model = fuel_model_mf->array(mfi);
        }

        // Bilinear reaches one atmospheric cell beyond the column holding the
        // fire cell, and the face average reaches one further in x or y. Both
        // land in the atmospheric fab's ghost region, which carries two cells
        // for z_phys_cc and more for the velocities; clamp anyway so a fab
        // without them degrades to one-sided rather than reading past the end.
        const Box& zbox = Box(z_phys_cc);
        const int ia_min = zbox.smallEnd(0);
        const int ia_max = zbox.bigEnd(0) - 1;
        const int ja_min = zbox.smallEnd(1);
        const int ja_max = zbox.bigEnd(1) - 1;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            const int i_f = iv[0];
            const int j_f = iv[1];

            // Phase 13A: Determine wind reference height for this cell
            Real z_ref_cell = z_ref;  // default: global fallback
            if (fuel_model_mf != nullptr && d_fcwh != nullptr) {
                // Per-fuel height lookup
                int fuel_code = static_cast<int>(fuel_model(i_f, j_f, 0));
                // Clamp to valid range [1, nfuelcats]
                fuel_code = amrex::max(1, amrex::min(fuel_code, nfuelcats));
                z_ref_cell = d_fcwh[fuel_code];
            }

            if (!bilinear) {
                // Nearest atmospheric column: every fire cell in a column shares
                // one wind vector, so the fire-grid wind is piecewise constant
                // on atmospheric cells.
                const int i_a = amrex::max(ia_min, amrex::min(i_f / C, ia_max));
                const int j_a = amrex::max(ja_min, amrex::min(j_f / C, ja_max));

                // Ground elevation of this fire cell's atmospheric column. This
                // is the terrain surface, not the first cell centre, which sits
                // half an atmospheric cell higher and would push the extraction
                // height up by an amount that varies with vertical resolution.
                const Real z_target = z_surf_arr(i_f, j_f, 0) + z_ref_cell;
                extract_z(i_f, j_f, 0) = z_target;

                Real u_out, v_out;
                column_wind_at_height(xvel, yvel, z_phys_cc, i_a, j_a, nz,
                                      z_target, u_out, v_out);
                fire_wind(i_f, j_f, 0, 0) = u_out;
                fire_wind(i_f, j_f, 0, 1) = v_out;
                return;
            }

            // Bilinear between the four surrounding atmospheric columns. The
            // fire cell centre sits at continuous atmospheric cell-centre
            // coordinate (i_f + 0.5)/C - 0.5.
            const Real gx = (Real(i_f) + 0.5) / Real(C) - 0.5;
            const Real gy = (Real(j_f) + 0.5) / Real(C) - 0.5;
            const int i0 = static_cast<int>(std::floor(gx));
            const int j0 = static_cast<int>(std::floor(gy));
            const Real wx = gx - Real(i0);
            const Real wy = gy - Real(j0);

            Real u_sum = 0.0;
            Real v_sum = 0.0;
            Real z_sum = 0.0;

            for (int c = 0; c < 4; ++c) {
                const int di = (c & 1);
                const int dj = ((c >> 1) & 1);
                const int ia = amrex::max(ia_min, amrex::min(i0 + di, ia_max));
                const int ja = amrex::max(ja_min, amrex::min(j0 + dj, ja_max));

                const Real w = (di ? wx : (1.0 - wx)) * (dj ? wy : (1.0 - wy));

                // Each column is sampled at the same height above *its own*
                // ground: near the surface the profile is anchored to the local
                // terrain, so blending at one absolute height would sample an
                // upslope neighbour far higher above ground than asked for.
                const Real z_target_c = z_col_arr(i_f, j_f, 0, c) + z_ref_cell;

                Real u_c, v_c;
                column_wind_at_height(xvel, yvel, z_phys_cc, ia, ja, nz,
                                      z_target_c, u_c, v_c);

                u_sum += w * u_c;
                v_sum += w * v_c;
                z_sum += w * z_target_c;
            }

            extract_z(i_f, j_f, 0) = z_sum;
            fire_wind(i_f, j_f, 0, 0) = u_sum;
            fire_wind(i_f, j_f, 0, 1) = v_sum;
        });
    }
}
