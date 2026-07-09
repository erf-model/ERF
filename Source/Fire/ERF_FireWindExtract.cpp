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

            // Classify terrain position and compute speed factor
            Real factor = 1.0;

            // Ridge: convex, windward slope
            if (curv_val > 0.01 && slope_mag > 0.05 && cos_upslope > 0.0) {
                factor = 1.0 + (k_ridge - 1.0) * amrex::min(slope_mag, 1.0_rt);
            }
            // Shelter: lee-side slope (wind blowing downslope)
            else if (cos_upslope < -0.5 && curv_val > 0.01) {
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

void fill_fire_wind_from_interpolation(
    MultiFab&       fire_wind_ref,
    MultiFab&       fire_wind_extract_z_mf,
    const MultiFab& xvel_mf,
    const MultiFab& yvel_mf,
    const MultiFab& z_phys_cc_mf,
    const FireGrid&        fg,
    Real            z_ref,
    int             nz,
    const MultiFab* fuel_model_mf,
    const Real*     d_fcwh,
    int             nfuelcats)
{
    // Direct vertical interpolation from atmospheric grid to fire grid
    int C = fg.C;

    for (MFIter mfi(fire_wind_ref, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> fire_wind = fire_wind_ref.array(mfi);
        Array4<Real> extract_z = fire_wind_extract_z_mf.array(mfi);
        Array4<const Real> xvel = xvel_mf.array(mfi);
        Array4<const Real> yvel = yvel_mf.array(mfi);
        Array4<const Real> z_phys_cc = z_phys_cc_mf.array(mfi);
        
        // Phase 13A: Per-fuel wind height support
        Array4<const Real> fuel_model;
        if (fuel_model_mf != nullptr) {
            fuel_model = fuel_model_mf->array(mfi);
        }

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            int i_f = iv[0];
            int j_f = iv[1];

            // Map to atmospheric column
            int i_a = i_f / C;
            int j_a = j_f / C;

            // Get surface height
            Real z_surf = z_phys_cc(i_a, j_a, 0);
            
            // Phase 13A: Determine wind reference height for this cell
            Real z_ref_cell = z_ref;  // default: global fallback
            if (fuel_model_mf != nullptr && d_fcwh != nullptr) {
                // Per-fuel height lookup
                int fuel_code = static_cast<int>(fuel_model(i_f, j_f, 0));
                // Clamp to valid range [1, nfuelcats]
                fuel_code = amrex::max(1, amrex::min(fuel_code, nfuelcats));
                z_ref_cell = d_fcwh[fuel_code];
            }
            
            // Compute target height
            Real z_target = z_surf + z_ref_cell;

            // Store the extraction height for plotfile output
            extract_z(i_f, j_f, 0) = z_target;

            // Find vertical level bracket.
            // Initialise k_lo to the top interval (nz-2) so that if z_target
            // is above all levels, the topmost wind values are used.
            int k_lo = nz - 2;
            for (int k = 0; k < nz - 1; ++k) {
                if (z_phys_cc(i_a, j_a, k) <= z_target && 
                    z_target < z_phys_cc(i_a, j_a, k + 1)) {
                    k_lo = k;
                    break;
                }
            }
            // k_lo is already guaranteed to be in [0, nz-2]; no additional clamp needed.

            // Compute interpolation weight
            Real z_lo = z_phys_cc(i_a, j_a, k_lo);
            Real z_hi = z_phys_cc(i_a, j_a, k_lo + 1);
            Real alpha = 0.0;
            if (z_hi > z_lo) {
                alpha = (z_target - z_lo) / (z_hi - z_lo);
                alpha = amrex::max(0.0, amrex::min(1.0, alpha));
            }

            int k_hi = k_lo + 1;

            // Average u/v from faces to cell centers
            Real u_cc_lo = 0.5 * (xvel(i_a, j_a, k_lo) + xvel(i_a + 1, j_a, k_lo));
            Real v_cc_lo = 0.5 * (yvel(i_a, j_a, k_lo) + yvel(i_a, j_a + 1, k_lo));

            Real u_cc_hi = 0.5 * (xvel(i_a, j_a, k_hi) + xvel(i_a + 1, j_a, k_hi));
            Real v_cc_hi = 0.5 * (yvel(i_a, j_a, k_hi) + yvel(i_a, j_a + 1, k_hi));

            // Interpolate to target height
            fire_wind(i_f, j_f, 0, 0) = u_cc_lo + alpha * (u_cc_hi - u_cc_lo);
            fire_wind(i_f, j_f, 0, 1) = v_cc_lo + alpha * (v_cc_hi - v_cc_lo);
        });
    }
}
