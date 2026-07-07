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
    MultiFab& fire_wind,
    const MultiFab& fire_slopes,
    const MultiFab& curvature,
    Real k_ridge, Real k_shelter,
    Real k_valley, Real k_deflect,
    Real min_curv)
{
    // Apply terrain corrections to wind field in-place
    // Ridge speed-up, lee sheltering, valley channeling, deflection

    for (MFIter mfi(fire_wind, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> wind = fire_wind.array(mfi);
        Array4<const Real> slopes = fire_slopes.array(mfi);
        Array4<const Real> curv = curvature.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            int i = iv[0];
            int j = iv[1];
            int k = 0;

            Real u = wind(i, j, k, 0);
            Real v = wind(i, j, k, 1);
            Real wind_mag = std::sqrt(u*u + v*v);

            if (wind_mag < 1.0e-6) {
                // No wind, no correction
                return;
            }

            // Terrain properties
            Real dzdx = slopes(i, j, k, 0);
            Real dzdy = slopes(i, j, k, 1);
            Real slope_mag = std::sqrt(dzdx*dzdx + dzdy*dzdy);
            Real terrain_curv = curv(i, j, k);

            // Ridge/valley detection
            Real speed_factor = 1.0;
            if (std::abs(terrain_curv) > min_curv) {
                if (terrain_curv > 0.0) {
                    // Ridge (convex)
                    speed_factor = 1.0 + (k_ridge - 1.0) * amrex::min(1.0, terrain_curv / 0.1);
                } else {
                    // Valley (concave)
                    speed_factor = 1.0 - (1.0 - k_valley) * amrex::min(1.0, std::abs(terrain_curv) / 0.1);
                }
            }

            // Sheltering on lee side
            // Wind direction
            Real wind_dir = std::atan2(v, u);
            // Slope aspect (direction of maximum slope)
            Real slope_aspect = std::atan2(-dzdy, -dzdx);
            // Angle between wind and slope
            Real angle_diff = wind_dir - slope_aspect;

            // Normalize to [-π, π]
            const Real pi = 3.14159265358979323846;
            while (angle_diff > pi) angle_diff -= 2.0*pi;
            while (angle_diff < -pi) angle_diff += 2.0*pi;

            // Lee side sheltering: reduce wind on downwind side
            Real cos_angle = std::cos(angle_diff);
            if (cos_angle < -0.5) {
                // Lee side
                speed_factor *= k_shelter;
            }

            // Deflection: wind turns on slopes (up to ±45°)
            Real max_deflect = 45.0 * pi / 180.0;  // Convert to radians
            Real deflect_rad = k_deflect * std::sin(angle_diff) * slope_mag;
            deflect_rad = amrex::max(-max_deflect, amrex::min(max_deflect, deflect_rad));

            // Apply speed factor
            Real new_mag = wind_mag * speed_factor;

            // Apply deflection by rotating wind
            Real new_dir = wind_dir + deflect_rad;
            wind(i, j, k, 0) = new_mag * std::cos(new_dir);
            wind(i, j, k, 1) = new_mag * std::sin(new_dir);
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
    int             nz)
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

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            int i_f = iv[0];
            int j_f = iv[1];

            // Map to atmospheric column
            int i_a = i_f / C;
            int j_a = j_f / C;

            // Get surface height and compute target height
            Real z_surf = z_phys_cc(i_a, j_a, 0);
            Real z_target = z_surf + z_ref;

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
