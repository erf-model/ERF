#include <ERF_FireWindExtract.H>


using namespace amrex;

void fill_fire_wind_from_most(
    MultiFab& fire_wind_ref,
    const MultiFab& ustar_mf,
    const MultiFab& z0_mf,
    const MultiFab& olen_mf,
    const MultiFab& uavg_mf,
    const MultiFab& vavg_mf,
    const FireGrid& fg,
    Real z_target)
{
    // Per-cell kernel: map atmospheric wind to fire grid
    int C = fg.C;

    for (MFIter mfi(fire_wind_ref, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> fire_wind = fire_wind_ref.array(mfi);
        Array4<const Real> ustar = ustar_mf.array(mfi);
        Array4<const Real> z0 = z0_mf.array(mfi);
        Array4<const Real> olen = olen_mf.array(mfi);
        Array4<const Real> uavg = uavg_mf.array(mfi);
        Array4<const Real> vavg = vavg_mf.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            int i_f = iv[0];
            int j_f = iv[1];

            // Map to atmospheric column via integer division
            int i_a = i_f / C;
            int j_a = j_f / C;

            // Read MOST diagnostics
            Real u_star = ustar(i_a, j_a, 0);
            Real z0_val = z0(i_a, j_a, 0);
            Real olen_val = olen(i_a, j_a, 0);

            // Reconstruct wind at z_target
            Real U_ref = most_wind_at_height(u_star, z0_val, olen_val, z_target);

            // Read wind direction from averages
            Real u_avg = uavg(i_a, j_a, 0);
            Real v_avg = vavg(i_a, j_a, 0);
            Real wind_mag = std::sqrt(u_avg*u_avg + v_avg*v_avg);

            // Set fire grid wind
            if (wind_mag > 1.0e-6) {
                fire_wind(i_f, j_f, 0, 0) = U_ref * (u_avg / wind_mag);
                fire_wind(i_f, j_f, 0, 1) = U_ref * (v_avg / wind_mag);
            } else {
                fire_wind(i_f, j_f, 0, 0) = 0.0;
                fire_wind(i_f, j_f, 0, 1) = 0.0;
            }
        });
    }
}

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
