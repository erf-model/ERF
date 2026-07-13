/**
 * @file ERF_DustWindExtract.cpp
 * @brief Implementation of wind and surface field extraction from ERF 3D solver to 2D dust grid.
 *
 * Extracts atmospheric wind at reference height and surface fields
 * from the 3D atmospheric solver onto the 2D dust grid each timestep.
 * The wind interpolation algorithm copies fill_fire_wind_from_interpolation
 * from Source/Fire/ERF_FireWindExtract.cpp, with DustGrid substituted for FireGrid.
 */

#include <ERF_DustWindExtract.H>

#ifdef ERF_USE_DUST

#include <AMReX_MFIter.H>

using namespace amrex;

void fill_dust_wind_from_interpolation(
    MultiFab&       dust_wind_ref,
    const MultiFab& xvel_mf,
    const MultiFab& yvel_mf,
    const MultiFab& z_phys_cc_mf,
    const DustGrid& dg,
    Real            zref,
    int             nz)
{
    // Direct vertical interpolation from atmospheric grid to dust grid
    int C = dg.grid_ratio;

    for (MFIter mfi(dust_wind_ref, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> dust_wind = dust_wind_ref.array(mfi);
        Array4<const Real> xvel = xvel_mf.array(mfi);
        Array4<const Real> yvel = yvel_mf.array(mfi);
        Array4<const Real> z_phys_cc = z_phys_cc_mf.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            int i_d = iv[0];  // dust grid index
            int j_d = iv[1];

            // Map to atmospheric column
            int i_a = i_d / C;
            int j_a = j_d / C;

            // Get surface height
            Real z_surf = z_phys_cc(i_a, j_a, 0);

            // Compute target height
            Real z_target = z_surf + zref;

            // Find vertical level bracket.
            // Initialize k_lo to the top interval (nz-2) so that if z_target
            // is above all levels, the topmost wind values are used.
            int k_lo = nz - 2;
            for (int k = 0; k < nz - 1; ++k) {
                if (z_phys_cc(i_a, j_a, k) <= z_target && 
                    z_target < z_phys_cc(i_a, j_a, k + 1)) {
                    k_lo = k;
                    break;
                }
            }

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
            dust_wind(i_d, j_d, 0, 0) = u_cc_lo + alpha * (u_cc_hi - u_cc_lo);
            dust_wind(i_d, j_d, 0, 1) = v_cc_lo + alpha * (v_cc_hi - v_cc_lo);
        });
    }
}

void fill_dust_ustar_from_surface_layer(
    MultiFab&       dust_ustar_in,
    const MultiFab& ustar_atm,
    const DustGrid& dg)
{
    const int C = dg.grid_ratio;
    for (MFIter mfi(dust_ustar_in, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto du = dust_ustar_in.array(mfi);
        auto ua = ustar_atm.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            du(i,j,k) = ua(i/C, j/C, 0);
        });
    }
}

void fill_dust_scalar_from_atm(
    MultiFab&       dust_field,
    const MultiFab& atm_field,
    const DustGrid& dg)
{
    const int C = dg.grid_ratio;
    for (MFIter mfi(dust_field, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto df = dust_field.array(mfi);
        auto af = atm_field.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            df(i,j,k) = af(i/C, j/C, 0);
        });
    }
}

#endif // ERF_USE_DUST
