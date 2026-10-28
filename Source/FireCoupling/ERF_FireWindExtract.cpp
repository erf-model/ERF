#include <ERF_FireWindExtract.H>

using namespace amrex;

/**
 * @brief Wind at a target height within one atmospheric column
 *
 * Brackets the target height by bisection on the column's terrain-following
 * cell-centre heights, averages the face velocities to cell centres at the two
 * bracketing levels, and interpolates linearly between them. A target above the
 * highest cell centre takes that level's wind. A target below the lowest cell
 * centre takes the lowest centre's wind as it is (below_first_log = false, which
 * with a 20 m first cell hands the fire the 10 m wind as its 6.1 m wind, 12 %
 * high at z0 = 0.1 m), or scaled down a neutral log profile of roughness z0_log
 * (below_first_log = true), U(z) = U(z_1) ln(z/z0) / ln(z_1/z0), the better form
 * for a first cell above the reference height. z_ground is the ground height the column heights are
 * measured from; z_target - z_ground is the height above ground.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void column_wind_at_height(const Array4<const Real>& xvel,
                           const Array4<const Real>& yvel,
                           const Array4<const Real>& z_cc,
                           int ia, int ja, int nz,
                           Real z_target, Real z_ground,
                           Real& u_out, Real& v_out,
                           bool below_first_log = false, Real z0_log = Real(0.1)) noexcept
{
    // Bracket by bisection; z_cc is monotonically increasing in k.
    int k_lo;
    // Log-profile factor for a target below the lowest centre, from floored
    // heights (both arguments valid whether or not the branch is taken, so
    // nothing speculated under the fpe traps can raise a flag).
    Real below_scale = Real(1.0);
    if (z_target <= z_cc(ia, ja, 0)) {
        k_lo = 0;
        if (below_first_log) {
            constexpr Real z_tiny = Real(1.0e-30);
            const Real z0  = std::abs(z0_log) + z_tiny;
            const Real z1  = amrex::max(z_cc(ia, ja, 0) - z_ground, z0 * Real(1.000001));
            const Real zt  = amrex::max(amrex::min(z_target - z_ground, z1), z0);
            below_scale = std::log(zt / z0) / std::log(z1 / z0);
        }
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
    Real alpha = Real(0.0);
    if (z_hi > z_lo) {
        alpha = (z_target - z_lo) / (z_hi - z_lo);
        alpha = amrex::max(Real(0.0), amrex::min(Real(1.0), alpha));
    }

    const Real u_lo = Real(0.5) * (xvel(ia, ja, k_lo) + xvel(ia + 1, ja, k_lo));
    const Real v_lo = Real(0.5) * (yvel(ia, ja, k_lo) + yvel(ia, ja + 1, k_lo));
    const Real u_hi = Real(0.5) * (xvel(ia, ja, k_hi) + xvel(ia + 1, ja, k_hi));
    const Real v_hi = Real(0.5) * (yvel(ia, ja, k_hi) + yvel(ia, ja + 1, k_hi));

    u_out = below_scale * (u_lo + alpha * (u_hi - u_lo));
    v_out = below_scale * (v_lo + alpha * (v_hi - v_lo));
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
    const MultiFab* z_ref_cell_mf,
    int             wind_interp,
    const MultiFab* col_open,
    const MultiFab* col_roof,
    Real            z_sample,
    Real            z0_log,
    bool            below_first_log)
{
    const int C = fg.C;
    // Fire cell (i_f, j_f) lies in atmospheric cell (i_f / C + lo_x, j_f / C + lo_y)
    const int lo_x = fg.atm_lo[0];
    const int lo_y = fg.atm_lo[1];
    const bool bilinear = (wind_interp == 1);
    const bool open_weights = bilinear && (col_open != nullptr) && (col_roof != nullptr);
    const bool resample = (z_sample > Real(0.0));

    for (MFIter mfi(fire_wind_ref, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> fire_wind = fire_wind_ref.array(mfi);
        Array4<Real> extract_z = fire_wind_extract_z_mf.array(mfi);
        Array4<const Real> xvel = xvel_mf.array(mfi);
        Array4<const Real> yvel = yvel_mf.array(mfi);
        Array4<const Real> z_phys_cc = z_phys_cc_mf.array(mfi);
        Array4<const Real> z_surf_arr = fire_surface_z_mf.const_array(mfi);
        Array4<const Real> z_col_arr  = fire_col_ground_mf.const_array(mfi);

        // Per-cell sampling height, when the engine supplies one
        Array4<const Real> z_ref_cell_arr;
        if (z_ref_cell_mf != nullptr) {
            z_ref_cell_arr = z_ref_cell_mf->const_array(mfi);
        }
        const bool per_cell_ht = (z_ref_cell_mf != nullptr);

        // Open-column weights: the 2D fields carry one ghost cell, which is
        // all the bilinear stencil reaches; clamp to their box regardless.
        Array4<const Real> open_arr, roof_arr;
        int io_min = 0, io_max = 0, jo_min = 0, jo_max = 0;
        if (open_weights) {
            open_arr = col_open->const_array(mfi);
            roof_arr = col_roof->const_array(mfi);
            const Box& obox = Box(open_arr);
            io_min = obox.smallEnd(0); io_max = obox.bigEnd(0);
            jo_min = obox.smallEnd(1); jo_max = obox.bigEnd(1);
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

            // Wind reference height for this cell: its own, or the global one
            Real z_ref_cell = z_ref;
            if (per_cell_ht) { z_ref_cell = z_ref_cell_arr(i_f, j_f, 0); }

            // Height actually sampled and the log-law factor that brings the
            // sampled wind down to the target height (1 without resampling).
            // The log ratio is written so that every argument is valid whether
            // or not it is selected: clang evaluates the ternary's arm
            // speculatively in Release, and with resampling off z_sample is 0,
            // so log(0) raised the divide-by-zero flag and tripped
            // amrex.fpe_trap_zero on the first step. abs and the tiny offset
            // are not selects, so nothing is folded back into the arm; with
            // resampling on the values are unchanged.
            constexpr Real z_tiny  = Real(1.0e-30);
            const Real z0_safe     = std::abs(z0_log) + z_tiny;
            const Real log_ref     = std::log((std::abs(z_ref_cell) + z_tiny) / z0_safe);
            const Real log_sample  = std::log((std::abs(z_sample)   + z_tiny) / z0_safe);
            const Real z_at  = resample ? z_sample : z_ref_cell;
            const Real scale = resample ? log_ref / (log_sample + z_tiny) : Real(1.0);

            if (!bilinear) {
                // Nearest atmospheric column: every fire cell in a column shares
                // one wind vector, so the fire-grid wind is piecewise constant
                // on atmospheric cells.
                const int i_a = amrex::max(ia_min, amrex::min(i_f / C + lo_x, ia_max));
                const int j_a = amrex::max(ja_min, amrex::min(j_f / C + lo_y, ja_max));

                // Ground elevation of this fire cell's atmospheric column. This
                // is the terrain surface, not the first cell centre, which sits
                // half an atmospheric cell higher and would push the extraction
                // height up by an amount that varies with vertical resolution.
                const Real z_target = z_surf_arr(i_f, j_f, 0) + z_at;
                extract_z(i_f, j_f, 0) = z_target;

                Real u_out, v_out;
                column_wind_at_height(xvel, yvel, z_phys_cc, i_a, j_a, nz,
                                      z_target, z_surf_arr(i_f, j_f, 0), u_out, v_out,
                                      below_first_log, z0_log);
                fire_wind(i_f, j_f, 0, 0) = scale * u_out;
                fire_wind(i_f, j_f, 0, 1) = scale * v_out;
                return;
            }

            // Bilinear between the four surrounding atmospheric columns. The
            // fire cell centre sits at continuous atmospheric cell-centre
            // coordinate (i_f + 0.5)/C - 0.5 relative to the region's corner; the
            // weights come from that, the column indices add the corner.
            const Real gx = (Real(i_f) + Real(0.5)) / Real(C) - Real(0.5);
            const Real gy = (Real(j_f) + Real(0.5)) / Real(C) - Real(0.5);
            const int i0_local = static_cast<int>(std::floor(gx));
            const int j0_local = static_cast<int>(std::floor(gy));
            const Real wx = gx - Real(i0_local);
            const Real wy = gy - Real(j0_local);
            const int i0 = i0_local + lo_x;
            const int j0 = j0_local + lo_y;

            Real u_sum = Real(0.0);
            Real v_sum = Real(0.0);
            Real z_sum = Real(0.0);

            // Open-column factors and their bilinear sum; fall back to the
            // plain weights when every column is closed.
            Real f_open[4] = {Real(1.0), Real(1.0), Real(1.0), Real(1.0)};
            Real w_norm = Real(1.0);
            if (open_weights) {
                Real wsum = Real(0.0);
                for (int c = 0; c < 4; ++c) {
                    const int di = (c & 1);
                    const int dj = ((c >> 1) & 1);
                    const int io = amrex::max(io_min, amrex::min(i0 + di, io_max));
                    const int jo = amrex::max(jo_min, amrex::min(j0 + dj, jo_max));
                    const Real w = (di ? wx : (Real(1.0) - wx)) * (dj ? wy : (Real(1.0) - wy));
                    f_open[c] = (roof_arr(io, jo, 0) > z_at) ? open_arr(io, jo, 0) : Real(1.0);
                    wsum += w * f_open[c];
                }
                if (wsum > Real(1.0e-6)) {
                    w_norm = Real(1.0) / wsum;
                } else {
                    for (int c = 0; c < 4; ++c) { f_open[c] = Real(1.0); }
                }
            }

            for (int c = 0; c < 4; ++c) {
                const int di = (c & 1);
                const int dj = ((c >> 1) & 1);
                const int ia = amrex::max(ia_min, amrex::min(i0 + di, ia_max));
                const int ja = amrex::max(ja_min, amrex::min(j0 + dj, ja_max));

                const Real w = (di ? wx : (Real(1.0) - wx)) * (dj ? wy : (Real(1.0) - wy))
                             * f_open[c] * w_norm;

                // Each column is sampled at the same height above *its own*
                // ground: near the surface the profile is anchored to the local
                // terrain, so blending at one absolute height would sample an
                // upslope neighbour far higher above ground than asked for.
                const Real z_target_c = z_col_arr(i_f, j_f, 0, c) + z_at;

                Real u_c, v_c;
                column_wind_at_height(xvel, yvel, z_phys_cc, ia, ja, nz,
                                      z_target_c, z_col_arr(i_f, j_f, 0, c), u_c, v_c,
                                      below_first_log, z0_log);

                u_sum += w * u_c;
                v_sum += w * v_c;
                z_sum += w * z_target_c;
            }

            extract_z(i_f, j_f, 0) = z_sum;
            fire_wind(i_f, j_f, 0, 0) = scale * u_sum;
            fire_wind(i_f, j_f, 0, 1) = scale * v_sum;
        });
    }
}

void extract_fire_wind_k0(
    amrex::MultiFab&       u_k0,
    amrex::MultiFab&       v_k0,
    amrex::MultiFab&       h_k0,
    const amrex::MultiFab& xvel,
    const amrex::MultiFab& yvel,
    const amrex::MultiFab& z_phys_cc,
    const amrex::MultiFab* z_phys_nd,
    const amrex::Geometry& geom)
{
    using namespace amrex;
    const int  k0      = geom.Domain().smallEnd(2);
    const Real z_floor = geom.ProbLo(2);
    const bool use_znd = (z_phys_nd != nullptr);
    for (MFIter mfi(u_k0); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto u   = xvel.const_array(mfi);
        const auto v   = yvel.const_array(mfi);
        const auto zcc = z_phys_cc.const_array(mfi);
        const Array4<Real const> znd = use_znd ? z_phys_nd->const_array(mfi) : Array4<Real const>{};
        const auto us = u_k0.array(mfi);
        const auto vs = v_k0.array(mfi);
        const auto hs = h_k0.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            us(i,j,0) = Real(0.5) * (u(i,j,k0) + u(i+1,j,k0));
            vs(i,j,0) = Real(0.5) * (v(i,j,k0) + v(i,j+1,k0));
            const Real zsfc = use_znd ? Real(0.25) * (znd(i,j,k0) + znd(i+1,j,k0) + znd(i,j+1,k0) + znd(i+1,j+1,k0))
                                      : z_floor;
            hs(i,j,0) = zcc(i,j,k0) - zsfc;
        });
    }
    u_k0.FillBoundary(geom.periodicity());
    v_k0.FillBoundary(geom.periodicity());
    h_k0.FillBoundary(geom.periodicity());
}
