#include <ERF_FireDustCoupling.H>
#include <ERF_HostFabView.H>

#if defined(ERF_ENABLE_FIRE) && defined(ERF_USE_DUST)

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>
#include <AMReX_Loop.H>

using namespace amrex;

void FireDustCoupling::apply_burned_area_to_crust(
    MultiFab&       dust_crust_index,
    const Geometry& geom_dust) const
{
    if (!enabled || fire_phi_mf == nullptr) {
        return;
    }

    // -----------------------------------------------------------------------
    // Build a host-side flat array of fire phi values so we can look up
    // any fire cell by (i_f, j_f) index without needing MFIter alignment.
    // This is safe because the fire sub-grid is small (e.g. 16x16x1).
    // -----------------------------------------------------------------------
    const Box& fire_domain = geom_fire.Domain();
    const int fi_lo = fire_domain.smallEnd(0);
    const int fj_lo = fire_domain.smallEnd(1);
    const int fi_hi = fire_domain.bigEnd(0);
    const int fj_hi = fire_domain.bigEnd(1);
    const int fnx   = fi_hi - fi_lo + 1;
    const int fny   = fj_hi - fj_lo + 1;

    // Copy fire phi into a CPU vector indexed as phi_host[j * fnx + i]
    amrex::Vector<amrex::Real> phi_host(fnx * fny, 0.0_rt);

    for (MFIter mfi_f(*fire_phi_mf); mfi_f.isValid(); ++mfi_f) {
        const Box& bx_f = mfi_f.validbox();
        // amrex::Loop runs on the host, and under CUDA the fire level set lives
        // in device memory, so stage the FAB before reading it. Gathering into a
        // flat host vector is the right shape here: the result feeds an MPI
        // reduce so every rank ends up with the whole fire field.
        const ERFHostFabView phi_view((*fire_phi_mf)[mfi_f]);
        auto phi_arr = phi_view.array();
        amrex::Loop(bx_f, [&](int i, int j, int /*k*/) {
            int ii = i - fi_lo;
            int jj = j - fj_lo;
            if (ii >= 0 && ii < fnx && jj >= 0 && jj < fny) {
                phi_host[jj * fnx + ii] = phi_arr(i, j, 0);
            }
        });
    }
    // Reduce across MPI ranks so every rank has the full fire phi array
    amrex::ParallelDescriptor::ReduceRealSum(
        phi_host.data(), static_cast<int>(phi_host.size()));

    // -----------------------------------------------------------------------
    // For each dust cell, map its physical centre to a fire cell index and
    // check whether fire_phi < 0 (burned). Apply crust reduction if so.
    // -----------------------------------------------------------------------
    const auto& fire_plo  = geom_fire.ProbLoArray();
    const auto  fire_dx   = geom_fire.CellSizeArray();
    const auto& dust_plo  = geom_dust.ProbLoArray();
    const auto  dust_dx   = geom_dust.CellSizeArray();

    const amrex::Real reduction = post_fire_crust_reduction;

    // Copy fire phi into a GPU-resident vector for device access
    amrex::Gpu::DeviceVector<amrex::Real> phi_device(fnx * fny);
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, phi_host.begin(), phi_host.end(), phi_device.begin());

    // Capture geometry data as integers to avoid floating-point issues
    const int fi_lo_int = fi_lo;
    const int fj_lo_int = fj_lo;
    const int fi_hi_int = fi_hi;
    const int fj_hi_int = fj_hi;
    const int fnx_int = fnx;

    for (MFIter mfi(dust_crust_index, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto crust = dust_crust_index.array(mfi);

        // Get raw device pointer for GPU access
        const amrex::Real* phi_ptr = phi_device.data();
        const amrex::Real fire_plo_0 = fire_plo[0];
        const amrex::Real fire_plo_1 = fire_plo[1];
        const amrex::Real fire_dx_0 = fire_dx[0];
        const amrex::Real fire_dx_1 = fire_dx[1];
        const amrex::Real dust_plo_0 = dust_plo[0];
        const amrex::Real dust_plo_1 = dust_plo[1];
        const amrex::Real dust_dx_0 = dust_dx[0];
        const amrex::Real dust_dx_1 = dust_dx[1];

        // Use ParallelFor for GPU/CPU compatibility
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Physical centre of this dust cell
            amrex::Real xc = dust_plo_0 + (i + 0.5_rt) * dust_dx_0;
            amrex::Real yc = dust_plo_1 + (j + 0.5_rt) * dust_dx_1;

            // Corresponding fire cell index (nearest cell containing xc, yc)
            int i_f = static_cast<int>((xc - fire_plo_0) / fire_dx_0);
            int j_f = static_cast<int>((yc - fire_plo_1) / fire_dx_1);

            // Clamp to valid fire domain
            i_f = amrex::max(fi_lo_int, amrex::min(fi_hi_int, i_f + fi_lo_int));
            j_f = amrex::max(fj_lo_int, amrex::min(fj_hi_int, j_f + fj_lo_int));

            int ii = i_f - fi_lo_int;
            int jj = j_f - fj_lo_int;

            // Apply crust reduction if fire_phi < 0 (burned)
            if (phi_ptr[jj * fnx_int + ii] < 0.0_rt) {
                crust(i, j, k) *= (1.0_rt - reduction);
                crust(i, j, k)  = amrex::max(crust(i, j, k), 0.0_rt);
            }
        });
    }

    // ParallelFor is asynchronous; phi_device frees its device allocation at the
    // end of this function. Let the kernels finish reading it.
    amrex::Gpu::streamSynchronize();

    // Report debug info
    amrex::Real crust_min = dust_crust_index.min(0);
    amrex::Real crust_max = dust_crust_index.max(0);
    amrex::Print() << "[DUST DEBUG] Fire-dust coupling: modified crust values "
                   << "crust_min=" << crust_min << ", crust_max=" << crust_max
                   << ", reduction=" << reduction << "\n";
}

void FireDustCoupling::apply_fire_wind_to_dust_ustar(
    amrex::MultiFab&       dust_ustar_in,
    const amrex::MultiFab& fire_wind_scratch,
    const amrex::Geometry& /*geom_dust*/,
    amrex::Real            z0,
    amrex::Real            zref,
    int                    C) const
{
    if (!enabled || !fire_wind_to_dust) return;

    // Clamp C to at least 1 (if fire and dust grids have same resolution)
    const int C_ratio = amrex::max(C, 1);

    // Log-law constant: u* = U * kappa / ln(zref/z0)
    constexpr amrex::Real kappa = 0.4_rt;
    const amrex::Real log_ratio = (zref > z0 && z0 > 0.0_rt)
        ? std::log(zref / z0) : 1.0_rt;

    for (amrex::MFIter mfi(dust_ustar_in, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto ustar = dust_ustar_in.array(mfi);
        auto wind  = fire_wind_scratch.const_array(mfi);

        amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                // Average fire wind over the C×C fire cells that cover
                // dust cell (i, j). When C=1 this is a direct cell lookup.
                amrex::Real u_avg = 0.0_rt, v_avg = 0.0_rt;
                for (int di = 0; di < C_ratio; ++di) {
                    for (int dj = 0; dj < C_ratio; ++dj) {
                        u_avg += wind(i*C_ratio + di, j*C_ratio + dj, 0, 0);
                        v_avg += wind(i*C_ratio + di, j*C_ratio + dj, 0, 1);
                    }
                }
                const amrex::Real inv = 1.0_rt / amrex::Real(C_ratio * C_ratio);
                u_avg *= inv;
                v_avg *= inv;

                // Derive u* from fire wind speed via log-law
                const amrex::Real spd        = std::sqrt(u_avg*u_avg + v_avg*v_avg);
                const amrex::Real ustar_fire = spd * kappa / log_ratio;

                // Dust emission driven by the larger of MRF u* and fire u*
                ustar(i, j, k) = amrex::max(ustar(i, j, k), ustar_fire);
            });
    }

    // Debug output — after all MPI collectives (LNG_MPI_SKILLS Rule B1)
    amrex::Real ustar_max = dust_ustar_in.max(0);
    amrex::Real wind_max  = fire_wind_scratch.max(0);
    amrex::Print() << "[DUST DEBUG] Phase 2 fire-wind coupling:"
                   << " fire_wind_max=" << wind_max << " m/s"
                   << " ustar_max_after=" << ustar_max << " m/s"
                   << " C=" << C_ratio << "\n";
}

#endif
