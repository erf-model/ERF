#include <ERF_FireDustCoupling.H>

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
        auto phi_arr = fire_phi_mf->const_array(mfi_f);
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

    // Track stats for debug output
    int n_burned_cells = 0;

    // Copy fire phi into a GPU-resident vector for device access
    amrex::Gpu::DeviceVector<amrex::Real> phi_device(fnx * fny);
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, phi_host.begin(), phi_host.end(), phi_device.begin());

    for (MFIter mfi(dust_crust_index, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto crust = dust_crust_index.array(mfi);

        // Get raw device pointer for GPU access
        const amrex::Real* phi_ptr = phi_device.data();
        const amrex::Real fi_lo_val = fi_lo;
        const amrex::Real fj_lo_val = fj_lo;
        const amrex::Real fnx_val = fnx;
        const amrex::Real fny_val = fny;
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
            i_f = amrex::max(static_cast<int>(fi_lo_val), amrex::min(static_cast<int>(fi_lo_val) + static_cast<int>(fnx_val) - 1, i_f + static_cast<int>(fi_lo_val)));
            j_f = amrex::max(static_cast<int>(fj_lo_val), amrex::min(static_cast<int>(fj_lo_val) + static_cast<int>(fny_val) - 1, j_f + static_cast<int>(fj_lo_val)));

            int ii = i_f - static_cast<int>(fi_lo_val);
            int jj = j_f - static_cast<int>(fj_lo_val);

            // Apply crust reduction if fire_phi < 0 (burned)
            if (phi_ptr[static_cast<int>(jj * fnx_val + ii)] < 0.0_rt) {
                crust(i, j, k) *= (1.0_rt - reduction);
                crust(i, j, k)  = amrex::max(crust(i, j, k), 0.0_rt);
            }
        });
    }

    // Report debug info
    amrex::Real crust_min = dust_crust_index.min(0);
    amrex::Real crust_max = dust_crust_index.max(0);
    amrex::Print() << "[DUST DEBUG] Fire-dust coupling: modified crust values "
                   << "crust_min=" << crust_min << ", crust_max=" << crust_max
                   << ", reduction=" << reduction << "\n";
}

#endif
