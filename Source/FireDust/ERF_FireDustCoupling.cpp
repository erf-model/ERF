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

    for (MFIter mfi(dust_crust_index, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto crust = dust_crust_index.array(mfi);

        // phi_host lives on CPU; loop on CPU via amrex::Loop (host path)
        amrex::Loop(bx, [&](int i, int j, int k) noexcept {
            // Physical centre of this dust cell
            amrex::Real xc = dust_plo[0] + (i + 0.5_rt) * dust_dx[0];
            amrex::Real yc = dust_plo[1] + (j + 0.5_rt) * dust_dx[1];

            // Corresponding fire cell index (nearest cell containing xc, yc)
            int i_f = static_cast<int>((xc - fire_plo[0]) / fire_dx[0]);
            int j_f = static_cast<int>((yc - fire_plo[1]) / fire_dx[1]);

            // Clamp to valid fire domain
            i_f = amrex::max(fi_lo, amrex::min(fi_hi, i_f + fi_lo));
            j_f = amrex::max(fj_lo, amrex::min(fj_hi, j_f + fj_lo));

            int ii = i_f - fi_lo;
            int jj = j_f - fj_lo;

            // Apply crust reduction if fire_phi < 0 (burned)
            if (phi_host[jj * fnx + ii] < 0.0_rt) {
                crust(i, j, k) *= (1.0_rt - reduction);
                crust(i, j, k)  = amrex::max(crust(i, j, k), 0.0_rt);
            }
        });
    }
}

#endif
