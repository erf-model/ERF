#include <ERF_FireDustCoupling.H>

#if defined(ERF_ENABLE_FIRE) && defined(ERF_USE_DUST)

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelFor.H>

using namespace amrex;

void FireDustCoupling::apply_burned_area_to_crust(
    MultiFab& dust_crust_index,
    const Geometry& geom_dust) const
{
    if (!enabled || fire_phi_mf == nullptr) {
        return;
    }

    // Iterate over dust grid cells and apply crust reduction where fire has burned
    for (MFIter mfi(dust_crust_index, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> crust = dust_crust_index.array(mfi);
        Array4<const Real> fire_phi = fire_phi_mf->const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            // Check if cell has been burned (fire_phi < 0)
            if (fire_phi(iv) < 0.0_rt) {
                // Apply crust reduction: crust *= (1 - reduction)
                crust(iv) *= (1.0_rt - post_fire_crust_reduction);
                crust(iv) = amrex::max(crust(iv), 0.0_rt);
            }
        });
    }
}

#endif
