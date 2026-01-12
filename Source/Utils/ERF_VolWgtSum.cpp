#include <iomanip>

#include "ERF.H"

using namespace amrex;

/**
 * Utility function for computing a volume weighted sum of MultiFab data for a single component
 *
 * @param lev Current level
 * @param mf_to_be_summed : MultiFab on which we do the volume weighted sum
 * @param dJ       : volume weighting due to metric terms
 * @param mfmx     : map factor in x-direction at cell centers
 * @param mfmy     : map factor in y-direction at cell centers
 * @param comp     : Index of the component we want to sum
 * @param finemask : If a finer level is available, determines whether we mask fine data
 * @param local    : Boolean sets whether or not to reduce the sum over the domain (false) or compute sums local to each MPI rank (true)
 */
Real
ERF::volWgtSumMF (int lev,
                  const MultiFab& mf_to_be_summed, int comp,
                  const MultiFab& dJ, const MultiFab& mfmx, const MultiFab& mfmy,
                  bool finemask,
                  bool local)
{
    BL_PROFILE("ERF::volWgtSumMF()");

    Real sum = 0.0;
    MultiFab tmp(mf_to_be_summed.boxArray(), mf_to_be_summed.DistributionMap(), 1, 0);

    // The quantity that is conserved is not (rho S), but rather (rho S / m^2) where
    // m is the map scale factor at cell centers
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx   = mfi.tilebox();
        const auto  dst_arr = tmp.array(mfi);
        const auto  src_arr = mf_to_be_summed.array(mfi);
        const auto& mfx_arr = mfmx.const_array(mfi);
        const auto& mfy_arr = mfmy.const_array(mfi);

        if (SolverChoice::terrain_type != TerrainType::EB) {
            if (SolverChoice::mesh_type == MeshType::ConstantDz) {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    dst_arr(i,j,k,0) = src_arr(i,j,k,comp) / (mfx_arr(i,j,0)*mfy_arr(i,j,0));
                });
            } else {
                const auto&  dJ_arr = dJ.const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    dst_arr(i,j,k,0) = src_arr(i,j,k,comp) * dJ_arr(i,j,k) / (mfx_arr(i,j,0)*mfy_arr(i,j,0));
                });
            }
        } else {
            const auto&  dJ_arr = dJ.const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                dst_arr(i,j,k,0) = src_arr(i,j,k,comp) * dJ_arr(i,j,k);
            });
        }

    } // mfi

    if (lev < finest_level && finemask) {
        MultiFab::Multiply(tmp, *fine_mask[lev+1].get(), 0, 0, 1, 0);
    }

    // If local = true  then "sum" will be the sum only over the FABs on each rank
    // If local = false then "sum" will be the sum over the whole MultiFab, and will be broadcast to all ranks
    sum = tmp.sum(0,local);

    auto const& dx = geom[lev].CellSizeArray();

    sum *= dx[0]*dx[1]*dx[2];

    return sum;
}

void
ERF::volWgtColumnSum (int lev, const MultiFab& mf_to_be_summed, int comp,
                      MultiFab& mf_2d, const MultiFab& dJ)
{
    BL_PROFILE("ERF::volWgtSumColumnMF()");

    mf_2d.setVal(0.);

    // The quantity that is conserved is not (rho S), but rather (rho S / m^2) where
    // m is the map scale factor at cell centers
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf_to_be_summed, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx   = mfi.tilebox();
        const auto  dst_arr = mf_2d.array(mfi);
        const auto  src_arr = mf_to_be_summed.array(mfi);
        if (SolverChoice::mesh_type == MeshType::ConstantDz) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                amrex::HostDevice::Atomic::Add(&dst_arr(i,j,0),src_arr(i,j,k,comp));
            });
        } else {
            const auto& dJ_arr = dJ.const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                amrex::HostDevice::Atomic::Add(&dst_arr(i,j,0),src_arr(i,j,k,comp)*dJ_arr(i,j,k));
            });
        }
    } // mfi

    auto const& dx = geom[lev].CellSizeArray();

    mf_2d.mult(dx[2]);
}

/**
 * Helper function for constructing a fine mask, that is, a MultiFab
 * masking coarser data at a lower level by zeroing out covered cells
 * in the fine mask MultiFab we compute.
 *
 * @param level Fine level index which masks underlying coarser data
 */
void
ERF::build_fine_mask (int level, MultiFab& fine_mask_lev)
{
    // Mask for zeroing covered cells
    AMREX_ASSERT(level > 0);

    BoxArray cba            = grids[level-1];
    DistributionMapping cdm = dmap[level-1];

    BoxArray fba            = fine_mask_lev.boxArray();

    iMultiFab ifine_mask_lev = makeFineMask(cba, cdm, fba, ref_ratio[level-1], 1, 0);

    const auto  fma =  fine_mask_lev.arrays();
    const auto ifma = ifine_mask_lev.arrays();
    ParallelFor(fine_mask_lev, [=] AMREX_GPU_DEVICE(int bno, int i, int j, int k) noexcept
    {
       fma[bno](i,j,k) = ifma[bno](i,j,k);
    });
}
