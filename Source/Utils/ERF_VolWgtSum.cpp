#include <iomanip>

#include "ERF.H"

using namespace amrex;

/**
 * Utility function for computing a volume weighted sum of MultiFab data for a single component
 *
 * @param lev Current level
 * @param mf MultiFab on which we do the volume weighted sum
 * @param comp Index of the component we want to sum
 * @param local Boolean sets whether or not to reduce the sum over the domain (false) or compute sums local to each MPI rank (true)
 * @param finemask If a finer level is available, determines whether we mask fine data
 */
Real
ERF::volWgtSumMF (int lev,
                  const MultiFab& mf,
                  int comp,
                  bool finemask)
{
    BL_PROFILE("ERF::volWgtSumMF()");

    Real sum = 0.0;
    MultiFab tmp(grids[lev], dmap[lev], 1, 0);
    MultiFab::Copy(tmp, mf, comp, 0, 1, 0);

    // The quantity that is conserved is not (rho S), but rather (rho S / m^2) where
    // m is the map scale factor at cell centers
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx   = mfi.tilebox();
        const auto  dst = tmp.array(mfi);
        const auto& mfx = mapfac[lev][MapFacType::m_x]->const_array(mfi);
        const auto& mfy = mapfac[lev][MapFacType::m_y]->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            dst(i,j,k) /= (mfx(i,j,0)*mfy(i,j,0));
        });
    } // mfi

    if (lev < finest_level && finemask) {
        const MultiFab& mask = build_fine_mask(lev+1);
        MultiFab::Multiply(tmp, mask, 0, 0, 1, 0);
    }

    // Get volume including terrain (consistent with volWgtSumMF routine)
    MultiFab volume(grids[lev], dmap[lev], 1, 0);
    auto const& dx = geom[lev].CellSizeArray();
    Real cell_vol  = dx[0]*dx[1]*dx[2];
    volume.setVal(cell_vol);
    if (SolverChoice::mesh_type != MeshType::ConstantDz) {
        MultiFab::Multiply(volume, *detJ_cc[lev], 0, 0, 1, 0);
    }

    //
    // Note that when we send in local = true, NO ParallelAllReduce::Sum
    //      is called inside the Dot product -- we will do that before we print
    //
    bool local = true;
    sum = MultiFab::Dot(tmp, 0, volume, 0, 1, 0, local);

    return sum;
}

/**
 * Helper function for constructing a fine mask, that is, a MultiFab
 * masking coarser data at a lower level by zeroing out covered cells
 * in the fine mask MultiFab we compute.
 *
 * @param level Fine level index which masks underlying coarser data
 */
MultiFab&
ERF::build_fine_mask (int level)
{
    // Mask for zeroing covered cells
    AMREX_ASSERT(level > 0);

    const BoxArray& cba = grids[level-1];
    const DistributionMapping& cdm = dmap[level-1];

    // TODO -- we should make a vector of these a member of ERF class
    fine_mask.define(cba, cdm, 1, 0, MFInfo());
    fine_mask.setVal(1.0);

    BoxArray fba = grids[level];
    iMultiFab ifine_mask = makeFineMask(cba, cdm, fba, ref_ratio[level-1], 1, 0);

    const auto  fma =  fine_mask.arrays();
    const auto ifma = ifine_mask.arrays();
    ParallelFor(fine_mask, [=] AMREX_GPU_DEVICE(int bno, int i, int j, int k) noexcept
    {
       fma[bno](i,j,k) = ifma[bno](i,j,k);
    });

    return fine_mask;
}
