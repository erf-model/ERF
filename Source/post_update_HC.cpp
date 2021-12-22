#include "ET_Integration.H"

using namespace amrex;

void rescale_state (MultiFab& state_mf)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
  for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();
    const auto ncomp = state_mf.nComp();

    const auto& state_fab = state_mf.array(mfi);

    // For each grid, loop over all the valid points
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
	// for example:
        

	});
  }
}
