#include <ERF.H>
#include <Utils.H>

using namespace amrex;

#if defined(ERF_USE_WARM_NO_PRECIP)
void
ERF::condensation_source (MultiFab& source, MultiFab& S_new, Real dt_lev)
{
    for ( MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();
        auto  s_arr = S_new.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // qv_array(i,j,k)  = qt_array(i,j,k) - qn_array(i,j,k);
        });

    } // mfi
}
#endif
