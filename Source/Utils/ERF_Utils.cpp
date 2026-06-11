#include "ERF_Utils.H"

using namespace amrex;

void
cons_to_prim(const MultiFab& cons_state, MultiFab& S_prim, int ng)
{
    BL_PROFILE("cons_to_prim()");

    int ncomp_prim = S_prim.nComp();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(cons_state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& gbx = mfi.growntilebox(ng);
        const Array4<const Real>& cons_arr     = cons_state.array(mfi);
        const Array4<      Real>& prim_arr     = S_prim.array(mfi);

        //
        // We may need > one ghost cells of prim in order to compute higher order advective terms
       //
       amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
           Real rho       = cons_arr(i,j,k,Rho_comp);
           Real rho_theta = cons_arr(i,j,k,RhoTheta_comp);
           prim_arr(i,j,k,PrimTheta_comp) = rho_theta / rho;
           for (int n = 1; n < ncomp_prim; ++n) {
               prim_arr(i,j,k,PrimTheta_comp + n) = cons_arr(i,j,k,RhoTheta_comp + n) / rho;
           }
       });
    } // mfi
};

void
make_qt(const MultiFab& cons_state, MultiFab& qt, int n_qstate_into_total)
{
    BL_PROFILE("make_qt()");

    // All moisture models are guaranteed to have RhoQ1_comp.
    MultiFab::Copy(qt, cons_state, RhoQ1_comp, 0, 1, qt.nGrowVect());

    for (int n = 1; n < n_qstate_into_total; ++n) {
        MultiFab::Add(qt, cons_state, RhoQ1_comp+n, 0, 1, qt.nGrowVect());
    }

    MultiFab::Divide(qt, cons_state, Rho_comp, 0, 1, qt.nGrowVect());
}
