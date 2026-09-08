#include "ERF_Utils.H"

using namespace amrex;

/**
 * Convert conservative variables to primitive variables.
 *
 * @param[in] cons_state MultiFab containing conservative state variables
 * @param[out] S_prim MultiFab to be filled with primitive state variables
 * @param[in] ng Number of ghost cells
 */
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
}

/**
 * Fill the ghost cells of the wall distance field.
 *
 * The wall distance is only computed on the valid region of each box (whether
 * from the analytic shortcuts, the Poisson solve, or the thin-body
 * correction), so without this the ghost cells still hold the large sentinel
 * value they were initialized with. Any operator that reads wall distance from
 * a grown box -- e.g. the RANS Dirichlet-k boundary condition in
 * SurfaceLayer::update_fluxes -- would then propagate that sentinel into the
 * state.
 *
 * Ghost cells shared with a neighboring box or across a periodic boundary are
 * filled by FillBoundary; ghost cells outside the domain in a non-periodic
 * direction are filled by zero-order extrapolation. Note that the analytic
 * wall distances vary only with k, so the extrapolation is exact in those
 * cases.
 *
 * @param[inout] wdist MultiFab holding the wall distance
 * @param[in] geom Geometry at this level
 */
void
fill_wall_dist_ghost_cells (MultiFab& wdist, const Geometry& geom)
{
    BL_PROFILE("fill_wall_dist_ghost_cells()");

    wdist.FillBoundary(geom.periodicity());

    const Box& domain = geom.Domain();
    const auto dom_lo = amrex::lbound(domain);
    const auto dom_hi = amrex::ubound(domain);

    const GpuArray<int,AMREX_SPACEDIM> is_per = {AMREX_D_DECL(geom.isPeriodic(0),
                                                              geom.isPeriodic(1),
                                                              geom.isPeriodic(2))};

    for (MFIter mfi(wdist); mfi.isValid(); ++mfi)
    {
        const Box& gbx = mfi.fabbox();
        const Array4<Real>& d_arr = wdist.array(mfi);

        //
        // Note that the cells we read from are always inside the domain in
        // every non-periodic direction, hence already filled, and are disjoint
        // from the cells we write to.
        //
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int ii = (is_per[0]) ? i : amrex::min(amrex::max(i,dom_lo.x),dom_hi.x);
            const int jj = (is_per[1]) ? j : amrex::min(amrex::max(j,dom_lo.y),dom_hi.y);
            const int kk = (is_per[2]) ? k : amrex::min(amrex::max(k,dom_lo.z),dom_hi.z);
            if (ii != i || jj != j || kk != k) {
                d_arr(i,j,k) = d_arr(ii,jj,kk);
            }
        });
    } // mfi
}

/**
 * Compute the total water mixing ratio from conservative variables.
 *
 * @param[in] cons_state MultiFab containing conservative state variables
 * @param[out] qt MultiFab to store the total water mixing ratio
 * @param[in] n_qstate_into_total Number of moisture components to include in the total
 */
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
