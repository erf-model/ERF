#include "ERF_Utils.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_Reduce.H>

using namespace amrex;

/**
 * Convert conservative variables to primitive variables.
 *
 * @param[in] cons_state MultiFab containing conservative state variables
 * @param[out] S_prim MultiFab to be filled with primitive state variables
 * @param[in] ng Number of ghost cells
 */
void
cons_to_prim (const MultiFab& cons_state, MultiFab& S_prim, int ng)
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
make_qt (const MultiFab& cons_state, MultiFab& qt, int n_qstate_into_total)
{
    BL_PROFILE("make_qt()");

    // All moisture models are guaranteed to have RhoQ1_comp.
    MultiFab::Copy(qt, cons_state, RhoQ1_comp, 0, 1, qt.nGrowVect());

    for (int n = 1; n < n_qstate_into_total; ++n) {
        MultiFab::Add(qt, cons_state, RhoQ1_comp+n, 0, 1, qt.nGrowVect());
    }

    MultiFab::Divide(qt, cons_state, Rho_comp, 0, 1, qt.nGrowVect());
}

/**
 * Fluid-mass-weighted average down over immersed-forcing levels (see ERF_Utils.H).
 */
void
if_average_down (const MultiFab& S_fine, MultiFab& S_crse,
                 const MultiFab& beta_fine, const MultiFab& beta_crse,
                 const MultiFab& vol_fine,
                 int scomp, int ncomp, const IntVect& ratio)
{
    BL_PROFILE("if_average_down()");
    AMREX_ALWAYS_ASSERT(S_fine.is_cell_centered() && S_crse.is_cell_centered());
    AMREX_ALWAYS_ASSERT(scomp >= 0 && ncomp > 0 && scomp + ncomp <= S_crse.nComp());
    AMREX_ALWAYS_ASSERT(S_fine.boxArray().coarsenable(ratio));
    AMREX_ALWAYS_ASSERT(beta_fine.boxArray() == S_fine.boxArray() && vol_fine.boxArray() == S_fine.boxArray());
    AMREX_ALWAYS_ASSERT(beta_fine.DistributionMap() == S_fine.DistributionMap() &&
                        vol_fine.DistributionMap() == S_fine.DistributionMap());
    AMREX_ALWAYS_ASSERT(beta_crse.boxArray() == S_crse.boxArray() &&
                        beta_crse.DistributionMap() == S_crse.DistributionMap());

    // Per coarse cell, on the coarsened fine layout: the fluid-mass-weighted sums of the
    // components, the fluid volume and the fluid mass
    BoxArray cba = S_fine.boxArray(); cba.coarsen(ratio);
    const int nsum = ncomp + 2;
    MultiFab csum(cba, S_fine.DistributionMap(), nsum, 0);
    const int rx = ratio[0], ry = ratio[1], rz = ratio[2];

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(csum, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const c = csum.array(mfi);
        auto const f = S_fine.const_array(mfi);
        auto const b = beta_fine.const_array(mfi);
        auto const v = vol_fine.const_array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real w = zero, wr = zero;
            for (int n = 0; n < ncomp; ++n) { c(i,j,k,n) = zero; }
            for (int kk = k*rz; kk < (k+1)*rz; ++kk) {
            for (int jj = j*ry; jj < (j+1)*ry; ++jj) {
            for (int ii = i*rx; ii < (i+1)*rx; ++ii) {
                // fluid volume of the fine cell (prvalues: the constants are not odr-used)
                const Real wf = (one - amrex::min(amrex::max(b(ii,jj,kk), Real(zero)), Real(one)))
                              * v(ii,jj,kk);
                w  += wf;
                wr += wf * f(ii,jj,kk,Rho_comp);
                for (int n = 0; n < ncomp; ++n) { c(i,j,k,n) += wf * f(ii,jj,kk,scomp+n); }
            }}}
            c(i,j,k,ncomp)   = w;
            c(i,j,k,ncomp+1) = wr;
        });
    }

    // The sums on the coarse layout: zero where the fine level does not cover the coarse level
    MultiFab csum_c(S_crse.boxArray(), S_crse.DistributionMap(), nsum, 0);
    csum_c.setVal(zero);
    csum_c.ParallelCopy(csum, 0, 0, nsum);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S_crse, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const s  = S_crse.array(mfi);
        auto const c  = csum_c.const_array(mfi);
        auto const bc = beta_crse.const_array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // uncovered coarse cells, coarse cells solid on their own level and coarse cells
            // without a fluid fine cell under them keep their value
            if (c(i,j,k,ncomp) > zero && c(i,j,k,ncomp+1) > zero && bc(i,j,k) < one) {
                // the coarse density: the fluid-volume mean of the fine density when it is in
                // the range, else the coarse cell's own
                const Real rho_c = (scomp == Rho_comp) ? c(i,j,k,ncomp+1) / c(i,j,k,ncomp)
                                                       : s(i,j,k,Rho_comp);
                const Real fac = rho_c / c(i,j,k,ncomp+1);
                for (int n = 0; n < ncomp; ++n) {
                    s(i,j,k,scomp+n) = (scomp+n == Rho_comp) ? rho_c : fac * c(i,j,k,n);
                }
            }
        });
    }
}

std::pair<int, Real>
if_cf_surface_mismatch (const MultiFab& fine_tb,
                        const BoxArray& fine_grids,
                        const BoxArray& crse_grids,
                        const DistributionMapping& crse_dmap,
                        const Geometry& crse_geom,
                        const IntVect& ratio,
                        Real warn)
{
    const IntVect rr = ratio;
    const iMultiFab covered = makeFineMask(crse_grids, crse_dmap, IntVect(1,1,0), fine_grids, rr,
                                           crse_geom.periodicity(), 0, 1);

    // spread of the fine fractions under each coarse cell, formed on the coarsened fine layout
    MultiFab ftmp(amrex::coarsen(fine_grids, rr), fine_tb.DistributionMap(), 1, 0);
    const int rx = rr[0], ry = rr[1], rz = rr[2];
    for (MFIter mfi(ftmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const s = ftmp.array(mfi);
        auto const f = fine_tb.const_array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real bmin = Real(one), bmax = Real(zero);
            for (int kk = k*rz; kk < (k+1)*rz; ++kk) {
            for (int jj = j*ry; jj < (j+1)*ry; ++jj) {
            for (int ii = i*rx; ii < (i+1)*rx; ++ii) {
                const Real b = amrex::min(amrex::max(f(ii,jj,kk), Real(zero)), Real(one));
                bmin = amrex::min(bmin, b);
                bmax = amrex::max(bmax, b);
            }}}
            s(i,j,k) = bmax - bmin;
        });
    }
    MultiFab spread(crse_grids, crse_dmap, 1, 0);
    spread.setVal(zero);
    spread.ParallelCopy(ftmp, 0, 0, 1);

    const Box& dom = crse_geom.Domain();
    const int per_x = crse_geom.isPeriodic(0), per_y = crse_geom.isPeriodic(1);
    const int dlo_x = dom.smallEnd(0), dhi_x = dom.bigEnd(0), dlo_y = dom.smallEnd(1), dhi_y = dom.bigEnd(1);
    ReduceOps<ReduceOpSum, ReduceOpMax> reduce_op;
    ReduceData<int, Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (MFIter mfi(spread, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const m  = covered.const_array(mfi);
        auto const sp = spread.const_array(mfi);
        reduce_op.eval(mfi.tilebox(), reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
        {
            if (m(i,j,k) != 1) { return {0, Real(zero)}; }
            const bool edge = ((per_x || i > dlo_x) && m(i-1,j,k) == 0) || ((per_x || i < dhi_x) && m(i+1,j,k) == 0)
                           || ((per_y || j > dlo_y) && m(i,j-1,k) == 0) || ((per_y || j < dhi_y) && m(i,j+1,k) == 0);
            if (!edge) { return {0, Real(zero)}; }
            return {(sp(i,j,k) >= warn) ? 1 : 0, sp(i,j,k)};
        });
    }
    ReduceTuple hv = reduce_data.value(reduce_op);
    int  nflagged = amrex::get<0>(hv);
    Real worst    = amrex::get<1>(hv);
    ParallelDescriptor::ReduceIntSum(nflagged);
    ParallelDescriptor::ReduceRealMax(worst);
    return {nflagged, worst};
}
