/**
 * \file ERF_CFFluxHistory.cpp
 *
 * The fluxes through the coarse/fine faces of a refined level, taken from the coarse level's
 * acoustic substeps.
 *
 * The density flux that each level applies over a time step is the time integral, over its
 * acoustic substeps, of the momentum on each face (the slow part, built from the stage
 * momenta, cancels against the stage momenta in the fast update).  The fine level's
 * coarse/fine faces are set from the coarse level at the substep times, interpolated linearly
 * in time between the coarse level's old and new states.  Unless the coarse momentum on a face
 * changes linearly over the coarse step -- and acoustic waves through a face make it anything
 * but -- the two integrals differ.  Density is not refluxed, so the difference is lost: mass
 * is not conserved, and at large time steps the run blows up (issue #4325).
 *
 * So the coarse level records, in thin slabs of faces under the coarse/fine faces of the fine
 * level, the flux each face applied over each substep of its last RK stage (the stage whose
 * fluxes are the step's), and the vertical quantity on its z-faces at the substep times:
 * (rho w) without terrain, Omega with it.  The fine level, in its own last
 * RK stage, takes on each coarse/fine face the average of the coarse flux over each of its
 * substeps, interpolated conservatively in space, and in z that vertical quantity at the end of
 * the substep.
 * Whatever the substeps of the two levels, its fluxes through the face then add up, over the
 * coarse step, to the coarse level's.
 *
 * With terrain (erf_substep_T) the fluxes are those the terrain integrator forms: h_zeta (rho u)/m
 * through x- and y-faces, whose area h_zeta is, and Omega through z-faces.  The fine level sets
 * the Omega of its coarse/fine z-faces, and so (rho w) there follows from its own (rho u) and
 * (rho v) and the terrain slopes.
 *
 * This is done for the fast integrators without terrain and with fixed terrain (erf_substep_NS,
 * erf_substep_T), not with moving terrain (erf_substep_MT).  It is done under one-way coupling
 * too, so that the fine level sees the same fluxes through its coarse/fine faces whether or
 * not it is averaged down and refluxed onto the coarse level.
 */
#include "ERF.H"
#include "ERF_Constants.H"

using namespace amrex;

/**
 * Whether the coarse/fine face fluxes of level crse_lev+1 come from the substeps of level
 * crse_lev
 *
 * @param crse_lev the coarse level
 */
bool
ERF::cf_flux_history_active (int crse_lev) const
{
    // erf.cf_flux_history = false keeps the coarse/fine faces interpolated linearly in time
    static const bool enabled = [] {
        bool on = true;
        ParmParse("erf").query("cf_flux_history", on);
        return on;
    }();
    if (!enabled) { return false; }
    if (crse_lev >= finest_level) { return false; }
    if (solverChoice.terrain_type == TerrainType::MovingFittedMesh ||
        solverChoice.terrain_type == TerrainType::EB) { return false; }
    for (int lev = crse_lev; lev <= crse_lev+1; ++lev) {
        if (solverChoice.anelastic[lev] == 1) { return false; }
        if (solverChoice.substepping_type[lev] == SubsteppingType::None) { return false; }
    }
    return true;
}

/**
 * The history of level lev's coarse/fine face fluxes, its layout made again when the level's
 * grids or distribution have changed
 *
 * @param lev the fine level (> 0)
 */
ERF::CFFluxHistory&
ERF::cf_flux_history (int lev)
{
    AMREX_ASSERT(lev > 0);
    if (static_cast<int>(m_cf_flux_hist.size()) <= lev) {
        m_cf_flux_hist.resize(lev+1);
    }
    CFFluxHistory& H = m_cf_flux_hist[lev];

    const BoxArray& fba = grids[lev];
    const DistributionMapping& fdm = dmap[lev];
    if (H.fba == fba && H.fdm == fdm) { return H; }

    H = CFFluxHistory{};
    H.fba = fba;
    H.fdm = fdm;

    const Box& fdom = geom[lev].Domain();
    const IntVect rr = refRatio(lev-1);

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        BoxList slabs;
        Vector<int> owner;
        for (int n = 0; n < fba.size(); ++n) {
            const Box& b = fba[n];
            for (int side = 0; side < 2; ++side) {
                // The layer of cells just outside this side of the box ...
                const int kc = (side == 0) ? b.smallEnd(dir) - 1 : b.bigEnd(dir) + 1;
                const int kf = (side == 0) ? b.smallEnd(dir)     : b.bigEnd(dir) + 1;
                int shift = 0;
                if (kc < fdom.smallEnd(dir) || kc > fdom.bigEnd(dir)) {
                    if (!geom[lev].isPeriodic(dir)) { continue; }
                    shift = (kc < fdom.smallEnd(dir)) ? fdom.length(dir) : -fdom.length(dir);
                }
                Box outside(b);
                outside.setRange(dir, kc + shift);

                // ... where it is not in the level: the box's coarse/fine faces
                for (Box f : fba.complementIn(outside)) {
                    f.shift(dir, -shift);
                    f.surroundingNodes(dir);
                    f.setRange(dir, kf);

                    // The coarse faces under them, with a face on each side for the slopes
                    Box slab = amrex::coarsen(f, rr);
                    for (int t = 0; t < AMREX_SPACEDIM; ++t) {
                        if (t != dir) { slab.grow(t, 1); }
                    }
                    slabs.push_back(slab);
                    owner.push_back(fdm[n]);
                    H.slab_box[dir].push_back(n);
                    H.slab_faces[dir].push_back(f);
                }
            }
        }
        if (!slabs.isEmpty()) {
            H.sba[dir] = BoxArray(slabs);
            H.sdm[dir] = DistributionMapping(owner);
        }
    }

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        H.tgt[dir].define(convert(fba, IntVect::TheDimensionVector(dir)), fdm, 2, 0);
    }
    H.tgt_w.define(convert(fba, IntVect::TheDimensionVector(2)), fdm, 1, 0);

    return H;
}

/**
 * Record the fluxes through the faces of level crse_lev over one of the acoustic substeps of
 * its last RK stage, under the coarse/fine faces of level crse_lev+1
 *
 * @param crse_lev  the coarse level
 * @param fast_step the substep, from 0
 * @param n_sub     the number of substeps in the stage
 * @param t_beg     the time at the start of the substep
 * @param t_end     the time at its end
 * @param rec       the flux through each face over the substep (CFFaceFluxes::rec)
 * @param rho_w_old (rho w) at the start of the substep, Omega with terrain
 * @param rho_w_new (rho w) at its end, Omega with terrain
 */
void
ERF::cf_flux_history_record (int crse_lev, int fast_step, int n_sub, Real t_beg, Real t_end,
                             const Array<MultiFab, AMREX_SPACEDIM>& rec,
                             const MultiFab& rho_w_old, const MultiFab& rho_w_new)
{
    CFFluxHistory& H = cf_flux_history(crse_lev+1);
    const Periodicity& period = geom[crse_lev].periodicity();

    if (fast_step == 0) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            H.flux[dir].clear();
            H.flux[dir].resize(n_sub);
        }
        H.rho_w.clear();
        H.rho_w.resize(n_sub+1);
        H.time.assign(n_sub+1, Real(0));
        H.time[0] = t_beg;
        if (!H.sba[2].empty()) {
            H.rho_w[0].define(H.sba[2], H.sdm[2], 1, 0);
            H.rho_w[0].setVal(zero);
            H.rho_w[0].ParallelCopy(rho_w_old, 0, 0, 1, IntVect(0), IntVect(0), period);
        }
    }
    H.complete = false;
    if (static_cast<int>(H.flux[0].size()) != n_sub) { return; }

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (H.sba[dir].empty()) { continue; }
        MultiFab& f = H.flux[dir][fast_step];
        f.define(H.sba[dir], H.sdm[dir], 1, 0);
        f.setVal(zero);
        f.ParallelCopy(rec[dir], 0, 0, 1, IntVect(0), IntVect(0), period);
    }
    if (!H.sba[2].empty()) {
        MultiFab& w = H.rho_w[fast_step+1];
        w.define(H.sba[2], H.sdm[2], 1, 0);
        w.setVal(zero);
        w.ParallelCopy(rho_w_new, 0, 0, 1, IntVect(0), IntVect(0), period);
    }
    H.time[fast_step+1] = t_end;
    H.complete = (fast_step == n_sub-1);
}

/**
 * Fill CFFluxHistory::tgt and tgt_w of level lev for one of its substeps: on its coarse/fine
 * faces, the coarse level's flux averaged over [t_beg, t_end], and (rho w) at t_end (Omega
 * with terrain)
 *
 * @param lev   the fine level (> 0)
 * @param t_beg the time at the start of the substep
 * @param t_end the time at its end
 */
void
ERF::cf_flux_history_targets (int lev, Real t_beg, Real t_end)
{
    CFFluxHistory& H = cf_flux_history(lev);
    if (!H.complete) {
        Abort("cf_flux_history_targets: the coarse level has not recorded its substeps");
    }

    // The weight of each coarse substep in this fine one: the fraction of the fine substep
    // it overlaps.  Over the fine substeps of a coarse step they add up to the coarse substeps.
    const int n_sub = static_cast<int>(H.time.size()) - 1;
    const Real dt_sub = t_end - t_beg;
    Vector<Real> wgt(n_sub, Real(0));
    Real wsum = Real(0);
    for (int m = 0; m < n_sub; ++m) {
        Real overlap = amrex::min(H.time[m+1], t_end) - amrex::max(H.time[m], t_beg);
        wgt[m] = amrex::max(overlap, Real(0)) / dt_sub;
        wsum += wgt[m];
    }
    if (std::abs(wsum - one) > Real(1.e-6)) {
        Abort("cf_flux_history_targets: the substep is not inside the coarse level's last stage");
    }

    // (rho w), or Omega with terrain, at t_end, linear between the coarse substep times
    int mw = 0;
    while (mw < n_sub-1 && H.time[mw+1] < t_end) { ++mw; }
    Real fw = (t_end - H.time[mw]) / (H.time[mw+1] - H.time[mw]);
    fw = amrex::min(amrex::max(fw, Real(0)), one);

    const IntVect rr = refRatio(lev-1);
    const Box& cdom = geom[lev-1].Domain();
    const GpuArray<int, AMREX_SPACEDIM> per {AMREX_D_DECL(geom[lev-1].isPeriodic(0),
                                                          geom[lev-1].isPeriodic(1),
                                                          geom[lev-1].isPeriodic(2))};

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        H.tgt[dir].setVal(zero);
    }
    H.tgt_w.setVal(zero);

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (H.sba[dir].empty()) { continue; }
        const bool do_w = (dir == 2);
        for (MFIter mfi(H.flux[dir][0]); mfi.isValid(); ++mfi) {
            const Box& sbx = mfi.validbox();
            FArrayBox avg_fab(sbx, do_w ? 2 : 1, The_Async_Arena());
            const Array4<Real> avg = avg_fab.array();
            ParallelFor(sbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                avg(i,j,k,0) = zero;
            });
            for (int m = 0; m < n_sub; ++m) {
                if (wgt[m] == Real(0)) { continue; }
                const Real w = wgt[m];
                const Array4<const Real> f = H.flux[dir][m].const_array(mfi);
                ParallelFor(sbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    avg(i,j,k,0) += w * f(i,j,k);
                });
            }
            if (do_w) {
                const Array4<const Real> w0 = H.rho_w[mw  ].const_array(mfi);
                const Array4<const Real> w1 = H.rho_w[mw+1].const_array(mfi);
                ParallelFor(sbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    avg(i,j,k,1) = (one - fw) * w0(i,j,k) + fw * w1(i,j,k);
                });
            }

            const int n = H.slab_box[dir][mfi.index()];
            const Array4<Real> tgt   = H.tgt[dir].array(n);
            const Array4<Real> tgt_w = H.tgt_w.array(n);
            const int ncomp = do_w ? 2 : 1;

            // Linear in the directions along the face, with monotonized central slopes: the
            // fine faces under a coarse face are placed symmetrically about it, so their
            // values average to the coarse one
            ParallelFor(H.slab_faces[dir][mfi.index()], [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const int fi[3] = {i, j, k};
                int ci[3];
                Real off[3];
                for (int t = 0; t < 3; ++t) {
                    ci[t]  = amrex::coarsen(fi[t], rr[t]);
                    off[t] = (t == dir) ? Real(0)
                           : (Real(fi[t] - ci[t]*rr[t]) + myhalf) / Real(rr[t]) - myhalf;
                }
                for (int c = 0; c < ncomp; ++c) {
                    const Real a0 = avg(ci[0],ci[1],ci[2],c);
                    Real val = a0;
                    for (int t = 0; t < 3; ++t) {
                        if (t == dir || rr[t] == 1) { continue; }
                        const bool lo_ok = per[t] || ci[t]-1 >= cdom.smallEnd(t);
                        const bool hi_ok = per[t] || ci[t]+1 <= cdom.bigEnd(t);
                        if (!lo_ok || !hi_ok) { continue; }
                        int cm[3] = {ci[0], ci[1], ci[2]}; cm[t] -= 1;
                        int cp[3] = {ci[0], ci[1], ci[2]}; cp[t] += 1;
                        const Real dl = a0 - avg(cm[0],cm[1],cm[2],c);
                        const Real dr = avg(cp[0],cp[1],cp[2],c) - a0;
                        Real slope = Real(0);
                        if (dl*dr > Real(0)) {
                            const Real dc = myhalf * (dl + dr);
                            slope = amrex::min(amrex::Math::abs(dc), Real(2)*amrex::Math::abs(dl), Real(2)*amrex::Math::abs(dr));
                            slope = (dc > Real(0)) ? slope : -slope;
                        }
                        val += slope * off[t];
                    }
                    if (c == 0) {
                        tgt(i,j,k,0) = val;
                        tgt(i,j,k,1) = one;
                    } else {
                        tgt_w(i,j,k) = val;
                    }
                }
            });
        }
    }
}
