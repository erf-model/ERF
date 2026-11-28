/**
 * \file ERF_NestedBaseState.cpp
 *
 * The base state of a level below the finest the run allows, as the average of the base
 * state of the mesh at max_level resolution, with the cell-center heights placed so that it
 * is in discrete hydrostatic balance.
 *
 * With the coarse metrics built from the finest mesh (metrics_from_finest), the volume of a
 * coarse cell is the sum of the volumes of the fine cells it holds.  Each level's base state
 * is in hydrostatic balance and satisfies the equation of state, but built separately on each
 * level the base-state mass of a coarse cell is not that of the fine cells it holds, so a
 * regrid changes the composite mass by the difference.  Here the coarse base state is instead
 *   1. the volume-weighted average of rho_0, rho_0 theta_0 and rho_0 qv_0 of the base state
 *      built (as the finest level builds its own) on the max_level mesh over the coarse columns,
 *   2. with p_0 and pi_0 from the equation of state, and
 *   3. with the cell-center heights z_phys_cc placed, column by column from the lowest cell
 *      up, so that it is in the discrete hydrostatic balance rebalance_columns enforces:
 *          p_0(k) - p_0(k-1) = -g (z_cc(k) - z_cc(k-1)) (rho_t(k) + rho_t(k-1)) / 2,
 *      rho_t = rho_0 (1 + qv_0).
 * The vertical pressure gradient, in the slow and the fast terms, uses that same distance
 * between the cell centers.  The volumes and areas (from the nodes) are not changed, and
 * z_cc(k) differs from the average of the nodes by a small fraction of the cell.
 *
 * A level's base state depends only on its own columns, not on where the finer levels are,
 * so a regrid changes no base state, and a fine cell's base state is the one its coarse
 * cell's was averaged from.
 */

#include <ERF.H>
#include <ERF_EOS.H>
#include <ERF_Constants.H>
#include <ERF_Utils.H>

#include <AMReX_MultiFabUtil.H>

using namespace amrex;

/**
 * Whether level lev's base state is nested (see the file comment): below max_level, with
 * gravity, on
 *   - a static terrain-fitted mesh with the metrics from the finest mesh (BTF, terrain from the
 *     problem), or
 *   - a flat mesh (ConstantDz or StretchedDz, no terrain), whose finer levels' heights follow
 *     from the levels alone,
 * and a base state built in hydrostatic balance from
 *   - the input sounding (not SoundingType::ConstantDensity), or
 *   - the problem's dry hydrostatic construction (InitType::Isentropic: erf_init_dens_hse_dry
 *     and erf_enforce_hse, as initHSE does).
 *
 * @param[in] lev level of refinement
 */
bool
ERF::nested_base_state (int lev)
{
    if (lev >= max_level) { return false; }
    if (!solverChoice.use_gravity) { return false; }

    // A StretchedDz mesh is flat whatever terrain_type says (erf.initial_dz and
    //    erf.terrain_z_levels set it to StaticFittedMesh); a ConstantDz mesh only without
    //    terrain (not EB or immersed forcing)
    const bool flat = (SolverChoice::mesh_type == MeshType::StretchedDz) ||
                      (SolverChoice::mesh_type == MeshType::ConstantDz &&
                       solverChoice.terrain_type == TerrainType::None);
    if (!flat && !metrics_from_finest(lev)) { return false; }

    if (solverChoice.init_type == InitType::Input_Sounding) {
        return (solverChoice.sounding_type != SoundingType::ConstantDensity);
    }
    return (solverChoice.init_type == InitType::Isentropic);
}

/**
 * Build the nested base state of level lev in bs, a base state on lev's grids, and place the
 * cell-center heights z_phys_cc[lev] to hold it in balance (see the file comment).
 *
 * The work is done on the columns of lev's grids, grown laterally by the ghost cells of the
 * base state and of z_phys_cc and spanning the full height of the domain, so that every
 * column is built from the ground up whatever the grids, and a ghost column is the column it
 * would be were it valid.  The nodes of each finer level over those columns are built from
 * the coarser level's as the levels themselves are (build_zphys_nd), from level 0 up.
 *
 * The lateral ghost cells outside the domain (not across a periodic boundary) of bs take the
 * values of the nearest cell inside it, and those of z_phys_cc keep their offset from the
 * nearest cell inside it, until the boundary conditions fill bs.
 *
 * @param[in]     lev level of refinement, below max_level
 * @param[in,out] bs  base state on lev's grids
 */
void
ERF::build_nested_base_state (int lev, MultiFab& bs)
{
    AMREX_ALWAYS_ASSERT(nested_base_state(lev));
    AMREX_ALWAYS_ASSERT(bs.boxArray() == z_phys_cc[lev]->boxArray() &&
                        bs.DistributionMap() == z_phys_cc[lev]->DistributionMap());

    const Box& domain = geom[lev].Domain();
    const IntVect ngb = bs.nGrowVect();
    const IntVect ngc = z_phys_cc[lev]->nGrowVect();
    const IntVect ngn = z_phys_nd[lev]->nGrowVect();
    const int grow    = std::max({ngb[0], ngb[1], ngc[0], ngc[1]});
    const int ngz     = ngc[2];
    const Real time   = t_new[lev];
    AMREX_ALWAYS_ASSERT(ngn[2] >= ngz);

    //
    // The columns: lev's grids grown laterally, inside the domain, over the full height
    //
    BoxList bl;
    for (int n = 0; n < static_cast<int>(bs.boxArray().size()); ++n) {
        Box b = bs.boxArray()[n];
        b.grow(0, grow); b.grow(1, grow);
        b.setRange(2, domain.smallEnd(2), domain.length(2));
        bl.push_back(b & domain);
    }
    BoxArray ba_col(std::move(bl));
    ba_col.removeOverlap();
    const DistributionMapping dm_col(ba_col);

    //
    // The nodes of each level over the columns, from level 0 up: level 0's own, then each finer
    //    level built from the one below it.  Those of lev are kept for its own cell centers.
    //
    auto ba_at = [&] (int l) {
        BoxArray ba_l(ba_col);
        if (l <= lev) {
            IntVect rr(1);
            for (int m = l; m < lev; ++m) { rr *= refRatio(m); }
            ba_l.coarsen(rr);
        } else {
            IntVect rr(1);
            for (int m = lev; m < l; ++m) { rr *= refRatio(m); }
            ba_l.refine(rr);
        }
        return ba_l;
    };

    std::unique_ptr<MultiFab> znd_lev, znd;
    if (SolverChoice::mesh_type == MeshType::VariableDz) {
        znd = std::make_unique<MultiFab>(convert(ba_at(0), IntVect(1)), dm_col, 1, ngn);
        // The ghost nodes first and then the valid ones, so a node some box holds as valid gets
        //    that value (ParallelCopy does not order overlapping sources)
        znd->ParallelCopy(*z_phys_nd[0], 0, 0, 1, z_phys_nd[0]->nGrowVect(), ngn);
        znd->ParallelCopy(*z_phys_nd[0], 0, 0, 1, IntVect(0), ngn);

        for (int l = 1; l <= max_level; ++l) {
            if (l == lev+1) { znd_lev = std::move(znd); }
            auto znd_l = std::make_unique<MultiFab>(convert(ba_at(l), IntVect(1)), dm_col, 1, ngn);
            build_zphys_nd(l, (l == lev+1) ? znd_lev.get() : znd.get(), *znd_l, time);
            znd_l->FillBoundary(geom[l].periodicity());
            znd = std::move(znd_l);
        }
    } else {
        // A flat mesh: the heights of each level are its levels (or the uniform spacing)
        znd_lev = std::make_unique<MultiFab>(convert(ba_at(lev), IntVect(1)), dm_col, 1, ngn);
        znd     = std::make_unique<MultiFab>(convert(ba_at(max_level), IntVect(1)), dm_col, 1, ngn);
        flat_heights(lev, *znd_lev);
        flat_heights(max_level, *znd);
    }
    AMREX_ALWAYS_ASSERT(znd_lev != nullptr);

    //
    // The base state at max_level over the columns, built as the finest level builds its own
    //
    const BoxArray ba_f = ba_at(max_level);
    auto zcc_f = std::make_unique<MultiFab>(ba_f, dm_col, 1, 1);
    make_zcc(geom[max_level], *znd, *zcc_f);

    // One ghost cell all round, which erf_enforce_hse writes
    MultiFab bs_f(ba_f, dm_col, BaseState::num_comps, 1);
    bs_f.setVal(zero);
    if (solverChoice.init_type == InitType::Input_Sounding) {
        build_base_state_from_input_sounding_on(max_level, bs_f, *znd, zcc_f.get());
    } else {
        // InitType::Isentropic: as initHSE builds a level's base state
        MultiFab r (bs_f, make_alias, BaseState::r0_comp,  1);
        MultiFab pr(bs_f, make_alias, BaseState::p0_comp,  1);
        MultiFab pi(bs_f, make_alias, BaseState::pi0_comp, 1);
        MultiFab th(bs_f, make_alias, BaseState::th0_comp, 1);
        MultiFab qv(bs_f, make_alias, BaseState::qv0_comp, 1);
        prob->erf_init_dens_hse_dry(r, znd, zcc_f, geom[max_level], stretched_dz_h[max_level],
                                    SolverChoice::mesh_type == MeshType::ConstantDz,
                                    SolverChoice::mesh_type == MeshType::StretchedDz);
        erf_enforce_hse(max_level, r, pr, pi, th, qv, zcc_f);
    }

    //
    // The conserved densities, J, J rho_0, J rho_0 theta_0 and J rho_0 qv_0, averaged down one
    //    level at a time to lev (J is then detJ_cc[lev], as make_metrics_from_finest forms it)
    //
    enum { cJ = 0, cR, cRT, cRQ, n_avg };
    auto U = std::make_unique<MultiFab>(ba_f, dm_col, n_avg, 0);
    {
        MultiFab J_f(ba_f, dm_col, 1, 0);
        make_J(geom[max_level], *znd, J_f);  // dz_k/dz on a flat mesh
        for (MFIter mfi(*U, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            const auto& u = U->array(mfi);
            const auto& b = bs_f.const_array(mfi);
            const auto& J = J_f.const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const Real Jr = J(i,j,k) * b(i,j,k,BaseState::r0_comp);
                u(i,j,k,cJ ) = J(i,j,k);
                u(i,j,k,cR ) = Jr;
                u(i,j,k,cRT) = Jr * b(i,j,k,BaseState::th0_comp);
                u(i,j,k,cRQ) = Jr * b(i,j,k,BaseState::qv0_comp);
            });
        }
    }
    znd.reset();

    for (int l = max_level; l > lev; --l) {
        auto U_c = std::make_unique<MultiFab>(ba_at(l-1), dm_col, n_avg, 0);
        average_down(*U, *U_c, 0, n_avg, refRatio(l-1));
        U = std::move(U_c);
    }

    //
    // rho_0, theta_0 and qv_0 of lev, and p_0 and pi_0 from the equation of state
    //
    MultiFab bs_col(ba_col, dm_col, BaseState::num_comps, 0);
    for (MFIter mfi(bs_col, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const auto& u = U->const_array(mfi);
        const auto& b = bs_col.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            b(i,j,k,BaseState::r0_comp ) = u(i,j,k,cR ) / u(i,j,k,cJ);
            b(i,j,k,BaseState::th0_comp) = u(i,j,k,cRT) / u(i,j,k,cR);
            b(i,j,k,BaseState::qv0_comp) = u(i,j,k,cRQ) / u(i,j,k,cR);
        });
    }
    U.reset();
    recompute_p_pi_hse(bs_col, solverChoice.rdOcp);

    //
    // The cell-center heights: the lowest cell keeps the average of its nodes, and each cell
    //    above it is placed so that the column is in balance; the ghost cells below and above
    //    keep their offsets from the lowest and highest cells.  Each height must lie between
    //    the lowest and the highest node of its cell.
    //
    MultiFab zcc_col(ba_col, dm_col, 1, IntVect(0,0,ngz));
    make_zcc(geom[lev], *znd_lev, zcc_col);

    const int kbot = domain.smallEnd(2);
    const int ktop = domain.bigEnd(2);
    const Real grav = CONST_GRAV;

    MultiFab bad(ba_col, dm_col, 1, 0);
    bad.setVal(zero);
    for (MFIter mfi(zcc_col); mfi.isValid(); ++mfi) {
        Box slab = mfi.validbox(); slab.makeSlab(2, kbot);
        const auto& z   = zcc_col.array(mfi);
        const auto& b   = bs_col.const_array(mfi);
        const auto& zn  = znd_lev->const_array(mfi);
        const auto& nb  = bad.array(mfi);
        ParallelFor(slab, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            const Real z_geo_bot = z(i,j,kbot);
            const Real z_geo_top = z(i,j,ktop);
            for (int k = kbot+1; k <= ktop; ++k) {
                const Real rt_lo = b(i,j,k-1,BaseState::r0_comp) * (one + b(i,j,k-1,BaseState::qv0_comp));
                const Real rt_hi = b(i,j,k  ,BaseState::r0_comp) * (one + b(i,j,k  ,BaseState::qv0_comp));
                z(i,j,k) = z(i,j,k-1) - (b(i,j,k,BaseState::p0_comp) - b(i,j,k-1,BaseState::p0_comp))
                                        / (grav * myhalf * (rt_lo + rt_hi));
            }
            for (int kk = 1; kk <= ngz; ++kk) {
                z(i,j,kbot-kk) += z(i,j,kbot) - z_geo_bot;
                z(i,j,ktop+kk) += z(i,j,ktop) - z_geo_top;
            }
            for (int k = kbot; k <= ktop; ++k) {
                const Real zlo = amrex::min(amrex::min(zn(i,j,k  ), zn(i+1,j,k  )), amrex::min(zn(i,j+1,k  ), zn(i+1,j+1,k  )));
                const Real zhi = amrex::max(amrex::max(zn(i,j,k+1), zn(i+1,j,k+1)), amrex::max(zn(i,j+1,k+1), zn(i+1,j+1,k+1)));
                if (!(z(i,j,k) > zlo && z(i,j,k) < zhi)) { nb(i,j,k) = one; }
            }
        });
    }
    znd_lev.reset();
    if (bad.max(0) > zero) {
        const IntVect loc = bad.maxIndex(0);
        Abort("build_nested_base_state: at level " + std::to_string(lev) +
              ", the height that holds the nested base state in balance lies outside its cell " +
              "at " + std::to_string(loc[0]) + " " + std::to_string(loc[1]) + " " + std::to_string(loc[2]));
    }

    //
    // Into bs and z_phys_cc[lev]: the valid cells and the lateral ghost cells inside the domain
    //    (and across a periodic boundary); below and above the domain, z_phys_cc's ghost cells
    //    come from the columns and bs's are left to the boundary conditions
    //
    // The average of the nodes, for the offsets of the ghost cells outside the domain
    MultiFab zcc_old(z_phys_cc[lev]->boxArray(), z_phys_cc[lev]->DistributionMap(), 1, ngc);
    make_zcc(geom[lev], *z_phys_nd[lev], zcc_old);

    bs.ParallelCopy(bs_col, 0, 0, BaseState::num_comps, IntVect(0), IntVect(ngb[0],ngb[1],0),
                    geom[lev].periodicity());
    z_phys_cc[lev]->ParallelCopy(zcc_col, 0, 0, 1, IntVect(0,0,ngz), ngc, geom[lev].periodicity());

    Box dom_g(domain);
    for (int dir = 0; dir < 2; ++dir) {
        if (geom[lev].isPeriodic(dir)) { dom_g.grow(dir, std::max(ngb[dir], ngc[dir])); }
    }
    const auto dlo = lbound(dom_g);
    const auto dhi = ubound(dom_g);
    const int nc = BaseState::num_comps;
    for (MFIter mfi(bs); mfi.isValid(); ++mfi) {
        const auto& d  = bs.array(mfi);
        const auto& z  = z_phys_cc[lev]->array(mfi);
        const auto& zo = zcc_old.const_array(mfi);
        Box gb = mfi.fabbox(); gb.setRange(2, domain.smallEnd(2), domain.length(2));
        ParallelFor(gb, nc, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int ii = amrex::min(amrex::max(i, dlo.x), dhi.x);
            const int jj = amrex::min(amrex::max(j, dlo.y), dhi.y);
            if (ii != i || jj != j) { d(i,j,k,n) = d(ii,jj,k,n); }
        });
        Box gz = mfi.validbox(); gz.grow(ngc);
        ParallelFor(gz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int ii = amrex::min(amrex::max(i, dlo.x), dhi.x);
            const int jj = amrex::min(amrex::max(j, dlo.y), dhi.y);
            if (ii != i || jj != j) { z(i,j,k) = z(ii,jj,k) + (zo(i,j,k) - zo(ii,jj,k)); }
        });
    }
    Gpu::streamSynchronize();

    if (physbcs_base[lev]) {
        (*physbcs_base[lev])(bs, 0, bs.nComp(), bs.nGrowVect());
    }

    // The vertical pressure gradient of this level uses these heights (already so on a mesh other
    //    than ConstantDz)
    if (static_cast<int>(solverChoice.zcc_vertical_metric_lev.size()) <= lev) {
        solverChoice.zcc_vertical_metric_lev.resize(lev+1, 0);
    }
    solverChoice.zcc_vertical_metric_lev[lev] = 1;
}

/**
 * The nodal heights of a flat mesh (ConstantDz or StretchedDz) at level lev's resolution, on
 * the boxes of z_nd and its ghost nodes: the uniform spacing from geometry.prob_lo[2], or the
 * stretched levels (extended past either end by the end spacing).
 *
 * @param[in]  lev  level whose resolution z_nd has
 * @param[out] z_nd nodal heights
 */
void
ERF::flat_heights (int lev, MultiFab& z_nd)
{
    const Real z_lo = geom[lev].ProbLo(2);
    const Real dz   = geom[lev].CellSize(2);
    const bool stretched = (SolverChoice::mesh_type == MeshType::StretchedDz);
    const int  nzs = stretched ? static_cast<int>(zlevels_stag[lev].size()) : 0;
    Gpu::DeviceVector<Real> d_zs(nzs);
    if (stretched) {
        AMREX_ALWAYS_ASSERT(nzs >= 2);
        Gpu::copy(Gpu::hostToDevice, zlevels_stag[lev].begin(), zlevels_stag[lev].end(), d_zs.begin());
    }
    const Real* zs = d_zs.data();
    for (MFIter mfi(z_nd); mfi.isValid(); ++mfi) {
        const auto& z = z_nd.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (stretched) {
                if (k < 0) {
                    z(i,j,k) = zs[0] + static_cast<Real>(k) * (zs[1] - zs[0]);
                } else if (k > nzs-1) {
                    z(i,j,k) = zs[nzs-1] + static_cast<Real>(k-nzs+1) * (zs[nzs-1] - zs[nzs-2]);
                } else {
                    z(i,j,k) = zs[k];
                }
            } else {
                z(i,j,k) = z_lo + static_cast<Real>(k) * dz;
            }
        });
    }
    Gpu::streamSynchronize();
}
