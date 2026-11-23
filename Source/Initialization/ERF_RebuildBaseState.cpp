/**
 * \file ERF_RebuildBaseState.cpp
 *
 * The base state of a level that a regrid made or remade, and the check that a level's base
 * state is in discrete hydrostatic balance (erf-model/ERF#5228, #5236).
 */

#include <ERF.H>
#include <iomanip>
#include <sstream>
#include <ERF_EOS.H>
#include <ERF_Constants.H>
#include <ERF_Utils.H>
#include <ERF_HSEUtils.H>
#include <ERF_ColumnBands.H>

using namespace amrex;

/**
 * Re-compute p_0 and pi_0 from rho_0, theta_0 and qv_0 over the valid region.
 *
 * @param[in,out] base_state_mf base state, with all BaseState components
 * @param[in]     l_rdOcp       R_d / c_p
 */
void
recompute_p_pi_hse (MultiFab& base_state_mf, const Real l_rdOcp)
{
    for (MFIter mfi(base_state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox();

        const Array4<Real>& bs = base_state_mf.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhotheta = bs(i,j,k,BaseState::r0_comp) * bs(i,j,k,BaseState::th0_comp);
            bs(i,j,k,BaseState::p0_comp)  = getPgivenRTh(rhotheta, bs(i,j,k,BaseState::qv0_comp));
            bs(i,j,k,BaseState::pi0_comp) = getExnergivenRTh(rhotheta, l_rdOcp, bs(i,j,k,BaseState::qv0_comp));
        });
    }
}

/**
 * A copy of bs, a base state of level lev, on the boxes of bs's BoxArray grown laterally by
 * grow, and level lev's heights on the same boxes.
 *
 * The grown boxes overlap where bs's boxes are neighbors; a column they share is built from
 * the same heights in each.  A base state built in their valid cells and copied back from
 * those alone, base_state.ParallelCopy(bs_g, 0, 0, ncomp, IntVect(0,0,ng[2]), bs.nGrowVect()),
 * therefore holds in each of bs's lateral ghost cells (grow = bs's lateral ghost cells) the
 * column it would hold were the cell valid, the same whichever box it is a ghost cell of.
 *
 * Where the heights of a whole column are known without z_phys_nd[lev] (MeshType::ConstantDz
 * and StretchedDz, which are flat), the grown boxes span the full height of the domain, so
 * that each column is built from the bottom of the domain up, as on level 0.  The base state
 * then does not depend on where the level's boxes start: a box that starts above the bottom
 * would otherwise start its column from the coarse interpolation below it, and boxes that
 * start at different heights would build the same height differently (erf-model/ERF#5308).
 * Over terrain (MeshType::VariableDz) the heights of the level are known only on its grids,
 * and the grown boxes keep the vertical range of their box.
 *
 * The lateral ghost cells of bs_g, z_nd_g and z_cc_g are there for constructions that also
 * write the ring around the cells they are given; what lands in them is not to be used.
 * bs_g's cells that bs does not cover take the values of the nearest cell it does; where
 * z_phys_nd[lev] does not reach z_nd_g's they start from the heights of the flat mesh.
 *
 * @param[in]  lev    level of bs
 * @param[in]  bs     base state to copy, valid and ghost cells
 * @param[in]  grow   lateral growth of the boxes, at most bs's lateral ghost cells
 * @param[in]  ng     ghost cells of bs_g, z_nd_g and z_cc_g; ng[2] at most bs's vertical ones
 * @param[out] bs_g   the copy
 * @param[out] z_nd_g z_phys_nd[lev] on the grown boxes (nodal)
 * @param[out] z_cc_g the cell-center heights on the grown boxes, as make_zcc forms them
 */
void
ERF::make_base_state_on_grown_boxes (int lev, const MultiFab& bs, int grow, const IntVect& ng,
                                     MultiFab& bs_g,
                                     std::unique_ptr<MultiFab>& z_nd_g,
                                     std::unique_ptr<MultiFab>& z_cc_g)
{
    const IntVect ngb = bs.nGrowVect();
    AMREX_ALWAYS_ASSERT(grow <= ngb[0] && grow <= ngb[1] && ng[2] <= ngb[2]);
    AMREX_ALWAYS_ASSERT(z_phys_nd[lev] != nullptr);
    AMREX_ALWAYS_ASSERT(z_phys_nd[lev]->nGrowVect().allGE(IntVect(grow, grow, ng[2])));

    // The boxes grown inside the domain only: outside it the heights need not be set yet and
    //    the boundary conditions fill the base state, and across a periodic boundary the
    //    columns are the periodic images of columns inside it (which not every construction
    //    integrates where it lies outside the domain)
    const Box& dom_g = geom[lev].Domain();
    const bool full_columns = (solverChoice.mesh_type != MeshType::VariableDz);
    BoxList bl;
    for (int n = 0; n < static_cast<int>(bs.boxArray().size()); ++n) {
        Box b = bs.boxArray()[n];
        b.grow(0, grow); b.grow(1, grow);
        if (full_columns) {
            b.setRange(2, dom_g.smallEnd(2), dom_g.length(2));
        } else {
            b.setRange(2, bs.boxArray()[n].smallEnd(2), bs.boxArray()[n].length(2));
        }
        bl.push_back(b & dom_g);
    }
    const BoxArray ba_g(std::move(bl));
    const DistributionMapping& dm = bs.DistributionMap();

    //
    // Each copy takes the source's ghost cells first and then its valid cells, so that a cell
    //    some box of the source holds as valid gets that value: ParallelCopy does not order
    //    overlapping sources, and a ghost cell (or a ghost node, where nodal boxes share a
    //    face) need not agree with the valid data it overlaps.
    //
    bs_g.define(ba_g, dm, bs.nComp(), ng);
    bs_g.ParallelCopy(bs, 0, 0, bs.nComp(), ngb, ng);
    bs_g.ParallelCopy(bs, 0, 0, bs.nComp(), IntVect(0), ng);

    // The lateral ghost cells of bs_g lie beyond bs's ghost cells, and in a full column the
    //    cells above and below the box's own vertical ghost cells are not in bs either: give
    //    them the values of the nearest cell that is, so that a construction that reads them
    //    (e.g. the cell below a box that starts above the bottom) sees a value of the right size
    for (MFIter mfi(bs_g); mfi.isValid(); ++mfi) {
        const Box vb = mfi.validbox();
        const Box sb = bs.boxArray()[mfi.index()];
        const auto lo = lbound(vb);
        const auto hi = ubound(vb);
        const int klo_s = amrex::max(sb.smallEnd(2) - ngb[2], vb.smallEnd(2) - ng[2]);
        const int khi_s = amrex::min(sb.bigEnd(2)   + ngb[2], vb.bigEnd(2)   + ng[2]);
        const auto& b = bs_g.array(mfi);
        ParallelFor(mfi.fabbox(), bs_g.nComp(), [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int ii = amrex::min(amrex::max(i, lo.x), hi.x);
            const int jj = amrex::min(amrex::max(j, lo.y), hi.y);
            const int kk = amrex::min(amrex::max(k, klo_s), khi_s);
            if (ii != i || jj != j || kk != k) { b(i,j,k,n) = b(ii,jj,kk,n); }
        });
    }

    BoxArray ba_nd(ba_g);
    ba_nd.surroundingNodes();
    z_nd_g = std::make_unique<MultiFab>(ba_nd, dm, 1, ng);
    {
        // The heights of the flat mesh: uniform, or the stretched levels (extended past
        //    either end by the end spacing)
        const Real z_lo = geom[lev].ProbLo(2);
        const Real dz   = geom[lev].CellSize(2);
        const bool stretched = (solverChoice.mesh_type == MeshType::StretchedDz);
        const int  nzs = stretched ? static_cast<int>(zlevels_stag[lev].size()) : 0;
        Gpu::DeviceVector<Real> d_zs(nzs);
        if (stretched) {
            AMREX_ALWAYS_ASSERT(nzs >= 2);
            Gpu::copy(Gpu::hostToDevice, zlevels_stag[lev].begin(), zlevels_stag[lev].end(), d_zs.begin());
        }
        const Real* zs = d_zs.data();
        for (MFIter mfi(*z_nd_g); mfi.isValid(); ++mfi) {
            const auto& z = z_nd_g->array(mfi);
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
    z_nd_g->ParallelCopy(*z_phys_nd[lev], 0, 0, 1, z_phys_nd[lev]->nGrowVect(), ng);
    z_nd_g->ParallelCopy(*z_phys_nd[lev], 0, 0, 1, IntVect(0), ng);

    z_cc_g = std::make_unique<MultiFab>(ba_g, dm, 1, ng);
    make_zcc(geom[lev], *z_nd_g, *z_cc_g);
}

/**
 * Copy a base state built by make_base_state_on_grown_boxes back into bs.
 *
 * Each box's valid cells (and the cells above and below them) come from that box's own grown
 * box, and the ghost cells that lie in another box's valid region from that box: grown boxes
 * that start at different heights integrate a column they share from different cells, so they
 * need not agree there.  The other lateral ghost cells inside the domain come from the grown
 * boxes (across a periodic boundary, from the periodic image), and those outside it from the
 * nearest cell inside it, until the boundary conditions fill them.
 *
 * @param[in]     lev  level of bs
 * @param[in,out] bs   base state, with the BoxArray and DistributionMapping bs_g was made from
 * @param[in]     bs_g the base state built on the grown boxes
 */
void
ERF::copy_base_state_from_grown_boxes (int lev, MultiFab& bs, const MultiFab& bs_g)
{
    const IntVect ngb = bs.nGrowVect();
    const int     ngz = std::min(ngb[2], bs_g.nGrowVect()[2]);
    const int     nc  = bs.nComp();

    bs.ParallelCopy(bs_g, 0, 0, nc, IntVect(0,0,ngz), ngb, geom[lev].periodicity());

    for (MFIter mfi(bs); mfi.isValid(); ++mfi) {
        const Box vb = amrex::grow(mfi.validbox(), IntVect(0,0,ngz));
        const auto& d = bs.array(mfi);
        const auto& g = bs_g.const_array(mfi);
        ParallelFor(vb, nc, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            d(i,j,k,n) = g(i,j,k,n);
        });
    }

    bs.FillBoundary(geom[lev].periodicity());

    Box dom_g(geom[lev].Domain());
    for (int dir = 0; dir < 2; ++dir) {
        if (geom[lev].isPeriodic(dir)) { dom_g.grow(dir, ngb[dir]); }
    }
    const auto dlo = lbound(dom_g);
    const auto dhi = ubound(dom_g);
    for (MFIter mfi(bs); mfi.isValid(); ++mfi) {
        const auto& d = bs.array(mfi);
        ParallelFor(mfi.fabbox(), nc, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int ii = amrex::min(amrex::max(i, dlo.x), dhi.x);
            const int jj = amrex::min(amrex::max(j, dlo.y), dhi.y);
            if (ii != i || jj != j) { d(i,j,k,n) = d(ii,jj,k,n); }
        });
    }
}

/**
 * Rebuild base_state[lev] after a regrid made or remade level lev (> 0).
 *
 * A regrid fills a new or remade level's base state by conservative interpolation from the
 * coarser level.  Next to the bottom and top walls that interpolation has no vertical slope,
 * so p_0 in the two cells nearest each wall is out of discrete hydrostatic balance by about
 * rho_0 g dz, and the c/f ghost cells of a level whose grids changed get those values too.
 * The remade level then no longer matches one built at the same place from scratch.  Here
 * the base state is built again on the new grids, by the same construction as at
 * initialization wherever the run can repeat it, and from the coarser level otherwise.
 *
 * Called after base_state[lev] has been defined on the new grids and the terrain arrays and
 * physbcs of the level have been made.
 *
 * @param[in] lev level that was made or remade
 */
void
ERF::rebuild_base_state (int lev)
{
    AMREX_ALWAYS_ASSERT(lev > 0);

    if (solverChoice.init_type == InitType::Input_Sounding &&
        solverChoice.sounding_type != SoundingType::ConstantDensity)
    {
        // The same steps as init_from_input_sounding for a fine level, so that a remade level
        // has the base state a level built from scratch on the same grids has
        base_state[lev-1].FillBoundary(geom[lev-1].periodicity());
        InterpFromCoarseLevel(base_state[lev], base_state[lev].nGrowVect(),
                              IntVect(0,0,0), // do not fill ghost cells outside the domain
                              base_state[lev-1], 0, 0, base_state[lev].nComp(),
                              geom[lev-1], geom[lev],
                              refRatio(lev-1), &cell_cons_interp,
                              domain_bcs_type, BaseBCVars::rho0_bc_comp);
        (*physbcs_base[lev])(base_state[lev],0,base_state[lev].nComp(),base_state[lev].nGrowVect());

        rebuild_base_state_from_input_sounding(lev, base_state[lev]);
        (*physbcs_base[lev])(base_state[lev],0,base_state[lev].nComp(),base_state[lev].nGrowVect());
    }
    else if (solverChoice.init_type == InitType::WRFInput)
    {
        // RemakeLevel and MakeNewLevelFromCoarse have already rebuilt it with
        // rebuild_base_state_from_wrfinput, or read it from this level's own file
    }
    else if (solverChoice.init_type == InitType::Metgrid ||
             solverChoice.init_type == InitType::NCFile)
    {
        // Only the coarse level holds the data these were initialized from
        rebuild_base_state_from_coarse(lev);
    }
    else
    {
        // Uniform, ConstantDensity, ConstantDensityLinearTheta, Isentropic, MoistBaseState,
        // HindCast, and an input sounding of SoundingType::ConstantDensity: the problem's own
        // hydrostatic construction, as at initialization
        initHSE(lev);
    }

    base_state[lev].FillBoundary(geom[lev].periodicity());

    check_base_state_hse(lev, "after a regrid");
}

/**
 * Rebuild base_state[lev] from the coarser level alone, in discrete hydrostatic balance.
 *
 * theta_0 and qv_0 are interpolated from the coarser level; they vary smoothly and are kept.
 * In each column of the valid region and of the lateral ghost cells, the lowest cell
 * of the level is then given the density that puts it in hydrostatic balance with the coarse
 * cell that contains it, integrating from the height of that coarse cell's center to its own,
 * and the column is integrated upward from there with rebalance_columns.  p_0 and pi_0 follow
 * from the equation of state.  Working in physical heights, this also holds where the fine
 * terrain differs from the coarse terrain.
 *
 * @param[in] lev level to rebuild (> 0)
 */
void
ERF::rebuild_base_state_from_coarse (int lev)
{
    AMREX_ALWAYS_ASSERT(lev > 0);

    // Everything the construction below does not overwrite: ghost cells below a level that
    // starts above the bottom, and outside the domain
    interp_base_state_from_coarse(lev);

    if (!solverChoice.use_gravity) { return; }

    const BoxArray& ba = base_state[lev].boxArray();
    const DistributionMapping& dm = base_state[lev].DistributionMap();
    const Box& domain = geom[lev].Domain();

    // The column integration continues a box from the box below it only through
    // rebalance_columns' bands, which do not extend to the lateral ghost cells built here
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(column_bands(ba).size() == 1,
        "Rebuilding the base state of a regridded level from the coarser level requires "
        "every grid of the level to span the level's full height (no split in z)");

    // The valid region and the lateral ghost cells inside the domain, as grids of their own,
    //    so that each ghost column is built as it would be were it valid
    const IntVect ngb = base_state[lev].nGrowVect();
    BoxList bl;
    for (int n = 0; n < static_cast<int>(ba.size()); ++n) {
        Box b = ba[n];
        b.grow(0,ngb[0]); b.grow(1,ngb[1]);
        b &= domain;
        bl.push_back(b);
    }
    const BoxArray gba(std::move(bl));
    BoxArray gba_nd(gba); gba_nd.surroundingNodes();

    const int ncomp = BaseState::num_comps;
    MultiFab work(gba, dm, ncomp, 0);
    work.ParallelCopy(base_state[lev], 0, 0, ncomp, base_state[lev].nGrowVect(), IntVect(0));

    MultiFab z_nd(gba_nd, dm, 1, 0);
    z_nd.ParallelCopy(*z_phys_nd[lev], 0, 0, 1, z_phys_nd[lev]->nGrowVect(), IntVect(0));

    // The coarse cell containing each fine cell: p_0, theta_0, qv_0 and the height of its center
    enum { Pc = 0, Thc, Qvc, Zc, n_crse };
    MultiFab crse(base_state[lev-1].boxArray(), base_state[lev-1].DistributionMap(), n_crse, 1);
    base_state[lev-1].FillBoundary(geom[lev-1].periodicity());
    MultiFab::Copy(crse, base_state[lev-1], BaseState::p0_comp,  Pc,  1, 1);
    MultiFab::Copy(crse, base_state[lev-1], BaseState::th0_comp, Thc, 1, 1);
    MultiFab::Copy(crse, base_state[lev-1], BaseState::qv0_comp, Qvc, 1, 1);
    MultiFab::Copy(crse, *z_phys_cc[lev-1], 0, Zc, 1, 1);

    MultiFab seed(gba, dm, n_crse, 0);
    Vector<BCRec> bcs(n_crse, domain_bcs_type[BCVars::cons_bc]);
    InterpFromCoarseLevel(seed, IntVect(0), IntVect(0), crse, 0, 0, n_crse,
                          geom[lev-1], geom[lev], refRatio(lev-1), &pc_interp, bcs, 0);

#ifdef AMREX_USE_FLOAT
    const Real tol = Real(1.0e-6);
#else
    const Real tol = Real(1.0e-10);
#endif
    const Real grav = CONST_GRAV;

    for (MFIter mfi(work); mfi.isValid(); ++mfi)
    {
        const int klo = mfi.validbox().smallEnd(2);
        Box slab = mfi.validbox(); slab.makeSlab(2, klo);

        const Array4<Real      >& w   = work.array(mfi);
        const Array4<Real const>& s   = seed.const_array(mfi);
        const Array4<Real const>& zna = z_nd.const_array(mfi);

        ParallelFor(slab, [=,RdoCp_d=RdoCp] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // From the coarse cell center ...
            const Real P_lo  = s(i,j,k,Pc);
            const Real Th_lo = s(i,j,k,Thc);
            const Real qv_lo = s(i,j,k,Qvc);
            const Real z_lo  = s(i,j,k,Zc);

            // ... to the center of this fine cell, at the height rebalance_columns uses
            const Real z_hi = Real(0.125) * (zna(i,j,k  ) + zna(i+1,j,k  ) + zna(i,j+1,k  ) + zna(i+1,j+1,k  )
                                            +zna(i,j,k+1) + zna(i+1,j,k+1) + zna(i,j+1,k+1) + zna(i+1,j+1,k+1));
            const Real dz = z_hi - z_lo;

            const Real R_lo = getRhogivenThetaPress(Th_lo, P_lo, RdoCp_d, qv_lo);
            const Real C    = -P_lo + myhalf * R_lo * (one + qv_lo) * grav * dz;

            const Real Th_hi = w(i,j,k,BaseState::th0_comp);
            const Real qv_hi = w(i,j,k,BaseState::qv0_comp);
            Real P_hi = P_lo;
            Real T_hi = getTgivenPandTh(P_hi, Th_hi, RdoCp_d);
            Real R_hi = getRhogivenThetaPress(Th_hi, P_hi, RdoCp_d, qv_hi);
            Real F    = P_hi + myhalf * R_hi * (one + qv_hi) * grav * dz + C;

            HSEutils::Newton_Raphson_hse(tol, RdoCp_d, dz, grav, C, Th_hi, T_hi,
                                         qv_hi, qv_hi, P_hi, R_hi, F, true);

            w(i,j,k,BaseState::r0_comp) = R_hi;
        });
    }

    // Integrate each column upward from that cell, holding theta_0 and qv_0
    MultiFab r_w (work, make_alias, BaseState::r0_comp,  1);
    MultiFab th_w(work, make_alias, BaseState::th0_comp, 1);
    MultiFab qv_w(work, make_alias, BaseState::qv0_comp, 1);
    const bool maintain_Th = true;
    rebalance_columns(r_w, th_w, qv_w, qv_w, &z_nd, geom[lev], maintain_Th);

    recompute_p_pi_hse(work, solverChoice.rdOcp);

    // Back into the valid region and the lateral ghost cells
    base_state[lev].ParallelCopy(work, 0, 0, ncomp, IntVect(0), IntVect(ngb[0],ngb[1],0));

    (*physbcs_base[lev])(base_state[lev],0,base_state[lev].nComp(),base_state[lev].nGrowVect());
}

/**
 * Check that base_state[lev] is in discrete hydrostatic balance, and abort (or warn) if it
 * is not.
 *
 * The residual is dp_0/dz + rho_0 (1 + qv_0) |g| between vertically adjacent cells of the
 * same grid, with dz the distance between the cell centers (from z_phys_nd), which is the balance
 * rebalance_columns and the initializations enforce.  It is checked over the valid region,
 * as compute_max_pressure_gradient_diagnostic does.  The tolerance is the one
 * compute_max_pressure_gradient_diagnostic uses.
 *
 * Not checked: runs without gravity, EB (covered cells hold no meaningful base state), and
 * the init types whose base state is not built in hydrostatic balance (Uniform and
 * ConstantDensityLinearTheta, which skip erf_enforce_hse).
 *
 * A restart only warns: a checkpoint written before the base state of a refined level was
 * rebuilt in balance (erf-model/ERF#5228, #5236) holds one that is not, and the run that
 * wrote it should still be able to continue.
 *
 * @param[in] lev           level to check
 * @param[in] when          when the check is made, for the message
 * @param[in] abort_if_fail abort if out of balance (else print a warning)
 */
void
ERF::check_base_state_hse (int lev, const std::string& when, bool abort_if_fail)
{
    if (!solverChoice.use_gravity) { return; }
    if (solverChoice.terrain_type == TerrainType::EB) { return; }
    if (solverChoice.init_type == InitType::Uniform ||
        solverChoice.init_type == InitType::ConstantDensityLinearTheta) { return; }

    const MultiFab& bs = base_state[lev];
    const MultiFab& znd = *z_phys_nd[lev];
    const Real grav    = CONST_GRAV;

    // The residual on each face, kept per cell so that its largest value can be located
    MultiFab res(bs.boxArray(), bs.DistributionMap(), 1, 0);
    res.setVal(zero);

    for (MFIter mfi(bs); mfi.isValid(); ++mfi)
    {
        // Faces between cells k-1 and k of this grid
        Box bx = mfi.validbox();
        bx.growLo(2,-1);
        if (!bx.ok()) { continue; }

        const Array4<Real const>& b = bs.const_array(mfi);
        const Array4<Real const>& zn = znd.const_array(mfi);
        const Array4<Real      >& r = res.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Cell-center heights as rebalance_columns forms them, from the nodal heights
            auto zcc = [&] (int kk) {
                return Real(0.125) * (zn(i,j,kk  ) + zn(i+1,j,kk  ) + zn(i,j+1,kk  ) + zn(i+1,j+1,kk  )
                                     +zn(i,j,kk+1) + zn(i+1,j,kk+1) + zn(i,j+1,kk+1) + zn(i+1,j+1,kk+1));
            };
            const Real dz  = zcc(k) - zcc(k-1);
            const Real rhi = b(i,j,k  ,BaseState::r0_comp) * (one + b(i,j,k  ,BaseState::qv0_comp));
            const Real rlo = b(i,j,k-1,BaseState::r0_comp) * (one + b(i,j,k-1,BaseState::qv0_comp));
            r(i,j,k) = std::abs((b(i,j,k,BaseState::p0_comp) - b(i,j,k-1,BaseState::p0_comp)) / dz
                                + myhalf * grav * (rhi + rlo));
        });
    }

    const Real max_res   = res.max(0, 0);
    const IntVect loc    = res.maxIndex(0);

    // As in compute_max_pressure_gradient_diagnostic: a fixed bound, or the roundoff in p_0
    // spread over one vertical cell where that is larger (single precision)
    const Real tol_fixed    = Real(1.e-8);
    const Real p_max        = bs.max(BaseState::p0_comp);
    const Real dz_for_tol   = (lev < static_cast<int>(dz_min.size()) && dz_min[lev] > Real(0))
                            ? dz_min[lev] : geom[lev].CellSize(2);
    const Real tol_roundoff = Real(16) * std::numeric_limits<Real>::epsilon() * p_max / dz_for_tol;
    const Real tol          = std::max(tol_fixed, tol_roundoff);

    if (verbose > 0) {
        Print() << "Base state at level " << lev << " " << when
                << ": max |dp0/dz + rho0*|g|| = " << max_res << " at face " << loc
                << " (tolerance " << tol << ")" << std::endl;
    }
    if (!(max_res <= tol)) {
        std::ostringstream msg;
        msg << std::setprecision(6) << "The base state at level " << lev
            << " is out of hydrostatic balance " << when << ": max |dp0/dz + rho0*|g|| = "
            << max_res << " at face " << loc << " exceeds " << tol;
        if (abort_if_fail) {
            Abort(msg.str());
        } else {
            Print() << "WARNING: " << msg.str() << ".  The checkpoint may have been written "
                    << "before refined levels' base states were rebuilt in balance; the level "
                    << "keeps this base state until a regrid remakes it." << std::endl;
        }
    }
}
