/**
 * \file ERF_MLTerrainPoisson.cpp
 */
#include "ERF_MLTerrainPoisson.H"
#include "ERF_Utils.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelReduce.H>

#include <limits>

using namespace amrex;


MLTerrainPoisson::MLTerrainPoisson (const Vector<Geometry>& a_geom,
                                    const Vector<BoxArray>& a_grids,
                                    const Vector<DistributionMapping>& a_dmap,
                                    const LPInfo& a_info,
                                    const MultiFab& z_phys_nd,
                                    const MultiFab& ax,
                                    const MultiFab& ay,
                                    const MultiFab& az,
                                    const MultiFab& dJ)
{
    BL_PROFILE("MLTerrainPoisson::MLTerrainPoisson()");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_geom.size() == 1 && a_grids.size() == 1 && a_dmap.size() == 1,
                                     "MLTerrainPoisson solves one AMR level at a time");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_info.hidden_direction != 2,
                                     "MLTerrainPoisson: the vertical direction cannot be hidden");

    define(a_geom, a_grids, a_dmap, a_info, {});

    define_levels(z_phys_nd, ax, ay, az, dJ);
}

int
MLTerrainPoisson::hidden_index (int mglev) const
{
    const int hd = hiddenDirection();
    return (hd >= 0) ? m_geom[0][mglev].Domain().smallEnd(hd) : 0;
}

Dim3
MLTerrainPoisson::periodic_lengths (int mglev) const
{
    const Geometry& geom = m_geom[0][mglev];
    const IntVect len = geom.Domain().length();
    Dim3 n{0,0,0};
    if (geom.isPeriodic(0)) { n.x = len[0]; }
    if (geom.isPeriodic(1)) { n.y = len[1]; }
    if (geom.isPeriodic(2)) { n.z = len[2]; }
    return n;
}

TerrainGhostRule
MLTerrainPoisson::ghost_rule (int mglev) const
{
    TerrainGhostRule r;
    const Box& dom = m_geom[0][mglev].Domain();
    r.dlo = lbound(dom);
    r.dhi = ubound(dom);
    // A box face inside the domain is a coarse/fine (or box-union) boundary
    r.cf_sign = (m_coarse_fine_bc_type == LinOpBCType::Dirichlet) ? -1 : 1;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        auto lo = m_lobc[0][dir];
        auto hi = m_hibc[0][dir];
        // Beyond a periodic domain face the only uncovered ghost cells are those of a
        // box union that does not span the direction, i.e. coarse/fine ones
        r.lo_sign[dir] = (lo == LinOpBCType::Dirichlet) ? -1 : ((lo == LinOpBCType::Neumann) ? 1 : r.cf_sign);
        r.hi_sign[dir] = (hi == LinOpBCType::Dirichlet) ? -1 : ((hi == LinOpBCType::Neumann) ? 1 : r.cf_sign);
    }
    return r;
}

IntVect
MLTerrainPoisson::coarsening_depth (const Box& b, int hidden_dir, int min_width)
{
    IntVect depth(0);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (d == hidden_dir) { continue; }
        IntVect ratio(1); ratio[d] = 2;
        IntVect minw(0);  minw[d]  = min_width;
        Box c = b;
        while (c.coarsenable(ratio, minw)) {
            c.coarsen(ratio);
            ++depth[d];
        }
    }
    return depth;
}

IntVect
MLTerrainPoisson::coarsening_depth (const BoxArray& ba, int hidden_dir, int min_width)
{
    IntVect depth(0);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (d == hidden_dir) { continue; }
        IntVect ratio(1); ratio[d] = 2;
        IntVect minw(0);  minw[d]  = min_width;
        BoxArray c = ba;
        while (c.coarsenable(ratio, minw)) {
            c.coarsen(ratio);
            ++depth[d];
        }
    }
    return depth;
}

BoxArray
MLTerrainPoisson::multigrid_grids (const Box& domain, int hidden_dir, const IntVect& depth, int max_size)
{
    const IntVect dom_depth = coarsening_depth(domain, hidden_dir);
    AMREX_ALWAYS_ASSERT(depth.allGE(IntVect(0)) && depth.allLE(dom_depth));
    AMREX_ALWAYS_ASSERT(max_size > 0);

    // Segment each direction into pieces of n * 2^depth cells
    Array<Vector<int>,AMREX_SPACEDIM> seg_lo, seg_len;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const int L = domain.length(d);
        if (d == hidden_dir) {
            seg_lo[d].push_back(domain.smallEnd(d));
            seg_len[d].push_back(L);
            continue;
        }
        const int unit = 1 << depth[d];
        const int q = L / unit;                               // exact by the assertion above
        // A segment must keep at least two cells after depth[d] coarsenings, i.e. two units
        int nseg = (L + max_size - 1) / max_size;
        nseg = amrex::max(1, amrex::min(nseg, q / 2));
        const int base = q / nseg;
        const int rem  = q % nseg;
        int lo = domain.smallEnd(d);
        for (int n = 0; n < nseg; ++n) {
            const int len = (base + ((n < rem) ? 1 : 0)) * unit;
            seg_lo[d].push_back(lo);
            seg_len[d].push_back(len);
            lo += len;
        }
    }

    BoxList bl;
    for (int kseg = 0; kseg < static_cast<int>(seg_lo[2].size()); ++kseg) {
    for (int jseg = 0; jseg < static_cast<int>(seg_lo[1].size()); ++jseg) {
    for (int iseg = 0; iseg < static_cast<int>(seg_lo[0].size()); ++iseg) {
        IntVect lo(seg_lo[0][iseg], seg_lo[1][jseg], seg_lo[2][kseg]);
        IntVect hi(lo[0] + seg_len[0][iseg] - 1, lo[1] + seg_len[1][jseg] - 1, lo[2] + seg_len[2][kseg] - 1);
        bl.push_back(Box(lo, hi));
    }}}
    return BoxArray(std::move(bl));
}

void
MLTerrainPoisson::subsample (const MultiFab& fine, MultiFab& crse, const IntVect& ratio)
{
    BoxArray fba_c = amrex::coarsen(fine.boxArray(), ratio);
    const bool direct = (fba_c == crse.boxArray()) && (fine.DistributionMap() == crse.DistributionMap());

    MultiFab tmp;
    MultiFab* dst = &crse;
    if (!direct) {
        tmp.define(fba_c, fine.DistributionMap(), 1, 0);
        dst = &tmp;
    }

    const Dim3 r = ratio.dim3();
    for (MFIter mfi(*dst); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        Array4<Real const> const& f = fine.const_array(mfi);
        Array4<Real      > const& c = dst->array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            c(i,j,k) = f(i*r.x, j*r.y, k*r.z);
        });
    }

    if (!direct) {
        crse.ParallelCopy(tmp, 0, 0, 1, IntVect(0), IntVect(0));
    }
}

void
MLTerrainPoisson::fill_uncovered_nodes (MultiFab& znd, const Geometry& geom)
{
    AMREX_ALWAYS_ASSERT(znd.ixType().nodeCentered());
    const IntVect ng = znd.nGrowVect();

    // 1 on the ghost nodes no box (or periodic image of a box) holds
    iMultiFab nmask(znd.boxArray(), znd.DistributionMap(), 1, ng);
    nmask.setVal(1);
    nmask.setVal(0, 0, 1, 0);
    nmask.FillBoundary(geom.periodicity());

    for (MFIter mfi(znd); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        const Box gbx = amrex::grow(vbx, ng);
        const auto lo = lbound(vbx);
        const auto hi = ubound(vbx);
        Array4<Real     > const& z = znd.array(mfi);
        Array4<int const> const& m = nmask.const_array(mfi);

        // Lateral x: copy the nearest node inside the box
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ((i < lo.x || i > hi.x) && j >= lo.y && j <= hi.y && k >= lo.z && k <= hi.z && m(i,j,k) == 1) {
                int ii = (i < lo.x) ? lo.x : hi.x;
                z(i,j,k) = z(ii,j,k);
            }
        });
        // Lateral y, over the x-extended slab
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ((j < lo.y || j > hi.y) && k >= lo.z && k <= hi.z && m(i,j,k) == 1) {
                int jj = (j < lo.y) ? lo.y : hi.y;
                z(i,j,k) = z(i,jj,k);
            }
        });
        // Vertical: linear extrapolation from the two nearest node layers
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ((k < lo.z || k > hi.z) && m(i,j,k) == 1) {
                if (k < lo.z) {
                    // ERF writes the node below the surface as 2 z0 - z1
                    z(i,j,k) = (k == lo.z-1) ? Real(2.0)*z(i,j,lo.z) - z(i,j,lo.z+1)
                             : z(i,j,lo.z) + static_cast<Real>(lo.z - k) * (z(i,j,lo.z) - z(i,j,lo.z+1));
                } else {
                    z(i,j,k) = z(i,j,hi.z) + static_cast<Real>(k - hi.z) * (z(i,j,hi.z) - z(i,j,hi.z-1));
                }
            }
        });
    }
}

void
MLTerrainPoisson::define_levels (const MultiFab& z_phys_nd,
                                 const MultiFab& ax, const MultiFab& ay,
                                 const MultiFab& az, const MultiFab& dJ)
{
    BL_PROFILE("MLTerrainPoisson::define_levels()");

    const int nmg = NMGLevels(0);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(z_phys_nd.nGrowVect().allGE(IntVect(1)),
                                     "MLTerrainPoisson: z_phys_nd needs at least one ghost node");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(z_phys_nd.boxArray() == amrex::convert(m_grids[0][0], IntVect(1)) &&
                                     z_phys_nd.DistributionMap() == m_dmap[0][0],
                                     "MLTerrainPoisson: z_phys_nd must live on the nodal version of the solve grids");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ax.boxArray() == amrex::convert(m_grids[0][0], IntVect(1,0,0)) &&
                                     ay.boxArray() == amrex::convert(m_grids[0][0], IntVect(0,1,0)) &&
                                     az.boxArray() == amrex::convert(m_grids[0][0], IntVect(0,0,1)) &&
                                     dJ.boxArray() == m_grids[0][0],
                                     "MLTerrainPoisson: ax, ay, az, dJ must live on the solve grids");

    m_zphys.resize(nmg);
    m_ax.resize(nmg);
    m_ay.resize(nmg);
    m_az.resize(nmg);
    m_dJ.resize(nmg);
    m_tri.resize(nmg);
    m_mask.resize(nmg);
    m_src.resize(nmg);
    m_sgn.resize(nmg);
    m_box_has_cf.resize(nmg);
    m_sx.resize(nmg);
    m_sy.resize(nmg);
    m_sz.resize(nmg);

    // The finest level is the one the projection hands us
    m_zphys[0] = MultiFab(z_phys_nd, make_alias, 0, 1);
    m_ax[0]    = MultiFab(ax, make_alias, 0, 1);
    m_ay[0]    = MultiFab(ay, make_alias, 0, 1);
    m_az[0]    = MultiFab(az, make_alias, 0, 1);
    m_dJ[0]    = MultiFab(dJ, make_alias, 0, 1);

    for (int mglev = 0; mglev < nmg; ++mglev) {
        m_mask[mglev].define(m_grids[0][mglev], m_dmap[0][mglev], 1, 1);
        m_mask[mglev].setVal(1);
        m_mask[mglev].setVal(0, 0, 1, 0);
        m_mask[mglev].FillBoundary(m_geom[0][mglev].periodicity());
    }
    // The ghost source maps depend on the boundary conditions, which are set after
    // construction, so prepareForSolve builds them

    if (nmg == 1) { return; }

    //
    // The projection scales ax, ay and az by map factors before it builds the operator.
    // Recover that scale relative to the plain metric areas so the coarse levels can
    // carry it too; when it is identically one (no map factors) skip it altogether.
    //
    {
        const Geometry& g0 = m_geom[0][0];
        MultiFab ax_raw(ax.boxArray(), m_dmap[0][0], 1, 1);
        MultiFab ay_raw(ay.boxArray(), m_dmap[0][0], 1, 1);
        MultiFab az_raw(az.boxArray(), m_dmap[0][0], 1, 1);
        ax_raw.setVal(1.0); ay_raw.setVal(1.0); az_raw.setVal(1.0);
        make_areas(g0, m_zphys[0], ax_raw, ay_raw, az_raw);

        m_sx[0].define(ax.boxArray(), m_dmap[0][0], 1, 0);
        m_sy[0].define(ay.boxArray(), m_dmap[0][0], 1, 0);
        m_sz[0].define(az.boxArray(), m_dmap[0][0], 1, 0);
        MultiFab::Copy(m_sx[0], ax, 0, 0, 1, 0); MultiFab::Divide(m_sx[0], ax_raw, 0, 0, 1, 0);
        MultiFab::Copy(m_sy[0], ay, 0, 0, 1, 0); MultiFab::Divide(m_sy[0], ay_raw, 0, 0, 1, 0);
        MultiFab::Copy(m_sz[0], az, 0, 0, 1, 0); MultiFab::Divide(m_sz[0], az_raw, 0, 0, 1, 0);

        Real dev = Real(0.0);
        for (MultiFab* s : {&m_sx[0], &m_sy[0], &m_sz[0]}) {
            MultiFab tmp(s->boxArray(), s->DistributionMap(), 1, 0);
            MultiFab::Copy(tmp, *s, 0, 0, 1, 0);
            tmp.plus(Real(-1.0), 0, 1, 0);
            dev = amrex::max(dev, tmp.norm0(0, 0, false));
        }
        // Map factors differ from one by far more than round-off when they are present
        m_has_scale = (dev > Real(100.0) * std::numeric_limits<Real>::epsilon());
    }

    for (int mglev = 1; mglev < nmg; ++mglev)
    {
        const Geometry& gc = m_geom[0][mglev];
        const Geometry& gf = m_geom[0][mglev-1];
        const IntVect ratio = gf.Domain().length() / gc.Domain().length();
        const BoxArray& cba = m_grids[0][mglev];
        const DistributionMapping& cdm = m_dmap[0][mglev];

        // The coarse surface is the fine surface at every ratio-th node
        m_zphys[mglev].define(amrex::convert(cba, IntVect(1)), cdm, 1, 1);
        subsample(m_zphys[mglev-1], m_zphys[mglev], ratio);
        m_zphys[mglev].FillBoundary(gc.periodicity());
        fill_uncovered_nodes(m_zphys[mglev], gc);

        m_dJ[mglev].define(cba, cdm, 1, 1);
        m_ax[mglev].define(amrex::convert(cba, IntVect(1,0,0)), cdm, 1, 1);
        m_ay[mglev].define(amrex::convert(cba, IntVect(0,1,0)), cdm, 1, 1);
        m_az[mglev].define(amrex::convert(cba, IntVect(0,0,1)), cdm, 1, 1);
        // The metric routines leave the ghost layer below the surface alone
        m_dJ[mglev].setVal(1.0);
        m_ax[mglev].setVal(1.0);
        m_ay[mglev].setVal(1.0);
        m_az[mglev].setVal(1.0);

        make_J(gc, m_zphys[mglev], m_dJ[mglev]);
        make_areas(gc, m_zphys[mglev], m_ax[mglev], m_ay[mglev], m_az[mglev]);

        if (m_has_scale) {
            m_sx[mglev].define(m_ax[mglev].boxArray(), cdm, 1, 0);
            m_sy[mglev].define(m_ay[mglev].boxArray(), cdm, 1, 0);
            m_sz[mglev].define(m_az[mglev].boxArray(), cdm, 1, 0);
            subsample(m_sx[mglev-1], m_sx[mglev], ratio);
            subsample(m_sy[mglev-1], m_sy[mglev], ratio);
            subsample(m_sz[mglev-1], m_sz[mglev], ratio);
            MultiFab::Multiply(m_ax[mglev], m_sx[mglev], 0, 0, 1, 0);
            MultiFab::Multiply(m_ay[mglev], m_sy[mglev], 0, 0, 1, 0);
            MultiFab::Multiply(m_az[mglev], m_sz[mglev], 0, 0, 1, 0);
            m_ax[mglev].FillBoundary(gc.periodicity());
            m_ay[mglev].FillBoundary(gc.periodicity());
            m_az[mglev].FillBoundary(gc.periodicity());
        }
    }
}

void
MLTerrainPoisson::define_ghost_fill (int mglev)
{
    BL_PROFILE("MLTerrainPoisson::define_ghost_fill()");

    const TerrainGhostRule rule = ghost_rule(mglev);
    const int hd = hiddenDirection();

    m_src[mglev].define(m_grids[0][mglev], m_dmap[0][mglev], 3, 1);
    m_sgn[mglev].define(m_grids[0][mglev], m_dmap[0][mglev], 1, 1);
    m_src[mglev].setVal(0);
    m_sgn[mglev].setVal(0.0);
    m_box_has_cf[mglev].define(m_grids[0][mglev], m_dmap[0][mglev]);

    const Dim3 dlo = lbound(m_geom[0][mglev].Domain());
    const Dim3 dhi = ubound(m_geom[0][mglev].Domain());
    const Dim3 nper = periodic_lengths(mglev);
    iMultiFab cf(m_grids[0][mglev], m_dmap[0][mglev], 1, 1);
    cf.setVal(0);

    for (MFIter mfi(m_mask[mglev]); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        const Dim3 blo = lbound(vbx);
        const Dim3 bhi = ubound(vbx);
        // The MLMG vectors have no ghost cells in a hidden direction
        IntVect ng(1);
        if (hd >= 0) { ng[hd] = 0; }
        const Box gbx = amrex::grow(vbx, ng);
        Array4<int  const> const& mask = m_mask[mglev].const_array(mfi);
        Array4<int       > const& src  = m_src[mglev].array(mfi);
        Array4<Real      > const& sgn  = m_sgn[mglev].array(mfi);
        Array4<int       > const& cfa  = cf.array(mfi);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            bool inside = (i >= blo.x && i <= bhi.x && j >= blo.y && j <= bhi.y && k >= blo.z && k <= bhi.z);
            if (!inside && mask(i,j,k) == 1) {
                int si, sj, sk;
                sgn(i,j,k) = terrain_ghost_source(i, j, k, blo, bhi, mask, rule, si, sj, sk);
                src(i,j,k,0) = si;
                src(i,j,k,1) = sj;
                src(i,j,k,2) = sk;
                // only the face neighbours bound a face of the box
                int nout = ((i < blo.x || i > bhi.x) ? 1 : 0) + ((j < blo.y || j > bhi.y) ? 1 : 0) + ((k < blo.z || k > bhi.z) ? 1 : 0);
                if (nout == 1 && terrain_is_cf_cell(i, j, k, mask, dlo, dhi, nper)) { cfa(i,j,k) = 1; }
            }
        });
        m_box_has_cf[mglev][mfi] = (cf[mfi].template sum<RunOn::Device>(gbx, 0) > 0) ? 1 : 0;
    }
}

void
MLTerrainPoisson::compute_tridiagonal (int mglev)
{
    BL_PROFILE("MLTerrainPoisson::compute_tridiagonal()");

    const Geometry& geom = m_geom[0][mglev];
    const auto dxinv = geom.InvCellSizeArray();
    const Dim3 dlo = lbound(geom.Domain());
    const Dim3 dhi = ubound(geom.Domain());
    const Dim3 nper = periodic_lengths(mglev);
    const int hd = hiddenDirection();
    const int hidx = hidden_index(mglev);

    m_tri[mglev].define(m_grids[0][mglev], m_dmap[0][mglev], 3, 0);

    for (MFIter mfi(m_tri[mglev]); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        const Dim3 blo = lbound(bx);
        const Dim3 bhi = ubound(bx);
        Array4<Real      > const& tri  = m_tri[mglev].array(mfi);
        Array4<int  const> const& mask = m_mask[mglev].const_array(mfi);
        Array4<int  const> const& src  = m_src[mglev].const_array(mfi);
        Array4<Real const> const& sgn  = m_sgn[mglev].const_array(mfi);
        Array4<Real const> const& axa  = m_ax[mglev].const_array(mfi);
        Array4<Real const> const& aya  = m_ay[mglev].const_array(mfi);
        Array4<Real const> const& aza  = m_az[mglev].const_array(mfi);
        Array4<Real const> const& dJa  = m_dJ[mglev].const_array(mfi);
        Array4<Real const> const& zpa  = m_zphys[mglev].const_array(mfi);
        const bool has_cf = (m_box_has_cf[mglev][mfi] != 0);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for (int d = -1; d <= 1; ++d) {
                TerrainUnitVector e{i, j, k+d, blo, bhi, mask, src, sgn, dlo, dhi, nper, hd, hidx};
                tri(i,j,k,d+1) = has_cf ? terrain_adotx_cf(i, j, k, e, mask, dlo, dhi, nper, hd, axa, aya, aza, dJa, zpa,
                                                            dxinv[0], dxinv[1], dxinv[2])
                                        : terrpoisson_adotx_value(i, j, k, e, axa, aya, aza, dJa, zpa,
                                                                  dxinv[0], dxinv[1], dxinv[2]);
            }
        });
    }
}

void
MLTerrainPoisson::prepareForSolve ()
{
    BL_PROFILE("MLTerrainPoisson::prepareForSolve()");

    MLCellLinOp::prepareForSolve();

    const int hd = hiddenDirection();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        for (auto bc : {m_lobc[0][dir], m_hibc[0][dir]}) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bc == LinOpBCType::Periodic ||
                                             bc == LinOpBCType::Neumann  ||
                                             bc == LinOpBCType::Dirichlet,
                "MLTerrainPoisson supports periodic, Neumann and Dirichlet boundaries only");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dir != hd || bc != LinOpBCType::Dirichlet,
                "MLTerrainPoisson: a hidden (one-cell) direction must be periodic or Neumann");
        }
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_coarse_fine_bc_type == LinOpBCType::Neumann ||
                                     m_coarse_fine_bc_type == LinOpBCType::Dirichlet,
                                     "MLTerrainPoisson: the coarse/fine boundary must be Neumann or Dirichlet");

    // Singular when no face pins the solution (same rule as MLPoisson)
    m_is_singular.clear();
    m_is_singular.resize(m_num_amr_levels, 0);
    bool any_dirichlet = false;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (m_lobc[0][dir] == LinOpBCType::Dirichlet || m_hibc[0][dir] == LinOpBCType::Dirichlet) {
            any_dirichlet = true;
        }
    }
    if (!any_dirichlet && m_domain_covered[0]) {
        m_is_singular[0] = 1;
    }
    if (!m_is_singular[0] && m_needs_coarse_data_for_bc && m_coarse_fine_bc_type == LinOpBCType::Neumann)
    {
        Box bbox = m_grids[0][0].minimalBox();
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            if (m_lobc[0][dir] == LinOpBCType::Dirichlet) { bbox.growLo(dir,1); }
            if (m_hibc[0][dir] == LinOpBCType::Dirichlet) { bbox.growHi(dir,1); }
        }
        if (m_geom[0][0].Domain().contains(bbox)) {
            m_is_singular[0] = 1;
        }
    }

    // Neither the boundary conditions nor the metric terms change between solves on
    // the same object
    if (!m_coeffs_ready) {
        for (int mglev = 0; mglev < NMGLevels(0); ++mglev) {
            define_ghost_fill(mglev);
            compute_tridiagonal(mglev);
        }
        m_coeffs_ready = true;
    }
}

void
MLTerrainPoisson::applyBC (int amrlev, int mglev, MF& in, BCMode /*bc_mode*/, StateMode /*s_mode*/,
                           const MLMGBndry* /*bndry*/, bool skip_fillboundary) const
{
    BL_PROFILE("MLTerrainPoisson::applyBC()");

    // The level boundary data are homogeneous (the projection sets none), so the
    // inhomogeneous and homogeneous fills coincide and bndry is not consulted.

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_coeffs_ready, "MLTerrainPoisson::applyBC before prepareForSolve");

    const Geometry& geom = m_geom[amrlev][mglev];
    const int hd = hiddenDirection();

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        AMREX_ALWAYS_ASSERT(dir == hd || in.nGrowVect()[dir] >= 1);
    }

    if (!skip_fillboundary) {
        // Full fill including the edge and corner ghost cells the stencil reads
        in.FillBoundary(0, 1, geom.periodicity());
    }

    // Every uncovered ghost cell takes its sign times its source, which is a valid or
    // covered cell and so already final: a single gather with no ordering between cells
    IntVect ng(1);
    if (hd >= 0) { ng[hd] = 0; }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(in); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        const Box gbx = amrex::grow(vbx, ng);
        const Dim3 blo = lbound(vbx);
        const Dim3 bhi = ubound(vbx);
        Array4<Real      > const& phi = in.array(mfi);
        Array4<int  const> const& mk  = m_mask[mglev].const_array(mfi);
        Array4<int  const> const& src = m_src[mglev].const_array(mfi);
        Array4<Real const> const& sgn = m_sgn[mglev].const_array(mfi);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            bool inside = (i >= blo.x && i <= bhi.x && j >= blo.y && j <= bhi.y && k >= blo.z && k <= bhi.z);
            if (!inside && mk(i,j,k) == 1) {
                phi(i,j,k) = sgn(i,j,k) * phi(src(i,j,k,0), src(i,j,k,1), src(i,j,k,2));
            }
        });
    }
}

void
MLTerrainPoisson::Fapply (int amrlev, int mglev, MF& out, const MF& in) const
{
    BL_PROFILE("MLTerrainPoisson::Fapply()");

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    const int hd = hiddenDirection();
    const int hidx = hidden_index(mglev);
    const Dim3 dlo = lbound(m_geom[amrlev][mglev].Domain());
    const Dim3 dhi = ubound(m_geom[amrlev][mglev].Domain());
    const Dim3 nper = periodic_lengths(mglev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(out, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real      > const& y    = out.array(mfi);
        Array4<Real const> const& axa  = m_ax[mglev].const_array(mfi);
        Array4<Real const> const& aya  = m_ay[mglev].const_array(mfi);
        Array4<Real const> const& aza  = m_az[mglev].const_array(mfi);
        Array4<Real const> const& dJa  = m_dJ[mglev].const_array(mfi);
        Array4<Real const> const& zpa  = m_zphys[mglev].const_array(mfi);
        Array4<int  const> const& mask = m_mask[mglev].const_array(mfi);
        TerrainStencilAccess x{in.const_array(mfi), hd, hidx};

        if (m_box_has_cf[mglev][mfi] != 0) {
            // zero correction flux through the coarse/fine faces of this box
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                y(i,j,k) = terrain_adotx_cf(i, j, k, x, mask, dlo, dhi, nper, hd, axa, aya, aza, dJa, zpa,
                                            dxinv[0], dxinv[1], dxinv[2]);
            });
        } else {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                terrpoisson_adotx(i, j, k, y, x, axa, aya, aza, dJa, zpa, dxinv[0], dxinv[1], dxinv[2]);
            });
        }
    }
}

void
MLTerrainPoisson::Fsmooth (int amrlev, int mglev, MF& sol, const MF& rhs, int redblack) const
{
    BL_PROFILE("MLTerrainPoisson::Fsmooth()");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_coeffs_ready, "MLTerrainPoisson::Fsmooth before prepareForSolve");

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    const int hd = hiddenDirection();
    const int hidx = hidden_index(mglev);
    const Dim3 dlo = lbound(m_geom[amrlev][mglev].Domain());
    const Dim3 dhi = ubound(m_geom[amrlev][mglev].Domain());
    const Dim3 nper = periodic_lengths(mglev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(sol); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        const int klo = bx.smallEnd(2);
        const int khi = bx.bigEnd(2);
        const Box slab = amrex::makeSlab(bx, 2, klo);

        FArrayBox cpf(bx, 1, The_Async_Arena());
        FArrayBox dpf(bx, 1, The_Async_Arena());

        Array4<Real      > const& s    = sol.array(mfi);
        Array4<Real const> const& b    = rhs.const_array(mfi);
        Array4<Real const> const& tri  = m_tri[mglev].const_array(mfi);
        Array4<Real const> const& axa  = m_ax[mglev].const_array(mfi);
        Array4<Real const> const& aya  = m_ay[mglev].const_array(mfi);
        Array4<Real const> const& aza  = m_az[mglev].const_array(mfi);
        Array4<Real const> const& dJa  = m_dJ[mglev].const_array(mfi);
        Array4<Real const> const& zpa  = m_zphys[mglev].const_array(mfi);
        Array4<Real      > const& cp   = cpf.array();
        Array4<Real      > const& dp   = dpf.array();
        Array4<int  const> const& mask = m_mask[mglev].const_array(mfi);
        TerrainStencilAccess x{sol.const_array(mfi), hd, hidx};
        const bool has_cf = (m_box_has_cf[mglev][mfi] != 0);

        ParallelFor(slab, [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/) noexcept
        {
            if (((i + j + redblack) & 1) == 0) {
                terrain_column_relax(i, j, klo, khi, s, x, b, tri, axa, aya, aza, dJa, zpa,
                                     dxinv[0], dxinv[1], dxinv[2], cp, dp,
                                     has_cf, mask, dlo, dhi, nper, hd);
            }
        });
    }
}

void
MLTerrainPoisson::FFlux (int amrlev, const MFIter& mfi,
                         const Array<FAB*,AMREX_SPACEDIM>& flux,
                         const FAB& sol, Location /*loc*/, int /*face_only*/) const
{
    BL_PROFILE("MLTerrainPoisson::FFlux()");

    const int mglev = 0;
    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    const int hd = hiddenDirection();
    const int hidx = hidden_index(mglev);

    Array4<Real const> const& zpa = m_zphys[mglev].const_array(mfi);
    TerrainStencilAccess x{sol.const_array(), hd, hidx};

    const Box& xbx = flux[0]->box();
    const Box& ybx = flux[1]->box();
    const Box& zbx = flux[2]->box();
    Array4<Real> const& fx = flux[0]->array();
    Array4<Real> const& fy = flux[1]->array();
    Array4<Real> const& fz = flux[2]->array();

    ParallelFor(xbx, ybx, zbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        fx(i,j,k) = terrpoisson_flux_x(i, j, k, x, zpa, dxinv[0]);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        fy(i,j,k) = terrpoisson_flux_y(i, j, k, x, zpa, dxinv[1]);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        fz(i,j,k) = terrpoisson_flux_z(i, j, k, x, zpa, dxinv[0], dxinv[1]);
    });

    if (m_box_has_cf[mglev][mfi] != 0) {
        // no correction flux through a coarse/fine face, as in Fapply
        const Dim3 dlo = lbound(m_geom[amrlev][mglev].Domain());
        const Dim3 dhi = ubound(m_geom[amrlev][mglev].Domain());
        const Dim3 nper = periodic_lengths(mglev);
        Array4<int const> const& mask = m_mask[mglev].const_array(mfi);
        ParallelFor(xbx, ybx, zbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (hd != 0 && (terrain_is_cf_cell(i-1,j,k, mask, dlo, dhi, nper) || terrain_is_cf_cell(i,j,k, mask, dlo, dhi, nper))) { fx(i,j,k) = Real(0.0); }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (hd != 1 && (terrain_is_cf_cell(i,j-1,k, mask, dlo, dhi, nper) || terrain_is_cf_cell(i,j,k, mask, dlo, dhi, nper))) { fy(i,j,k) = Real(0.0); }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (terrain_is_cf_cell(i,j,k-1, mask, dlo, dhi, nper) || terrain_is_cf_cell(i,j,k, mask, dlo, dhi, nper)) { fz(i,j,k) = Real(0.0); }
        });
    }
}

void
MLTerrainPoisson::getFluxes (const Vector<Array<MF*,AMREX_SPACEDIM>>& a_flux,
                             const Vector<MF*>& a_sol, Location a_loc) const
{
    BL_PROFILE("MLTerrainPoisson::getFluxes()");
    for (int alev = 0; alev < m_num_amr_levels; ++alev) {
        compFlux(alev, a_flux[alev], *a_sol[alev], a_loc);
    }
}

void
MLTerrainPoisson::restriction (int amrlev, int cmglev, MF& crse, MF& fine) const
{
    if (!m_weighted_restriction) {
        MLCellLinOp::restriction(amrlev, cmglev, crse, fine);
        return;
    }

    BL_PROFILE("MLTerrainPoisson::restriction()");

    AMREX_ALWAYS_ASSERT(amrlev == 0);
    IntVect ratio = mg_coarsen_ratio_vec[cmglev-1];
    if (hasHiddenDimension()) { ratio[hiddenDirection()] = 1; }

    const MultiFab& wgt = m_dJ[cmglev-1];

    BoxArray fba_c = amrex::coarsen(fine.boxArray(), ratio);
    const bool direct = (fba_c == crse.boxArray()) && (fine.DistributionMap() == crse.DistributionMap());

    MultiFab tmp;
    MultiFab* dst = &crse;
    if (!direct) {
        tmp.define(fba_c, fine.DistributionMap(), 1, 0);
        dst = &tmp;
    }

    const Dim3 r = ratio.dim3();
    for (MFIter mfi(*dst); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        Array4<Real const> const& f = fine.const_array(mfi);
        Array4<Real const> const& w = wgt.const_array(mfi);
        Array4<Real      > const& c = dst->array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int I, int J, int K) noexcept
        {
            Real num = Real(0.0);
            Real den = Real(0.0);
            for (int kk = 0; kk < r.z; ++kk) {
            for (int jj = 0; jj < r.y; ++jj) {
            for (int ii = 0; ii < r.x; ++ii) {
                int i = I*r.x + ii;
                int j = J*r.y + jj;
                int k = K*r.z + kk;
                num += w(i,j,k) * f(i,j,k);
                den += w(i,j,k);
            }}}
            c(I,J,K) = num / den;
        });
    }

    if (!direct) {
        crse.ParallelCopy(tmp, 0, 0, 1, IntVect(0), IntVect(0));
    }
}

Vector<Real>
MLTerrainPoisson::getSolvabilityOffset (int /*amrlev*/, int mglev, MF const& rhs) const
{
    // The operator is (1/dJ) div(...), so its left null vector is dJ: the compatible
    // right-hand side has zero dJ-weighted mean, not zero plain mean.
    Real sums[2];
    sums[0] = MultiFab::Dot(rhs, 0, m_dJ[mglev], 0, 1, 0, true);
    sums[1] = m_dJ[mglev].sum(0, true);
    ParallelAllReduce::Sum(sums, 2, ParallelContext::CommunicatorSub());
    return Vector<Real>{sums[0] / sums[1]};
}

void
MLTerrainPoisson::fixSolvabilityByOffset (int /*amrlev*/, int /*mglev*/, MF& rhs, Vector<Real> const& offset) const
{
    rhs.plus(-offset[0], 0, 1, 0);
}

