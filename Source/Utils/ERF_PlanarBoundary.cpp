/**
 * \file ERF_PlanarBoundary.cpp
 */
#include <ERF_PlanarBoundary.H>

#include <algorithm>
#include <utility>

using namespace amrex;

/**
 * Record the surface copies of a planar BoxArray (see ERF_PlanarBoundary.H).
 *
 * @param[in] ba3d          3D BoxArray the planar BoxArray was collapsed from
 * @param[in] ba2d          planar BoxArray, one box per box of ba3d and in the same order
 * @param[in] dm            DistributionMapping shared by ba3d and ba2d
 * @param[in] surface_index index of the surface cell in the normal direction
 * @param[in] is_low        whether the surface is the low side of the domain
 * @param[in] normal_dir    normal direction of the surface
 */
void
PlanarBoundary::define (const BoxArray& ba3d,
                        const BoxArray& ba2d,
                        const DistributionMapping& dm,
                        int surface_index,
                        bool is_low,
                        int normal_dir)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ba2d.size() == ba3d.size(),
        "PlanarBoundary::define: the planar BoxArray must hold one box per 3D box");

    m_nplanar = static_cast<int>(ba2d.size());
    m_src_index.clear();
    m_buffers.clear();
    m_strips.clear();

    // The planar boxes are taken as they are, whatever index they were collapsed to
    BoxList bl_sfc(IndexType::TheCellType());
    Vector<int> pmap;
    for (int ib = 0; ib < m_nplanar; ++ib) {
        const bool touches_surface = is_low
            ? (ba3d[ib].smallEnd(normal_dir) == surface_index)
            : (ba3d[ib].bigEnd(normal_dir) == surface_index);
        if (touches_surface) {
            bl_sfc.push_back(enclosedCells(ba2d[ib]));
            pmap.push_back(dm[ib]);
            m_src_index.push_back(ib);
        }
    }
    m_ba_sfc = BoxArray(std::move(bl_sfc));
    m_dm_sfc = DistributionMapping(std::move(pmap));

    // Exactly one computed copy per column: the surface boxes must not overlap
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_ba_sfc.isDisjoint(),
        "PlanarBoundary::define: the boxes touching the surface overlap in the plane, so a "
        "planar field would have more than one computed copy per column");
}

/**
 * Fill every copy of a planar MultiFab, valid region and ghost cells, from the surface
 * copies (see ERF_PlanarBoundary.H).
 *
 * @param[in,out] mf     planar MultiFab on the planar BoxArray given to define
 * @param[in]     period periodicity of the level
 */
void
PlanarBoundary::fill (MultiFab& mf, const Periodicity& period)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(static_cast<int>(mf.size()) == m_nplanar,
        "PlanarBoundary::fill: the MultiFab is not on the planar BoxArray given to define");

    const int nsfc = static_cast<int>(m_src_index.size());

    // Every planar box is a surface box: no duplicates, FillBoundary is well defined
    if (nsfc == m_nplanar) {
        mf.FillBoundary(period);
        return;
    }

    // No box on this level reaches the surface, so no copy holds computed data to fill
    // from; leave the MultiFab as it is
    if (nsfc == 0) { return; }

    const int ncomp = mf.nComp();
    MultiFab& buf = buffer(mf);
    gather_surface(mf, buf, 0, 0, ncomp);

    // A face-centered buffer's boxes share a face with their neighbors, and the gather above
    // takes each box's face from its own surface copy, so two boxes can hold different values
    // there.  Give every shared face one value before the ParallelCopy -- OverrideSync's
    // precedence is the lowest global box index, so the result does not depend on the
    // decomposition or the rank count -- and every copy of mf then ends up with the same
    // value.  Returns immediately for cell-centered data, which has no shared faces.
    buf.OverrideSync(period);

    // Fill every copy, valid region and ghost cells, from the computed surface copies
    mf.ParallelCopy(buf, 0, 0, ncomp, IntVect(0), mf.nGrowVect(), period);
}

/**
 * Give every copy of a ghost cell outside a non-periodic lateral domain edge the value of
 * the surface copy that holds the column it faces (see ERF_PlanarBoundary.H).
 *
 * @param[in,out] mf   planar MultiFab on the planar BoxArray given to define
 * @param[in]     geom geometry of the level
 */
void
PlanarBoundary::reconcile_outside_domain (MultiFab& mf, const Geometry& geom)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(static_cast<int>(mf.size()) == m_nplanar,
        "PlanarBoundary::reconcile_outside_domain: the MultiFab is not on the planar "
        "BoxArray given to define");

    if (!mf.ixType().cellCentered() || m_src_index.empty()) { return; }

    EdgeStrips& strips = edge_strips(mf, geom);
    if (strips.ba.empty()) { return; }

    // Copy each strip from its owner into a buffer, kept for each number of components ...
    const int ncomp = mf.nComp();
    MultiFab* bufp = nullptr;
    for (auto& b : strips.bufs) {
        if (b->nComp() == ncomp) { bufp = b.get(); break; }
    }
    if (bufp == nullptr) {
        strips.bufs.push_back(std::make_unique<MultiFab>(strips.ba, strips.dm, ncomp, 0));
        bufp = strips.bufs.back().get();
    }
    MultiFab& buf = *bufp;

    // Each strip is on the rank of its owner, so taking it from the owner is a local copy
    for (MFIter mfi(buf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const int owner = strips.owner[mfi.index()];
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mf.DistributionMap()[owner] == ParallelDescriptor::MyProc(),
            "PlanarBoundary::reconcile_outside_domain: an edge strip is not on the rank of its owner");
        const auto src = mf.const_array(owner);
        const auto dst = buf.array(mfi);
        const int si = strips.shift[mfi.index()][0];
        const int sj = strips.shift[mfi.index()][1];
        ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            dst(i,j,k,n) = src(i+si,j+sj,k,n);
        });
    }

    // ... and then into every copy of those cells, the owner's included
    mf.ParallelCopy(buf, 0, 0, ncomp, IntVect(0), mf.nGrowVect(), geom.periodicity());
}

/**
 * The strips beyond the non-periodic domain edges owned by each surface copy of a planar
 * layout, built on first use.
 *
 * A surface copy that touches a non-periodic edge owns the ghost cells beyond it that face
 * its valid region, out to the ghost width of the field, and at a domain corner the cells
 * beyond both edges.  Since the surface copies do not overlap, neither do these strips.
 *
 * The remaining cells beyond a non-periodic edge that lie in the ghost region of a surface
 * copy face a column no grid of the level covers.  Each goes to the surface copy with the
 * lowest global index that holds it, in the image inside the domain in the periodic
 * directions, so that the ParallelCopy from the strips reaches all its images.
 *
 * @param[in] mf   planar MultiFab
 * @param[in] geom geometry of the level
 */
PlanarBoundary::EdgeStrips&
PlanarBoundary::edge_strips (const MultiFab& mf, const Geometry& geom)
{
    const IntVect ng = mf.nGrowVect();
    const BoxArray& pba = mf.boxArray();
    const DistributionMapping& pdm = mf.DistributionMap();
    for (auto& s : m_strips) {
        if (s.ngrow == ng && s.planar_ba == pba && s.planar_dm == pdm) { return s; }
    }

    const Box& dom = geom.Domain();
    BoxList bl;
    Vector<int> pmap;
    Vector<int> owner;
    Vector<IntVect> shift;

    // The cells that face the valid region of a surface copy
    for (int ib : m_src_index) {
        const Box vbx = pba[ib];
        Box grown = vbx;
        for (int d = 0; d < 2; ++d) {
            if (geom.isPeriodic(d) || ng[d] == 0) { continue; }
            if (vbx.smallEnd(d) == dom.smallEnd(d)) { grown.growLo(d, ng[d]); }
            if (vbx.bigEnd(d)   == dom.bigEnd(d))   { grown.growHi(d, ng[d]); }
        }
        if (grown == vbx) { continue; }
        for (const Box& b : boxDiff(grown, vbx)) {
            bl.push_back(b);
            pmap.push_back(pdm[ib]);
            owner.push_back(ib);
            shift.emplace_back(0);
        }
    }

    // The cells facing no grid, to the first surface copy that holds them.  taken lists every
    // cell assigned so far; the surface copies are in increasing global index.
    Vector<Box> taken(bl.begin(), bl.end());
    for (int ib : m_src_index) {
        const Box vbx = pba[ib];
        Box ghost = vbx;
        Box inside = vbx;
        for (int d = 0; d < 2; ++d) {
            ghost.grow(d, ng[d]);
            inside.grow(d, ng[d]);
            if (!geom.isPeriodic(d)) {
                inside.setSmall(d, std::max(inside.smallEnd(d), dom.smallEnd(d)));
                inside.setBig  (d, std::min(inside.bigEnd(d),   dom.bigEnd(d)));
            }
        }
        if (inside == ghost) { continue; }

        for (const Box& beyond : boxDiff(ghost, inside)) {
            // Split along the periodic directions into the parts below, inside and above the
            // domain, and move each into the domain; the owner's cell is the image's + shift
            Vector<std::pair<Box,IntVect>> pieces{{beyond, IntVect(0)}};
            for (int d = 0; d < 2; ++d) {
                if (!geom.isPeriodic(d)) { continue; }
                const int len = dom.length(d);
                Vector<std::pair<Box,IntVect>> split;
                for (const auto& [b, sh] : pieces) {
                    Box lo = b; lo.setBig(d, std::min(b.bigEnd(d), dom.smallEnd(d)-1));
                    Box mid = b; mid.setSmall(d, std::max(b.smallEnd(d), dom.smallEnd(d)));
                                 mid.setBig  (d, std::min(b.bigEnd(d),   dom.bigEnd(d)));
                    Box hi = b; hi.setSmall(d, std::max(b.smallEnd(d), dom.bigEnd(d)+1));
                    IntVect shlo = sh; shlo[d] -= len;
                    IntVect shhi = sh; shhi[d] += len;
                    if (lo.ok())  { split.emplace_back(amrex::shift(lo, d,  len), shlo); }
                    if (mid.ok()) { split.emplace_back(mid, sh); }
                    if (hi.ok())  { split.emplace_back(amrex::shift(hi, d, -len), shhi); }
                }
                pieces = std::move(split);
            }

            for (const auto& [b, sh] : pieces) {
                BoxList rest(b);
                for (const Box& t : taken) {
                    if (!t.intersects(b)) { continue; }
                    BoxList next;
                    for (const Box& r : rest) {
                        if (r.intersects(t)) { next.join(boxDiff(r, t)); }
                        else                 { next.push_back(r); }
                    }
                    rest = std::move(next);
                    if (rest.isEmpty()) { break; }
                }
                for (const Box& r : rest) {
                    bl.push_back(r);
                    pmap.push_back(pdm[ib]);
                    owner.push_back(ib);
                    shift.push_back(sh);
                    taken.push_back(r);
                }
            }
        }
    }

    EdgeStrips s;
    s.planar_ba = pba;
    s.planar_dm = pdm;
    s.ngrow = ng;
    if (!bl.isEmpty()) {
        s.ba = BoxArray(std::move(bl));
        s.dm = DistributionMapping(std::move(pmap));
        s.owner = std::move(owner);
        s.shift = std::move(shift);
    }
    m_strips.push_back(std::move(s));
    return m_strips.back();
}

/**
 * Copy the computed surface copies of a planar MultiFab into a MultiFab without duplicates
 * (see ERF_PlanarBoundary.H).
 *
 * @param[in]  mf    planar MultiFab on the planar BoxArray given to define
 * @param[out] dst   MultiFab on the surface boxes, in the index type of mf
 * @param[in]  scomp first component to read from mf
 * @param[in]  dcomp first component to write in dst
 * @param[in]  ncomp number of components
 */
void
PlanarBoundary::gather_surface (const MultiFab& mf, MultiFab& dst,
                                int scomp, int dcomp, int ncomp) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(static_cast<int>(mf.size()) == m_nplanar,
        "PlanarBoundary::gather_surface: the MultiFab is not on the planar BoxArray given to define");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dst.size() == static_cast<Long>(m_src_index.size()),
        "PlanarBoundary::gather_surface: the destination is not on the surface boxes");

    for (MFIter mfi(dst); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const int src = m_src_index[mfi.index()];
        // The surface copy must be the planar box with the same footprint, on this rank
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bx == mf.boxArray()[src] &&
                                         mf.DistributionMap()[src] == ParallelDescriptor::MyProc(),
            "PlanarBoundary::gather_surface: surface box does not match its planar box");
        dst[mfi].copy<RunOn::Device>(mf[src], bx, scomp, bx, dcomp, ncomp);
    }
}

/**
 * The ground plane of a field on a level's 3D grids, on one planar box per 3D box (see
 * ERF_PlanarBoundary.H).
 *
 * @param[in] src field on the level's 3D grids
 * @param[in] klo index of the lowest cell of the domain
 */
std::unique_ptr<MultiFab>
copy_ground_plane (const MultiFab& src, int klo)
{
    const IndexType ixtype = src.ixType();
    const BoxArray& ba = src.boxArray();
    const int ncomp = src.nComp();

    BoxList bl(ixtype);
    for (int ib = 0; ib < static_cast<int>(ba.size()); ++ib) {
        Box b = enclosedCells(ba[ib]);
        b.setRange(2, klo);
        bl.push_back(convert(b, ixtype));
    }
    auto plane = std::make_unique<MultiFab>(BoxArray(std::move(bl)), src.DistributionMap(), ncomp, 0);

    // Every box first holds the lowest plane of its own 3D box, so that it holds finite
    // values even over columns in which no box reaches the ground
    for (MFIter mfi(*plane); mfi.isValid(); ++mfi) {
        const Box& pbx = mfi.validbox();
        const int dk = enclosedCells(ba[mfi.index()]).smallEnd(2) - klo;
        (*plane)[mfi].copy<RunOn::Device>(src[mfi], amrex::shift(pbx, 2, dk), 0, pbx, 0, ncomp);
    }

    // Then the ground plane, from the boxes that hold it
    plane->ParallelCopy(src, 0, 0, ncomp);

    return plane;
}

/**
 * Finish a column sum computed box by box on a planar MultiFab (see ERF_PlanarBoundary.H).
 *
 * @param[in,out] mf    planar MultiFab, one box per 3D box
 * @param[in]     comp  first component holding a column sum
 * @param[in]     ncomp number of components
 */
void
sum_over_stacked_boxes (MultiFab& mf, int comp, int ncomp)
{
    if (mf.boxArray().isDisjoint()) { return; }

    // The columns covered by the planar boxes, each in exactly one box
    BoxArray ba_cols(mf.boxArray().boxList());
    ba_cols.removeOverlap();
    MultiFab sums(ba_cols, DistributionMapping(ba_cols), ncomp, 0);
    sums.setVal(0.0);

    // Every planar box adds its part into the columns it covers, then each takes the totals back
    sums.ParallelAdd(mf, comp, 0, ncomp);
    mf.ParallelCopy(sums, 0, comp, ncomp);
}

/**
 * Gather buffer for one target layout, index type, and number of components, allocated on
 * first use.  The source MultiFab list avoids rebuilding the layout key on every fill.
 *
 * @param[in] mf planar MultiFab to buffer
 */
MultiFab&
PlanarBoundary::buffer (const MultiFab& mf)
{
    const IndexType ixtype = mf.ixType();
    const int ncomp = mf.nComp();

    // define() clears m_buffers whenever the underlying layout changes, so a
    // source pointer is a stable per-field cache key for the lifetime of this
    // PlanarBoundary definition.
    for (auto& b : m_buffers) {
        if (b.ixtype == ixtype && b.ncomp == ncomp &&
            std::find(b.sources.begin(), b.sources.end(), &mf) != b.sources.end()) {
            return *b.mf;
        }
    }

    // Only an unfamiliar source field needs the more expensive derived-layout
    // construction below.  Fields with the same layout continue to share one
    // gather buffer.
    BoxList bl_sfc(ixtype);
    Vector<int> pmap;
    for (int src : m_src_index) {
        bl_sfc.push_back(mf.boxArray()[src]);
        pmap.push_back(mf.DistributionMap()[src]);
    }
    BoxArray ba(std::move(bl_sfc));
    DistributionMapping dm(std::move(pmap));

    for (auto& b : m_buffers) {
        if (b.ixtype == ixtype && b.ncomp == ncomp && b.ba == ba && b.dm == dm) {
            b.sources.push_back(&mf);
            return *b.mf;
        }
    }
    m_buffers.push_back(Buffer{ixtype, ncomp, ba, dm, {&mf},
                               std::make_unique<MultiFab>(ba, dm, ncomp, 0)});
    return *m_buffers.back().mf;
}
