/**
 * \file ERF_PlanarBoundary.cpp
 */
#include <ERF_PlanarBoundary.H>

#include <algorithm>

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

    // A face-centered buffer's boxes share a face with their neighbours, and the gather above
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
