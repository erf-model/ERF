/**
 * \file ERF_PlanarBoundary.cpp
 */
#include <ERF_PlanarBoundary.H>

using namespace amrex;

/**
 * Record the surface copies of a planar BoxArray (see ERF_PlanarBoundary.H).
 *
 * @param[in] ba3d 3D BoxArray the planar BoxArray was collapsed from
 * @param[in] ba2d planar BoxArray, one box per box of ba3d and in the same order
 * @param[in] dm   DistributionMapping shared by ba3d and ba2d
 * @param[in] klo  k index of the lowest cell in the domain
 */
void
PlanarBoundary::define (const BoxArray& ba3d,
                        const BoxArray& ba2d,
                        const DistributionMapping& dm,
                        int klo)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ba2d.size() == ba3d.size(),
        "PlanarBoundary::define: the planar BoxArray must hold one box per 3D box");

    m_nplanar = static_cast<int>(ba2d.size());
    m_src_index.clear();
    m_buffers.clear();

    // The planar boxes are taken as they are, whatever k they were collapsed to
    BoxList bl_sfc(IndexType::TheCellType());
    Vector<int> pmap;
    for (int ib = 0; ib < m_nplanar; ++ib) {
        if (ba3d[ib].smallEnd(2) == klo) {
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
    MultiFab& buf = buffer(mf.ixType(), ncomp);
    for (MFIter mfi(buf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const int src = m_src_index[mfi.index()];
        // The surface copy must be the planar box with the same footprint, on this rank
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bx == mf.boxArray()[src] &&
                                         mf.DistributionMap()[src] == ParallelDescriptor::MyProc(),
            "PlanarBoundary::fill: surface box does not match its planar box");
        buf[mfi].copy<RunOn::Device>(mf[src], bx, 0, bx, 0, ncomp);
    }

    // Fill every copy, valid region and ghost cells, from the computed surface copies
    mf.ParallelCopy(buf, 0, 0, ncomp, IntVect(0), mf.nGrowVect(), period);
}

/**
 * Gather buffer for one index type and number of components, allocated on first use.
 *
 * @param[in] ixtype index type of the planar MultiFab
 * @param[in] ncomp  number of components of the planar MultiFab
 */
MultiFab&
PlanarBoundary::buffer (IndexType ixtype, int ncomp)
{
    for (auto& b : m_buffers) {
        if (b.ixtype == ixtype && b.ncomp == ncomp) { return *b.mf; }
    }
    m_buffers.push_back(Buffer{ixtype, ncomp,
                               std::make_unique<MultiFab>(convert(m_ba_sfc, ixtype), m_dm_sfc, ncomp, 0)});
    return *m_buffers.back().mf;
}
