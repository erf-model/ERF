#include <algorithm>
#include <utility>

#include "ERF_MicrophysicsColumns.H"

using namespace amrex;

void
MicrophysicsColumns::define (const BoxArray& ba, const Box& domain)
{
    m_boxes_span_domain = true;
    for (int i = 0; i < ba.size(); ++i) {
        if (ba[i].smallEnd(2) != domain.smallEnd(2) || ba[i].bigEnd(2) != domain.bigEnd(2)) {
            m_boxes_span_domain = false;
            break;
        }
    }

    // The z-flattened footprint of every box. A column is split between boxes exactly
    // when two footprints overlap, since the boxes themselves are disjoint.
    BoxList fl;
    for (int i = 0; i < ba.size(); ++i) {
        Box f(ba[i]); f.setRange(2,0);
        fl.push_back(f);
    }
    const BoxArray fba(std::move(fl));

    m_needs_gather = false;
    for (int i = 0; i < fba.size() && !m_needs_gather; ++i) {
        m_needs_gather = (fba.intersections(fba[i]).size() > 1);
    }

    m_col_ba = BoxArray();
    m_col_dm = DistributionMapping();
    if (!m_needs_gather) { return; }

    // Cut each footprint along the x and y edges of every footprint that overlaps it.
    // Each resulting piece lies either wholly inside or wholly outside every footprint,
    // so all of its columns are covered by the same boxes. A piece is emitted only by
    // the lowest-numbered footprint containing it, so the pieces are disjoint.
    BoxList col;
    for (int i = 0; i < fba.size(); ++i) {
        const Box& F = fba[i];
        const auto overlaps = fba.intersections(F);

        Vector<int> xs, ys;
        for (const auto& [j, bx] : overlaps) {
            xs.push_back(bx.smallEnd(0)); xs.push_back(bx.bigEnd(0)+1);
            ys.push_back(bx.smallEnd(1)); ys.push_back(bx.bigEnd(1)+1);
        }
        std::sort(xs.begin(), xs.end()); xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
        std::sort(ys.begin(), ys.end()); ys.erase(std::unique(ys.begin(), ys.end()), ys.end());

        for (int a = 0; a+1 < xs.size(); ++a) {
            for (int b = 0; b+1 < ys.size(); ++b) {
                const Box p(IntVect(xs[a], ys[b], 0), IntVect(xs[a+1]-1, ys[b+1]-1, 0));

                Vector<std::pair<int,int>> zruns;
                bool owner = true;
                for (const auto& [j, bx] : overlaps) {
                    if (fba[j].contains(p)) {
                        if (j < i) { owner = false; break; }
                        zruns.emplace_back(ba[j].smallEnd(2), ba[j].bigEnd(2));
                    }
                }
                if (!owner) { continue; }

                // Merge the z ranges of the boxes over this piece into contiguous runs
                std::sort(zruns.begin(), zruns.end());
                int lo = zruns[0].first, hi = zruns[0].second;
                for (int n = 1; n <= zruns.size(); ++n) {
                    if (n < zruns.size() && zruns[n].first == hi+1) {
                        hi = zruns[n].second;
                        continue;
                    }
                    col.push_back(Box(IntVect(p.smallEnd(0), p.smallEnd(1), lo),
                                      IntVect(p.bigEnd(0),   p.bigEnd(1),   hi)));
                    if (n < zruns.size()) { lo = zruns[n].first; hi = zruns[n].second; }
                }
            }
        }
    }

    // Merge the pieces, then chop them in x and y (never in z) to no more than the
    // largest box of the original grids, so the layout still spreads over the ranks.
    col.simplify();
    IntVect chunk(1, 1, domain.length(2));
    for (int i = 0; i < ba.size(); ++i) {
        chunk[0] = std::max(chunk[0], ba[i].length(0));
        chunk[1] = std::max(chunk[1], ba[i].length(1));
    }
    col.maxSize(chunk);

    m_col_ba = BoxArray(std::move(col));
    AMREX_ALWAYS_ASSERT(m_col_ba.numPts() == ba.numPts());
    m_col_dm = DistributionMapping(m_col_ba);
}

MicrophysicsColumns::FabPtr
MicrophysicsColumns::gather_one (const MultiFab& src) const
{
    // The kernels stage arrays over the FAB box, so give the copy the same ghost cells and
    // zero them first: ghost cells no source box covers would otherwise be uninitialized.
    // Copy only from the source's valid cells. Its ghost cells may be stale (the
    // accumulators' are never refreshed), and where a source box's ghost region overlaps
    // another source box's valid cells, a copy that read ghost cells could take either.
    auto dst = std::make_shared<MultiFab>(convert(m_col_ba, src.ixType()), m_col_dm,
                                          src.nComp(), src.nGrowVect(),
                                          MFInfo().SetArena(src.arena()));
    dst->setVal(0.0);
    dst->ParallelCopy(src, 0, 0, src.nComp(), IntVect(0), src.nGrowVect());
    return dst;
}

std::unique_ptr<MultiFab>
MicrophysicsColumns::gather_nodal (const MultiFab* src) const
{
    if (!src) { return nullptr; }
    auto dst = std::make_unique<MultiFab>(convert(m_col_ba, src->ixType()), m_col_dm,
                                          src->nComp(), src->nGrowVect(),
                                          MFInfo().SetArena(src->arena()));
    dst->setVal(0.0);
    dst->ParallelCopy(*src, 0, 0, src->nComp(), IntVect(0), src->nGrowVect());
    return dst;
}

std::unique_ptr<iMultiFab>
MicrophysicsColumns::gather_2d (const iMultiFab* src) const
{
    if (!src) { return nullptr; }
    BoxList bl = m_col_ba.boxList();
    for (auto& b : bl) { b.setRange(2,0); }
    auto dst = std::make_unique<iMultiFab>(BoxArray(std::move(bl)), m_col_dm,
                                           src->nComp(), src->nGrowVect());
    dst->setVal(0);
    dst->ParallelCopy(*src, 0, 0, src->nComp(), IntVect(0), src->nGrowVect());
    return dst;
}
