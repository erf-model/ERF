#include <set>

#include <AMReX_BoxList.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MFIter.H>
#include <AMReX_Periodicity.H>

#include "ERF_ColumnBands.H"

using namespace amrex;

Vector<int>
column_bands (const BoxArray& ba)
{
    std::set<int> klo;
    for (int ib = 0; ib < static_cast<int>(ba.size()); ++ib) {
        klo.insert(ba[ib].smallEnd(2));
    }
    return Vector<int>(klo.begin(), klo.end());
}

void
fill_below_band (MultiFab& mf, int icomp, int ncomp, int klo_band,
                 const IntVect& lateral_ng, const Geometry& geom)
{
    AMREX_ALWAYS_ASSERT(lateral_ng[2] == 0);
    AMREX_ALWAYS_ASSERT(mf.nGrowVect().allGE(lateral_ng + IntVect(0,0,1)));

    const BoxArray& ba = mf.boxArray();
    const DistributionMapping& dm = mf.DistributionMap();

    // The slab just below each box of the band, on the rank that owns that box
    BoxList slabs;
    Vector<int> ranks;
    Vector<int> owner;
    for (int ib = 0; ib < static_cast<int>(ba.size()); ++ib) {
        if (ba[ib].smallEnd(2) == klo_band) {
            Box slab = amrex::grow(ba[ib], lateral_ng);
            slab.setRange(2, klo_band-1);
            slabs.push_back(slab);
            ranks.push_back(dm[ib]);
            owner.push_back(ib);
        }
    }
    if (slabs.isEmpty()) { return; }

    MultiFab below(BoxArray(std::move(slabs)), DistributionMapping(std::move(ranks)), ncomp, 0);

    // A cell with no box of this level below it keeps the value it has
    for (MFIter mfi(below); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        below[mfi].copy<RunOn::Device>(mf[owner[mfi.index()]], bx, icomp, bx, 0, ncomp);
    }

    //
    // Only a box below the band reaches a slab: a box of the band or above starts at klo_band or
    //    higher, and neither its lateral ghost cells nor its lateral periodic images extend below
    //    that.  So both copies read columns that are already integrated.  The lateral ghost cells
    //    of the boxes below go first, so that a cell inside a box below then takes that box's own
    //    value.  The periodicity is lateral only: a periodic image in z would bring a box from the
    //    top of the domain below the band.
    //
    const Box& domain = geom.Domain();
    const Periodicity lateral_period(IntVect(geom.isPeriodic(0) ? domain.length(0) : 0,
                                             geom.isPeriodic(1) ? domain.length(1) : 0,
                                             0));
    below.ParallelCopy(mf, icomp, 0, ncomp, lateral_ng, IntVect(0), lateral_period);
    below.ParallelCopy(mf, icomp, 0, ncomp, IntVect(0), IntVect(0), lateral_period);

    for (MFIter mfi(below); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        mf[owner[mfi.index()]].copy<RunOn::Device>(below[mfi], bx, 0, bx, icomp, ncomp);
    }
}
