#ifndef ERF_GTEST_MOST_AVERAGE_PARALLEL_COMMON_H_
#define ERF_GTEST_MOST_AVERAGE_PARALLEL_COMMON_H_

#include "ERF_GTestMOSTAverageCommon.H"

#include <AMReX_FabArrayBase.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <limits>
#include <memory>

namespace erf_most_average_test {

inline Geometry
make_parallel_geometry ()
{
    const Box domain(IntVect(0), IntVect(31, 31, 3));
    const RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                           {AMREX_D_DECL(32.0, 32.0, 4.0)});
    const Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
    return Geometry(domain, &real_box, 0, is_periodic.data());
}

inline BoxArray
make_parallel_box_array (const Box& domain)
{
    BoxArray ba(domain);
    // Keep each box column-complete for TileNoZ(), while providing multiple
    // boxes and multiple lateral tiles for MPI/decomposition coverage.
    ba.maxSize(IntVect(AMREX_D_DECL(16, 16, 1024)));
    return ba;
}

class ScopedMFIterTileSize
{
public:
    explicit ScopedMFIterTileSize (const IntVect& tile_size)
        : m_saved(FabArrayBase::mfiter_tile_size)
    {
        // TileNoZ() reads this value when its MFIter is constructed. Set it
        // identically on every rank and restore it at the end of the test.
        FabArrayBase::mfiter_tile_size = tile_size;
    }

    ~ScopedMFIterTileSize ()
    {
        FabArrayBase::mfiter_tile_size = m_saved;
    }

private:
    IntVect m_saved;
};

struct MOSTAverageFields
{
    Geometry geom;
    Box domain;
    BoxArray ba;
    DistributionMapping dm;
    MultiFab cons;
    MultiFab xvel;
    MultiFab yvel;
    MultiFab zvel;
    std::unique_ptr<MultiFab> theta;
    std::unique_ptr<MultiFab> qv;
    std::unique_ptr<MultiFab> qr;
    std::unique_ptr<MultiFab> z_phys_nd;
    Vector<MultiFab*> vars_old;

    MOSTAverageFields ()
        : MOSTAverageFields(true)
    {}

    // Exercise both one-large-FAB and distributed BoxArray ownership with
    // the same field values and BoxArray index ordering.
    explicit MOSTAverageFields (const bool distributed)
        : MOSTAverageFields(make_parallel_geometry(), distributed)
    {}

    MOSTAverageFields (const Geometry& geometry, const bool distributed)
        : geom(geometry),
          domain(geom.Domain()),
          ba(distributed ? make_parallel_box_array(domain) : BoxArray(domain)),
          dm(ba),
          cons(ba, dm, 1, 1),
          xvel(make_face_box_array(ba, 0), dm, 1, 1),
          yvel(make_face_box_array(ba, 1), dm, 1, 1),
          zvel(make_face_box_array(ba, 2), dm, 1, 1),
          theta(std::make_unique<MultiFab>(ba, dm, 1, 1)),
          qv(std::make_unique<MultiFab>(ba, dm, 1, 1))
    {
        cons.setVal(Real(1.0));
        xvel.setVal(test_u);
        yvel.setVal(test_v);
        zvel.setVal(test_w);
        theta->setVal(test_theta);
        qv->setVal(test_qv);
        vars_old = {&cons, &xvel, &yvel, &zvel};
    }

private:
    static BoxArray
    make_face_box_array (const BoxArray& cell_ba, const int dir)
    {
        BoxArray face_ba(cell_ba);
        face_ba.surroundingNodes(dir);
        return face_ba;
    }
};

struct GlobalRange
{
    Real lo{0};
    Real hi{0};
    Long count{0};
};

inline GlobalRange
selected_average_range (const MultiFab& average,
                        const BoxArray& source_ba,
                        const Box& domain,
                        const Orientation face)
{
    const int dir = face.coordDir();
    const int face_index = face.isLow()
        ? domain.smallEnd(dir)
        : domain.bigEnd(dir) + (source_ba.ixType().nodeCentered(dir) ? 1 : 0);
    Real local_lo = std::numeric_limits<Real>::max();
    Real local_hi = std::numeric_limits<Real>::lowest();
    Long local_count = 0;

    for (MFIter mfi(average, false); mfi.isValid(); ++mfi) {
        const Box& source = source_ba[mfi.index()];
        const bool owns_face = face.isLow()
            ? source.smallEnd(dir) == face_index
            : source.bigEnd(dir) == face_index;
        if (!owns_face) { continue; }

        const Box& valid = mfi.validbox();
        const auto& fab = average[mfi];
        local_count += static_cast<Long>(valid.numPts());
        local_lo = std::min(local_lo, fab.min(valid, 0));
        local_hi = std::max(local_hi, fab.max(valid, 0));
    }

    ParallelDescriptor::ReduceLongSum(local_count);
    ParallelDescriptor::ReduceRealMin(local_lo);
    ParallelDescriptor::ReduceRealMax(local_hi);
    return {local_lo, local_hi, local_count};
}

struct RegionRange
{
    Real selected_lo{0};
    Real selected_hi{0};
    Long selected_count{0};
    Long nonselected_count{0};
    Long nonselected_bad{0};
};

inline RegionRange
region_range (const MultiFab& average,
              const BoxArray& source_ba,
              const Box& domain,
              const Orientation face,
              const Real nonselected_value)
{
    const int dir = face.coordDir();
    const int face_index = face.isLow()
        ? domain.smallEnd(dir)
        : domain.bigEnd(dir) + (source_ba.ixType().nodeCentered(dir) ? 1 : 0);
    Real selected_lo = std::numeric_limits<Real>::max();
    Real selected_hi = std::numeric_limits<Real>::lowest();
    Long selected_count = 0;
    Long nonselected_count = 0;
    Long nonselected_bad = 0;

    for (MFIter mfi(average, false); mfi.isValid(); ++mfi) {
        const Box& source = source_ba[mfi.index()];
        const bool selected = face.isLow()
            ? source.smallEnd(dir) == face_index
            : source.bigEnd(dir) == face_index;
        const Box& valid = mfi.validbox();
        const auto& fab = average[mfi];
        const Long points = static_cast<Long>(valid.numPts());
        if (selected) {
            selected_count += points;
            selected_lo = std::min(selected_lo, fab.min(valid, 0));
            selected_hi = std::max(selected_hi, fab.max(valid, 0));
        } else {
            nonselected_count += points;
            const Real lo = fab.min(valid, 0);
            const Real hi = fab.max(valid, 0);
            if (lo != nonselected_value || hi != nonselected_value) {
                nonselected_bad += points;
            }
        }
    }

    ParallelDescriptor::ReduceLongSum(selected_count);
    ParallelDescriptor::ReduceLongSum(nonselected_count);
    ParallelDescriptor::ReduceLongSum(nonselected_bad);
    ParallelDescriptor::ReduceRealMin(selected_lo);
    ParallelDescriptor::ReduceRealMax(selected_hi);
    return {selected_lo, selected_hi, selected_count,
            nonselected_count, nonselected_bad};
}

} // namespace erf_most_average_test

#endif
