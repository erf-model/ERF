#include <cmath>
#include <memory>
#include <string>
#include <utility>

#include <AMReX_BoxArray.H>
#include <AMReX_BoxList.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_Math.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

#include <gtest/gtest.h>

#include "ERF_ColumnBands.H"
#include "../../../Exec/ERF_Prob.H"

namespace {

// Largest |a| over the given region of every box of mf (valid cells grown by ng)
amrex::Real max_abs (const amrex::MultiFab& mf, const amrex::IntVect& ng)
{
    amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const auto a = mf.const_array(mfi);
        reduce_op.eval(amrex::grow(mfi.validbox(), ng), reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Real> {
                return { std::abs(a(i,j,k)) };
            });
    }
    amrex::Real m = amrex::get<0>(reduce_data.value());
    amrex::ParallelDescriptor::ReduceRealMax(m);
    return m;
}

// ---------------------------------------------------------------------------------------------
// fill_below_band on its own
// ---------------------------------------------------------------------------------------------

//
// NOTE: the loops that launch device kernels live here rather than in the TEST bodies below.
//       nvcc rejects an extended __device__ lambda in a function with private or protected
//       access within its class, and gtest makes each TEST body a private TestBody() member.
//

// Valid cells of the left and right boxes hold 1 and 5, their ghost cells 2; the upper box
// holds 3 everywhere
void set_precedence_values (amrex::MultiFab& mf, const amrex::Box& left, const amrex::Box& upper)
{
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const auto a = mf.array(mfi);
        const amrex::Box vbx = mfi.validbox();
        const bool is_upper = (vbx == upper);
        const amrex::Real valid = is_upper ? 3.0 : ((vbx == left) ? 1.0 : 5.0);
        const amrex::Real ghost = is_upper ? 3.0 : 2.0;
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i,j,k) = vbx.contains(i,j,k) ? valid : ghost;
        });
    }
}

// Subtract from mf what fill_below_band should have left, and store the difference in err.
// Below the upper box: a cell inside a box below takes that box's value even where it is also
// a ghost cell of its neighbour; a ghost cell of a box below takes the ghost value; a cell with
// no box below keeps the upper box's value; the upper box's own cells and the boxes below are
// untouched.
void subtract_precedence_expected (const amrex::MultiFab& mf, amrex::MultiFab& err,
                                   const amrex::Box& left, const amrex::Box& upper)
{
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const auto a = mf.const_array(mfi);
        const auto e = err.array(mfi);
        const amrex::Box vbx = mfi.validbox();
        const bool is_upper = (vbx == upper);
        const amrex::Real valid = is_upper ? 3.0 : ((vbx == left) ? 1.0 : 5.0);
        const amrex::Real ghost = is_upper ? 3.0 : 2.0;
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            amrex::Real expected = vbx.contains(i,j,k) ? valid : ghost;
            if (is_upper && k == 3) {
                const bool yin = (j >= 0 && j <= 7);
                if      (i >= 0 && i <= 3)   { expected = yin ? 1.0 : 2.0; }
                else if (i >= 4 && i <= 5)   { expected = yin ? 5.0 : 2.0; }
                else if (i == -1 || i == 6)  { expected = 2.0; }
                else                         { expected = 3.0; }
            }
            e(i,j,k) = a(i,j,k) - expected;
        });
    }
}

// Two boxes side by side at the bottom, one wider box stacked on them; not periodic
TEST(ColumnBandsParallel, FillBelowBandPrecedence)
{
    const amrex::Box domain(amrex::IntVect(0,0,0), amrex::IntVect(7,7,7));
    const amrex::RealBox rb(0.0, 0.0, 0.0, 8.0, 8.0, 8.0);
    const amrex::Geometry geom(domain, rb, amrex::CoordSys::cartesian, {0, 0, 0});

    const amrex::Box left (amrex::IntVect(0,0,0), amrex::IntVect(3,7,3));
    const amrex::Box right(amrex::IntVect(4,0,0), amrex::IntVect(5,7,3));
    const amrex::Box upper(amrex::IntVect(0,0,4), amrex::IntVect(7,7,7));

    amrex::BoxList bl;
    bl.push_back(left); bl.push_back(right); bl.push_back(upper);
    const amrex::BoxArray ba(std::move(bl));
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab mf(ba, dm, 1, 1);

    EXPECT_EQ(column_bands(ba), amrex::Vector<int>({0, 4}));

    set_precedence_values(mf, left, upper);

    fill_below_band(mf, 0, 1, 4, amrex::IntVect(1,1,0), geom);

    amrex::MultiFab err(ba, dm, 1, 1);
    err.setVal(0.0);
    subtract_precedence_expected(mf, err, left, upper);
    EXPECT_EQ(max_abs(err, amrex::IntVect(1)), 0.0);
}

// ---------------------------------------------------------------------------------------------
// The hydrostatic density over terrain, erf_init_dens_hse_dry
// ---------------------------------------------------------------------------------------------

// An 8 x 8 x 24 cell domain, 800 m x 800 m x 1200 m, periodic laterally
constexpr int nx = 8;
constexpr int ny = 8;
constexpr int nz = 24;

amrex::Geometry make_geom ()
{
    const amrex::Box domain(amrex::IntVect(0,0,0), amrex::IntVect(nx-1,ny-1,nz-1));
    const amrex::RealBox rb(0.0, 0.0, 0.0, 800.0, 800.0, 1200.0);
    return amrex::Geometry(domain, rb, amrex::CoordSys::cartesian, {1, 1, 0});
}

// i taken into [0, n), so that a ghost cell holds exactly the value of its periodic image, as it
// would after FillBoundary
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int periodic_index (int i, int n) noexcept
{
    return ((i % n) + n) % n;
}

// Terrain-following node heights over a periodic 60 m bump
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real z_node (int i, int j, int k) noexcept
{
    const amrex::Real twopi = amrex::Real(2.0) * amrex::Math::pi<amrex::Real>();
    const amrex::Real x = amrex::Real(periodic_index(i, nx)) / amrex::Real(nx);
    const amrex::Real y = amrex::Real(periodic_index(j, ny)) / amrex::Real(ny);
    const amrex::Real h = amrex::Real(15.0) * (amrex::Real(1.0) - std::cos(twopi*x))
                                            * (amrex::Real(1.0) - std::cos(twopi*y));
    return h + (amrex::Real(1200.0) - h) * amrex::Real(k) / amrex::Real(nz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real z_cell (int i, int j, int k) noexcept
{
    return amrex::Real(0.125) *
        (z_node(i,j,k  ) + z_node(i+1,j,k  ) + z_node(i,j+1,k  ) + z_node(i+1,j+1,k  ) +
         z_node(i,j,k+1) + z_node(i+1,j,k+1) + z_node(i,j+1,k+1) + z_node(i+1,j+1,k+1));
}

struct Column
{
    amrex::MultiFab rho;
    std::unique_ptr<amrex::MultiFab> z_cc;
};

// Cell-centre heights, and a density that stands in for the one interpolated from a coarser
// level: out of balance and varying in every direction, so where an integration starts matters
Column make_column (const amrex::BoxArray& ba)
{
    const amrex::DistributionMapping dm(ba);
    Column c;
    c.rho.define(ba, dm, 1, 1);
    c.z_cc = std::make_unique<amrex::MultiFab>(ba, dm, 1, 1);
    for (amrex::MFIter mfi(c.rho); mfi.isValid(); ++mfi) {
        const auto rho = c.rho.array(mfi);
        const auto zcc = c.z_cc->array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real twopi = amrex::Real(2.0) * amrex::Math::pi<amrex::Real>();
            const amrex::Real x = amrex::Real(periodic_index(i, nx)) / amrex::Real(nx);
            const amrex::Real y = amrex::Real(periodic_index(j, ny)) / amrex::Real(ny);
            zcc(i,j,k) = z_cell(i,j,k);
            rho(i,j,k) = amrex::Real(1.1) * std::exp(-z_cell(i,j,k) / amrex::Real(8000.0))
                       * (amrex::Real(1.0) + amrex::Real(0.03) * std::sin(twopi*x + amrex::Real(0.5)*amrex::Real(k))
                                                               * std::cos(twopi*y));
        });
    }
    return c;
}

void init_dens_hse (Column& c)
{
    const amrex::Geometry geom = make_geom();
    Problem prob(geom.ProbLo(), geom.ProbHi());
    std::unique_ptr<amrex::MultiFab> z_nd; // unused on this path
    const amrex::Vector<amrex::Real> stretched_dz;
    prob.erf_init_dens_hse_dry(c.rho, z_nd, c.z_cc, geom, stretched_dz, false, false);
}

// Largest |other - ref| over the valid and lateral ghost cells of other; ref is one box whose
// ghost cells cover all of them.  Ghost cells in z are not compared: a box's top ghost cell is
// not written by the integration.
amrex::Real max_difference (const amrex::MultiFab& ref, const amrex::MultiFab& other)
{
    const amrex::IntVect ng(1,1,0);
    amrex::MultiFab diff(other.boxArray(), other.DistributionMap(), 1, ng);
    diff.ParallelCopy(ref, 0, 0, 1, ref.nGrowVect(), ng);
    amrex::MultiFab::Subtract(diff, other, 0, 0, 1, ng);
    return max_abs(diff, ng);
}

amrex::BoxArray chopped (const amrex::Box& region, const amrex::IntVect& max_size)
{
    amrex::BoxArray ba(region);
    ba.maxSize(max_size);
    return ba;
}

// Two boxes below, three above with different footprints and heights
amrex::BoxArray staggered_footprints (const amrex::Box& r)
{
    const int xm = (r.smallEnd(0) + r.bigEnd(0) + 1) / 2;
    const int ym = (r.smallEnd(1) + r.bigEnd(1) + 1) / 2;
    const int z1 = r.smallEnd(2) + (r.length(2) * 5) / 12;
    const int z2 = r.smallEnd(2) + (r.length(2) * 9) / 12;
    const auto lo = r.smallEnd();
    const auto hi = r.bigEnd();
    amrex::BoxList bl;
    bl.push_back(amrex::Box(amrex::IntVect(lo[0],lo[1],lo[2]), amrex::IntVect(xm-1, hi[1], z1-1)));
    bl.push_back(amrex::Box(amrex::IntVect(xm,   lo[1],lo[2]), amrex::IntVect(hi[0],hi[1], z1-1)));
    bl.push_back(amrex::Box(amrex::IntVect(lo[0],lo[1],z1   ), amrex::IntVect(hi[0],ym-1,  hi[2])));
    bl.push_back(amrex::Box(amrex::IntVect(lo[0],ym,   z1   ), amrex::IntVect(hi[0],hi[1], z2-1)));
    bl.push_back(amrex::Box(amrex::IntVect(lo[0],ym,   z2   ), amrex::IntVect(hi[0],hi[1], hi[2])));
    return amrex::BoxArray(std::move(bl));
}

void expect_layouts_match_one_box (const amrex::Box& region)
{
    Column ref = make_column(amrex::BoxArray(region));
    init_dens_hse(ref);

    // The integration must have changed the density, or the comparison proves nothing
    const Column untouched = make_column(amrex::BoxArray(region));
    ASSERT_GT(max_difference(untouched.rho, ref.rho), amrex::Real(1.0e-3));

    const amrex::IntVect half(region.length(0)/2, region.length(1)/2, nz);
    const std::pair<const char*, amrex::BoxArray> layouts[] = {
        {"columns",                    chopped(region, half)},
        {"3 boxes per column",         chopped(region, amrex::IntVect(nx, ny, 8))},
        {"columns of 3 boxes",         chopped(region, amrex::IntVect(half[0], half[1], 8))},
        {"columns of 5 uneven boxes",  chopped(region, amrex::IntVect(half[0], half[1], 5))},
        {"staggered footprints",       staggered_footprints(region)},
    };

    // The same arithmetic runs on the same numbers whatever the boxes, so the densities are equal
    for (const auto& layout : layouts) {
        SCOPED_TRACE(layout.first);
        Column split = make_column(layout.second);
        init_dens_hse(split);
        EXPECT_EQ(max_difference(ref.rho, split.rho), amrex::Real(0.0));
    }
}

} // namespace

// A column split into boxes stacked in z continues its integration across the split, so the
// density is the one a single box gives, in the valid cells and the lateral ghost cells; each box
// used to start afresh from the ghost cell below it
TEST(ColumnBandsParallel, InitDensHSESplitInZMatchesOneBox)
{
    expect_layouts_match_one_box(make_geom().Domain());
}

// A refined patch, whose lateral ghost cells lie outside the fine grids
TEST(ColumnBandsParallel, InitDensHSEPatchMatchesOneBox)
{
    expect_layouts_match_one_box(amrex::Box(amrex::IntVect(2,2,0), amrex::IntVect(5,5,nz-1)));
}

// A refined patch aloft: the lowest boxes start from the ghost cell below them, and the boxes
// stacked on them continue from there
TEST(ColumnBandsParallel, InitDensHSEPatchAloftMatchesOneBox)
{
    expect_layouts_match_one_box(amrex::Box(amrex::IntVect(2,2,6), amrex::IntVect(5,5,nz-1)));
}
