#include <limits>
#include <utility>

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>
#include <AMReX_iMultiFab.H>

#include <gtest/gtest.h>

#include <ERF_SEBSkinInit.H>

// The first value of the force-restore balance's skin temperature from the skin
// temperature of the initial file (ERF_SEBSkinInit.H, erf.radiation.seb_init_skin_from_file).
// An 8 x 4 level in two boxes, i = 0..3 and i = 4..7. The balance's skin starts at 300 K in
// every column and ghost column. The file's skin is 280 + 0.5 i + 0.25 j [K], and -1 K in
// its ghost columns. Water is i >= 2 with j < 2 (in both boxes); the column (0, 3) carries
// the mask of a building's column (2), which is land.

namespace {

using namespace amrex;

constexpr Real t_default = Real(300.0);
constexpr Real t_lo = Real(200.0);
constexpr Real t_hi = Real(340.0);

Real file_value (int i, int j) { return Real(280.0) + Real(0.5) * i + Real(0.25) * j; }
bool is_water (int i, int j) { return i >= 2 && j < 2; }

// The value of one column, valid or ghost, summed over the boxes that hold it.
std::pair<Real, int> value_and_count (const MultiFab& mf, const IntVect& point)
{
    ReduceOps<ReduceOpSum, ReduceOpSum> reduce_op;
    ReduceData<Real, int> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    const Box point_box(point, point);
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box overlap = point_box & mfi.fabbox();
        if (overlap.isEmpty()) { continue; }
        const auto array = mf.const_array(mfi);
        reduce_op.eval(overlap, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple { return {array(i, j, k), 1}; });
    }
    Gpu::streamSynchronize();
    Real value = amrex::get<0>(reduce_data.value());
    int count = amrex::get<1>(reduce_data.value());
    ParallelDescriptor::ReduceRealSum(value);
    ParallelDescriptor::ReduceIntSum(count);
    return {value, count};
}

// A valid column lies in one box's valid region; it may also be a ghost column of the
// other box. Both copies must agree for the caller's check, so the mean is returned.
Real value_at (const MultiFab& mf, int i, int j)
{
    const auto vc = value_and_count(mf, IntVect(i, j, 0));
    return vc.first / static_cast<Real>(vc.second);
}

void set_column (MultiFab& mf, int ic, int jc, Real value)
{
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const auto a = mf.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (i == ic && j == jc) { a(i, j, k) = value; }
        });
    }
    Gpu::streamSynchronize();
}

struct SkinCase
{
    Box domain{IntVect(0, 0, 0), IntVect(7, 3, 0)};
    BoxArray ba;
    DistributionMapping dm;
    MultiFab skin, file;
    iMultiFab mask;

    SkinCase ()
    {
        ba.define(domain);
        ba.maxSize(IntVect(4, 4, 1));
        dm.define(ba);
        const IntVect ng(1, 1, 0);
        skin.define(ba, dm, 1, ng);
        file.define(ba, dm, 1, ng);
        mask.define(ba, dm, 1, ng);
        skin.setVal(t_default);
        file.setVal(Real(-1.0));
        // The mask's ghost columns read as water: only the valid columns may be asked.
        mask.setVal(0);
        fill(file, mask);
    }

    static void fill (MultiFab& file_mf, iMultiFab& mask_mf)
    {
        for (MFIter mfi(file_mf); mfi.isValid(); ++mfi) {
            const auto f = file_mf.array(mfi);
            const auto m = mask_mf.array(mfi);
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                f(i, j, k) = Real(280.0) + Real(0.5) * i + Real(0.25) * j;
                m(i, j, k) = (i >= 2 && j < 2) ? 0 : ((i == 0 && j == 3) ? 2 : 1);
            });
        }
        Gpu::streamSynchronize();
    }
};

SEBSkinInitReport over_ranks (SEBSkinInitReport r)
{
    ParallelDescriptor::ReduceLongSum(r.n_set);
    ParallelDescriptor::ReduceLongSum(r.n_refused);
    ParallelDescriptor::ReduceRealMin(r.t_min);
    ParallelDescriptor::ReduceRealMax(r.t_max);
    return r;
}

} // namespace

TEST(SEBSkinInit, UsableIsFiniteAndInsideTheBalanceRange)
{
    EXPECT_TRUE(seb_skin_init_usable(Real(288.0), t_lo, t_hi));
    EXPECT_TRUE(seb_skin_init_usable(t_lo, t_lo, t_hi));
    EXPECT_TRUE(seb_skin_init_usable(t_hi, t_lo, t_hi));
    EXPECT_FALSE(seb_skin_init_usable(Real(199.99), t_lo, t_hi));
    EXPECT_FALSE(seb_skin_init_usable(Real(340.01), t_lo, t_hi));
    EXPECT_FALSE(seb_skin_init_usable(Real(0.0), t_lo, t_hi));
    EXPECT_FALSE(seb_skin_init_usable(Real(-999.0), t_lo, t_hi));
    EXPECT_FALSE(seb_skin_init_usable(std::numeric_limits<Real>::quiet_NaN(), t_lo, t_hi));
    EXPECT_FALSE(seb_skin_init_usable(std::numeric_limits<Real>::infinity(), t_lo, t_hi));
}

// Land columns take the file's value; water columns and every ghost column keep 300 K.
TEST(SEBSkinInit, LandTakesTheFileWaterKeepsTheDefault)
{
    SkinCase c;
    const SEBSkinInitReport r = over_ranks(seb_init_skin_from_file(c.skin, c.file, &c.mask, t_lo, t_hi));

    // 32 columns, 12 of them water (i = 2..7, j = 0..1).
    EXPECT_EQ(r.n_set, 20);
    EXPECT_EQ(r.n_refused, 0);
    EXPECT_EQ(r.t_min, file_value(0, 0));   // 280
    EXPECT_EQ(r.t_max, file_value(7, 3));   // 284.25
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 8; ++i) {
            const auto vc = value_and_count(c.skin, IntVect(i, j, 0));
            // The columns next to the box seam, i = 3 and i = 4, are also ghost columns of
            // the other box, and those were not written: the caller fills them.
            const bool seam = (i == 3 || i == 4);
            ASSERT_EQ(vc.second, seam ? 2 : 1) << "column " << i << "," << j;
            const Real own = is_water(i, j) ? t_default : file_value(i, j);
            EXPECT_EQ(vc.first, seam ? own + t_default : own) << "column " << i << "," << j;
        }
    }
    // The building's column is land.
    EXPECT_EQ(value_at(c.skin, 0, 3), file_value(0, 3));
    // Ghost columns beyond the domain keep the default.
    EXPECT_EQ(value_at(c.skin, -1, 2), t_default);
    EXPECT_EQ(value_at(c.skin, 8, 3), t_default);
    EXPECT_EQ(value_at(c.skin, 5, 4), t_default);
}

// With no mask every column is land.
TEST(SEBSkinInit, NoMaskMeansAllLand)
{
    SkinCase c;
    const SEBSkinInitReport r = over_ranks(seb_init_skin_from_file(c.skin, c.file, nullptr, t_lo, t_hi));
    EXPECT_EQ(r.n_set, 32);
    EXPECT_EQ(r.n_refused, 0);
    EXPECT_EQ(value_at(c.skin, 6, 0), file_value(6, 0));
    EXPECT_EQ(r.t_max, file_value(7, 3));
}

// A land value that is not usable is counted and the column is left as it was. A water
// column's value is not read. The bounds themselves are usable.
TEST(SEBSkinInit, UnusableLandValuesAreCountedAndNotCopied)
{
    SkinCase c;
    set_column(c.file, 0, 0, std::numeric_limits<Real>::quiet_NaN());
    set_column(c.file, 5, 2, Real(150.0));     // below the range, in the second box
    set_column(c.file, 1, 3, Real(340.5));     // above the range
    set_column(c.file, 1, 1, t_lo);            // the lower bound
    set_column(c.file, 6, 3, t_hi);            // the upper bound
    set_column(c.file, 6, 1, std::numeric_limits<Real>::quiet_NaN());   // water
    set_column(c.file, 3, 0, Real(-999.0));                             // water

    const SEBSkinInitReport r = over_ranks(seb_init_skin_from_file(c.skin, c.file, &c.mask, t_lo, t_hi));
    EXPECT_EQ(r.n_refused, 3);
    EXPECT_EQ(r.n_set, 17);
    EXPECT_EQ(r.t_min, t_lo);
    EXPECT_EQ(r.t_max, t_hi);
    EXPECT_EQ(value_at(c.skin, 0, 0), t_default);
    EXPECT_EQ(value_at(c.skin, 5, 2), t_default);
    EXPECT_EQ(value_at(c.skin, 1, 3), t_default);
    EXPECT_EQ(value_at(c.skin, 1, 1), t_lo);
    EXPECT_EQ(value_at(c.skin, 6, 3), t_hi);
    EXPECT_EQ(value_at(c.skin, 6, 1), t_default);
    EXPECT_EQ(value_at(c.skin, 2, 2), file_value(2, 2));
}

// A level of water only: nothing is set and nothing is refused, whatever the file holds.
TEST(SEBSkinInit, AllWaterSetsNothing)
{
    SkinCase c;
    c.mask.setVal(0);
    c.file.setVal(Real(-999.0));
    const SEBSkinInitReport r = over_ranks(seb_init_skin_from_file(c.skin, c.file, &c.mask, t_lo, t_hi));
    EXPECT_EQ(r.n_set, 0);
    EXPECT_EQ(r.n_refused, 0);
    EXPECT_EQ(value_at(c.skin, 0, 0), t_default);
    EXPECT_EQ(value_at(c.skin, 7, 3), t_default);
}
