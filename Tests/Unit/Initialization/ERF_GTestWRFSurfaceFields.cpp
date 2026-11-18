#include <cmath>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Periodicity.H>
#include <AMReX_Utility.H>

#include <gtest/gtest.h>

#include <ERF_WRFSurfaceFields.H>

// The holder of the wrfinput land and soil fields that no land model takes
// (WRFSurfaceFields). A level of 6 x 5 columns in boxes of at most 4 x 3, one ghost
// column, not periodic. The source carries a different value in every column and soil
// layer, so a swapped index, a reversed layer order or a missed clamp changes the answer.

namespace {

using namespace amrex;
using Field = WRFSurfaceFields::Field;
using Kind  = WRFSurfaceFields::Kind;

constexpr int nx = 6;
constexpr int ny = 5;

struct ScopedThrowOnAbort
{
    ScopedThrowOnAbort () : m_saved(amrex::system::throw_exception) { amrex::system::throw_exception = true; }
    ~ScopedThrowOnAbort () { amrex::system::throw_exception = m_saved; }
    bool m_saved;
};

struct Level
{
    Box domain{IntVect(0, 0, 0), IntVect(nx - 1, ny - 1, 0)};
    BoxArray ba2d;
    DistributionMapping dm;
    IntVect ng{1, 1, 0};
    explicit Level (const IntVect& max_size = IntVect(4, 3, 1))
    {
        ba2d = BoxArray(domain);
        ba2d.maxSize(max_size);
        dm = DistributionMapping(ba2d);
    }
};

// Multiples of 1/64: exact in single and double precision, so the comparisons below are exact
Real albedo_at (int i, int j) { return Real(0.015625) * Real(i) + Real(0.125) * Real(j); }
Real soil_at (int i, int j, int n) { return Real(280.0) + Real(n) + Real(0.015625) * Real(i) + Real(0.125) * Real(j); }

// Columns ilo..ihi x jlo..jhi with nlay values in z from k0, as a rank reads them.
template <typename F>
FArrayBox make_source (int ilo, int ihi, int jlo, int jhi, int nlay, int k0, F&& value, bool nodal_z = false)
{
    Box bx(IntVect(ilo, jlo, k0), IntVect(ihi, jhi, k0 + nlay - 1));
    if (nodal_z) { bx.setType(IndexType(IntVect(0, 0, 1))); }
    FArrayBox src(bx, 1, The_Pinned_Arena());
    const auto a = src.array();
    for (int k = 0; k < nlay; ++k) {
        for (int j = jlo; j <= jhi; ++j) {
            for (int i = ilo; i <= ihi; ++i) { a(i, j, k0 + k) = value(i, j, k); }
        }
    }
    return src;
}

// The whole level from one file
template <typename F>
FArrayBox make_source (int nlay, F&& value)
{
    return make_source(0, nx - 1, 0, ny - 1, nlay, 0, std::forward<F>(value));
}

Box columns (int ilo, int ihi, int jlo, int jhi)
{
    return Box(IntVect(ilo, jlo, 0), IntVect(ihi, jhi, 0));
}

Box all_columns () { return columns(0, nx - 1, 0, ny - 1); }

// A soil layer column of the file: one (i,j), nodal in z as the reader builds it
FArrayBox make_column (std::initializer_list<Real> vals)
{
    Box bx(IntVect(0, 0, 0), IntVect(0, 0, static_cast<int>(vals.size()) - 1));
    bx.setType(IndexType(IntVect(0, 0, 1)));
    FArrayBox col(bx, 1, The_Pinned_Arena());
    const auto a = col.array();
    int k = 0;
    for (const Real v : vals) { a(0, 0, k++) = v; }
    return col;
}

// Visit every cell of the field, ghost cells included, on the host.
template <typename F>
void for_each_cell (const MultiFab& mf, F&& visit)
{
    for (MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const Box gbx = mfi.fabbox();
        FArrayBox host(gbx, mf.nComp(), The_Pinned_Arena());
        Gpu::dtoh_memcpy(host.dataPtr(), mf[mfi].dataPtr(), host.nBytes());
        const auto a = host.const_array();
        const Dim3 lo = lbound(gbx);
        const Dim3 hi = ubound(gbx);
        for (int n = 0; n < mf.nComp(); ++n) {
            for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) { visit(i, j, n, a(i, j, lo.z, n)); }
            }
        }
    }
}

int clamp (int v, int lo, int hi) { return std::min(std::max(v, lo), hi); }

// Every cell, ghost cells included, holds the value of the nearest column of the level
template <typename F>
void expect_field (const MultiFab& mf, F&& value)
{
    int ncell = 0;
    for_each_cell(mf, [&] (int i, int j, int n, Real v) {
        EXPECT_EQ(v, value(clamp(i, 0, nx - 1), clamp(j, 0, ny - 1), n)) << "i=" << i << " j=" << j << " n=" << n;
        ++ncell;
    });
    EXPECT_GT(ncell, 0);
}

Real fine_albedo_at (int i, int j) { return Real(0.25) + Real(0.0078125) * Real(i) + Real(0.015625) * Real(j); }
Real fine_soil_at (int i, int j, int n) { return Real(300.0) + Real(n) + Real(0.015625) * Real(i) + Real(0.125) * Real(j); }

// Level 0 of 6 x 5 columns with ALBEDO, LAI, TSLB and ZS; level 1, refined by 2, on the
// fine columns 2..5 x 2..5 with ALBEDO, TSLB and ZS of its own (other values than level 0)
struct TwoLevels
{
    static constexpr int nlay = 3;
    Level L;
    IntVect ng{1, 1, 0};
    Box fine_box{IntVect(2, 2, 0), IntVect(5, 5, 0)};
    BoxArray fine_ba;
    DistributionMapping fine_dm;
    WRFSurfaceFields f;
    TwoLevels ()
    {
        fine_ba = BoxArray(fine_box);
        fine_ba.maxSize(IntVect(2, 4, 1));
        fine_dm = DistributionMapping(fine_ba);
        f.resize(2);
        const FArrayBox alb  = make_source(1, [] (int i, int j, int) { return albedo_at(i, j); });
        const FArrayBox lai  = make_source(1, [] (int i, int j, int) { return Real(2.0) + albedo_at(i, j); });
        const FArrayBox soil = make_source(nlay, [] (int i, int j, int n) { return soil_at(i, j, n); });
        f.fill(0, Field::ALBEDO, L.ba2d, L.dm, ng, 1,    alb,  all_columns());
        f.fill(0, Field::LAI,    L.ba2d, L.dm, ng, 1,    lai,  all_columns());
        f.fill(0, Field::TSLB,   L.ba2d, L.dm, ng, nlay, soil, all_columns());
        f.set_soil_layers(0, Field::ZS, make_column({Real(0.05), Real(0.25), Real(0.7)}));
        f.finish_level(0, Periodicity::NonPeriodic());
        const FArrayBox falb  = make_source(2, 5, 2, 5, 1, 0, [] (int i, int j, int) { return fine_albedo_at(i, j); });
        const FArrayBox fsoil = make_source(2, 5, 2, 5, nlay, 0, [] (int i, int j, int n) { return fine_soil_at(i, j, n); });
        f.fill(1, Field::ALBEDO, fine_ba, fine_dm, ng, 1,    falb,  fine_box);
        f.fill(1, Field::TSLB,   fine_ba, fine_dm, ng, nlay, fsoil, fine_box);
        f.set_soil_layers(1, Field::ZS, make_column({Real(0.125), Real(0.25), Real(0.75)}));
        f.finish_level(1, Periodicity::NonPeriodic());
    }
};

// One file that covers the level, then the close of the level
void fill_whole (WRFSurfaceFields& f, const Level& L, Field id, const FArrayBox& src, int nlay = 1)
{
    f.fill(0, id, L.ba2d, L.dm, L.ng, nlay, src, all_columns());
    f.finish_level(0, Periodicity::NonPeriodic());
}

} // namespace

TEST(WRFSurfaceFields, CatalogNamesKindsAndBounds)
{
    const auto& cat = WRFSurfaceFields::catalog();
    constexpr int nf = 13;
    const char* names[nf] = {"LAI", "VEGFRA", "TMN", "ZNT", "ALBEDO", "ALBBCK", "SEAICE", "SNOW", "SNOWH",
                             "TSLB", "SMOIS", "ZS", "DZS"};
    const char* units[nf] = {"m2/m2", "%", "K", "m", "1", "1", "1", "kg/m2", "m", "K", "m3/m3", "m", "m"};
    const Real  hi[nf]    = {Real(100.0), Real(100.0), Real(1000.0), Real(100.0), Real(1.0), Real(1.0), Real(1.0),
                             Real(1.0e6), Real(1.0e4), Real(1000.0), Real(1.0), Real(1.0e4), Real(1.0e4)};
    ASSERT_EQ(static_cast<int>(cat.size()), nf);
    ASSERT_EQ(WRFSurfaceFields::num_fields, nf);
    for (int n = 0; n < nf; ++n) {
        EXPECT_STREQ(cat[n].nc_name, names[n]);
        EXPECT_STREQ(cat[n].units, units[n]);
        EXPECT_EQ(static_cast<int>(cat[n].id), n);
        EXPECT_EQ(cat[n].lo, Real(0.0));
        EXPECT_EQ(cat[n].hi, hi[n]) << names[n];
        const Kind kind = (n < 9) ? Kind::Surface : ((n < 11) ? Kind::SoilField : Kind::LayerColumn);
        EXPECT_EQ(cat[n].kind, kind) << names[n];
        EXPECT_EQ(WRFSurfaceFields::find(names[n]), &cat[n]);
        EXPECT_EQ(WRFSurfaceFields::has_soil_layers(names[n]), n >= 9) << names[n];
    }
    EXPECT_EQ(WRFSurfaceFields::find("TSK"), nullptr);
    EXPECT_EQ(WRFSurfaceFields::find("albedo"), nullptr);
    EXPECT_FALSE(WRFSurfaceFields::has_soil_layers("SH2O"));   // the land model's own
}

TEST(WRFSurfaceFields, ColumnFieldKeepsEveryColumnAndClampsTheGhosts)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(2);
    const FArrayBox src = make_source(1, [] (int i, int j, int) { return albedo_at(i, j); });
    fill_whole(f, L, Field::ALBEDO, src);
    const MultiFab* mf = f.get(0, Field::ALBEDO);
    ASSERT_NE(mf, nullptr);
    EXPECT_EQ(mf->nComp(), 1);
    EXPECT_EQ(mf->nGrowVect(), L.ng);
    expect_field(*mf, [] (int i, int j, int) { return albedo_at(i, j); });
    EXPECT_TRUE(f.has(0, Field::ALBEDO));
    EXPECT_FALSE(f.has(0, Field::LAI));
    EXPECT_FALSE(f.has(1, Field::ALBEDO));
    EXPECT_EQ(f.get(-1, Field::ALBEDO), nullptr);
    EXPECT_EQ(f.get(2, Field::ALBEDO), nullptr);
    EXPECT_TRUE(f.holds_any(0));
    EXPECT_FALSE(f.holds_any(1));
}

// The source's z index starts at 3 and is nodal in z, as the reader's fab of a soil
// variable is; the first layer of the file is component 0.
TEST(WRFSurfaceFields, SoilFieldHasOneComponentPerLayerTopFirst)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(1);
    const int nlay = 4;
    const FArrayBox src = make_source(0, nx - 1, 0, ny - 1, nlay, 3,
                                      [] (int i, int j, int n) { return soil_at(i, j, n); }, true);
    fill_whole(f, L, Field::TSLB, src, nlay);
    const MultiFab* mf = f.get(0, Field::TSLB);
    ASSERT_NE(mf, nullptr);
    ASSERT_EQ(mf->nComp(), nlay);
    expect_field(*mf, [] (int i, int j, int n) { return soil_at(i, j, n); });
}

// A surface field takes the first z index only, whatever the source holds above it
TEST(WRFSurfaceFields, ColumnFieldIgnoresHigherLayersOfTheSource)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(1);
    const FArrayBox src = make_source(3, [] (int i, int j, int n) { return soil_at(i, j, n); });
    fill_whole(f, L, Field::TMN, src, 3);
    const MultiFab* mf = f.get(0, Field::TMN);
    ASSERT_NE(mf, nullptr);
    ASSERT_EQ(mf->nComp(), 1);
    expect_field(*mf, [] (int i, int j, int) { return soil_at(i, j, 0); });
}

// A level put together from two files: each call gives the columns of one file, from a
// source that holds only those columns.
TEST(WRFSurfaceFields, TwoFilesEachGiveTheirPart)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(1);
    const int nlay = 3;
    const auto val = [] (int i, int j, int n) { return soil_at(i, j, n); };
    const FArrayBox west = make_source(0, 2, 0, ny - 1, nlay, 0, val);
    const FArrayBox east = make_source(3, nx - 1, 0, ny - 1, nlay, 0, val);
    f.fill(0, Field::SMOIS, L.ba2d, L.dm, L.ng, nlay, west, columns(0, 2, 0, ny - 1));
    f.fill(0, Field::SMOIS, L.ba2d, L.dm, L.ng, nlay, east, columns(3, nx - 1, 0, ny - 1));
    // SMOIS is bounded by 1: these values (about 280) are outside, so the field is dropped
    f.finish_level(0, Periodicity::NonPeriodic());
    EXPECT_FALSE(f.has(0, Field::SMOIS));

    f.fill(0, Field::TSLB, L.ba2d, L.dm, L.ng, nlay, west, columns(0, 2, 0, ny - 1));
    f.fill(0, Field::TSLB, L.ba2d, L.dm, L.ng, nlay, east, columns(3, nx - 1, 0, ny - 1));
    f.finish_level(0, Periodicity::NonPeriodic());
    const MultiFab* mf = f.get(0, Field::TSLB);
    ASSERT_NE(mf, nullptr);
    ASSERT_EQ(mf->nComp(), nlay);
    // The ghost columns across the seam hold the other file's values
    expect_field(*mf, val);
}

// The files give only part of the level: the field is reported and not held; the run goes on
TEST(WRFSurfaceFields, FilesThatDoNotCoverTheLevelLeaveTheFieldOut)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(1);
    ScopedThrowOnAbort guard;
    const FArrayBox west = make_source(0, 2, 0, ny - 1, 1, 0, [] (int i, int j, int) { return albedo_at(i, j); });
    EXPECT_NO_THROW(f.fill(0, Field::ALBEDO, L.ba2d, L.dm, L.ng, 1, west, columns(0, 2, 0, ny - 1)));
    EXPECT_NO_THROW(f.finish_level(0, Periodicity::NonPeriodic()));
    EXPECT_FALSE(f.has(0, Field::ALBEDO));
    EXPECT_FALSE(f.holds_any(0));
}

// The file says it covers the level, and the region read does not reach: not read there,
// reported, not held
TEST(WRFSurfaceFields, SourceThatDoesNotReachLeavesTheFieldOut)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(1);
    ScopedThrowOnAbort guard;
    const FArrayBox src = make_source(0, nx - 2, 0, ny - 1, 1, 0, [] (int i, int j, int) { return albedo_at(i, j); });
    EXPECT_NO_THROW(f.fill(0, Field::LAI, L.ba2d, L.dm, L.ng, 1, src, all_columns()));
    EXPECT_NO_THROW(f.finish_level(0, Periodicity::NonPeriodic()));
    EXPECT_FALSE(f.has(0, Field::LAI));
    // A soil source with fewer layers than the file declares
    const FArrayBox thin = make_source(2, [] (int i, int j, int n) { return soil_at(i, j, n); });
    EXPECT_NO_THROW(f.fill(0, Field::TSLB, L.ba2d, L.dm, L.ng, 3, thin, all_columns()));
    EXPECT_NO_THROW(f.finish_level(0, Periodicity::NonPeriodic()));
    EXPECT_FALSE(f.has(0, Field::TSLB));
}

TEST(WRFSurfaceFields, TwoFilesThatDisagreeOnTheLayerCountLeaveTheFieldOut)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(1);
    const auto val = [] (int i, int j, int n) { return soil_at(i, j, n); };
    const FArrayBox west = make_source(0, 2, 0, ny - 1, 3, 0, val);
    const FArrayBox east = make_source(3, nx - 1, 0, ny - 1, 4, 0, val);
    f.fill(0, Field::TSLB, L.ba2d, L.dm, L.ng, 3, west, columns(0, 2, 0, ny - 1));
    f.fill(0, Field::TSLB, L.ba2d, L.dm, L.ng, 4, east, columns(3, nx - 1, 0, ny - 1));
    testing::internal::CaptureStdout();
    f.finish_level(0, Periodicity::NonPeriodic());
    const std::string report = testing::internal::GetCapturedStdout();
    EXPECT_FALSE(f.has(0, Field::TSLB));
    // The reason is the disagreement, not a count of values: the second file is not copied
    EXPECT_NE(report.find("disagree on its number of soil layers"), std::string::npos) << report;
}

// Grids that do not fill their bounding box: a ghost column that no box covers, inside the
// bounding box, takes the nearest column of its own box. One layout leaves the gap in y,
// the other in x.
TEST(WRFSurfaceFields, GhostColumnsNoBoxCoversTakeTheNearestOfTheirBox)
{
    const auto check = [] (const Box& small, const Box& other, int gap_dir) {
        BoxList bl;
        bl.push_back(small);
        bl.push_back(other);
        const BoxArray ba(bl);
        const DistributionMapping dm(ba);
        MultiFab mf(ba, dm, 1, IntVect(1, 1, 0));
        mf.setVal(Real(-1.0));
        for (MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
            const Box vbx = mfi.validbox();
            FArrayBox host(vbx, 1, The_Pinned_Arena());
            const auto h = host.array();
            for (int j = vbx.smallEnd(1); j <= vbx.bigEnd(1); ++j) {
                for (int i = vbx.smallEnd(0); i <= vbx.bigEnd(0); ++i) { h(i, j, 0) = albedo_at(i, j); }
            }
            mf[mfi].copy<RunOn::Host>(host, vbx, 0, vbx, 0, 1);
        }
        WRFSurfaceFields::fill_ghost_cells(mf, Periodicity::NonPeriodic());
        int n_checked = 0;
        for (MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
            if (mfi.validbox() != small) { continue; }
            const Box gbx = mfi.fabbox();
            FArrayBox host(gbx, 1, The_Pinned_Arena());
            Gpu::dtoh_memcpy(host.dataPtr(), mf[mfi].dataPtr(), host.nBytes());
            const auto a = host.const_array();
            // The row (or line) of ghost columns just beyond the small box, in the gap
            const int g = small.bigEnd(gap_dir) + 1;
            const int t = 1 - gap_dir;
            for (int s = small.smallEnd(t); s <= small.bigEnd(t); ++s) {
                const int i = (gap_dir == 0) ? g : s;
                const int j = (gap_dir == 0) ? s : g;
                const int ci = (gap_dir == 0) ? g - 1 : s;
                const int cj = (gap_dir == 0) ? s : g - 1;
                EXPECT_EQ(a(i, j, 0), albedo_at(ci, cj)) << "i=" << i << " j=" << j;
                ++n_checked;
            }
        }
        EXPECT_GT(n_checked, 0);
    };
    // Gap in y: the columns 0..3 x 3..4 belong to no box
    check(columns(0, 3, 0, 2), columns(4, nx - 1, 0, ny - 1), 1);
    // Gap in x: the columns 3..5 x 0..2 belong to no box
    check(columns(0, 2, 0, 2), columns(0, nx - 1, 3, ny - 1), 0);
}

// A value that cannot be the quantity leaves the field out, and does not stop the run
TEST(WRFSurfaceFields, ValuesOutsideTheBoundsLeaveTheFieldOut)
{
    Level L;
    ScopedThrowOnAbort guard;
    // One column of the level holds the value; the others hold a value inside the bounds
    const auto held = [&L] (Field id, Real bad) {
        WRFSurfaceFields f;
        f.resize(1);
        const FArrayBox src = make_source(1, [bad] (int i, int j, int) {
            return (i == 4 && j == 3) ? bad : Real(0.5);
        });
        EXPECT_NO_THROW(fill_whole(f, L, id, src));
        return f.has(0, id);
    };
    EXPECT_FALSE(held(Field::ALBEDO, Real(1.5)));
    EXPECT_FALSE(held(Field::ALBEDO, Real(-0.01)));
    EXPECT_FALSE(held(Field::VEGFRA, Real(100.5)));
    EXPECT_FALSE(held(Field::SEAICE, Real(2.0)));
    EXPECT_FALSE(held(Field::SEAICE, Real(-1.0e-6)));
    EXPECT_FALSE(held(Field::SMOIS, Real(1.2)));
    EXPECT_FALSE(held(Field::TMN, Real(-1.0)));
    EXPECT_FALSE(held(Field::ALBBCK, Real(1.5)));
    EXPECT_FALSE(held(Field::SNOWH, Real(-999.0)));
    EXPECT_FALSE(held(Field::ZNT, std::numeric_limits<Real>::quiet_NaN()));
    EXPECT_FALSE(held(Field::SNOW, std::numeric_limits<Real>::infinity()));
    EXPECT_FALSE(held(Field::SNOWH, -std::numeric_limits<Real>::infinity()));
    // The next number above a bound
    EXPECT_FALSE(held(Field::ALBEDO, std::nextafter(Real(1.0), Real(2.0))));
    EXPECT_FALSE(held(Field::VEGFRA, std::nextafter(Real(100.0), Real(200.0))));
    // A NetCDF fill value is finite: the upper bounds leave it out
    for (const Field id : {Field::LAI, Field::TMN, Field::ZNT, Field::SNOW, Field::SNOWH, Field::TSLB}) {
        EXPECT_FALSE(held(id, Real(9.96921e36)));
    }
    // The bounds themselves pass
    EXPECT_TRUE(held(Field::ALBEDO, Real(1.0)));
    EXPECT_TRUE(held(Field::ALBEDO, Real(0.0)));
    EXPECT_TRUE(held(Field::VEGFRA, Real(100.0)));
    EXPECT_TRUE(held(Field::SNOWH, Real(1.0e4)));
}

// A field left out does not take the others with it
TEST(WRFSurfaceFields, OneFieldLeftOutLeavesTheOthersHeld)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(1);
    const FArrayBox good = make_source(1, [] (int i, int j, int) { return albedo_at(i, j); });
    const FArrayBox bad  = make_source(1, [] (int, int, int) { return Real(-999.0); });
    f.fill(0, Field::ALBBCK, L.ba2d, L.dm, L.ng, 1, good, all_columns());
    f.fill(0, Field::SNOWH,  L.ba2d, L.dm, L.ng, 1, bad,  all_columns());
    f.finish_level(0, Periodicity::NonPeriodic());
    EXPECT_FALSE(f.has(0, Field::SNOWH));
    const MultiFab* mf = f.get(0, Field::ALBBCK);
    ASSERT_NE(mf, nullptr);
    expect_field(*mf, [] (int i, int j, int) { return albedo_at(i, j); });
}

TEST(WRFSurfaceFields, PeriodicGhostsComeFromTheFarSide)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(1);
    const Periodicity period(IntVect(nx, 0, 0));
    const FArrayBox src = make_source(1, [] (int i, int j, int) { return albedo_at(i, j); });
    f.fill(0, Field::ALBEDO, L.ba2d, L.dm, L.ng, 1, src, all_columns());
    f.finish_level(0, period);
    const MultiFab* mf = f.get(0, Field::ALBEDO);
    ASSERT_NE(mf, nullptr);
    for_each_cell(*mf, [] (int i, int j, int, Real v) {
        const int pi = (i + nx) % nx;
        EXPECT_EQ(v, albedo_at(pi, clamp(j, 0, ny - 1))) << "i=" << i << " j=" << j;
    });
}

TEST(WRFSurfaceFields, SoilLayersComeFromTheColumn)
{
    WRFSurfaceFields f;
    f.resize(2);
    f.set_soil_layers(0, Field::ZS,  make_column({Real(0.05), Real(0.25), Real(0.7), Real(1.5)}));
    f.set_soil_layers(0, Field::DZS, make_column({Real(0.1), Real(0.3), Real(0.6), Real(1.0)}));
    ASSERT_EQ(f.soil_depth(0).size(), 4);
    ASSERT_EQ(f.soil_thickness(0).size(), 4);
    EXPECT_EQ(f.soil_depth(0)[0], Real(0.05));
    EXPECT_EQ(f.soil_depth(0)[3], Real(1.5));
    EXPECT_EQ(f.soil_thickness(0)[1], Real(0.3));
    EXPECT_TRUE(f.soil_depth(1).empty());
    EXPECT_TRUE(f.holds_any(0));
    EXPECT_FALSE(f.holds_any(1));
    // Not fields: get() stays null for them
    EXPECT_EQ(f.get(0, Field::ZS), nullptr);
    EXPECT_EQ(f.get(0, Field::DZS), nullptr);
}

// A land model with its first level at the surface gives ZS = 0 there
TEST(WRFSurfaceFields, FirstSoilDepthMayBeZero)
{
    WRFSurfaceFields f;
    f.resize(1);
    f.set_soil_layers(0, Field::ZS, make_column({Real(0.0), Real(0.05), Real(0.2), Real(0.6)}));
    ASSERT_EQ(f.soil_depth(0).size(), 4);
    EXPECT_EQ(f.soil_depth(0)[0], Real(0.0));
    EXPECT_EQ(f.soil_depth(0)[1], Real(0.05));
}

TEST(WRFSurfaceFields, SoilLayerColumnsThatBreakTheRulesAreLeftOut)
{
    ScopedThrowOnAbort guard;
    const auto depth_held = [] (std::initializer_list<Real> v) {
        WRFSurfaceFields f;
        f.resize(1);
        EXPECT_NO_THROW(f.set_soil_layers(0, Field::ZS, make_column(v)));
        return !f.soil_depth(0).empty();
    };
    const auto thick_held = [] (std::initializer_list<Real> v) {
        WRFSurfaceFields f;
        f.resize(1);
        EXPECT_NO_THROW(f.set_soil_layers(0, Field::DZS, make_column(v)));
        return !f.soil_thickness(0).empty();
    };
    EXPECT_TRUE (depth_held({Real(0.05), Real(0.25), Real(0.7)}));
    EXPECT_FALSE(depth_held({Real(0.0), Real(0.0), Real(0.0)}));        // not given
    EXPECT_FALSE(depth_held({Real(0.05), Real(0.05), Real(0.7)}));      // does not increase
    EXPECT_FALSE(depth_held({Real(0.7), Real(0.25), Real(0.05)}));      // bottom layer first
    EXPECT_FALSE(depth_held({Real(-0.05), Real(0.25), Real(0.7)}));     // negative
    EXPECT_FALSE(depth_held({Real(0.0), Real(0.0), Real(0.7)}));        // zero twice
    EXPECT_FALSE(depth_held({Real(0.05), std::numeric_limits<Real>::quiet_NaN(), Real(0.7)}));
    EXPECT_FALSE(depth_held({Real(0.05), Real(0.25), Real(9.96921e36)}));
    EXPECT_TRUE (thick_held({Real(0.1), Real(0.3), Real(0.6)}));
    EXPECT_FALSE(thick_held({Real(0.0), Real(0.0), Real(0.0)}));        // not given
    EXPECT_FALSE(thick_held({Real(0.0), Real(0.3), Real(0.6)}));        // a zero thickness
    EXPECT_FALSE(thick_held({Real(0.1), Real(-0.3), Real(0.6)}));
    EXPECT_FALSE(thick_held({Real(0.1), std::numeric_limits<Real>::infinity(), Real(0.6)}));

    // Two files of a level that disagree: left out; a later file does not bring it back
    WRFSurfaceFields f;
    f.resize(1);
    f.set_soil_layers(0, Field::ZS, make_column({Real(0.05), Real(0.25), Real(0.7)}));
    f.set_soil_layers(0, Field::ZS, make_column({Real(0.05), Real(0.25), Real(0.7)}));
    EXPECT_EQ(f.soil_depth(0).size(), 3);
    f.set_soil_layers(0, Field::ZS, make_column({Real(0.05), Real(0.3), Real(0.7)}));
    EXPECT_TRUE(f.soil_depth(0).empty());
    f.set_soil_layers(0, Field::ZS, make_column({Real(0.05), Real(0.25), Real(0.7)}));
    EXPECT_TRUE(f.soil_depth(0).empty());

    EXPECT_THROW(amrex::ignore_unused(f.soil_depth(1)), std::runtime_error);
}

TEST(WRFSurfaceFields, SoilDataThatDisagreeOnTheLayerCountAreLeftOut)
{
    Level L;
    ScopedThrowOnAbort guard;
    const FArrayBox soil3 = make_source(3, [] (int i, int j, int n) { return soil_at(i, j, n); });
    const FArrayBox alb   = make_source(1, [] (int i, int j, int) { return albedo_at(i, j); });
    {
        // Agreement: all held
        WRFSurfaceFields f;
        f.resize(1);
        f.set_soil_layers(0, Field::ZS,  make_column({Real(0.05), Real(0.25), Real(0.7)}));
        f.set_soil_layers(0, Field::DZS, make_column({Real(0.1), Real(0.3), Real(0.6)}));
        f.fill(0, Field::TSLB, L.ba2d, L.dm, L.ng, 3, soil3, all_columns());
        f.finish_level(0, Periodicity::NonPeriodic());
        EXPECT_TRUE(f.has(0, Field::TSLB));
        EXPECT_EQ(f.soil_depth(0).size(), 3);
        EXPECT_EQ(f.soil_thickness(0).size(), 3);
    }
    {
        // Four depths, three layers of temperature: the four soil data go, the albedo stays
        WRFSurfaceFields f;
        f.resize(1);
        f.set_soil_layers(0, Field::ZS,  make_column({Real(0.05), Real(0.25), Real(0.7), Real(1.5)}));
        f.set_soil_layers(0, Field::DZS, make_column({Real(0.1), Real(0.3), Real(0.6)}));
        f.fill(0, Field::TSLB,   L.ba2d, L.dm, L.ng, 3, soil3, all_columns());
        f.fill(0, Field::ALBBCK, L.ba2d, L.dm, L.ng, 1, alb,   all_columns());
        EXPECT_NO_THROW(f.finish_level(0, Periodicity::NonPeriodic()));
        EXPECT_FALSE(f.has(0, Field::TSLB));
        EXPECT_TRUE(f.soil_depth(0).empty());
        EXPECT_TRUE(f.soil_thickness(0).empty());
        EXPECT_TRUE(f.has(0, Field::ALBBCK));
    }
}

TEST(WRFSurfaceFields, ClearDropsOneLevel)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(2);
    const FArrayBox src = make_source(1, [] (int i, int j, int) { return albedo_at(i, j); });
    for (int lev = 0; lev < 2; ++lev) {
        f.fill(lev, Field::ALBEDO, L.ba2d, L.dm, L.ng, 1, src, all_columns());
        f.finish_level(lev, Periodicity::NonPeriodic());
        f.set_soil_layers(lev, Field::ZS, make_column({Real(0.05), Real(0.25)}));
    }
    f.clear(0);
    EXPECT_FALSE(f.has(0, Field::ALBEDO));
    EXPECT_TRUE(f.soil_depth(0).empty());
    EXPECT_FALSE(f.holds_any(0));
    EXPECT_TRUE(f.has(1, Field::ALBEDO));
    EXPECT_EQ(f.soil_depth(1).size(), 2);
    f.clear(5);    // out of range: nothing
    f.clear(-1);
    EXPECT_TRUE(f.has(1, Field::ALBEDO));
}

// A regrid or a rebalance: the same columns in other boxes keep every value
TEST(WRFSurfaceFields, RemakeCarriesTheFieldsOntoNewGrids)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(1);
    const int nlay = 3;
    const FArrayBox alb  = make_source(1, [] (int i, int j, int) { return albedo_at(i, j); });
    const FArrayBox soil = make_source(nlay, [] (int i, int j, int n) { return soil_at(i, j, n); });
    f.fill(0, Field::ALBEDO, L.ba2d, L.dm, L.ng, 1,    alb,  all_columns());
    f.fill(0, Field::TSLB,   L.ba2d, L.dm, L.ng, nlay, soil, all_columns());
    f.set_soil_layers(0, Field::ZS, make_column({Real(0.05), Real(0.25), Real(0.7)}));
    f.finish_level(0, Periodicity::NonPeriodic());

    const Level N(IntVect(2, 5, 1));
    ASSERT_FALSE(N.ba2d == L.ba2d);
    f.remake(0, N.ba2d, N.dm, Periodicity::NonPeriodic(), IntVect(1), false);
    const MultiFab* a = f.get(0, Field::ALBEDO);
    const MultiFab* t = f.get(0, Field::TSLB);
    ASSERT_NE(a, nullptr);
    ASSERT_NE(t, nullptr);
    EXPECT_TRUE(a->boxArray() == N.ba2d);
    EXPECT_TRUE(t->boxArray() == N.ba2d);
    EXPECT_EQ(a->nGrowVect(), L.ng);
    ASSERT_EQ(t->nComp(), nlay);
    expect_field(*a, [] (int i, int j, int) { return albedo_at(i, j); });
    expect_field(*t, [] (int i, int j, int n) { return soil_at(i, j, n); });
    EXPECT_EQ(f.soil_depth(0).size(), 3);
}

// Level 0 has no coarser level: new grids that reach beyond the grids that held a field
// drop it, with one line; the run goes on and the soil layer columns are kept
TEST(WRFSurfaceFields, RemakeDropsTheFieldsBeyondTheirRegion)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(1);
    ScopedThrowOnAbort guard;
    const FArrayBox alb = make_source(1, [] (int i, int j, int) { return albedo_at(i, j); });
    fill_whole(f, L, Field::ALBEDO, alb);
    f.set_soil_layers(0, Field::ZS, make_column({Real(0.05), Real(0.25)}));
    BoxArray wider(Box(IntVect(0, 0, 0), IntVect(nx + 1, ny - 1, 0)));
    const DistributionMapping dm(wider);
    testing::internal::CaptureStdout();
    EXPECT_NO_THROW(f.remake(0, wider, dm, Periodicity::NonPeriodic(), IntVect(1), true));
    const std::string report = testing::internal::GetCapturedStdout();
    EXPECT_FALSE(f.has(0, Field::ALBEDO));
    EXPECT_NE(report.find("Warning: wrfinput ALBEDO at level 0 is dropped"), std::string::npos) << report;
    ASSERT_EQ(f.soil_depth(0).size(), 2);
    EXPECT_EQ(f.soil_depth(0)[1], Real(0.25));
}

// A regrid of a fine level that moves and grows: the columns its old grids held keep the
// level's own values, every other column takes the coarse column it lies in
TEST(WRFSurfaceFields, RemakeTakesTheNewColumnsFromTheCoarserLevel)
{
    TwoLevels T;
    const Box moved = columns(4, 9, 2, 7);   // old grids: 2..5 x 2..5
    BoxArray nba(moved);
    nba.maxSize(IntVect(4, 3, 1));
    const DistributionMapping ndm(nba);
    testing::internal::CaptureStdout();
    // The ratio in z is ignored: the fields have no vertical extent
    T.f.remake(1, nba, ndm, Periodicity::NonPeriodic(), IntVect(2, 2, 2), false, 1);
    const std::string report = testing::internal::GetCapturedStdout();

    const MultiFab* a = T.f.get(1, Field::ALBEDO);
    const MultiFab* t = T.f.get(1, Field::TSLB);
    ASSERT_NE(a, nullptr);
    ASSERT_NE(t, nullptr);
    EXPECT_TRUE(a->boxArray() == nba);
    EXPECT_EQ(a->nGrowVect(), T.ng);
    ASSERT_EQ(t->nComp(), TwoLevels::nlay);
    // Sets: a ghost cell of one box is a valid cell of its neighbor
    std::set<std::pair<int,int>> n_own;
    std::set<std::pair<int,int>> ncrse;
    for_each_cell(*a, [&] (int i, int j, int, Real v) {
        if (!moved.contains(IntVect(i, j, 0))) { return; }
        const bool own = T.fine_box.contains(IntVect(i, j, 0));
        EXPECT_EQ(v, own ? fine_albedo_at(i, j) : albedo_at(i/2, j/2)) << "i=" << i << " j=" << j;
        (own ? n_own : ncrse).insert({i, j});
    });
    EXPECT_EQ(n_own.size(), 8);
    EXPECT_EQ(ncrse.size(), 28);
    for_each_cell(*t, [&] (int i, int j, int n, Real v) {
        if (!moved.contains(IntVect(i, j, 0))) { return; }
        const bool own = T.fine_box.contains(IntVect(i, j, 0));
        EXPECT_EQ(v, own ? fine_soil_at(i, j, n) : soil_at(i/2, j/2, n)) << "i=" << i << " j=" << j << " n=" << n;
    });
    // The level's own soil layer columns are kept; LAI, which only level 0 holds, is not taken
    ASSERT_EQ(T.f.soil_depth(1).size(), TwoLevels::nlay);
    EXPECT_EQ(T.f.soil_depth(1)[0], Real(0.125));
    EXPECT_FALSE(T.f.has(1, Field::LAI));
    EXPECT_NE(report.find("wrfinput surface field ALBEDO at level 1 on the new grids: 1 component, "
                          "28 of 36 columns from level 0"), std::string::npos) << report;
    EXPECT_EQ(report.find("Warning"), std::string::npos) << report;
}

// A level without a file of its own (made by tagging): every field of the coarser level,
// with its soil layers and soil layer columns. A level with a file takes none it lacks.
TEST(WRFSurfaceFields, ALevelWithoutAFileTakesEveryFieldOfTheCoarserLevel)
{
    for (const bool take : {true, false}) {
        TwoLevels T;
        T.f.clear(1);
        ASSERT_FALSE(T.f.holds_any(1));
        const Box where = columns(6, 11, 0, 5);
        BoxArray nba(where);
        nba.maxSize(IntVect(4, 3, 1));
        const DistributionMapping ndm(nba);
        T.f.remake(1, nba, ndm, Periodicity::NonPeriodic(), IntVect(2, 2, 1), take);
        if (!take) {
            EXPECT_FALSE(T.f.holds_any(1));
            continue;
        }
        const MultiFab* a = T.f.get(1, Field::ALBEDO);
        const MultiFab* l = T.f.get(1, Field::LAI);
        const MultiFab* t = T.f.get(1, Field::TSLB);
        ASSERT_NE(a, nullptr);
        ASSERT_NE(l, nullptr);
        ASSERT_NE(t, nullptr);
        EXPECT_FALSE(T.f.has(1, Field::SMOIS));
        EXPECT_EQ(a->nGrowVect(), T.ng);
        ASSERT_EQ(t->nComp(), TwoLevels::nlay);
        int ncell = 0;
        for_each_cell(*t, [&] (int i, int j, int n, Real v) {
            // Ghost cells included: beyond the level's box they take its nearest column
            const int ci = clamp(i, 6, 11) / 2;
            const int cj = clamp(j, 0, 5) / 2;
            EXPECT_EQ(v, soil_at(ci, cj, n)) << "i=" << i << " j=" << j << " n=" << n;
            ++ncell;
        });
        EXPECT_GT(ncell, 0);
        for_each_cell(*a, [&] (int i, int j, int, Real v) {
            EXPECT_EQ(v, albedo_at(clamp(i, 6, 11) / 2, clamp(j, 0, 5) / 2)) << "i=" << i << " j=" << j;
        });
        ASSERT_EQ(T.f.soil_depth(1).size(), TwoLevels::nlay);
        EXPECT_EQ(T.f.soil_depth(1)[2], Real(0.7));
    }
}

// Grids that grow where the coarser level cannot give a field (it does not hold it, or
// holds another number of soil layers): that field is dropped, the others are carried
TEST(WRFSurfaceFields, RemakeDropsOnlyTheFieldsTheCoarserLevelCannotGive)
{
    TwoLevels T;
    // At level 1 only: SNOW. With two layers where level 0 has three: SMOIS.
    const auto flat = [] (int, int, int) { return Real(0.25); };
    const FArrayBox snow  = make_source(2, 5, 2, 5, 1, 0, flat);
    const FArrayBox smois = make_source(2, 5, 2, 5, 2, 0, flat);
    const FArrayBox smois0 = make_source(TwoLevels::nlay, flat);
    WRFSurfaceFields& f = T.f;
    f.clear(1);
    f.fill(0, Field::SMOIS, T.L.ba2d, T.L.dm, T.ng, TwoLevels::nlay, smois0, all_columns());
    f.finish_level(0, Periodicity::NonPeriodic());
    const FArrayBox alb = make_source(2, 5, 2, 5, 1, 0, [] (int i, int j, int) { return fine_albedo_at(i, j); });
    f.fill(1, Field::ALBEDO, T.fine_ba, T.fine_dm, T.ng, 1, alb,   T.fine_box);
    f.fill(1, Field::SNOW,   T.fine_ba, T.fine_dm, T.ng, 1, snow,  T.fine_box);
    f.fill(1, Field::SMOIS,  T.fine_ba, T.fine_dm, T.ng, 2, smois, T.fine_box);
    f.finish_level(1, Periodicity::NonPeriodic());
    ASSERT_TRUE(f.has(1, Field::SNOW));
    ASSERT_TRUE(f.has(1, Field::SMOIS));

    // Inside the old grids: nothing is needed from level 0, every field is carried
    BoxArray inner(columns(2, 5, 2, 3));
    const DistributionMapping idm(inner);
    f.remake(1, inner, idm, Periodicity::NonPeriodic(), IntVect(2, 2, 1), false);
    EXPECT_TRUE(f.has(1, Field::ALBEDO));
    EXPECT_TRUE(f.has(1, Field::SNOW));
    ASSERT_TRUE(f.has(1, Field::SMOIS));
    EXPECT_EQ(f.get(1, Field::SMOIS)->nComp(), 2);

    // Beyond them
    BoxArray grown(columns(2, 7, 2, 3));
    const DistributionMapping gdm(grown);
    testing::internal::CaptureStdout();
    f.remake(1, grown, gdm, Periodicity::NonPeriodic(), IntVect(2, 2, 1), false);
    const std::string report = testing::internal::GetCapturedStdout();
    EXPECT_TRUE(f.has(1, Field::ALBEDO));
    EXPECT_FALSE(f.has(1, Field::SNOW));
    EXPECT_FALSE(f.has(1, Field::SMOIS));
    EXPECT_NE(report.find("Warning: wrfinput SNOW at level 1 is dropped"), std::string::npos) << report;
    EXPECT_NE(report.find("Warning: wrfinput SMOIS at level 1 is dropped"), std::string::npos) << report;
    EXPECT_EQ(report.find("ALBEDO"), std::string::npos) << report;
}

TEST(WRFSurfaceFields, CheckpointHeaderRoundTrip)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(2);
    const int nlay = 3;
    const FArrayBox alb  = make_source(1, [] (int i, int j, int) { return albedo_at(i, j); });
    const FArrayBox soil = make_source(nlay, [] (int i, int j, int n) { return soil_at(i, j, n); });
    f.fill(0, Field::ALBBCK, L.ba2d, L.dm, L.ng, 1,    alb,  all_columns());
    f.fill(0, Field::TSLB,   L.ba2d, L.dm, L.ng, nlay, soil, all_columns());
    // Not multiples of a power of two: the text must give them back bit for bit
    const Real third = Real(1.0) / Real(3.0);
    f.set_soil_layers(0, Field::ZS,  make_column({Real(0.05), third, Real(0.7)}));
    f.set_soil_layers(0, Field::DZS, make_column({Real(0.1), Real(0.3), third}));
    f.finish_level(0, Periodicity::NonPeriodic());

    std::ostringstream os;
    f.write_header(os, 2);
    std::istringstream is(os.str());
    Vector<WRFSurfaceFields::HeaderLevel> levels;
    ASSERT_TRUE(WRFSurfaceFields::parse_header(is, levels));
    ASSERT_EQ(levels.size(), 2);
    for (int n = 0; n < WRFSurfaceFields::num_fields; ++n) {
        const Field id = static_cast<Field>(n);
        const int expect = (id == Field::ALBBCK) ? 1
                         : ((id == Field::TSLB || id == Field::ZS || id == Field::DZS) ? nlay : 0);
        EXPECT_EQ(levels[0].ncomp[n], expect) << n;
        EXPECT_EQ(levels[1].ncomp[n], 0) << n;
    }
    ASSERT_EQ(levels[0].soil_depth.size(), nlay);
    ASSERT_EQ(levels[0].soil_thickness.size(), nlay);
    EXPECT_EQ(levels[0].soil_depth[1], third);
    EXPECT_EQ(levels[0].soil_thickness[2], third);
    EXPECT_EQ(levels[0].soil_depth[0], Real(0.05));
    EXPECT_TRUE(levels[1].soil_depth.empty());

    // Not a header of this class
    std::istringstream other("HyperCLaw-V1.1\n3\n");
    EXPECT_FALSE(WRFSurfaceFields::parse_header(other, levels));
    EXPECT_TRUE(levels.empty());
    std::istringstream cut("WRFSurfaceFields 1\n1\nLevel 0 2\nALBBCK 1\n");
    EXPECT_FALSE(WRFSurfaceFields::parse_header(cut, levels));
    std::istringstream newer("WRFSurfaceFields 2\n0\n");
    EXPECT_FALSE(WRFSurfaceFields::parse_header(newer, levels));
}

// Write the held fields, read them back on other boxes: every value and ghost cell returns
TEST(WRFSurfaceFields, CheckpointRoundTripOntoOtherGrids)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(1);
    const int nlay = 3;
    const FArrayBox alb  = make_source(1, [] (int i, int j, int) { return albedo_at(i, j); });
    const FArrayBox soil = make_source(nlay, [] (int i, int j, int n) { return soil_at(i, j, n); });
    f.fill(0, Field::ALBBCK, L.ba2d, L.dm, L.ng, 1,    alb,  all_columns());
    f.fill(0, Field::TSLB,   L.ba2d, L.dm, L.ng, nlay, soil, all_columns());
    f.set_soil_layers(0, Field::ZS, make_column({Real(0.0), Real(0.25), Real(0.7)}));
    f.finish_level(0, Periodicity::NonPeriodic());

    const std::string dir = "wrfsfc_gtest_chk";
    amrex::UtilCreateCleanDirectory(dir, false);
    amrex::UtilCreateDirectory(dir + "/Level_0", 0755);
    f.write_checkpoint(dir, 1);

    const Level N(IntVect(3, 2, 1));
    WRFSurfaceFields g;
    g.resize(1);
    ASSERT_TRUE(g.read_checkpoint(dir, 1, {N.ba2d}, {N.dm}, N.ng, {Periodicity::NonPeriodic()}, {}));
    const MultiFab* a = g.get(0, Field::ALBBCK);
    const MultiFab* t = g.get(0, Field::TSLB);
    ASSERT_NE(a, nullptr);
    ASSERT_NE(t, nullptr);
    EXPECT_TRUE(a->boxArray() == N.ba2d);
    ASSERT_EQ(t->nComp(), nlay);
    expect_field(*a, [] (int i, int j, int) { return albedo_at(i, j); });
    expect_field(*t, [] (int i, int j, int n) { return soil_at(i, j, n); });
    EXPECT_FALSE(g.has(0, Field::ALBEDO));
    ASSERT_EQ(g.soil_depth(0).size(), nlay);
    EXPECT_EQ(g.soil_depth(0)[0], Real(0.0));
    EXPECT_EQ(g.soil_depth(0)[2], Real(0.7));
    EXPECT_TRUE(g.soil_thickness(0).empty());

    // A checkpoint written before the fields were stored: nothing is restored
    WRFSurfaceFields h;
    h.resize(1);
    amrex::FileSystem::Remove(dir + "/" + WRFSurfaceFields::checkpoint_header_name);
    testing::internal::CaptureStdout();
    EXPECT_FALSE(h.read_checkpoint(dir, 1, {N.ba2d}, {N.dm}, N.ng, {Periodicity::NonPeriodic()}, {}));
    const std::string report = testing::internal::GetCapturedStdout();
    EXPECT_FALSE(h.holds_any(0));
    EXPECT_NE(report.find("WRFSurfaceFieldsHeader is not in the checkpoint"), std::string::npos) << report;
    EXPECT_NE(report.find("none is restored"), std::string::npos) << report;

    amrex::FileSystem::RemoveAll(dir);
}

// A restart onto grids that reach beyond those of the checkpoint follows the rule of a
// regrid: a fine level takes the missing columns from the level below; level 0, which has
// none below, does not restore the field. No cell is left without a value.
TEST(WRFSurfaceFields, CheckpointOntoGridsBeyondThoseOfTheCheckpoint)
{
    TwoLevels T;
    const std::string dir = "wrfsfc_gtest_chk_beyond";
    amrex::UtilCreateCleanDirectory(dir, false);
    amrex::UtilCreateDirectory(dir + "/Level_0", 0755);
    amrex::UtilCreateDirectory(dir + "/Level_1", 0755);
    T.f.write_checkpoint(dir, 2);
    const Vector<Periodicity> period(2, Periodicity::NonPeriodic());
    const Vector<IntVect> ratio(1, IntVect(2, 2, 1));

    // Level 1 on 0..7 x 2..5 where the checkpoint has 2..5 x 2..5
    {
        const Box wide = columns(0, 7, 2, 5);
        BoxArray wba(wide);
        wba.maxSize(IntVect(3, 4, 1));
        const DistributionMapping wdm(wba);
        WRFSurfaceFields g;
        g.resize(2);
        testing::internal::CaptureStdout();
        ASSERT_TRUE(g.read_checkpoint(dir, 2, {T.L.ba2d, wba}, {T.L.dm, wdm}, T.ng, period, ratio));
        const std::string report = testing::internal::GetCapturedStdout();
        EXPECT_EQ(report.find("Warning"), std::string::npos) << report;
        const MultiFab* t = g.get(1, Field::TSLB);
        ASSERT_NE(t, nullptr);
        std::set<std::pair<int,int>> n_own;
        std::set<std::pair<int,int>> ncrse;
        for_each_cell(*t, [&] (int i, int j, int n, Real v) {
            const int ci = clamp(i, 0, 7);   // ghost cells take the nearest column of the level
            const int cj = clamp(j, 2, 5);
            const bool own = T.fine_box.contains(IntVect(ci, cj, 0));
            EXPECT_EQ(v, own ? fine_soil_at(ci, cj, n) : soil_at(ci/2, cj/2, n)) << "i=" << i << " j=" << j << " n=" << n;
            if (wide.contains(IntVect(i, j, 0))) { (own ? n_own : ncrse).insert({i, j}); }
        });
        EXPECT_EQ(n_own.size(), 16);
        EXPECT_EQ(ncrse.size(), 16);
        EXPECT_EQ(g.soil_depth(1)[0], Real(0.125));
    }

    // Level 0 one column wider than the checkpoint
    {
        BoxArray wider(Box(IntVect(0, 0, 0), IntVect(nx, ny - 1, 0)));
        const DistributionMapping dm(wider);
        WRFSurfaceFields g;
        g.resize(1);
        testing::internal::CaptureStdout();
        ASSERT_TRUE(g.read_checkpoint(dir, 1, {wider}, {dm}, T.ng, {Periodicity::NonPeriodic()}, {}));
        const std::string report = testing::internal::GetCapturedStdout();
        EXPECT_FALSE(g.has(0, Field::ALBEDO));
        EXPECT_FALSE(g.has(0, Field::TSLB));
        EXPECT_NE(report.find("Warning: wrfinput ALBEDO at level 0 is not restored: this run's grids "
                              "reach beyond the grids of the checkpoint"), std::string::npos) << report;
        EXPECT_EQ(g.soil_depth(0).size(), TwoLevels::nlay);
    }

    amrex::FileSystem::RemoveAll(dir);
}
