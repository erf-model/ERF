#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Periodicity.H>

#include <gtest/gtest.h>

#include <ERF_WRFSurfaceFields.H>

// The holder of the wrfinput land and soil fields that no land model takes
// (WRFSurfaceFields). A level of 6 x 5 columns in boxes of at most 4 x 3, one ghost
// column, not periodic. The source carries a different value in every column and soil
// layer, so a swapped index, a reversed layer order or a missed clamp changes the answer.

namespace {

using namespace amrex;
using Field = WRFSurfaceFields::Field;

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
    Level ()
    {
        ba2d = BoxArray(domain);
        ba2d.maxSize(IntVect(4, 3, 1));
        dm = DistributionMapping(ba2d);
    }
};

// Multiples of 1/64: exact in single and double precision, so the comparisons below are exact
Real albedo_at (int i, int j) { return Real(0.015625) * Real(i) + Real(0.125) * Real(j); }
Real soil_at (int i, int j, int n) { return Real(280.0) + Real(n) + Real(0.015625) * Real(i) + Real(0.125) * Real(j); }

// The columns the file covers, with nlay values in z, filled on the host.
template <typename F>
FArrayBox make_source (int nlay, F&& value)
{
    const Box bx(IntVect(0, 0, 0), IntVect(nx - 1, ny - 1, nlay - 1));
    FArrayBox src(bx, 1, The_Pinned_Arena());
    const auto a = src.array();
    for (int k = 0; k < nlay; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) { a(i, j, k) = value(i, j, k); }
        }
    }
    return src;
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

} // namespace

TEST(WRFSurfaceFields, CatalogNamesAndBounds)
{
    const auto& cat = WRFSurfaceFields::catalog();
    constexpr int nf = 11;
    const char* names[nf] = {"LAI", "VEGFRA", "TMN", "ZNT", "ALBEDO", "ALBBCK", "SEAICE", "SNOW", "SNOWH", "TSLB", "SMOIS"};
    const char* units[nf] = {"m2/m2", "%", "K", "m", "1", "1", "1", "kg/m2", "m", "K", "m3/m3"};
    const Real  hi[nf]    = {Real(100.0), Real(100.0), Real(1000.0), Real(100.0), Real(1.0), Real(1.0), Real(1.0),
                             Real(1.0e6), Real(1.0e4), Real(1000.0), Real(1.0)};
    ASSERT_EQ(static_cast<int>(cat.size()), nf);
    ASSERT_EQ(WRFSurfaceFields::num_fields, nf);
    for (int n = 0; n < nf; ++n) {
        EXPECT_STREQ(cat[n].nc_name, names[n]);
        EXPECT_STREQ(cat[n].units, units[n]);
        EXPECT_EQ(static_cast<int>(cat[n].id), n);
        EXPECT_EQ(cat[n].lo, Real(0.0));
        EXPECT_EQ(cat[n].hi, hi[n]) << names[n];
        EXPECT_EQ(cat[n].soil, n >= 9);
        EXPECT_EQ(WRFSurfaceFields::find(names[n]), &cat[n]);
    }
    EXPECT_EQ(WRFSurfaceFields::find("TSK"), nullptr);
    EXPECT_EQ(WRFSurfaceFields::find("ZS"), nullptr);
    EXPECT_EQ(WRFSurfaceFields::find("albedo"), nullptr);
}

// A nested level: the file covers columns 2..5 x 1..3 of a 6 x 5 level, the source's z
// index starts at 3, and the source holds more columns than the file covers. A ghost or
// valid column outside the file's range takes the nearest column inside it, not the
// nearest column of the source.
TEST(WRFSurfaceFields, OffsetFileRangeAndOffsetLayerIndex)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(1);
    const int nlay = 3;
    const int k0 = 3;
    const Box bx(IntVect(-1, -1, k0), IntVect(nx, ny, k0 + nlay - 1));
    FArrayBox src(bx, 1, The_Pinned_Arena());
    const auto a = src.array();
    for (int k = k0; k < k0 + nlay; ++k) {
        for (int j = -1; j <= ny; ++j) {
            for (int i = -1; i <= nx; ++i) { a(i, j, k) = soil_at(i, j, k - k0); }
        }
    }
    const int i_lo = 2, i_hi = 5, j_lo = 1, j_hi = 3;
    f.fill(0, Field::TSLB, L.ba2d, L.dm, L.ng, Periodicity::NonPeriodic(), src, i_lo, i_hi, j_lo, j_hi);
    const MultiFab* mf = f.get(0, Field::TSLB);
    ASSERT_NE(mf, nullptr);
    ASSERT_EQ(mf->nComp(), nlay);
    for_each_cell(*mf, [&] (int i, int j, int n, Real v) {
        EXPECT_EQ(v, soil_at(clamp(i, i_lo, i_hi), clamp(j, j_lo, j_hi), n)) << i << " " << j << " " << n;
    });
}

// On a periodic level the ghost columns across the seam hold the far side's values.
TEST(WRFSurfaceFields, PeriodicGhostsComeFromTheFarSide)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(1);
    const FArrayBox src = make_source(1, [] (int i, int j, int) { return albedo_at(i, j); });
    const Periodicity period(IntVect(nx, ny, 0));
    f.fill(0, Field::ALBEDO, L.ba2d, L.dm, L.ng, period, src, 0, nx - 1, 0, ny - 1);
    const MultiFab* mf = f.get(0, Field::ALBEDO);
    ASSERT_NE(mf, nullptr);
    for_each_cell(*mf, [&] (int i, int j, int, Real v) {
        EXPECT_EQ(v, albedo_at((i + nx) % nx, (j + ny) % ny)) << i << " " << j;
    });
}

TEST(WRFSurfaceFields, ColumnFieldKeepsEveryColumnAndClampsTheGhosts)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(2);
    EXPECT_EQ(f.get(0, Field::ALBEDO), nullptr);
    EXPECT_FALSE(f.has(0, Field::ALBEDO));

    const FArrayBox src = make_source(1, [] (int i, int j, int) { return albedo_at(i, j); });
    f.fill(0, Field::ALBEDO, L.ba2d, L.dm, L.ng, Periodicity::NonPeriodic(), src, 0, nx - 1, 0, ny - 1);

    const MultiFab* mf = f.get(0, Field::ALBEDO);
    ASSERT_NE(mf, nullptr);
    EXPECT_TRUE(f.has(0, Field::ALBEDO));
    EXPECT_EQ(mf->nComp(), 1);
    EXPECT_EQ(mf->nGrowVect(), L.ng);
    EXPECT_TRUE(mf->boxArray() == L.ba2d);
    // The other fields and the other level are untouched
    EXPECT_EQ(f.get(0, Field::LAI), nullptr);
    EXPECT_EQ(f.get(1, Field::ALBEDO), nullptr);
    EXPECT_EQ(f.get(2, Field::ALBEDO), nullptr);   // beyond the levels held

    int cells = 0;
    for_each_cell(*mf, [&] (int i, int j, int, Real v) {
        EXPECT_EQ(v, albedo_at(clamp(i, 0, nx - 1), clamp(j, 0, ny - 1))) << i << " " << j;
        ++cells;
    });
    if (ParallelDescriptor::NProcs() == 1) { EXPECT_GT(cells, nx * ny); }
}

TEST(WRFSurfaceFields, SoilFieldHasOneComponentPerLayerTopFirst)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(1);
    const int nlay = 4;
    const FArrayBox src = make_source(nlay, [] (int i, int j, int k) { return soil_at(i, j, k); });
    f.fill(0, Field::TSLB, L.ba2d, L.dm, L.ng, Periodicity::NonPeriodic(), src, 0, nx - 1, 0, ny - 1);

    const MultiFab* mf = f.get(0, Field::TSLB);
    ASSERT_NE(mf, nullptr);
    ASSERT_EQ(mf->nComp(), nlay);
    for_each_cell(*mf, [&] (int i, int j, int n, Real v) {
        EXPECT_EQ(v, soil_at(clamp(i, 0, nx - 1), clamp(j, 0, ny - 1), n)) << i << " " << j << " " << n;
    });
}

// A field with one value per column takes the source's lowest z index only.
TEST(WRFSurfaceFields, ColumnFieldIgnoresHigherLayersOfTheSource)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(1);
    const FArrayBox src = make_source(3, [] (int i, int j, int k) { return soil_at(i, j, k); });
    f.fill(0, Field::TMN, L.ba2d, L.dm, L.ng, Periodicity::NonPeriodic(), src, 0, nx - 1, 0, ny - 1);
    const MultiFab* mf = f.get(0, Field::TMN);
    ASSERT_NE(mf, nullptr);
    ASSERT_EQ(mf->nComp(), 1);
    for_each_cell(*mf, [&] (int i, int j, int, Real v) {
        EXPECT_EQ(v, soil_at(clamp(i, 0, nx - 1), clamp(j, 0, ny - 1), 0));
    });
}

TEST(WRFSurfaceFields, ValuesOutsideTheBoundsStopTheRun)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(1);
    ScopedThrowOnAbort guard;
    const auto fill = [&] (Field id, Real value_at_one_column) {
        const FArrayBox src = make_source(1, [=] (int i, int j, int) {
            return (i == 3 && j == 2) ? value_at_one_column : Real(0.5);
        });
        f.fill(0, id, L.ba2d, L.dm, L.ng, Periodicity::NonPeriodic(), src, 0, nx - 1, 0, ny - 1);
    };
    // At a bound, and a field left at zero: accepted
    EXPECT_NO_THROW(fill(Field::ALBEDO, Real(1.0)));
    EXPECT_NO_THROW(fill(Field::ALBEDO, Real(0.0)));
    EXPECT_NO_THROW(fill(Field::VEGFRA, Real(100.0)));
    EXPECT_NO_THROW(fill(Field::SNOW, Real(5000.0)));
    EXPECT_NO_THROW(fill(Field::TMN, Real(1000.0)));
    EXPECT_NO_THROW(fill(Field::ALBBCK, Real(1.0)));
    // One column outside
    EXPECT_THROW(fill(Field::ALBEDO, Real(1.5)), std::runtime_error);
    EXPECT_THROW(fill(Field::ALBEDO, Real(-0.01)), std::runtime_error);
    EXPECT_THROW(fill(Field::VEGFRA, Real(100.5)), std::runtime_error);
    EXPECT_THROW(fill(Field::SEAICE, Real(2.0)), std::runtime_error);
    EXPECT_THROW(fill(Field::SMOIS, Real(1.2)), std::runtime_error);
    EXPECT_THROW(fill(Field::TMN, Real(-1.0)), std::runtime_error);
    EXPECT_THROW(fill(Field::ALBBCK, Real(1.5)), std::runtime_error);
    EXPECT_THROW(fill(Field::ZNT, std::numeric_limits<Real>::quiet_NaN()), std::runtime_error);
    EXPECT_THROW(fill(Field::SNOW, std::numeric_limits<Real>::infinity()), std::runtime_error);
    EXPECT_THROW(fill(Field::SNOWH, -std::numeric_limits<Real>::infinity()), std::runtime_error);
    // The next value above a bound
    EXPECT_THROW(fill(Field::ALBEDO, std::nextafter(Real(1.0), Real(2.0))), std::runtime_error);
    EXPECT_THROW(fill(Field::VEGFRA, std::nextafter(Real(100.0), Real(200.0))), std::runtime_error);
    // The NetCDF fill value is finite: the upper bounds stop it
    for (const Field id : {Field::LAI, Field::TMN, Field::ZNT, Field::SNOW, Field::SNOWH, Field::TSLB}) {
        EXPECT_THROW(fill(id, Real(9.96921e36)), std::runtime_error);
    }
}

// The source must cover the rank's columns; a region that is one column short stops the
// run. The count of uncovered cells is reduced before the stop, so every rank throws.
TEST(WRFSurfaceFields, SourceThatDoesNotCoverTheColumnsStopsTheRun)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(1);
    ScopedThrowOnAbort guard;
    const Box short_box(IntVect(0, 0, 0), IntVect(nx - 2, ny - 1, 0));
    FArrayBox src(short_box, 1, The_Pinned_Arena());
    src.setVal<RunOn::Host>(Real(0.5));
    EXPECT_THROW(f.fill(0, Field::LAI, L.ba2d, L.dm, L.ng, Periodicity::NonPeriodic(), src, 0, nx - 1, 0, ny - 1),
                 std::runtime_error);
}

TEST(WRFSurfaceFields, ClearDropsOneLevel)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(2);
    const FArrayBox src = make_source(1, [] (int i, int j, int) { return albedo_at(i, j); });
    f.fill(0, Field::ALBEDO, L.ba2d, L.dm, L.ng, Periodicity::NonPeriodic(), src, 0, nx - 1, 0, ny - 1);
    f.fill(1, Field::ALBEDO, L.ba2d, L.dm, L.ng, Periodicity::NonPeriodic(), src, 0, nx - 1, 0, ny - 1);
    f.clear(0);
    EXPECT_EQ(f.get(0, Field::ALBEDO), nullptr);
    EXPECT_NE(f.get(1, Field::ALBEDO), nullptr);
    f.clear(7);   // a level that is not held: nothing happens
    EXPECT_NE(f.get(1, Field::ALBEDO), nullptr);
}

TEST(WRFSurfaceFields, SoilLayersComeFromTheColumn)
{
    WRFSurfaceFields f;
    f.resize(1);
    EXPECT_TRUE(f.soil_depth(0).empty());
    EXPECT_TRUE(f.soil_thickness(0).empty());

    // The Noah layers: thickness 0.1, 0.3, 0.6, 1.0 m; centres 0.05, 0.25, 0.7, 1.5 m.
    // The column sits away from the origin, as it does on a nested level.
    const Real dzs[4] = {Real(0.1), Real(0.3), Real(0.6), Real(1.0)};
    const Real zs[4]  = {Real(0.05), Real(0.25), Real(0.7), Real(1.5)};
    const Box col(IntVect(12, 7, 0), IntVect(12, 7, 3));
    FArrayBox thick(col, 1, The_Pinned_Arena());
    FArrayBox depth(col, 1, The_Pinned_Arena());
    for (int k = 0; k < 4; ++k) {
        thick.array()(12, 7, k) = dzs[k];
        depth.array()(12, 7, k) = zs[k];
    }
    f.set_soil_layers(0, true, thick);
    f.set_soil_layers(0, false, depth);
    ASSERT_EQ(f.soil_thickness(0).size(), 4);
    ASSERT_EQ(f.soil_depth(0).size(), 4);
    for (int k = 0; k < 4; ++k) {
        EXPECT_EQ(f.soil_thickness(0)[k], dzs[k]);
        EXPECT_EQ(f.soil_depth(0)[k], zs[k]);
    }

    EXPECT_NO_THROW(f.check_soil_layers(0));

    ScopedThrowOnAbort guard;
    // One layer at zero among positive ones
    thick.array()(12, 7, 2) = Real(0.0);
    EXPECT_THROW(f.set_soil_layers(0, true, thick), std::runtime_error);
    // Depths that do not increase
    depth.array()(12, 7, 2) = Real(0.2);
    EXPECT_THROW(f.set_soil_layers(0, false, depth), std::runtime_error);
    // Two layers at the same depth
    depth.array()(12, 7, 2) = zs[1];
    EXPECT_THROW(f.set_soil_layers(0, false, depth), std::runtime_error);
    depth.array()(12, 7, 1) = std::numeric_limits<Real>::quiet_NaN();
    EXPECT_THROW(f.set_soil_layers(0, false, depth), std::runtime_error);
    // Every layer at zero: taken as not given
    for (int k = 0; k < 4; ++k) { thick.array()(12, 7, k) = Real(0.0); }
    EXPECT_NO_THROW(f.set_soil_layers(0, true, thick));
    EXPECT_TRUE(f.soil_thickness(0).empty());
    // A level that is not held
    EXPECT_THROW(amrex::ignore_unused(f.soil_depth(1)), std::runtime_error);

    f.clear(0);
    EXPECT_TRUE(f.soil_depth(0).empty());
    EXPECT_TRUE(f.soil_thickness(0).empty());
}

// TSLB, SMOIS, ZS and DZS, whichever are read, must agree on the number of layers.
TEST(WRFSurfaceFields, SoilDataMustAgreeOnTheLayerCount)
{
    Level L;
    WRFSurfaceFields f;
    f.resize(1);
    EXPECT_NO_THROW(f.check_soil_layers(0));   // nothing read

    const auto column = [] (int nlay) {
        FArrayBox c(Box(IntVect(0, 0, 0), IntVect(0, 0, nlay - 1)), 1, The_Pinned_Arena());
        for (int k = 0; k < nlay; ++k) { c.array()(0, 0, k) = Real(0.25) * Real(k + 1); }
        return c;
    };
    const auto soil = [] (int nlay) {
        return make_source(nlay, [] (int, int, int) { return Real(0.25); });
    };
    f.fill(0, Field::SMOIS, L.ba2d, L.dm, L.ng, Periodicity::NonPeriodic(), soil(4), 0, nx - 1, 0, ny - 1);
    EXPECT_NO_THROW(f.check_soil_layers(0));   // one of the four alone
    f.set_soil_layers(0, false, column(4));
    f.set_soil_layers(0, true, column(4));
    f.fill(0, Field::TSLB, L.ba2d, L.dm, L.ng, Periodicity::NonPeriodic(), soil(4), 0, nx - 1, 0, ny - 1);
    EXPECT_NO_THROW(f.check_soil_layers(0));

    ScopedThrowOnAbort guard;
    f.set_soil_layers(0, true, column(3));
    EXPECT_THROW(f.check_soil_layers(0), std::runtime_error);
    f.set_soil_layers(0, true, column(4));
    f.fill(0, Field::TSLB, L.ba2d, L.dm, L.ng, Periodicity::NonPeriodic(), soil(3), 0, nx - 1, 0, ny - 1);
    EXPECT_THROW(f.check_soil_layers(0), std::runtime_error);
}
