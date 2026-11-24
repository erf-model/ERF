// Unit tests for the flood coupling (Source/FloodCoupling), without an ERF object.
//
// The atmospheric level is a box of 8 x 6 x 3 cells, 100 m on a side and 50 m
// deep. Its boxes cover the columns i = 2..5, j = 1..4 (an offset refined region),
// in two boxes, so the flood grid has atm_lo = (2, 1). The surface precipitation
// accumulators are MultiFabs on those boxes that the tests set by hand; a test
// engine records what the coupler hands it. Every value is a small integer or a
// quarter of one, so every expected value is exact in single and double precision.

#include <ERF_FloodCoupler.H>
#include <ERF_FloodCouplingParams.H>
#include <ERF_FloodEngineRegistry.H>
#include <ERF_FloodGrid.H>
#include <ERF_SurfacePrecipitation.H>
#include "../MovingBodies/ERF_GTestThrowOnAbort.H"

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include <gtest/gtest.h>

#include <AMReX_Math.H>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

using namespace amrex;

namespace {

constexpr int nx = 8, ny = 6, nz = 3;
constexpr int C  = 2;

Geometry level_geom ()
{
    const Box dom(IntVect(0, 0, 0), IntVect(nx - 1, ny - 1, nz - 1));
    const RealBox rb({Real(0.0), Real(0.0), Real(0.0)},
                     {Real(100.0) * nx, Real(100.0) * ny, Real(50.0) * nz});
    return Geometry(dom, rb, CoordSys::cartesian, {0, 0, 0});
}

// Two boxes over the columns i = 2..5, j = 1..4, each spanning the height
BoxArray level_boxes ()
{
    BoxList bl;
    bl.push_back(Box(IntVect(2, 1, 0), IntVect(3, 4, nz - 1)));
    bl.push_back(Box(IntVect(4, 1, 0), IntVect(5, 4, nz - 1)));
    return BoxArray(std::move(bl));
}

// Set component 0 of mf to f(i, j) at k = 0 and to a value no step should read above it
template <typename F>
void set_k0 (MultiFab& mf, F f)
{
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const auto a = mf.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i, j, k, 0) = (k == 0) ? f(i, j) : Real(-999.0);
        });
    }
}

// Largest |mf comp - f(i, j)| over the valid cells of every rank
template <typename F>
Real max_error (const MultiFab& mf, int comp, F f)
{
    MultiFab err(mf.boxArray(), mf.DistributionMap(), 1, 0);
    for (MFIter mfi(err); mfi.isValid(); ++mfi) {
        const auto e = err.array(mfi);
        const auto a = mf.const_array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            e(i, j, k) = amrex::Math::abs(a(i, j, k, comp) - f(i, j, k));
        });
    }
    return err.norm0(0);
}

/// Records what the coupler hands it; no ERF calls
class RecordingEngine : public FloodEngine
{
public:
    explicit RecordingEngine (std::vector<std::string> names = {}, bool wants_precip = true)
        : m_names(std::move(names)), m_wants(wants_precip) {}
    const char* name () const override { return "gtest_recording"; }
    void verify_prerequisites (const ERF&, int) override {}
    void initialize (const ERF&, int, const FloodGrid&, const MultiFab*) override {}
    void advance (const FloodForcing& f) override
    {
        ++n_advance;
        last_time = f.time;
        last_dt   = f.dt;
        had_precip = (f.precip_k0 != nullptr);
        if (f.precip_k0) {
            precip = std::make_unique<MultiFab>(f.precip_k0->boxArray(), f.precip_k0->DistributionMap(), 2, 0);
            MultiFab::Copy(*precip, *f.precip_k0, 0, 0, 2, 0);
        }
    }
    bool wants_surface_precip () const override { return m_wants; }
    std::vector<std::string> plot_var_names () const override { return m_names; }
    const MultiFab* plot_var (const std::string& n) const override
    {
        return (!m_names.empty() && n == m_names.front()) ? field : nullptr;
    }

    int  n_advance = 0;
    Real last_time = Real(-1.0);
    Real last_dt   = Real(-1.0);
    bool had_precip = false;
    std::unique_ptr<MultiFab> precip;
    const MultiFab* field = nullptr;

private:
    std::vector<std::string> m_names;
    bool m_wants;
};

ERF_REGISTER_FLOOD_ENGINE(RecordingEngine, "gtest_recording")

// A coupler on the offset region with a recording engine; the engine stays reachable
struct Rig
{
    Geometry geom = level_geom();
    BoxArray ba   = level_boxes();
    DistributionMapping dm{ba};
    RecordingEngine* eng = nullptr;
    std::unique_ptr<FloodCoupler> cpl;
    FloodGrid fg;

    explicit Rig (bool wants_precip = true, std::vector<std::string> names = {})
    {
        FloodCouplingParams p;
        p.grid_ratio = C;
        auto e = std::make_unique<RecordingEngine>(std::move(names), wants_precip);
        eng = e.get();
        cpl = std::make_unique<FloodCoupler>(p, std::move(e));
        fg = create_flood_grid(ba, dm, geom, C);
        fg.lev = 1;
        cpl->set_flood_grid(1, fg, ba, dm, geom);
    }
    void advance (double t, double dt, const SurfacePrecipAccumulationSources& s)
    {
        MultiFab u(convert(ba, IntVect(1, 0, 0)), dm, 1, 0); u.setVal(0.0_rt);
        MultiFab v(convert(ba, IntVect(0, 1, 0)), dm, 1, 0); v.setVal(0.0_rt);
        MultiFab z(ba, dm, 1, 0); z.setVal(0.0_rt);
        cpl->advance(t, dt, u, v, z, s);
    }
};

/// Sets erf.flood inputs for one test and removes them afterwards
class ScopedFloodInputs
{
public:
    ScopedFloodInputs () : m_pp("erf.flood") {}
    ~ScopedFloodInputs () { for (const auto& k : m_keys) { m_pp.remove(k); } }
    template <typename T> void add (const std::string& k, const T& v) { m_pp.add(k.c_str(), v); m_keys.push_back(k); }
    ScopedFloodInputs (const ScopedFloodInputs&) = delete;
    ScopedFloodInputs& operator= (const ScopedFloodInputs&) = delete;
private:
    ParmParse m_pp;
    std::vector<std::string> m_keys;
};

} // namespace

// The flood grid over an offset region: atm_lo is the region's corner, the grid
// spans C times its columns from index 0, and its physical extent is the region's.
// A field copied from the atmosphere and averaged back is unchanged.
TEST(FloodGrid, OffsetRegionMapsBothWays)
{
    const Geometry geom = level_geom();
    const BoxArray ba = level_boxes();
    const DistributionMapping dm(ba);
    const FloodGrid fg = create_flood_grid(ba, dm, geom, C);

    EXPECT_EQ(fg.atm_lo[0], 2);
    EXPECT_EQ(fg.atm_lo[1], 1);
    EXPECT_EQ(fg.geom.Domain().smallEnd(0), 0);
    EXPECT_EQ(fg.geom.Domain().length(0), C * 4);
    EXPECT_EQ(fg.geom.Domain().length(1), C * 4);
    EXPECT_EQ(fg.geom.ProbLo(0), Real(200.0));
    EXPECT_EQ(fg.geom.ProbHi(0), Real(600.0));
    EXPECT_EQ(fg.geom.ProbLo(1), Real(100.0));
    EXPECT_EQ(fg.geom.ProbHi(1), Real(500.0));
    EXPECT_EQ(fg.geom.CellSize(0), Real(50.0));
    ASSERT_EQ(fg.ba.size(), ba.size());

    BoxArray ba_k0 = ba; ba_k0.coarsen(IntVect(1, 1, nz));
    MultiFab atm(ba_k0, dm, 1, 0);
    set_k0(atm, [] AMREX_GPU_DEVICE (int i, int j) { return Real(i + 10 * j); });
    MultiFab fl(fg.ba, fg.dm, 1, 0);
    fill_flood_from_atm_k0(fl, atm, fg);
    // flood cell (i, j) lies in column (i / C + 2, j / C + 1)
    EXPECT_EQ(max_error(fl, 0, [] AMREX_GPU_DEVICE (int i, int j, int) {
        return Real(i / C + 2 + 10 * (j / C + 1)); }), Real(0.0));

    MultiFab back(ba_k0, dm, 1, 0);
    coarsen_flood_to_atm(back, fl, fg);
    EXPECT_EQ(max_error(back, 0, [] AMREX_GPU_DEVICE (int i, int j, int) {
        return Real(i + 10 * j); }), Real(0.0));
}

// Each step hands the engine the change of the accumulation since the step
// before (or since the reference), split into liquid and frozen, per column.
// A non-zero reference, as on a restart, is not handed over again.
TEST(FloodCoupler, HandsTheStepPrecipitationLiquidAndFrozen)
{
    Rig r;
    MultiFab rain(r.ba, r.dm, 1, 0), snow(r.ba, r.dm, 1, 0);
    SurfacePrecipAccumulationSources s;
    s.rain = {&rain, Real(1.0)};
    s.snow = {&snow, Real(0.5)};   // native units: half a kg/m2 each

    // Reference: what a checkpoint would hold
    set_k0(rain, [] AMREX_GPU_DEVICE (int i, int j) { return Real(100 + i + j); });
    set_k0(snow, [] AMREX_GPU_DEVICE (int, int)     { return Real(8.0); });
    r.cpl->set_precip_reference(s);

    // Step 1: rain grows by i, snow by 2 native (1 kg/m2)
    set_k0(rain, [] AMREX_GPU_DEVICE (int i, int j) { return Real(100 + 2 * i + j); });
    set_k0(snow, [] AMREX_GPU_DEVICE (int, int)     { return Real(10.0); });
    r.advance(30.0, 5.0, s);
    ASSERT_TRUE(r.eng->had_precip);
    EXPECT_EQ(r.eng->n_advance, 1);
    EXPECT_EQ(r.eng->last_time, Real(30.0));
    EXPECT_EQ(r.eng->last_dt, Real(5.0));
    EXPECT_DOUBLE_EQ(r.cpl->fields_time(), 35.0);
    EXPECT_EQ(max_error(*r.eng->precip, 0, [] AMREX_GPU_DEVICE (int i, int, int) { return Real(i); }), Real(0.0));
    EXPECT_EQ(max_error(*r.eng->precip, 1, [] AMREX_GPU_DEVICE (int, int, int)   { return Real(1.0); }), Real(0.0));

    // Step 2: rain grows by j only; no snow
    set_k0(rain, [] AMREX_GPU_DEVICE (int i, int j) { return Real(100 + 2 * i + 2 * j); });
    r.advance(35.0, 5.0, s);
    EXPECT_EQ(max_error(*r.eng->precip, 0, [] AMREX_GPU_DEVICE (int, int j, int) { return Real(j); }), Real(0.0));
    EXPECT_EQ(max_error(*r.eng->precip, 1, [] AMREX_GPU_DEVICE (int, int, int)   { return Real(0.0); }), Real(0.0));
}

// A scheme with a total and frozen species: liquid is the total less the frozen part
TEST(FloodCoupler, LiquidIsTotalLessFrozen)
{
    Rig r;
    MultiFab total(r.ba, r.dm, 1, 0), graupel(r.ba, r.dm, 1, 0);
    SurfacePrecipAccumulationSources s;
    s.total   = {&total, Real(0.25)};
    s.graupel = {&graupel, Real(1.0)};
    set_k0(total,   [] AMREX_GPU_DEVICE (int, int) { return Real(0.0); });
    set_k0(graupel, [] AMREX_GPU_DEVICE (int, int) { return Real(0.0); });
    r.cpl->set_precip_reference(s);

    set_k0(total,   [] AMREX_GPU_DEVICE (int i, int) { return Real(4 * (i + 1)); });   // i + 1 kg/m2
    set_k0(graupel, [] AMREX_GPU_DEVICE (int, int)   { return Real(0.25); });
    r.advance(0.0, 1.0, s);
    EXPECT_EQ(max_error(*r.eng->precip, 0, [] AMREX_GPU_DEVICE (int i, int, int) { return Real(i) + Real(0.75); }), Real(0.0));
    EXPECT_EQ(max_error(*r.eng->precip, 1, [] AMREX_GPU_DEVICE (int, int, int)   { return Real(0.25); }), Real(0.0));
}

// An engine that does not read the precipitation gets none and needs no reference
TEST(FloodCoupler, NoPrecipitationForAnEngineThatDoesNotWantIt)
{
    Rig r(false);
    r.advance(0.0, 1.0, SurfacePrecipAccumulationSources{});
    EXPECT_EQ(r.eng->n_advance, 1);
    EXPECT_FALSE(r.eng->had_precip);
}

// An engine that reads the precipitation needs the reference first
TEST(FloodCoupler, AdvanceBeforeTheReferenceStops)
{
    Rig r;
    const std::string msg = erf_gtest::abort_message([&] { r.advance(0.0, 1.0, SurfacePrecipAccumulationSources{}); });
    EXPECT_NE(msg.find("set_precip_reference() before the first advance()"), std::string::npos) << msg;
    EXPECT_EQ(r.eng->n_advance, 0);
}

// The ERF plotfile gets the mean of the engine's field over each column at k = 0,
// zero above and on every other level
TEST(FloodCoupler, PlotFieldIsTheColumnMeanAtTheSurface)
{
    Rig r(true, {"flood_test"});
    MultiFab field(r.fg.ba, r.fg.dm, 1, 0);
    for (MFIter mfi(field); mfi.isValid(); ++mfi) {
        const auto a = field.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { a(i, j, k) = Real(i + 10 * j); });
    }
    r.eng->field = &field;

    MultiFab dst(r.ba, r.dm, 2, 0);
    dst.setVal(Real(7.0));
    r.cpl->fill_plot_var("flood_test", dst, 1, 1);
    // column (I, J) holds flood cells 2(I-2)+{0,1}, 2(J-1)+{0,1}: mean 2(I-2)+0.5 + 10 (2(J-1)+0.5)
    EXPECT_EQ(max_error(dst, 1, [] AMREX_GPU_DEVICE (int i, int j, int k) {
        return (k == 0) ? Real(2 * (i - 2)) + Real(0.5) + Real(10) * (Real(2 * (j - 1)) + Real(0.5)) : Real(0.0); }),
        Real(0.0));
    EXPECT_EQ(max_error(dst, 0, [] AMREX_GPU_DEVICE (int, int, int) { return Real(7.0); }), Real(0.0))
        << "another component was written";

    r.cpl->fill_plot_var("flood_test", dst, 1, 0);
    EXPECT_EQ(dst.norm0(1), Real(0.0)) << "a level without the flood grid holds the field";
    r.cpl->fill_plot_var("flood_other", dst, 1, 1);
    EXPECT_EQ(dst.norm0(1), Real(0.0)) << "a field the engine does not have was written as non-zero";
}

// Plot field names must start with "flood_" and be listed once
TEST(FloodCoupler, PlotFieldNamesAreChecked)
{
    for (const auto& bad : {std::vector<std::string>{"depth"},
                            std::vector<std::string>{"flood_"},
                            std::vector<std::string>{"flood_depth", "flood_depth"}}) {
        const std::string msg = erf_gtest::abort_message([&] {
            FloodCoupler cpl(FloodCouplingParams{}, std::make_unique<RecordingEngine>(bad));
        });
        EXPECT_NE(msg.find("[FLOOD]"), std::string::npos) << "names " << bad.front() << ": " << msg;
        EXPECT_NE(msg.find(bad.size() > 1 ? "twice" : "must start with \"flood_\""), std::string::npos) << msg;
    }
    EXPECT_TRUE(erf_gtest::abort_message([&] {
        FloodCoupler cpl(FloodCouplingParams{}, std::make_unique<RecordingEngine>(std::vector<std::string>{"flood_h", "flood_u"}));
    }).empty());
}

// A checkpoint's flood grid record must match this run's level, ratio and columns
TEST(FloodCoupler, CheckpointMismatchNamesTheInputs)
{
    Rig r;
    const Box cols(IntVect(2, 1, 0), IntVect(5, 4, 0));
    EXPECT_TRUE(r.cpl->checkpoint_mismatch(1, C, cols).empty());
    EXPECT_NE(r.cpl->checkpoint_mismatch(0, C, cols).find("erf.flood.anchor_level"), std::string::npos);
    EXPECT_NE(r.cpl->checkpoint_mismatch(1, C + 1, cols).find("erf.flood.grid_ratio = 3"), std::string::npos);
    EXPECT_FALSE(r.cpl->checkpoint_mismatch(1, C, Box(IntVect(2, 1, 0), IntVect(6, 4, 0))).empty());
}

// The registry: sorted names, creation by name, and the two registration errors
TEST(FloodEngineRegistry, NamesCreationAndErrors)
{
    const std::vector<std::string> names = flood_engine_names();
    EXPECT_TRUE(std::is_sorted(names.begin(), names.end()));
    ASSERT_NE(std::find(names.begin(), names.end(), "gtest_recording"), names.end());
    auto e = make_flood_engine("gtest_recording", FloodCouplingParams{});
    ASSERT_NE(e, nullptr);
    EXPECT_STREQ(e->name(), "gtest_recording");
    EXPECT_EQ(make_flood_engine("no_such_engine", FloodCouplingParams{}), nullptr);

    auto creator = [](const FloodCouplingParams&) -> std::unique_ptr<FloodEngine> { return nullptr; };
    EXPECT_NE(erf_gtest::abort_message([&] { register_flood_engine("", creator); }).find("empty name"),
              std::string::npos);
    EXPECT_NE(erf_gtest::abort_message([&] { register_flood_engine("gtest_recording", creator); })
                  .find("two engines registered the name \"gtest_recording\""), std::string::npos);
    EXPECT_EQ(flood_engine_names(), names) << "a refused registration changed the registry";
}

// The coupling inputs: range checks abort naming the key; with several engines
// compiled in, the deck names one
TEST(FloodCouplingParams, RangeChecksAndEngineChoice)
{
    {
        ScopedFloodInputs in;
        in.add("grid_ratio", 0);
        EXPECT_NE(erf_gtest::abort_message([] { FloodCouplingParams::read(); })
                      .find("erf.flood.grid_ratio must be >= 1, got 0"), std::string::npos);
    }
    {
        ScopedFloodInputs in;
        in.add("anchor_level", -2);
        EXPECT_NE(erf_gtest::abort_message([] { FloodCouplingParams::read(); })
                      .find("erf.flood.anchor_level must be -1"), std::string::npos);
    }
    {
        ScopedFloodInputs in;
        in.add("engine", std::string("no_such_engine"));
        EXPECT_NE(erf_gtest::abort_message([] { FloodCouplingParams::read(); })
                      .find("erf.flood.engine = \"no_such_engine\" is not one of"), std::string::npos);
    }
    {
        ScopedFloodInputs in;
        in.add("enable", true);
        in.add("engine", std::string("gtest_recording"));
        in.add("grid_ratio", 3);
        FloodCouplingParams p;
        EXPECT_TRUE(erf_gtest::abort_message([&] { p = FloodCouplingParams::read(); }).empty());
        EXPECT_TRUE(p.enable);
        EXPECT_EQ(p.engine, "gtest_recording");
        EXPECT_EQ(p.grid_ratio, 3);
        EXPECT_EQ(p.anchor_level, -1);
    }
    if (flood_engine_names().size() > 1) {
        ScopedFloodInputs in;
        in.add("enable", true);
        EXPECT_NE(erf_gtest::abort_message([] { FloodCouplingParams::read(); })
                      .find("several flood engines are compiled in"), std::string::npos);
    }
}
