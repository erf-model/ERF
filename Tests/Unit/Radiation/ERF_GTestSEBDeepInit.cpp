#include <cmath>
#include <limits>
#include <string>
#include <utility>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>
#include <AMReX_Utility.H>

#include <gtest/gtest.h>

#include <ERF_Constants.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_LandSurface.H>
#include <ERF_RadStruct.H>
#include <ERF_TwoStreamRadiation.H>

// The deep (restoring) temperature of the force-restore balance as carried state, with
// erf.radiation.seb_init_t_deep_from_file. ERF::init_seb_t_deep_from_file writes the land
// columns from the initial file's TMN once, with the copy that ERF_GTestSEBSkinInit.cpp
// tests; what is tested here is that the balance then keeps the value: it is not reset to
// erf.radiation.seb_t_deep_default each step, the skin restores toward it, a new level
// takes it from its parent (seb_prognostic_fields) and the checkpoint holds it.
//
// A 2 x 1 level of 100 m columns. The balance runs on the scalar flux defaults, so the two
// columns differ only in their deep temperature.

namespace {

using namespace amrex;

constexpr Real t_skin0 = Real(300.0);
constexpr Real t_deep_default = Real(300.0);
constexpr Real t_deep_west = Real(300.0);   // column 0: the default's value
constexpr Real t_deep_east = Real(295.0);   // column 1
constexpr Real dt = Real(1000.0);

BoxArray collapse_z (const BoxArray& ba)
{
    BoxList boxes = ba.boxList();
    for (auto& box : boxes) { box.setRange(2, 0); }
    return BoxArray(std::move(boxes));
}

Real value_at (const MultiFab& mf, int ic)
{
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    const Box point_box(IntVect(ic, 0, 0), IntVect(ic, 0, 0));
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box overlap = point_box & mfi.validbox();
        if (overlap.isEmpty()) { continue; }
        const auto array = mf.const_array(mfi);
        reduce_op.eval(overlap, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> GpuTuple<Real> { return {array(i, j, k)}; });
    }
    Gpu::streamSynchronize();
    Real value = amrex::get<0>(reduce_data.value());
    ParallelDescriptor::ReduceRealSum(value);
    return value;
}

// One column's value of a one-component field on the valid region.
void set_column (MultiFab& mf, int ic, Real value)
{
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const auto a = mf.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (i == ic) { a(i, j, k) = value; }
        });
    }
    Gpu::streamSynchronize();
}

RadChoice balance (bool t_deep_from_file)
{
    RadChoice rad;
    rad.enabled = true;
    rad.sw_enabled = false;
    rad.lw_enabled = true;
    rad.fixed_solar_zenith_angle = Real(0.5);
    rad.fixed_total_solar_irradiance = Real(1360.9);
    rad.rad_t_sfc = t_skin0;
    rad.surface_emissivity_lw = Real(1.0);
    rad.seb_enable = true;
    rad.seb_use_radiation_fluxes = false;
    rad.seb_prognostic_enable = true;
    rad.seb_sw_flux_default = Real(200.0);
    rad.seb_lw_flux_default = Real(0.0);
    rad.seb_hfx_default = Real(0.0);
    rad.seb_lh_default = Real(0.0);
    rad.seb_grdflx_default = Real(0.0);
    rad.seb_t_deep_default = t_deep_default;
    rad.seb_surface_heat_capacity = Real(20000.0);
    rad.seb_turbulent_flux_source = SEBTurbulentFluxSource::Defaults;
    rad.seb_init_t_deep_from_file = t_deep_from_file;
    return rad;
}

struct Level
{
    Box domain{IntVect(0, 0, 0), IntVect(1, 0, 3)};
    RealBox real_box{{0.0, 0.0, 0.0}, {200.0, 100.0, 400.0}};
    int is_periodic[3] = {1, 1, 0};
    Geometry geom{domain, &real_box, 0, is_periodic};
    BoxArray ba{domain};
    BoxArray ba2d{collapse_z(ba)};
    DistributionMapping dm{ba};
};

struct Deep { Real west, east; };
struct StepResult { Deep skin, deep; };

// The deep temperature the test writes in place of ERF::init_seb_t_deep_from_file. With
// the input off the class hands out no writable deep temperature, so the test takes the
// read-only one: the balance must then overwrite what is written here.
void write_deep (TwoStreamRadiation& radiation)
{
    MultiFab* state_deep = radiation.seb_t_deep_state(0);
    MultiFab* deep = state_deep ? state_deep : const_cast<MultiFab*>(radiation.seb_t_deep(0));
    AMREX_ALWAYS_ASSERT(deep != nullptr);
    set_column(*deep, 0, t_deep_west);
    set_column(*deep, 1, t_deep_east);
}

// One step of the balance from a skin of 300 K with the deep temperatures of write_deep.
StepResult one_step (bool t_deep_from_file)
{
    const Level L{};
    const RadChoice rad = balance(t_deep_from_file);

    MultiFab state(L.ba, L.dm, RhoQ2_comp + 1, 0);
    state.setVal(Real(0.0));
    state.setVal(Real(1.0), Rho_comp, 1);
    state.setVal(getThgivenRandT(Real(1.0), Real(290.0), RdoCp), RhoTheta_comp, 1);

    TwoStreamRadiation radiation;
    radiation.resize(1);
    radiation.define_level(0, rad, RdoCp, L.ba2d, L.dm, L.ba, L.domain);
    write_deep(radiation);

    LandSurface lsm;
    lsm.ReSize(1);
    lsm.SetModel<NullSurf>();
    MultiFab qheating(L.ba, L.dm, 2, 0);
    MultiFab rad_fluxes(convert(L.ba, IntVect(0, 0, 1)), L.dm, 4, 0);
    const Vector<const MultiFab*> radiation_inputs {};

    radiation.advance(0, 1, Real(0.0), dt, "pre_dycore",
                      state, nullptr, L.geom, lsm, radiation_inputs, false,
                      &qheating, &rad_fluxes, nullptr, nullptr, nullptr, nullptr, nullptr,
                      0.0, false);
    radiation.advance(0, 1, dt, dt, "post_dycore",
                      state, nullptr, L.geom, lsm, radiation_inputs, false,
                      nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                      0.0, false);

    const MultiFab* skin = radiation.prognostic_surface_temperature_state(0);
    const MultiFab* deep = radiation.seb_t_deep(0);
    AMREX_ALWAYS_ASSERT(skin != nullptr && deep != nullptr);
    return {{value_at(*skin, 0), value_at(*skin, 1)}, {value_at(*deep, 0), value_at(*deep, 1)}};
}

} // namespace

TEST(SEBDeepInit, DefaultIsOff)
{
    const RadChoice rad{};
    EXPECT_FALSE(rad.seb_init_t_deep_from_file);
}

// With the input the deep temperature is one of the fields a new level takes from its
// parent, the two-way coupling averages down and RemakeLevel carries over.
TEST(SEBDeepInit, DeepTemperatureIsCarriedStateOnlyWithTheInput)
{
    const Level L{};
    for (const bool on : {false, true}) {
        const RadChoice rad = balance(on);
        TwoStreamRadiation radiation;
        radiation.resize(1);
        radiation.define_level(0, rad, RdoCp, L.ba2d, L.dm, L.ba, L.domain);
        EXPECT_EQ(radiation.seb_t_deep_is_state(), on);
        EXPECT_EQ(radiation.seb_t_deep_state(0) != nullptr, on);
        const Vector<MultiFab*> fields = radiation.seb_prognostic_fields(0);
        // The skin temperature and moisture, and with the input the deep temperature.
        ASSERT_EQ(static_cast<int>(fields.size()), on ? 3 : 2);
        if (on) { EXPECT_EQ(fields[2], radiation.seb_t_deep_state(0)); }
        // define_level leaves the default in every column either way.
        EXPECT_EQ(radiation.seb_t_deep(0)->min(0), t_deep_default);
        EXPECT_EQ(radiation.seb_t_deep(0)->max(0), t_deep_default);
    }
}

// Without the input the balance sets the deep temperature to the default before each
// step, as it did. With it the columns keep their own and the skin restores toward it:
// the balance is explicit, T_s += dt (R / C_s - (2 pi / tau)(T_s - T_deep)), and the two
// columns share R, C_s and the first T_s, so their skins differ by
// dt (2 pi / tau)(T_deep,west - T_deep,east) after one step.
TEST(SEBDeepInit, BalanceKeepsThePerColumnDeepTemperatureAndRestoresTowardIt)
{
    const StepResult off = one_step(false);
    EXPECT_EQ(off.deep.west, t_deep_default);
    EXPECT_EQ(off.deep.east, t_deep_default);
    EXPECT_EQ(off.skin.west, off.skin.east);

    const StepResult on = one_step(true);
    EXPECT_EQ(on.deep.west, t_deep_west);
    EXPECT_EQ(on.deep.east, t_deep_east);
    // The west column holds the default's value: the run without the input, bit for bit.
    EXPECT_EQ(on.skin.west, off.skin.west);
    const Real tau = RadChoice().seb_restore_timescale_s;
    const Real expected = dt * (Real(2.0) * PI / tau) * (t_deep_west - t_deep_east);
    EXPECT_GT(expected, Real(0.36));   // 0.3636 K: far above the tolerance below
    // The two skins are near 300 K: the tolerance is a few of their rounding errors.
    const Real tol = Real(16.0) * std::numeric_limits<Real>::epsilon() * t_skin0;
    EXPECT_NEAR(on.skin.west - on.skin.east, expected, tol);
}

// The checkpoint holds the deep temperature with the input and not without it. A
// checkpoint without it leaves the default and says so; one with it, read without the
// input, is not read and the run says so.
TEST(SEBDeepInit, CheckpointHoldsTheDeepTemperatureWithTheInput)
{
    const Level L{};
    const Periodicity period = L.geom.periodicity();
    const std::string with = "unit_seb_deep_init_chk_with";
    const std::string without = "unit_seb_deep_init_chk_without";
    for (const auto& chk : {with, without}) {
        amrex::UtilCreateCleanDirectory(chk, true);
        amrex::UtilCreateCleanDirectory(chk + "/Level_0", true);
    }
    const std::string deep_file = "/Level_0/TwoStream_TDeep_H";

    for (const bool on : {true, false}) {
        const RadChoice rad = balance(on);
        TwoStreamRadiation written;
        written.resize(1);
        written.define_level(0, rad, RdoCp, L.ba2d, L.dm, L.ba, L.domain);
        write_deep(written);
        written.write_checkpoint(0, on ? with : without);
    }
    EXPECT_TRUE(amrex::FileExists(with + deep_file));
    EXPECT_FALSE(amrex::FileExists(without + deep_file));
    // The deep moisture belongs to erf.radiation.seb_deep_prognostic_enable alone.
    EXPECT_FALSE(amrex::FileExists(with + "/Level_0/TwoStream_QDeep_H"));

    struct Case { bool on; const std::string* chk; Real west, east; const char* note; };
    const Case cases[] = {
        {true,  &with,    t_deep_west,    t_deep_east,    nullptr},
        {true,  &without, t_deep_default, t_deep_default,
         "NOTE: the checkpoint has no TwoStream_TDeep: erf.radiation.seb_init_t_deep_from_file "
         "is not used on this restart"},
        {false, &with,    t_deep_default, t_deep_default,
         "NOTE: the checkpoint holds TwoStream_TDeep, which is not read"},
        {false, &without, t_deep_default, t_deep_default, nullptr},
    };
    for (const Case& c : cases) {
        const RadChoice rad = balance(c.on);
        TwoStreamRadiation radiation;
        radiation.resize(1);
        radiation.define_level(0, rad, RdoCp, L.ba2d, L.dm, L.ba, L.domain);
        testing::internal::CaptureStdout();
        radiation.read_checkpoint(0, *c.chk, period);
        const std::string report = testing::internal::GetCapturedStdout();
        EXPECT_EQ(value_at(*radiation.seb_t_deep(0), 0), c.west) << *c.chk << " " << c.on;
        EXPECT_EQ(value_at(*radiation.seb_t_deep(0), 1), c.east) << *c.chk << " " << c.on;
        // Only the I/O rank prints.
        if (!ParallelDescriptor::IOProcessor()) { continue; }
        if (c.note) {
            EXPECT_NE(report.find(c.note), std::string::npos) << report;
        } else {
            EXPECT_EQ(report.find("TwoStream_TDeep"), std::string::npos) << report;
        }
    }
    if (ParallelDescriptor::IOProcessor()) {
        for (const auto& chk : {with, without}) { amrex::FileSystem::RemoveAll(chk); }
    }
}
