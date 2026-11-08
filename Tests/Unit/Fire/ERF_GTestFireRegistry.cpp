/**
 * @file ERF_GTestFireRegistry.cpp
 * @brief The fire-engine registry: an engine that registers itself with
 *        ERF_REGISTER_FIRE_ENGINE is listed and created by name
 *
 * The registry replaces a factory ERF had to be edited for. These tests
 * register an engine from this file alone, as an engine directory outside ERF
 * does, and check that erf.fire.engine finds it, that an unknown name finds
 * nothing, and that the defaults of the optional methods hold (an engine
 * without restart support stops on a restart).
 */
#include <gtest/gtest.h>

#include <ERF_FireEngineRegistry.H>
#include <AMReX_MultiFab.H>

#include <algorithm>

using namespace amrex;

namespace {

/// The required methods only: every optional one keeps its default.
class RegistryTestEngine : public FireEngine
{
public:
    explicit RegistryTestEngine (const FireCouplingParams& p) : grid_ratio(p.grid_ratio) {}
    const char* name () const override { return "registry-test"; }
    void verify_prerequisites (const ERF&, int) override {}
    void initialize (const ERF&, int, const FireGrid&, const SurfaceLayer*, const MultiFab*) override {}
    void advance (const FireForcing&) override {}
    void fill_fluxes (FireFluxes&) override {}
    Real smoke_heat_per_kg () const override { return 1.0; }
    int grid_ratio;
};

/// Default-constructed: the macro accepts an engine without a params constructor.
class RegistryTestEngineDefault : public RegistryTestEngine
{
public:
    RegistryTestEngineDefault () : RegistryTestEngine(FireCouplingParams{}) {}
    const char* name () const override { return "registry-test-default"; }
};

} // namespace

ERF_REGISTER_FIRE_ENGINE(RegistryTestEngine, "registry-test")
ERF_REGISTER_FIRE_ENGINE(RegistryTestEngineDefault, "registry-test-default")

TEST(FireRegistry, ListsSelfRegisteredEngines)
{
    const auto names = fire_engine_names();
    EXPECT_NE(std::find(names.begin(), names.end(), "registry-test"), names.end());
    EXPECT_NE(std::find(names.begin(), names.end(), "registry-test-default"), names.end());
    EXPECT_TRUE(std::is_sorted(names.begin(), names.end()));
}

TEST(FireRegistry, CreatesByNameWithTheCouplingInputs)
{
    FireCouplingParams p;
    p.grid_ratio = 7;
    auto engine = make_fire_engine("registry-test", p);
    ASSERT_NE(engine, nullptr);
    EXPECT_STREQ(engine->name(), "registry-test");
    EXPECT_EQ(dynamic_cast<RegistryTestEngine*>(engine.get())->grid_ratio, 7);

    auto dflt = make_fire_engine("registry-test-default", p);
    ASSERT_NE(dflt, nullptr);
    EXPECT_STREQ(dflt->name(), "registry-test-default");
}

TEST(FireRegistry, UnknownNameCreatesNothing)
{
    EXPECT_EQ(make_fire_engine("no-such-engine", FireCouplingParams{}), nullptr);
}

TEST(FireRegistry, OptionalMethodDefaults)
{
    auto engine = make_fire_engine("registry-test", FireCouplingParams{});
    ASSERT_NE(engine, nullptr);
    EXPECT_TRUE(engine->plot_var_names().empty());
    EXPECT_EQ(engine->plot_var("fire_phi"), nullptr);
    EXPECT_LT(engine->burned_area(), 0.0);
    EXPECT_FALSE(engine->wants_surface_precip());
    EXPECT_TRUE(engine->needs_surface_layer());
    engine->start_run(1.0, 10, 1);              // no-op
    engine->write_plotfile("unused", 0.0, 0);   // no-op
}

// amrex::Abort throws instead of ending the process while this is alive
struct AbortThrows
{
    bool saved = amrex::system::throw_exception;
    AbortThrows () { amrex::system::throw_exception = true; }
    ~AbortThrows () { amrex::system::throw_exception = saved; }
};

TEST(FireRegistry, ReadCheckpointDefaultStops)
{
    auto engine = make_fire_engine("registry-test", FireCouplingParams{});
    ASSERT_NE(engine, nullptr);
    AbortThrows guard;
    EXPECT_ANY_THROW(engine->read_checkpoint("chk00000", 0, 0));
}

TEST(FireRegistry, DuplicateNameStops)
{
    AbortThrows guard;
    EXPECT_ANY_THROW(register_fire_engine("registry-test",
                                          [](const FireCouplingParams& p) -> std::unique_ptr<FireEngine> {
                                              return std::make_unique<RegistryTestEngine>(p);
                                          }));
    EXPECT_ANY_THROW(register_fire_engine("", nullptr));
}
