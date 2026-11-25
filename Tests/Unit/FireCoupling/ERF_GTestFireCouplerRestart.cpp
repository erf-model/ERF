/**
 * @file ERF_GTestFireCouplerRestart.cpp
 * @brief FireCoupler with a stand-in engine: the start time the engine is
 * handed, and the levels a restart searches for the lagged fire fluxes
 *
 * - ERF keeps its time in double; the engine gets the start of the step in
 *   double too, so a time that single precision cannot hold reaches it
 *   unchanged (in a double-precision build Real is double and the check
 *   holds either way; the single-precision build is the one that tells).
 * - A restart searches every level the checkpoint holds for the lagged
 *   fluxes, not only the levels of this run's amr.max_level. A checkpoint
 *   with the fluxes on level 1, read by a run with max_level = 0, is refused
 *   with the mismatch instead of starting the fire with no lagged fluxes, and
 *   the engine is asked to search the same levels.
 */

#include <ERF_FireCoupler.H>
#include <ERF_FireCouplingParams.H>
#include <ERF_FireEngine.H>
#include <ERF_FireGrid.H>
#include <ERF_IndexDefines.H>
#include <ERF_EOS.H>
#include "../ERF_GTestTempDir.H"
#include "../MovingBodies/ERF_GTestThrowOnAbort.H"

#include <AMReX_FileSystem.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

#include <gtest/gtest.h>

#include <memory>
#include <string>
#include <type_traits>

using namespace amrex;

namespace {

constexpr int NX = 4, NY = 4, NZ = 8, C = 2;

/// Records what the coupler hands it; reports no flux.
class StandInEngine : public FireEngine
{
public:
    const char* name () const override { return "stand-in"; }
    void verify_prerequisites (const ERF&, int) override {}
    void initialize (const ERF&, int, const FireGrid&, const SurfaceLayer*, const MultiFab*) override {}
    void advance (const FireForcing& forcing) override { start_time = forcing.time; }
    void fill_fluxes (FireFluxes& fluxes) override { fluxes.sensible = nullptr; fluxes.latent = nullptr; }
    Real smoke_heat_per_kg () const override { return Real(1.0e7); }
    void read_checkpoint (const std::string&, int max_level, int) override { searched_to = max_level; }

    double start_time = -1.0;
    int searched_to = -1;
};

/// A periodic NX x NY x NZ level in two boxes, its k = 0 slabs and its state
struct Level
{
    BoxArray ba, ba2d;
    DistributionMapping dm;
    Geometry geom;
    MultiFab S, xvel, yvel, zcc;
};

std::unique_ptr<Level> make_level ()
{
    auto L = std::make_unique<Level>();
    const Box domain(IntVect(0, 0, 0), IntVect(NX - 1, NY - 1, NZ - 1));
    BoxList bl;
    bl.push_back(Box(IntVect(0, 0, 0), IntVect(1, NY - 1, NZ - 1)));
    bl.push_back(Box(IntVect(2, 0, 0), IntVect(NX - 1, NY - 1, NZ - 1)));
    L->ba = BoxArray(bl);
    L->dm = DistributionMapping(L->ba);
    L->geom = Geometry(domain, RealBox(0.0, 0.0, 0.0, NX * 100.0, NY * 100.0, NZ * 10.0),
                       CoordSys::cartesian, {1, 1, 0});
    L->ba2d = L->ba;
    L->ba2d.coarsen(IntVect(1, 1, NZ));
    L->S.define(L->ba, L->dm, NVAR_max, 1);
    L->S.setVal(Real(0.0));
    L->S.setVal(Real(1.15), Rho_comp, 1, 1);
    L->S.setVal(getRhoThetagivenP(Real(9.0e4)), RhoTheta_comp, 1, 1);
    L->xvel.define(convert(L->ba, IntVect(1, 0, 0)), L->dm, 1, 1); L->xvel.setVal(Real(1.0));
    L->yvel.define(convert(L->ba, IntVect(0, 1, 0)), L->dm, 1, 1); L->yvel.setVal(Real(1.0));
    L->zcc.define(L->ba, L->dm, 1, 1);
    L->zcc.setVal(Real(5.0));
    return L;
}

std::unique_ptr<FireCoupler> make_coupler (const Level& lev, StandInEngine*& engine)
{
    FireCouplingParams p;
    p.grid_ratio    = C;
    p.coupling_type = "lagged";
    auto eng = std::make_unique<StandInEngine>();
    engine = eng.get();
    auto cpl = std::make_unique<FireCoupler>(p, std::move(eng));
    cpl->set_fire_grid(0, create_fire_grid(lev.ba, lev.dm, lev.geom, C), lev.ba, lev.dm, lev.geom);
    return cpl;
}

} // namespace

// The start of the step reaches the engine in double: 1e8 + 1/8 s is held
// exactly by a double and rounds to 1e8 s in single precision
TEST(FireCoupler, TheEngineGetsTheStartTimeInDouble)
{
    static_assert(std::is_same_v<decltype(FireForcing::time), double>,
                  "FireForcing::time is the double ERF keeps its time in");
    auto lev = make_level();
    StandInEngine* eng = nullptr;
    auto cpl = make_coupler(*lev, eng);
    const double t = 1.0e8 + 0.125;
    ASSERT_NE(static_cast<double>(static_cast<float>(t)), t) << "the test time must not fit a float";
    cpl->advance(t, 0.5, lev->S, lev->xvel, lev->yvel, lev->xvel, lev->yvel, lev->zcc,
                 lev->ba2d, lev->geom, SurfacePrecipAccumulationSources{});
    EXPECT_EQ(eng->start_time, t);
}

// A checkpoint with the lagged fluxes on level 1, read by a run with
// amr.max_level = 0 and the fire on level 0, is refused with the mismatch
TEST(FireCoupler, ACheckpointWithTheFireAboveThisRunsLevelsIsRefused)
{
    auto lev = make_level();
    const std::string chk = erf_gtest_temp_path("fire_coupler_chk_finer").string();
    if (ParallelDescriptor::IOProcessor()) {
        UtilCreateDirectory(chk + "/Level_0", 0755);
        UtilCreateDirectory(chk + "/Level_1", 0755);
    }
    ParallelDescriptor::Barrier();
    // the lagged fluxes of a fire on level 1 (refinement ratio 2, the whole domain)
    BoxArray ba1 = lev->ba2d;
    ba1.refine(IntVect(2, 2, 1));
    MultiFab q1(ba1, lev->dm, 1, 0);
    q1.setVal(Real(1.0));
    VisMF::Write(q1, MultiFabFileFullPrefix(1, chk, "Level_", "FireQAtmPrev"));
    ParallelDescriptor::Barrier();

    StandInEngine* eng = nullptr;
    auto cpl = make_coupler(*lev, eng);
    const std::string msg = erf_gtest::abort_message([&] { cpl->read_checkpoint(chk, 0, 0); });
    EXPECT_NE(msg.find("holds the lagged fire fluxes on level 1"), std::string::npos) << msg;
    EXPECT_EQ(eng->searched_to, 1) << "the engine searches every level the checkpoint holds";

    // the same checkpoint without Level_1: nothing on this run's levels, the fluxes stay zero
    if (ParallelDescriptor::IOProcessor()) { FileSystem::RemoveAll(chk + "/Level_1"); }
    ParallelDescriptor::Barrier();
    StandInEngine* eng2 = nullptr;
    auto cpl2 = make_coupler(*lev, eng2);
    EXPECT_EQ(erf_gtest::abort_message([&] { cpl2->read_checkpoint(chk, 0, 0); }), "");
    EXPECT_EQ(eng2->searched_to, 0);
    EXPECT_EQ(cpl2->get_Q_atm_prev()->norm0(0, 1, cpl2->get_Q_atm_prev()->nGrowVect()), Real(0.0));

    if (ParallelDescriptor::IOProcessor()) { FileSystem::RemoveAll(chk); }
    ParallelDescriptor::Barrier();
}
