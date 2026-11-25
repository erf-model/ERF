// Unit tests for the surface-coupling layer the fire and the flood couplings
// share (Source/SurfaceCoupling), without an ERF object.
//
// The atmospheric level 0 is 12 x 9 x 3 cells in four boxes; the surface grid
// refines it by C = 3. The tests check that the coarsening gives the bits of
// amrex::average_down() on level 0 (the fire coupling used average_down()
// there before it moved to this layer), and that the checks on the inputs of
// the maps abort with the coupling's own names.

#include <ERF_SurfaceCouplingNames.H>
#include <ERF_SurfaceEngineRegistry.H>
#include <ERF_SurfaceGrid.H>
#include "../MovingBodies/ERF_GTestThrowOnAbort.H"

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <string>

using namespace amrex;

namespace {

constexpr int nx = 12, ny = 9, nz = 3;
constexpr int C  = 3;
// With C = 1 every summation order gives the same bits and the first test
// would prove nothing
static_assert(C >= 2, "the coarsening test needs C >= 2");

const SurfaceCouplingNames& test_names ()
{
    static const SurfaceCouplingNames names{"fire", "Fire", "FIRE"};
    return names;
}

Geometry level_geom ()
{
    const Box dom(IntVect(0, 0, 0), IntVect(nx - 1, ny - 1, nz - 1));
    const RealBox rb({Real(0.0), Real(0.0), Real(0.0)},
                     {Real(30.0) * nx, Real(30.0) * ny, Real(20.0) * nz});
    return Geometry(dom, rb, CoordSys::cartesian, {1, 0, 0});
}

// Four boxes of unequal size, each spanning the height
BoxArray level_boxes ()
{
    BoxList bl;
    bl.push_back(Box(IntVect(0, 0, 0), IntVect(4, 3, nz - 1)));
    bl.push_back(Box(IntVect(5, 0, 0), IntVect(11, 3, nz - 1)));
    bl.push_back(Box(IntVect(0, 4, 0), IntVect(6, 8, nz - 1)));
    bl.push_back(Box(IntVect(7, 4, 0), IntVect(11, 8, nz - 1)));
    return BoxArray(std::move(bl));
}

// The k = 0 slabs of the level's boxes, in their order
BoxArray k0_slabs (const BoxArray& ba)
{
    BoxList bl;
    for (int n = 0; n < ba.size(); ++n) {
        Box b = ba[n];
        b.setSmall(2, 0);
        b.setBig(2, 0);
        bl.push_back(b);
    }
    return BoxArray(std::move(bl));
}

// Values whose sums round differently in a different order
void set_irregular (MultiFab& mf)
{
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const auto a = mf.array(mfi);
        ParallelFor(mfi.validbox(), mf.nComp(), [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            a(i, j, k, n) = Real(1.0) / Real(3.0)
                          + Real(1.0e3) * std::sin(Real(1.3) * Real(i) + Real(0.7) * Real(j) + Real(n))
                          + Real(1.0e-3) * Real(i * j + k);
        });
    }
}

struct DummyEngine { int value = 0; };

} // namespace

// On level 0 the region's corner is the origin, and surface_coarsen_to_atm()
// sums in the order of amrex::average_down() with ratio (C, C, 1): the fire
// coupling's level-0 answers stay the same to the last bit.
TEST(SurfaceGrid, CoarsenGivesTheBitsOfAverageDownOnLevelZero)
{
    const Geometry geom = level_geom();
    const BoxArray ba   = level_boxes();
    const DistributionMapping dm(ba);
    const SurfaceGrid sg = create_surface_grid(ba, dm, geom, C, test_names());
    ASSERT_EQ(sg.atm_lo, IntVect(0, 0, 0));

    MultiFab surf(sg.ba, sg.dm, 2, 0);
    set_irregular(surf);

    const BoxArray slabs = k0_slabs(ba);
    MultiFab mine(slabs, dm, 2, 0);
    MultiFab ref(slabs, dm, 2, 0);
    mine.setVal(Real(-1.0));
    ref.setVal(Real(-2.0));
    surface_coarsen_to_atm(mine, surf, sg, test_names());
    amrex::average_down(surf, ref, 0, 2, IntVect(C, C, 1));

    // Device-side count of the cells that differ, summed over the ranks, so
    // the test also runs on a GPU build
    const auto m = mine.const_arrays();
    const auto r = ref.const_arrays();
    const int ndiff = ParReduce(TypeList<ReduceOpSum>{}, TypeList<int>{}, mine, IntVect(0), 2,
        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n) noexcept -> GpuTuple<int>
        {
            return { m[b](i, j, k, n) != r[b](i, j, k, n) ? 1 : 0 };
        });
    int ndiff_all = ndiff;
    ParallelDescriptor::ReduceIntSum(ndiff_all);
    EXPECT_EQ(ndiff_all, 0);
}

// Every check on the inputs of the maps names the coupling and its grid
TEST(SurfaceGrid, ChecksNameTheCoupling)
{
    const Geometry geom = level_geom();
    const BoxArray ba   = level_boxes();
    const DistributionMapping dm(ba);

    std::string msg = erf_gtest::abort_message([&] {
        create_surface_grid(ba, dm, geom, 0, test_names());
    });
    EXPECT_NE(msg.find("[FIRE] Building the fire grid: erf.fire.grid_ratio must be >= 1"), std::string::npos) << msg;

    const SurfaceGrid sg = create_surface_grid(ba, dm, geom, C, test_names());
    MultiFab surf(sg.ba, sg.dm, 1, 0);
    surf.setVal(Real(1.0));

    // Same number of boxes and the same owners, but box n of the slabs is not
    // the box surface box n refines
    BoxList swapped;
    const BoxArray slabs = k0_slabs(ba);
    swapped.push_back(slabs[1]);
    swapped.push_back(slabs[0]);
    swapped.push_back(slabs[2]);
    swapped.push_back(slabs[3]);
    MultiFab atm_bad(BoxArray(std::move(swapped)), dm, 1, 0);
    atm_bad.setVal(Real(0.0));

    msg = erf_gtest::abort_message([&] { surface_coarsen_to_atm(atm_bad, surf, sg, test_names()); });
    EXPECT_NE(msg.find("[FIRE] Averaging the fire grid onto the atmospheric columns: the two must pair box for box"),
              std::string::npos) << msg;

    msg = erf_gtest::abort_message([&] { surface_fill_from_atm_k0(surf, atm_bad, sg, test_names()); });
    EXPECT_NE(msg.find("[FIRE] Copying the atmospheric columns to the fire grid: the two must pair box for box"),
              std::string::npos) << msg;

    // The paired slabs pass both maps
    MultiFab atm_ok(slabs, dm, 1, 0);
    atm_ok.setVal(Real(0.0));
    EXPECT_TRUE(surface_pairs_with_atm(surf, atm_ok, sg));
    msg = erf_gtest::abort_message([&] {
        surface_coarsen_to_atm(atm_ok, surf, sg, test_names());
        surface_fill_from_atm_k0(surf, atm_ok, sg, test_names());
    });
    EXPECT_EQ(msg, "");

    MultiFab nodal(convert(sg.ba, IntVect(1, 0, 0)), sg.dm, 1, 1);
    EXPECT_FALSE(surface_pairs_with_atm(nodal, atm_ok, sg));
    msg = erf_gtest::abort_message([&] { surface_fill_boundary(nodal, sg.geom, test_names()); });
    EXPECT_NE(msg.find("[FIRE] Filling the ghost cells of the fire grid: the MultiFab must be cell-centred"),
              std::string::npos) << msg;
}

// The registry behind the registration macros: sorted names, creation, and the
// two registration errors with the coupling's macro name. The registry is
// process-wide, so the two engines are added only on the first run of the test
// (gtest_repeat runs it again in the same process).
TEST(SurfaceEngineRegistry, NamesCreationAndErrors)
{
    using Registry = SurfaceEngineRegistry<DummyEngine, int>;
    const SurfaceCouplingNames& names = test_names();
    if (Registry::names().empty()) {
        Registry::add(names, "zeta", [](const int& v) { auto e = std::make_unique<DummyEngine>(); e->value = v; return e; });
        Registry::add(names, "alpha", [](const int&) { return std::make_unique<DummyEngine>(); });
    }

    const auto list = Registry::names();
    ASSERT_EQ(list.size(), 2u);
    EXPECT_EQ(list[0], "alpha");
    EXPECT_EQ(list[1], "zeta");
    auto e = Registry::make("zeta", 7);
    ASSERT_NE(e, nullptr);
    EXPECT_EQ(e->value, 7);
    EXPECT_EQ(Registry::make("none", 0), nullptr);

    std::string msg = erf_gtest::abort_message([&] {
        Registry::add(names, "zeta", [](const int&) { return std::make_unique<DummyEngine>(); });
    });
    EXPECT_NE(msg.find("[FIRE] ERF_REGISTER_FIRE_ENGINE: two engines registered the name \"zeta\""),
              std::string::npos) << msg;
    msg = erf_gtest::abort_message([&] {
        Registry::add(names, "", [](const int&) { return std::make_unique<DummyEngine>(); });
    });
    EXPECT_NE(msg.find("[FIRE] ERF_REGISTER_FIRE_ENGINE: an engine registered an empty name"),
              std::string::npos) << msg;
}
