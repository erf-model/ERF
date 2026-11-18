#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_Gpu.H>

#include <ERF_EBIFBuildings.H>
#include <ERF_EBIFTerrain.H>

#include <gtest/gtest.h>

#include <limits>
#include <type_traits>
#include <utility>

namespace {

using namespace amrex;
using TerrainBuildingsIF = decltype(EB2::makeUnion(
    std::declval<TerrainIF>(), std::declval<BuildingsLevelsetIF>()));

// GeometryShop dispatches on this trait. A host-only shop would dereference
// device memory in both the buildings array and the terrain/dz arrays (#5187).
static_assert(EB2::IsGPUable<BuildingsLevelsetIF>::value);
static_assert(EB2::IsGPUable<TerrainBuildingsIF>::value);
static_assert(std::is_trivially_copyable_v<TerrainBuildingsIF>);

template <typename F>
void check_geometry (F const& implicit_function, Geometry const& geom, bool terrain)
{
    auto shop = EB2::makeShop(implicit_function);
    Box const nodes = surroundingNodes(geom.Domain());
    FArrayBox levelset(nodes, 1);
    // Exercise GeometryShop's dispatch, not just a direct call to the functor.
    shop.fillFab(levelset, geom, RunOn::Gpu, nodes);
    FArrayBox host_levelset(nodes, 1, The_Pinned_Arena());
    Gpu::dtoh_memcpy(host_levelset.dataPtr(), levelset.dataPtr(),
                     levelset.size()*sizeof(Real));
    auto const values = host_levelset.const_array();

    // The fixture is [0,1]^3. These nodes are away from cube edges, where
    // trilinear interpolation reproduces the planar signed distance exactly.
    // Allow roundoff from STL distances, interpolation, and GeometryShop.
    Real const tol = Real(128)*std::numeric_limits<Real>::epsilon();
    EXPECT_NEAR(values(18,8,8), Real(0.25), tol);  // ( .25, .5, .5), solid
    EXPECT_NEAR(values(14,8,8), Real(-0.25), tol); // (-.25, .5, .5), fluid
    EXPECT_NEAR(values(14,8,5), terrain ? Real(0.125) : Real(-0.25), tol);

    // Keep every backing allocation alive through both EB construction paths.
    EB2::Build(shop, geom, 0, 0, 4, false);
    auto* index_space = &EB2::IndexSpace::top();
    EXPECT_TRUE(index_space->getLevel(geom).isOK());
    EXPECT_FALSE(index_space->getLevel(geom).isAllRegular());
    EB2::BuildFC();
    Gpu::streamSynchronize();
    EB2::IndexSpace::erase(index_space);
}

void run_buildings_test (bool include_terrain)
{
    // More than 32 cells in x makes the owner's sampling span multiple boxes
    // and, in the MPI suite, multiple ranks before replication.
    Box const domain(IntVect(0), IntVect(39,15,15));
    RealBox const bounds({Real(-2), Real(-0.5), Real(-0.5)},
                         {Real(3), Real(1.5), Real(1.5)});
    int const periodic[3] = {0,0,0};
    Geometry const geom(domain, &bounds, 0, periodic);
    BuildingsIF owner(ERF_BUILDINGS_TEST_STL, Real(1), {Real(0),Real(0),Real(0)}, 0, geom);

    if (include_terrain) {
        FArrayBox heights(makeSlab(surroundingNodes(domain), 2, 0), 1);
        heights.setVal<RunOn::Gpu>(Real(0.25));
        Gpu::DeviceVector<Real> dz(domain.length(2), Real(0.125));
        TerrainIF const terrain(heights, geom, dz);
        check_geometry(EB2::makeUnion(terrain, owner.getIF()), geom, true);
    } else {
        check_geometry(owner.getIF(), geom, false);
    }
}

TEST(BuildingsIF, BuildingsOnlyGeometryShop) { run_buildings_test(false); }
TEST(BuildingsIF, TerrainUnionGeometryShop) { run_buildings_test(true); }

} // namespace
