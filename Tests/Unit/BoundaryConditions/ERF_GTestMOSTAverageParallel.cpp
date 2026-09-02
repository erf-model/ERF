#include "ERF_GTestMOSTAverageParallelCommon.H"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <string>

using namespace amrex;
using namespace erf_most_average_test;

namespace {

std::array<Orientation, 6> all_faces ()
{
    return {{
        Orientation(Direction::x, Orientation::low),
        Orientation(Direction::x, Orientation::high),
        Orientation(Direction::y, Orientation::low),
        Orientation(Direction::y, Orientation::high),
        Orientation(Direction::z, Orientation::low),
        Orientation(Direction::z, Orientation::high)
    }};
}

void expect_reference_metadata (const MOSTAverage& averages,
                                const Geometry& geom,
                                const Orientation face)
{
    const int dir = face.coordDir();
    const Box& domain = geom.Domain();
    const int expected_index = face.isLow()
        ? domain.smallEnd(dir) : domain.bigEnd(dir);
    const Real expected_zref = (dir == 2 && face.isLow())
        ? geom.ProbLo(dir) + Real(0.5) * geom.CellSize(dir)
        : Real(0.5) * geom.CellSize(dir);
    EXPECT_EQ(averages.get_k_indices(0)->min(0), expected_index);
    EXPECT_EQ(averages.get_k_indices(0)->max(0), expected_index);
    EXPECT_NEAR(averages.get_zref(0)->min(0), expected_zref,
                tolerance(expected_zref));
    EXPECT_NEAR(averages.get_zref(0)->max(0), expected_zref,
                tolerance(expected_zref));
}

void expect_all_plane_values (const Vector<Real>& actual,
                              const Orientation face,
                              const bool distributed)
{
    const auto expected = expected_average_values();
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t comp = 0; comp < 7; ++comp) {
        EXPECT_NEAR(actual[comp], expected[comp], tolerance(expected[comp]))
            << "direction=" << face.coordDir()
            << ", high=" << !face.isLow()
            << ", distributed=" << distributed
            << ", component=" << comp;
    }
    if (face.coordDir() < 2) {
        EXPECT_NEAR(actual[7], expected[7], tolerance(expected[7]))
            << "direction=" << face.coordDir()
            << ", high=" << !face.isLow()
            << ", distributed=" << distributed << ", component=7";
        EXPECT_NEAR(actual[8], expected[8], tolerance(expected[8]))
            << "direction=" << face.coordDir()
            << ", high=" << !face.isLow()
            << ", distributed=" << distributed << ", component=8";
    }
}

// Motivation: plane averages must be independent of MPI ownership and the
// lateral tile boundaries used by TileNoZ(). Metadata is part of this oracle
// because it selects the source plane used by every average component.
TEST(MOSTAverageParallel, DistributedAndTiledPlaneAverageMatchesReference)
{
    ScopedMFIterTileSize tile_size(IntVect(AMREX_D_DECL(1024, 4, 1024)));

    for (const bool distributed : {false, true}) {
        for (const auto& face : all_faces()) {
            MOSTAverageFields fields(distributed);
            const std::string prefix =
                "unit_most_parallel_plane_" + std::to_string(distributed) + "_" +
                std::to_string(face.coordDir()) + (face.isLow() ? "_lo" : "_hi");
            ScopedMOSTParams params(prefix.c_str());
            MOSTAverage averages(face, Vector<Geometry>{fields.geom}, false,
                                  prefix, MeshType::ConstantDz, TerrainType::None);
            averages.make_MOSTAverage_at_level(
                0, fields.vars_old, fields.theta, fields.qv,
                fields.qr, fields.z_phys_nd);
            expect_reference_metadata(averages, fields.geom, face);
            averages.compute_averages(0);
            Gpu::streamSynchronize();
            expect_all_plane_values(averages.get_plane_average(0), face, distributed);
        }
    }
}

// Motivation: both averaging policies must honor explicit zref and k_arr_in
// on every distributed low/high Cartesian face.
TEST(MOSTAverageParallel, BothPoliciesSupportKAndZrefOnEveryWall)
{
    const Real requested_distance = Real(4.0);
    for (const int policy : {0, 1}) {
        for (const bool use_zref : {true, false}) {
            for (const auto& face : all_faces()) {
                MOSTAverageFields fields(true);
                const int dir = face.coordDir();
                const Real dz = fields.geom.CellSize(dir);
                const bool zlo_absolute = dir == 2 && face.isLow();
                const Real requested_zref = zlo_absolute
                    ? fields.geom.ProbLo(dir) + requested_distance
                    : requested_distance;
                const int requested_k = fields.domain.smallEnd(dir) + 1;
                const std::string prefix =
                    "unit_most_parallel_policy_" + std::to_string(policy) + "_" +
                    (use_zref ? "zref_" : "k_") + std::to_string(dir) +
                    (face.isLow() ? "_lo" : "_hi");
                ScopedMOSTParams params(prefix.c_str(), policy, 0);
                if (use_zref) { params.add_zref(requested_zref); }
                else { params.add_k_index(requested_k); }

                MOSTAverage averages(face, Vector<Geometry>{fields.geom}, false,
                                      prefix, MeshType::ConstantDz,
                                      TerrainType::None);
                averages.make_MOSTAverage_at_level(
                    0, fields.vars_old, fields.theta, fields.qv,
                    fields.qr, fields.z_phys_nd);

                int wall_offset = 0;
                int expected_index = requested_k;
                if (use_zref) {
                    wall_offset = static_cast<int>(std::floor(
                        (zlo_absolute ? requested_zref - fields.geom.ProbLo(dir)
                                      : requested_zref) / dz - myhalf));
                    expected_index = face.isLow()
                        ? fields.domain.smallEnd(dir) + wall_offset
                        : fields.domain.bigEnd(dir) - wall_offset;
                } else {
                    wall_offset = face.isLow()
                        ? requested_k - fields.domain.smallEnd(dir)
                        : fields.domain.bigEnd(dir) - requested_k;
                }
                const Real expected_wall_distance =
                    (static_cast<Real>(wall_offset) + myhalf) * dz;
                const Real expected_zref = zlo_absolute
                    ? fields.geom.ProbLo(dir) + expected_wall_distance
                    : expected_wall_distance;
                EXPECT_EQ(averages.get_k_indices(0)->min(0), expected_index);
                EXPECT_EQ(averages.get_k_indices(0)->max(0), expected_index);
                EXPECT_NEAR(averages.get_zref(0)->min(0), expected_zref,
                            tolerance(expected_zref));
                EXPECT_NEAR(averages.get_zref(0)->max(0), expected_zref,
                            tolerance(expected_zref));

                averages.compute_averages(0);
                Gpu::streamSynchronize();
                expect_all_plane_values(averages.get_plane_average(0), face, true);
            }
        }
    }
}

// Motivation: fitted terrain selects its reference index from nodal heights,
// and both averaging policies must consume that index. The distributed box
// layout deliberately leaves the vertical extent unsplit.
TEST(MOSTAverageParallel, DistributedTerrainKIndexDrivesBothPolicies)
{
    const Real requested_distance = Real(1.5);
    const Geometry geom = make_parallel_geometry();
    for (const auto& face : all_faces()) {
        const int dir = face.coordDir();
        const int wall_offset = static_cast<int>(std::floor(
            requested_distance / geom.CellSize(dir) - myhalf));
        const int expected_index = face.isLow()
            ? geom.Domain().smallEnd(dir) + wall_offset
            : geom.Domain().bigEnd(dir) - wall_offset;

        for (const int policy : {0, 1}) {
            MOSTAverageFields fields(geom, true);
            for (MFIter mfi(*fields.theta, false); mfi.isValid(); ++mfi) {
                auto theta_arr = fields.theta->array(mfi);
                ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const int normal_index = dir == 0 ? i : (dir == 1 ? j : k);
                    theta_arr(i,j,k) = Real(100.0) + Real(10.0) * normal_index;
                });
            }
            BoxArray ba_nd(fields.ba);
            ba_nd.surroundingNodes();
            auto z_phys_nd = std::make_unique<MultiFab>(ba_nd, fields.dm, 1, 1);
            initialize_uniform_terrain_height(*z_phys_nd, fields.geom);

            const std::string prefix = "unit_most_parallel_terrain_" +
                std::to_string(dir) + (face.isLow() ? "_lo_" : "_hi_") +
                std::to_string(policy);
            ScopedMOSTParams params(prefix.c_str(), policy, 0);
            params.add_zref(requested_distance);
            MOSTAverage averages(face, Vector<Geometry>{fields.geom}, false,
                                  prefix, MeshType::ConstantDz,
                                  TerrainType::StaticFittedMesh);
            averages.make_MOSTAverage_at_level(
                0, fields.vars_old, fields.theta, fields.qv,
                fields.qr, z_phys_nd);
            Gpu::streamSynchronize();
            EXPECT_EQ(averages.get_k_indices(0)->min(0), expected_index);
            EXPECT_EQ(averages.get_k_indices(0)->max(0), expected_index);

            averages.compute_averages(0);
            Gpu::streamSynchronize();
            const Real expected_theta = Real(100.0) + Real(10.0) * expected_index;
            if (policy == 0) {
                ASSERT_GE(averages.get_plane_average(0).size(), std::size_t(4));
                EXPECT_NEAR(averages.get_plane_average(0)[3], expected_theta,
                            tolerance(expected_theta));
            } else {
                EXPECT_NEAR(averages.get_average(0, 3)->min(0), expected_theta,
                            tolerance(expected_theta));
            }
        }
    }
}

// Motivation: region averages must write only selected boundary boxes. The
// untouched-box oracle also catches a write that accidentally uses a rank's
// local FAB rather than the selected global face.
TEST(MOSTAverageParallel, DistributedAndTiledRegionAverageMatchesReference)
{
    ScopedMFIterTileSize tile_size(IntVect(AMREX_D_DECL(1024, 4, 1024)));
    const auto expected = expected_average_values();

    for (const bool distributed : {false, true}) {
        for (const auto& face : all_faces()) {
            MOSTAverageFields fields(distributed);
            const std::string prefix = "unit_most_parallel_region_" +
                std::to_string(distributed) + "_" + std::to_string(face.coordDir()) +
                (face.isLow() ? "_lo" : "_hi");
            ScopedMOSTParams params(prefix.c_str(), 1, 0);
            MOSTAverage averages(face, Vector<Geometry>{fields.geom}, false,
                                  prefix, MeshType::ConstantDz, TerrainType::None);
            averages.make_MOSTAverage_at_level(
                0, fields.vars_old, fields.theta, fields.qv,
                fields.qr, fields.z_phys_nd);
            averages.compute_averages(0);
            Gpu::streamSynchronize();

            const std::array<const BoxArray*, 9> source_ba{{
                &fields.xvel.boxArray(), &fields.yvel.boxArray(), &fields.zvel.boxArray(),
                &fields.ba, &fields.ba, &fields.ba, &fields.ba, &fields.ba, &fields.ba
            }};
            const std::size_t nchecked = face.coordDir() < 2 ? expected.size() : 7;
            for (std::size_t comp = 0; comp < nchecked; ++comp) {
                const auto* average = averages.get_average(0, static_cast<int>(comp));
                const Real initial = comp == 2
                    ? Real(1.e34)
                    : (comp == 4 ? Real(0.0) : bogus_large_value);
                const auto range = region_range(
                    *average, *source_ba[comp], fields.domain, face, initial);
                EXPECT_GT(range.selected_count, 0);
                EXPECT_NEAR(range.selected_lo, expected[comp], tolerance(expected[comp]));
                EXPECT_NEAR(range.selected_hi, expected[comp], tolerance(expected[comp]));
                if (distributed) {
                    EXPECT_GT(range.nonselected_count, 0);
                } else {
                    EXPECT_EQ(range.nonselected_count, 0);
                }
                EXPECT_EQ(range.nonselected_bad, 0)
                    << "direction=" << face.coordDir()
                    << ", high=" << !face.isLow()
                    << ", distributed=" << distributed
                    << ", component=" << comp;
            }
        }
    }
}

} // namespace
