#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

#include <ERF_Constants.H>
#include <ERF_MOSTAverage.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cmath>
#include <limits>
#include <memory>
#include <string>

using namespace amrex;

namespace {

class ScopedMOSTParams
{
public:
    explicit ScopedMOSTParams (const char* prefix,
                               const int policy = 0,
                               const int radius = 0)
        : m_pp(prefix)
    {
        m_pp.add("most.average_policy", policy);
        m_pp.add("most.radius", radius);
    }

    void add_zref (const Real zref)
    {
        m_pp.add("most.zref", zref);
    }

    void add_k_index (const int k_index)
    {
        m_pp.addarr("most.k_arr_in", Vector<int>{k_index});
    }

    ~ScopedMOSTParams ()
    {
        for (const auto* name : m_names) {
            m_pp.remove(name);
        }
    }

private:
    ParmParse m_pp;
    const std::array<const char*, 10> m_names{{
        "most.radius",
        "most.time_average",
        "most.terrain_rotate",
        "most.use_interpolation",
        "most.use_normal_vector",
        "most.time_window",
        "most.include_subgrid_vel",
        "most.average_policy",
        "most.zref",
        "most.k_arr_in"
    }};
};

Geometry
make_geometry ()
{
    const Box domain(IntVect(3, 5, 7), IntVect(5, 8, 11));
    const RealBox real_box({AMREX_D_DECL(10.0, 20.0, 30.0)},
                           {AMREX_D_DECL(16.0, 32.0, 50.0)});
    const Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
    return Geometry(domain, &real_box, 0, is_periodic.data());
}

Geometry
make_terrain_geometry ()
{
    const Box domain(IntVect(0), IntVect(2, 3, 4));
    const RealBox real_box({AMREX_D_DECL(10.0, 20.0, 30.0)},
                           {AMREX_D_DECL(13.0, 24.0, 35.0)});
    const Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
    return Geometry(domain, &real_box, 0, is_periodic.data());
}

void
initialize_uniform_terrain_height (MultiFab& z_phys_nd, const Geometry& geom)
{
    const Real zlo = geom.ProbLo(2);
    const Real dz = geom.CellSize(2);
    for (MFIter mfi(z_phys_nd, false); mfi.isValid(); ++mfi) {
        Box zbox = z_phys_nd.boxArray()[mfi.index()];
        zbox.grow(z_phys_nd.nGrowVect());
        auto z_arr = z_phys_nd.array(mfi);
        ParallelFor(zbox, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            z_arr(i,j,k) = zlo + static_cast<Real>(k) * dz;
        });
    }
    Gpu::streamSynchronize();
}

std::array<Orientation, 6>
all_faces ()
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

Real
tolerance (const Real expected)
{
    return Real(128.0) * std::numeric_limits<Real>::epsilon() *
           std::max(Real(1.0), std::abs(expected));
}

const Real most_test_u = Real(3.0);
const Real most_test_v = Real(4.0);
const Real most_test_w = Real(12.0);
const Real most_test_theta = Real(300.0);
const Real most_test_qv = Real(0.01);

std::array<Real, 7>
expected_average_values ()
{
    return {{
        most_test_u,
        most_test_v,
        most_test_w,
        most_test_theta,
        most_test_qv,
        most_test_theta * (one + epsv * most_test_qv),
        std::sqrt(most_test_u*most_test_u + most_test_v*most_test_v)
    }};
}

struct MOSTAverageFields
{
    Geometry geom;
    Box domain;
    BoxArray ba;
    DistributionMapping dm;
    MultiFab cons;
    MultiFab xvel;
    MultiFab yvel;
    MultiFab zvel;
    std::unique_ptr<MultiFab> theta;
    std::unique_ptr<MultiFab> qv;
    std::unique_ptr<MultiFab> qr;
    std::unique_ptr<MultiFab> z_phys_nd;
    Vector<MultiFab*> vars_old;

    explicit MOSTAverageFields (const Geometry& geometry = make_geometry())
        : geom(geometry),
          domain(geom.Domain()),
          ba(domain),
          dm(ba),
          cons(ba, dm, 1, 1),
          xvel(BoxArray(surroundingNodes(domain, 0)), dm, 1, 1),
          yvel(BoxArray(surroundingNodes(domain, 1)), dm, 1, 1),
          zvel(BoxArray(surroundingNodes(domain, 2)), dm, 1, 1),
          theta(std::make_unique<MultiFab>(ba, dm, 1, 1)),
          qv(std::make_unique<MultiFab>(ba, dm, 1, 1))
    {
        cons.setVal(Real(1.0));
        xvel.setVal(most_test_u);
        yvel.setVal(most_test_v);
        zvel.setVal(most_test_w);
        theta->setVal(most_test_theta);
        qv->setVal(most_test_qv);

        vars_old = {&cons, &xvel, &yvel, &zvel};
    }
};

} // namespace

// Motivation: Plane MOST averages must use the interior high cell for the
// default reference. The z-low face keeps its absolute-coordinate reference
// convention, while all other faces use wall-relative distance. Tangential and
// cell-centered fields must still select the high cell plane, while only a
// normal nodal field extends to the exterior high node.
TEST(MOSTAverage, PlanarDefaultReferenceWorksOnEveryWall)
{
    MOSTAverageFields fields;
    const auto& geom = fields.geom;
    const auto& domain = fields.domain;

    ScopedMOSTParams params("unit_most_planar_default_reference");
    const auto faces = all_faces();

    for (const auto& face : faces) {
        MOSTAverage averages(face, Vector<Geometry>{geom}, false,
                            "unit_most_planar_default_reference",
                            MeshType::ConstantDz, TerrainType::None);
        averages.make_MOSTAverage_at_level(
            0, fields.vars_old, fields.theta, fields.qv,
            fields.qr, fields.z_phys_nd);

        const int dir = face.coordDir();
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

        averages.compute_averages(0);
        Gpu::streamSynchronize();
        const auto plane_average = averages.get_plane_average(0);
        ASSERT_EQ(plane_average.size(), static_cast<std::size_t>(averages.get_navg()));

        const auto expected = expected_average_values();
        for (std::size_t comp = 0; comp < expected.size(); ++comp) {
            EXPECT_NEAR(plane_average[comp], expected[comp], tolerance(expected[comp]))
                << "direction=" << dir << ", high=" << !face.isLow()
                << ", component=" << comp;
        }

        if (dir < 2) {
            const Real expected_xz =
                std::sqrt(most_test_u*most_test_u + most_test_w*most_test_w);
            const Real expected_yz =
                std::sqrt(most_test_v*most_test_v + most_test_w*most_test_w);
            EXPECT_NEAR(plane_average[7], expected_xz, tolerance(expected_xz));
            EXPECT_NEAR(plane_average[8], expected_yz, tolerance(expected_yz));
        }
    }
}

// Motivation: both averaging policies must use the same face-aware reference
// setup. Exercise the explicit zref and k_arr_in branches for every Cartesian
// low/high face.
TEST(MOSTAverage, BothPoliciesSupportKAndZrefOnEveryWall)
{
    const Real requested_distance = Real(4.0);
    const auto faces = all_faces();

    for (const int policy : {0, 1}) {
        for (const bool use_zref : {true, false}) {
            for (std::size_t face_index = 0; face_index < faces.size(); ++face_index) {
                MOSTAverageFields fields;
                const auto& geom = fields.geom;
                const auto& domain = fields.domain;
                const auto& face = faces[face_index];
                const int dir = face.coordDir();
                const Real dz = geom.CellSize(dir);
                const bool zlo_absolute = (dir == 2 && face.isLow());
                const Real requested_zref = zlo_absolute
                    ? geom.ProbLo(dir) + requested_distance
                    : requested_distance;
                const int requested_k = domain.smallEnd(dir) + 1;
                const std::string prefix =
                    "unit_most_policy_" + std::to_string(policy) + "_" +
                    (use_zref ? "zref_" : "k_") +
                    std::to_string(face_index);
                ScopedMOSTParams params(prefix.c_str(), policy, 0);
                if (use_zref) {
                    params.add_zref(requested_zref);
                } else {
                    params.add_k_index(requested_k);
                }

                MOSTAverage averages(face, Vector<Geometry>{geom}, false,
                                      prefix, MeshType::ConstantDz,
                                      TerrainType::None);
                averages.make_MOSTAverage_at_level(
                    0, fields.vars_old, fields.theta, fields.qv,
                    fields.qr, fields.z_phys_nd);

                int wall_offset;
                int expected_index;
                if (use_zref) {
                    wall_offset = static_cast<int>(
                        std::floor((zlo_absolute
                            ? requested_zref - geom.ProbLo(dir)
                            : requested_zref) / dz - myhalf));
                    expected_index = face.isLow()
                        ? domain.smallEnd(dir) + wall_offset
                        : domain.bigEnd(dir) - wall_offset;
                } else {
                    expected_index = requested_k;
                    wall_offset = face.isLow()
                        ? requested_k - domain.smallEnd(dir)
                        : domain.bigEnd(dir) - requested_k;
                }
                const Real expected_wall_distance =
                    (static_cast<Real>(wall_offset) + myhalf) * dz;
                const Real expected_zref = zlo_absolute
                    ? geom.ProbLo(dir) + expected_wall_distance
                    : expected_wall_distance;
                EXPECT_EQ(averages.get_k_indices(0)->min(0), expected_index)
                    << "policy=" << policy << ", zref=" << use_zref
                    << ", direction=" << dir << ", high=" << !face.isLow();
                EXPECT_EQ(averages.get_k_indices(0)->max(0), expected_index)
                    << "policy=" << policy << ", zref=" << use_zref
                    << ", direction=" << dir << ", high=" << !face.isLow();
                EXPECT_NEAR(averages.get_zref(0)->min(0), expected_zref,
                            tolerance(expected_zref))
                    << "policy=" << policy << ", zref=" << use_zref
                    << ", direction=" << dir << ", high=" << !face.isLow();
                EXPECT_NEAR(averages.get_zref(0)->max(0), expected_zref,
                            tolerance(expected_zref))
                    << "policy=" << policy << ", zref=" << use_zref
                    << ", direction=" << dir << ", high=" << !face.isLow();

                averages.compute_averages(0);
                Gpu::streamSynchronize();
                const auto expected = expected_average_values();
                for (std::size_t comp = 0; comp < expected.size(); ++comp) {
                    const Real actual =
                        averages.get_average(0, static_cast<int>(comp))->min(0);
                    EXPECT_NEAR(actual, expected[comp], tolerance(expected[comp]))
                        << "policy=" << policy << ", zref=" << use_zref
                        << ", direction=" << dir << ", high=" << !face.isLow()
                        << ", component=" << comp;
                }

                if (dir < 2) {
                    const Real expected_xz =
                        std::sqrt(most_test_u*most_test_u + most_test_w*most_test_w);
                    const Real expected_yz =
                        std::sqrt(most_test_v*most_test_v + most_test_w*most_test_w);
                    EXPECT_NEAR(averages.get_average(0, 7)->min(0),
                                expected_xz, tolerance(expected_xz))
                        << "policy=" << policy << ", zref=" << use_zref
                        << ", direction=" << dir << ", high=" << !face.isLow();
                    EXPECT_NEAR(averages.get_average(0, 8)->min(0),
                                expected_yz, tolerance(expected_yz))
                        << "policy=" << policy << ", zref=" << use_zref
                        << ", direction=" << dir << ", high=" << !face.isLow();
                }
            }
        }
    }
}

// Motivation: set_k_indices_N is the fallback path for non-terrain meshes and
// must apply the same low/high wall conventions on every Cartesian face.
TEST(MOSTAverage, NonTerrainKIndicesAreSetOnEveryWall)
{
    MOSTAverageFields fields;
    const auto& geom = fields.geom;
    const auto& domain = fields.domain;
    const Real requested_distance = Real(4.0);
    const auto faces = all_faces();

    for (std::size_t face_index = 0; face_index < faces.size(); ++face_index) {
        const auto& face = faces[face_index];
        const int dir = face.coordDir();
        const Real dz = geom.CellSize(dir);
        const bool zlo_absolute = (dir == 2 && face.isLow());
        const Real requested_zref = zlo_absolute
            ? geom.ProbLo(dir) + requested_distance
            : requested_distance;
        const std::string prefix =
            "unit_most_nonterrain_k_indices_" + std::to_string(face_index);
        ScopedMOSTParams params(prefix.c_str());
        params.add_zref(requested_zref);

        MOSTAverage averages(face, Vector<Geometry>{geom}, false, prefix,
                             MeshType::ConstantDz, TerrainType::None);
        averages.make_MOSTAverage_at_level(
            0, fields.vars_old, fields.theta, fields.qv,
            fields.qr, fields.z_phys_nd);
        Gpu::streamSynchronize();

        const int wall_offset = static_cast<int>(
            std::floor(requested_distance / dz - myhalf));
        const int expected_index = face.isLow()
            ? domain.smallEnd(dir) + wall_offset
            : domain.bigEnd(dir) - wall_offset;
        const Real expected_wall_distance =
            (static_cast<Real>(wall_offset) + myhalf) * dz;
        const Real expected_zref = zlo_absolute
            ? geom.ProbLo(dir) + expected_wall_distance
            : expected_wall_distance;

        const auto* k_indices = averages.get_k_indices(0);
        EXPECT_EQ(k_indices->min(0), expected_index)
            << "direction=" << dir << ", high=" << !face.isLow();
        EXPECT_EQ(k_indices->max(0), expected_index)
            << "direction=" << dir << ", high=" << !face.isLow();

        const auto* zref = averages.get_zref(0);
        EXPECT_NEAR(zref->min(0), expected_zref, tolerance(expected_zref))
            << "direction=" << dir << ", high=" << !face.isLow();
        EXPECT_NEAR(zref->max(0), expected_zref, tolerance(expected_zref))
            << "direction=" << dir << ", high=" << !face.isLow();
    }
}

// Motivation: set_k_indices_T must select the cell containing the requested
// wall-relative reference point on every Cartesian face.  Lateral terrain
// faces use the wall-normal path, while the z faces search the nodal terrain
// heights directly.
TEST(MOSTAverage, TerrainKIndicesAreSetOnEveryWall)
{
    const Geometry geom = make_terrain_geometry();
    const Box domain = geom.Domain();
    const Real requested_distance = Real(1.5);
    const auto faces = all_faces();

    for (std::size_t face_index = 0; face_index < faces.size(); ++face_index) {
        MOSTAverageFields fields(geom);
        const auto& face = faces[face_index];
        const std::string prefix =
            "unit_most_terrain_k_indices_" + std::to_string(face_index);
        ScopedMOSTParams params(prefix.c_str(), 0, 0);
        params.add_zref(requested_distance);

        BoxArray ba_nd(fields.ba);
        ba_nd.surroundingNodes();
        auto z_phys_nd = std::make_unique<MultiFab>(ba_nd, fields.dm, 1, 1);
        initialize_uniform_terrain_height(*z_phys_nd, geom);

        MOSTAverage averages(face, Vector<Geometry>{geom}, false, prefix,
                              MeshType::ConstantDz,
                              TerrainType::StaticFittedMesh);
        averages.make_MOSTAverage_at_level(
            0, fields.vars_old, fields.theta, fields.qv,
            fields.qr, z_phys_nd);
        Gpu::streamSynchronize();

        const int dir = face.coordDir();
        const int wall_offset = static_cast<int>(
            std::floor(requested_distance / geom.CellSize(dir) - myhalf));
        const int expected_index = face.isLow()
            ? domain.smallEnd(dir) + wall_offset
            : domain.bigEnd(dir) - wall_offset;
        const auto* k_indices = averages.get_k_indices(0);
        EXPECT_EQ(k_indices->min(0), expected_index)
            << "direction=" << dir << ", high=" << !face.isLow();
        EXPECT_EQ(k_indices->max(0), expected_index)
            << "direction=" << dir << ", high=" << !face.isLow();

        const auto* zref = averages.get_zref(0);
        EXPECT_NEAR(zref->min(0), requested_distance,
                    tolerance(requested_distance))
            << "direction=" << dir << ", high=" << !face.isLow();
        EXPECT_NEAR(zref->max(0), requested_distance,
                    tolerance(requested_distance))
            << "direction=" << dir << ", high=" << !face.isLow();
    }
}

// Motivation: terrain k-indices are consumed by both averaging policies on
// every face.  Vary the scalar along the face normal so sampling a neighboring
// cell is observable for x, y, and z low/high faces.
TEST(MOSTAverage, TerrainKIndexDrivesBothAveragingPoliciesOnEveryWall)
{
    const Geometry geom = make_terrain_geometry();
    const Box domain = geom.Domain();
    const Real requested_distance = Real(1.5);
    const auto faces = all_faces();

    for (std::size_t face_index = 0; face_index < faces.size(); ++face_index) {
        const auto& face = faces[face_index];
        const int dir = face.coordDir();
        const int wall_offset = static_cast<int>(
            std::floor(requested_distance / geom.CellSize(dir) - myhalf));
        const int expected_index = face.isLow()
            ? domain.smallEnd(dir) + wall_offset
            : domain.bigEnd(dir) - wall_offset;
        const Real expected_theta = Real(100.0) + Real(10.0) * expected_index;

        for (const int policy : {0, 1}) {
            MOSTAverageFields fields(geom);
            for (MFIter mfi(*fields.theta, false); mfi.isValid(); ++mfi) {
                auto theta_arr = fields.theta->array(mfi);
                ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const int normal_index = (dir == 0) ? i : ((dir == 1) ? j : k);
                    theta_arr(i,j,k) = Real(100.0) + Real(10.0) * normal_index;
                });
            }
            Gpu::streamSynchronize();

            const std::string prefix =
                "unit_most_terrain_average_" + std::to_string(face_index) +
                "_" + std::to_string(policy);
            ScopedMOSTParams params(prefix.c_str(), policy, 0);
            params.add_zref(requested_distance);

            BoxArray ba_nd(fields.ba);
            ba_nd.surroundingNodes();
            auto z_phys_nd = std::make_unique<MultiFab>(ba_nd, fields.dm, 1, 1);
            initialize_uniform_terrain_height(*z_phys_nd, geom);

            MOSTAverage averages(face, Vector<Geometry>{geom}, false, prefix,
                                  MeshType::ConstantDz,
                                  TerrainType::StaticFittedMesh);
            averages.make_MOSTAverage_at_level(
                0, fields.vars_old, fields.theta, fields.qv,
                fields.qr, z_phys_nd);
            Gpu::streamSynchronize();

            const auto* k_indices = averages.get_k_indices(0);
            EXPECT_EQ(k_indices->min(0), expected_index)
                << "policy=" << policy << ", direction=" << dir
                << ", high=" << !face.isLow();
            EXPECT_EQ(k_indices->max(0), expected_index)
                << "policy=" << policy << ", direction=" << dir
                << ", high=" << !face.isLow();

            averages.compute_averages(0);
            Gpu::streamSynchronize();
            if (policy == 0) {
                const auto plane_average = averages.get_plane_average(0);
                ASSERT_GE(plane_average.size(), static_cast<std::size_t>(4));
                EXPECT_NEAR(plane_average[3], expected_theta,
                            tolerance(expected_theta))
                    << "direction=" << dir << ", high=" << !face.isLow();
            } else {
                EXPECT_NEAR(averages.get_average(0, 3)->min(0), expected_theta,
                            tolerance(expected_theta))
                    << "direction=" << dir << ", high=" << !face.isLow();
            }
        }
    }
}
