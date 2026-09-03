#include "ERF_GTestSurfaceLayerParallelCommon.H"

#include "../ERF_GTestAssertions.H"

#include <gtest/gtest.h>

#include <array>
#include <string>

using namespace amrex;
using namespace erf_surface_layer_test;

namespace {

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

struct StressComponents
{
    MultiFab* required;
    MultiFab* transpose;
};

StressComponents
stress_components (SurfaceLayerFields& fields,
                   const int normal_dir,
                   const int tangential_dir)
{
    if (normal_dir == 0 && tangential_dir == 1) {
        return {fields.tau[TauType::tau21].get(), fields.tau[TauType::tau12].get()};
    } else if (normal_dir == 0 && tangential_dir == 2) {
        return {fields.tau[TauType::tau31].get(), fields.tau[TauType::tau13].get()};
    } else if (normal_dir == 1 && tangential_dir == 0) {
        return {fields.tau[TauType::tau12].get(), fields.tau[TauType::tau21].get()};
    } else if (normal_dir == 1 && tangential_dir == 2) {
        return {fields.tau[TauType::tau32].get(), fields.tau[TauType::tau23].get()};
    } else if (normal_dir == 2 && tangential_dir == 0) {
        return {fields.tau[TauType::tau13].get(), fields.tau[TauType::tau31].get()};
    }
    return {fields.tau[TauType::tau23].get(), fields.tau[TauType::tau32].get()};
}

void
expect_lateral_stresses (SurfaceLayerFields& fields, const Orientation face)
{
    const int dir = face.coordDir();
    const std::array<int, 2> tangential_dirs = dir == 0
        ? std::array<int, 2>{1, 2}
        : std::array<int, 2>{0, 2};
    for (const int tangential_dir : tangential_dirs) {
        const auto stresses = stress_components(fields, dir, tangential_dir);
        const auto required_range =
            face_range(*stresses.required, fields.domain, face);
        const auto transpose_range =
            face_range(*stresses.transpose, fields.domain, face);
        EXPECT_GT(required_range.count, 0);
        EXPECT_TRUE(is_changed(required_range.lo));
        EXPECT_TRUE(is_changed(required_range.hi));
        EXPECT_TRUE(is_changed(transpose_range.lo));
        EXPECT_TRUE(is_changed(transpose_range.hi));
        EXPECT_EQ(required_range.sentinel_count,
                  required_range.total_count - required_range.count);
        EXPECT_EQ(transpose_range.sentinel_count,
                  transpose_range.total_count - transpose_range.count);
        ERF_EXPECT_NEAR(required_range.lo, transpose_range.lo,
                        halo_tolerance(required_range.lo));
        ERF_EXPECT_NEAR(required_range.hi, transpose_range.hi,
                        halo_tolerance(required_range.hi));
    }
}

TEST(SurfaceLayerParallel, DistributedFaceStressIsFaceOwnedAndMatchesSerialReference)
{
    ScopedMFIterTileSize tile_size(IntVect(AMREX_D_DECL(4, 4, 1024)));
    ScopedSurfaceLayerParams params("unit_surface_layer_parallel");

    for (const auto& face : all_faces()) {
        SurfaceLayerFields fields;
        SurfaceLayerFields reference(true);
        auto layer = fields.prepare_layer(
            face, active_face(face), "unit_surface_layer_parallel",
            false, true, false);
        auto reference_layer = reference.prepare_layer(
            face, active_face(face), "unit_surface_layer_parallel",
            false, true, false);
        fields.impose(*layer);
        reference.impose(*reference_layer);

        const MultiFab* required_a = nullptr;
        const MultiFab* required_b = nullptr;
        const MultiFab* transpose_a = nullptr;
        const MultiFab* transpose_b = nullptr;
        const int dir = face.coordDir();
        if (dir == 0) {
            required_a = fields.tau[TauType::tau21].get();
            required_b = fields.tau[TauType::tau31].get();
            transpose_a = fields.tau[TauType::tau12].get();
            transpose_b = fields.tau[TauType::tau13].get();
        } else if (dir == 1) {
            required_a = fields.tau[TauType::tau12].get();
            required_b = fields.tau[TauType::tau32].get();
            transpose_a = fields.tau[TauType::tau21].get();
            transpose_b = fields.tau[TauType::tau23].get();
        } else {
            required_a = fields.tau[TauType::tau13].get();
            required_b = fields.tau[TauType::tau23].get();
            transpose_a = fields.tau[TauType::tau31].get();
            transpose_b = fields.tau[TauType::tau32].get();
        }

        const auto required_a_range =
            face_owned_range(*required_a, fields.ba, fields.domain, face);
        const auto required_b_range =
            face_owned_range(*required_b, fields.ba, fields.domain, face);
        const auto transpose_a_range =
            face_owned_range(*transpose_a, fields.ba, fields.domain, face);
        const auto transpose_b_range =
            face_owned_range(*transpose_b, fields.ba, fields.domain, face);

        EXPECT_GT(required_a_range.count, 0);
        EXPECT_GT(required_b_range.count, 0);
        EXPECT_TRUE(is_changed(required_a_range.lo));
        EXPECT_TRUE(is_changed(required_a_range.hi));
        EXPECT_TRUE(is_changed(required_b_range.lo));
        EXPECT_TRUE(is_changed(required_b_range.hi));
        const auto reference_a_range = face_owned_range(
            *reference.tau[dir == 0 ? TauType::tau21 :
                           (dir == 1 ? TauType::tau12 : TauType::tau13)],
            reference.ba, reference.domain, face);
        const auto reference_b_range = face_owned_range(
            *reference.tau[dir == 0 ? TauType::tau31 :
                           (dir == 1 ? TauType::tau32 : TauType::tau23)],
            reference.ba, reference.domain, face);
        // Nodal stress FABs count shared interface nodes once per FAB in a
        // multi-box layout. Compare values and ownership below, but do not
        // compare raw point counts with the single-box reference.
        EXPECT_NEAR(required_a_range.lo, reference_a_range.lo, Real(1.e-10));
        EXPECT_NEAR(required_a_range.hi, reference_a_range.hi, Real(1.e-10));
        EXPECT_NEAR(required_b_range.lo, reference_b_range.lo, Real(1.e-10));
        EXPECT_NEAR(required_b_range.hi, reference_b_range.hi, Real(1.e-10));
        EXPECT_EQ(required_a_range.sentinel_count,
                  required_a_range.total_count - required_a_range.count);
        EXPECT_EQ(required_b_range.sentinel_count,
                  required_b_range.total_count - required_b_range.count);
        EXPECT_EQ(transpose_a_range.sentinel_count,
                  transpose_a_range.total_count - transpose_a_range.count);
        EXPECT_EQ(transpose_b_range.sentinel_count,
                  transpose_b_range.total_count - transpose_b_range.count);
        EXPECT_NEAR(required_a_range.lo, transpose_a_range.lo, Real(1.e-10));
        EXPECT_NEAR(required_a_range.hi, transpose_a_range.hi, Real(1.e-10));
        EXPECT_NEAR(required_b_range.lo, transpose_b_range.lo, Real(1.e-10));
        EXPECT_NEAR(required_b_range.hi, transpose_b_range.hi, Real(1.e-10));

        const auto reference_transpose_a_range = face_owned_range(
            *reference.tau[dir == 0 ? TauType::tau12 :
                           (dir == 1 ? TauType::tau21 : TauType::tau31)],
            reference.ba, reference.domain, face);
        const auto reference_transpose_b_range = face_owned_range(
            *reference.tau[dir == 0 ? TauType::tau13 :
                           (dir == 1 ? TauType::tau23 : TauType::tau32)],
            reference.ba, reference.domain, face);
        EXPECT_NEAR(transpose_a_range.lo, reference_transpose_a_range.lo,
                    Real(1.e-10));
        EXPECT_NEAR(transpose_a_range.hi, reference_transpose_a_range.hi,
                    Real(1.e-10));
        EXPECT_NEAR(transpose_b_range.lo, reference_transpose_b_range.lo,
                    Real(1.e-10));
        EXPECT_NEAR(transpose_b_range.hi, reference_transpose_b_range.hi,
                    Real(1.e-10));
    }
}

// Motivation: lateral surface parameters are computed only on face-owned
// grids. Their tangential halos must nevertheless be communicated without
// importing the sentinel values from interior grids that share the collapsed
// wall plane.
TEST(SurfaceLayerParallel, LateralSurfaceParameterGhostsAreFaceOwned)
{
    ScopedMFIterTileSize tile_size(IntVect(AMREX_D_DECL(4, 4, 1024)));
    ScopedSurfaceLayerParams params("unit_surface_layer_parallel_halos");

    for (const auto face : {
             Orientation(Direction::x, Orientation::low),
             Orientation(Direction::x, Orientation::high),
             Orientation(Direction::y, Orientation::low),
             Orientation(Direction::y, Orientation::high)}) {
        SCOPED_TRACE(std::string("direction=") +
                     std::to_string(face.coordDir()) +
                     ", high=" + std::to_string(!face.isLow()));
        SurfaceLayerFields fields;
        auto layer = fields.prepare_layer(
            face, active_face(face), "unit_surface_layer_parallel_halos",
            true, true, true);

        const auto pairs = tangential_halo_pairs(
            *layer->get_u_star(0), fields.ba, fields.domain, face);
        EXPECT_GT(pairs.size(), std::size_t(0));
        if (pairs.empty()) { continue; }
        const auto& pmap = layer->get_u_star(0)->DistributionMap().ProcessorMap();
        bool has_cross_rank_pair = false;
        for (const auto& pair : pairs) {
            if (pmap[pair.target_fab] != pmap[pair.neighbor_fab]) {
                has_cross_rank_pair = true;
                break;
            }
        }
        EXPECT_TRUE(has_cross_rank_pair);

        const std::array<const MultiFab*, 4> parameters{{
            layer->get_u_star(0), layer->get_t_star(0),
            layer->get_q_star(0), layer->get_olen(0)}};
        for (std::size_t parameter = 0; parameter < parameters.size(); ++parameter) {
            for (const auto& pair : pairs) {
                const Real halo = global_fab_value(
                    *parameters[parameter], pair.target_fab, pair.point, true);
                const Real reference = global_fab_value(
                    *parameters[parameter], pair.neighbor_fab, pair.point, false);
                ERF_EXPECT_NEAR(halo, reference, halo_tolerance(reference))
                    << "parameter=" << parameter
                    << ", target_fab=" << pair.target_fab
                    << ", neighbor_fab=" << pair.neighbor_fab;
            }
        }

        const int dir = face.coordDir();
        for (int ibox = 0; ibox < fields.ba.size(); ++ibox) {
            const Box& source = fields.ba[ibox];
            const bool selected = face.isLow()
                ? source.smallEnd(dir) == fields.domain.smallEnd(dir)
                : source.bigEnd(dir) == fields.domain.bigEnd(dir);
            if (selected) { continue; }

            const Box target_box = layer->get_u_star(0)->boxArray()[ibox];
            const IntVect point = target_box.smallEnd();
            EXPECT_EQ(global_fab_value(*layer->get_u_star(0), ibox, point, false),
                      seeded_value(u_star_seed, ibox));
            EXPECT_EQ(global_fab_value(*layer->get_t_star(0), ibox, point, false),
                      seeded_value(t_star_seed, ibox));
            EXPECT_EQ(global_fab_value(*layer->get_q_star(0), ibox, point, false),
                      seeded_value(q_star_seed, ibox));
            EXPECT_EQ(global_fab_value(*layer->get_olen(0), ibox, point, false),
                      seeded_value(olen_seed, ibox));
        }

        fields.impose(*layer);
        expect_lateral_stresses(fields, face);
    }
}

// Motivation: qsat must be written only by tiles that contain the selected
// boundary plane.  The distributed layout deliberately leaves other FABs on
// each rank, so a valid result cannot be established from one local FAB.
TEST(SurfaceLayerParallel, DistributedQsurfUpdatesSelectedFace)
{
    ScopedMFIterTileSize tile_size(IntVect(AMREX_D_DECL(4, 4, 1024)));
    ScopedSurfaceLayerParams params("unit_surface_layer_parallel_qsurf");

    for (const auto& face : all_faces()) {
        SurfaceLayerFields fields;
        fields.lmask[0]->setVal(0);
        auto layer = fields.prepare_layer(
            face, active_face(face), "unit_surface_layer_parallel_qsurf", true);
        const auto* qsurf = layer->get_q_surf(0);
        const auto counts = value_counts(*qsurf);

        Long expected_finite = 0;
        const int dir = face.coordDir();
        for (int ibox = 0; ibox < fields.ba.size(); ++ibox) {
            const Box& source = fields.ba[ibox];
            const bool selected = face.isLow()
                ? source.smallEnd(dir) == fields.domain.smallEnd(dir)
                : source.bigEnd(dir) == fields.domain.bigEnd(dir);
            if (selected) {
                expected_finite += static_cast<Long>(
                    qsurf->boxArray()[ibox].numPts());
            }
        }

        EXPECT_GT(counts.finite, 0);
        EXPECT_EQ(counts.finite, expected_finite)
            << "direction=" << dir << ", high=" << !face.isLow();
        EXPECT_EQ(counts.sentinel, counts.total - expected_finite)
            << "direction=" << dir << ", high=" << !face.isLow();
        EXPECT_EQ(counts.finite + counts.sentinel, counts.total)
            << "direction=" << dir << ", high=" << !face.isLow();

        const auto selected = face_owned_range(
            *qsurf, fields.ba, fields.domain, face);
        EXPECT_EQ(selected.count, counts.finite)
            << "direction=" << dir << ", high=" << !face.isLow();
        const Real expected = expected_qsat(fields.geom);
        EXPECT_NEAR(selected.lo, expected, qsat_tolerance(expected))
            << "direction=" << dir << ", high=" << !face.isLow();
        EXPECT_NEAR(selected.hi, expected, qsat_tolerance(expected))
            << "direction=" << dir << ", high=" << !face.isLow();
    }
}

// Motivation: shared corners are processed by more than one active face. The
// second face may add its transpose value, but must not erase the first face's
// required normal stress. This is the distributed version of the serial
// corner contract.
TEST(SurfaceLayerParallel, DistributedMixedFaceCornersPreserveBothStresses)
{
    ScopedSurfaceLayerParams params("unit_surface_layer_parallel_corners");
    const std::array<std::array<Orientation, 2>, 3> direction_pairs{{
        {{Orientation(Direction::x, Orientation::low),
          Orientation(Direction::y, Orientation::low)}},
        {{Orientation(Direction::x, Orientation::high),
          Orientation(Direction::z, Orientation::high)}},
        {{Orientation(Direction::y, Orientation::high),
          Orientation(Direction::z, Orientation::low)}}
    }};

    for (const auto& pair : direction_pairs) {
        SurfaceLayerFields fields;
        GpuArray<int, AMREX_SPACEDIM*2> active{};
        active[static_cast<int>(pair[0])] = 1;
        active[static_cast<int>(pair[1])] = 1;
        auto first = fields.prepare_layer(
            pair[0], active, "unit_surface_layer_parallel_corners");
        fields.impose(*first);

        const int first_dir = pair[0].coordDir();
        const int second_dir = pair[1].coordDir();
        IntVect corner(fields.domain.smallEnd());
        corner[first_dir] = pair[0].isLow()
            ? fields.domain.smallEnd(first_dir) : fields.domain.bigEnd(first_dir) + 1;
        corner[second_dir] = pair[1].isLow()
            ? fields.domain.smallEnd(second_dir) : fields.domain.bigEnd(second_dir) + 1;
        const int remaining_dir = 3 - first_dir - second_dir;
        corner[remaining_dir] = fields.domain.smallEnd(remaining_dir) + 1;

        const auto first_stress = stress_components(fields, first_dir, second_dir);
        const auto second_stress = stress_components(fields, second_dir, first_dir);
        const Real first_value = mf_value(*first_stress.required, corner);
        EXPECT_TRUE(is_changed(first_value));
        EXPECT_EQ(mf_value(*first_stress.transpose, corner), tau_sentinel);

        auto second = fields.prepare_layer(
            pair[1], active, "unit_surface_layer_parallel_corners");
        fields.impose(*second);
        EXPECT_TRUE(is_changed(mf_value(*second_stress.required, corner)));
        EXPECT_NEAR(mf_value(*second_stress.transpose, corner), first_value,
                    Real(1.e-10));
        EXPECT_NEAR(mf_value(*first_stress.required, corner), first_value,
                    Real(1.e-10));
    }
}

} // namespace
