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

namespace {

// f(i,j) = 1 + i + 2 j in every cell of every planar box, valid region and ghost cells, so every
// duplicate copy of a surface cell holds the same value, as after fill_planar_boundary
void set_planar_pattern (MultiFab& mf)
{
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto arr = mf.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            arr(i,j,k) = Real(1.0) + static_cast<Real>(i) + Real(2.0) * static_cast<Real>(j);
        });
    }
    Gpu::streamSynchronize();
}

} // namespace

// On grids split in z the planar surface-layer arrays hold one duplicate box per stacked 3D box;
// surface_sum must count each surface cell once (the surface history averages were multiplied by
// the number of stacked boxes).
TEST(SurfaceLayerParallel, SurfaceSumCountsEachSurfaceCellOnceOnZSplitGrids)
{
    ScopedSurfaceLayerParams params("unit_surface_layer_zsplit_sum");
    const Orientation zlo(Direction::z, Orientation::low);

    // 32 x 32 x 4 cells in 16 x 16 x 2 boxes: 4 columns of 2 stacked boxes
    SurfaceLayerFields fields(false, IntVect(AMREX_D_DECL(16, 16, 2)));
    ASSERT_EQ(static_cast<int>(fields.ba.size()), 8);
    auto layer = fields.prepare_layer(zlo, active_face(zlo), "unit_surface_layer_zsplit_sum",
                                      false, false, false, false);
    MultiFab& ustar = *layer->get_u_star(0);
    ASSERT_EQ(static_cast<int>(ustar.boxArray().size()), 8);

    const Real ncell = Real(32 * 32);
    // Sum over i, j in [0,31] of 1 + i + 2 j
    const Real pattern_sum = ncell + Real(32 * 496) + Real(2 * 32 * 496);

    ustar.setVal(Real(1.0));
    // The plain sum counts every surface cell twice, so the layout does duplicate the planar boxes
    ERF_EXPECT_NEAR(ustar.sum(0), Real(2.0) * ncell, Real(1.0e-10) * ncell);
    ERF_EXPECT_NEAR(layer->surface_sum(0, ustar), ncell, Real(1.0e-10) * ncell);

    set_planar_pattern(ustar);
    ERF_EXPECT_NEAR(layer->surface_sum(0, ustar), pattern_sum, Real(1.0e-10) * pattern_sum);

    // Without a split the planar boxes are not duplicated and the two sums agree
    SurfaceLayerFields unsplit;
    auto unsplit_layer = unsplit.prepare_layer(zlo, active_face(zlo), "unit_surface_layer_zsplit_sum",
                                               false, false, false, false);
    MultiFab& ustar_unsplit = *unsplit_layer->get_u_star(0);
    set_planar_pattern(ustar_unsplit);
    ERF_EXPECT_NEAR(unsplit_layer->surface_sum(0, ustar_unsplit), pattern_sum,
                    Real(1.0e-10) * pattern_sum);
    ERF_EXPECT_NEAR(ustar_unsplit.sum(0), pattern_sum, Real(1.0e-10) * pattern_sum);
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
        const bool expect_cross_rank_pair =
            ParallelDescriptor::NProcs() > 1;
        EXPECT_EQ(has_cross_rank_pair, expect_cross_rank_pair);

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
        auto layer = fields.prepare_layer(
            face, active_face(face), "unit_surface_layer_parallel_qsurf",
            true, false, false, false);
        fields.lmask[0]->setVal(0);
        layer->get_t_surf(0)->setVal(test_surface_temperature);
        std::unique_ptr<MultiFab> z_phys_nd;
        layer->fill_qsurf_with_qsat(0, fields.cons, z_phys_nd);
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

namespace {

// theta = 300 K in the two lowest cells and 302 K above, rho = 1, no TKE, over every cell of every
// box including ghost cells.  With dz = 1 m the MYNN25 estimator puts the PBL height where theta_v
// reaches min + 1.25 K between the cell centres at 1.5 m and 2.5 m: 1.5 + 1.25/2 = 2.125 m.
//
// k_jump moves the jump: theta = 302 K from cell k_jump up, and the height is k_jump + 0.125 m.
void set_theta_jump (MultiFab& cons, const int k_jump = 2)
{
    for (MFIter mfi(cons); mfi.isValid(); ++mfi) {
        auto arr = cons.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            arr(i,j,k,Rho_comp)      = Real(1.0);
            arr(i,j,k,RhoTheta_comp) = (k >= k_jump) ? Real(302.0) : Real(300.0);
            arr(i,j,k,RhoKE_comp)    = Real(0.0);
        });
    }
    Gpu::streamSynchronize();
}

// The same jump of 2 K in theta_v made by water vapour alone: theta = 300 K everywhere, and qv
// steps from zero to 2 / (300 epsv) at k = 2.  The estimator only sees it if it is handed qv.
void set_vapour_jump (MultiFab& cons)
{
    const Real qv_above = Real(2.0) / (Real(300.0) * epsv);
    for (MFIter mfi(cons); mfi.isValid(); ++mfi) {
        auto arr = cons.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            arr(i,j,k,Rho_comp)      = Real(1.0);
            arr(i,j,k,RhoTheta_comp) = Real(300.0);
            arr(i,j,k,RhoKE_comp)    = Real(0.0);
            arr(i,j,k,RhoQ1_comp)    = (k >= 2) ? qv_above : Real(0.0);
        });
    }
    Gpu::streamSynchronize();
}

std::pair<Real,Real> valid_min_max (const MultiFab& mf)
{
    ReduceOps<ReduceOpMin, ReduceOpMax> reduce_op;
    ReduceData<Real, Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const auto arr = mf.const_array(mfi);
        reduce_op.eval(mfi.validbox(), reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                return {arr(i,j,k), arr(i,j,k)};
            });
    }
    auto result = reduce_data.value(reduce_op);
    Real lo = get<0>(result);
    Real hi = get<1>(result);
    ParallelDescriptor::ReduceRealMin(lo);
    ParallelDescriptor::ReduceRealMax(hi);
    return {lo, hi};
}

Real pblh_on_layout (const bool single_box, const IntVect& max_size, int& nplanar,
                     std::pair<Real,Real>& range, const bool vapour_jump = false)
{
    const Orientation zlo(Direction::z, Orientation::low);
    SurfaceLayerFields fields(single_box, max_size);
    MoistureComponentIndices moisture_indices;
    if (vapour_jump) {
        set_vapour_jump(fields.cons);
        moisture_indices.qv = RhoQ1_comp;
    } else {
        set_theta_jump(fields.cons);
    }
    auto layer = fields.prepare_layer(zlo, active_face(zlo), "unit_surface_layer_zsplit_pblh",
                                      false, false, false, false);

    Vector<Vector<MultiFab>> vars(1);
    vars[0].resize(Vars::NumTypes);
    vars[0][Vars::cons] = MultiFab(fields.cons, amrex::make_alias, 0, fields.cons.nComp());
    layer->update_pblh(0, vars, nullptr, moisture_indices);

    const MultiFab& pblh = *layer->get_pblh(0);
    nplanar = static_cast<int>(pblh.boxArray().size());
    range = valid_min_max(pblh);
    return range.first;
}

// Largest difference, over the valid cells of every planar copy of pblh, from a field that is
// `found` inside the box `where` but outside the box `hole`, and zero elsewhere
Real pblh_error_on_grids (const BoxArray& grids, const int k_jump, const Box& where,
                          const Box& hole, const Real found, int& nplanar)
{
    const Orientation zlo(Direction::z, Orientation::low);
    SurfaceLayerFields fields(grids);
    set_theta_jump(fields.cons, k_jump);
    auto layer = fields.prepare_layer(zlo, active_face(zlo), "unit_surface_layer_zsplit_pblh",
                                      false, false, false, false);

    Vector<Vector<MultiFab>> vars(1);
    vars[0].resize(Vars::NumTypes);
    vars[0][Vars::cons] = MultiFab(fields.cons, amrex::make_alias, 0, fields.cons.nComp());
    layer->update_pblh(0, vars, nullptr, MoistureComponentIndices{});

    const MultiFab& pblh = *layer->get_pblh(0);
    nplanar = static_cast<int>(pblh.boxArray().size());

    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (MFIter mfi(pblh); mfi.isValid(); ++mfi) {
        const auto arr = pblh.const_array(mfi);
        reduce_op.eval(mfi.validbox(), reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                const IntVect iv(AMREX_D_DECL(i, j, 0));
                const Real expected = (where.contains(iv) && !hole.contains(iv)) ? found : Real(0.0);
                return {std::abs(arr(i,j,k) - expected)};
            });
    }
    Real err = get<0>(reduce_data.value(reduce_op));
    ParallelDescriptor::ReduceRealMax(err);
    return err;
}

// The four lateral quarters of the 32 x 32 domain over the cells k_lo to k_hi
BoxList quarter_boxes (const int k_lo, const int k_hi)
{
    BoxList bl;
    for (int jq = 0; jq < 2; ++jq) {
        for (int iq = 0; iq < 2; ++iq) {
            bl.push_back(Box(IntVect(AMREX_D_DECL(16*iq,    16*jq,    k_lo)),
                             IntVect(AMREX_D_DECL(16*iq+15, 16*jq+15, k_hi))));
        }
    }
    return bl;
}

} // namespace

// The MYNN25 PBL-height estimator scans whole columns.  On grids split in z it must still see the
// whole column, and every duplicate planar copy of pblh must hold the same, exact height.
TEST(SurfaceLayerParallel, PBLHeightScansWholeColumnsOnZSplitGrids)
{
    ScopedSurfaceLayerParams params("unit_surface_layer_zsplit_pblh");
    ParmParse pp("unit_surface_layer_zsplit_pblh");
    pp.add("most.pblh_calc", std::string("MYNN25"));

    const Real expected = Real(2.125);
    const Real tol = Real(1.0e-10);

    int nplanar = 0;
    std::pair<Real,Real> range;

    pblh_on_layout(true, IntVect(AMREX_D_DECL(16, 16, 1024)), nplanar, range);
    EXPECT_EQ(nplanar, 1);
    ERF_EXPECT_NEAR(range.first,  expected, tol);
    ERF_EXPECT_NEAR(range.second, expected, tol);

    // 4 columns of 2 stacked boxes: the upper box of each column holds no cell below k = 2
    pblh_on_layout(false, IntVect(AMREX_D_DECL(16, 16, 2)), nplanar, range);
    EXPECT_EQ(nplanar, 8);
    ERF_EXPECT_NEAR(range.first,  expected, tol);
    ERF_EXPECT_NEAR(range.second, expected, tol);

    pp.remove("most.pblh_calc");
}

// Grids that stop below the top of the domain, as on a refined level of partial height: the scan
// must stop with the box.  The cell above the box is a ghost cell, which holds data, so a jump
// between the top cell and that ghost cell is found and a jump above it is not.
TEST(SurfaceLayerParallel, PBLHeightStopsAtTheTopOfPartialHeightGrids)
{
    ScopedSurfaceLayerParams params("unit_surface_layer_zsplit_pblh");
    ParmParse pp("unit_surface_layer_zsplit_pblh");
    pp.add("most.pblh_calc", std::string("MYNN25"));

    const Real tol = Real(1.0e-10);
    const Box everywhere(IntVect(AMREX_D_DECL(0, 0, 0)), IntVect(AMREX_D_DECL(31, 31, 0)));
    const Box nowhere;
    const BoxArray lower_half(quarter_boxes(0, 1));
    int nplanar = 0;

    ERF_EXPECT_NEAR(pblh_error_on_grids(lower_half, 2, everywhere, nowhere, Real(2.125), nplanar),
                    Real(0.0), tol);
    EXPECT_EQ(nplanar, 4);

    // The jump sits between k = 2 and k = 3, which these grids do not hold
    ERF_EXPECT_NEAR(pblh_error_on_grids(lower_half, 3, nowhere, nowhere, Real(0.0), nplanar),
                    Real(0.0), tol);

    pp.remove("most.pblh_calc");
}

// Boxes stacked in z with different footprints, and one that is partly aloft: three quarters of
// the domain hold k = 0, 1 and a box over 8 <= i <= 23 holds k = 2, 3, also over the fourth
// quarter, where nothing lies below it.  The jump between k = 2 and k = 3 is found only in the
// columns that reach the ground and hold those cells; every planar copy agrees, and the columns
// over the fourth quarter report zero, the value for a height that was not found.
TEST(SurfaceLayerParallel, PBLHeightOnUnevenStacksAndBoxesAloft)
{
    ScopedSurfaceLayerParams params("unit_surface_layer_zsplit_pblh");
    ParmParse pp("unit_surface_layer_zsplit_pblh");
    pp.add("most.pblh_calc", std::string("MYNN25"));

    const Real tol = Real(1.0e-10);
    BoxList bl = quarter_boxes(0, 1);
    const Box fourth_quarter = bl.data().back();
    bl.data().pop_back();
    bl.push_back(Box(IntVect(AMREX_D_DECL(8, 0, 2)), IntVect(AMREX_D_DECL(23, 31, 3))));

    const Box under_the_upper_box(IntVect(AMREX_D_DECL(8, 0, 0)), IntVect(AMREX_D_DECL(23, 31, 0)));
    Box hole(fourth_quarter);
    hole.setRange(2, 0);

    int nplanar = 0;
    ERF_EXPECT_NEAR(pblh_error_on_grids(BoxArray(bl), 3, under_the_upper_box, hole, Real(3.125),
                                        nplanar),
                    Real(0.0), tol);
    EXPECT_EQ(nplanar, 4);

    pp.remove("most.pblh_calc");
}

// A level none of whose boxes reaches the ground has no column to scan: the PBL height is zero
// on every planar box, not the value the field was allocated with.
TEST(SurfaceLayerParallel, PBLHeightIsZeroOnGridsEntirelyAloft)
{
    ScopedSurfaceLayerParams params("unit_surface_layer_zsplit_pblh");
    ParmParse pp("unit_surface_layer_zsplit_pblh");
    pp.add("most.pblh_calc", std::string("MYNN25"));

    const Box nowhere;
    int nplanar = 0;
    ERF_EXPECT_NEAR(pblh_error_on_grids(BoxArray(quarter_boxes(2, 3)), 3, nowhere, nowhere,
                                        Real(0.0), nplanar),
                    Real(0.0), Real(1.0e-10));
    EXPECT_EQ(nplanar, 4);

    pp.remove("most.pblh_calc");
}

// The columns hold only the components the estimator reads.  A jump in theta_v that water
// vapour alone makes is found on grids split in z only if the moisture species are among them.
TEST(SurfaceLayerParallel, PBLHeightSeesMoistureOnZSplitGrids)
{
    ScopedSurfaceLayerParams params("unit_surface_layer_zsplit_pblh");
    ParmParse pp("unit_surface_layer_zsplit_pblh");
    pp.add("most.pblh_calc", std::string("MYNN25"));

    // 1.5 + 1.25 / (300 epsv qv) goes through a division that is not exact
    const Real expected = Real(2.125);
    const Real tol = Real(1.0e3) * std::numeric_limits<Real>::epsilon();

    int nplanar = 0;
    std::pair<Real,Real> range;

    pblh_on_layout(true, IntVect(AMREX_D_DECL(16, 16, 1024)), nplanar, range, true);
    EXPECT_EQ(nplanar, 1);
    ERF_EXPECT_NEAR(range.first,  expected, tol);
    ERF_EXPECT_NEAR(range.second, expected, tol);

    pblh_on_layout(false, IntVect(AMREX_D_DECL(16, 16, 2)), nplanar, range, true);
    EXPECT_EQ(nplanar, 8);
    ERF_EXPECT_NEAR(range.first,  expected, tol);
    ERF_EXPECT_NEAR(range.second, expected, tol);

    pp.remove("most.pblh_calc");
}
