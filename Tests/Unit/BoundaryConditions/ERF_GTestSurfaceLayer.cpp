#include "AMReX_Gpu.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFab.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Reduce.H"

#include "ERF_SurfaceLayer.H"
#include "ERF_GTestSurfaceLayerCommon.H"

#include <gtest/gtest.h>

#include "../ERF_GTestAssertions.H"

#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <memory>
#include <string>

using namespace amrex;
using erf_surface_layer_test::expected_qsat;
using erf_surface_layer_test::qsat_tolerance;
using erf_surface_layer_test::stress_has_expected_sign;
using erf_surface_layer_test::stress_is_antisymmetric;
using erf_surface_layer_test::tau_sentinel;
using erf_surface_layer_test::test_rho;
using erf_surface_layer_test::test_rho_theta;
using erf_surface_layer_test::test_qv;
using erf_surface_layer_test::test_primitive_z_ng;
using erf_surface_layer_test::test_surface_temperature;
using erf_surface_layer_test::test_state_ng;
using erf_surface_layer_test::test_u;
using erf_surface_layer_test::test_v;
using erf_surface_layer_test::test_w;
using erf_surface_layer_test::test_velocity_ng;

namespace {

Geometry
make_geometry ()
{
    const Box domain(IntVect(2, 2, 2), IntVect(4, 4, 4));
    const RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                           {AMREX_D_DECL(3.0, 3.0, 3.0)});
    const Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
    return Geometry(domain, &real_box, 0, is_periodic.data());
}

Geometry
make_qsurf_geometry ()
{
    const Box domain(IntVect(0), IntVect(2));
    const RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                           {AMREX_D_DECL(3.0, 3.0, 3.0)});
    const Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
    return Geometry(domain, &real_box, 0, is_periodic.data());
}

BoxArray
make_qsurf_box_array (const Box& domain)
{
    BoxArray ba(domain);
    // Keep complete vertical columns while creating non-face FABs in both
    // lateral directions for the ownership contract.
    ba.maxSize(IntVect(AMREX_D_DECL(1, 1, 1024)));
    return ba;
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

BoxArray
collapse_z (const BoxArray& ba)
{
    BoxList boxes = ba.boxList();
    for (auto& box : boxes) {
        box.setRange(2, 0);
    }
    return BoxArray(std::move(boxes));
}

Real
mf_value (const MultiFab& mf, const IntVect& point)
{
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    const Box point_box(point, point, mf.boxArray().ixType());
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box overlap = point_box & mfi.validbox();
        if (overlap.isEmpty()) { continue; }
        const auto array = mf.const_array(mfi);
        reduce_op.eval(overlap, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> GpuTuple<Real>
            {
                return {array(i,j,k)};
            });
    }
    Gpu::streamSynchronize();
    return get<0>(reduce_data.value());
}

Real
single_value (const MultiFab& mf, const Box& box, int comp)
{
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box overlap = box & mfi.validbox();
        if (overlap.isEmpty()) { continue; }
        const auto array = mf.const_array(mfi);
        reduce_op.eval(overlap, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> GpuTuple<Real>
            {
                return {array(i,j,k,comp)};
            });
    }
    Gpu::streamSynchronize();
    return get<0>(reduce_data.value());
}

bool
is_changed (const Real value)
{
    return std::isfinite(value) && value != tau_sentinel;
}

class ScopedSurfaceLayerParams
{
public:
    explicit ScopedSurfaceLayerParams (const char* prefix)
        : m_pp(prefix)
    {
        m_pp.add("most.average_policy", 0);
        m_pp.add("most.z0", Real(0.1));
    }

    ~ScopedSurfaceLayerParams ()
    {
        m_pp.remove("most.average_policy");
        m_pp.remove("most.z0");
    }

private:
    ParmParse m_pp;
};

struct SurfaceLayerFields
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

    Vector<std::unique_ptr<iMultiFab>> lmask;
    std::unique_ptr<MultiFab> no_walldist;

    Vector<std::unique_ptr<MultiFab>> tau;
    MultiFab xheat_flux;
    MultiFab yheat_flux;
    MultiFab zheat_flux;

    Vector<MultiFab*> state;

    explicit SurfaceLayerFields (const Geometry& geometry = make_geometry(),
                                 const bool split_lateral = false)
        : geom(geometry),
          domain(geom.Domain()),
          ba(split_lateral ? make_qsurf_box_array(domain) : BoxArray(domain)),
          dm(ba),
          cons(ba, dm, RhoQ1_comp + 1, test_state_ng),
          xvel(convert(ba, IntVect(AMREX_D_DECL(1, 0, 0))), dm, 1,
               test_velocity_ng),
          yvel(convert(ba, IntVect(AMREX_D_DECL(0, 1, 0))), dm, 1,
               test_velocity_ng),
          zvel(convert(ba, IntVect(AMREX_D_DECL(0, 0, 1))), dm, 1,
               test_velocity_ng),
          theta(std::make_unique<MultiFab>(
              ba, dm, 1,
              IntVect(AMREX_D_DECL(test_state_ng, test_state_ng,
                                   test_primitive_z_ng)))),
          xheat_flux(convert(ba, IntVect(AMREX_D_DECL(1, 0, 0))), dm, 1, 1),
          yheat_flux(convert(ba, IntVect(AMREX_D_DECL(0, 1, 0))), dm, 1, 1),
          zheat_flux(convert(ba, IntVect(AMREX_D_DECL(0, 0, 1))), dm, 1, 1)
    {
        cons.setVal(Real(0.0));
        cons.setVal(test_rho, Rho_comp, 1);
        cons.setVal(test_rho_theta, RhoTheta_comp, 1);
        cons.setVal(test_rho * test_qv, RhoQ1_comp, 1);
        xvel.setVal(test_u);
        yvel.setVal(test_v);
        zvel.setVal(test_w);
        theta->setVal(test_surface_temperature);

        lmask.emplace_back(std::make_unique<iMultiFab>(
            collapse_z(ba), dm, 1,
            IntVect(AMREX_D_DECL(test_state_ng, test_state_ng, 0))));
        lmask[0]->setVal(1);

        tau.resize(9);
        const BoxArray ba12 = convert(ba, IntVect(1, 1, 0));
        const BoxArray ba13 = convert(ba, IntVect(1, 0, 1));
        const BoxArray ba23 = convert(ba, IntVect(0, 1, 1));
        tau[TauType::tau11] = std::make_unique<MultiFab>(ba, dm, 1, 1);
        tau[TauType::tau22] = std::make_unique<MultiFab>(ba, dm, 1, 1);
        tau[TauType::tau33] = std::make_unique<MultiFab>(ba, dm, 1, 1);
        tau[TauType::tau12] = std::make_unique<MultiFab>(ba12, dm, 1, 1);
        tau[TauType::tau13] = std::make_unique<MultiFab>(ba13, dm, 1, 1);
        tau[TauType::tau23] = std::make_unique<MultiFab>(ba23, dm, 1, 1);
        tau[TauType::tau21] = std::make_unique<MultiFab>(ba12, dm, 1, 1);
        tau[TauType::tau31] = std::make_unique<MultiFab>(ba13, dm, 1, 1);
        tau[TauType::tau32] = std::make_unique<MultiFab>(ba23, dm, 1, 1);

        state = {&cons, &xvel, &yvel, &zvel};
        reset_outputs();
    }

    void reset_outputs ()
    {
        for (auto& stress : tau) {
            stress->setVal(tau_sentinel);
        }
        xheat_flux.setVal(tau_sentinel);
        yheat_flux.setVal(tau_sentinel);
        zheat_flux.setVal(tau_sentinel);
        Gpu::streamSynchronize();
    }

    void set_varying_surface_temperature ()
    {
        for (MFIter mfi(*theta, false); mfi.isValid(); ++mfi) {
            const Box box = mfi.fabbox();
            auto theta_arr = theta->array(mfi);
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                theta_arr(i,j,k) = test_surface_temperature
                    + Real(0.25) * static_cast<Real>(i)
                    + Real(0.5) * static_cast<Real>(j);
            });
        }
        Gpu::streamSynchronize();
    }

    std::unique_ptr<SurfaceLayer>
    prepare_layer (const Orientation face,
                   const GpuArray<int, AMREX_SPACEDIM*2>& active_faces,
                   const std::string& prefix,
                   const bool with_moisture = false)
    {
        bool rotate = false;
        Vector<Geometry> geoms{geom};
        Vector<std::unique_ptr<MultiFab>> qv_prim(1);
        Vector<std::unique_ptr<MultiFab>> z_phys_nd(1);
        if (with_moisture) {
            qv_prim[0] = std::make_unique<MultiFab>(
                ba, dm, 1,
                IntVect(AMREX_D_DECL(test_state_ng, test_state_ng,
                                     test_primitive_z_ng)));
            qv_prim[0]->setVal(test_qv);
        }
        auto layer = std::make_unique<SurfaceLayer>(
            face, geoms, rotate, prefix, qv_prim, z_phys_nd,
            MeshType::ConstantDz, TerrainType::None, TurbChoice{},
            0.0, 0.0);
        layer->set_surface_layer_faces(active_faces);

        std::unique_ptr<MultiFab> qr_prim;
        Vector<MultiFab*> empty_mfs;
        Vector<std::string> empty_names;
        Vector<std::unique_ptr<MultiFab>> sst;
        Vector<std::unique_ptr<MultiFab>> tsk;
        layer->make_SurfaceLayer_at_level(
            0, 1, state, theta, qv_prim[0], qr_prim, z_phys_nd[0],
            nullptr, nullptr, nullptr, empty_mfs, empty_names,
            empty_mfs, empty_names, sst, tsk, lmask);
        if (with_moisture) {
            layer->get_q_surf(0)->setVal(tau_sentinel);
        }
        layer->update_fluxes(0, 0.0, 0.0, cons, z_phys_nd[0],
                             no_walldist, 20);
        return layer;
    }

    void impose (SurfaceLayer& layer)
    {
        Vector<const MultiFab*> const_state{&cons, &xvel, &yvel, &zvel};
        layer.impose_SurfaceLayer_bcs(
            0, const_state, tau, &xheat_flux, &yheat_flux, &zheat_flux,
            nullptr, nullptr, nullptr, nullptr);
        Gpu::streamSynchronize();
    }
};

GpuArray<int, AMREX_SPACEDIM*2>
active_faces (std::initializer_list<Orientation> faces)
{
    GpuArray<int, AMREX_SPACEDIM*2> result{};
    for (const auto face : faces) {
        result[static_cast<int>(face)] = 1;
    }
    return result;
}

IntVect
face_point (const Box& domain, const Orientation face)
{
    IntVect point(domain.smallEnd());
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (dir == face.coordDir()) {
            point[dir] = face.isLow() ? domain.smallEnd(dir)
                                      : domain.bigEnd(dir) + 1;
        } else {
            point[dir] = domain.smallEnd(dir) + 1;
        }
    }
    return point;
}

} // namespace

// Motivation: the Moeng stress functor has separate x-, y-, and z-wall
// interpolation paths.  Exercise each path at both wall orientations so the
// SurfaceLayer caller can rely on the normal high-face index being mapped back
// to the adjacent interior cell.
TEST(SurfaceLayer, MoengDirectionalFluxesAreFiniteOnEveryWall)
{
    const Geometry geom = make_geometry();
    const Box domain = geom.Domain();
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);
    constexpr int ng = 2;

    MultiFab cons(ba, dm, 3, ng);
    MultiFab xvel(BoxArray(surroundingNodes(domain, 0)), dm, 1, ng);
    MultiFab yvel(BoxArray(surroundingNodes(domain, 1)), dm, 1, ng);
    MultiFab zvel(BoxArray(surroundingNodes(domain, 2)), dm, 1, ng);
    MultiFab um(ba, dm, 1, ng);
    MultiFab vm(ba, dm, 1, ng);
    MultiFab wm(ba, dm, 1, ng);
    MultiFab umm(ba, dm, 1, ng);
    MultiFab ustar(ba, dm, 1, ng);
    MultiFab stress(ba, dm, 2, 0);

    cons.setVal(Real(1.0));
    xvel.setVal(Real(3.0));
    yvel.setVal(Real(4.0));
    zvel.setVal(Real(5.0));
    um.setVal(Real(3.0));
    vm.setVal(Real(4.0));
    wm.setVal(Real(5.0));
    umm.setVal(Real(6.4));
    ustar.setVal(Real(0.8));

    const auto faces = all_faces();
    std::array<Real, AMREX_SPACEDIM> low_u{};
    std::array<Real, AMREX_SPACEDIM> high_u{};
    std::array<Real, AMREX_SPACEDIM> low_v{};
    std::array<Real, AMREX_SPACEDIM> high_v{};

    for (const auto& face : faces) {
        const moeng_flux flux(Real(0.1), face.isLow(),
                              domain.smallEnd(2), domain.bigEnd(2));
        const int dir = face.coordDir();
        int i = domain.smallEnd(0) + 1;
        int j = domain.smallEnd(1) + 1;
        int k = domain.smallEnd(2) + 1;
        const int normal_index = face.isLow()
            ? domain.smallEnd(dir) : domain.bigEnd(dir) + 1;
        if (dir == 0) {
            i = normal_index;
        } else if (dir == 1) {
            j = normal_index;
        } else {
            k = normal_index;
        }
        const auto cons_arr = cons[0].const_array();
        const auto xvel_arr = xvel[0].const_array();
        const auto yvel_arr = yvel[0].const_array();
        const auto zvel_arr = zvel[0].const_array();
        const auto um_arr = um[0].const_array();
        const auto vm_arr = vm[0].const_array();
        const auto wm_arr = wm[0].const_array();
        const auto umm_arr = umm[0].const_array();
        const auto ustar_arr = ustar[0].const_array();
        auto stress_arr = stress[0].array();
        const Box output_box(domain.smallEnd(), domain.smallEnd());
        ParallelFor(output_box, [=] AMREX_GPU_DEVICE (int oi, int oj, int ok)
        {
            stress_arr(oi,oj,ok,0) = flux.compute_u_flux(
                i, j, k, dir, cons_arr, xvel_arr, yvel_arr, zvel_arr,
                umm_arr, um_arr, vm_arr, wm_arr, ustar_arr);
            stress_arr(oi,oj,ok,1) = flux.compute_v_flux(
                i, j, k, dir, cons_arr, xvel_arr, yvel_arr, zvel_arr,
                umm_arr, um_arr, vm_arr, wm_arr, ustar_arr);
        });
        Gpu::streamSynchronize();
        const Real stress_u = single_value(stress, output_box, 0);
        const Real stress_v = single_value(stress, output_box, 1);

        EXPECT_TRUE(std::isfinite(stress_u))
            << "direction=" << dir << ", high=" << !face.isLow();
        EXPECT_TRUE(std::isfinite(stress_v))
            << "direction=" << dir << ", high=" << !face.isLow();

        if (face.isLow()) {
            low_u[dir] = stress_u;
            low_v[dir] = stress_v;
        } else {
            high_u[dir] = stress_u;
            high_v[dir] = stress_v;
        }
    }

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        EXPECT_TRUE(stress_has_expected_sign(low_u[dir], true))
            << "direction=" << dir;
        EXPECT_TRUE(stress_has_expected_sign(low_v[dir], true))
            << "direction=" << dir;
        EXPECT_TRUE(stress_is_antisymmetric(low_u[dir], high_u[dir], Real(1.e-12)))
            << "direction=" << dir;
        EXPECT_TRUE(stress_is_antisymmetric(low_v[dir], high_v[dir], Real(1.e-12)))
            << "direction=" << dir;
    }
}

// Motivation: each Cartesian face owns two normal momentum-stress components,
// while an isolated face must still populate the equivalent transpose
// components with the same computed Moeng flux.
TEST(SurfaceLayer, FaceStressIsConsistentForNonconstantInputs)
{
    ScopedSurfaceLayerParams params("unit_surface_layer_faces");
    const auto faces = all_faces();

    for (const auto face : faces) {
        SurfaceLayerFields fields;
        fields.set_varying_surface_temperature();
        auto layer = fields.prepare_layer(
            face, active_faces({face}), "unit_surface_layer_faces");
        fields.impose(*layer);

        const IntVect point = face_point(fields.domain, face);
        const int dir = face.coordDir();
        const MultiFab* required_a = nullptr;
        const MultiFab* required_b = nullptr;
        const MultiFab* transpose_a = nullptr;
        const MultiFab* transpose_b = nullptr;
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

        const Real required_a_value = mf_value(*required_a, point);
        const Real required_b_value = mf_value(*required_b, point);
        const Real transpose_a_value = mf_value(*transpose_a, point);
        const Real transpose_b_value = mf_value(*transpose_b, point);
        EXPECT_TRUE(is_changed(required_a_value))
            << "required component A, direction=" << dir
            << ", high=" << !face.isLow();
        EXPECT_TRUE(is_changed(required_b_value))
            << "required component B, direction=" << dir
            << ", high=" << !face.isLow();
        EXPECT_NEAR(required_a_value, transpose_a_value, Real(1.e-10));
        EXPECT_NEAR(required_b_value, transpose_b_value, Real(1.e-10));
    }
}

// Motivation: the water path must match the saturation oracle on the
// selected face and preserve the land sentinel on every other output box.
// This is the serial counterpart of the distributed qsurf ownership test.
TEST(SurfaceLayer, QsurfMatchesReferenceOnSelectedFace)
{
    ScopedSurfaceLayerParams params("unit_surface_layer_qsurf_serial");

    for (const auto& face : all_faces()) {
        SCOPED_TRACE(std::string("direction=") +
                     std::to_string(face.coordDir()) +
                     ", high=" + std::to_string(!face.isLow()));
        SurfaceLayerFields fields(make_qsurf_geometry(), true);
        fields.lmask[0]->setVal(0);
        auto layer = fields.prepare_layer(
            face, active_faces({face}), "unit_surface_layer_qsurf_serial", true);
        const MultiFab* qsurf = layer->get_q_surf(0);
        const Real expected = expected_qsat(fields.geom);
        Long selected_count = 0;
        for (MFIter mfi(*qsurf, false); mfi.isValid(); ++mfi) {
            const Box& source = fields.ba[mfi.index()];
            const int dir = face.coordDir();
            const bool selected = face.isLow()
                ? source.smallEnd(dir) == fields.domain.smallEnd(dir)
                : source.bigEnd(dir) == fields.domain.bigEnd(dir);
            const Box& valid = mfi.validbox();
            const auto qsurf_arr = qsurf->const_array(mfi);
            ReduceOps<ReduceOpMin, ReduceOpMax> reduce_op;
            ReduceData<Real, Real> reduce_data(reduce_op);
            reduce_op.eval(valid, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    -> GpuTuple<Real, Real>
                {
                    const Real value = qsurf_arr(i,j,k);
                    return {value, value};
                });
            Gpu::streamSynchronize();
            const auto qsurf_range = reduce_data.value();
            const Real qsurf_min = get<0>(qsurf_range);
            const Real qsurf_max = get<1>(qsurf_range);
            if (selected) {
                selected_count += static_cast<Long>(valid.numPts());
                ERF_EXPECT_NEAR(qsurf_min, expected,
                                qsat_tolerance(expected));
                ERF_EXPECT_NEAR(qsurf_max, expected,
                                qsat_tolerance(expected));
            } else {
                EXPECT_EQ(qsurf_min, tau_sentinel);
                EXPECT_EQ(qsurf_max, tau_sentinel);
            }
        }
        EXPECT_GT(selected_count, 0);
    }
}

// Motivation: tau31 and tau32 are optional away from their corresponding
// lateral faces.  The optional paths must not dereference absent outputs.
TEST(SurfaceLayer, OptionalTransposeStressCanBeAbsent)
{
    ScopedSurfaceLayerParams params("unit_surface_layer_missing_tau");

    const auto exercise = [] (const Orientation face,
                              const bool omit_tau31,
                              const bool omit_tau32) {
        SCOPED_TRACE(std::string("direction=") +
                     std::to_string(face.coordDir()) +
                     ", high=" + std::to_string(!face.isLow()));
        SurfaceLayerFields fields;
        auto layer = fields.prepare_layer(
            face, active_faces({face}), "unit_surface_layer_missing_tau");
        if (omit_tau31) { fields.tau[TauType::tau31].reset(); }
        if (omit_tau32) { fields.tau[TauType::tau32].reset(); }
        fields.impose(*layer);

        const IntVect point = face_point(fields.domain, face);
        const MultiFab* required_a = nullptr;
        const MultiFab* required_b = nullptr;
        if (face.coordDir() == 0) {
            required_a = fields.tau[TauType::tau21].get();
            required_b = fields.tau[TauType::tau31].get();
        } else if (face.coordDir() == 1) {
            required_a = fields.tau[TauType::tau12].get();
            required_b = fields.tau[TauType::tau32].get();
        } else {
            required_a = fields.tau[TauType::tau13].get();
            required_b = fields.tau[TauType::tau23].get();
        }
        ASSERT_NE(required_a, nullptr);
        ASSERT_NE(required_b, nullptr);
        EXPECT_TRUE(is_changed(mf_value(*required_a, point)));
        EXPECT_TRUE(is_changed(mf_value(*required_b, point)));
    };

    exercise(Orientation(Direction::x, Orientation::low), false, true);
    exercise(Orientation(Direction::y, Orientation::low), true, false);
    exercise(Orientation(Direction::z, Orientation::low), true, true);
}

// Motivation: when multiple surface-layer faces meet, transpose writes at
// shared corners must not overwrite the stress component owned by the
// neighboring face; both normal components must remain available.
TEST(SurfaceLayer, MoengPreservesRequiredStressAtMixedFaceCorners)
{
    ScopedSurfaceLayerParams params("unit_surface_layer_corners");
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
        const auto active = active_faces({pair[0], pair[1]});
        auto first = fields.prepare_layer(
            pair[0], active, "unit_surface_layer_corners");
        fields.impose(*first);

        const Orientation first_face = pair[0];
        const Orientation second_face = pair[1];
        const int first_dir = first_face.coordDir();
        const int second_dir = second_face.coordDir();
        const int corner_dir0 = std::min(first_dir, second_dir);
        const int corner_dir1 = std::max(first_dir, second_dir);
        IntVect corner(fields.domain.smallEnd());
        corner[corner_dir0] = first_dir == corner_dir0
            ? (first_face.isLow() ? fields.domain.smallEnd(corner_dir0)
                                  : fields.domain.bigEnd(corner_dir0) + 1)
            : (second_face.isLow() ? fields.domain.smallEnd(corner_dir0)
                                   : fields.domain.bigEnd(corner_dir0) + 1);
        corner[corner_dir1] = first_dir == corner_dir1
            ? (first_face.isLow() ? fields.domain.smallEnd(corner_dir1)
                                  : fields.domain.bigEnd(corner_dir1) + 1)
            : (second_face.isLow() ? fields.domain.smallEnd(corner_dir1)
                                   : fields.domain.bigEnd(corner_dir1) + 1);
        const int cell_dir = 3 - corner_dir0 - corner_dir1;
        corner[cell_dir] = fields.domain.smallEnd(cell_dir) + 1;

        const MultiFab* first_required = nullptr;
        const MultiFab* first_transpose = nullptr;
        const MultiFab* second_required = nullptr;
        const MultiFab* second_transpose = nullptr;
        if (first_dir == 0) {
            first_required = (second_dir == 1)
                ? fields.tau[TauType::tau21].get()
                : fields.tau[TauType::tau31].get();
            first_transpose = (second_dir == 1)
                ? fields.tau[TauType::tau12].get()
                : fields.tau[TauType::tau13].get();
        } else if (first_dir == 1) {
            first_required = (second_dir == 0)
                ? fields.tau[TauType::tau12].get()
                : fields.tau[TauType::tau32].get();
            first_transpose = (second_dir == 0)
                ? fields.tau[TauType::tau21].get()
                : fields.tau[TauType::tau23].get();
        } else {
            first_required = (second_dir == 0)
                ? fields.tau[TauType::tau13].get()
                : fields.tau[TauType::tau23].get();
            first_transpose = (second_dir == 0)
                ? fields.tau[TauType::tau31].get()
                : fields.tau[TauType::tau32].get();
        }

        if (second_dir == 0) {
            second_required = (first_dir == 1)
                ? fields.tau[TauType::tau21].get()
                : fields.tau[TauType::tau31].get();
            second_transpose = (first_dir == 1)
                ? fields.tau[TauType::tau12].get()
                : fields.tau[TauType::tau13].get();
        } else if (second_dir == 1) {
            second_required = (first_dir == 0)
                ? fields.tau[TauType::tau12].get()
                : fields.tau[TauType::tau32].get();
            second_transpose = (first_dir == 0)
                ? fields.tau[TauType::tau21].get()
                : fields.tau[TauType::tau23].get();
        } else {
            second_required = (first_dir == 0)
                ? fields.tau[TauType::tau13].get()
                : fields.tau[TauType::tau23].get();
            second_transpose = (first_dir == 0)
                ? fields.tau[TauType::tau31].get()
                : fields.tau[TauType::tau32].get();
        }

        const Real first_value = mf_value(*first_required, corner);
        const Real first_transpose_before = mf_value(*first_transpose, corner);
        EXPECT_TRUE(is_changed(first_value));
        EXPECT_EQ(first_transpose_before, tau_sentinel);

        auto second = fields.prepare_layer(
            second_face, active, "unit_surface_layer_corners");
        fields.impose(*second);

        const Real second_value = mf_value(*second_required, corner);
        const Real second_transpose_after = mf_value(*second_transpose, corner);
        EXPECT_TRUE(is_changed(second_value));
        EXPECT_NEAR(second_transpose_after, first_value, Real(1.e-10));
        EXPECT_NEAR(first_value, mf_value(*first_required, corner), Real(1.e-10));
    }
}
