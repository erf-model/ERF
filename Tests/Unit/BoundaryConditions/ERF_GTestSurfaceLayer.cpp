#include "AMReX_Gpu.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFab.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Reduce.H"

#include "ERF_SurfaceLayer.H"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <memory>
#include <string>

using namespace amrex;

namespace {

constexpr Real tau_sentinel = Real(-987654.0);

Geometry
make_geometry ()
{
    const Box domain(IntVect(2, 2, 2), IntVect(4, 4, 4));
    const RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                           {AMREX_D_DECL(3.0, 3.0, 3.0)});
    const Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
    return Geometry(domain, &real_box, 0, is_periodic.data());
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

    SurfaceLayerFields ()
        : geom(make_geometry()),
          domain(geom.Domain()),
          ba(domain),
          dm(ba),
          cons(ba, dm, 3, 2),
          xvel(BoxArray(surroundingNodes(domain, 0)), dm, 1, 2),
          yvel(BoxArray(surroundingNodes(domain, 1)), dm, 1, 2),
          zvel(BoxArray(surroundingNodes(domain, 2)), dm, 1, 2),
          theta(std::make_unique<MultiFab>(ba, dm, 1, 2)),
          xheat_flux(BoxArray(surroundingNodes(domain, 0)), dm, 1, 1),
          yheat_flux(BoxArray(surroundingNodes(domain, 1)), dm, 1, 1),
          zheat_flux(BoxArray(surroundingNodes(domain, 2)), dm, 1, 1)
    {
        cons.setVal(Real(1.0));
        xvel.setVal(Real(3.0));
        yvel.setVal(Real(4.0));
        zvel.setVal(Real(5.0));
        theta->setVal(Real(300.0));

        lmask.emplace_back(std::make_unique<iMultiFab>(
            collapse_z(ba), dm, 1, IntVect(2, 2, 0)));
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

    std::unique_ptr<SurfaceLayer>
    prepare_layer (const Orientation face,
                   const GpuArray<int, AMREX_SPACEDIM*2>& active_faces,
                   const std::string& prefix)
    {
        bool rotate = false;
        Vector<Geometry> geoms{geom};
        Vector<std::unique_ptr<MultiFab>> qv_prim(1);
        Vector<std::unique_ptr<MultiFab>> z_phys_nd(1);
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
        const moeng_flux flux(Real(0.1), face.isLow());
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
        EXPECT_LT(low_u[dir], Real(0.0)) << "direction=" << dir;
        EXPECT_LT(low_v[dir], Real(0.0)) << "direction=" << dir;
        EXPECT_NEAR(high_u[dir], -low_u[dir], Real(1.e-12))
            << "direction=" << dir;
        EXPECT_NEAR(high_v[dir], -low_v[dir], Real(1.e-12))
            << "direction=" << dir;
    }
}

// Motivation: each Cartesian face owns two normal momentum-stress components,
// while an isolated face must still populate the equivalent transpose
// components with the same computed Moeng flux.
TEST(SurfaceLayer, MoengWritesRequiredComponentsAndPreservesIsolatedSymmetry)
{
    ScopedSurfaceLayerParams params("unit_surface_layer_faces");
    const auto faces = all_faces();

    for (const auto face : faces) {
        SurfaceLayerFields fields;
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
