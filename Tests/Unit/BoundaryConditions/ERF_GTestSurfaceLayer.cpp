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
#include <cstdio>
#include <filesystem>
#include <initializer_list>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <utility>

using namespace amrex;
using erf_surface_layer_test::expected_qsat;
using erf_surface_layer_test::expected_surface_pressure;
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

struct ScopedTestFile
{
    explicit ScopedTestFile (std::filesystem::path path_in)
        : path(std::move(path_in)) {}

    ~ScopedTestFile ()
    {
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }

    std::filesystem::path path;
};

bool
is_changed (const Real value)
{
    return std::isfinite(value) && value != tau_sentinel;
}

struct MoengStressValues
{
    std::array<Real, AMREX_SPACEDIM> low_u{};
    std::array<Real, AMREX_SPACEDIM> high_u{};
    std::array<Real, AMREX_SPACEDIM> low_v{};
    std::array<Real, AMREX_SPACEDIM> high_v{};
};

MoengStressValues compute_moeng_stress_values ()
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

    MoengStressValues values;
    for (const auto& face : all_faces()) {
        const moeng_flux flux(Real(0.1), face.isLow(),
                              domain.smallEnd(2), domain.bigEnd(2));
        const int dir = face.coordDir();
        int i = domain.smallEnd(0) + 1;
        int j = domain.smallEnd(1) + 1;
        int k = domain.smallEnd(2) + 1;
        const int normal_index = face.isLow()
            ? domain.smallEnd(dir) : domain.bigEnd(dir) + 1;
        if (dir == 0) { i = normal_index; }
        else if (dir == 1) { j = normal_index; }
        else { k = normal_index; }

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
        if (face.isLow()) {
            values.low_u[dir] = stress_u;
            values.low_v[dir] = stress_v;
        } else {
            values.high_u[dir] = stress_u;
            values.high_v[dir] = stress_v;
        }
    }
    return values;
}

class ScopedSurfaceLayerParams
{
public:
    explicit ScopedSurfaceLayerParams (const char* prefix)
        : m_pp(prefix)
    {
        m_pp.add("most.average_policy", 0);
        m_pp.add("most.z0", Real(0.1));
        m_pp.add("most.surf_temp", Real(300.0));
    }

    ~ScopedSurfaceLayerParams ()
    {
        m_pp.remove("most.average_policy");
        m_pp.remove("most.z0");
        m_pp.remove("most.surf_temp");
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
    std::unique_ptr<MultiFab> lsm_tsurf;
    std::unique_ptr<MultiFab> coupled_sst;
    std::unique_ptr<iMultiFab> coupled_valid;

    Vector<std::unique_ptr<iMultiFab>> lmask;
    std::unique_ptr<MultiFab> no_walldist;

    Vector<std::unique_ptr<MultiFab>> tau;
    MultiFab xheat_flux;
    MultiFab yheat_flux;
    MultiFab zheat_flux;

    Vector<MultiFab*> state;

    explicit SurfaceLayerFields (const Geometry& geometry = make_geometry(),
                                 const bool split_lateral = false,
                                 const int ncons = RhoQ1_comp + 1)
        : geom(geometry),
          domain(geom.Domain()),
          ba(split_lateral ? make_qsurf_box_array(domain) : BoxArray(domain)),
          dm(ba),
          cons(ba, dm, ncons, test_state_ng),
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
        if (cons.nComp() > RhoQ1_comp) {
            cons.setVal(test_rho * test_qv, RhoQ1_comp, 1);
        }
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

    void set_surface_cell_pressure (const Real pressure, const Real qv = Real(0.0))
    {
        const Real rho_theta = (p_0 / R_d) *
            std::pow(pressure / p_0, Real(1.0) / Gamma) /
            (Real(1.0) + RvoRd * qv);
        for (MFIter mfi(cons, false); mfi.isValid(); ++mfi) {
            auto& fab = cons[mfi];
            fab.setVal<RunOn::Device>(test_rho, fab.box(), Rho_comp, 1);
            fab.setVal<RunOn::Device>(rho_theta, fab.box(), RhoTheta_comp, 1);
            if (fab.nComp() > RhoQ1_comp) {
                fab.setVal<RunOn::Device>(test_rho * qv, fab.box(), RhoQ1_comp, 1);
            }
        }
        Gpu::streamSynchronize();
    }

    std::unique_ptr<SurfaceLayer>
    prepare_layer (const Orientation face,
                   const GpuArray<int, AMREX_SPACEDIM*2>& active_faces,
                   const std::string& prefix,
                   const bool with_moisture = false,
                   const bool update_fluxes = true,
                   const std::string& lsm_name = "",
                   const Real lsm_value = Real(0.0),
                   const bool coupled_active = false,
                   const Real rdOcp = RdoCp)
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
            Vector<Vector<Real>>{},
            MeshType::ConstantDz, TerrainType::None, TurbChoice{}, rdOcp,
            0.0, 0.0);
        layer->set_surface_layer_faces(active_faces);
        if (coupled_active) {
            layer->set_coupled_sst_active(true);
            coupled_sst = std::make_unique<MultiFab>(
                collapse_z(ba), dm, 1,
                IntVect(AMREX_D_DECL(0, 0, 0)));
            coupled_sst->setVal(Real(290.0));
            coupled_valid = std::make_unique<iMultiFab>(
                collapse_z(ba), dm, 1,
                IntVect(AMREX_D_DECL(0, 0, 0)));
            coupled_valid->setVal(0);
        }

        std::unique_ptr<MultiFab> qr_prim;
        Vector<MultiFab*> empty_mfs;
        Vector<std::string> empty_names;
        Vector<MultiFab*> lsm_data;
        Vector<std::string> lsm_names;
        if (!lsm_name.empty()) {
            lsm_tsurf = std::make_unique<MultiFab>(
                collapse_z(ba), dm, 1,
                IntVect(AMREX_D_DECL(test_state_ng, test_state_ng, 0)));
            lsm_tsurf->setVal(lsm_value);
            lsm_data.push_back(lsm_tsurf.get());
            lsm_names.push_back(lsm_name);
        }
        Vector<std::unique_ptr<MultiFab>> sst;
        Vector<std::unique_ptr<MultiFab>> tsk;
        layer->make_SurfaceLayer_at_level(
            0, 1, state, theta, qv_prim[0], qr_prim, z_phys_nd[0],
            nullptr, nullptr, nullptr, lsm_data, lsm_names,
            empty_mfs, empty_names, sst, tsk, lmask);
        if (coupled_active) {
            layer->update_coupled_sst_ptr(0, coupled_sst.get(), coupled_valid.get());
        }
        if (with_moisture) {
            layer->get_q_surf(0)->setVal(tau_sentinel);
        }
        if (update_fluxes) {
            layer->update_fluxes(0, 0.0, 0.0, cons, z_phys_nd[0],
                                 no_walldist, 20);
        }
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

Long check_qsurf_values (const SurfaceLayerFields& fields,
                         const MultiFab& qsurf,
                         const Orientation face,
                         const Real expected)
{
    Long selected_count = 0;
    for (MFIter mfi(qsurf, false); mfi.isValid(); ++mfi) {
        const Box& source = fields.ba[mfi.index()];
        const int dir = face.coordDir();
        const bool selected = face.isLow()
            ? source.smallEnd(dir) == fields.domain.smallEnd(dir)
            : source.bigEnd(dir) == fields.domain.bigEnd(dir);
        const Box& valid = mfi.validbox();
        const auto qsurf_arr = qsurf.const_array(mfi);
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
            ERF_EXPECT_NEAR(qsurf_min, expected, qsat_tolerance(expected));
            ERF_EXPECT_NEAR(qsurf_max, expected, qsat_tolerance(expected));
        } else {
            EXPECT_EQ(qsurf_min, tau_sentinel);
            EXPECT_EQ(qsurf_max, tau_sentinel);
        }
    }
    return selected_count;
}

// nvcc rejects an extended __device__ lambda whose enclosing function has
// private access, and gtest generates TestBody() as a private member, so the
// device fills below live here rather than inside the TEST bodies.
void set_quadratic_node_heights (MultiFab& z_phys_nd)
{
    for (MFIter mfi(z_phys_nd, false); mfi.isValid(); ++mfi) {
        auto z_arr = z_phys_nd.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            z_arr(i,j,k) = static_cast<Real>(k * k);
        });
    }
    Gpu::streamSynchronize();
}

} // namespace

// Motivation: the Moeng stress functor has separate x-, y-, and z-wall
// interpolation paths.  Exercise each path at both wall orientations so the
// SurfaceLayer caller can rely on the normal high-face index being mapped back
// to the adjacent interior cell.
TEST(SurfaceLayer, MoengDirectionalFluxesAreFiniteOnEveryWall)
{
    const auto values = compute_moeng_stress_values();
    const auto& low_u = values.low_u;
    const auto& high_u = values.high_u;
    const auto& low_v = values.low_v;
    const auto& high_v = values.high_v;

    for (const auto& face : all_faces()) {
        const int dir = face.coordDir();
        EXPECT_TRUE(std::isfinite(face.isLow() ? low_u[dir] : high_u[dir]))
            << "direction=" << dir << ", high=" << !face.isLow();
        EXPECT_TRUE(std::isfinite(face.isLow() ? low_v[dir] : high_v[dir]))
            << "direction=" << dir << ", high=" << !face.isLow();
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
    const std::string prefix = "unit_surface_layer_qsurf_serial";
    ScopedSurfaceLayerParams params(prefix.c_str());
    ParmParse pp(prefix);
    pp.add("most.roughness_type_sea", std::string("constant"));

    for (const auto& face : all_faces()) {
        SCOPED_TRACE(std::string("direction=") +
                     std::to_string(face.coordDir()) +
                     ", high=" + std::to_string(!face.isLow()));
        SurfaceLayerFields fields(make_qsurf_geometry(), true);
        auto layer = fields.prepare_layer(
            face, active_faces({face}), "unit_surface_layer_qsurf_serial",
            true, false);
        fields.lmask[0]->setVal(0);
        const Real pressure = expected_surface_pressure(fields.geom, face);
        const Real surface_theta = test_surface_temperature *
            std::pow(p_0 / pressure, RdoCp);
        layer->get_t_surf(0)->setVal(surface_theta);
        std::unique_ptr<MultiFab> z_phys_nd;
        // The qsurf boundary fill intentionally visits vertical state ghosts;
        // initialize those cells so this test exercises qsat rather than an
        // unrelated invalid-density path.
        for (MFIter mfi(fields.cons, false); mfi.isValid(); ++mfi) {
            auto& fab = fields.cons[mfi];
            fab.setVal<RunOn::Device>(test_rho, fab.box(), Rho_comp, 1);
            fab.setVal<RunOn::Device>(test_rho_theta, fab.box(), RhoTheta_comp, 1);
            fab.setVal<RunOn::Device>(test_rho * test_qv, fab.box(), RhoQ1_comp, 1);
        }
        Gpu::streamSynchronize();
        layer->fill_qsurf_with_qsat(0, fields.cons, z_phys_nd);
        const MultiFab* qsurf = layer->get_q_surf(0);
        const Real expected = expected_qsat(fields.geom, face);
        const Long selected_count = check_qsurf_values(
            fields, *qsurf, face, expected);
        EXPECT_GT(selected_count, 0);
    }
    pp.remove("most.roughness_type_sea");
}

// Motivation: a terrain-following z-high boundary must use the local upper
// W-face height with a negative signed offset. This independent oracle makes
// the upper-face pressure differ materially from both p_cc and the old
// ground-relative Compute_Zrel_AtCellCenter path.
TEST(SurfaceLayer, QsurfUsesLocalSignedZHighFacePressure)
{
    const std::string prefix = "unit_surface_layer_qsurf_zhigh_geometry";
    ScopedSurfaceLayerParams params(prefix.c_str());
    ParmParse pp(prefix);
    pp.add("most.roughness_type_sea", std::string("constant"));
    const Geometry geom = make_qsurf_geometry();
    const Orientation face(Direction::z, Orientation::high);
    SurfaceLayerFields fields(geom, true);
    fields.lmask[0]->setVal(0);
    auto layer = fields.prepare_layer(
        face, active_faces({face}), "unit_surface_layer_qsurf_zhigh_geometry",
        true, false);

    BoxArray node_ba(fields.ba);
    node_ba.convert(IntVect::TheNodeVector());
    auto z_phys_nd = std::make_unique<MultiFab>(node_ba, fields.dm, 1, 0);
    set_quadratic_node_heights(*z_phys_nd);

    // For the three-cell fixture, z_cc(k=2)-z_upper_face(k=3) = 6.5-9 = -2.5.
    constexpr Real local_delta_z = Real(-2.5);
    const Real pressure = expected_surface_pressure(geom, face, local_delta_z);
    const Real surface_theta = test_surface_temperature *
        std::pow(p_0 / pressure, RdoCp);
    layer->get_t_surf(0)->setVal(surface_theta);
    layer->fill_qsurf_with_qsat(0, fields.cons, z_phys_nd);

    const Real expected = expected_qsat(geom, face, local_delta_z);
    EXPECT_GT(check_qsurf_values(
        fields, *layer->get_q_surf(0), face, expected), Long(0));
    pp.remove("most.roughness_type_sea");
}

// Motivation: lateral and tangential ghost cells are not authoritative
// physical surface state. An invalid halo must be ignored without raising the
// collective fatal flag, while the valid physical z-high slab still receives
// qsat values.
TEST(SurfaceLayer, QsurfInvalidHaloIsNonfatal)
{
    const std::string prefix = "unit_surface_layer_qsurf_invalid_halo";
    ScopedSurfaceLayerParams params(prefix.c_str());
    ParmParse pp(prefix);
    pp.add("most.roughness_type_sea", std::string("constant"));
    const Geometry geom = make_qsurf_geometry();
    const Orientation face(Direction::z, Orientation::high);
    SurfaceLayerFields fields(geom, true);
    fields.lmask[0]->setVal(0);
    auto layer = fields.prepare_layer(
        face, active_faces({face}), "unit_surface_layer_qsurf_invalid_halo",
        true, false);
    const Real pressure = expected_surface_pressure(geom, face);
    const Real surface_theta = test_surface_temperature *
        std::pow(p_0 / pressure, RdoCp);
    layer->get_t_surf(0)->setVal(surface_theta);

    bool injected = false;
    for (MFIter mfi(fields.cons, false); mfi.isValid() && !injected; ++mfi) {
        Box ghost = mfi.validbox();
        ghost.grow(1);
        const IntVect point(ghost.bigEnd(0), ghost.smallEnd(1), mfi.validbox().bigEnd(2));
        if (mfi.fabbox().contains(point) && !mfi.validbox().contains(point)) {
            fields.cons[mfi].setVal<RunOn::Device>(
                Real(0.0), Box(point, point), Rho_comp, 1);
            injected = true;
        }
    }
    ASSERT_TRUE(injected);
    layer->fill_qsurf_with_qsat(0, fields.cons, nullptr);

    const Real expected = expected_qsat(geom, face);
    EXPECT_GT(check_qsurf_values(
        fields, *layer->get_q_surf(0), face, expected), Long(0));
    pp.remove("most.roughness_type_sea");
}

// Motivation: lateral-face qsat uses a grown tangential halo, so an invalid
// state value just outside the physical face must be ignored rather than
// treated as an authoritative conversion failure. This direct x-low case
// complements the z-high halo regression and verifies that the physical slab
// still receives the independent p_cc pressure oracle.
TEST(SurfaceLayer, QsurfLateralInvalidHaloIsNonfatal)
{
    const std::string prefix = "unit_surface_layer_qsurf_lateral_invalid_halo";
    ScopedSurfaceLayerParams params(prefix.c_str());
    ParmParse pp(prefix);
    pp.add("most.roughness_type_sea", std::string("constant"));
    const Geometry geom = make_qsurf_geometry();
    const Orientation face(Direction::x, Orientation::low);
    SurfaceLayerFields fields(geom, true);
    fields.lmask[0]->setVal(0);
    auto layer = fields.prepare_layer(
        face, active_faces({face}), prefix, true, false);
    const Real pressure = expected_surface_pressure(geom, face);
    const Real surface_theta = test_surface_temperature *
        std::pow(p_0 / pressure, RdoCp);
    layer->get_t_surf(0)->setVal(surface_theta);

    bool injected = false;
    for (MFIter mfi(fields.cons, false); mfi.isValid() && !injected; ++mfi) {
        const Box& source = mfi.validbox();
        if (source.smallEnd(0) != geom.Domain().smallEnd(0)) { continue; }
        const IntVect point(source.smallEnd(0), source.smallEnd(1) - 1,
                            source.smallEnd(2));
        if (mfi.fabbox().contains(point) && !source.contains(point)) {
            fields.cons[mfi].setVal<RunOn::Device>(
                Real(0.0), Box(point, point), Rho_comp, 1);
            injected = true;
        }
    }
    ASSERT_TRUE(injected);

    layer->fill_qsurf_with_qsat(0, fields.cons, nullptr);
    EXPECT_GT(check_qsurf_values(
        fields, *layer->get_q_surf(0), face, expected_qsat(geom, face)), Long(0));
    pp.remove("most.roughness_type_sea");
}

// Motivation: MOST's prescribed surface value is already potential
// temperature. A pressure-dependent conversion at this producer would change
// the long-standing MOST flux contract.
TEST(SurfaceLayer, PrescribedMostSurfaceTemperatureRemainsPotentialTemperature)
{
    ScopedSurfaceLayerParams params("unit_surface_layer_most_theta");
    SurfaceLayerFields fields;
    const Orientation face(Direction::z, Orientation::low);
    auto layer = fields.prepare_layer(
        face, active_faces({face}), "unit_surface_layer_most_theta");

    const IntVect point = face_point(fields.domain, face);
    EXPECT_EQ(mf_value(*layer->get_t_surf(0), point), test_surface_temperature);
}

// Motivation: the prescribed MOST heating rate is a tendency of the stored
// potential temperature, so it must not be Exner-converted either.
TEST(SurfaceLayer, PrescribedMostHeatingRateRemainsPotentialTemperatureTendency)
{
    const std::string prefix = "unit_surface_layer_most_heating";
    ScopedSurfaceLayerParams params(prefix.c_str());
    ParmParse pp(prefix);
    pp.add("most.surf_heating_rate", Real(3600.0));

    SurfaceLayerFields fields;
    const Orientation face(Direction::z, Orientation::low);
    auto layer = fields.prepare_layer(face, active_faces({face}), prefix);
    layer->update_fluxes(0, 2.0, 2.0, fields.cons, nullptr,
                         fields.no_walldist, 20);

    const IntVect point = face_point(fields.domain, face);
    EXPECT_EQ(mf_value(*layer->get_t_surf(0), point), test_surface_temperature + Real(2.0));
    pp.remove("most.surf_heating_rate");
}

// Motivation: Noah-MP's t_sfc is an absolute radiative temperature owned by
// the radiation path, not a SurfaceLayer theta-like input. It must not be
// adopted as MOST's thermal boundary; the configured MOST theta remains the
// fallback when no recognized theta-like LSM field is present.
TEST(SurfaceLayer, NoahRadiativeSurfaceTemperatureIsNotAdoptedBySurfaceLayer)
{
    const Orientation face(Direction::z, Orientation::low);
    ScopedSurfaceLayerParams params("unit_surface_layer_noah_t_sfc");
    SurfaceLayerFields fields;
    auto layer = fields.prepare_layer(
        face, active_faces({face}), "unit_surface_layer_noah_t_sfc",
        false, true, "t_sfc", Real(290.0));

    const IntVect point = face_point(fields.domain, face);
    EXPECT_EQ(mf_value(*layer->get_t_surf(0), point), test_surface_temperature);
}

// Motivation: SLM exposes its surface field as theta. The SurfaceLayer must
// preserve that canonical value, rather than treating every LSM temperature
// field as an absolute-temperature Noah-MP field.
TEST(SurfaceLayer, PotentialTemperatureLsmFieldRemainsUnchanged)
{
    const Geometry geom = make_qsurf_geometry();
    const Orientation face(Direction::z, Orientation::low);
    ScopedSurfaceLayerParams params("unit_surface_layer_lsm_theta");
    SurfaceLayerFields fields(geom);
    fields.set_surface_cell_pressure(Real(0.9) * p_0);
    auto layer = fields.prepare_layer(
        face, active_faces({face}), "unit_surface_layer_lsm_theta",
        false, true, "theta", Real(275.0));

    const IntVect point = face_point(fields.domain, face);
    EXPECT_EQ(mf_value(*layer->get_t_surf(0), point), Real(275.0));
}

// Motivation: the text-file SST contract is absolute temperature, while the
// SurfaceLayer field consumed by MOST is theta. This exercises the actual
// z-low forcing path at reduced pressure rather than only testing the EOS
// helper in isolation.
TEST(SurfaceLayer, TextSstIsConvertedToPotentialTemperature)
{
    const std::string prefix = "unit_surface_layer_text_sst";
    const auto file = std::filesystem::current_path() /
        ("erf_surface_layer_text_sst_" + std::to_string(sizeof(Real)) + ".txt");
    ScopedTestFile cleanup(file);
    {
        std::ofstream out(file);
        ASSERT_TRUE(out.good());
        out << "day sst(K)\n0.0 290.0\n1.0 290.0\n";
    }

    ScopedSurfaceLayerParams params(prefix.c_str());
    ParmParse pp(prefix);
    pp.add("most.use_sfc_sst", true);
    pp.add("most.sfc_file", file.string());

    const Geometry geom = make_qsurf_geometry();
    const Orientation face(Direction::z, Orientation::low);
    SurfaceLayerFields fields(geom);
    fields.lmask[0]->setVal(0);
    const Real pressure = Real(0.9) * p_0;
    fields.set_surface_cell_pressure(
        pressure - test_rho * CONST_GRAV * myhalf * geom.CellSize(2));
    auto layer = fields.prepare_layer(face, active_faces({face}), prefix);

    const Real expected_theta = Real(290.0) * std::pow(p_0 / pressure, RdoCp);
    const IntVect point = face_point(fields.domain, face);
    EXPECT_NEAR(mf_value(*layer->get_t_surf(0), point), expected_theta,
                Real(64.0) * std::numeric_limits<Real>::epsilon() * expected_theta);

    pp.remove("most.use_sfc_sst");
    pp.remove("most.sfc_file");
}

// Motivation: production dry conserved state has four components and no
// moisture density. The text-SST conversion must therefore obtain qv=0
// without reading the out-of-range RhoQ1_comp slot.
TEST(SurfaceLayer, DryTextSstDoesNotReadMoistureComponent)
{
    const std::string prefix = "unit_surface_layer_dry_text_sst";
    const auto file = std::filesystem::current_path() /
        ("erf_surface_layer_dry_text_sst_" + std::to_string(sizeof(Real)) + ".txt");
    ScopedTestFile cleanup(file);
    {
        std::ofstream out(file);
        ASSERT_TRUE(out.good());
        out << "day sst(K)\n0.0 290.0\n1.0 290.0\n";
    }

    ScopedSurfaceLayerParams params(prefix.c_str());
    ParmParse pp(prefix);
    pp.add("most.use_sfc_sst", true);
    pp.add("most.sfc_file", file.string());

    const Geometry geom = make_qsurf_geometry();
    const Orientation face(Direction::z, Orientation::low);
    SurfaceLayerFields fields(geom, false, RhoQ1_comp);
    ASSERT_EQ(fields.cons.nComp(), RhoQ1_comp);
    fields.lmask[0]->setVal(0);
    const Real pressure = Real(0.9) * p_0;
    fields.set_surface_cell_pressure(
        pressure - test_rho * CONST_GRAV * myhalf * geom.CellSize(2));
    auto layer = fields.prepare_layer(face, active_faces({face}), prefix);

    const Real expected_theta = Real(290.0) * std::pow(p_0 / pressure, RdoCp);
    const IntVect point = face_point(fields.domain, face);
    EXPECT_NEAR(mf_value(*layer->get_t_surf(0), point), expected_theta,
                Real(64.0) * std::numeric_limits<Real>::epsilon() * expected_theta);

    pp.remove("most.use_sfc_sst");
    pp.remove("most.sfc_file");
}

// Motivation: coupled SST is an absolute-temperature producer with partial
// water coverage. Only covered cells should be converted into theta; an
// uncovered water cell must retain the existing SurfaceLayer fallback.
TEST(SurfaceLayer, CoupledSstConvertsOnlyCoveredWaterCells)
{
    const std::string prefix = "unit_surface_layer_coupled_sst";
    ScopedSurfaceLayerParams params(prefix.c_str());
    const Geometry geom = make_qsurf_geometry();
    const Orientation face(Direction::z, Orientation::low);
    SurfaceLayerFields fields(geom);
    fields.lmask[0]->setVal(0);
    const Real pressure = Real(0.9) * p_0;
    fields.set_surface_cell_pressure(
        pressure - test_rho * CONST_GRAV * myhalf * geom.CellSize(2));
    auto layer = fields.prepare_layer(
        face, active_faces({face}), prefix, false, false, "", Real(0.0), true);

    const IntVect covered = face_point(fields.domain, face);
    fields.coupled_valid->setVal(1, Box(covered, covered), 0, 1);
    layer->update_fluxes(0, 0.0, 0.0, fields.cons, nullptr,
                         fields.no_walldist, 20);

    const Real expected_theta = Real(290.0) * std::pow(p_0 / pressure, RdoCp);
    EXPECT_NEAR(mf_value(*layer->get_t_surf(0), covered), expected_theta,
                Real(64.0) * std::numeric_limits<Real>::epsilon() * expected_theta);
    const IntVect uncovered(0, 0, 0);
    EXPECT_EQ(mf_value(*layer->get_t_surf(0), uncovered), test_surface_temperature);
}

// Motivation: production dry coupled-SST state has no RhoQ1_comp component.
// The covered conversion must use qv=0 without reading beyond the four
// conserved components while still applying the pressure-dependent oracle.
TEST(SurfaceLayer, DryCoupledSstDoesNotReadMoistureComponent)
{
    const std::string prefix = "unit_surface_layer_dry_coupled_sst";
    ScopedSurfaceLayerParams params(prefix.c_str());
    const Geometry geom = make_qsurf_geometry();
    const Orientation face(Direction::z, Orientation::low);
    SurfaceLayerFields fields(geom, false, RhoQ1_comp);
    ASSERT_EQ(fields.cons.nComp(), RhoQ1_comp);
    fields.lmask[0]->setVal(0);
    const Real pressure = Real(0.9) * p_0;
    fields.set_surface_cell_pressure(
        pressure - test_rho * CONST_GRAV * myhalf * geom.CellSize(2));
    auto layer = fields.prepare_layer(
        face, active_faces({face}), prefix, false, false, "", Real(0.0), true);
    fields.coupled_valid->setVal(1);
    layer->update_fluxes(0, 0.0, 0.0, fields.cons, nullptr,
                         fields.no_walldist, 20);

    const Real expected_theta = Real(290.0) * std::pow(p_0 / pressure, RdoCp);
    const IntVect point = face_point(fields.domain, face);
    EXPECT_NEAR(mf_value(*layer->get_t_surf(0), point), expected_theta,
                Real(64.0) * std::numeric_limits<Real>::epsilon() * expected_theta);
}

// Motivation: production coupled-SST donors have no lateral ghost cells, but
// the SurfaceLayer destination carries a one-cell grown halo. A covered
// nonperiodic edge must therefore sample the clamped physical donor rather
// than be dropped when the target is grown. Repeating the update without
// coverage also verifies that the edge retains the existing fallback.
TEST(SurfaceLayer, CoupledSstZeroGhostDonorClampsNonperiodicEdge)
{
    const std::string prefix = "unit_surface_layer_coupled_sst_zero_ghost_edge";
    ScopedSurfaceLayerParams params(prefix.c_str());
    const Geometry geom = make_qsurf_geometry();
    const Orientation face(Direction::z, Orientation::low);
    SurfaceLayerFields fields(geom);
    fields.lmask[0]->setVal(0);
    const Real pressure = Real(0.9) * p_0;
    fields.set_surface_cell_pressure(
        pressure - test_rho * CONST_GRAV * myhalf * geom.CellSize(2));
    auto layer = fields.prepare_layer(
        face, active_faces({face}), prefix, false, false, "", Real(0.0), true);
    fields.coupled_valid->setVal(1);
    layer->update_fluxes(0, 0.0, 0.0, fields.cons, nullptr,
                         fields.no_walldist, 20);

    const Real expected_theta = Real(290.0) * std::pow(p_0 / pressure, RdoCp);
    const MultiFab* t_surf = layer->get_t_surf(0);
    bool checked_edge = false;
    for (int ibox = 0; ibox < fields.ba.size(); ++ibox) {
        const Box& source = fields.ba[ibox];
        if (source.smallEnd(0) != geom.Domain().smallEnd(0)) { continue; }
        const Box target = t_surf->boxArray()[ibox];
        const IntVect interior = target.smallEnd();
        IntVect edge = interior;
        edge[0] -= 1;
        ASSERT_TRUE((*t_surf)[ibox].box().contains(edge));
        const auto t_arr = (*t_surf)[ibox].const_array();
        EXPECT_NEAR(t_arr(interior[0], interior[1], interior[2]),
                    expected_theta, Real(64.0) * std::numeric_limits<Real>::epsilon() *
                    expected_theta);
        EXPECT_NEAR(t_arr(edge[0], edge[1], edge[2]),
                    expected_theta, Real(64.0) * std::numeric_limits<Real>::epsilon() *
                    expected_theta);
        checked_edge = true;
    }
    ASSERT_TRUE(checked_edge);

    fields.coupled_valid->setVal(0);
    layer->get_t_surf(0)->setVal(test_surface_temperature);
    layer->update_fluxes(0, 0.0, 0.0, fields.cons, nullptr,
                         fields.no_walldist, 20);
    for (int ibox = 0; ibox < fields.ba.size(); ++ibox) {
        const Box& source = fields.ba[ibox];
        if (source.smallEnd(0) != geom.Domain().smallEnd(0)) { continue; }
        const Box target = t_surf->boxArray()[ibox];
        const IntVect edge(target.smallEnd(0) - 1,
                           target.smallEnd(1), target.smallEnd(2));
        const auto t_arr = (*t_surf)[ibox].const_array();
        const Real fallback = t_arr(edge[0], edge[1], edge[2]);
        EXPECT_EQ(fallback, test_surface_temperature);
    }
}

// Motivation: a coupled SST pointer without a coverage mask carries no
// evidence that any cell has an ocean donor. The SurfaceLayer fallback must
// therefore remain untouched instead of treating the null mask as all-valid.
TEST(SurfaceLayer, CoupledSstWithoutCoverageMaskLeavesFallbackUnchanged)
{
    const std::string prefix = "unit_surface_layer_coupled_sst_no_mask";
    ScopedSurfaceLayerParams params(prefix.c_str());
    const Geometry geom = make_qsurf_geometry();
    const Orientation face(Direction::z, Orientation::low);
    SurfaceLayerFields fields(geom);
    fields.lmask[0]->setVal(0);
    auto layer = fields.prepare_layer(
        face, active_faces({face}), prefix, false, false, "", Real(0.0), true);

    layer->update_coupled_sst_ptr(0, fields.coupled_sst.get(), nullptr);
    layer->update_fluxes(0, 0.0, 0.0, fields.cons, nullptr,
                         fields.no_walldist, 20);

    const IntVect point = face_point(fields.domain, face);
    EXPECT_EQ(mf_value(*layer->get_t_surf(0), point), test_surface_temperature);
}

// Motivation: the atmospheric thermodynamic exponent is solver configuration,
// so the actual text-SST producer path must use a non-default cp rather than
// silently reverting to the hard-coded dry-air exponent.
TEST(SurfaceLayer, ConfiguredRdOcpControlsTextSstConversion)
{
    const std::string prefix = "unit_surface_layer_text_sst_custom_rdOcp";
    const auto file = std::filesystem::current_path() /
        ("erf_surface_layer_text_sst_custom_rdOcp_" + std::to_string(sizeof(Real)) + ".txt");
    ScopedTestFile cleanup(file);
    {
        std::ofstream out(file);
        ASSERT_TRUE(out.good());
        out << "day sst(K)\n0.0 290.0\n1.0 290.0\n";
    }

    ScopedSurfaceLayerParams params(prefix.c_str());
    ParmParse pp(prefix);
    pp.add("most.use_sfc_sst", true);
    pp.add("most.sfc_file", file.string());

    const Geometry geom = make_qsurf_geometry();
    const Orientation face(Direction::z, Orientation::low);
    SurfaceLayerFields fields(geom);
    fields.lmask[0]->setVal(0);
    const Real pressure = Real(0.9) * p_0;
    fields.set_surface_cell_pressure(
        pressure - test_rho * CONST_GRAV * myhalf * geom.CellSize(2));
    const Real rdOcp = R_d / Real(1100.0);
    auto layer = fields.prepare_layer(
        face, active_faces({face}), prefix, false, true, "", Real(0.0), false, rdOcp);

    const Real expected_theta = Real(290.0) * std::pow(pressure / p_0, -rdOcp);
    const IntVect point = face_point(fields.domain, face);
    EXPECT_NEAR(mf_value(*layer->get_t_surf(0), point), expected_theta,
                Real(64.0) * std::numeric_limits<Real>::epsilon() * expected_theta);
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
