#ifndef ERF_GTEST_SURFACE_LAYER_PARALLEL_COMMON_H_
#define ERF_GTEST_SURFACE_LAYER_PARALLEL_COMMON_H_

#include "ERF_GTestSurfaceLayerCommon.H"

#include <AMReX_FabArrayBase.H>
#include <AMReX_Gpu.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>
#include <AMReX_iMultiFab.H>

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace erf_surface_layer_test {

inline Geometry
make_parallel_geometry ()
{
    const Box domain(IntVect(0), IntVect(31, 31, 3));
    const RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                           {AMREX_D_DECL(32.0, 32.0, 4.0)});
    const Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
    return Geometry(domain, &real_box, 0, is_periodic.data());
}

inline BoxArray
make_parallel_box_array (const Box& domain, const bool single_box = false)
{
    BoxArray ba(domain);
    if (single_box) { return ba; }
    // Lateral boxes span the full vertical domain, as required by the x/y
    // surface-layer kernels, while still distributing the face ownership.
    ba.maxSize(IntVect(AMREX_D_DECL(16, 16, 1024)));
    return ba;
}

inline DistributionMapping
make_parallel_distribution_mapping (const BoxArray& ba)
{
    Vector<int> pmap(ba.size(), 0);
    const int nprocs = ParallelDescriptor::NProcs();
    if (nprocs > 1) {
        for (int ibox = 0; ibox < ba.size(); ++ibox) {
            if (nprocs == 2) {
                // The two lateral neighbors of every selected face must be
                // assigned to different ranks.  This checkerboard mapping
                // achieves that for both x and y faces.
                pmap[ibox] = (ibox == 0 || ibox == 3) ? 0 : 1;
            } else {
                // With at least three ranks, assign the four-box layout to
                // distinct ranks where possible.  This also works when more
                // ranks are present than boxes; unused ranks are harmless.
                pmap[ibox] = ibox % nprocs;
            }
        }
    }
    return DistributionMapping(std::move(pmap));
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

class ScopedMFIterTileSize
{
public:
    explicit ScopedMFIterTileSize (const IntVect& tile_size)
        : m_saved(FabArrayBase::mfiter_tile_size)
    {
        FabArrayBase::mfiter_tile_size = tile_size;
    }

    ~ScopedMFIterTileSize ()
    {
        FabArrayBase::mfiter_tile_size = m_saved;
    }

private:
    IntVect m_saved;
};

inline BoxArray
collapse_z (const BoxArray& ba)
{
    BoxList boxes = ba.boxList();
    for (auto& box : boxes) {
        box.setRange(2, 0);
    }
    return BoxArray(std::move(boxes));
}

constexpr Real u_star_seed = Real(710000.0);
constexpr Real t_star_seed = Real(720000.0);
constexpr Real q_star_seed = Real(730000.0);
constexpr Real olen_seed = Real(740000.0);

inline Real
seeded_value (const Real base, const int box_index)
{
    return base + static_cast<Real>(box_index);
}

inline Real
halo_tolerance (const Real expected)
{
    // Halo exchange copies the same computed value; allow a small
    // backend-independent roundoff margin while retaining a scaled oracle.
    return Real(256.0) * std::numeric_limits<Real>::epsilon() *
           std::max(Real(1.0), std::abs(expected));
}

inline void
set_fab_values (MultiFab& mf, const Real base)
{
    for (MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const Real value = seeded_value(base, mfi.index());
        auto array = mf.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            array(i,j,k) = value;
        });
    }
    Gpu::streamSynchronize();
}

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
    Vector<std::unique_ptr<MultiFab>> qv_prim;
    Vector<std::unique_ptr<iMultiFab>> lmask;
    std::unique_ptr<MultiFab> no_walldist;
    Vector<std::unique_ptr<MultiFab>> tau;
    MultiFab xheat_flux;
    MultiFab yheat_flux;
    MultiFab zheat_flux;
    MultiFab xqv_flux;
    MultiFab yqv_flux;
    MultiFab zqv_flux;
    Vector<MultiFab*> state;

    explicit SurfaceLayerFields (const bool single_box = false)
        : geom(make_parallel_geometry()),
          domain(geom.Domain()),
          ba(make_parallel_box_array(domain, single_box)),
          dm(make_parallel_distribution_mapping(ba)),
          cons(ba, dm, RhoQ1_comp + 1,
               IntVect(AMREX_D_DECL(test_state_ng, test_state_ng, 1))),
          xvel(make_face_box_array(ba, 0), dm, 1, test_velocity_ng),
          yvel(make_face_box_array(ba, 1), dm, 1, test_velocity_ng),
          zvel(make_face_box_array(ba, 2), dm, 1, test_velocity_ng),
          theta(std::make_unique<MultiFab>(
              ba, dm, 1,
              IntVect(AMREX_D_DECL(test_state_ng, test_state_ng,
                                   test_primitive_z_ng)))),
          qv_prim(1),
          xheat_flux(make_face_box_array(ba, 0), dm, 1, 1),
          yheat_flux(make_face_box_array(ba, 1), dm, 1, 1),
          zheat_flux(make_face_box_array(ba, 2), dm, 1, 1),
          xqv_flux(make_face_box_array(ba, 0), dm, 1, 1),
          yqv_flux(make_face_box_array(ba, 1), dm, 1, 1),
          zqv_flux(make_face_box_array(ba, 2), dm, 1, 1)
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
        BoxArray ba12(ba); ba12.convert(IntVect(1, 1, 0));
        BoxArray ba13(ba); ba13.convert(IntVect(1, 0, 1));
        BoxArray ba23(ba); ba23.convert(IntVect(0, 1, 1));
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
        for (auto& stress : tau) { stress->setVal(tau_sentinel); }
        xheat_flux.setVal(tau_sentinel);
        yheat_flux.setVal(tau_sentinel);
        zheat_flux.setVal(tau_sentinel);
        xqv_flux.setVal(tau_sentinel);
        yqv_flux.setVal(tau_sentinel);
        zqv_flux.setVal(tau_sentinel);
        Gpu::streamSynchronize();
    }

    void set_varying_surface_temperature ()
    {
        for (MFIter mfi(*theta, false); mfi.isValid(); ++mfi) {
            const Box box = mfi.fabbox();
            auto theta_arr = theta->array(mfi);
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                theta_arr(i,j,k) = test_surface_temperature +
                    Real(0.25) * static_cast<Real>(i) +
                    Real(0.5) * static_cast<Real>(j);
            });
        }
        Gpu::streamSynchronize();
    }

    void set_varying_surface_humidity ()
    {
        if (!qv_prim[0]) { return; }
        for (MFIter mfi(*qv_prim[0], false); mfi.isValid(); ++mfi) {
            const Real box_offset = Real(0.001) * static_cast<Real>(mfi.index());
            auto qv_arr = qv_prim[0]->array(mfi);
            ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                qv_arr(i,j,k) = test_qv + box_offset;
            });
        }
        Gpu::streamSynchronize();
    }

    std::unique_ptr<SurfaceLayer>
    prepare_layer (const Orientation face,
                   const GpuArray<int, AMREX_SPACEDIM*2>& active,
                   const std::string& prefix,
                   const bool with_moisture = false,
                   const bool varying_inputs = false,
                   const bool seed_outputs = false)
    {
        Vector<Geometry> geoms{geom};
        Vector<std::unique_ptr<MultiFab>> z_phys_nd(1);
        qv_prim[0].reset();
        if (with_moisture) {
            qv_prim[0] = std::make_unique<MultiFab>(
                ba, dm, 1,
                IntVect(AMREX_D_DECL(test_state_ng, test_state_ng,
                                     test_primitive_z_ng)));
            qv_prim[0]->setVal(test_qv);
        }
        if (varying_inputs) {
            set_varying_surface_temperature();
            set_varying_surface_humidity();
        }
        bool rotate = false;
        auto layer = std::make_unique<SurfaceLayer>(
            face, geoms, rotate, prefix, qv_prim, z_phys_nd,
            MeshType::ConstantDz, TerrainType::None, TurbChoice{}, 0.0, 0.0);
        layer->set_surface_layer_faces(active);

        std::unique_ptr<MultiFab> qr_prim;
        Vector<MultiFab*> empty_mfs;
        Vector<std::string> empty_names;
        Vector<std::unique_ptr<MultiFab>> sst;
        Vector<std::unique_ptr<MultiFab>> tsk;
        layer->make_SurfaceLayer_at_level(
            0, 1, state, theta, qv_prim[0], qr_prim, z_phys_nd[0],
            nullptr, nullptr, nullptr, empty_mfs, empty_names,
            empty_mfs, empty_names, sst, tsk, lmask);
        if (seed_outputs) {
            set_fab_values(*layer->get_u_star(0), u_star_seed);
            set_fab_values(*layer->get_t_star(0), t_star_seed);
            set_fab_values(*layer->get_q_star(0), q_star_seed);
            set_fab_values(*layer->get_olen(0), olen_seed);
        }
        if (with_moisture) { layer->get_q_surf(0)->setVal(tau_sentinel); }
        layer->update_fluxes(0, 0.0, 0.0, cons, z_phys_nd[0],
                             no_walldist, 20);
        return layer;
    }

    void impose (SurfaceLayer& layer)
    {
        Vector<const MultiFab*> const_state{&cons, &xvel, &yvel, &zvel};
        layer.impose_SurfaceLayer_bcs(
            0, const_state, tau, &xheat_flux, &yheat_flux, &zheat_flux,
            &xqv_flux, &yqv_flux, &zqv_flux, nullptr);
        Gpu::streamSynchronize();
    }

private:
    static BoxArray make_face_box_array (const BoxArray& cell_ba, const int dir)
    {
        BoxArray face_ba(cell_ba);
        face_ba.surroundingNodes(dir);
        return face_ba;
    }
};

struct GlobalRange
{
    Real lo{0};
    Real hi{0};
    Long count{0};
    Long total_count{0};
    Long sentinel_count{0};
};

inline Real
mf_value (const MultiFab& mf, const IntVect& point)
{
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    const Box point_box(point, point, mf.boxArray().ixType());
    for (MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const Box overlap = point_box & mfi.validbox();
        if (overlap.isEmpty()) { continue; }
        const auto array = mf.const_array(mfi);
        reduce_op.eval(overlap, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> GpuTuple<Real>
            { return {array(i,j,k)}; });
    }
    Gpu::streamSynchronize();
    Real value = get<0>(reduce_data.value());
    ParallelDescriptor::ReduceRealSum(value);
    return value;
}

inline Real
global_fab_value (const MultiFab& mf,
                  const int fab_index,
                  const IntVect& point,
                  const bool include_ghost)
{
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    const Box point_box(point, point, mf.boxArray().ixType());
    for (MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        if (mfi.index() != fab_index) { continue; }
        const Box& available = include_ghost ? mfi.fabbox() : mfi.validbox();
        const Box overlap = point_box & available;
        if (overlap.isEmpty()) { continue; }
        const auto array = mf.const_array(mfi);
        reduce_op.eval(overlap, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> GpuTuple<Real>
            { return {array(i,j,k)}; });
    }
    Gpu::streamSynchronize();
    Real value = get<0>(reduce_data.value());
    ParallelDescriptor::ReduceRealSum(value);
    return value;
}

struct TangentialHaloPair
{
    int target_fab;
    int neighbor_fab;
    IntVect point;
};

inline std::vector<TangentialHaloPair>
tangential_halo_pairs (const MultiFab& mf,
                       const BoxArray& source_ba,
                       const Box& domain,
                       const Orientation face)
{
    std::vector<TangentialHaloPair> result;
    const int dir = face.coordDir();
    if (dir == 2) { return result; }
    const int tangent_dir = 1 - dir;
    const BoxArray& target_ba = mf.boxArray();

    for (int ibox = 0; ibox < target_ba.size(); ++ibox) {
        const Box& source = source_ba[ibox];
        const bool selected = face.isLow()
            ? source.smallEnd(dir) == domain.smallEnd(dir)
            : source.bigEnd(dir) == domain.bigEnd(dir);
        if (!selected) { continue; }

        const Box& target = target_ba[ibox];
        if (target.smallEnd(tangent_dir) > domain.smallEnd(tangent_dir)) {
            IntVect point = target.smallEnd();
            point[tangent_dir] -= 1;
            for (int jbox = 0; jbox < target_ba.size(); ++jbox) {
                const Box& neighbor_source = source_ba[jbox];
                const bool neighbor_selected = face.isLow()
                    ? neighbor_source.smallEnd(dir) == domain.smallEnd(dir)
                    : neighbor_source.bigEnd(dir) == domain.bigEnd(dir);
                if (jbox != ibox && neighbor_selected &&
                    target_ba[jbox].contains(point)) {
                    result.push_back({ibox, jbox, point});
                    break;
                }
            }
        }
        if (target.bigEnd(tangent_dir) < domain.bigEnd(tangent_dir)) {
            IntVect point = target.bigEnd();
            point[tangent_dir] += 1;
            for (int jbox = 0; jbox < target_ba.size(); ++jbox) {
                const Box& neighbor_source = source_ba[jbox];
                const bool neighbor_selected = face.isLow()
                    ? neighbor_source.smallEnd(dir) == domain.smallEnd(dir)
                    : neighbor_source.bigEnd(dir) == domain.bigEnd(dir);
                if (jbox != ibox && neighbor_selected &&
                    target_ba[jbox].contains(point)) {
                    result.push_back({ibox, jbox, point});
                    break;
                }
            }
        }
    }
    return result;
}

inline GlobalRange
face_range (const MultiFab& mf, const Box& domain, const Orientation face)
{
    Box plane(domain);
    plane.convert(mf.boxArray().ixType());
    const int dir = face.coordDir();
    const int index = face.isLow() ? plane.smallEnd(dir) : plane.bigEnd(dir);
    plane.setSmall(dir, index);
    plane.setBig(dir, index);

    Long local_count = 0;
    Long local_total_count = 0;
    Long local_sentinel_count = 0;
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Long> reduce_data(reduce_op);
    for (MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const Box& valid = mfi.validbox();
        local_total_count += static_cast<Long>(valid.numPts());
        const auto array = mf.const_array(mfi);
        reduce_op.eval(valid, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> GpuTuple<Long>
            { return {array(i,j,k) == tau_sentinel ? Long(1) : Long(0)}; });
        local_count += static_cast<Long>((plane & valid).numPts());
    }

    Real local_lo = std::numeric_limits<Real>::max();
    Real local_hi = std::numeric_limits<Real>::lowest();
    if (local_count > 0) {
        local_lo = mf.min(plane, 0, 0, true);
        local_hi = mf.max(plane, 0, 0, true);
    }
    Gpu::streamSynchronize();
    local_sentinel_count = get<0>(reduce_data.value());
    ParallelDescriptor::ReduceLongSum(local_count);
    ParallelDescriptor::ReduceLongSum(local_total_count);
    ParallelDescriptor::ReduceLongSum(local_sentinel_count);
    ParallelDescriptor::ReduceRealMin(local_lo);
    ParallelDescriptor::ReduceRealMax(local_hi);
    return {local_lo, local_hi, local_count, local_total_count,
            local_sentinel_count};
}

// The lateral surface fields collapse the face-normal coordinate.  Therefore
// a geometric plane intersects non-face source boxes as well; use the
// original source ownership to identify the selected FABs.
inline GlobalRange
face_owned_range (const MultiFab& mf,
                  const BoxArray& source_ba,
                  const Box& domain,
                  const Orientation face)
{
    const int dir = face.coordDir();
    const int source_face_index = face.isLow()
        ? domain.smallEnd(dir) : domain.bigEnd(dir);

    Box face_plane(domain);
    face_plane.convert(mf.boxArray().ixType());
    const int field_face_index = face.isLow()
        ? face_plane.smallEnd(dir) : face_plane.bigEnd(dir);
    face_plane.setRange(dir, field_face_index);

    Long local_count = 0;
    Long local_total_count = 0;
    Long local_sentinel_count = 0;
    Real local_lo = std::numeric_limits<Real>::max();
    Real local_hi = std::numeric_limits<Real>::lowest();

    ReduceOps<ReduceOpMin, ReduceOpMax> range_op;
    ReduceData<Real, Real> range_data(range_op);
    ReduceOps<ReduceOpSum> sentinel_op;
    ReduceData<Long> sentinel_data(sentinel_op);

    for (MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const Box& valid = mfi.validbox();
        const auto array = mf.const_array(mfi);
        local_total_count += static_cast<Long>(valid.numPts());
        sentinel_op.eval(valid, sentinel_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> GpuTuple<Long>
            { return {array(i,j,k) == tau_sentinel ? Long(1) : Long(0)}; });

        const Box& source = source_ba[mfi.index()];
        const bool selected = face.isLow()
            ? source.smallEnd(dir) == source_face_index
            : source.bigEnd(dir) == source_face_index;
        if (!selected) { continue; }

        // Stress fields retain their three-dimensional valid boxes.  Only
        // the selected wall-normal plane is populated by the surface-layer
        // kernel; the remaining points in a face-owned FAB are intentionally
        // left at the sentinel value.
        const Box selected_plane = valid & face_plane;
        if (selected_plane.isEmpty()) { continue; }

        local_count += static_cast<Long>(selected_plane.numPts());
        range_op.eval(selected_plane, range_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
                -> GpuTuple<Real, Real>
            {
                const Real value = array(i,j,k);
                return {value, value};
            });
    }

    Gpu::streamSynchronize();
    const auto range = range_data.value();
    local_lo = get<0>(range);
    local_hi = get<1>(range);
    local_sentinel_count = get<0>(sentinel_data.value());
    ParallelDescriptor::ReduceLongSum(local_count);
    ParallelDescriptor::ReduceLongSum(local_total_count);
    ParallelDescriptor::ReduceLongSum(local_sentinel_count);
    ParallelDescriptor::ReduceRealMin(local_lo);
    ParallelDescriptor::ReduceRealMax(local_hi);
    return {local_lo, local_hi, local_count, local_total_count,
            local_sentinel_count};
}

struct GlobalValueCounts { Long total{0}; Long finite{0}; Long sentinel{0}; };

inline GlobalValueCounts
value_counts (const MultiFab& mf)
{
    ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_op;
    ReduceData<Long, Long, Long> reduce_data(reduce_op);
    for (MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const auto array = mf.const_array(mfi);
        reduce_op.eval(mfi.validbox(), reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
                -> GpuTuple<Long, Long, Long>
            {
                const Real value = array(i,j,k);
                return {Long(1),
                        (std::isfinite(value) && value != tau_sentinel)
                            ? Long(1) : Long(0),
                        value == tau_sentinel ? Long(1) : Long(0)};
            });
    }
    Gpu::streamSynchronize();
    auto counts = reduce_data.value();
    Long total = get<0>(counts), finite = get<1>(counts), sentinel = get<2>(counts);
    ParallelDescriptor::ReduceLongSum(total);
    ParallelDescriptor::ReduceLongSum(finite);
    ParallelDescriptor::ReduceLongSum(sentinel);
    return {total, finite, sentinel};
}

inline bool is_changed (const Real value)
{
    return std::isfinite(value) && value != tau_sentinel;
}

inline GpuArray<int, AMREX_SPACEDIM*2>
active_face (const Orientation face)
{
    GpuArray<int, AMREX_SPACEDIM*2> result{};
    result[static_cast<int>(face)] = 1;
    return result;
}

} // namespace erf_surface_layer_test

#endif
