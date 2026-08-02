#ifndef ERF_CLOUD_CHAMBER_WATER_LEDGER_H_
#define ERF_CLOUD_CHAMBER_WATER_LEDGER_H_

#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParReduce.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Array4.H>

#include <array>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>

#include <ERF_IndexDefines.H>

/**
 * Opt-in diagnostic ledger for the two non-precipitating water components.
 *
 * The ledger is deliberately separate from CloudChamberBudget.  The budget
 * answers whether the final state closes against resolved wall fluxes; this
 * ledger records the intermediate state transitions and the amount removed
 * by the qv/qc positivity safeguard.  It is disabled unless explicitly
 * requested with erf.cloud_chamber_water_ledger = true.
 */
class CloudChamberWaterLedger
{
public:
    enum Water : int { Qv = 0, Qc = 1, NumWater = 2 };

    enum DiagnosticComponent : int {
        DeltaQv = 0,
        DeltaQc,
        AbsDeltaQv,
        AbsDeltaQc,
        RawXQv,
        RawXQc,
        RawYQv,
        RawYQc,
        RawZQv,
        RawZQc,
        WeightedXQv,
        WeightedXQc,
        WeightedYQv,
        WeightedYQc,
        WeightedZQv,
        WeightedZQc,
        BoundaryXQv,
        BoundaryXQc,
        BoundaryYQv,
        BoundaryYQc,
        BoundaryZQv,
        BoundaryZQc,
        AbsRawQv,
        AbsRawQc,
        NumDiagnosticComponents
    };

    enum FaceAuditComponent : int {
        FaceSumQv = 0,
        FaceSumQc,
        FaceMinQv,
        FaceMinQc,
        FaceMaxQv,
        FaceMaxQc,
        FaceCount,
        NumFaceAuditComponents
    };

    explicit CloudChamberWaterLedger (bool enabled) : m_enabled(enabled) {}

    bool enabled () const noexcept { return m_enabled; }

    void begin_step (const amrex::MultiFab& cons,
                     const amrex::Geometry& geom,
                     double time,
                     double dt)
    {
        if (!m_enabled) { return; }
        reset_device_values();
        if (!m_advection_diagnostics) {
            m_advection_diagnostics = std::make_unique<amrex::MultiFab>(
                cons.boxArray(), cons.DistributionMap(), NumDiagnosticComponents, 0);
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                m_face_audit[dir] = std::make_unique<amrex::MultiFab>(
                    amrex::convert(cons.boxArray(), amrex::IntVect::TheDimensionVector(dir)),
                    cons.DistributionMap(), NumFaceAuditComponents, 0);
            }
        }
        m_advection_diagnostics->setVal(amrex::Real(0.0));
        reset_face_audit();
        m_time = time;
        m_dt = dt;
        m_start = state_integrals(cons, geom);
        m_before_slow = m_start;
        m_after_dycore = {};
        m_after_microphysics = {};
        m_end = {};
    }

    void begin_advection_stage ()
    {
        if (m_enabled) { reset_face_audit(); }
    }

    void record_after_dycore (const amrex::MultiFab& cons,
                              const amrex::Geometry& geom) {
        if (m_enabled) { m_after_dycore = state_integrals(cons, geom); }
    }

    void record_after_microphysics (const amrex::MultiFab& cons,
                                    const amrex::Geometry& geom) {
        if (m_enabled) { m_after_microphysics = state_integrals(cons, geom); }
    }

    void finish_step (const amrex::MultiFab& cons,
                      const amrex::Geometry& geom,
                      double time)
    {
        if (!m_enabled) { return; }
        m_end = state_integrals(cons, geom);

        const amrex::Real volume = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
        constexpr int base_columns = 7 * NumWater + 6 * NumWater;
        constexpr int independent_columns = 13 * NumWater;
        constexpr int face_columns = 5 * AMREX_SPACEDIM;
        std::array<amrex::Real, base_columns + independent_columns + face_columns> packed{};
        int n = 0;
        for (const auto value : m_start) { packed[n++] = value; }
        for (const auto value : m_before_slow) { packed[n++] = value; }
        packed[n++] = m_candidate_qv.dataValue() * volume;
        packed[n++] = m_candidate_qc.dataValue() * volume;
        packed[n++] = m_after_positivity_qv.dataValue() * volume;
        packed[n++] = m_after_positivity_qc.dataValue() * volume;
        for (const auto value : m_after_dycore) { packed[n++] = value; }
        for (const auto value : m_after_microphysics) { packed[n++] = value; }
        for (const auto value : m_end) { packed[n++] = value; }
        packed[n++] = m_pre_advection_rhs_qv.dataValue() * volume * m_dt;
        packed[n++] = m_pre_advection_rhs_qc.dataValue() * volume * m_dt;
        packed[n++] = m_advective_boundary_qv.dataValue() * m_dt;
        packed[n++] = m_advective_boundary_qc.dataValue() * m_dt;
        packed[n++] = m_advection_rhs_qv.dataValue() * volume * m_dt;
        packed[n++] = m_advection_rhs_qc.dataValue() * volume * m_dt;
        packed[n++] = m_raw_advection_divergence_qv.dataValue() * volume * m_dt;
        packed[n++] = m_raw_advection_divergence_qc.dataValue() * volume * m_dt;
        packed[n++] = m_weighted_advection_divergence_qv.dataValue() * volume * m_dt;
        packed[n++] = m_weighted_advection_divergence_qc.dataValue() * volume * m_dt;
        packed[n++] = m_clip_qv.dataValue() * volume;
        packed[n++] = m_clip_qc.dataValue() * volume;

        packed[n++] = m_advection_delta_qv.dataValue() * volume * m_dt;
        packed[n++] = m_advection_delta_qc.dataValue() * volume * m_dt;

        using SumOps = amrex::TypeList<
            amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum,
            amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum,
            amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum,
            amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum,
            amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum,
            amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum>;
        using SumTypes = amrex::TypeList<
            amrex::Real, amrex::Real, amrex::Real, amrex::Real,
            amrex::Real, amrex::Real, amrex::Real, amrex::Real,
            amrex::Real, amrex::Real, amrex::Real, amrex::Real,
            amrex::Real, amrex::Real, amrex::Real, amrex::Real,
            amrex::Real, amrex::Real, amrex::Real, amrex::Real,
            amrex::Real, amrex::Real, amrex::Real, amrex::Real>;
        const auto& diagnostic_arrays = m_advection_diagnostics->const_arrays();
        const auto stable = amrex::ParReduce(
            SumOps{}, SumTypes{}, *m_advection_diagnostics, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept {
                const auto& a = diagnostic_arrays[box_no];
                return amrex::GpuTuple<amrex::Real, amrex::Real, amrex::Real, amrex::Real,
                                       amrex::Real, amrex::Real, amrex::Real, amrex::Real,
                                       amrex::Real, amrex::Real, amrex::Real, amrex::Real,
                                       amrex::Real, amrex::Real, amrex::Real, amrex::Real,
                                       amrex::Real, amrex::Real, amrex::Real, amrex::Real,
                                       amrex::Real, amrex::Real, amrex::Real, amrex::Real>{
                    a(i,j,k,DeltaQv), a(i,j,k,DeltaQc),
                    a(i,j,k,AbsDeltaQv), a(i,j,k,AbsDeltaQc),
                    a(i,j,k,RawXQv), a(i,j,k,RawXQc),
                    a(i,j,k,RawYQv), a(i,j,k,RawYQc),
                    a(i,j,k,RawZQv), a(i,j,k,RawZQc),
                    a(i,j,k,WeightedXQv), a(i,j,k,WeightedXQc),
                    a(i,j,k,WeightedYQv), a(i,j,k,WeightedYQc),
                    a(i,j,k,WeightedZQv), a(i,j,k,WeightedZQc),
                    a(i,j,k,BoundaryXQv), a(i,j,k,BoundaryXQc),
                    a(i,j,k,BoundaryYQv), a(i,j,k,BoundaryYQc),
                    a(i,j,k,BoundaryZQv), a(i,j,k,BoundaryZQc),
                    a(i,j,k,AbsRawQv), a(i,j,k,AbsRawQc)};
            });
        const amrex::Real rhs_scale = volume * m_dt;
        const amrex::Real boundary_scale = m_dt;
        packed[n++] = amrex::get<0>(stable) * rhs_scale;
        packed[n++] = amrex::get<1>(stable) * rhs_scale;
        packed[n++] = amrex::get<2>(stable) * rhs_scale;
        packed[n++] = amrex::get<3>(stable) * rhs_scale;
        packed[n++] = amrex::get<4>(stable) * rhs_scale;
        packed[n++] = amrex::get<5>(stable) * rhs_scale;
        packed[n++] = amrex::get<6>(stable) * rhs_scale;
        packed[n++] = amrex::get<7>(stable) * rhs_scale;
        packed[n++] = amrex::get<8>(stable) * rhs_scale;
        packed[n++] = amrex::get<9>(stable) * rhs_scale;
        packed[n++] = amrex::get<10>(stable) * rhs_scale;
        packed[n++] = amrex::get<11>(stable) * rhs_scale;
        packed[n++] = amrex::get<12>(stable) * rhs_scale;
        packed[n++] = amrex::get<13>(stable) * rhs_scale;
        packed[n++] = amrex::get<14>(stable) * rhs_scale;
        packed[n++] = amrex::get<15>(stable) * rhs_scale;
        packed[n++] = amrex::get<16>(stable) * boundary_scale;
        packed[n++] = amrex::get<17>(stable) * boundary_scale;
        packed[n++] = amrex::get<18>(stable) * boundary_scale;
        packed[n++] = amrex::get<19>(stable) * boundary_scale;
        packed[n++] = amrex::get<20>(stable) * boundary_scale;
        packed[n++] = amrex::get<21>(stable) * boundary_scale;
        packed[n++] = amrex::get<22>(stable) * rhs_scale;
        packed[n++] = amrex::get<23>(stable) * rhs_scale;

        const int sum_columns = n;
        const auto face_summary = face_audit_summary();
        for (const auto value : face_summary) { packed[n++] = value; }
        amrex::ParallelDescriptor::ReduceRealSum(packed.data(), sum_columns);

        const amrex::Real min_qv = global_min(m_min_candidate_qv.dataValue());
        const amrex::Real min_qc = global_min(m_min_candidate_qc.dataValue());
        const amrex::Real largest_qv = global_max(m_largest_deficit_qv.dataValue());
        const amrex::Real largest_qc = global_max(m_largest_deficit_qc.dataValue());
        const amrex::Long negative_qv = global_count(m_negative_qv.dataValue());
        const amrex::Long negative_qc = global_count(m_negative_qc.dataValue());

        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::ofstream out("cloud_chamber_water_ledger.dat", std::ios::app);
            if (out.tellp() == std::streampos(0)) {
                out << "# time dt "
                       "W_timestep_start_qv W_timestep_start_qc "
                       "W_before_slow_moisture_update_qv W_before_slow_moisture_update_qc "
                       "W_candidate_before_positivity_qv W_candidate_before_positivity_qc "
                       "W_after_positivity_qv W_after_positivity_qc "
                       "W_before_microphysics_qv W_before_microphysics_qc "
                       "W_after_microphysics_qv W_after_microphysics_qc "
                       "W_timestep_end_qv W_timestep_end_qc "
                       "pre_advection_rhs_integral_qv pre_advection_rhs_integral_qc "
                       "advective_boundary_integral_qv advective_boundary_integral_qc "
                       "advection_rhs_integral_qv advection_rhs_integral_qc "
                       "raw_advection_divergence_integral_qv raw_advection_divergence_integral_qc "
                       "weighted_advection_divergence_integral_qv weighted_advection_divergence_integral_qc "
                       "qv_clipping_increment qc_clipping_increment "
                       "advection_rhs_delta_atomic_qv advection_rhs_delta_atomic_qc "
                       "advection_rhs_delta_stable_qv advection_rhs_delta_stable_qc "
                       "advection_rhs_delta_abs_qv advection_rhs_delta_abs_qc "
                       "stable_raw_divergence_x_qv stable_raw_divergence_x_qc "
                       "stable_raw_divergence_y_qv stable_raw_divergence_y_qc "
                       "stable_raw_divergence_z_qv stable_raw_divergence_z_qc "
                       "stable_weighted_divergence_x_qv stable_weighted_divergence_x_qc "
                       "stable_weighted_divergence_y_qv stable_weighted_divergence_y_qc "
                       "stable_weighted_divergence_z_qv stable_weighted_divergence_z_qc "
                       "stable_boundary_flux_x_qv stable_boundary_flux_x_qc "
                       "stable_boundary_flux_y_qv stable_boundary_flux_y_qc "
                       "stable_boundary_flux_z_qv stable_boundary_flux_z_qc "
                       "stable_raw_divergence_abs_qv stable_raw_divergence_abs_qc "
                       "tile_face_mismatch_max_x_qv tile_face_mismatch_max_x_qc "
                       "tile_face_mismatch_abs_sum_x_qv tile_face_mismatch_abs_sum_x_qc "
                       "tile_face_mismatch_count_x "
                       "tile_face_mismatch_max_y_qv tile_face_mismatch_max_y_qc "
                       "tile_face_mismatch_abs_sum_y_qv tile_face_mismatch_abs_sum_y_qc "
                       "tile_face_mismatch_count_y "
                       "tile_face_mismatch_max_z_qv tile_face_mismatch_max_z_qc "
                       "tile_face_mismatch_abs_sum_z_qv tile_face_mismatch_abs_sum_z_qc "
                       "tile_face_mismatch_count_z "
                       "qv_negative_cell_count qc_negative_cell_count "
                       "qv_largest_deficit qc_largest_deficit "
                       "qv_min_pre_positivity qc_min_pre_positivity\n";
            }
            out << std::setprecision(17) << time << ' ' << m_dt;
            for (const auto value : packed) { out << ' ' << value; }
            out << ' ' << negative_qv << ' ' << negative_qc
                << ' ' << largest_qv << ' ' << largest_qc
                << ' ' << min_qv << ' ' << min_qc << '\n';
        }

        (void) time;
    }

    amrex::Real* candidate_ptr (int water) noexcept {
        return water == Qv ? m_candidate_qv.dataPtr() : m_candidate_qc.dataPtr();
    }

    amrex::Real* after_positivity_ptr (int water) noexcept {
        return water == Qv ? m_after_positivity_qv.dataPtr() : m_after_positivity_qc.dataPtr();
    }

    amrex::Real* clipping_ptr (int water) noexcept {
        return water == Qv ? m_clip_qv.dataPtr() : m_clip_qc.dataPtr();
    }

    amrex::Real* advective_boundary_ptr (int water) noexcept {
        return water == Qv ? m_advective_boundary_qv.dataPtr() : m_advective_boundary_qc.dataPtr();
    }

    amrex::Real* advection_rhs_ptr (int water) noexcept {
        return water == Qv ? m_advection_rhs_qv.dataPtr() : m_advection_rhs_qc.dataPtr();
    }

    amrex::Real* pre_advection_rhs_ptr (int water) noexcept {
        return water == Qv ? m_pre_advection_rhs_qv.dataPtr() : m_pre_advection_rhs_qc.dataPtr();
    }

    amrex::Real* raw_advection_divergence_ptr (int water) noexcept {
        return water == Qv ? m_raw_advection_divergence_qv.dataPtr() : m_raw_advection_divergence_qc.dataPtr();
    }

    amrex::Real* weighted_advection_divergence_ptr (int water) noexcept {
        return water == Qv ? m_weighted_advection_divergence_qv.dataPtr() : m_weighted_advection_divergence_qc.dataPtr();
    }

    amrex::Real* advection_delta_ptr (int water) noexcept {
        return water == Qv ? m_advection_delta_qv.dataPtr() : m_advection_delta_qc.dataPtr();
    }

    amrex::Real* raw_advection_divergence_abs_ptr (int water) noexcept {
        return water == Qv ? m_abs_raw_advection_qv.dataPtr() : m_abs_raw_advection_qc.dataPtr();
    }

    amrex::Array4<amrex::Real> advection_diagnostics_array (const amrex::MFIter& mfi) noexcept {
        return m_advection_diagnostics->array(mfi);
    }

    amrex::Array4<amrex::Real> face_audit_array (int dir, const amrex::MFIter& mfi) noexcept {
        return m_face_audit[dir]->array(mfi);
    }

    amrex::Long* negative_count_ptr (int water) noexcept {
        return water == Qv ? m_negative_qv.dataPtr() : m_negative_qc.dataPtr();
    }

    amrex::Real* minimum_candidate_ptr (int water) noexcept {
        return water == Qv ? m_min_candidate_qv.dataPtr() : m_min_candidate_qc.dataPtr();
    }

    amrex::Real* largest_deficit_ptr (int water) noexcept {
        return water == Qv ? m_largest_deficit_qv.dataPtr() : m_largest_deficit_qc.dataPtr();
    }

private:
    using State = std::array<amrex::Real, NumWater>;

    static State state_integrals (const amrex::MultiFab& cons,
                                  const amrex::Geometry& geom)
    {
        const amrex::Real volume = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
        return {
            (RhoQ1_comp < cons.nComp()) ? cons.sum(RhoQ1_comp, true) * volume : amrex::Real(0.0),
            (RhoQ2_comp < cons.nComp()) ? cons.sum(RhoQ2_comp, true) * volume : amrex::Real(0.0)};
    }

    template <typename T>
    static void set_device_value (amrex::Gpu::DeviceScalar<T>& value, T host_value)
    {
        amrex::Gpu::htod_memcpy(value.dataPtr(), &host_value, sizeof(T));
    }

    void reset_device_values ()
    {
        set_device_value(m_candidate_qv, amrex::Real(0.0));
        set_device_value(m_candidate_qc, amrex::Real(0.0));
        set_device_value(m_after_positivity_qv, amrex::Real(0.0));
        set_device_value(m_after_positivity_qc, amrex::Real(0.0));
        set_device_value(m_clip_qv, amrex::Real(0.0));
        set_device_value(m_clip_qc, amrex::Real(0.0));
        set_device_value(m_advective_boundary_qv, amrex::Real(0.0));
        set_device_value(m_advective_boundary_qc, amrex::Real(0.0));
        set_device_value(m_advection_rhs_qv, amrex::Real(0.0));
        set_device_value(m_advection_rhs_qc, amrex::Real(0.0));
        set_device_value(m_pre_advection_rhs_qv, amrex::Real(0.0));
        set_device_value(m_pre_advection_rhs_qc, amrex::Real(0.0));
        set_device_value(m_raw_advection_divergence_qv, amrex::Real(0.0));
        set_device_value(m_raw_advection_divergence_qc, amrex::Real(0.0));
        set_device_value(m_weighted_advection_divergence_qv, amrex::Real(0.0));
        set_device_value(m_weighted_advection_divergence_qc, amrex::Real(0.0));
        set_device_value(m_advection_delta_qv, amrex::Real(0.0));
        set_device_value(m_advection_delta_qc, amrex::Real(0.0));
        set_device_value(m_abs_raw_advection_qv, amrex::Real(0.0));
        set_device_value(m_abs_raw_advection_qc, amrex::Real(0.0));
        set_device_value(m_negative_qv, amrex::Long(0));
        set_device_value(m_negative_qc, amrex::Long(0));
        set_device_value(m_min_candidate_qv, std::numeric_limits<amrex::Real>::max());
        set_device_value(m_min_candidate_qc, std::numeric_limits<amrex::Real>::max());
        set_device_value(m_largest_deficit_qv, amrex::Real(0.0));
        set_device_value(m_largest_deficit_qc, amrex::Real(0.0));
    }

    void reset_face_audit ()
    {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            auto& face = *m_face_audit[dir];
            face.setVal(amrex::Real(0.0), FaceSumQv, 2, 0);
            face.setVal(std::numeric_limits<amrex::Real>::max(), FaceMinQv, 2, 0);
            face.setVal(-std::numeric_limits<amrex::Real>::max(), FaceMaxQv, 2, 0);
            face.setVal(amrex::Real(0.0), FaceCount, 1, 0);
        }
    }

    std::array<amrex::Real, 5 * AMREX_SPACEDIM> face_audit_summary () const
    {
        std::array<amrex::Real, 5 * AMREX_SPACEDIM> result{};
        using Ops = amrex::TypeList<
            amrex::ReduceOpMax, amrex::ReduceOpMax,
            amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum>;
        using Types = amrex::TypeList<
            amrex::Real, amrex::Real, amrex::Real, amrex::Real, amrex::Real>;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            const auto& arrays = m_face_audit[dir]->const_arrays();
            const auto local = amrex::ParReduce(
                Ops{}, Types{}, *m_face_audit[dir], amrex::IntVect(0),
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept {
                    const auto& a = arrays[box_no];
                    const amrex::Real count = a(i,j,k,FaceCount);
                    const bool shared = count > amrex::Real(1.5);
                    const amrex::Real qv_mismatch = shared ?
                        a(i,j,k,FaceMaxQv) - a(i,j,k,FaceMinQv) : amrex::Real(0.0);
                    const amrex::Real qc_mismatch = shared ?
                        a(i,j,k,FaceMaxQc) - a(i,j,k,FaceMinQc) : amrex::Real(0.0);
                    return amrex::GpuTuple<amrex::Real, amrex::Real,
                                           amrex::Real, amrex::Real, amrex::Real>{
                        qv_mismatch, qc_mismatch,
                        qv_mismatch, qc_mismatch,
                        shared ? amrex::Real(1.0) : amrex::Real(0.0)};
                });
            const int offset = 5 * dir;
            result[offset + 0] = amrex::get<0>(local);
            result[offset + 1] = amrex::get<1>(local);
            result[offset + 2] = amrex::get<2>(local);
            result[offset + 3] = amrex::get<3>(local);
            result[offset + 4] = amrex::get<4>(local);
            amrex::ParallelDescriptor::ReduceRealMax(result.data() + offset, 2);
            amrex::ParallelDescriptor::ReduceRealSum(result.data() + offset + 2, 3);
        }
        return result;
    }

    static amrex::Real global_min (amrex::Real value)
    {
        amrex::ParallelDescriptor::ReduceRealMin(value);
        return value;
    }

    static amrex::Real global_max (amrex::Real value)
    {
        amrex::ParallelDescriptor::ReduceRealMax(value);
        return value;
    }

    static amrex::Long global_count (amrex::Long value)
    {
        amrex::ParallelDescriptor::ReduceLongSum(value);
        return value;
    }

    bool m_enabled = false;
    double m_time = 0.0;
    double m_dt = 0.0;
    State m_start{};
    State m_before_slow{};
    State m_after_dycore{};
    State m_after_microphysics{};
    State m_end{};

    amrex::Gpu::DeviceScalar<amrex::Real> m_candidate_qv{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_candidate_qc{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_after_positivity_qv{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_after_positivity_qc{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_clip_qv{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_clip_qc{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_advective_boundary_qv{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_advective_boundary_qc{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_advection_rhs_qv{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_advection_rhs_qc{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_pre_advection_rhs_qv{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_pre_advection_rhs_qc{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_raw_advection_divergence_qv{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_raw_advection_divergence_qc{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_weighted_advection_divergence_qv{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_weighted_advection_divergence_qc{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_advection_delta_qv{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_advection_delta_qc{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_abs_raw_advection_qv{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_abs_raw_advection_qc{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Long> m_negative_qv{amrex::Long(0)};
    amrex::Gpu::DeviceScalar<amrex::Long> m_negative_qc{amrex::Long(0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_min_candidate_qv{std::numeric_limits<amrex::Real>::max()};
    amrex::Gpu::DeviceScalar<amrex::Real> m_min_candidate_qc{std::numeric_limits<amrex::Real>::max()};
    amrex::Gpu::DeviceScalar<amrex::Real> m_largest_deficit_qv{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_largest_deficit_qc{amrex::Real(0.0)};
    std::unique_ptr<amrex::MultiFab> m_advection_diagnostics;
    std::array<std::unique_ptr<amrex::MultiFab>, AMREX_SPACEDIM> m_face_audit;
};

#endif
