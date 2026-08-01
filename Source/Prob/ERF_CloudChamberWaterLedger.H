#ifndef ERF_CLOUD_CHAMBER_WATER_LEDGER_H_
#define ERF_CLOUD_CHAMBER_WATER_LEDGER_H_

#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <array>
#include <fstream>
#include <iomanip>
#include <limits>

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

    explicit CloudChamberWaterLedger (bool enabled) : m_enabled(enabled) {}

    bool enabled () const noexcept { return m_enabled; }

    void begin_step (const amrex::MultiFab& cons,
                     const amrex::Geometry& geom,
                     double time,
                     double dt)
    {
        if (!m_enabled) { return; }
        reset_device_values();
        m_time = time;
        m_dt = dt;
        m_start = state_integrals(cons, geom);
        m_before_slow = m_start;
        m_after_dycore = {};
        m_after_microphysics = {};
        m_end = {};
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
        std::array<amrex::Real, 7 * NumWater + 6 * NumWater> packed{};
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
        amrex::ParallelDescriptor::ReduceRealSum(packed.data(), static_cast<int>(packed.size()));

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
        set_device_value(m_negative_qv, amrex::Long(0));
        set_device_value(m_negative_qc, amrex::Long(0));
        set_device_value(m_min_candidate_qv, std::numeric_limits<amrex::Real>::max());
        set_device_value(m_min_candidate_qc, std::numeric_limits<amrex::Real>::max());
        set_device_value(m_largest_deficit_qv, amrex::Real(0.0));
        set_device_value(m_largest_deficit_qc, amrex::Real(0.0));
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
    amrex::Gpu::DeviceScalar<amrex::Long> m_negative_qv{amrex::Long(0)};
    amrex::Gpu::DeviceScalar<amrex::Long> m_negative_qc{amrex::Long(0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_min_candidate_qv{std::numeric_limits<amrex::Real>::max()};
    amrex::Gpu::DeviceScalar<amrex::Real> m_min_candidate_qc{std::numeric_limits<amrex::Real>::max()};
    amrex::Gpu::DeviceScalar<amrex::Real> m_largest_deficit_qv{amrex::Real(0.0)};
    amrex::Gpu::DeviceScalar<amrex::Real> m_largest_deficit_qc{amrex::Real(0.0)};
};

#endif
