#ifndef ERF_CLOUD_CHAMBER_WATER_LEDGER_H_
#define ERF_CLOUD_CHAMBER_WATER_LEDGER_H_

#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Geometry.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_ParReduce.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Array4.H>
#include <AMReX_FabArrayUtility.H>

#include <algorithm>
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

    struct FaceMismatch
    {
        amrex::Real signed_mismatch = amrex::Real(0.0);
        amrex::Real absolute_mismatch = amrex::Real(0.0);
        amrex::Real tolerance = amrex::Real(0.0);
        bool mismatched = false;
    };

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

    AMREX_GPU_HOST_DEVICE static FaceMismatch compare_face_values (
        amrex::Real reference, amrex::Real duplicate,
        amrex::Real flux_floor = amrex::Real(1.0e-30)) noexcept
    {
        const amrex::Real scale = amrex::max(
            amrex::max(amrex::Math::abs(reference), amrex::Math::abs(duplicate)),
            flux_floor);
        const amrex::Real tolerance = amrex::Real(64.0) *
            std::numeric_limits<amrex::Real>::epsilon() * scale;
        const amrex::Real signed_mismatch = duplicate - reference;
        const amrex::Real absolute_mismatch = amrex::Math::abs(signed_mismatch);
        return {signed_mismatch, absolute_mismatch, tolerance,
                absolute_mismatch > tolerance};
    }

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
        m_face_summary.fill(amrex::Real(0.0));
        m_interface_summary.fill(amrex::Real(0.0));
        m_number_of_grids = 0;
        m_number_of_mfiter_tiles = 0;
        m_maximum_tiles_per_grid = 0;
        m_moisture_tiling_enabled = false;
        m_mpi_rank_count = 1;
        m_openmp_thread_count = 1;
        m_advection_stage_count = 0;
        m_last_advection_stage = -1;
    }

    void begin_advection_stage ()
    {
        if (m_enabled) {
            reset_face_audit();
            ++m_advection_stage_count;
        }
    }

    void record_advection_stage (int nrk)
    {
        if (!m_enabled) { return; }
        m_last_advection_stage = nrk;
        const auto stage_summary = face_audit_summary();
        const auto interface_summary = face_interface_summary();
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            const int offset = 7 * dir;
            m_face_summary[offset + 0] += stage_summary[offset + 0];
            m_face_summary[offset + 1] += stage_summary[offset + 1];
            m_face_summary[offset + 2] += stage_summary[offset + 2];
            m_face_summary[offset + 3] = amrex::max(
                m_face_summary[offset + 3], stage_summary[offset + 3]);
            m_face_summary[offset + 4] = amrex::max(
                m_face_summary[offset + 4], stage_summary[offset + 4]);
            m_face_summary[offset + 5] += stage_summary[offset + 5];
            m_face_summary[offset + 6] += stage_summary[offset + 6];
        }
        for (std::size_t i = 0; i < m_interface_summary.size(); ++i) {
            const bool is_maximum = (i % 6 == 1) || (i % 6 == 4);
            if (is_maximum) {
                m_interface_summary[i] = amrex::max(
                    m_interface_summary[i], interface_summary[i]);
            } else {
                m_interface_summary[i] += interface_summary[i];
            }
        }
    }

    void record_topology (int number_of_grids, int number_of_mfiter_tiles,
                          int maximum_tiles_per_grid,
                          bool moisture_tiling_enabled,
                          int mpi_rank_count, int openmp_thread_count)
    {
        if (!m_enabled) { return; }
        m_number_of_grids = number_of_grids;
        m_number_of_mfiter_tiles = std::max(m_number_of_mfiter_tiles,
                                            number_of_mfiter_tiles);
        m_maximum_tiles_per_grid = std::max(m_maximum_tiles_per_grid,
                                            maximum_tiles_per_grid);
        m_moisture_tiling_enabled = m_moisture_tiling_enabled ||
                                     moisture_tiling_enabled;
        m_mpi_rank_count = mpi_rank_count;
        m_openmp_thread_count = openmp_thread_count;
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
        constexpr int face_columns = 19 * AMREX_SPACEDIM;
        constexpr int topology_columns = 9;
        std::array<amrex::Real,
                   base_columns + independent_columns + face_columns + topology_columns> packed{};
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
        for (const auto value : m_face_summary) { packed[n++] = value; }
        // Keep the on-disk order aligned with the header: all grid-interface
        // qv/qc summaries by direction, followed by all MPI-interface
        // summaries by direction.  The in-memory array is grouped by
        // direction so each stage can update one face orientation together.
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            for (int component = 0; component < 6; ++component) {
                packed[n++] = m_interface_summary[12 * dir + component];
            }
        }
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            for (int component = 6; component < 12; ++component) {
                packed[n++] = m_interface_summary[12 * dir + component];
            }
        }
        packed[n++] = static_cast<amrex::Real>(m_number_of_grids);
        packed[n++] = static_cast<amrex::Real>(m_number_of_mfiter_tiles);
        packed[n++] = static_cast<amrex::Real>(m_maximum_tiles_per_grid);
        packed[n++] = m_moisture_tiling_enabled ? amrex::Real(1.0) : amrex::Real(0.0);
        packed[n++] = static_cast<amrex::Real>(m_mpi_rank_count);
        packed[n++] = static_cast<amrex::Real>(m_openmp_thread_count);
        packed[n++] = (m_mpi_rank_count == 1) ? amrex::Real(1.0) : amrex::Real(0.0);
        packed[n++] = static_cast<amrex::Real>(m_advection_stage_count);
        packed[n++] = static_cast<amrex::Real>(m_last_advection_stage);
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
                       "tile_face_shared_count_x "
                       "tile_face_mismatched_count_x_qv tile_face_mismatched_count_x_qc "
                       "tile_face_mismatch_max_x_qv tile_face_mismatch_max_x_qc "
                       "tile_face_mismatch_abs_sum_x_qv tile_face_mismatch_abs_sum_x_qc "
                       "tile_face_shared_count_y "
                       "tile_face_mismatched_count_y_qv tile_face_mismatched_count_y_qc "
                       "tile_face_mismatch_max_y_qv tile_face_mismatch_max_y_qc "
                       "tile_face_mismatch_abs_sum_y_qv tile_face_mismatch_abs_sum_y_qc "
                       "tile_face_shared_count_z "
                       "tile_face_mismatched_count_z_qv tile_face_mismatched_count_z_qc "
                       "tile_face_mismatch_max_z_qv tile_face_mismatch_max_z_qc "
                       "tile_face_mismatch_abs_sum_z_qv tile_face_mismatch_abs_sum_z_qc "
                       "grid_face_mismatched_count_x_qv grid_face_mismatch_max_x_qv grid_face_mismatch_abs_sum_x_qv "
                       "grid_face_mismatched_count_x_qc grid_face_mismatch_max_x_qc grid_face_mismatch_abs_sum_x_qc "
                       "grid_face_mismatched_count_y_qv grid_face_mismatch_max_y_qv grid_face_mismatch_abs_sum_y_qv "
                       "grid_face_mismatched_count_y_qc grid_face_mismatch_max_y_qc grid_face_mismatch_abs_sum_y_qc "
                       "grid_face_mismatched_count_z_qv grid_face_mismatch_max_z_qv grid_face_mismatch_abs_sum_z_qv "
                       "grid_face_mismatched_count_z_qc grid_face_mismatch_max_z_qc grid_face_mismatch_abs_sum_z_qc "
                       "mpi_face_mismatched_count_x_qv mpi_face_mismatch_max_x_qv mpi_face_mismatch_abs_sum_x_qv "
                       "mpi_face_mismatched_count_x_qc mpi_face_mismatch_max_x_qc mpi_face_mismatch_abs_sum_x_qc "
                       "mpi_face_mismatched_count_y_qv mpi_face_mismatch_max_y_qv mpi_face_mismatch_abs_sum_y_qv "
                       "mpi_face_mismatched_count_y_qc mpi_face_mismatch_max_y_qc mpi_face_mismatch_abs_sum_y_qc "
                       "mpi_face_mismatched_count_z_qv mpi_face_mismatch_max_z_qv mpi_face_mismatch_abs_sum_z_qv "
                       "mpi_face_mismatched_count_z_qc mpi_face_mismatch_max_z_qc mpi_face_mismatch_abs_sum_z_qc "
                       "number_of_grids number_of_mfiter_tiles maximum_tiles_per_grid "
                       "moisture_tiling_enabled mpi_rank_count openmp_thread_count "
                       "mpi_interface_audit_enabled "
                       "advection_stage_count last_advection_stage "
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

    std::array<amrex::Real, 7 * AMREX_SPACEDIM> face_audit_summary () const
    {
        std::array<amrex::Real, 7 * AMREX_SPACEDIM> result{};
        using Ops = amrex::TypeList<
            amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum,
            amrex::ReduceOpMax, amrex::ReduceOpMax,
            amrex::ReduceOpSum, amrex::ReduceOpSum>;
        using Types = amrex::TypeList<
            amrex::Real, amrex::Real, amrex::Real, amrex::Real,
            amrex::Real, amrex::Real, amrex::Real>;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            const auto& arrays = m_face_audit[dir]->const_arrays();
            const auto local = amrex::ParReduce(
                Ops{}, Types{}, *m_face_audit[dir], amrex::IntVect(0),
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept {
                    const auto& a = arrays[box_no];
                    const amrex::Real count = a(i,j,k,FaceCount);
                    const bool shared = count > amrex::Real(1.5);
                    const amrex::Real qv_spread = shared ?
                        a(i,j,k,FaceMaxQv) - a(i,j,k,FaceMinQv) : amrex::Real(0.0);
                    const amrex::Real qc_spread = shared ?
                        a(i,j,k,FaceMaxQc) - a(i,j,k,FaceMinQc) : amrex::Real(0.0);
                    const auto qv_check = compare_face_values(
                        a(i,j,k,FaceMinQv), a(i,j,k,FaceMaxQv));
                    const auto qc_check = compare_face_values(
                        a(i,j,k,FaceMinQc), a(i,j,k,FaceMaxQc));
                    return amrex::GpuTuple<amrex::Real, amrex::Real,
                                           amrex::Real, amrex::Real,
                                           amrex::Real, amrex::Real, amrex::Real>{
                        shared ? amrex::Real(1.0) : amrex::Real(0.0),
                        qv_check.mismatched ? amrex::Real(1.0) : amrex::Real(0.0),
                        qc_check.mismatched ? amrex::Real(1.0) : amrex::Real(0.0),
                        qv_spread, qc_spread, qv_spread, qc_spread};
                });
            const int offset = 7 * dir;
            result[offset + 0] = amrex::get<0>(local);
            result[offset + 1] = amrex::get<1>(local);
            result[offset + 2] = amrex::get<2>(local);
            result[offset + 3] = amrex::get<3>(local);
            result[offset + 4] = amrex::get<4>(local);
            result[offset + 5] = amrex::get<5>(local);
            result[offset + 6] = amrex::get<6>(local);
            amrex::ParallelDescriptor::ReduceRealSum(result.data() + offset, 3);
            amrex::ParallelDescriptor::ReduceRealMax(result.data() + offset + 3, 2);
            amrex::ParallelDescriptor::ReduceRealSum(result.data() + offset + 5, 2);
        }
        return result;
    }

    // The tile-local audit above sees duplicate writes inside one FAB.  This
    // opt-in owner comparison is separate: it preserves each grid's local
    // face value, synchronizes a copy from AMReX's lowest-numbered owner,
    // and classifies the remaining duplicates by owner rank.  It therefore
    // covers grid interfaces and, on supported distributed configurations,
    // MPI interfaces without changing production fluxes.
    std::array<amrex::Real, 12 * AMREX_SPACEDIM> face_interface_summary () const
    {
        std::array<amrex::Real, 12 * AMREX_SPACEDIM> result{};
        // AMReX 26.06's distributed nodal OverrideSync path is not safe for
        // this diagnostic MultiFab layout on the current Spack MPI build.
        // Do not let an optional ledger abort a valid multi-rank run; the
        // ledger explicitly marks this coverage unavailable in its schema.
        if (amrex::ParallelDescriptor::NProcs() > 1) {
            return result;
        }
        using Ops = amrex::TypeList<
            amrex::ReduceOpSum, amrex::ReduceOpMax, amrex::ReduceOpSum,
            amrex::ReduceOpSum, amrex::ReduceOpMax, amrex::ReduceOpSum,
            amrex::ReduceOpSum, amrex::ReduceOpMax, amrex::ReduceOpSum,
            amrex::ReduceOpSum, amrex::ReduceOpMax, amrex::ReduceOpSum>;
        using Types = amrex::TypeList<
            amrex::Real, amrex::Real, amrex::Real,
            amrex::Real, amrex::Real, amrex::Real,
            amrex::Real, amrex::Real, amrex::Real,
            amrex::Real, amrex::Real, amrex::Real>;

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            const auto& local_face = *m_face_audit[dir];
            amrex::MultiFab owner_values(local_face.boxArray(),
                                          local_face.DistributionMap(), 2, 0);
            amrex::MultiFab::Copy(owner_values, local_face, 0, FaceMinQv, 2, 0);

            auto owner_mask = owner_values.OwnerMask();
            owner_values.OverrideSync(*owner_mask);

            amrex::MultiFab local_rank(local_face.boxArray(),
                                       local_face.DistributionMap(), 1, 0);
            for (amrex::MFIter mfi(local_rank); mfi.isValid(); ++mfi) {
                local_rank[mfi].setVal(
                    static_cast<amrex::Real>(local_face.DistributionMap()[mfi.index()]));
            }
            local_rank.OverrideSync(*owner_mask,
                                    amrex::Periodicity::NonPeriodic());

            const auto& local_arrays = local_face.const_arrays();
            const auto& owner_arrays = owner_values.const_arrays();
            const auto& mask_arrays = owner_mask->const_arrays();
            const auto& rank_arrays = local_rank.const_arrays();
            const int local_rank_id = amrex::ParallelDescriptor::MyProc();
            const auto local = amrex::ParReduce(
                Ops{}, Types{}, owner_values, amrex::IntVect(0),
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept {
                    if (mask_arrays[box_no](i,j,k) != 0) {
                        return amrex::GpuTuple<
                            amrex::Real, amrex::Real, amrex::Real,
                            amrex::Real, amrex::Real, amrex::Real,
                            amrex::Real, amrex::Real, amrex::Real,
                            amrex::Real, amrex::Real, amrex::Real>{
                            amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0),
                            amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0),
                            amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0),
                            amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0)};
                    }
                    const auto qv_check = compare_face_values(
                        local_arrays[box_no](i,j,k,FaceMinQv),
                        owner_arrays[box_no](i,j,k,0));
                    const auto qc_check = compare_face_values(
                        local_arrays[box_no](i,j,k,FaceMinQc),
                        owner_arrays[box_no](i,j,k,1));
                    const bool mpi_duplicate =
                        rank_arrays[box_no](i,j,k) != local_rank_id;
                    const amrex::Real qv_count =
                        qv_check.mismatched ? amrex::Real(1.0) : amrex::Real(0.0);
                    const amrex::Real qc_count =
                        qc_check.mismatched ? amrex::Real(1.0) : amrex::Real(0.0);
                    const amrex::Real grid_qv_count = mpi_duplicate ? amrex::Real(0.0) : qv_count;
                    const amrex::Real grid_qc_count = mpi_duplicate ? amrex::Real(0.0) : qc_count;
                    const amrex::Real mpi_qv_count = mpi_duplicate ? qv_count : amrex::Real(0.0);
                    const amrex::Real mpi_qc_count = mpi_duplicate ? qc_count : amrex::Real(0.0);
                    return amrex::GpuTuple<
                        amrex::Real, amrex::Real, amrex::Real,
                        amrex::Real, amrex::Real, amrex::Real,
                        amrex::Real, amrex::Real, amrex::Real,
                        amrex::Real, amrex::Real, amrex::Real>{
                        grid_qv_count,
                        mpi_duplicate ? amrex::Real(0.0) : qv_check.absolute_mismatch,
                        mpi_duplicate ? amrex::Real(0.0) : qv_check.absolute_mismatch,
                        grid_qc_count,
                        mpi_duplicate ? amrex::Real(0.0) : qc_check.absolute_mismatch,
                        mpi_duplicate ? amrex::Real(0.0) : qc_check.absolute_mismatch,
                        mpi_qv_count,
                        mpi_duplicate ? qv_check.absolute_mismatch : amrex::Real(0.0),
                        mpi_duplicate ? qv_check.absolute_mismatch : amrex::Real(0.0),
                        mpi_qc_count,
                        mpi_duplicate ? qc_check.absolute_mismatch : amrex::Real(0.0),
                        mpi_duplicate ? qc_check.absolute_mismatch : amrex::Real(0.0)};
                });
            const int offset = 12 * dir;
            result[offset + 0] = amrex::get<0>(local);
            result[offset + 1] = amrex::get<1>(local);
            result[offset + 2] = amrex::get<2>(local);
            result[offset + 3] = amrex::get<3>(local);
            result[offset + 4] = amrex::get<4>(local);
            result[offset + 5] = amrex::get<5>(local);
            result[offset + 6] = amrex::get<6>(local);
            result[offset + 7] = amrex::get<7>(local);
            result[offset + 8] = amrex::get<8>(local);
            result[offset + 9] = amrex::get<9>(local);
            result[offset + 10] = amrex::get<10>(local);
            result[offset + 11] = amrex::get<11>(local);
            amrex::ParallelDescriptor::ReduceRealSum(result.data() + offset + 0, 1);
            amrex::ParallelDescriptor::ReduceRealMax(result.data() + offset + 1, 1);
            amrex::ParallelDescriptor::ReduceRealSum(result.data() + offset + 2, 1);
            amrex::ParallelDescriptor::ReduceRealSum(result.data() + offset + 3, 1);
            amrex::ParallelDescriptor::ReduceRealMax(result.data() + offset + 4, 1);
            amrex::ParallelDescriptor::ReduceRealSum(result.data() + offset + 5, 1);
            amrex::ParallelDescriptor::ReduceRealSum(result.data() + offset + 6, 1);
            amrex::ParallelDescriptor::ReduceRealMax(result.data() + offset + 7, 1);
            amrex::ParallelDescriptor::ReduceRealSum(result.data() + offset + 8, 1);
            amrex::ParallelDescriptor::ReduceRealSum(result.data() + offset + 9, 1);
            amrex::ParallelDescriptor::ReduceRealMax(result.data() + offset + 10, 1);
            amrex::ParallelDescriptor::ReduceRealSum(result.data() + offset + 11, 1);
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
    std::array<amrex::Real, 7 * AMREX_SPACEDIM> m_face_summary{};
    std::array<amrex::Real, 12 * AMREX_SPACEDIM> m_interface_summary{};
    int m_number_of_grids = 0;
    int m_number_of_mfiter_tiles = 0;
    int m_maximum_tiles_per_grid = 0;
    bool m_moisture_tiling_enabled = false;
    int m_mpi_rank_count = 1;
    int m_openmp_thread_count = 1;
    int m_advection_stage_count = 0;
    int m_last_advection_stage = -1;
};

#endif
