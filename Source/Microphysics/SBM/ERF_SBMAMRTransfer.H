#ifndef ERF_SBM_AMR_TRANSFER_H_
#define ERF_SBM_AMR_TRANSFER_H_

#include "ERF_SBMLayout.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <string>
#include <cstdint>
#include <optional>
#include <vector>

namespace erf_sbm {

class SBMStateAdmission;

/** Optional per-call measurement of current AMReX FAB payload above the
 * caller-established baseline. Samples are taken at actual operator
 * allocation high-water points; this is local-rank FAB payload, not RSS.
 */
struct SBMTransferAllocationCensus {
    std::uint64_t baseline_bytes{0};
    std::uint64_t peak_bytes{0};
    std::uint64_t support_loss_collective_votes{0};
    // Optional test trace: first global component of each processed complete
    // semantic chunk. Empty in production unless explicitly supplied by a test.
    std::vector<int> processed_chunk_first_components;
    void begin () noexcept;
    void sample () noexcept;
    [[nodiscard]] std::uint64_t peak_extra_bytes () const noexcept
    { return peak_bytes > baseline_bytes ? peak_bytes - baseline_bytes : 0; }
};

/** A spectrum, carrier, and mapped measure tied to explicit semantic times. */
struct SBMAMRStateView {
    const amrex::MultiFab* spectrum{nullptr};
    double spectrum_time{0.0};
    const amrex::MultiFab* dry_air_density{nullptr};
    int density_component{-1};
    double density_time{0.0};
    const amrex::MultiFab* mapped_measure{nullptr};
    int measure_component{-1};
    double measure_time{0.0};
};

/**
 * Form an H-conservative coarse candidate from matching-time spectral views.
 * Cross-level endpoint labels may use a timestep-bounded comparison when the
 * caller supplies the interval's actual minimum physical timestep; each input
 * tuple's spectrum/density/measure labels remain exactly identical.
 * On success, overwrite the destination's valid cells with the complete,
 * admitted result. On failure (including a rejection on another MPI rank),
 * leave every destination valid cell and all source views unchanged. This
 * operator does not itself commit authoritative SBM state.
 */
bool RestrictMappedSpectrum (const SBMLayout& layout,
                             const SBMAMRStateView& fine,
                             const SBMAMRStateView& coarse,
                             const amrex::IntVect& ratio,
                             int coarse_level,
                             amrex::MultiFab& coarse_candidate,
                             std::string& diagnostic,
                             int max_groups_per_chunk = 16,
                             bool reverse_chunk_order = false,
                             const SBMStateAdmission* admission = nullptr,
                             SBMTransferAllocationCensus* allocation_census = nullptr,
                             std::optional<double> cross_level_dt_min = std::nullopt);

/**
 * Private-candidate variant for coarse synchronization. The caller supplies
 * an unpublished full-spectrum buffer; it may be overwritten on failure.
 * The routine reuses it for the restricted mapped amount and final spectrum.
 */
bool RestrictMappedSpectrumPrivate (const SBMLayout& layout,
                                    const SBMAMRStateView& fine,
                                    const SBMAMRStateView& coarse,
                                    const amrex::IntVect& ratio,
                                    int coarse_level,
                                    amrex::MultiFab& private_candidate,
                                    std::string& diagnostic,
                                    int max_groups_per_chunk = 16,
                                    bool reverse_chunk_order = false,
                                    const SBMStateAdmission* admission = nullptr,
                                    SBMTransferAllocationCensus* allocation_census = nullptr,
                                    std::optional<double> cross_level_dt_min = std::nullopt,
                                    const amrex::MultiFab* covered_parent_mask = nullptr);

/**
 * Form a same-time piecewise-constant carrier-relative fine candidate. On
 * success, overwrite the destination's valid cells with the complete, admitted
 * result. On failure (including a rejection on another MPI rank), leave every
 * destination valid cell and all source views unchanged. This operator does
 * not itself commit authoritative SBM state.
 */
bool ProlongCarrierRelativeSpectrum (const SBMLayout& layout,
                                     const SBMAMRStateView& coarse,
                                     const SBMAMRStateView& fine_target,
                                     const amrex::Geometry& coarse_geometry,
                                     const amrex::Geometry& fine_geometry,
                                     const amrex::IntVect& ratio,
                                     int fine_level,
                                     amrex::MultiFab& fine_candidate,
                                     std::string& diagnostic,
                                     int max_groups_per_chunk = 16,
                                     bool reverse_chunk_order = false,
                                     const SBMStateAdmission* admission = nullptr,
                                     SBMTransferAllocationCensus* allocation_census = nullptr);

/**
 * Private-candidate prolongation for a new level. The caller supplies the one
 * unpublished full fine spectrum; it may be overwritten on failure. Intensive
 * interpolation scratch remains bounded by complete constraint-group chunks.
 */
bool ProlongCarrierRelativeSpectrumPrivate (
    const SBMLayout& layout,
    const SBMAMRStateView& coarse,
    const SBMAMRStateView& fine_target,
    const amrex::Geometry& coarse_geometry,
    const amrex::Geometry& fine_geometry,
    const amrex::IntVect& ratio,
    int fine_level,
    amrex::MultiFab& private_candidate,
    std::string& diagnostic,
    int max_groups_per_chunk = 16,
    bool reverse_chunk_order = false,
    const SBMStateAdmission* admission = nullptr,
    SBMTransferAllocationCensus* allocation_census = nullptr);

/** Fine cells of a remade level, classified from grid metadata alone. */
struct SBMRegridCellCounts {
    amrex::Long cells{0};      // valid cells on the new grids
    amrex::Long retained{0};   // covered by both old and new grids
    amrex::Long prolonged{0};  // newly covered; filled from the parent
    amrex::Long released{0};   // covered only by the old grids
};

[[nodiscard]] SBMRegridCellCounts CountRegridCells (
    const amrex::BoxArray& old_grids, const amrex::BoxArray& new_grids);

/**
 * Private candidate for a fine level remade on new grids at a synchronized
 * instant. Cells the retained old-grid spectrum covers keep its values
 * bitwise; every other cell takes the same-time carrier-relative
 * prolongation of its parent. The new-grid density and measure, the retained
 * spectrum, and the coarse tuple must all carry semantic_time. Every parent
 * of the new grids must satisfy the mapped carrier identity and the complete
 * candidate must be admissible. Inputs are never modified; the candidate may
 * be overwritten on failure.
 */
bool RemapRegriddedSpectrumPrivate (
    const SBMLayout& layout,
    const SBMAMRStateView& coarse,
    const amrex::MultiFab& retained_spectrum,
    const amrex::MultiFab& fine_density, int density_component,
    const amrex::MultiFab& fine_measure, int measure_component,
    double semantic_time,
    const amrex::Geometry& coarse_geometry,
    const amrex::Geometry& fine_geometry,
    const amrex::IntVect& ratio,
    int fine_level,
    amrex::MultiFab& private_candidate,
    std::string& diagnostic,
    int max_groups_per_chunk = 16,
    const SBMStateAdmission* coarse_admission = nullptr,
    const SBMStateAdmission* fine_admission = nullptr,
    SBMTransferAllocationCensus* allocation_census = nullptr);

/** Add one bounded YA correction chunk in mapped H coordinates only where
 * the scalar covered-parent mask proves a coarse/fine face interface. The
 * candidate and invalid mask are private scratch; invalid is accumulated
 * across chunks and checked collectively after the final chunk.
 */
[[nodiscard]] bool ApplyMappedRefluxCorrectionChunk (
    amrex::MultiFab& private_candidate,
    const amrex::MultiFab& mapped_correction,
    const amrex::MultiFab& mapped_measure,
    const amrex::MultiFab& covered_parent_mask,
    const std::vector<int>& global_components,
    amrex::MultiFab& invalid_mask, std::string& diagnostic);

[[nodiscard]] bool MappedRefluxCorrectionMaskIsValid (
    const amrex::MultiFab& invalid_mask, std::string& diagnostic);

} // namespace erf_sbm

#endif
