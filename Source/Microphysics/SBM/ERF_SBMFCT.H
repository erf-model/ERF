#ifndef ERF_SBM_FCT_H_
#define ERF_SBM_FCT_H_

#include "ERF_SBMConstraintGroups.H"
#include "ERF_AuxiliaryStageContext.H"

#include <AMReX_REAL.H>

#include <string>
#include <vector>

namespace erf_sbm {

enum class TransportMethod { DonorCell, GroupedFCT_WENOZ3 };

/** A face is oriented from left to right and owns exactly one transfer. */
struct FCTFaceTransfer {
    int left_cell{0};
    int right_cell{0};
    amrex::Real left_volume{1.0};
    amrex::Real right_volume{1.0};
    std::vector<amrex::Real> low;
    std::vector<amrex::Real> high;
};

struct FCTResult {
    std::vector<FCTFaceTransfer> accepted_faces;
    std::vector<amrex::Real> limiter;
    std::vector<amrex::Real> updated_state;
    amrex::Real minimum_margin{0.0};
};

/**
 * Apply one invariant-domain limiter to each shared face.  `low_state` is the
 * already-updated donor/diffusive state and face vectors are instantaneous
 * fluxes multiplied by the stage interval (integrated transfers).  A single
 * lambda is selected per face over every affected group constraint.
 */
[[nodiscard]] FCTResult
limit_grouped(const std::vector<amrex::Real>& low_state,
              int ncell, int ncomp,
              const std::vector<FCTFaceTransfer>& faces,
              const std::vector<ConstraintGroup>& groups,
              int chunk_size = 0);

/** Difference of two accepted stage transfers, preserving component order. */
[[nodiscard]] std::vector<amrex::Real>
accepted_stage_correction(const std::vector<amrex::Real>& low_transfer,
                          const std::vector<amrex::Real>& high_transfer,
                          amrex::Real stage_interval);

struct StageWeightContract {
    amrex::Real baseline_weight{1.0};
    amrex::Real low_order_weight{1.0};
    amrex::Real correction_weight{1.0};
    amrex::Real completed_ledger_weight{0.0};
};

/** Exact weights used by the two supported ERF scalar recurrences. */
[[nodiscard]] StageWeightContract
stage_weight_contract(::erf_auxiliary::IntegrationMethod method, int stage,
                      bool completes_step) noexcept;

} // namespace erf_sbm

#endif
