#ifndef ERF_SBM_DIFFUSION_H_
#define ERF_SBM_DIFFUSION_H_

#include "ERF_SBMConstraintGroups.H"

#include <AMReX_REAL.H>

#include <vector>

namespace erf_sbm {

/** One physical face for the explicit, density-weighted Fick operator. */
struct DiffusionFace {
    int left_cell{0};
    int right_cell{0};
    amrex::Real area{1.0};
    amrex::Real distance{1.0};
    amrex::Real left_volume{1.0};
    amrex::Real right_volume{1.0};
    amrex::Real rho_left{1.0};
    amrex::Real rho_right{1.0};
    amrex::Real rho_face{1.0};
    amrex::Real coefficient{0.0};
};

struct DiffusionResult {
    std::vector<amrex::Real> integrated_transfers;
    std::vector<amrex::Real> updated_state;
    amrex::Real maximum_outgoing_demand{0.0};
};

/**
 * Explicit two-point diffusion of every extensive component.
 *
 * The transported intensive variable is X/rho and the physical transfer is
 * I = -A dt rho_face K (X_R/rho_R-X_L/rho_L)/distance.  The caller supplies
 * the face coefficient K and physical geometry; no unit-area or unit-time
 * factors are inserted by this routine.
 */
[[nodiscard]] DiffusionResult
explicit_two_point_diffusion(const std::vector<amrex::Real>& state,
                             int ncell, int ncomp,
                             const std::vector<amrex::Real>& volumes,
                             const std::vector<DiffusionFace>& faces,
                             amrex::Real dt);

/** Conservative explicit bound based on total positive outgoing demand. */
[[nodiscard]] amrex::Real
admissible_explicit_timestep(const std::vector<amrex::Real>& state,
                             int ncell, int ncomp,
                             const std::vector<amrex::Real>& volumes,
                             const std::vector<DiffusionFace>& faces,
                             amrex::Real coefficient_scale = 1.0);

} // namespace erf_sbm

#endif
