#ifndef ERF_SBM_AMR_H_
#define ERF_SBM_AMR_H_

#include "ERF_SBMConstraintGroups.H"

#include <AMReX_REAL.H>

#include <string>
#include <vector>

namespace erf_sbm {

/** Volume-weighted restriction for extensive auxiliary densities. */
[[nodiscard]] std::vector<amrex::Real>
volume_weighted_restrict(const std::vector<amrex::Real>& fine_state,
                         int nfine, int ncoarse, int ncomp,
                         const std::vector<int>& parent,
                         const std::vector<amrex::Real>& fine_volumes,
                         const std::vector<amrex::Real>& coarse_volumes);

/** Conservative piecewise-constant prolongation/injection. */
[[nodiscard]] std::vector<amrex::Real>
piecewise_constant_prolong(const std::vector<amrex::Real>& coarse_state,
                           int ncoarse, int nfine, int ncomp,
                           const std::vector<int>& parent);

/** YAFluxRegister expects a face flux per area and multiplies by dt/dx. */
[[nodiscard]] amrex::Real
register_flux_from_integrated_transfer(amrex::Real integrated_transfer,
                                       amrex::Real face_area,
                                       amrex::Real dt);

struct PostRefluxCheck {
    bool admissible{false};
    int level{0};
    int cell{0};
    std::string group;
    std::string constraint;
    amrex::Real pre_value{0.0};
    amrex::Real correction{0.0};
    amrex::Real post_value{0.0};
    amrex::Real margin{0.0};
};

/** Fail closed after reflux; this routine never clips or repairs a state. */
[[nodiscard]] PostRefluxCheck
validate_post_reflux(const std::vector<amrex::Real>& pre_state,
                     const std::vector<amrex::Real>& correction,
                     const std::vector<amrex::Real>& post_state,
                     int level, int ncell, int ncomp,
                     const std::vector<ConstraintGroup>& groups);

} // namespace erf_sbm

#endif
