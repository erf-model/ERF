#ifndef ERF_SBM_BOUNDARY_H_
#define ERF_SBM_BOUNDARY_H_

#include "ERF_SBMConstraintGroups.H"

#include <AMReX_REAL.H>

#include <string>
#include <vector>

namespace erf_sbm {

enum class BoundaryKind { Periodic, PrescribedSpectralInflow, AdvectiveOutflow, ImpermeableWall };

struct BoundaryDescriptor {
    BoundaryKind kind{BoundaryKind::Periodic};
    int face_id{0};
    std::string semantic_id;
    std::vector<amrex::Real> prescribed_state;
};

struct BoundaryTransferBudget {
    std::vector<amrex::Real> auxiliary;
    amrex::Real qc{0.0};
    amrex::Real qr{0.0};
};

[[nodiscard]] bool validate_prescribed_inflow(const BoundaryDescriptor& boundary,
                                              const std::vector<ConstraintGroup>& groups,
                                              std::string* failure = nullptr);

/** Convert a boundary face flux into the same signed integrated-transfer ledger. */
[[nodiscard]] BoundaryTransferBudget
make_boundary_budget(const BoundaryDescriptor& boundary,
                     const std::vector<amrex::Real>& face_flux,
                     const std::vector<ConstraintGroup>& groups,
                     const std::vector<amrex::Real>& qc_projection,
                     const std::vector<amrex::Real>& qr_projection,
                     amrex::Real area, amrex::Real dt);

} // namespace erf_sbm

#endif
