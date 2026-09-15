#include "ERF_SBMBoundary.H"

#include <cmath>
#include <stdexcept>

namespace erf_sbm {

bool validate_prescribed_inflow(const BoundaryDescriptor& boundary,
                                const std::vector<ConstraintGroup>& groups,
                                std::string* failure)
{
    if (boundary.kind != BoundaryKind::PrescribedSpectralInflow) {
        if (failure) *failure = "boundary is not prescribed spectral inflow";
        return false;
    }
    for (const auto& group : groups) {
        amrex::Real margin = 0.0;
        std::string failed;
        if (!group.admissible(boundary.prescribed_state, &margin, &failed)) {
            if (failure) *failure = group.semantic_id + ":" + failed;
            return false;
        }
    }
    return true;
}

BoundaryTransferBudget make_boundary_budget(const BoundaryDescriptor& boundary,
                                            const std::vector<amrex::Real>& face_flux,
                                            const std::vector<ConstraintGroup>& groups,
                                            const std::vector<amrex::Real>& qc_projection,
                                            const std::vector<amrex::Real>& qr_projection,
                                            const amrex::Real area, const amrex::Real dt)
{
    if (!std::isfinite(area) || area < 0.0 || !std::isfinite(dt) || dt < 0.0 ||
        qc_projection.size() != qr_projection.size()) {
        throw std::invalid_argument("invalid SBM boundary transfer geometry");
    }
    if (boundary.kind == BoundaryKind::PrescribedSpectralInflow) {
        std::string failure;
        if (!validate_prescribed_inflow(boundary, groups, &failure)) {
            throw std::domain_error("invalid prescribed SBM inflow: " + failure);
        }
    }
    BoundaryTransferBudget result;
    result.auxiliary.resize(face_flux.size(), 0.0);
    if (boundary.kind == BoundaryKind::ImpermeableWall || boundary.kind == BoundaryKind::Periodic) return result;
    for (std::size_t i = 0; i < face_flux.size(); ++i) {
        if (!std::isfinite(face_flux[i])) throw std::invalid_argument("nonfinite SBM boundary flux");
        result.auxiliary[i] = area * dt * face_flux[i];
    }
    if (qc_projection.size() != face_flux.size() || qr_projection.size() != face_flux.size()) {
        throw std::invalid_argument("boundary bulk projection has the wrong number of components");
    }
    for (std::size_t i = 0; i < face_flux.size(); ++i) {
        result.qc += area * dt * qc_projection[i];
        result.qr += area * dt * qr_projection[i];
    }
    return result;
}

} // namespace erf_sbm

