#include "ERF_SBMDiffusion.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace erf_sbm {

namespace {

void validate_common(const std::vector<amrex::Real>& state, const int ncell,
                     const int ncomp, const std::vector<amrex::Real>& volumes,
                     const std::vector<DiffusionFace>& faces)
{
    if (ncell <= 0 || ncomp <= 0 || state.size() != static_cast<std::size_t>(ncell*ncomp) ||
        volumes.size() != static_cast<std::size_t>(ncell)) {
        throw std::invalid_argument("invalid SBM diffusion state or volume array");
    }
    for (const auto volume : volumes) if (!(volume > 0.0) || !std::isfinite(volume)) {
        throw std::invalid_argument("SBM diffusion volumes must be finite and positive");
    }
    for (const auto& face : faces) {
        if (face.left_cell < 0 || face.left_cell >= ncell || face.right_cell < 0 || face.right_cell >= ncell ||
            face.left_cell == face.right_cell || !(face.area > 0.0) || !(face.distance > 0.0) ||
            !(face.left_volume > 0.0) || !(face.right_volume > 0.0) ||
            !(face.rho_left > 0.0) || !(face.rho_right > 0.0) || !(face.rho_face >= 0.0) ||
            !std::isfinite(face.area) || !std::isfinite(face.distance) ||
            !std::isfinite(face.coefficient)) {
            throw std::invalid_argument("invalid SBM diffusion face geometry or density");
        }
    }
}

} // namespace

DiffusionResult explicit_two_point_diffusion(const std::vector<amrex::Real>& state,
                                             const int ncell, const int ncomp,
                                             const std::vector<amrex::Real>& volumes,
                                             const std::vector<DiffusionFace>& faces,
                                             const amrex::Real dt)
{
    validate_common(state, ncell, ncomp, volumes, faces);
    if (!std::isfinite(dt) || dt < 0.0) throw std::invalid_argument("invalid SBM diffusion timestep");
    DiffusionResult result;
    result.integrated_transfers.assign(faces.size()*static_cast<std::size_t>(ncomp), 0.0);
    result.updated_state = state;
    for (std::size_t fi = 0; fi < faces.size(); ++fi) {
        const auto& face = faces[fi];
        const auto& left_volume = volumes[static_cast<std::size_t>(face.left_cell)];
        const auto& right_volume = volumes[static_cast<std::size_t>(face.right_cell)];
        const amrex::Real conductance = face.area * dt * face.rho_face * face.coefficient / face.distance;
        for (int c = 0; c < ncomp; ++c) {
            const auto left = state[static_cast<std::size_t>(face.left_cell*ncomp+c)];
            const auto right = state[static_cast<std::size_t>(face.right_cell*ncomp+c)];
            if (!std::isfinite(left) || !std::isfinite(right)) throw std::domain_error("nonfinite SBM diffusion state");
            const auto intensive_difference = right / face.rho_right - left / face.rho_left;
            const auto transfer = -conductance * intensive_difference;
            result.integrated_transfers[fi*static_cast<std::size_t>(ncomp)+static_cast<std::size_t>(c)] = transfer;
            result.updated_state[static_cast<std::size_t>(face.left_cell*ncomp+c)] -= transfer / left_volume;
            result.updated_state[static_cast<std::size_t>(face.right_cell*ncomp+c)] += transfer / right_volume;
            result.maximum_outgoing_demand = std::max(result.maximum_outgoing_demand,
                                                       std::max(amrex::Real(0.0), -transfer / left_volume));
            result.maximum_outgoing_demand = std::max(result.maximum_outgoing_demand,
                                                       std::max(amrex::Real(0.0), transfer / right_volume));
        }
    }
    return result;
}

amrex::Real admissible_explicit_timestep(const std::vector<amrex::Real>& state,
                                         const int ncell, const int ncomp,
                                         const std::vector<amrex::Real>& volumes,
                                         const std::vector<DiffusionFace>& faces,
                                         const amrex::Real coefficient_scale)
{
    if (!std::isfinite(coefficient_scale) || coefficient_scale < 0.0) {
        throw std::invalid_argument("invalid SBM diffusion coefficient scale");
    }
    validate_common(state, ncell, ncomp, volumes, faces);
    amrex::Real rate = 0.0;
    for (const auto& face : faces) {
        const auto conductance = face.area * face.rho_face * face.coefficient * coefficient_scale / face.distance;
        rate = std::max(rate, conductance / volumes[static_cast<std::size_t>(face.left_cell)] / face.rho_left);
        rate = std::max(rate, conductance / volumes[static_cast<std::size_t>(face.right_cell)] / face.rho_right);
    }
    return rate > 0.0 ? amrex::Real(0.5) / rate : std::numeric_limits<amrex::Real>::infinity();
}

} // namespace erf_sbm
