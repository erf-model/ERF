#include "ERF_SBMAMR.H"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace erf_sbm {

std::vector<amrex::Real> volume_weighted_restrict(const std::vector<amrex::Real>& fine_state,
                                                  const int nfine, const int ncoarse, const int ncomp,
                                                  const std::vector<int>& parent,
                                                  const std::vector<amrex::Real>& fine_volumes,
                                                  const std::vector<amrex::Real>& coarse_volumes)
{
    if (nfine <= 0 || ncoarse <= 0 || ncomp <= 0 || fine_state.size() != static_cast<std::size_t>(nfine*ncomp) ||
        parent.size() != static_cast<std::size_t>(nfine) || fine_volumes.size() != static_cast<std::size_t>(nfine) ||
        coarse_volumes.size() != static_cast<std::size_t>(ncoarse)) throw std::invalid_argument("invalid SBM restriction arrays");
    std::vector<amrex::Real> result(static_cast<std::size_t>(ncoarse*ncomp), 0.0);
    std::vector<amrex::Real> covered(static_cast<std::size_t>(ncoarse), 0.0);
    for (int f = 0; f < nfine; ++f) {
        const int c = parent[static_cast<std::size_t>(f)];
        if (c < 0 || c >= ncoarse || !(fine_volumes[static_cast<std::size_t>(f)] > 0.0)) throw std::invalid_argument("invalid SBM parent map");
        covered[static_cast<std::size_t>(c)] += fine_volumes[static_cast<std::size_t>(f)];
        for (int q = 0; q < ncomp; ++q) result[static_cast<std::size_t>(c*ncomp+q)] += fine_volumes[static_cast<std::size_t>(f)] * fine_state[static_cast<std::size_t>(f*ncomp+q)];
    }
    for (int c = 0; c < ncoarse; ++c) {
        if (!(coarse_volumes[static_cast<std::size_t>(c)] > 0.0) || covered[static_cast<std::size_t>(c)] > coarse_volumes[static_cast<std::size_t>(c)]*(1.0+1.0e-12)) throw std::invalid_argument("fine volume exceeds coarse volume");
        for (int q = 0; q < ncomp; ++q) result[static_cast<std::size_t>(c*ncomp+q)] /= coarse_volumes[static_cast<std::size_t>(c)];
    }
    return result;
}

std::vector<amrex::Real> piecewise_constant_prolong(const std::vector<amrex::Real>& coarse_state,
                                                    const int ncoarse, const int nfine, const int ncomp,
                                                    const std::vector<int>& parent)
{
    if (ncoarse <= 0 || nfine <= 0 || ncomp <= 0 || coarse_state.size() != static_cast<std::size_t>(ncoarse*ncomp) || parent.size() != static_cast<std::size_t>(nfine)) throw std::invalid_argument("invalid SBM prolongation arrays");
    std::vector<amrex::Real> result(static_cast<std::size_t>(nfine*ncomp));
    for (int f = 0; f < nfine; ++f) {
        const int c = parent[static_cast<std::size_t>(f)];
        if (c < 0 || c >= ncoarse) throw std::invalid_argument("invalid SBM prolongation parent");
        for (int q = 0; q < ncomp; ++q) result[static_cast<std::size_t>(f*ncomp+q)] = coarse_state[static_cast<std::size_t>(c*ncomp+q)];
    }
    return result;
}

amrex::Real register_flux_from_integrated_transfer(const amrex::Real integrated_transfer,
                                                   const amrex::Real face_area,
                                                   const amrex::Real dt)
{
    if (!std::isfinite(integrated_transfer) || !std::isfinite(face_area) || !std::isfinite(dt) || face_area <= 0.0 || dt <= 0.0) throw std::invalid_argument("physical transfer adapter requires positive face area and dt");
    return integrated_transfer / (face_area * dt);
}

PostRefluxCheck validate_post_reflux(const std::vector<amrex::Real>& pre_state,
                                     const std::vector<amrex::Real>& correction,
                                     const std::vector<amrex::Real>& post_state,
                                     const int level, const int ncell, const int ncomp,
                                     const std::vector<ConstraintGroup>& groups)
{
    if (level < 0 || ncell <= 0 || ncomp <= 0 || pre_state.size() != static_cast<std::size_t>(ncell*ncomp) || correction.size() != pre_state.size() || post_state.size() != pre_state.size()) throw std::invalid_argument("invalid post-reflux state");
    PostRefluxCheck result;
    result.level = level;
    for (int c = 0; c < ncell; ++c) {
        std::vector<amrex::Real> state(static_cast<std::size_t>(ncomp));
        for (int q = 0; q < ncomp; ++q) state[static_cast<std::size_t>(q)] = pre_state[static_cast<std::size_t>(c*ncomp+q)] + correction[static_cast<std::size_t>(c*ncomp+q)];
        for (const auto& group : groups) {
            amrex::Real margin = 0.0;
            std::string failed;
            if (!group.admissible(state, &margin, &failed)) {
                result.admissible = false; result.cell = c; result.group = group.semantic_id; result.constraint = failed; result.margin = margin;
                for (const auto& constraint : group.constraints) if (constraint.semantic_id == failed) {
                    result.pre_value = group.evaluate(constraint, std::vector<amrex::Real>(pre_state.begin()+c*ncomp, pre_state.begin()+(c+1)*ncomp));
                    result.post_value = group.evaluate(constraint, std::vector<amrex::Real>(post_state.begin()+c*ncomp, post_state.begin()+(c+1)*ncomp));
                    result.correction = result.post_value - result.pre_value;
                    break;
                }
                return result;
            }
            result.margin = std::min(result.margin == 0.0 ? margin : result.margin, margin);
        }
    }
    result.admissible = true;
    return result;
}

} // namespace erf_sbm
