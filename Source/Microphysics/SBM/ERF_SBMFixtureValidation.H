#ifndef ERF_SBM_FIXTURE_VALIDATION_H_
#define ERF_SBM_FIXTURE_VALIDATION_H_

#include "ERF_SBMRemapping.H"

#include <vector>

namespace erf_sbm {

/** Result of the canonical-state check used for user/fixture initialization. */
struct FixtureInitialStateCanonicality {
    bool canonical{false};
    int population_id{-1};
    int bin{-1};
    MomentMode moment_mode{MomentMode::OneMoment};
};

/**
 * Apply the persisted-state admission rule to each bin before ERF accepts a
 * user/fixture initial spectrum. Kept as the constructor's shared validation
 * seam so unit tests can exercise the exact candidate-layout check.
 */
inline FixtureInitialStateCanonicality
validate_fixture_initial_state_canonicality (
    const SBMLayout& layout, const std::vector<amrex::Real>& state) noexcept
{
    FixtureInitialStateCanonicality result;
    if (state.size() != static_cast<std::size_t>(layout.ncomp())) return result;

    result.canonical = true;
    const int state_count = static_cast<int>(state.size());
    for (const auto& population : layout.populations()) {
        const auto view = population_remap_view(layout, population.population_id);
        for (int bin = 0; bin < population.grid.nbins(); ++bin) {
            if (!remap_detail::canonical_persisted_bin_state(
                    view, bin, state.data(), state_count)) {
                result.canonical = false;
                result.population_id = population.population_id;
                result.bin = bin;
                result.moment_mode = population.moment_mode;
                return result;
            }
        }
    }
    return result;
}

} // namespace erf_sbm

#endif // ERF_SBM_FIXTURE_VALIDATION_H_
