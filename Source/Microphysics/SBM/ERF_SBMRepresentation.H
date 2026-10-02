#ifndef ERF_SBM_REPRESENTATION_H_
#define ERF_SBM_REPRESENTATION_H_

#include "ERF_SpectralGrid.H"

namespace erf_sbm {

/** Stable scientific identities selected from a population's stored moments. */
struct RepresentationIdentity {
    const char* representation;
    const char* reconstruction;
    const char* packet_remap;
};

[[nodiscard]] constexpr RepresentationIdentity
representation_identity (const MomentMode mode) noexcept
{
    switch (mode) {
    case MomentMode::OneMoment:
        return {"fixed-pivot-1m-v1", "fixed-pivot-delta-v1", "fixed-pivot-two-center-v1"};
    case MomentMode::TwoMoment:
        return {"interval-2m-v1", "mean-delta-2m-v1", "actual-mass-interval-v1"};
    }
    return {"invalid-representation", "invalid-reconstruction", "invalid-packet-remap"};
}

} // namespace erf_sbm

#endif
