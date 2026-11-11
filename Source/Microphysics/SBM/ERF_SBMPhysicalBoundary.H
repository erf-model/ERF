#ifndef ERF_SBM_PHYSICAL_BOUNDARY_H_
#define ERF_SBM_PHYSICAL_BOUNDARY_H_

#include "ERF_SBMAdvectionBoundary.H"
#include "ERF_SBMLayout.H"
#include "ERF_IndexDefines.H"

#include <AMReX_Orientation.H>
#include <AMReX_REAL.H>

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace erf_sbm {

inline constexpr int physical_face_count = 2 * AMREX_SPACEDIM;

/**
 * Spectral advection category and constant inflow spectrum of each domain
 * face, indexed by amrex::Orientation. Inflow values are dry-air-relative
 * z = U / rho_d for every layout component.
 */
struct SBMPhysicalBoundaryConfig {
    std::array<AdvectionFaceKind, physical_face_count> kind{};
    std::array<std::vector<amrex::Real>, physical_face_count> inflow_z;

    [[nodiscard]] bool all_periodic () const noexcept;
    [[nodiscard]] bool periodic_in (int dir) const noexcept;
    /** Restart-schema line; an all-periodic configuration keeps the M4b identity. */
    [[nodiscard]] std::string identity () const;
};

[[nodiscard]] const char* face_name (int face) noexcept;
[[nodiscard]] const char* face_kind_name (AdvectionFaceKind kind) noexcept;

/** Map one ERF physical boundary type, or explain why SBM rejects it. */
[[nodiscard]] bool map_erf_boundary (ERF_BC bc, AdvectionFaceKind& kind,
                                     std::string& reason);

/** Require a complete, canonical, group-admissible inflow spectrum. */
[[nodiscard]] bool validate_inflow_state (const SBMLayout& layout,
                                          const std::vector<amrex::Real>& z,
                                          std::string& reason);

/**
 * Cumulative composite outward transfer through each physical face: coarse
 * faces under a finer level are excluded, so the sum closes the leaf-domain
 * inventory. Gross is the sum of absolute face transfers, a roundoff scale.
 */
struct SBMBoundaryLedger {
    int ncomp{0};
    std::vector<double> outward;
    std::vector<double> gross;
    std::uint64_t face_stage_terms{0};

    void reset (int components);
    [[nodiscard]] double total_outward (int component) const;
    [[nodiscard]] double total_gross (int component) const;
    [[nodiscard]] std::string serialize () const;
    [[nodiscard]] static bool parse (const std::string& text, int components,
                                     SBMBoundaryLedger& ledger,
                                     std::string& reason);
};

} // namespace erf_sbm

#endif
