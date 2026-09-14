#ifndef ERF_SBM_TRANSFER_CLOSURE_H_
#define ERF_SBM_TRANSFER_CLOSURE_H_

#include "ERF_AuxiliaryFaceTransfer.H"
#include "ERF_SBMLayout.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace erf_sbm {

/**
 * @brief Local closure residuals for one accepted auxiliary transport step.
 *
 * The compact inputs are ERF state MultiFabs and their component indices are
 * explicit because the production baseline is stored as a compact two-field
 * snapshot while the accepted state remains in the full ERF conserved bundle.
 */
struct AcceptedTransferClosure {
    amrex::Real spectral_max{0.0};
    amrex::Real qc_max{0.0};
    amrex::Real qr_max{0.0};
    amrex::Real spectral_tolerance{0.0};
    amrex::Real qc_tolerance{0.0};
    amrex::Real qr_tolerance{0.0};

    [[nodiscard]] bool passes() const noexcept {
        return spectral_max <= spectral_tolerance &&
               qc_max <= qc_tolerance && qr_max <= qr_tolerance;
    }
};

/**
 * Evaluate X_new - X_old + div(I) for every spectral component and for the
 * accepted compact cloud/rain face transfers.  The face arrays are the
 * already-accepted production objects; no flux is reconstructed here.
 */
AcceptedTransferClosure evaluate_accepted_transfer_closure(
    const SBMLayout& layout,
    const amrex::MultiFab& old_spectral,
    const amrex::MultiFab& new_spectral,
    const ::erf_auxiliary::AuxiliaryFaceTransfer& accepted_spectral,
    const amrex::MultiFab& old_compact,
    int old_qc_comp,
    int old_qr_comp,
    const amrex::MultiFab& new_compact,
    int new_qc_comp,
    int new_qr_comp,
    const ::erf_auxiliary::AuxiliaryFaceTransfer& accepted_compact,
    const amrex::Geometry& geometry);

} // namespace erf_sbm

#endif
