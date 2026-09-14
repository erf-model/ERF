#include "ERF_SBMTransferClosure.H"

#include <AMReX_MFParallelFor.H>

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace erf_sbm {

namespace {

struct ResidualMetric {
    amrex::Real maximum{0.0};
    amrex::Real scale{1.0};
};

amrex::Real scaled_tolerance(const amrex::Real scale)
{
    // The residual contains one state subtraction and six face terms.  A
    // 256-epsilon safety factor covers their accumulation while remaining a
    // machine-precision-scaled check rather than a broad physical tolerance.
    constexpr amrex::Real safety_factor = amrex::Real(256.0);
    return safety_factor * std::numeric_limits<amrex::Real>::epsilon() *
           std::max(amrex::Real(1.0), scale);
}

ResidualMetric measure_component(const amrex::MultiFab& old_state,
                                 const int old_comp,
                                 const amrex::MultiFab& new_state,
                                 const int new_comp,
                                 const ::erf_auxiliary::AuxiliaryFaceTransfer& transfer,
                                 const int transfer_comp,
                                 const amrex::Geometry& geometry)
{
    if (old_comp < 0 || new_comp < 0 || transfer_comp < 0 ||
        old_comp >= old_state.nComp() || new_comp >= new_state.nComp() ||
        transfer_comp >= transfer.ncomp()) {
        throw std::invalid_argument("accepted-transfer closure component is out of range");
    }
    amrex::MultiFab metric(old_state.boxArray(), old_state.DistributionMap(), 2, 0);
    metric.setVal(amrex::Real(0.0));
    const amrex::Real dxi = static_cast<amrex::Real>(geometry.InvCellSize(0));
    const amrex::Real dyi = static_cast<amrex::Real>(geometry.InvCellSize(1));
    const amrex::Real dzi = static_cast<amrex::Real>(geometry.InvCellSize(2));
    for (amrex::MFIter mfi(metric); mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.validbox();
        const auto old_arr = old_state.const_array(mfi);
        const auto new_arr = new_state.const_array(mfi);
        const auto fx = transfer.x().const_array(mfi);
        const auto fy = transfer.y().const_array(mfi);
        const auto fz = transfer.z().const_array(mfi);
        const auto out = metric.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const amrex::Real divergence =
                (fx(i+1,j,k,transfer_comp) - fx(i,j,k,transfer_comp)) * dxi +
                (fy(i,j+1,k,transfer_comp) - fy(i,j,k,transfer_comp)) * dyi +
                (fz(i,j,k+1,transfer_comp) - fz(i,j,k,transfer_comp)) * dzi;
            const amrex::Real residual = new_arr(i,j,k,new_comp) -
                old_arr(i,j,k,old_comp) + divergence;
            const amrex::Real scale =
                amrex::Math::abs(new_arr(i,j,k,new_comp)) +
                amrex::Math::abs(old_arr(i,j,k,old_comp)) +
                (amrex::Math::abs(fx(i+1,j,k,transfer_comp)) +
                 amrex::Math::abs(fx(i,j,k,transfer_comp))) * dxi +
                (amrex::Math::abs(fy(i,j+1,k,transfer_comp)) +
                 amrex::Math::abs(fy(i,j,k,transfer_comp))) * dyi +
                (amrex::Math::abs(fz(i,j,k+1,transfer_comp)) +
                 amrex::Math::abs(fz(i,j,k,transfer_comp))) * dzi +
                amrex::Real(1.0);
            out(i,j,k,0) = amrex::Math::abs(residual);
            out(i,j,k,1) = scale;
        });
    }
    return {metric.max(0), metric.max(1)};
}

ResidualMetric measure_spectral(const SBMLayout& layout,
                                const amrex::MultiFab& old_spectral,
                                const amrex::MultiFab& new_spectral,
                                const ::erf_auxiliary::AuxiliaryFaceTransfer& transfer,
                                const amrex::Geometry& geometry)
{
    if (old_spectral.nComp() < layout.ncomp() || new_spectral.nComp() < layout.ncomp() ||
        transfer.ncomp() < layout.ncomp()) {
        throw std::invalid_argument("accepted-transfer closure spectral layout is incompatible");
    }
    ResidualMetric result;
    for (int comp = 0; comp < layout.ncomp(); ++comp) {
        const auto component = measure_component(old_spectral, comp, new_spectral, comp,
                                                 transfer, comp, geometry);
        result.maximum = std::max(result.maximum, component.maximum);
        result.scale = std::max(result.scale, component.scale);
    }
    return result;
}

} // namespace

AcceptedTransferClosure evaluate_accepted_transfer_closure(
    const SBMLayout& layout,
    const amrex::MultiFab& old_spectral,
    const amrex::MultiFab& new_spectral,
    const ::erf_auxiliary::AuxiliaryFaceTransfer& accepted_spectral,
    const amrex::MultiFab& old_compact,
    const int old_qc_comp,
    const int old_qr_comp,
    const amrex::MultiFab& new_compact,
    const int new_qc_comp,
    const int new_qr_comp,
    const ::erf_auxiliary::AuxiliaryFaceTransfer& accepted_compact,
    const amrex::Geometry& geometry)
{
    const auto spectral = measure_spectral(layout, old_spectral, new_spectral,
                                            accepted_spectral, geometry);
    const auto qc = measure_component(old_compact, old_qc_comp, new_compact, new_qc_comp,
                                      accepted_compact, 0, geometry);
    const auto qr = measure_component(old_compact, old_qr_comp, new_compact, new_qr_comp,
                                      accepted_compact, 1, geometry);
    return {spectral.maximum, qc.maximum, qr.maximum,
            scaled_tolerance(spectral.scale), scaled_tolerance(qc.scale),
            scaled_tolerance(qr.scale)};
}

} // namespace erf_sbm
