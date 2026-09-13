#include "ERF_SBMTransportPrototype.H"

#include "ERF_IndexDefines.H"
#include "ERF_SBMBulkProjection.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MFParallelFor.H>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace erf_sbm {

namespace {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real donor_ratio (const amrex::Array4<const amrex::Real>& spectral,
                         const amrex::Array4<const amrex::Real>& density,
                         const int i, const int j, const int k, const int comp) noexcept
{
    const amrex::Real rho = density(i,j,k);
    return rho > amrex::Real(0.0) ? spectral(i,j,k,comp) / rho : amrex::Real(0.0);
}

} // namespace

void update_stage(const ::erf_auxiliary::StageContext& context,
                  const HostState& old_state,
                  const HostState& evaluation_state,
                  const HostState& rhs,
                  HostState& output_state)
{
    if (old_state.size() != rhs.size() || evaluation_state.size() != old_state.size()) {
        throw std::invalid_argument("auxiliary stage vectors have inconsistent sizes");
    }
    output_state.resize(old_state.size());
    if (context.method == ::erf_auxiliary::IntegrationMethod::CompressibleRK3) {
        // ERF uses the old full-step baseline for all three callback updates.
        for (std::size_t n = 0; n < old_state.size(); ++n) {
            output_state[n] = old_state[n] + static_cast<amrex::Real>(context.stage_interval) * rhs[n];
        }
    } else if (context.stage_index == 0) {
        for (std::size_t n = 0; n < old_state.size(); ++n) {
            output_state[n] = old_state[n] + static_cast<amrex::Real>(context.stage_interval) * rhs[n];
        }
    } else {
        // Anelastic Heun: old + 1/2[(predictor-old) + h R1].
        for (std::size_t n = 0; n < old_state.size(); ++n) {
            output_state[n] = old_state[n] + static_cast<amrex::Real>(0.5) *
                ((evaluation_state[n] - old_state[n]) +
                 static_cast<amrex::Real>(context.stage_interval) * rhs[n]);
        }
    }
}

HostState accepted_ledger(const ::erf_auxiliary::StageContext& context,
                          const std::vector<HostState>& stage_fluxes,
                          const double full_step)
{
    if ((context.method == ::erf_auxiliary::IntegrationMethod::CompressibleRK3 && stage_fluxes.size() != 3) ||
        (context.method == ::erf_auxiliary::IntegrationMethod::AnelasticHeun && stage_fluxes.size() != 2)) {
        throw std::invalid_argument("wrong number of stage fluxes for auxiliary ledger");
    }
    if (stage_fluxes.empty()) return {};
    HostState result(stage_fluxes.front().size(), amrex::Real(0.0));
    if (context.method == ::erf_auxiliary::IntegrationMethod::CompressibleRK3) {
        for (std::size_t n = 0; n < result.size(); ++n) {
            result[n] = context.completes_level_step ? static_cast<amrex::Real>(full_step) * stage_fluxes[2][n] : amrex::Real(0.0);
        }
    } else {
        for (std::size_t n = 0; n < result.size(); ++n) {
            result[n] = static_cast<amrex::Real>(0.5 * full_step) * (stage_fluxes[0][n] + stage_fluxes[1][n]);
        }
    }
    return result;
}

void advance_stage(::erf_auxiliary::AuxiliaryStateManager& manager,
                   const SBMLayout& layout,
                   const ::erf_auxiliary::StageContext& context,
                   const amrex::MultiFab& rho_evaluation,
                   amrex::MultiFab& core_state,
                   const amrex::MultiFab& carrier_x,
                   const amrex::MultiFab& carrier_y,
                   const amrex::MultiFab& carrier_z,
                   const amrex::Geometry& geometry)
{
    if (layout.populations().size() != 1) {
        throw std::invalid_argument("P1 transport supports exactly one spectral population");
    }
    if (geometry.isAllPeriodic() == false || geometry.Domain().length(0) <= 0) {
        throw std::invalid_argument("P1 transport requires a periodic Cartesian geometry");
    }
    const auto& population = layout.populations().front();
    const int first = population.liquid_mass_offset;
    const int nbins = population.grid.nbins();
    const int cloud_count = population.grid.cloud_bin_count();
    const auto& evaluation = (context.stage_index == 0) ? manager.old(0) : manager.evaluation(0);
    const auto& old = manager.old(0);
    const auto& predictor = manager.evaluation(0);
    auto& output = manager.output(0);
    const amrex::Real dt = static_cast<amrex::Real>(context.stage_interval);
    const amrex::Real dxi = static_cast<amrex::Real>(geometry.InvCellSize(0));
    const amrex::Real dyi = static_cast<amrex::Real>(geometry.InvCellSize(1));
    const amrex::Real dzi = static_cast<amrex::Real>(geometry.InvCellSize(2));

    // The accepted stage is formed in one kernel from the same face-transfer
    // expression used for every liquid bin.  No independently advected qc/qr
    // flux is constructed.
    for (amrex::MFIter mfi(output); mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.validbox();
        const auto eval = evaluation.const_array(mfi);
        const auto old_arr = old.const_array(mfi);
        const auto pred = predictor.const_array(mfi);
        const auto rho = rho_evaluation.const_array(mfi);
        const auto cx = carrier_x.const_array(mfi);
        const auto cy = carrier_y.const_array(mfi);
        const auto cz = carrier_z.const_array(mfi);
        const auto out = output.array(mfi);
        const auto core = core_state.array(mfi);

        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int b = 0; b < nbins; ++b) {
                const int n = first + b;
                const amrex::Real fx_lo = cx(i,j,k) * ((cx(i,j,k) >= amrex::Real(0.0)) ?
                    donor_ratio(eval, rho, i-1,j,k,n) : donor_ratio(eval, rho, i,j,k,n));
                const amrex::Real fx_hi = cx(i+1,j,k) * ((cx(i+1,j,k) >= amrex::Real(0.0)) ?
                    donor_ratio(eval, rho, i,j,k,n) : donor_ratio(eval, rho, i+1,j,k,n));
                const amrex::Real fy_lo = cy(i,j,k) * ((cy(i,j,k) >= amrex::Real(0.0)) ?
                    donor_ratio(eval, rho, i,j-1,k,n) : donor_ratio(eval, rho, i,j,k,n));
                const amrex::Real fy_hi = cy(i,j+1,k) * ((cy(i,j+1,k) >= amrex::Real(0.0)) ?
                    donor_ratio(eval, rho, i,j,k,n) : donor_ratio(eval, rho, i,j+1,k,n));
                const amrex::Real fz_lo = cz(i,j,k) * ((cz(i,j,k) >= amrex::Real(0.0)) ?
                    donor_ratio(eval, rho, i,j,k-1,n) : donor_ratio(eval, rho, i,j,k,n));
                const amrex::Real fz_hi = cz(i,j,k+1) * ((cz(i,j,k+1) >= amrex::Real(0.0)) ?
                    donor_ratio(eval, rho, i,j,k,n) : donor_ratio(eval, rho, i,j,k+1,n));
                const amrex::Real rhs = -((fx_hi-fx_lo)*dxi + (fy_hi-fy_lo)*dyi + (fz_hi-fz_lo)*dzi);
                if (context.method == ::erf_auxiliary::IntegrationMethod::CompressibleRK3 || context.stage_index == 0) {
                    out(i,j,k,n) = old_arr(i,j,k,n) + dt*rhs;
                } else {
                    out(i,j,k,n) = old_arr(i,j,k,n) + amrex::Real(0.5) *
                        ((pred(i,j,k,n)-old_arr(i,j,k,n)) + dt*rhs);
                }
            }
            amrex::Real qc = amrex::Real(0.0);
            amrex::Real qr = amrex::Real(0.0);
            for (int b = 0; b < cloud_count; ++b) qc += out(i,j,k,first+b);
            for (int b = cloud_count; b < nbins; ++b) qr += out(i,j,k,first+b);
            core(i,j,k,RhoQ2_comp) = qc;
            core(i,j,k,RhoQ3_comp) = qr;
        });
    }
    output.FillBoundary(geometry.periodicity());

    // P1 is deliberately fail-closed.  A material negative auxiliary state
    // must stop the run; clipping here would hide a donor/geometry/stage bug
    // and would make the compact projection disagree with its authoritative
    // spectral state.  The tolerance only covers roundoff at zero.
    for (int comp = 0; comp < layout.ncomp(); ++comp) {
        const amrex::Real minimum = output.min(comp);
        const amrex::Real scale = amrex::max(amrex::Real(1.0),
                                              amrex::Math::abs(minimum));
        const amrex::Real tolerance = amrex::Real(128.0) *
            std::numeric_limits<amrex::Real>::epsilon() * scale;
        if (!std::isfinite(minimum) || minimum < -tolerance) {
            throw std::runtime_error("SBM P1 auxiliary state violated the nonnegative invariant; no clipping is permitted");
        }
    }
    manager.accept_stage(0);
    core_state.FillBoundary(geometry.periodicity());
}

} // namespace erf_sbm
