#include "ERF_SBMTransportPrototype.H"

#include "ERF_IndexDefines.H"
#include "ERF_SBMBulkProjection.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MFParallelFor.H>

#include <cmath>
#include <limits>
#include <string>
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

void validate_nonnegative_state(const amrex::MultiFab& state, const int ncomp,
                                const char* context)
{
    if (ncomp <= 0 || ncomp > state.nComp()) {
        throw std::invalid_argument("invalid component count for auxiliary finite-state check");
    }
    if (!state.is_finite(0, ncomp, state.nGrowVect())) {
        throw std::runtime_error(std::string(context) + " contains NaN or infinite values");
    }
    for (int comp = 0; comp < ncomp; ++comp) {
        const amrex::Real minimum = state.min(comp);
        const amrex::Real scale = amrex::max(amrex::Real(1.0),
                                              amrex::Math::abs(minimum));
        const amrex::Real tolerance = amrex::Real(128.0) *
            std::numeric_limits<amrex::Real>::epsilon() * scale;
        if (minimum < -tolerance) {
            throw std::runtime_error(std::string(context) +
                                     " contains materially negative values; no clipping is permitted");
        }
    }
}

void advance_stage(::erf_auxiliary::AuxiliaryStateManager& manager,
                   const SBMLayout& layout,
                   const ::erf_auxiliary::StageContext& context,
                   const amrex::MultiFab& rho_evaluation,
                   amrex::MultiFab& core_state,
                   const amrex::MultiFab& carrier_x,
                   const amrex::MultiFab& carrier_y,
                   const amrex::MultiFab& carrier_z,
                   const amrex::Geometry& geometry,
                   ::erf_auxiliary::AuxiliaryFaceTransfer& stage_flux)
{
    if (layout.populations().size() != 1 ||
        layout.populations().front().moment_mode != MomentMode::OneMoment ||
        layout.ncomp() != layout.populations().front().grid.nbins()) {
        throw std::invalid_argument("P1 transport supports exactly one one-moment spectral population");
    }
    if (geometry.isAllPeriodic() == false || geometry.Domain().length(0) <= 0) {
        throw std::invalid_argument("P1 transport requires a periodic Cartesian geometry");
    }
    if (!stage_flux.defined() || stage_flux.ncomp() != layout.ncomp()) {
        throw std::invalid_argument("P1 stage face-transfer storage does not match the spectral layout");
    }
    const auto& population = layout.populations().front();
    const int first = population.mass_offset;
    const int nbins = population.grid.nbins();
    const auto& evaluation = (context.stage_index == 0) ? manager.old(0) : manager.evaluation(0);
    const auto& old = manager.old(0);
    const auto& predictor = manager.evaluation(0);
    auto& output = manager.output(0);
    const amrex::Real dt = static_cast<amrex::Real>(context.stage_interval);
    const amrex::Real dxi = static_cast<amrex::Real>(geometry.InvCellSize(0));
    const amrex::Real dyi = static_cast<amrex::Real>(geometry.InvCellSize(1));
    const amrex::Real dzi = static_cast<amrex::Real>(geometry.InvCellSize(2));

    stage_flux.setVal(amrex::Real(0.0));

    // Construct each numerical face flux exactly once on its face-centered
    // MultiFab.  The same arrays feed the divergence below and the accepted
    // full-step ledger; no cell-centered or independently reconstructed bulk
    // flux is involved.
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const auto* carrier = dir == 0 ? &carrier_x : (dir == 1 ? &carrier_y : &carrier_z);
        auto& flux = stage_flux.direction(dir);
        for (amrex::MFIter mfi(flux); mfi.isValid(); ++mfi) {
            const amrex::Box box = mfi.validbox();
            const auto carrier_arr = carrier->const_array(mfi);
            const auto eval = evaluation.const_array(mfi);
            const auto rho = rho_evaluation.const_array(mfi);
            const auto out = flux.array(mfi);
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                for (int b = 0; b < nbins; ++b) {
                    const amrex::Real face_mass_flux = carrier_arr(i,j,k);
                    const int donor_i = (dir == 0 && face_mass_flux >= amrex::Real(0.0)) ? i-1 : i;
                    const int donor_j = (dir == 1 && face_mass_flux >= amrex::Real(0.0)) ? j-1 : j;
                    const int donor_k = (dir == 2 && face_mass_flux >= amrex::Real(0.0)) ? k-1 : k;
                    out(i,j,k,first+b) = face_mass_flux *
                        donor_ratio(eval, rho, donor_i, donor_j, donor_k, first+b);
                }
            });
        }
    }

    for (amrex::MFIter mfi(output); mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.validbox();
        const auto old_arr = old.const_array(mfi);
        const auto pred = predictor.const_array(mfi);
        const auto fx = stage_flux.x().const_array(mfi);
        const auto fy = stage_flux.y().const_array(mfi);
        const auto fz = stage_flux.z().const_array(mfi);
        const auto out = output.array(mfi);

        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int b = 0; b < nbins; ++b) {
                const int n = first + b;
                const amrex::Real rhs = -((fx(i+1,j,k,n)-fx(i,j,k,n))*dxi +
                                           (fy(i,j+1,k,n)-fy(i,j,k,n))*dyi +
                                           (fz(i,j,k+1,n)-fz(i,j,k,n))*dzi);
                if (context.method == ::erf_auxiliary::IntegrationMethod::CompressibleRK3 || context.stage_index == 0) {
                    out(i,j,k,n) = old_arr(i,j,k,n) + dt*rhs;
                } else {
                    out(i,j,k,n) = old_arr(i,j,k,n) + amrex::Real(0.5) *
                        ((pred(i,j,k,n)-old_arr(i,j,k,n)) + dt*rhs);
                }
            }
        });
    }
    output.FillBoundary(geometry.periodicity());

    // The compact host fields are projections of the updated authoritative
    // spectral state.  This is deliberately separate from face-flux
    // construction: qc/qr never get their own numerical transport path.
    const SBMBulkProjection bulk_projection(layout);
    for (amrex::MFIter mfi(output); mfi.isValid(); ++mfi) {
        bulk_projection.apply_to_core(mfi.validbox(), output.const_array(mfi),
                                      core_state.array(mfi));
    }

    validate_nonnegative_state(output, layout.ncomp());
    manager.accept_stage(0);
    manager.record_stage_face_transfer(0, context, stage_flux);
    core_state.FillBoundary(geometry.periodicity());
}

} // namespace erf_sbm
