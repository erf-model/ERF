#include "ERF_SBMTransportPrototype.H"

#include "ERF_IndexDefines.H"
#include "ERF_SBMBulkProjection.H"
#include "ERF_Interpolation_WENO_Z.H"

#include <AMReX_ParReduce.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MFParallelFor.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_GpuUtility.H>
#include <AMReX_GpuContainers.H>

#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <stdexcept>

namespace erf_sbm {

using amrex::Real;

namespace {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real donor_ratio (const amrex::Array4<const amrex::Real>& spectral,
                         const amrex::Array4<const amrex::Real>& density,
                         const int i, const int j, const int k, const int comp) noexcept
{
    const amrex::Real rho = density(i,j,k);
    return rho > amrex::Real(0.0) ? spectral(i,j,k,comp) / rho : amrex::Real(0.0);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void restrict_constraint (const amrex::Real delta_form,
                          const amrex::Real left_margin,
                          const amrex::Real right_margin,
                          const amrex::Real scale,
                          amrex::Real& lambda) noexcept
{
    const amrex::Real left_demand = -scale * delta_form;
    const amrex::Real right_demand = scale * delta_form;
    if (left_demand > amrex::Real(0.0)) {
        lambda = amrex::min(lambda, amrex::max(amrex::Real(0.0), left_margin / left_demand));
    }
    if (right_demand > amrex::Real(0.0)) {
        lambda = amrex::min(lambda, amrex::max(amrex::Real(0.0), right_margin / right_demand));
    }
}

struct ProductionProperty {
    int component{0};
    int kind{0};
    int has_upper{0};
    amrex::Real support_min{0.0};
    amrex::Real support_max{0.0};
    amrex::Real pivot{1.0};
};

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

    // One local GPU-capable traversal covers every requested component and
    // returns the complete validity tuple.  ParReduce is local in AMReX, so
    // the three fixed-size reductions below are the only MPI collectives;
    // their count is independent of ncomp.
    const auto& arrays = state.const_arrays();
    const auto local = amrex::ParReduce(
        amrex::TypeList<amrex::ReduceOpLogicalOr, amrex::ReduceOpMin,
                        amrex::ReduceOpMax>{},
        amrex::TypeList<int, amrex::Real, amrex::Real>{},
        state, state.nGrowVect(), ncomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int comp)
            -> amrex::GpuTuple<int, amrex::Real, amrex::Real> {
            const amrex::Real value = arrays[box_no](i,j,k,comp);
            const int nonfinite = (amrex::isnan(value) || amrex::isinf(value)) ? 1 : 0;
            return {nonfinite, value, amrex::Math::abs(value)};
        });
    int has_nonfinite = amrex::get<0>(local);
    amrex::Real minimum = amrex::get<1>(local);
    amrex::Real maximum_absolute = amrex::max(amrex::Real(0.0), amrex::get<2>(local));
    amrex::ParallelDescriptor::ReduceIntMax(has_nonfinite);
    amrex::ParallelDescriptor::ReduceRealMin(minimum);
    amrex::ParallelDescriptor::ReduceRealMax(maximum_absolute);
    if (has_nonfinite != 0) {
        throw std::runtime_error(std::string(context) + " contains NaN or infinite values");
    }

    // Scale only with the global maximum magnitude.  There is deliberately
    // no order-one floor: an all-zero state has tau_neg == 0.
    const amrex::Real tolerance = maximum_absolute == amrex::Real(0.0) ?
        amrex::Real(0.0) : amrex::Real(128.0) *
        std::numeric_limits<amrex::Real>::epsilon() * maximum_absolute;
    if (minimum < -tolerance) {
        std::ostringstream message;
        message << context << " has material negative value: minimum=" << minimum
                << ", global_max_abs=" << maximum_absolute
                << ", tolerance=" << tolerance;
        throw std::runtime_error(message.str());
    }
}

void validate_admissible_state(::erf_auxiliary::AuxiliaryStateManager& manager,
                               const SBMLayout& layout, const int level)
{
    if (!manager.has_level(level)) throw std::invalid_argument("SBM admissibility check references an undefined level");
    validate_nonnegative_state(manager.output(level), layout.ncomp(), "SBM post-reflux physical state");
    if (layout.populations().size() != 1) {
        throw std::invalid_argument("SBM admissibility check requires the current single-population P2 production contract");
    }
    const auto& population = layout.populations().front();
    if (population.moment_mode != MomentMode::TwoMoment) return;
    auto& scratch = manager.scratch(level);
    const auto& source = manager.output(level);
    amrex::MultiFab::Copy(scratch, source, 0, 0, layout.ncomp(), source.nGrowVect());
    const int first = population.mass_offset;
    const int nbins = population.grid.nbins();
    for (amrex::MFIter mfi(scratch); mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.validbox();
        const auto state = source.const_array(mfi);
        const auto endpoints = scratch.array(mfi);
        for (int b = 0; b < nbins; ++b) {
            const int mass = first + b;
            const int number = population.number_offset + b;
            const Real lower = population.grid.edges()[static_cast<std::size_t>(b)];
            const Real upper = population.grid.edges()[static_cast<std::size_t>(b+1)];
            const Real denominator = upper - lower;
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const Real M = state(i,j,k,mass);
                const Real C = state(i,j,k,number);
                const Real scale = amrex::Math::abs(M) + upper*amrex::Math::abs(C) + lower*amrex::Math::abs(C);
                const Real tolerance = Real(128.0) * std::numeric_limits<Real>::epsilon() * scale;
                Real L = (upper*C-M)/denominator;
                Real H = (M-lower*C)/denominator;
                if (L < Real(0.0) && L >= -tolerance) L = Real(0.0);
                if (H < Real(0.0) && H >= -tolerance) H = Real(0.0);
                endpoints(i,j,k,mass) = L;
                endpoints(i,j,k,number) = H;
                if (L < Real(0.0) || H < Real(0.0)) {
                    endpoints(i,j,k,mass) = -Real(1.0);
                    endpoints(i,j,k,number) = -Real(1.0);
                }
            });
        }
    }
    validate_nonnegative_state(scratch, layout.ncomp(), "SBM post-reflux endpoint state");
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
                   ::erf_auxiliary::AuxiliaryFaceTransfer& stage_flux,
                   const TransportMethod method,
                   const int level,
                   const amrex::Real diffusion_coefficient)
{
    if (layout.populations().size() != 1 || !manager.has_level(level)) {
        throw std::invalid_argument("ERF SBM transport requires one initialized runtime population and level");
    }
    if (geometry.isAllPeriodic() == false || geometry.Domain().length(0) <= 0) {
        throw std::invalid_argument("P1 transport requires a periodic Cartesian geometry");
    }
    if (!stage_flux.defined() || stage_flux.ncomp() != layout.ncomp()) {
        throw std::invalid_argument("P1 stage face-transfer storage does not match the spectral layout");
    }
    if (!std::isfinite(diffusion_coefficient) || diffusion_coefficient < amrex::Real(0.0)) {
        throw std::invalid_argument("SBM explicit diffusion coefficient must be finite and nonnegative");
    }
    if (method == TransportMethod::GroupedFCT_WENOZ3 && manager.scratch(level).nGrowVect().min() < 2) {
        throw std::invalid_argument("GroupedFCT_WENOZ3 requires two ghost cells for ERF WENO_Z3");
    }
    const auto& population = layout.populations().front();
    if (layout.ncomp() < population.component_count) {
        throw std::invalid_argument("ERF SBM transport layout is smaller than its population storage");
    }
    const int first = population.mass_offset;
    const int nbins = population.grid.nbins();
    const bool two_moment = population.moment_mode == MomentMode::TwoMoment;
    const int ncomp = layout.ncomp();
    const auto& evaluation = (context.stage_index == 0) ? manager.old(level) : manager.evaluation(level);
    const auto& old = manager.old(level);
    const auto& predictor = manager.evaluation(level);
    auto& output = manager.output(level);
    auto& transport_scratch = manager.scratch(level);
    const amrex::Real dt = static_cast<amrex::Real>(context.stage_interval);
    const amrex::Real dxi = static_cast<amrex::Real>(geometry.InvCellSize(0));
    const amrex::Real dyi = static_cast<amrex::Real>(geometry.InvCellSize(1));
    const amrex::Real dzi = static_cast<amrex::Real>(geometry.InvCellSize(2));
    if (!std::isfinite(dt) || dt < Real(0.0)) {
        throw std::invalid_argument("SBM stage interval must be finite and nonnegative");
    }
    const Real diffusion_bound = dt * diffusion_coefficient * (dxi*dxi + dyi*dyi + dzi*dzi);
    if (diffusion_bound > Real(0.5)) {
        std::ostringstream message;
        message << "SBM explicit diffusion timestep exceeds admissible bound: dt=" << dt
                << " bound=" << (diffusion_coefficient > Real(0.0) ?
                    Real(0.5) / (diffusion_coefficient * (dxi*dxi + dyi*dyi + dzi*dzi)) :
                    std::numeric_limits<Real>::infinity())
                << " diffusion_coefficient=" << diffusion_coefficient;
        throw std::invalid_argument(message.str());
    }

    stage_flux.setVal(amrex::Real(0.0));

    // Two-moment storage is (M,C), while transport uses nonnegative endpoint
    // variables (L,H).  The existing per-level scratch FAB holds endpoint
    // ratios, so the temporary is bounded by the local tile and never grows
    // with a compile-time MAX_BINS constant.
    if (two_moment) {
        for (amrex::MFIter mfi(transport_scratch); mfi.isValid(); ++mfi) {
            const amrex::Box box = mfi.validbox();
            const auto source = evaluation.const_array(mfi);
            const auto rho = rho_evaluation.const_array(mfi);
            const auto scratch = transport_scratch.array(mfi);
            ParallelFor(box, layout.ncomp(), [=] AMREX_GPU_DEVICE (int i, int j, int k, int c) noexcept {
                const Real density = rho(i,j,k);
                scratch(i,j,k,c) = density > Real(0.0) ? source(i,j,k,c)/density : Real(0.0);
            });
            for (int b = 0; b < nbins; ++b) {
                const int mass = first + b;
                const int number = population.number_offset + b;
                const Real lower = population.grid.edges()[static_cast<std::size_t>(b)];
                const Real upper = population.grid.edges()[static_cast<std::size_t>(b+1)];
                const Real denominator = upper - lower;
                ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    const Real M = source(i,j,k,mass);
                    const Real C = source(i,j,k,number);
                    const Real scale = amrex::Math::abs(M) + upper*amrex::Math::abs(C) + lower*amrex::Math::abs(C);
                    const Real tolerance = Real(128.0) * std::numeric_limits<Real>::epsilon() * scale;
                    Real L = (upper*C - M) / denominator;
                    Real H = (M - lower*C) / denominator;
                    if (L < Real(0.0) && L >= -tolerance) L = Real(0.0);
                    if (H < Real(0.0) && H >= -tolerance) H = Real(0.0);
                    const Real density = rho(i,j,k);
                    scratch(i,j,k,mass) = density > Real(0.0) ? L/density : Real(0.0);
                    scratch(i,j,k,number) = density > Real(0.0) ? H/density : Real(0.0);
                });
            }
        }
        transport_scratch.FillBoundary(geometry.periodicity());
    } else if (method == TransportMethod::GroupedFCT_WENOZ3) {
        // WENO reconstructs the intensive ratio X/rho.  Keep that ratio in
        // the same bounded scratch FAB used by the endpoint path.
        for (amrex::MFIter mfi(transport_scratch); mfi.isValid(); ++mfi) {
            const amrex::Box box = mfi.validbox();
            const auto source = evaluation.const_array(mfi);
            const auto rho = rho_evaluation.const_array(mfi);
            const auto scratch = transport_scratch.array(mfi);
            ParallelFor(box, layout.ncomp(), [=] AMREX_GPU_DEVICE (int i, int j, int k, int c) noexcept {
                const Real density = rho(i,j,k);
                scratch(i,j,k,c) = density > Real(0.0) ? source(i,j,k,c)/density : Real(0.0);
            });
        }
        transport_scratch.FillBoundary(geometry.periodicity());
    }

    std::unique_ptr<::erf_auxiliary::AuxiliaryFaceTransfer> high_flux;
    if (method == TransportMethod::GroupedFCT_WENOZ3) {
        high_flux = std::make_unique<::erf_auxiliary::AuxiliaryFaceTransfer>();
        high_flux->define(stage_flux.x().boxArray(), stage_flux.x().DistributionMap(), layout.ncomp(), 0);
    }

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
            const auto eval = two_moment ? transport_scratch.const_array(mfi) : evaluation.const_array(mfi);
            const auto physical = evaluation.const_array(mfi);
            const auto rho = rho_evaluation.const_array(mfi);
            const auto out = flux.array(mfi);
            const Real inverse_distance = geometry.InvCellSize(dir);
            for (int b = 0; b < nbins; ++b) {
                const int mass = first + b;
                const int number = two_moment ? population.number_offset + b : -1;
                const Real lower = two_moment ? population.grid.edges()[static_cast<std::size_t>(b)] : Real(0.0);
                const Real upper = two_moment ? population.grid.edges()[static_cast<std::size_t>(b+1)] : Real(0.0);
                ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    const amrex::Real face_mass_flux = carrier_arr(i,j,k);
                    const int donor_i = (dir == 0 && face_mass_flux >= amrex::Real(0.0)) ? i-1 : i;
                    const int donor_j = (dir == 1 && face_mass_flux >= amrex::Real(0.0)) ? j-1 : j;
                    const int donor_k = (dir == 2 && face_mass_flux >= amrex::Real(0.0)) ? k-1 : k;
                    if (two_moment) {
                        const Real left_endpoint = face_mass_flux * eval(donor_i,donor_j,donor_k,mass);
                        const Real right_endpoint = face_mass_flux * eval(donor_i,donor_j,donor_k,number);
                        out(i,j,k,mass) = lower*left_endpoint + upper*right_endpoint;
                        out(i,j,k,number) = left_endpoint + right_endpoint;
                    } else {
                        out(i,j,k,mass) = face_mass_flux *
                            donor_ratio(eval, rho, donor_i, donor_j, donor_k, mass);
                    }
                    if (diffusion_coefficient > Real(0.0)) {
                        const int left_i = dir == 0 ? i-1 : i;
                        const int left_j = dir == 1 ? j-1 : j;
                        const int left_k = dir == 2 ? k-1 : k;
                        const int right_i = i, right_j = j, right_k = k;
                        const Real rho_left = rho(left_i,left_j,left_k);
                        const Real rho_right = rho(right_i,right_j,right_k);
                        const Real rho_face = Real(0.5) * (rho_left + rho_right);
                        if (rho_left > Real(0.0) && rho_right > Real(0.0) && rho_face > Real(0.0)) {
                            if (two_moment) {
                                const Real ldiff = -rho_face * diffusion_coefficient *
                                    (eval(right_i,right_j,right_k,mass) - eval(left_i,left_j,left_k,mass)) * inverse_distance;
                                const Real hdiff = -rho_face * diffusion_coefficient *
                                    (eval(right_i,right_j,right_k,number) - eval(left_i,left_j,left_k,number)) * inverse_distance;
                                out(i,j,k,mass) += lower*ldiff + upper*hdiff;
                                out(i,j,k,number) += ldiff + hdiff;
                            } else {
                                out(i,j,k,mass) += -rho_face * diffusion_coefficient *
                                    (donor_ratio(physical, rho, right_i,right_j,right_k,mass) -
                                     donor_ratio(physical, rho, left_i,left_j,left_k,mass)) * inverse_distance;
                            }
                        }
                    }
                });
            }
            // Attached carrier-bin properties share the same donor mass flux
            // and density-weighted intensive transport as the carrier state.
            for (int c = 0; c < ncomp; ++c) {
                const bool is_mass = c >= first && c < first + nbins;
                const bool is_number = two_moment && c >= population.number_offset &&
                                       c < population.number_offset + nbins;
                if (is_mass || is_number) continue;
                ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    const Real face_mass_flux = carrier_arr(i,j,k);
                    const int donor_i = (dir == 0 && face_mass_flux >= Real(0.0)) ? i-1 : i;
                    const int donor_j = (dir == 1 && face_mass_flux >= Real(0.0)) ? j-1 : j;
                    const int donor_k = (dir == 2 && face_mass_flux >= Real(0.0)) ? k-1 : k;
                    const Real rho_left = rho(dir == 0 ? i-1 : i,
                                              dir == 1 ? j-1 : j,
                                              dir == 2 ? k-1 : k);
                    const Real rho_right = rho(i,j,k);
                    out(i,j,k,c) = face_mass_flux * donor_ratio(physical, rho, donor_i, donor_j, donor_k, c);
                    if (diffusion_coefficient > Real(0.0) && rho_left > Real(0.0) && rho_right > Real(0.0)) {
                        const Real rho_face = Real(0.5) * (rho_left + rho_right);
                        const int left_i = dir == 0 ? i-1 : i;
                        const int left_j = dir == 1 ? j-1 : j;
                        const int left_k = dir == 2 ? k-1 : k;
                        out(i,j,k,c) += -rho_face * diffusion_coefficient *
                            (physical(i,j,k,c)/rho_right - physical(left_i,left_j,left_k,c)/rho_left) *
                            inverse_distance;
                    }
                });
            }
        }
    }

    if (high_flux) {
        // The high-order candidate shares the exact carrier mass flux and the
        // same reusable WENO_Z3 helper used by ERF's scalar advection path.
        // Only the candidate differs from the donor low-order transfer.
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            const auto* carrier = dir == 0 ? &carrier_x : (dir == 1 ? &carrier_y : &carrier_z);
            auto& flux = high_flux->direction(dir);
            for (amrex::MFIter mfi(flux); mfi.isValid(); ++mfi) {
                const amrex::Box box = mfi.validbox();
                const auto carrier_arr = carrier->const_array(mfi);
                const auto ratio = transport_scratch.const_array(mfi);
                const auto out = flux.array(mfi);
                const auto density = rho_evaluation.const_array(mfi);
                const auto physical = evaluation.const_array(mfi);
                const Real inverse_distance = geometry.InvCellSize(dir);
                for (int b = 0; b < nbins; ++b) {
                    const int mass = first + b;
                    const int number = two_moment ? population.number_offset + b : -1;
                    const Real lower = two_moment ? population.grid.edges()[static_cast<std::size_t>(b)] : Real(0.0);
                    const Real upper = two_moment ? population.grid.edges()[static_cast<std::size_t>(b+1)] : Real(0.0);
                    ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        const Real mass_flux = carrier_arr(i,j,k);
                        WENO_Z3 weno(ratio, Real(0.0));
                        Real lo = 0.0;
                        if (dir == 0) weno.InterpolateInX(i,j,k,mass,lo,mass_flux);
                        else if (dir == 1) weno.InterpolateInY(i,j,k,mass,lo,mass_flux);
                        else weno.InterpolateInZ(i,j,k,mass,lo,mass_flux);
                        if (two_moment) {
                            Real hi = 0.0;
                            if (dir == 0) weno.InterpolateInX(i,j,k,number,hi,mass_flux);
                            else if (dir == 1) weno.InterpolateInY(i,j,k,number,hi,mass_flux);
                            else weno.InterpolateInZ(i,j,k,number,hi,mass_flux);
                            const Real left_flux = mass_flux * lo;
                            const Real right_flux = mass_flux * hi;
                            out(i,j,k,mass) = lower*left_flux + upper*right_flux;
                            out(i,j,k,number) = left_flux + right_flux;
                        } else {
                            out(i,j,k,mass) = mass_flux * lo;
                        }
                        if (diffusion_coefficient > Real(0.0)) {
                            const int left_i = dir == 0 ? i-1 : i;
                            const int left_j = dir == 1 ? j-1 : j;
                            const int left_k = dir == 2 ? k-1 : k;
                            const Real rho_left = density(left_i,left_j,left_k);
                            const Real rho_right = density(i,j,k);
                            const Real rho_face = Real(0.5) * (rho_left + rho_right);
                            if (rho_left > Real(0.0) && rho_right > Real(0.0) && rho_face > Real(0.0)) {
                                if (two_moment) {
                                    // The endpoint reconstructed values are
                                    // in lo/hi; reconstruct the low-side
                                    // diffusion gradient from the endpoint
                                    // ratios in the neighboring cells.
                                    const Real ldiff = -rho_face * diffusion_coefficient *
                                        (ratio(i,j,k,mass) - ratio(left_i,left_j,left_k,mass)) * inverse_distance;
                                    const Real hdiff = -rho_face * diffusion_coefficient *
                                        (ratio(i,j,k,number) - ratio(left_i,left_j,left_k,number)) * inverse_distance;
                                    out(i,j,k,mass) += lower*ldiff + upper*hdiff;
                                    out(i,j,k,number) += ldiff + hdiff;
                                } else {
                                    out(i,j,k,mass) += -rho_face * diffusion_coefficient *
                                        (physical(i,j,k,mass)/rho_right - physical(left_i,left_j,left_k,mass)/rho_left) * inverse_distance;
                                }
                            }
                        }
                    });
                }
                for (int c = 0; c < ncomp; ++c) {
                    const bool is_mass = c >= first && c < first + nbins;
                    const bool is_number = two_moment && c >= population.number_offset &&
                                           c < population.number_offset + nbins;
                    if (is_mass || is_number) continue;
                    ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        const Real mass_flux = carrier_arr(i,j,k);
                        WENO_Z3 weno(ratio, Real(0.0));
                        Real value = Real(0.0);
                        if (dir == 0) weno.InterpolateInX(i,j,k,c,value,mass_flux);
                        else if (dir == 1) weno.InterpolateInY(i,j,k,c,value,mass_flux);
                        else weno.InterpolateInZ(i,j,k,c,value,mass_flux);
                        out(i,j,k,c) = mass_flux * value;
                        if (diffusion_coefficient > Real(0.0)) {
                            const int left_i = dir == 0 ? i-1 : i;
                            const int left_j = dir == 1 ? j-1 : j;
                            const int left_k = dir == 2 ? k-1 : k;
                            const Real rho_left = density(left_i,left_j,left_k);
                            const Real rho_right = density(i,j,k);
                            const Real rho_face = Real(0.5) * (rho_left + rho_right);
                            if (rho_left > Real(0.0) && rho_right > Real(0.0) && rho_face > Real(0.0)) {
                                out(i,j,k,c) += -rho_face * diffusion_coefficient *
                                    (physical(i,j,k,c)/rho_right - physical(left_i,left_j,left_k,c)/rho_left) * inverse_distance;
                            }
                        }
                    });
                }
            }
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
                if (two_moment) {
                    const int number = population.number_offset + b;
                    const amrex::Real rhs_number = -((fx(i+1,j,k,number)-fx(i,j,k,number))*dxi +
                                                      (fy(i,j+1,k,number)-fy(i,j,k,number))*dyi +
                                                      (fz(i,j,k+1,number)-fz(i,j,k,number))*dzi);
                    if (context.method == ::erf_auxiliary::IntegrationMethod::CompressibleRK3 || context.stage_index == 0) {
                        out(i,j,k,number) = old_arr(i,j,k,number) + dt*rhs_number;
                    } else {
                        out(i,j,k,number) = old_arr(i,j,k,number) + amrex::Real(0.5) *
                            ((pred(i,j,k,number)-old_arr(i,j,k,number)) + dt*rhs_number);
                    }
                }
            }
        });
    }
    output.FillBoundary(geometry.periodicity());

    if (high_flux) {
        // Apply one lambda per shared face to the complete population/bin
        // group.  All group constraints are evaluated in physical storage;
        // two-moment endpoint constraints are expanded algebraically so the
        // device path never captures runtime STL containers.
        const Real anelastic_weight =
            (context.method == ::erf_auxiliary::IntegrationMethod::AnelasticHeun && context.stage_index > 0)
            ? Real(0.5) : Real(1.0);
        amrex::Gpu::ManagedVector<ProductionProperty> properties;
        for (std::size_t p = 0; p < layout.attached_properties().size(); ++p) {
            const auto& property = layout.attached_properties()[p];
            if (property.carrier_population != population.population_id) continue;
            ProductionProperty meta;
            meta.component = layout.property_offset(static_cast<int>(p));
            meta.kind = static_cast<int>(property.kind);
            meta.has_upper = std::isfinite(property.support_max) ? 1 : 0;
            meta.support_min = property.support_min;
            meta.support_max = property.support_max;
            properties.push_back(meta);
        }
        const ProductionProperty* property_data = properties.data();
        const int nproperties = static_cast<int>(properties.size());
        const int ncomp = layout.ncomp();
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            auto& low = stage_flux.direction(dir);
            const auto& candidate = high_flux->direction(dir);
            const Real inverse_length = geometry.InvCellSize(dir);
            for (amrex::MFIter mfi(low); mfi.isValid(); ++mfi) {
                const amrex::Box box = mfi.validbox();
                const auto low_arr = low.const_array(mfi);
                const auto high_arr = candidate.const_array(mfi);
                const auto state = output.const_array(mfi);
                const auto result = high_flux->direction(dir).array(mfi);
                const Real scale = dt * anelastic_weight * inverse_length;
                for (int b = 0; b < nbins; ++b) {
                    const int mass = first + b;
                    const int number = two_moment ? population.number_offset + b : -1;
                    const Real lower = two_moment ? population.grid.edges()[static_cast<std::size_t>(b)] : Real(0.0);
                    const Real upper = two_moment ? population.grid.edges()[static_cast<std::size_t>(b+1)] : Real(0.0);
                    const Real denominator = two_moment ? upper - lower : Real(1.0);
                    const Real pivot = two_moment ? Real(1.0) : population.grid.pivot(b);
                    ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        const int li = dir == 0 ? i-1 : i;
                        const int lj = dir == 1 ? j-1 : j;
                        const int lk = dir == 2 ? k-1 : k;
                        const int ri = i, rj = j, rk = k;
                        Real lambda = Real(1.0);
                        const Real delta_mass = high_arr(i,j,k,mass) - low_arr(i,j,k,mass);
                        const Real delta_number = two_moment ?
                            high_arr(i,j,k,number) - low_arr(i,j,k,number) : Real(0.0);
                        if (two_moment) {
                            const Real delta_low = (upper*delta_number - delta_mass) / denominator;
                            const Real delta_high = (delta_mass - lower*delta_number) / denominator;
                            const Real left_low = (upper*state(li,lj,lk,number)-state(li,lj,lk,mass))/denominator;
                            const Real left_high = (state(li,lj,lk,mass)-lower*state(li,lj,lk,number))/denominator;
                            const Real right_low = (upper*state(ri,rj,rk,number)-state(ri,rj,rk,mass))/denominator;
                            const Real right_high = (state(ri,rj,rk,mass)-lower*state(ri,rj,rk,number))/denominator;
                            restrict_constraint(delta_low, left_low, right_low, scale, lambda);
                            restrict_constraint(delta_high, left_high, right_high, scale, lambda);
                        } else {
                            restrict_constraint(delta_mass, state(li,lj,lk,mass), state(ri,rj,rk,mass), scale, lambda);
                        }
                        for (int p = 0; p < nproperties; ++p) {
                            const auto property = property_data[p];
                            const int component = property.component + b;
                            const Real delta_property = high_arr(i,j,k,component) - low_arr(i,j,k,component);
                            restrict_constraint(delta_property, state(li,lj,lk,component),
                                                state(ri,rj,rk,component), scale, lambda);
                            if (property.has_upper) {
                                if (two_moment) {
                                    const Real delta_upper = property.support_max*delta_number - delta_property;
                                    const Real left_upper = property.support_max*state(li,lj,lk,number) - state(li,lj,lk,component);
                                    const Real right_upper = property.support_max*state(ri,rj,rk,number) - state(ri,rj,rk,component);
                                    restrict_constraint(delta_upper, left_upper, right_upper, scale, lambda);
                                    if (property.support_min > Real(0.0)) {
                                        const Real delta_lower = delta_property - property.support_min*delta_number;
                                        const Real left_lower = state(li,lj,lk,component) - property.support_min*state(li,lj,lk,number);
                                        const Real right_lower = state(ri,rj,rk,component) - property.support_min*state(ri,rj,rk,number);
                                        restrict_constraint(delta_lower, left_lower, right_lower, scale, lambda);
                                    }
                                } else {
                                    const Real delta_upper = property.support_max/pivot*delta_mass - delta_property;
                                    const Real left_upper = property.support_max/pivot*state(li,lj,lk,mass) - state(li,lj,lk,component);
                                    const Real right_upper = property.support_max/pivot*state(ri,rj,rk,mass) - state(ri,rj,rk,component);
                                    restrict_constraint(delta_upper, left_upper, right_upper, scale, lambda);
                                    if (property.support_min > Real(0.0)) {
                                        const Real delta_lower = delta_property - property.support_min/pivot*delta_mass;
                                        const Real left_lower = state(li,lj,lk,component) - property.support_min/pivot*state(li,lj,lk,mass);
                                        const Real right_lower = state(ri,rj,rk,component) - property.support_min/pivot*state(ri,rj,rk,mass);
                                        restrict_constraint(delta_lower, left_lower, right_lower, scale, lambda);
                                    }
                                }
                            }
                            if (property.kind == static_cast<int>(PropertyKind::MassBoundedSubset)) {
                                const Real delta_bound = delta_mass - delta_property;
                                const Real left_bound = state(li,lj,lk,mass) - state(li,lj,lk,component);
                                const Real right_bound = state(ri,rj,rk,mass) - state(ri,rj,rk,component);
                                restrict_constraint(delta_bound, left_bound, right_bound, scale, lambda);
                            }
                        }
                        lambda = amrex::max(Real(0.0), amrex::min(Real(1.0), lambda));
                        result(i,j,k,mass) = low_arr(i,j,k,mass) + lambda*delta_mass;
                        if (two_moment) {
                            result(i,j,k,number) = low_arr(i,j,k,number) + lambda*delta_number;
                        }
                        for (int p = 0; p < nproperties; ++p) {
                            const int component = property_data[p].component + b;
                            result(i,j,k,component) = low_arr(i,j,k,component) + lambda*(high_arr(i,j,k,component)-low_arr(i,j,k,component));
                        }
                    });
                }
            }
        }
        // The accepted candidate replaces the high-order candidate in the
        // divergence and in the ledger.  Its correction is scaled by the same
        // stage recurrence coefficient as the low-order update.
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            const auto& accepted = high_flux->direction(dir);
            const auto& low = stage_flux.direction(dir);
            const Real inverse_length = geometry.InvCellSize(dir);
            const int ncomp_local = layout.ncomp();
            for (amrex::MFIter mfi(output); mfi.isValid(); ++mfi) {
                const amrex::Box box = mfi.validbox();
                const auto out = output.array(mfi);
                const auto corr = accepted.const_array(mfi);
                const auto original = low.const_array(mfi);
                ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    for (int c = 0; c < ncomp_local; ++c) {
                        const Real delta = corr(i+ (dir==0), j+(dir==1), k+(dir==2), c) - original(i+ (dir==0), j+(dir==1), k+(dir==2), c);
                        const Real delta_lo = corr(i,j,k,c) - original(i,j,k,c);
                        out(i,j,k,c) -= dt * anelastic_weight * inverse_length * (delta - delta_lo);
                    }
                });
            }
        }
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            amrex::MultiFab::Copy(stage_flux.direction(dir), high_flux->direction(dir),
                                  0, 0, ncomp, 0);
        }
        amrex::Gpu::synchronize();
    }

    if (two_moment) {
        for (amrex::MFIter mfi(output); mfi.isValid(); ++mfi) {
            const amrex::Box box = mfi.validbox();
            const auto source = output.const_array(mfi);
            const auto scratch = transport_scratch.array(mfi);
            for (int b = 0; b < nbins; ++b) {
                const int mass = first + b;
                const int number = population.number_offset + b;
                const Real lower = population.grid.edges()[static_cast<std::size_t>(b)];
                const Real upper = population.grid.edges()[static_cast<std::size_t>(b+1)];
                const Real denominator = upper - lower;
                ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    const Real M = source(i,j,k,mass);
                    const Real C = source(i,j,k,number);
                    const Real scale = amrex::Math::abs(M) + upper*amrex::Math::abs(C) + lower*amrex::Math::abs(C);
                    const Real tolerance = Real(128.0) * std::numeric_limits<Real>::epsilon() * scale;
                    Real L = (upper*C - M) / denominator;
                    Real H = (M - lower*C) / denominator;
                    if (L < Real(0.0) && L >= -tolerance) L = Real(0.0);
                    if (H < Real(0.0) && H >= -tolerance) H = Real(0.0);
                    if (C < Real(0.0) || L < Real(0.0) || H < Real(0.0)) { L = -Real(1.0); H = -Real(1.0); }
                    scratch(i,j,k,mass) = L;
                    scratch(i,j,k,number) = H;
                });
            }
        }
        validate_nonnegative_state(transport_scratch, layout.ncomp(), "SBM two-moment endpoint state");
    }

    // The compact host fields are projections of the updated authoritative
    // spectral state.  This is deliberately separate from face-flux
    // construction: qc/qr never get their own numerical transport path.
    const SBMBulkProjection bulk_projection(layout);
    for (amrex::MFIter mfi(output); mfi.isValid(); ++mfi) {
        bulk_projection.apply_to_core(mfi.validbox(), output.const_array(mfi),
                                      core_state.array(mfi));
    }

    validate_nonnegative_state(output, layout.ncomp());
    manager.accept_stage(level);
    manager.record_stage_face_transfer(level, context, stage_flux);
    core_state.FillBoundary(geometry.periodicity());
}

} // namespace erf_sbm
