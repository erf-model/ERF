#include "ERF_SBMAMRTransfer.H"

#include "AuxiliaryState/ERF_AuxiliaryMappedTransport.H"
#include "ERF_SBMAMRTime.H"
#include "ERF_SBMConstraintGroups.H"
#include "ERF_SBMRemapping.H"
#include "ERF_SBMRestart.H"

#include <AMReX_Arena.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

namespace erf_sbm {

void SBMTransferAllocationCensus::begin () noexcept
{
    baseline_bytes = static_cast<std::uint64_t>(amrex::TotalBytesAllocatedInFabs());
    peak_bytes = baseline_bytes;
    support_loss_collective_votes = 0;
    processed_chunk_first_components.clear();
}

void SBMTransferAllocationCensus::sample () noexcept
{
    peak_bytes = std::max(peak_bytes,
        static_cast<std::uint64_t>(amrex::TotalBytesAllocatedInFabs()));
}

namespace {

using amrex::MultiFab;
using amrex::Real;

bool collective_all_true (const bool local_ok,
                          const std::string& local_diagnostic,
                          const char* remote_failure_message,
                          std::string& diagnostic)
{
    int any_bad = local_ok ? 0 : 1;
    amrex::ParallelDescriptor::ReduceIntMax(any_bad);
    if (any_bad == 0) {
        diagnostic.clear();
        return true;
    }
    diagnostic = local_ok || local_diagnostic.empty()
                     ? remote_failure_message
                     : local_diagnostic;
    return false;
}

bool validate_view_structure (const SBMAMRStateView& view,
                              const SBMLayout& layout, const char* label,
                              std::string& diagnostic)
{
    diagnostic.clear();
    if (view.spectrum == nullptr || view.dry_air_density == nullptr ||
        view.mapped_measure == nullptr) {
        diagnostic = std::string(label) + " SBM transfer tuple is incomplete";
        return false;
    }
    if (!std::isfinite(view.spectrum_time) ||
        !std::isfinite(view.density_time) ||
        !std::isfinite(view.measure_time) ||
        !semantic_time_matches_exact(view.spectrum_time, view.density_time) ||
        !semantic_time_matches_exact(view.spectrum_time, view.measure_time)) {
        diagnostic = std::string(label) +
                     " spectrum, dry-air density, and mapped measure times do not match";
        return false;
    }
    if (view.spectrum->nComp() != layout.ncomp() ||
        view.density_component < 0 ||
        view.density_component >= view.dry_air_density->nComp() ||
        view.measure_component < 0 ||
        view.measure_component >= view.mapped_measure->nComp() ||
        !erf_auxiliary::SameCellLayout(*view.spectrum, *view.dry_air_density) ||
        !erf_auxiliary::SameCellLayout(*view.spectrum, *view.mapped_measure)) {
        diagnostic = std::string(label) +
                     " spectrum, density, and measure have incompatible component or cell layouts";
        return false;
    }
    return true;
}

bool validate_view_carriers (const SBMAMRStateView& view, const char* label,
                             std::string& diagnostic)
{
    std::string density_diagnostic;
    std::string measure_diagnostic;
    // Both checks contain AMReX reductions. Run them in the same order on all
    // ranks instead of short-circuiting after a rank-local result.
    const bool density_ok = erf_auxiliary::ValidatePositiveFiniteComponent(
        *view.dry_air_density, view.density_component, density_diagnostic);
    const bool measure_ok = erf_auxiliary::ValidatePositiveFiniteComponent(
        *view.mapped_measure, view.measure_component, measure_diagnostic);
    if (!density_ok) {
        diagnostic = std::string(label) + " dry-air density: " + density_diagnostic;
        return false;
    }
    if (!measure_ok) {
        diagnostic = std::string(label) + " mapped measure: " + measure_diagnostic;
        return false;
    }
    diagnostic.clear();
    return true;
}

bool validate_view_spectrum (const SBMAMRStateView& view,
                             const SBMLayout& layout, const int level,
                             const SBMStateAdmission* admission,
                             const char* label, std::string& diagnostic)
{
    const bool admissible = admission != nullptr
        ? admission->admissible(*view.spectrum, level, diagnostic)
        : authoritative_state_admissible(*view.spectrum, layout, level,
                                         &diagnostic);
    if (!admissible) {
        diagnostic = std::string(label) + " authoritative spectrum: " + diagnostic;
        return false;
    }
    diagnostic.clear();
    return true;
}

bool valid_ratio (const amrex::IntVect& ratio, std::string& diagnostic)
{
    for (int direction = 0; direction < AMREX_SPACEDIM; ++direction) {
        if (ratio[direction] <= 0) {
            diagnostic = "SBM transfer refinement ratio must be positive in every direction";
            return false;
        }
    }
    return true;
}

bool form_mapped_chunk (const MultiFab& spectrum, const MultiFab& measure,
                        const int measure_component,
                        const std::vector<int>& components,
                        const amrex::Gpu::DeviceVector<int>& device_components,
                        MultiFab& mapped, MultiFab& support, MultiFab& invalid)
{
    amrex::ignore_unused(components);
    const int ncomp = static_cast<int>(components.size());
    const int* component_map = device_components.data();
    for (amrex::MFIter mfi(spectrum, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.tilebox();
        const auto state = spectrum.const_array(mfi);
        const auto omega = measure.const_array(mfi);
        const auto output = mapped.array(mfi);
        const auto positive = support.array(mfi);
        const auto bad = invalid.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int cell_bad = bad(i, j, k, 0) > Real(0.0) ? 1 : 0;
            const Real scale = omega(i, j, k, measure_component);
            for (int local = 0; local < ncomp; ++local) {
                const int component = component_map[local];
                Real amount = Real(0.0);
                const auto status = remap_detail::checked_product(
                    scale, state(i, j, k, component), amount);
                output(i, j, k, local) = amount;
                positive(i, j, k, local) = amount > Real(0.0)
                    ? Real(1.0) : Real(0.0);
                if (status != remap_detail::ProductStatus::Ok) cell_bad = 1;
            }
            bad(i, j, k, 0) = static_cast<Real>(cell_bad);
        });
    }
    return true;
}

void write_restricted_chunk (const MultiFab& mapped,
                             const MultiFab& coarse_support,
                             const MultiFab& measure,
                             const int measure_component,
                             const MultiFab& coverage,
                             const std::vector<int>& components,
                             const amrex::Gpu::DeviceVector<int>& device_components,
                             MultiFab& candidate, MultiFab& lost_support,
                             MultiFab& invalid)
{
    const int ncomp = static_cast<int>(components.size());
    const int* component_map = device_components.data();
    for (amrex::MFIter mfi(candidate, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.tilebox();
        const auto amount = mapped.const_array(mfi);
        const auto support = coarse_support.const_array(mfi);
        const auto omega = measure.const_array(mfi);
        const auto mask = coverage.const_array(mfi);
        const auto output = candidate.array(mfi);
        const auto bad = invalid.array(mfi);
        const auto lost = lost_support.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (mask(i, j, k, 0) <= Real(0.0)) return;
            int cell_bad = bad(i, j, k, 0) > Real(0.0) ? 1 : 0;
            int cell_lost = lost(i, j, k, 0) > Real(0.0) ? 1 : 0;
            const Real scale = omega(i, j, k, measure_component);
            for (int local = 0; local < ncomp; ++local) {
                const int component = component_map[local];
                const Real h = amount(i, j, k, local);
                if (support(i, j, k, local) > Real(0.0) && h == Real(0.0)) {
                    cell_lost = 1;
                }
                Real value = Real(0.0);
                const auto status = remap_detail::checked_quotient(h, scale, value);
                output(i, j, k, component) = value;
                if (status != remap_detail::QuotientStatus::Ok) cell_bad = 1;
            }
            bad(i, j, k, 0) = static_cast<Real>(cell_bad);
            lost(i, j, k, 0) = static_cast<Real>(cell_lost);
        });
    }
}

bool form_carrier_relative_chunk (
    const MultiFab& spectrum, const MultiFab& density,
    const int density_component, const std::vector<int>& components,
    const amrex::Gpu::DeviceVector<int>& device_components,
    MultiFab& intensive, MultiFab& invalid, const MultiFab* used_mask)
{
    const int ncomp = static_cast<int>(components.size());
    const int* component_map = device_components.data();
    const bool has_mask = used_mask != nullptr;
    for (amrex::MFIter mfi(spectrum, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.tilebox();
        const auto state = spectrum.const_array(mfi);
        const auto rho = density.const_array(mfi);
        const auto output = intensive.array(mfi);
        const auto bad = invalid.array(mfi);
        const auto used = has_mask ? used_mask->const_array(mfi)
                                   : amrex::Array4<const Real>{};
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int cell_bad = bad(i, j, k, 0) > Real(0.0) ? 1 : 0;
            const Real carrier = rho(i, j, k, density_component);
            for (int local = 0; local < ncomp; ++local) {
                Real value = Real(0.0);
                const auto status = remap_detail::checked_quotient(
                    state(i, j, k, component_map[local]), carrier, value);
                output(i, j, k, local) = value;
                if (status != remap_detail::QuotientStatus::Ok) cell_bad = 1;
            }
            if (has_mask && used(i, j, k, 0) <= Real(0.0)) cell_bad = 0;
            bad(i, j, k, 0) = static_cast<Real>(cell_bad);
        });
    }
    return true;
}

bool reconstruct_carrier_relative_chunk (
    const MultiFab& intensive, const MultiFab& density,
    const int density_component, const MultiFab& measure,
    const int measure_component, const std::vector<int>& components,
    const amrex::Gpu::DeviceVector<int>& device_components,
    MultiFab& candidate, MultiFab& invalid, const MultiFab* used_mask)
{
    const int ncomp = static_cast<int>(components.size());
    const int* component_map = device_components.data();
    const bool has_mask = used_mask != nullptr;
    for (amrex::MFIter mfi(candidate, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.tilebox();
        const auto z = intensive.const_array(mfi);
        const auto rho = density.const_array(mfi);
        const auto omega = measure.const_array(mfi);
        const auto output = candidate.array(mfi);
        const auto bad = invalid.array(mfi);
        const auto used = has_mask ? used_mask->const_array(mfi)
                                   : amrex::Array4<const Real>{};
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (has_mask && used(i, j, k, 0) <= Real(0.0)) return;
            int cell_bad = bad(i, j, k, 0) > Real(0.0) ? 1 : 0;
            const Real carrier = rho(i, j, k, density_component);
            const Real scale = omega(i, j, k, measure_component);
            for (int local = 0; local < ncomp; ++local) {
                Real value = Real(0.0);
                Real amount = Real(0.0);
                const auto value_status = remap_detail::checked_product(
                    carrier, z(i, j, k, local), value);
                const auto amount_status = remap_detail::checked_product(
                    scale, value, amount);
                output(i, j, k, component_map[local]) = value;
                if (value_status != remap_detail::ProductStatus::Ok ||
                    amount_status != remap_detail::ProductStatus::Ok) {
                    cell_bad = 1;
                }
            }
            bad(i, j, k, 0) = static_cast<Real>(cell_bad);
        });
    }
    return true;
}

bool validate_mapped_carrier_inventory (const SBMAMRStateView& coarse,
                                        const SBMAMRStateView& fine,
                                        const amrex::IntVect& ratio,
                                        SBMTransferAllocationCensus* allocation_census,
                                        std::string& diagnostic,
                                        const MultiFab* prolonged_mask = nullptr)
{
    std::uint64_t child_count = 1;
    for (int direction = 0; direction < AMREX_SPACEDIM; ++direction) {
        const auto factor = static_cast<std::uint64_t>(ratio[direction]);
        if (factor == 0 ||
            child_count > std::numeric_limits<std::uint64_t>::max() / factor) {
            return collective_all_true(false,
                "mapped carrier refinement child count overflows",
                "mapped carrier refinement child count is invalid on another MPI rank",
                diagnostic);
        }
        child_count *= factor;
    }
    const Real unit_roundoff = std::numeric_limits<Real>::epsilon() / Real(2.0);
    const Real fine_ops = static_cast<Real>(child_count) + Real(3.0);
    const Real parent_ops = Real(3.0);
    if (!std::isfinite(fine_ops) || !(fine_ops * unit_roundoff < Real(1.0)) ||
        !(parent_ops * unit_roundoff < Real(1.0))) {
        diagnostic = "mapped carrier roundoff envelope is not defined for this refinement ratio and precision";
        return collective_all_true(false, diagnostic,
            "mapped carrier roundoff envelope is invalid on another MPI rank",
            diagnostic);
    }
    const Real gamma_fine = fine_ops * unit_roundoff /
                            (Real(1.0) - fine_ops * unit_roundoff);
    const Real gamma_parent = parent_ops * unit_roundoff /
                              (Real(1.0) - parent_ops * unit_roundoff);

    MultiFab fine_carrier(fine.spectrum->boxArray(),
                          fine.spectrum->DistributionMap(), 1, 0);
    // Only parents of prolonged cells are checked; without a mask, all are.
    MultiFab fine_coverage;
    if (prolonged_mask == nullptr) {
        fine_coverage.define(fine.spectrum->boxArray(),
                             fine.spectrum->DistributionMap(), 1, 0);
        fine_coverage.setVal(Real(1.0));
    }
    const MultiFab& checked_children =
        prolonged_mask != nullptr ? *prolonged_mask : fine_coverage;
    MultiFab coarse_average(coarse.spectrum->boxArray(),
                            coarse.spectrum->DistributionMap(), 1, 0);
    MultiFab coarse_coverage(coarse.spectrum->boxArray(),
                             coarse.spectrum->DistributionMap(), 1, 0);
    MultiFab fine_invalid(fine.spectrum->boxArray(),
                          fine.spectrum->DistributionMap(), 1, 0);
    MultiFab coarse_invalid(coarse.spectrum->boxArray(),
                            coarse.spectrum->DistributionMap(), 1, 0);
    fine_invalid.setVal(Real(0.0));
    coarse_invalid.setVal(Real(0.0));
    if (allocation_census != nullptr) allocation_census->sample();

    for (amrex::MFIter mfi(*fine.spectrum, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.tilebox();
        const auto rho = fine.dry_air_density->const_array(mfi);
        const auto omega = fine.mapped_measure->const_array(mfi);
        const auto checked = checked_children.const_array(mfi);
        const auto output = fine_carrier.array(mfi);
        const auto bad = fine_invalid.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            Real product = Real(0.0);
            const auto status = remap_detail::checked_product(
                omega(i, j, k, fine.measure_component),
                rho(i, j, k, fine.density_component), product);
            output(i, j, k, 0) = product;
            bad(i, j, k, 0) =
                status == remap_detail::ProductStatus::Ok ||
                        checked(i, j, k, 0) <= Real(0.0)
                    ? Real(0.0) : Real(1.0);
        });
    }

    coarse_coverage.setVal(Real(0.0));
    amrex::average_down(checked_children, coarse_coverage, 0, 1, ratio);
    coarse_average.setVal(Real(0.0));
    amrex::average_down(fine_carrier, coarse_average, 0, 1, ratio);

    MultiFab carrier_error(coarse.spectrum->boxArray(),
                           coarse.spectrum->DistributionMap(), 1, 0);
    MultiFab carrier_tolerance(coarse.spectrum->boxArray(),
                               coarse.spectrum->DistributionMap(), 1, 0);
    carrier_error.setVal(Real(0.0));
    carrier_tolerance.setVal(Real(0.0));
    for (amrex::MFIter mfi(*coarse.spectrum, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.tilebox();
        const auto rho = coarse.dry_air_density->const_array(mfi);
        const auto omega = coarse.mapped_measure->const_array(mfi);
        const auto fine_mean = coarse_average.const_array(mfi);
        const auto coverage = coarse_coverage.const_array(mfi);
        const auto bad = coarse_invalid.array(mfi);
        const auto error = carrier_error.array(mfi);
        const auto tolerance = carrier_tolerance.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (coverage(i, j, k, 0) <= Real(0.0)) return;
            Real coarse_product = Real(0.0);
            const auto status = remap_detail::checked_product(
                omega(i, j, k, coarse.measure_component),
                rho(i, j, k, coarse.density_component), coarse_product);
            const Real fine_value = fine_mean(i, j, k, 0);
            const Real gamma_f = Real(1.0) * gamma_fine;
            const Real gamma_c = Real(1.0) * gamma_parent;
            const Real allowed =
                gamma_f * amrex::Math::abs(fine_value) /
                    (Real(1.0) - gamma_f) +
                gamma_c * amrex::Math::abs(coarse_product);
            const Real residual = amrex::Math::abs(fine_value - coarse_product);
            const bool finite = amrex::Math::isfinite(fine_value) &&
                                amrex::Math::isfinite(coarse_product) &&
                                amrex::Math::isfinite(allowed) &&
                                amrex::Math::isfinite(residual);
            error(i, j, k, 0) = finite ? residual : Real(0.0);
            tolerance(i, j, k, 0) = finite ? allowed : Real(0.0);
            bad(i, j, k, 0) =
                status == remap_detail::ProductStatus::Ok && finite &&
                        residual <= allowed
                    ? Real(0.0) : Real(1.0);
        });
    }

    const bool fine_products_ok = fine_invalid.max(0, 0, true) <= Real(0.0);
    const bool parent_identity_ok = coarse_invalid.max(0, 0, true) <= Real(0.0);
    const Real max_residual = carrier_error.max(0, 0, true);
    const Real max_roundoff_bound = carrier_tolerance.max(0, 0, true);
    const bool local_ok = fine_products_ok && parent_identity_ok;
    std::ostringstream local_diagnostic;
    if (!fine_products_ok) {
        local_diagnostic << "mapped dry-carrier product is not representable";
    } else if (!parent_identity_ok) {
        local_diagnostic << "mapped dry-carrier child inventory differs from its parent beyond the operation-count rounding envelope"
                         << ": max_abs_residual="
                         << std::setprecision(std::numeric_limits<Real>::max_digits10)
                         << max_residual << ", max_allowed_roundoff="
                         << max_roundoff_bound << ", refinement_children="
                         << child_count << ", gamma_fine=" << gamma_fine
                         << ", gamma_parent=" << gamma_parent;
    }
    if (!collective_all_true(
            local_ok,
            local_diagnostic.str(),
            "mapped dry-carrier child inventory differs on another MPI rank",
            diagnostic)) {
        return false;
    }
    diagnostic.clear();
    return true;
}

bool candidate_is_separate (const MultiFab& candidate,
                            const SBMAMRStateView& view)
{
    using AddressRange = std::pair<std::uintptr_t, std::uintptr_t>;
    const auto storage_ranges = [] (const MultiFab& field) {
        std::vector<AddressRange> ranges;
        ranges.reserve(static_cast<std::size_t>(field.size()));
        for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi) {
            const auto& fab = field[mfi];
            const auto begin = reinterpret_cast<std::uintptr_t>(fab.dataPtr());
            const auto size = static_cast<std::uintptr_t>(fab.box().numPts()) *
                              static_cast<std::uintptr_t>(fab.nComp()) * sizeof(Real);
            ranges.emplace_back(begin, begin + size);
        }
        return ranges;
    };
    const auto candidate_ranges = storage_ranges(candidate);
    const auto overlaps = [&candidate, &candidate_ranges, &storage_ranges]
                          (const MultiFab& field) {
        if (&candidate == &field) { return true; }
        const auto field_ranges = storage_ranges(field);
        for (const auto& candidate_range : candidate_ranges) {
            for (const auto& field_range : field_ranges) {
                if (candidate_range.first < field_range.second &&
                    field_range.first < candidate_range.second) {
                    return true;
                }
            }
        }
        return false;
    };
    return view.spectrum != nullptr && view.dry_air_density != nullptr &&
           view.mapped_measure != nullptr && !overlaps(*view.spectrum) &&
           !overlaps(*view.dry_air_density) && !overlaps(*view.mapped_measure);
}

} // namespace

bool ApplyMappedRefluxCorrectionChunk (
    MultiFab& private_candidate, const MultiFab& mapped_correction,
    const MultiFab& mapped_measure, const MultiFab& covered_parent_mask,
    const std::vector<int>& global_components, MultiFab& invalid_mask,
    std::string& diagnostic)
{
    diagnostic.clear();
    bool local_layout_ok = !global_components.empty() &&
        mapped_correction.nComp() >= static_cast<int>(global_components.size()) &&
        invalid_mask.nComp() == 1 && mapped_measure.nComp() > 0 &&
        covered_parent_mask.nComp() == 1 &&
        private_candidate.boxArray() == mapped_correction.boxArray() &&
        private_candidate.boxArray() == mapped_measure.boxArray() &&
        private_candidate.boxArray() == covered_parent_mask.boxArray() &&
        private_candidate.boxArray() == invalid_mask.boxArray() &&
        private_candidate.DistributionMap() == mapped_correction.DistributionMap() &&
        private_candidate.DistributionMap() == mapped_measure.DistributionMap() &&
        private_candidate.DistributionMap() == covered_parent_mask.DistributionMap() &&
        private_candidate.DistributionMap() == invalid_mask.DistributionMap();
    for (int direction = 0; direction < AMREX_SPACEDIM; ++direction) {
        local_layout_ok = local_layout_ok &&
            covered_parent_mask.nGrowVect()[direction] >= 1;
    }
    for (const int component : global_components) {
        local_layout_ok = local_layout_ok && component >= 0 &&
            component < private_candidate.nComp();
    }
    if (!collective_all_true(
            local_layout_ok,
            "localized SBM YA correction fields or component map are incompatible",
            "localized SBM YA correction fields are incompatible on another MPI rank",
            diagnostic)) {
        return false;
    }

    amrex::Gpu::DeviceVector<int> device_components(global_components.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, global_components.begin(),
                     global_components.end(), device_components.begin());
    const int* components = device_components.data();
    const int ncomponents = static_cast<int>(global_components.size());
    for (amrex::MFIter mfi(private_candidate, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.tilebox();
        const auto values = private_candidate.array(mfi);
        const auto delta = mapped_correction.const_array(mfi);
        const auto omega = mapped_measure.const_array(mfi);
        const auto is_covered = covered_parent_mask.const_array(mfi);
        const auto bad = invalid_mask.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int cell_bad = bad(i, j, k, 0) > Real(0.0) ? 1 : 0;
            bool has_covered_face_neighbor = false;
            for (int direction = 0; direction < AMREX_SPACEDIM; ++direction) {
                int lo_i = i, lo_j = j, lo_k = k;
                int hi_i = i, hi_j = j, hi_k = k;
                if (direction == 0) { --lo_i; ++hi_i; }
                if (direction == 1) { --lo_j; ++hi_j; }
                if (direction == 2) { --lo_k; ++hi_k; }
                has_covered_face_neighbor = has_covered_face_neighbor ||
                    is_covered(lo_i, lo_j, lo_k, 0) > Real(0.5) ||
                    is_covered(hi_i, hi_j, hi_k, 0) > Real(0.5);
            }
            const Real scale = omega(i, j, k, 0);
            for (int local = 0; local < ncomponents; ++local) {
                const Real h_correction = delta(i, j, k, local);
                if (h_correction == Real(0.0)) continue;
                if (is_covered(i, j, k, 0) > Real(0.5) ||
                    !has_covered_face_neighbor) {
                    cell_bad = 1;
                    continue;
                }
                Real mapped_old = Real(0.0);
                const auto map_status = remap_detail::checked_product(
                    scale, values(i, j, k, components[local]), mapped_old);
                const Real mapped_new = mapped_old + h_correction;
                Real physical_new = Real(0.0);
                const auto quotient_status = remap_detail::checked_quotient(
                    mapped_new, scale, physical_new);
                if (map_status != remap_detail::ProductStatus::Ok ||
                    quotient_status != remap_detail::QuotientStatus::Ok ||
                    !amrex::Math::isfinite(h_correction) ||
                    !amrex::Math::isfinite(mapped_new)) {
                    cell_bad = 1;
                    continue;
                }
                values(i, j, k, components[local]) = physical_new;
            }
            bad(i, j, k, 0) = static_cast<Real>(cell_bad);
        });
    }
    return true;
}

bool MappedRefluxCorrectionMaskIsValid (const MultiFab& invalid_mask,
                                        std::string& diagnostic)
{
    const bool local_ok = invalid_mask.max(0, 0, true) <= Real(0.0);
    return collective_all_true(
        local_ok,
        "YA reflux produced a correction outside the uncovered coarse/fine interface mask or its mapped update was unrepresentable",
        "YA reflux correction support or mapped update was invalid on another MPI rank",
        diagnostic);
}

static bool restrict_mapped_spectrum_impl (
    const SBMLayout& layout, const SBMAMRStateView& fine,
    const SBMAMRStateView& coarse, const amrex::IntVect& ratio,
    const int coarse_level, MultiFab& coarse_candidate,
    const bool candidate_is_private, const int max_groups_per_chunk,
    const bool reverse_chunk_order, const SBMStateAdmission* admission,
    SBMTransferAllocationCensus* allocation_census,
    const std::optional<double> cross_level_dt_min,
    const MultiFab* covered_parent_mask,
    std::string& diagnostic)
{
    diagnostic.clear();
    if (allocation_census != nullptr) allocation_census->begin();
    std::string ratio_diagnostic;
    const bool ratio_ok = valid_ratio(ratio, ratio_diagnostic);
    const bool level_ok = coarse_level >= 0;
    std::string fine_structure_diagnostic;
    std::string coarse_structure_diagnostic;
    const bool fine_structure_ok = validate_view_structure(
        fine, layout, "fine", fine_structure_diagnostic);
    const bool coarse_structure_ok = validate_view_structure(
        coarse, layout, "coarse", coarse_structure_diagnostic);
    const bool same_time = cross_level_dt_min.has_value()
        ? semantic_time_matches(fine.spectrum_time, coarse.spectrum_time,
                                *cross_level_dt_min)
        : semantic_time_matches_exact(fine.spectrum_time,
                                      coarse.spectrum_time);
    // average_down requires complete fine children. Do not query coarsenability
    // until the ratio is positive and the fine spectrum pointer is valid.
    const bool fine_boxes_coarsenable =
        ratio_ok && fine_structure_ok &&
        fine.spectrum->boxArray().coarsenable(ratio);
    const bool local_structure_ok = ratio_ok && level_ok && fine_structure_ok &&
                                    coarse_structure_ok && same_time &&
                                    fine_boxes_coarsenable;
    const std::string structure_diagnostic =
        !ratio_ok ? ratio_diagnostic
        : !level_ok ? "restriction coarse level must be nonnegative"
        : !fine_structure_ok ? fine_structure_diagnostic
        : !coarse_structure_ok ? coarse_structure_diagnostic
        : !same_time ? "restriction requires fine and coarse state tuples at the same semantic time"
        : !fine_boxes_coarsenable
              ? "restriction fine BoxArray is not coarsenable by the refinement ratio; fine boxes must align to complete coarse cells"
                     : std::string{};
    if (!collective_all_true(
            local_structure_ok, structure_diagnostic,
            "SBM restriction input tuple is invalid on another MPI rank",
            diagnostic)) {
        return false;
    }

    const bool candidate_layout_ok =
        coarse_candidate.nComp() == layout.ncomp() &&
        erf_auxiliary::SameCellLayout(coarse_candidate, *coarse.spectrum);
    const bool coverage_layout_ok = covered_parent_mask == nullptr ||
        (covered_parent_mask->nComp() >= 1 &&
         covered_parent_mask->boxArray() == coarse.spectrum->boxArray() &&
         covered_parent_mask->DistributionMap() ==
             coarse.spectrum->DistributionMap());
    const bool separate_from_fine = candidate_is_separate(coarse_candidate, fine);
    const bool separate_from_coarse = candidate_is_separate(coarse_candidate, coarse);
    const bool local_candidate_ok = candidate_layout_ok && coverage_layout_ok &&
                                    separate_from_fine && separate_from_coarse;
    if (!collective_all_true(
            local_candidate_ok,
            "restriction candidate must be separate from both input tuples and match the coarse SBM layout, as must any covered-parent mask",
            "SBM restriction candidate aliases an input on another MPI rank",
            diagnostic)) {
        return false;
    }

    std::string fine_carrier_diagnostic;
    std::string coarse_carrier_diagnostic;
    const bool fine_carriers_ok = validate_view_carriers(
        fine, "fine", fine_carrier_diagnostic);
    const bool coarse_carriers_ok = validate_view_carriers(
        coarse, "coarse", coarse_carrier_diagnostic);
    std::string fine_spectrum_diagnostic;
    std::string coarse_spectrum_diagnostic;
    const bool fine_spectrum_ok = validate_view_spectrum(
        fine, layout, coarse_level + 1, admission, "fine",
        fine_spectrum_diagnostic);
    const bool coarse_spectrum_ok = validate_view_spectrum(
        coarse, layout, coarse_level, admission, "coarse",
        coarse_spectrum_diagnostic);
    const bool local_sources_ok = fine_carriers_ok && coarse_carriers_ok &&
                                  fine_spectrum_ok && coarse_spectrum_ok;
    const std::string source_diagnostic =
        !fine_carriers_ok ? fine_carrier_diagnostic
        : !coarse_carriers_ok ? coarse_carrier_diagnostic
        : !fine_spectrum_ok ? fine_spectrum_diagnostic
        : !coarse_spectrum_ok ? coarse_spectrum_diagnostic
                              : std::string{};
    if (!collective_all_true(
            local_sources_ok, source_diagnostic,
            "SBM restriction source is inadmissible on another MPI rank",
            diagnostic)) {
        return false;
    }

    if (max_groups_per_chunk <= 0) {
        diagnostic = "SBM restriction chunk group limit must be positive";
        return collective_all_true(false, diagnostic,
            "SBM restriction chunk group limit is invalid on another MPI rank",
            diagnostic);
    }
    const auto chunks = make_constraint_closure_chunks(
        layout, max_groups_per_chunk);
    int max_chunk_components = 0;
    for (std::size_t pass = 0; pass < chunks.size(); ++pass) {
        const std::size_t chunk_index = reverse_chunk_order
            ? chunks.size() - 1U - pass : pass;
        const auto& chunk = chunks[chunk_index];
        max_chunk_components = std::max(
            max_chunk_components, static_cast<int>(chunk.components.size()));
    }
    if (!collective_all_true(
            max_chunk_components > 0,
            "SBM restriction layout has no complete semantic component groups",
            "SBM restriction layout has no complete semantic component groups on another MPI rank",
            diagnostic)) {
        return false;
    }

    MultiFab fine_h_chunk(fine.spectrum->boxArray(),
                          fine.spectrum->DistributionMap(),
                          max_chunk_components, 0);
    MultiFab fine_support_chunk(fine.spectrum->boxArray(),
                                fine.spectrum->DistributionMap(),
                                max_chunk_components, 0);
    MultiFab coarse_h_chunk(coarse.spectrum->boxArray(),
                            coarse.spectrum->DistributionMap(),
                            max_chunk_components, 0);
    MultiFab coarse_support_chunk(coarse.spectrum->boxArray(),
                                  coarse.spectrum->DistributionMap(),
                                  max_chunk_components, 0);
    MultiFab fine_coverage;
    MultiFab coarse_coverage_scratch;
    if (covered_parent_mask == nullptr) {
        fine_coverage.define(fine.spectrum->boxArray(),
                             fine.spectrum->DistributionMap(), 1, 0);
        coarse_coverage_scratch.define(coarse.spectrum->boxArray(),
                                       coarse.spectrum->DistributionMap(), 1, 0);
    }
    MultiFab lost_positive_support(coarse.spectrum->boxArray(),
                                   coarse.spectrum->DistributionMap(), 1, 0);
    MultiFab fine_invalid(fine.spectrum->boxArray(),
                          fine.spectrum->DistributionMap(), 1, 0);
    MultiFab coarse_invalid(coarse.spectrum->boxArray(),
                            coarse.spectrum->DistributionMap(), 1, 0);
    MultiFab owned_candidate_scratch;
    MultiFab* work = &coarse_candidate;
    if (!candidate_is_private) {
        owned_candidate_scratch.define(coarse.spectrum->boxArray(),
                                       coarse.spectrum->DistributionMap(),
                                       layout.ncomp(), 0);
        work = &owned_candidate_scratch;
    }
    if (allocation_census != nullptr) allocation_census->sample();
    // One full private coarse candidate starts as accepted coarse U. Only
    // covered parents are replaced by H-restricted values; uncovered U never
    // takes a mapped round trip.
    MultiFab::Copy(*work, *coarse.spectrum, 0, 0, layout.ncomp(), 0);
    fine_invalid.setVal(Real(0.0));
    coarse_invalid.setVal(Real(0.0));
    lost_positive_support.setVal(Real(0.0));
    if (covered_parent_mask == nullptr) {
        fine_coverage.setVal(Real(1.0));
        coarse_coverage_scratch.setVal(Real(0.0));
        amrex::average_down(fine_coverage, coarse_coverage_scratch, 0, 1, ratio);
    }
    const MultiFab& coarse_coverage = covered_parent_mask != nullptr
        ? *covered_parent_mask : coarse_coverage_scratch;
    amrex::Gpu::DeviceVector<int> device_components(
        static_cast<std::size_t>(max_chunk_components));
    for (std::size_t pass = 0; pass < chunks.size(); ++pass) {
        const std::size_t chunk_index = reverse_chunk_order
            ? chunks.size() - 1U - pass : pass;
        const auto& chunk = chunks[chunk_index];
        if (allocation_census != nullptr && !chunk.components.empty()) {
            allocation_census->processed_chunk_first_components.push_back(
                chunk.components.front());
        }
        const int chunk_components = static_cast<int>(chunk.components.size());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, chunk.components.begin(),
                         chunk.components.end(), device_components.begin());
        form_mapped_chunk(*fine.spectrum, *fine.mapped_measure,
                          fine.measure_component, chunk.components,
                          device_components, fine_h_chunk,
                          fine_support_chunk, fine_invalid);
        amrex::average_down(fine_h_chunk, coarse_h_chunk, 0,
                            chunk_components, ratio);
        amrex::average_down(fine_support_chunk, coarse_support_chunk, 0,
                            chunk_components, ratio);
        write_restricted_chunk(coarse_h_chunk, coarse_support_chunk,
                               *coarse.mapped_measure,
                               coarse.measure_component, coarse_coverage,
                               chunk.components, device_components, *work,
                               lost_positive_support, coarse_invalid);
    }

    // Support loss is tracked for every component and decided collectively
    // once. average_down retains its required AMReX data communication.
    const bool fine_mapping_valid = fine_invalid.max(0, 0, true) <= Real(0.0);
    const bool coarse_quotient_valid =
        coarse_invalid.max(0, 0, true) <= Real(0.0);
    const bool support_preserved =
        lost_positive_support.max(0, 0, true) <= Real(0.0);
    const bool local_support_preserved = fine_mapping_valid &&
                                         coarse_quotient_valid &&
                                         support_preserved;
    const std::string local_transfer_diagnostic =
        !fine_mapping_valid
            ? "forming fine mapped amount H=omega*U overflowed or underflowed"
        : !coarse_quotient_valid
            ? "recovering coarse spectrum U=H/omega overflowed or underflowed"
        : !support_preserved
            ? "restriction averaging underflow erased positive mapped support to zero"
            : std::string{};
    if (allocation_census != nullptr) {
        ++allocation_census->support_loss_collective_votes;
    }
    if (!collective_all_true(local_support_preserved,
                             local_transfer_diagnostic,
                             "mapped restriction overflowed, underflowed, or erased positive support on another MPI rank",
                             diagnostic)) {
        return false;
    }

    std::string candidate_admission_diagnostic;
    const bool candidate_admissible = admission != nullptr
        ? admission->admissible(*work, coarse_level,
                                candidate_admission_diagnostic)
        : authoritative_state_admissible(*work, layout, coarse_level,
                                         &candidate_admission_diagnostic);
    const bool local_candidate_admissible = candidate_admissible;
    const std::string candidate_diagnostic = !candidate_admissible
        ? "restricted SBM candidate is inadmissible: " +
              candidate_admission_diagnostic
        : std::string{};
    if (!collective_all_true(
            local_candidate_admissible, candidate_diagnostic,
            "restricted SBM candidate is inadmissible on another MPI rank",
            diagnostic)) {
        return false;
    }
    if (!candidate_is_private) {
        MultiFab::Copy(coarse_candidate, *work, 0, 0, layout.ncomp(), 0);
    }
    return true;
}

bool RestrictMappedSpectrum (const SBMLayout& layout,
                             const SBMAMRStateView& fine,
                             const SBMAMRStateView& coarse,
                             const amrex::IntVect& ratio,
                             const int coarse_level,
                             MultiFab& coarse_candidate,
                             std::string& diagnostic,
                             const int max_groups_per_chunk,
                             const bool reverse_chunk_order,
                             const SBMStateAdmission* admission,
                             SBMTransferAllocationCensus* allocation_census,
                             const std::optional<double> cross_level_dt_min)
{
    return restrict_mapped_spectrum_impl(layout, fine, coarse, ratio,
                                         coarse_level, coarse_candidate, false,
                                         max_groups_per_chunk,
                                         reverse_chunk_order, admission,
                                         allocation_census,
                                         cross_level_dt_min, nullptr,
                                         diagnostic);
}

bool RestrictMappedSpectrumPrivate (const SBMLayout& layout,
                                    const SBMAMRStateView& fine,
                                    const SBMAMRStateView& coarse,
                                    const amrex::IntVect& ratio,
                                    const int coarse_level,
                                    MultiFab& private_candidate,
                                    std::string& diagnostic,
                                    const int max_groups_per_chunk,
                                    const bool reverse_chunk_order,
                                    const SBMStateAdmission* admission,
                                    SBMTransferAllocationCensus* allocation_census,
                                    const std::optional<double> cross_level_dt_min,
                                    const MultiFab* covered_parent_mask)
{
    return restrict_mapped_spectrum_impl(layout, fine, coarse, ratio,
                                         coarse_level, private_candidate, true,
                                         max_groups_per_chunk,
                                         reverse_chunk_order, admission,
                                         allocation_census,
                                         cross_level_dt_min,
                                         covered_parent_mask, diagnostic);
}

// Write the same-time piecewise-constant carrier-relative prolongation of the
// coarse spectrum into every valid cell of work. The fine tuple supplies only
// the target density and measure; work defines the target layout.
static bool prolong_carrier_relative_chunks (
    const SBMLayout& layout, const SBMAMRStateView& coarse,
    const SBMAMRStateView& fine_target,
    const amrex::Geometry& coarse_geometry,
    const amrex::Geometry& fine_geometry,
    const amrex::IntVect& ratio, MultiFab& work,
    const int max_groups_per_chunk, const bool reverse_chunk_order,
    SBMTransferAllocationCensus* allocation_census,
    std::string& diagnostic,
    const MultiFab* prolonged_mask = nullptr)
{
    if (max_groups_per_chunk <= 0) {
        diagnostic = "SBM prolongation chunk group limit must be positive";
        return collective_all_true(false, diagnostic,
            "SBM prolongation chunk group limit is invalid on another MPI rank",
            diagnostic);
    }
    const auto chunks = make_constraint_closure_chunks(
        layout, max_groups_per_chunk);
    int max_chunk_components = 0;
    for (const auto& chunk : chunks) {
        max_chunk_components = std::max(
            max_chunk_components, static_cast<int>(chunk.components.size()));
    }
    if (!collective_all_true(
            max_chunk_components > 0,
            "SBM prolongation layout has no complete semantic component groups",
            "SBM prolongation layout has no complete semantic component groups on another MPI rank",
            diagnostic)) {
        return false;
    }

    MultiFab coarse_z_chunk(coarse.spectrum->boxArray(),
                            coarse.spectrum->DistributionMap(),
                            max_chunk_components, 0);
    MultiFab fine_z_chunk(work.boxArray(), work.DistributionMap(),
                          max_chunk_components, 0);
    MultiFab coarse_invalid(coarse.spectrum->boxArray(),
                            coarse.spectrum->DistributionMap(), 1, 0);
    MultiFab fine_invalid(work.boxArray(), work.DistributionMap(), 1, 0);
    coarse_invalid.setVal(Real(0.0));
    fine_invalid.setVal(Real(0.0));
    // With a mask, only the parents of prolonged cells must form a ratio.
    MultiFab parent_mask;
    if (prolonged_mask != nullptr) {
        parent_mask.define(coarse.spectrum->boxArray(),
                           coarse.spectrum->DistributionMap(), 1, 0);
        parent_mask.setVal(Real(0.0));
        amrex::average_down(*prolonged_mask, parent_mask, 0, 1, ratio);
    }
    if (allocation_census != nullptr) allocation_census->sample();

    amrex::Gpu::DeviceVector<int> device_components(
        static_cast<std::size_t>(max_chunk_components));
    amrex::PhysBCFunctNoOp no_physical_boundary;
    for (std::size_t pass = 0; pass < chunks.size(); ++pass) {
        const std::size_t chunk_index = reverse_chunk_order
            ? chunks.size() - 1U - pass : pass;
        const auto& chunk = chunks[chunk_index];
        if (allocation_census != nullptr && !chunk.components.empty()) {
            allocation_census->processed_chunk_first_components.push_back(
                chunk.components.front());
        }
        const int chunk_components = static_cast<int>(chunk.components.size());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, chunk.components.begin(),
                         chunk.components.end(), device_components.begin());
        form_carrier_relative_chunk(
            *coarse.spectrum, *coarse.dry_air_density,
            coarse.density_component, chunk.components, device_components,
            coarse_z_chunk, coarse_invalid,
            prolonged_mask != nullptr ? &parent_mask : nullptr);

        amrex::Vector<amrex::BCRec> bcs(
            static_cast<std::size_t>(chunk_components));
        // Valid fine cells only: piecewise-constant prolongation reads the
        // in-domain parent and never these records or a physical ghost.
        for (auto& bc : bcs) {
            for (int direction = 0; direction < AMREX_SPACEDIM; ++direction) {
                const int type = fine_geometry.isPeriodic(direction)
                    ? amrex::BCType::int_dir : amrex::BCType::foextrap;
                bc.setLo(direction, type);
                bc.setHi(direction, type);
            }
        }
        amrex::InterpFromCoarseLevel(
            fine_z_chunk, amrex::IntVect(0),
            static_cast<Real>(coarse.spectrum_time), coarse_z_chunk, 0, 0,
            chunk_components, coarse_geometry, fine_geometry,
            no_physical_boundary, 0, no_physical_boundary, 0, ratio,
            &amrex::pc_interp, bcs, 0);
        reconstruct_carrier_relative_chunk(
            fine_z_chunk, *fine_target.dry_air_density,
            fine_target.density_component, *fine_target.mapped_measure,
            fine_target.measure_component, chunk.components, device_components,
            work, fine_invalid, prolonged_mask);
    }

    const bool reconstructed_ok =
        coarse_invalid.max(0, 0, true) <= Real(0.0) &&
        fine_invalid.max(0, 0, true) <= Real(0.0);
    return collective_all_true(
        reconstructed_ok,
        "forming or reconstructing a chunked carrier-relative SBM state overflowed or underflowed",
        "forming or reconstructing a chunked carrier-relative SBM state failed on another MPI rank",
        diagnostic);
}

static bool prolong_carrier_relative_spectrum_impl (
    const SBMLayout& layout, const SBMAMRStateView& coarse,
    const SBMAMRStateView& fine_target,
    const amrex::Geometry& coarse_geometry,
    const amrex::Geometry& fine_geometry,
    const amrex::IntVect& ratio, const int fine_level,
    MultiFab& fine_candidate, const bool candidate_is_private,
    const int max_groups_per_chunk, const bool reverse_chunk_order,
    const SBMStateAdmission* admission,
    SBMTransferAllocationCensus* allocation_census,
    std::string& diagnostic)
{
    diagnostic.clear();
    if (allocation_census != nullptr) allocation_census->begin();
    std::string ratio_diagnostic;
    const bool ratio_ok = valid_ratio(ratio, ratio_diagnostic);
    const bool level_ok = fine_level > 0;

    std::string coarse_structure_diagnostic;
    std::string fine_structure_diagnostic;
    const bool coarse_structure_ok = validate_view_structure(
        coarse, layout, "coarse", coarse_structure_diagnostic);
    const bool fine_structure_ok = validate_view_structure(
        fine_target, layout, "fine target", fine_structure_diagnostic);
    const bool same_time = semantic_time_matches_exact(coarse.spectrum_time,
                                                        fine_target.spectrum_time);
    const bool local_structure_ok = ratio_ok && level_ok && coarse_structure_ok &&
                                    fine_structure_ok && same_time;
    const std::string structure_diagnostic =
        !ratio_ok ? ratio_diagnostic
        : !level_ok ? "prolongation fine level must be greater than zero"
        : !coarse_structure_ok ? coarse_structure_diagnostic
        : !fine_structure_ok ? fine_structure_diagnostic
        : !same_time ? "prolongation requires same-time source and target tuples"
                     : std::string{};
    if (!collective_all_true(
            local_structure_ok, structure_diagnostic,
            "SBM prolongation input tuple is invalid on another MPI rank",
            diagnostic)) {
        return false;
    }

    const bool candidate_layout_ok =
        fine_candidate.nComp() == layout.ncomp() &&
        erf_auxiliary::SameCellLayout(fine_candidate, *fine_target.spectrum);
    const bool separate_from_coarse = candidate_is_separate(fine_candidate, coarse);
    const bool separate_from_fine_target =
        candidate_is_separate(fine_candidate, fine_target);
    bool local_candidate_ok = candidate_layout_ok && separate_from_coarse &&
                              separate_from_fine_target;
    std::string local_candidate_diagnostic =
        "prolongation candidate must be separate from both input tuples and match the fine SBM layout";
    if (local_candidate_ok) {
        for (int direction = 0; direction < AMREX_SPACEDIM; ++direction) {
            if (coarse_geometry.isPeriodic(direction) !=
                fine_geometry.isPeriodic(direction)) {
                local_candidate_ok = false;
                local_candidate_diagnostic =
                    "carrier-relative prolongation requires matching coarse and "
                    "fine periodicity";
                break;
            }
        }
    }
    if (local_candidate_ok) {
        amrex::Box refined_coarse_domain = coarse_geometry.Domain();
        refined_coarse_domain.refine(ratio);
        if (refined_coarse_domain != fine_geometry.Domain()) {
            local_candidate_ok = false;
            local_candidate_diagnostic =
                "prolongation refinement ratio does not map the coarse domain "
                "to the fine domain";
        }
    }
    if (!collective_all_true(
            local_candidate_ok, local_candidate_diagnostic,
            "SBM prolongation candidate or geometry is invalid on another MPI rank",
            diagnostic)) {
        return false;
    }

    std::string coarse_carrier_diagnostic;
    std::string fine_carrier_diagnostic;
    const bool coarse_carriers_ok = validate_view_carriers(
        coarse, "coarse", coarse_carrier_diagnostic);
    const bool fine_carriers_ok = validate_view_carriers(
        fine_target, "fine target", fine_carrier_diagnostic);
    std::string coarse_spectrum_diagnostic;
    const bool coarse_spectrum_ok = validate_view_spectrum(
        coarse, layout, fine_level - 1, admission, "coarse",
        coarse_spectrum_diagnostic);
    const bool local_sources_ok = coarse_carriers_ok && fine_carriers_ok &&
                                  coarse_spectrum_ok;
    const std::string source_diagnostic =
        !coarse_carriers_ok ? coarse_carrier_diagnostic
        : !fine_carriers_ok ? fine_carrier_diagnostic
        : !coarse_spectrum_ok ? coarse_spectrum_diagnostic
                              : std::string{};
    if (!collective_all_true(
            local_sources_ok, source_diagnostic,
            "SBM prolongation source is inadmissible on another MPI rank",
            diagnostic)) {
        return false;
    }

    if (!validate_mapped_carrier_inventory(coarse, fine_target, ratio,
                                           allocation_census,
                                           diagnostic)) {
        return false;
    }

    MultiFab private_candidate;
    MultiFab* work = &fine_candidate;
    if (!candidate_is_private) {
        private_candidate.define(fine_target.spectrum->boxArray(),
                                 fine_target.spectrum->DistributionMap(),
                                 layout.ncomp(), 0);
        work = &private_candidate;
    }
    if (!prolong_carrier_relative_chunks(
            layout, coarse, fine_target, coarse_geometry, fine_geometry,
            ratio, *work, max_groups_per_chunk, reverse_chunk_order,
            allocation_census, diagnostic)) {
        return false;
    }

    std::string candidate_admission_diagnostic;
    const bool candidate_admissible = admission != nullptr
        ? admission->admissible(*work, fine_level,
                                candidate_admission_diagnostic)
        : authoritative_state_admissible(*work, layout, fine_level,
                                         &candidate_admission_diagnostic);
    const bool local_candidate_admissible = candidate_admissible;
    const std::string candidate_diagnostic = !candidate_admissible
        ? "prolonged SBM candidate is inadmissible: " +
              candidate_admission_diagnostic
        : std::string{};
    if (!collective_all_true(
            local_candidate_admissible, candidate_diagnostic,
            "prolonged SBM candidate is inadmissible on another MPI rank",
            diagnostic)) {
        return false;
    }
    if (!candidate_is_private) {
        MultiFab::Copy(fine_candidate, *work, 0, 0, layout.ncomp(), 0);
    }
    return true;
}

bool ProlongCarrierRelativeSpectrum (const SBMLayout& layout,
                                     const SBMAMRStateView& coarse,
                                     const SBMAMRStateView& fine_target,
                                     const amrex::Geometry& coarse_geometry,
                                     const amrex::Geometry& fine_geometry,
                                     const amrex::IntVect& ratio,
                                     const int fine_level,
                                     MultiFab& fine_candidate,
                                     std::string& diagnostic,
                                     const int max_groups_per_chunk,
                                     const bool reverse_chunk_order,
                                     const SBMStateAdmission* admission,
                                     SBMTransferAllocationCensus* allocation_census)
{
    return prolong_carrier_relative_spectrum_impl(
        layout, coarse, fine_target, coarse_geometry, fine_geometry, ratio,
        fine_level, fine_candidate, false, max_groups_per_chunk,
        reverse_chunk_order, admission, allocation_census, diagnostic);
}

bool ProlongCarrierRelativeSpectrumPrivate (
    const SBMLayout& layout, const SBMAMRStateView& coarse,
    const SBMAMRStateView& fine_target,
    const amrex::Geometry& coarse_geometry,
    const amrex::Geometry& fine_geometry,
    const amrex::IntVect& ratio, const int fine_level,
    MultiFab& private_candidate, std::string& diagnostic,
    const int max_groups_per_chunk, const bool reverse_chunk_order,
    const SBMStateAdmission* admission,
    SBMTransferAllocationCensus* allocation_census)
{
    return prolong_carrier_relative_spectrum_impl(
        layout, coarse, fine_target, coarse_geometry, fine_geometry, ratio,
        fine_level, private_candidate, true, max_groups_per_chunk,
        reverse_chunk_order, admission, allocation_census, diagnostic);
}

SBMRegridCellCounts CountRegridCells (const amrex::BoxArray& old_grids,
                                      const amrex::BoxArray& new_grids)
{
    SBMRegridCellCounts counts;
    counts.cells = new_grids.numPts();
    for (int box = 0; !old_grids.empty() && box < new_grids.size(); ++box) {
        for (const auto& overlap : old_grids.intersections(new_grids[box])) {
            counts.retained += overlap.second.numPts();
        }
    }
    counts.prolonged = counts.cells - counts.retained;
    counts.released = old_grids.numPts() - counts.retained;
    return counts;
}

bool RemapRegriddedSpectrumPrivate (
    const SBMLayout& layout, const SBMAMRStateView& coarse,
    const MultiFab& retained_spectrum,
    const MultiFab& fine_density, const int density_component,
    const MultiFab& fine_measure, const int measure_component,
    const double semantic_time,
    const amrex::Geometry& coarse_geometry,
    const amrex::Geometry& fine_geometry,
    const amrex::IntVect& ratio, const int fine_level,
    MultiFab& private_candidate, std::string& diagnostic,
    const int max_groups_per_chunk,
    const SBMStateAdmission* coarse_admission,
    const SBMStateAdmission* fine_admission,
    SBMTransferAllocationCensus* allocation_census)
{
    diagnostic.clear();
    if (allocation_census != nullptr) allocation_census->begin();
    // The candidate supplies the new-grid layout; it is never read.
    const SBMAMRStateView fine_target{
        &private_candidate, semantic_time, &fine_density, density_component,
        semantic_time, &fine_measure, measure_component, semantic_time};
    const SBMAMRStateView retained{
        &retained_spectrum, semantic_time, &fine_density, density_component,
        semantic_time, &fine_measure, measure_component, semantic_time};

    std::string ratio_diagnostic;
    const bool ratio_ok = valid_ratio(ratio, ratio_diagnostic);
    const bool level_ok = fine_level > 0;
    std::string coarse_structure_diagnostic;
    std::string fine_structure_diagnostic;
    const bool coarse_structure_ok = validate_view_structure(
        coarse, layout, "coarse", coarse_structure_diagnostic);
    const bool fine_structure_ok = validate_view_structure(
        fine_target, layout, "regrid target", fine_structure_diagnostic);
    const bool same_time = semantic_time_matches_exact(coarse.spectrum_time,
                                                        semantic_time);
    const bool retained_ok = retained_spectrum.ok() &&
                             retained_spectrum.nComp() == layout.ncomp();
    const bool local_structure_ok = ratio_ok && level_ok && coarse_structure_ok &&
                                    fine_structure_ok && same_time && retained_ok;
    const std::string structure_diagnostic =
        !ratio_ok ? ratio_diagnostic
        : !level_ok ? "regrid fine level must be greater than zero"
        : !coarse_structure_ok ? coarse_structure_diagnostic
        : !fine_structure_ok ? fine_structure_diagnostic
        : !same_time ? "regrid requires same-time coarse and fine tuples"
        : !retained_ok ? "retained regrid spectrum does not match the SBM layout"
                       : std::string{};
    if (!collective_all_true(
            local_structure_ok, structure_diagnostic,
            "SBM regrid input tuple is invalid on another MPI rank",
            diagnostic)) {
        return false;
    }

    bool local_candidate_ok = candidate_is_separate(private_candidate, coarse) &&
                              candidate_is_separate(private_candidate, retained);
    std::string local_candidate_diagnostic =
        "regrid candidate must be separate from the coarse, retained, and carrier inputs";
    if (local_candidate_ok) {
        for (int direction = 0; direction < AMREX_SPACEDIM; ++direction) {
            if (coarse_geometry.isPeriodic(direction) !=
                fine_geometry.isPeriodic(direction)) {
                local_candidate_ok = false;
                local_candidate_diagnostic =
                    "carrier-relative regrid requires matching coarse and fine "
                    "periodicity";
                break;
            }
        }
    }
    if (local_candidate_ok) {
        amrex::Box refined_coarse_domain = coarse_geometry.Domain();
        refined_coarse_domain.refine(ratio);
        if (refined_coarse_domain != fine_geometry.Domain() ||
            (!retained_spectrum.boxArray().empty() &&
             !fine_geometry.Domain().contains(
                 retained_spectrum.boxArray().minimalBox()))) {
            local_candidate_ok = false;
            local_candidate_diagnostic =
                "regrid refinement ratio or retained grids do not match the fine domain";
        }
    }
    if (!collective_all_true(
            local_candidate_ok, local_candidate_diagnostic,
            "SBM regrid candidate or geometry is invalid on another MPI rank",
            diagnostic)) {
        return false;
    }

    std::string coarse_carrier_diagnostic;
    std::string fine_carrier_diagnostic;
    const bool coarse_carriers_ok = validate_view_carriers(
        coarse, "coarse", coarse_carrier_diagnostic);
    const bool fine_carriers_ok = validate_view_carriers(
        fine_target, "regrid target", fine_carrier_diagnostic);
    std::string coarse_spectrum_diagnostic;
    const bool coarse_spectrum_ok = validate_view_spectrum(
        coarse, layout, fine_level - 1, coarse_admission, "coarse",
        coarse_spectrum_diagnostic);
    const bool local_sources_ok = coarse_carriers_ok && fine_carriers_ok &&
                                  coarse_spectrum_ok;
    const std::string source_diagnostic =
        !coarse_carriers_ok ? coarse_carrier_diagnostic
        : !fine_carriers_ok ? fine_carrier_diagnostic
        : !coarse_spectrum_ok ? coarse_spectrum_diagnostic
                              : std::string{};
    if (!collective_all_true(
            local_sources_ok, source_diagnostic,
            "SBM regrid source is inadmissible on another MPI rank",
            diagnostic)) {
        return false;
    }

    // Retained cells are copied, so only newly covered cells and their parents
    // are prolonged and checked.
    MultiFab prolonged_mask(private_candidate.boxArray(),
                            private_candidate.DistributionMap(), 1, 0);
    prolonged_mask.setVal(Real(1.0));
    const amrex::BoxArray& retained_grids = retained_spectrum.boxArray();
    if (!retained_grids.empty()) {
        for (amrex::MFIter mfi(prolonged_mask); mfi.isValid(); ++mfi) {
            const auto mask = prolonged_mask.array(mfi);
            for (const auto& overlap : retained_grids.intersections(mfi.validbox())) {
                amrex::ParallelFor(overlap.second,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        mask(i, j, k, 0) = Real(0.0);
                    });
            }
        }
    }
    if (!validate_mapped_carrier_inventory(coarse, fine_target, ratio,
                                           allocation_census, diagnostic,
                                           &prolonged_mask)) {
        return false;
    }
    if (!prolong_carrier_relative_chunks(
            layout, coarse, fine_target, coarse_geometry, fine_geometry,
            ratio, private_candidate, max_groups_per_chunk, false,
            allocation_census, diagnostic, &prolonged_mask)) {
        return false;
    }
    private_candidate.ParallelCopy(retained_spectrum, 0, 0, layout.ncomp());

    std::string candidate_admission_diagnostic;
    const bool candidate_admissible = fine_admission != nullptr
        ? fine_admission->admissible(private_candidate, fine_level,
                                     candidate_admission_diagnostic)
        : authoritative_state_admissible(private_candidate, layout, fine_level,
                                         &candidate_admission_diagnostic);
    return collective_all_true(
        candidate_admissible,
        candidate_admissible ? std::string{}
            : "regridded SBM candidate is inadmissible: " +
                  candidate_admission_diagnostic,
        "regridded SBM candidate is inadmissible on another MPI rank",
        diagnostic);
}

} // namespace erf_sbm
