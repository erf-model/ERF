#include "ERF_SBMAMRTransfer.H"

#include "AuxiliaryState/ERF_AuxiliaryMappedTransport.H"
#include "ERF_SBMRemapping.H"
#include "ERF_SBMRestart.H"

#include <AMReX_FillPatchUtil.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

namespace erf_sbm {
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
        view.spectrum_time != view.density_time ||
        view.spectrum_time != view.measure_time) {
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
                             const char* label, std::string& diagnostic)
{
    if (!authoritative_state_admissible(*view.spectrum, layout, level,
                                        &diagnostic)) {
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

bool form_mapped_state (const MultiFab& spectrum, const MultiFab& measure,
                        const int measure_component, MultiFab& mapped,
                        std::string& diagnostic)
{
    MultiFab invalid(spectrum.boxArray(), spectrum.DistributionMap(), 1, 0);
    invalid.setVal(Real(0.0));
    const int ncomp = spectrum.nComp();
    for (amrex::MFIter mfi(spectrum, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.tilebox();
        const auto state = spectrum.const_array(mfi);
        const auto omega = measure.const_array(mfi);
        const auto output = mapped.array(mfi);
        const auto bad = invalid.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int cell_bad = 0;
            const Real scale = omega(i, j, k, measure_component);
            for (int component = 0; component < ncomp; ++component) {
                Real amount = Real(0.0);
                const auto status = remap_detail::checked_product(
                    scale, state(i, j, k, component), amount);
                output(i, j, k, component) = amount;
                if (status != remap_detail::ProductStatus::Ok) cell_bad = 1;
            }
            bad(i, j, k, 0) = static_cast<Real>(cell_bad);
        });
    }
    const bool locally_valid = invalid.max(0, 0, true) <= Real(0.0);
    if (!locally_valid) {
        diagnostic = "forming mapped spectral amount H=omega*U overflowed or underflowed";
        return false;
    }
    return true;
}

bool
divide_mapped_state (const MultiFab& mapped,
                     const MultiFab& measure,
                     const int measure_component,
                     const MultiFab& covered_coarse,
                     const int coverage_component,
                     MultiFab& candidate,
                     std::string& diagnostic)
{
    MultiFab invalid(mapped.boxArray(), mapped.DistributionMap(), 1, 0);
    invalid.setVal(Real(0.0));
    const int ncomp = mapped.nComp();
    for (amrex::MFIter mfi(mapped, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.tilebox();
        const auto amount = mapped.const_array(mfi);
        const auto omega = measure.const_array(mfi);
        const auto coverage = covered_coarse.const_array(mfi);
        const auto output = candidate.array(mfi);
        const auto bad = invalid.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int cell_bad = 0;
            if (coverage(i, j, k, coverage_component) > Real(0.0)) {
                const Real scale = omega(i, j, k, measure_component);
                for (int component = 0; component < ncomp; ++component) {
                    Real value = Real(0.0);
                    const auto status = remap_detail::checked_quotient(
                        amount(i, j, k, component), scale, value);
                    output(i, j, k, component) = value;
                    if (status != remap_detail::QuotientStatus::Ok)
                        cell_bad = 1;
                }
            }
            bad(i, j, k, 0) = static_cast<Real>(cell_bad);
        });
    }
    const bool locally_valid = invalid.max(0, 0, true) <= Real(0.0);
    if (!locally_valid) {
        diagnostic = "recovering coarse spectrum U=H/omega overflowed or underflowed";
        return false;
    }
    return true;
}

bool form_carrier_relative_state (const MultiFab& spectrum,
                                  const MultiFab& density,
                                  const int density_component,
                                  MultiFab& intensive,
                                  std::string& diagnostic)
{
    MultiFab invalid(spectrum.boxArray(), spectrum.DistributionMap(), 1, 0);
    invalid.setVal(Real(0.0));
    const int ncomp = spectrum.nComp();
    for (amrex::MFIter mfi(spectrum, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.tilebox();
        const auto state = spectrum.const_array(mfi);
        const auto rho = density.const_array(mfi);
        const auto output = intensive.array(mfi);
        const auto bad = invalid.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int cell_bad = 0;
            const Real carrier = rho(i, j, k, density_component);
            for (int component = 0; component < ncomp; ++component) {
                Real value = Real(0.0);
                const auto status = remap_detail::checked_quotient(
                    state(i, j, k, component), carrier, value);
                output(i, j, k, component) = value;
                if (status != remap_detail::QuotientStatus::Ok) cell_bad = 1;
            }
            bad(i, j, k, 0) = static_cast<Real>(cell_bad);
        });
    }
    const bool locally_valid = invalid.max(0, 0, true) <= Real(0.0);
    if (!locally_valid) {
        diagnostic = "forming carrier-relative spectrum z=U/rho_d overflowed or underflowed";
        return false;
    }
    return true;
}

bool reconstruct_carrier_relative_state (const MultiFab& intensive,
                                         const MultiFab& density,
                                         const int density_component,
                                         const MultiFab& measure,
                                         const int measure_component,
                                         MultiFab& candidate,
                                         std::string& diagnostic)
{
    MultiFab invalid(candidate.boxArray(), candidate.DistributionMap(), 1, 0);
    invalid.setVal(Real(0.0));
    const int ncomp = candidate.nComp();
    for (amrex::MFIter mfi(candidate, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.tilebox();
        const auto z = intensive.const_array(mfi);
        const auto rho = density.const_array(mfi);
        const auto omega = measure.const_array(mfi);
        const auto output = candidate.array(mfi);
        const auto bad = invalid.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int cell_bad = 0;
            const Real carrier = rho(i, j, k, density_component);
            const Real scale = omega(i, j, k, measure_component);
            for (int component = 0; component < ncomp; ++component) {
                Real value = Real(0.0);
                Real amount = Real(0.0);
                const auto value_status = remap_detail::checked_product(
                    carrier, z(i, j, k, component), value);
                const auto amount_status = remap_detail::checked_product(
                    scale, value, amount);
                output(i, j, k, component) = value;
                if (value_status != remap_detail::ProductStatus::Ok ||
                    amount_status != remap_detail::ProductStatus::Ok) {
                    cell_bad = 1;
                }
            }
            bad(i, j, k, 0) = static_cast<Real>(cell_bad);
        });
    }
    const bool locally_valid = invalid.max(0, 0, true) <= Real(0.0);
    if (!locally_valid) {
        diagnostic = "reconstructing U=rho_d*z or H=omega*U overflowed or underflowed";
        return false;
    }
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

bool RestrictMappedSpectrum (const SBMLayout& layout,
                             const SBMAMRStateView& fine,
                             const SBMAMRStateView& coarse,
                             const amrex::IntVect& ratio,
                             const int coarse_level,
                             MultiFab& coarse_candidate,
                             std::string& diagnostic)
{
    diagnostic.clear();
    std::string ratio_diagnostic;
    const bool ratio_ok = valid_ratio(ratio, ratio_diagnostic);
    const bool level_ok = coarse_level >= 0;
    std::string fine_structure_diagnostic;
    std::string coarse_structure_diagnostic;
    const bool fine_structure_ok = validate_view_structure(
        fine, layout, "fine", fine_structure_diagnostic);
    const bool coarse_structure_ok = validate_view_structure(
        coarse, layout, "coarse", coarse_structure_diagnostic);
    const bool same_time = fine.spectrum_time == coarse.spectrum_time;
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
    const bool separate_from_fine = candidate_is_separate(coarse_candidate, fine);
    const bool separate_from_coarse = candidate_is_separate(coarse_candidate, coarse);
    const bool local_candidate_ok = candidate_layout_ok && separate_from_fine &&
                                    separate_from_coarse;
    if (!collective_all_true(
            local_candidate_ok,
            "restriction candidate must be separate from both input tuples and match the coarse SBM layout",
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
        fine, layout, coarse_level + 1, "fine", fine_spectrum_diagnostic);
    const bool coarse_spectrum_ok = validate_view_spectrum(
        coarse, layout, coarse_level, "coarse", coarse_spectrum_diagnostic);
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

    MultiFab fine_h(fine.spectrum->boxArray(), fine.spectrum->DistributionMap(),
                    layout.ncomp(), 0);
    MultiFab coarse_h(coarse.spectrum->boxArray(), coarse.spectrum->DistributionMap(),
                      layout.ncomp(), 0);
    const int coverage_component = layout.ncomp();
    const int ncomp = layout.ncomp();
    MultiFab fine_support(fine.spectrum->boxArray(),
                          fine.spectrum->DistributionMap(), layout.ncomp() + 1,
                          0);
    MultiFab coarse_support(coarse.spectrum->boxArray(),
                            coarse.spectrum->DistributionMap(),
                            layout.ncomp() + 1, 0);
    std::string fine_mapped_diagnostic;
    const bool fine_mapped_ok = form_mapped_state(
        *fine.spectrum, *fine.mapped_measure, fine.measure_component, fine_h,
        fine_mapped_diagnostic);
    if (!collective_all_true(
            fine_mapped_ok, fine_mapped_diagnostic,
            "SBM restriction mapped-state formation failed on another MPI rank",
            diagnostic)) {
        return false;
    }

    fine_support.setVal(Real(0.0));
    for (amrex::MFIter mfi(fine_h, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const amrex::Box box = mfi.tilebox();
        const auto mapped = fine_h.const_array(mfi);
        const auto support = fine_support.array(mfi);
        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                for (int component = 0; component < ncomp; ++component) {
                    support(i, j, k, component) =
                        mapped(i, j, k, component) > Real(0.0) ? Real(1.0)
                                                               : Real(0.0);
                }
                support(i, j, k, coverage_component) = Real(1.0);
            });
    }

    // Both averages use the same fine coverage. The extra support components
    // remember whether any child had positive mapped inventory, even if the
    // inventory average itself rounds to zero.
    coarse_h.setVal(Real(0.0));
    coarse_support.setVal(Real(0.0));
    amrex::average_down(fine_h, coarse_h, 0, layout.ncomp(), ratio);
    amrex::average_down(fine_support, coarse_support, 0, layout.ncomp() + 1,
                        ratio);

    MultiFab lost_positive_support(coarse_h.boxArray(),
                                   coarse_h.DistributionMap(), 1, 0);
    lost_positive_support.setVal(Real(0.0));
    for (amrex::MFIter mfi(coarse_h, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const amrex::Box box = mfi.tilebox();
        const auto mapped = coarse_h.const_array(mfi);
        const auto support = coarse_support.const_array(mfi);
        const auto bad = lost_positive_support.array(mfi);
        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                int cell_bad = 0;
                if (support(i, j, k, coverage_component) > Real(0.0)) {
                    for (int component = 0; component < ncomp; ++component) {
                        if (support(i, j, k, component) > Real(0.0) &&
                            mapped(i, j, k, component) == Real(0.0)) {
                            cell_bad = 1;
                        }
                    }
                }
                bad(i, j, k, 0) = static_cast<Real>(cell_bad);
            });
    }
    // The nonnegative bad mask has max zero when clear; an empty rank's local
    // MultiFab maximum is negative, so both cases satisfy this predicate.
    const bool local_support_preserved =
        lost_positive_support.max(0, 0, true) <= Real(0.0);
    if (!collective_all_true(local_support_preserved,
                             "restriction averaging underflow erased positive "
                             "mapped support to zero",
                             "restriction averaging underflow erased positive "
                             "mapped support on another MPI rank",
                             diagnostic)) {
        return false;
    }

    // Build the complete candidate independently. Only cells represented by
    // fine coverage are reconstructed, so uncovered coarse values never take
    // an unnecessary multiply/divide round trip. The public destination is
    // untouched until every rank has admitted this scratch result.
    MultiFab coarse_candidate_scratch(coarse.spectrum->boxArray(),
                                      coarse.spectrum->DistributionMap(),
                                      layout.ncomp(), 0);
    MultiFab::Copy(coarse_candidate_scratch, *coarse.spectrum, 0, 0,
                   layout.ncomp(), 0);
    const bool divided_ok = divide_mapped_state(
        coarse_h, *coarse.mapped_measure, coarse.measure_component,
        coarse_support, coverage_component, coarse_candidate_scratch,
        diagnostic);
    std::string candidate_admission_diagnostic;
    const bool candidate_admissible = authoritative_state_admissible(
        coarse_candidate_scratch, layout, coarse_level,
        &candidate_admission_diagnostic);
    const bool local_candidate_admissible = divided_ok && candidate_admissible;
    const std::string candidate_diagnostic =
        !divided_ok ? diagnostic
                    : !candidate_admissible
                          ? "restricted SBM candidate is inadmissible: " +
                                candidate_admission_diagnostic
                          : std::string{};
    if (!collective_all_true(
            local_candidate_admissible, candidate_diagnostic,
            "restricted SBM candidate is inadmissible on another MPI rank",
            diagnostic)) {
        return false;
    }
    MultiFab::Copy(coarse_candidate, coarse_candidate_scratch, 0, 0,
                   layout.ncomp(), 0);
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
                                     std::string& diagnostic)
{
    diagnostic.clear();
    std::string ratio_diagnostic;
    const bool ratio_ok = valid_ratio(ratio, ratio_diagnostic);
    const bool level_ok = fine_level > 0;

    std::string coarse_structure_diagnostic;
    std::string fine_structure_diagnostic;
    const bool coarse_structure_ok = validate_view_structure(
        coarse, layout, "coarse", coarse_structure_diagnostic);
    const bool fine_structure_ok = validate_view_structure(
        fine_target, layout, "fine target", fine_structure_diagnostic);
    const bool same_time = coarse.spectrum_time == fine_target.spectrum_time;
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
            if (!coarse_geometry.isPeriodic(direction) ||
                !fine_geometry.isPeriodic(direction)) {
                local_candidate_ok = false;
                local_candidate_diagnostic =
                    "M4a carrier-relative prolongation requires periodic "
                    "geometry";
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
        coarse, layout, fine_level - 1, "coarse", coarse_spectrum_diagnostic);
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

    MultiFab coarse_z(coarse.spectrum->boxArray(),
                      coarse.spectrum->DistributionMap(), layout.ncomp(), 0);
    const bool coarse_z_ok = form_carrier_relative_state(
        *coarse.spectrum, *coarse.dry_air_density,
        coarse.density_component, coarse_z, diagnostic);
    if (!collective_all_true(
            coarse_z_ok, diagnostic,
            "forming coarse carrier-relative SBM state failed on another MPI rank",
            diagnostic)) {
        return false;
    }

    MultiFab fine_z(fine_target.spectrum->boxArray(),
                    fine_target.spectrum->DistributionMap(), layout.ncomp(), 0);
    amrex::PhysBCFunctNoOp no_physical_boundary;
    amrex::Vector<amrex::BCRec> bcs(static_cast<std::size_t>(layout.ncomp()));
    // PCInterp currently ignores these records under the periodic-only M4a
    // gate. Keep every face explicit; this is not a nonperiodic BC policy.
    for (auto& bc : bcs) {
        for (int direction = 0; direction < AMREX_SPACEDIM; ++direction) {
            bc.setLo(direction, amrex::BCType::int_dir);
            bc.setHi(direction, amrex::BCType::int_dir);
        }
    }
    amrex::InterpFromCoarseLevel(
        fine_z, amrex::IntVect(0), static_cast<Real>(coarse.spectrum_time),
        coarse_z, 0, 0, layout.ncomp(), coarse_geometry, fine_geometry,
        no_physical_boundary, 0, no_physical_boundary, 0, ratio,
        &amrex::pc_interp, bcs, 0);

    MultiFab fine_candidate_scratch(fine_target.spectrum->boxArray(),
                                    fine_target.spectrum->DistributionMap(),
                                    layout.ncomp(), 0);
    const bool reconstructed_ok = reconstruct_carrier_relative_state(
        fine_z, *fine_target.dry_air_density, fine_target.density_component,
        *fine_target.mapped_measure, fine_target.measure_component,
        fine_candidate_scratch, diagnostic);
    std::string candidate_admission_diagnostic;
    const bool candidate_admissible = authoritative_state_admissible(
        fine_candidate_scratch, layout, fine_level,
        &candidate_admission_diagnostic);
    const bool local_candidate_admissible = reconstructed_ok &&
                                            candidate_admissible;
    const std::string candidate_diagnostic =
        !reconstructed_ok ? diagnostic
                          : !candidate_admissible
                                ? "prolonged SBM candidate is inadmissible: " +
                                      candidate_admission_diagnostic
                                : std::string{};
    if (!collective_all_true(
            local_candidate_admissible, candidate_diagnostic,
            "prolonged SBM candidate is inadmissible on another MPI rank",
            diagnostic)) {
        return false;
    }
    MultiFab::Copy(fine_candidate, fine_candidate_scratch, 0, 0, layout.ncomp(),
                   0);
    return true;
}

} // namespace erf_sbm
