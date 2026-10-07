#include "ERF_SBMAMRTransfer.H"

#include "AuxiliaryState/ERF_AuxiliaryMappedTransport.H"
#include "ERF_SBMRemapping.H"
#include "ERF_SBMRestart.H"

#include <AMReX_FillPatchUtil.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFabUtil.H>

#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

namespace erf_sbm {
namespace {

using amrex::MultiFab;
using amrex::Real;

bool validate_view (const SBMAMRStateView& view, const SBMLayout& layout,
                    const int level, const char* label,
                    const bool validate_spectrum, std::string& diagnostic)
{
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
    if (!erf_auxiliary::ValidatePositiveFiniteComponent(
            *view.dry_air_density, view.density_component, diagnostic)) {
        diagnostic = std::string(label) + " dry-air density: " + diagnostic;
        return false;
    }
    if (!erf_auxiliary::ValidatePositiveFiniteComponent(
            *view.mapped_measure, view.measure_component, diagnostic)) {
        diagnostic = std::string(label) + " mapped measure: " + diagnostic;
        return false;
    }
    if (validate_spectrum &&
        !authoritative_state_admissible(*view.spectrum, layout, level,
                                        &diagnostic)) {
        diagnostic = std::string(label) + " authoritative spectrum: " + diagnostic;
        return false;
    }
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
    if (invalid.max(0) != Real(0.0)) {
        diagnostic = "forming mapped spectral amount H=omega*U overflowed or underflowed";
        return false;
    }
    return true;
}

bool divide_mapped_state (const MultiFab& mapped, const MultiFab& measure,
                         const int measure_component, MultiFab& candidate,
                         std::string& diagnostic)
{
    MultiFab invalid(mapped.boxArray(), mapped.DistributionMap(), 1, 0);
    invalid.setVal(Real(0.0));
    const int ncomp = mapped.nComp();
    for (amrex::MFIter mfi(mapped, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.tilebox();
        const auto amount = mapped.const_array(mfi);
        const auto omega = measure.const_array(mfi);
        const auto output = candidate.array(mfi);
        const auto bad = invalid.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int cell_bad = 0;
            const Real scale = omega(i, j, k, measure_component);
            for (int component = 0; component < ncomp; ++component) {
                Real value = Real(0.0);
                const auto status = remap_detail::checked_quotient(
                    amount(i, j, k, component), scale, value);
                output(i, j, k, component) = value;
                if (status != remap_detail::QuotientStatus::Ok) cell_bad = 1;
            }
            bad(i, j, k, 0) = static_cast<Real>(cell_bad);
        });
    }
    if (invalid.max(0) != Real(0.0)) {
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
    if (invalid.max(0) != Real(0.0)) {
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
    MultiFab mapped(candidate.boxArray(), candidate.DistributionMap(),
                    candidate.nComp(), 0);
    MultiFab invalid(candidate.boxArray(), candidate.DistributionMap(), 1, 0);
    invalid.setVal(Real(0.0));
    const int ncomp = candidate.nComp();
    for (amrex::MFIter mfi(candidate, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.tilebox();
        const auto z = intensive.const_array(mfi);
        const auto rho = density.const_array(mfi);
        const auto omega = measure.const_array(mfi);
        const auto output = candidate.array(mfi);
        const auto mapped_output = mapped.array(mfi);
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
                mapped_output(i, j, k, component) = amount;
                if (value_status != remap_detail::ProductStatus::Ok ||
                    amount_status != remap_detail::ProductStatus::Ok) {
                    cell_bad = 1;
                }
            }
            bad(i, j, k, 0) = static_cast<Real>(cell_bad);
        });
    }
    if (invalid.max(0) != Real(0.0)) {
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
    if (!valid_ratio(ratio, diagnostic)) return false;
    if (!validate_view(fine, layout, coarse_level + 1, "fine", true, diagnostic) ||
        !validate_view(coarse, layout, coarse_level, "coarse", true, diagnostic)) {
        return false;
    }
    if (fine.spectrum_time != coarse.spectrum_time) {
        diagnostic = "restriction requires fine and coarse state tuples at the same semantic time";
        return false;
    }
    if (coarse_level < 0 || coarse_candidate.nComp() != layout.ncomp() ||
        !erf_auxiliary::SameCellLayout(coarse_candidate, *coarse.spectrum) ||
        !candidate_is_separate(coarse_candidate, fine) ||
        !candidate_is_separate(coarse_candidate, coarse)) {
        diagnostic = "restriction candidate must be a separate coarse field with the SBM layout";
        return false;
    }

    MultiFab fine_h(fine.spectrum->boxArray(), fine.spectrum->DistributionMap(),
                    layout.ncomp(), 0);
    MultiFab coarse_h(coarse.spectrum->boxArray(), coarse.spectrum->DistributionMap(),
                      layout.ncomp(), 0);
    if (!form_mapped_state(*fine.spectrum, *fine.mapped_measure,
                           fine.measure_component, fine_h, diagnostic) ||
        !form_mapped_state(*coarse.spectrum, *coarse.mapped_measure,
                           coarse.measure_component, coarse_h, diagnostic)) {
        return false;
    }

    // Coarse H starts from the existing coarse state so cells not covered by
    // fine boxes retain their original values. average_down touches only the
    // covered coarse region.
    amrex::average_down(fine_h, coarse_h, 0, layout.ncomp(), ratio);
    if (!divide_mapped_state(coarse_h, *coarse.mapped_measure,
                             coarse.measure_component, coarse_candidate,
                             diagnostic)) {
        return false;
    }
    if (!authoritative_state_admissible(coarse_candidate, layout, coarse_level,
                                        &diagnostic)) {
        diagnostic = "restricted SBM candidate is inadmissible: " + diagnostic;
        return false;
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
                                     std::string& diagnostic)
{
    diagnostic.clear();
    if (!valid_ratio(ratio, diagnostic)) return false;
    if (!validate_view(coarse, layout, fine_level - 1, "coarse", true, diagnostic) ||
        !validate_view(fine_target, layout, fine_level, "fine target", false, diagnostic)) {
        return false;
    }
    if (fine_level <= 0 || coarse.spectrum_time != fine_target.spectrum_time) {
        diagnostic = "prolongation requires a valid fine level and same-time source/target tuples";
        return false;
    }
    if (fine_candidate.nComp() != layout.ncomp() ||
        !erf_auxiliary::SameCellLayout(fine_candidate, *fine_target.spectrum) ||
        !candidate_is_separate(fine_candidate, coarse) ||
        !candidate_is_separate(fine_candidate, fine_target)) {
        diagnostic = "prolongation candidate must be a separate fine field with the SBM layout";
        return false;
    }
    for (int direction = 0; direction < AMREX_SPACEDIM; ++direction) {
        if (!coarse_geometry.isPeriodic(direction) ||
            !fine_geometry.isPeriodic(direction) ||
            coarse_geometry.isPeriodic(direction) != fine_geometry.isPeriodic(direction)) {
            diagnostic = "M4a carrier-relative prolongation requires matching periodic geometry";
            return false;
        }
    }
    amrex::Box refined_coarse_domain = coarse_geometry.Domain();
    refined_coarse_domain.refine(ratio);
    if (refined_coarse_domain != fine_geometry.Domain()) {
        diagnostic = "prolongation refinement ratio does not map the coarse domain to the fine domain";
        return false;
    }

    MultiFab coarse_z(coarse.spectrum->boxArray(),
                      coarse.spectrum->DistributionMap(), layout.ncomp(), 0);
    if (!form_carrier_relative_state(*coarse.spectrum,
                                     *coarse.dry_air_density,
                                     coarse.density_component, coarse_z,
                                     diagnostic)) {
        return false;
    }

    MultiFab fine_z(fine_target.spectrum->boxArray(),
                    fine_target.spectrum->DistributionMap(), layout.ncomp(), 0);
    amrex::PhysBCFunctNoOp no_physical_boundary;
    amrex::Vector<amrex::BCRec> bcs(static_cast<std::size_t>(layout.ncomp()));
    amrex::InterpFromCoarseLevel(
        fine_z, amrex::IntVect(0), static_cast<Real>(coarse.spectrum_time),
        coarse_z, 0, 0, layout.ncomp(), coarse_geometry, fine_geometry,
        no_physical_boundary, 0, no_physical_boundary, 0, ratio,
        &amrex::pc_interp, bcs, 0);

    if (!reconstruct_carrier_relative_state(
            fine_z, *fine_target.dry_air_density,
            fine_target.density_component, *fine_target.mapped_measure,
            fine_target.measure_component, fine_candidate, diagnostic)) {
        return false;
    }
    if (!authoritative_state_admissible(fine_candidate, layout, fine_level,
                                        &diagnostic)) {
        diagnostic = "prolonged SBM candidate is inadmissible: " + diagnostic;
        return false;
    }
    return true;
}

} // namespace erf_sbm
