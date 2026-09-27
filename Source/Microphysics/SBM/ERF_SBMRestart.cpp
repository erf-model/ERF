#include "ERF_SBMRestart.H"

#include "ERF_SBMConstraintGroups.H"
#include <AMReX_Arena.H>
#include <AMReX_BoxIterator.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFabUtil.H>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace erf_sbm {

std::string restart_schema(const SBMLayout& layout)
{
    return std::string("ERF-SBM-RESTART-M1-v1\n") +
           "layout=" + layout.schema_identity() + "\n" +
           "constraint-policy=nonnegative-bin-mass-and-moments-v1\n" +
           "projection=liquid-mass-sum-to-qc-qr-v1\n" +
           "representation=bin-mass-density-v1\n" +
           "transport=zero-transport-fixture-v1\n";
}

bool restart_schema_matches(const SBMLayout& layout, const std::string& persisted)
{
    return restart_schema(layout) == persisted;
}

bool authoritative_state_admissible(const amrex::MultiFab& spectrum,
                                    const SBMLayout& layout,
                                    const int level,
                                    std::string* diagnostic)
{
    auto reject = [diagnostic](std::string message) {
        if (diagnostic) *diagnostic = std::move(message);
        return false;
    };
    if (diagnostic) diagnostic->clear();
    if (level < 0) return reject("SBM authoritative restart state has an invalid level");
    if (spectrum.nComp() != layout.ncomp()) {
        return reject("SBM authoritative restart state component count does not match its layout");
    }

    const auto groups = make_constraint_groups(layout);
    std::vector<amrex::Real> state(static_cast<std::size_t>(layout.ncomp()));
    const int precision = std::numeric_limits<amrex::Real>::max_digits10;

    for (amrex::MFIter mfi(spectrum); mfi.isValid(); ++mfi) {
        const amrex::FArrayBox& source_fab = spectrum[mfi];
        const amrex::FArrayBox* host_source = &source_fab;
        std::unique_ptr<amrex::FArrayBox> host_fab;
#ifdef AMREX_USE_GPU
        if (source_fab.arena()->isManaged() || source_fab.arena()->isDevice()) {
            host_fab = std::make_unique<amrex::FArrayBox>(
                source_fab.box(), source_fab.nComp(), amrex::The_Pinned_Arena());
            amrex::Gpu::dtoh_memcpy_async(
                host_fab->dataPtr(), source_fab.dataPtr(),
                static_cast<std::size_t>(source_fab.size()) * sizeof(amrex::Real));
            amrex::Gpu::streamSynchronize();
            host_source = host_fab.get();
        }
#endif
        const auto values = host_source->const_array();
        for (amrex::BoxIterator iterator(mfi.validbox()); iterator.ok(); ++iterator) {
            const amrex::IntVect cell = iterator();
            for (int component = 0; component < layout.ncomp(); ++component) {
                const amrex::Real value = values(cell, component);
                if (!std::isfinite(value)) {
                    std::ostringstream message;
                    message << "SBM authoritative restart state is inadmissible"
                            << ": level=" << level << ", cell=(";
                    for (int direction = 0; direction < AMREX_SPACEDIM; ++direction) {
                        if (direction != 0) message << ',';
                        message << cell[direction];
                    }
                    message << "), component=" << component
                            << ", constraint=finite, value="
                            << std::setprecision(precision) << value;
                    const auto group = std::find_if(groups.begin(), groups.end(),
                        [component](const ConstraintGroup& candidate) {
                            return candidate.contains(component);
                        });
                    if (group != groups.end()) {
                        message << ", population=" << group->population_id
                                << " (" << group->semantic_id << ")"
                                << ", bin=" << group->bin;
                    }
                    return reject(message.str());
                }
                state[static_cast<std::size_t>(component)] = value;
            }

            for (const auto& group : groups) {
                amrex::Real margin = amrex::Real(0.0);
                std::string failed_constraint;
                if (group.admissible(state, &margin, &failed_constraint)) continue;

                std::ostringstream message;
                message << "SBM authoritative restart state is inadmissible"
                        << ": level=" << level << ", cell=(";
                for (int direction = 0; direction < AMREX_SPACEDIM; ++direction) {
                    if (direction != 0) message << ',';
                    message << cell[direction];
                }
                message << "), population=" << group.population_id
                        << " (" << group.semantic_id << ")"
                        << ", bin=" << group.bin
                        << ", constraint=" << failed_constraint
                        << ", margin=" << std::setprecision(precision) << margin
                        << ", components={";
                for (std::size_t member = 0; member < group.members.size(); ++member) {
                    if (member != 0) message << ',';
                    const int component = group.members[member];
                    message << component << ':'
                            << std::setprecision(precision)
                            << state[static_cast<std::size_t>(component)];
                }
                message << '}';
                return reject(message.str());
            }
        }
    }
    return true;
}

bool restart_projection_matches(const amrex::MultiFab& spectrum,
                                const amrex::MultiFab& persisted_core,
                                const SBMBulkProjection& projection,
                                const int qc_component, const int qr_component,
                                const amrex::Real tolerance_scale)
{
    if (spectrum.boxArray() != persisted_core.boxArray() ||
        spectrum.DistributionMap() != persisted_core.DistributionMap() ||
        qc_component < 0 || qr_component < 0 ||
        qc_component >= persisted_core.nComp() || qr_component >= persisted_core.nComp() ||
        !std::isfinite(tolerance_scale) || !(tolerance_scale >= amrex::Real(0.0)) ||
        !spectrum.is_finite(0, spectrum.nComp(), 0) ||
        !persisted_core.is_finite(qc_component, 1, 0) ||
        !persisted_core.is_finite(qr_component, 1, 0)) {
        return false;
    }

    amrex::MultiFab expected(persisted_core.boxArray(), persisted_core.DistributionMap(), 2, 0);
    amrex::MultiFab stored(persisted_core.boxArray(), persisted_core.DistributionMap(), 2, 0);
    expected.setVal(amrex::Real(0.0));
    for (amrex::MFIter mfi(expected); mfi.isValid(); ++mfi) {
        projection.apply_to_core(mfi.validbox(), spectrum.const_array(mfi),
                                 expected.array(mfi), 0, 1);
    }
    if (!expected.is_finite(0, 2, 0)) return false;
    amrex::MultiFab::Copy(stored, persisted_core, qc_component, 0, 1, 0);
    amrex::MultiFab::Copy(stored, persisted_core, qr_component, 1, 1, 0);
    if (!stored.is_finite(0, 2, 0)) return false;

    const amrex::Real scale = std::max({expected.norm0(0), expected.norm0(1),
                                        stored.norm0(0), stored.norm0(1)});
    amrex::MultiFab::Subtract(expected, stored, 0, 0, 2, 0);
    if (!expected.is_finite(0, 2, 0)) return false;
    const amrex::Real difference = std::max(expected.norm0(0), expected.norm0(1));
    const amrex::Real tolerance = tolerance_scale *
        std::numeric_limits<amrex::Real>::epsilon() *
        std::max(scale, std::numeric_limits<amrex::Real>::min());
    return std::isfinite(scale) && std::isfinite(difference) &&
           std::isfinite(tolerance) && difference <= tolerance;
}

} // namespace erf_sbm
