#include "ERF_SBMRestart.H"

#include "ERF_SBMConstraintGroups.H"
#include "ERF_SBMRemapping.H"
#include <AMReX_Arena.H>
#include <AMReX_BoxIterator.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace erf_sbm {

struct SBMStateAdmission::Impl {
    explicit Impl (const SBMLayout& layout)
        : layout(layout), ncomp(layout.ncomp()), groups(make_constraint_groups(layout)),
          flattened(make_constraint_descriptors(layout))
    {
        constraints.resize(flattened.constraints.size());
        terms.resize(flattened.terms.size());
        if (!flattened.constraints.empty()) {
            amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                             flattened.constraints.begin(),
                             flattened.constraints.end(), constraints.begin());
        }
        if (!flattened.terms.empty()) {
            amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                             flattened.terms.begin(), flattened.terms.end(),
                             terms.begin());
        }

        std::vector<amrex::Real> host_edges;
        std::vector<amrex::Real> host_pivots;
        std::vector<int> host_properties;
        std::vector<int> edge_offsets, pivot_offsets, property_offsets;
        for (const auto& population : layout.populations()) {
            edge_offsets.push_back(static_cast<int>(host_edges.size()));
            host_edges.insert(host_edges.end(), population.grid.edges().begin(),
                              population.grid.edges().end());
            pivot_offsets.push_back(static_cast<int>(host_pivots.size()));
            host_pivots.insert(host_pivots.end(), population.grid.pivots().begin(),
                               population.grid.pivots().end());
            property_offsets.push_back(static_cast<int>(host_properties.size()));
            host_properties.insert(host_properties.end(),
                population.property_component_offsets.begin(),
                population.property_component_offsets.end());
        }
        edges.resize(host_edges.size());
        pivots.resize(host_pivots.size());
        properties.resize(host_properties.size());
        if (!host_edges.empty()) {
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, host_edges.begin(),
                             host_edges.end(), edges.begin());
        }
        if (!host_pivots.empty()) {
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, host_pivots.begin(),
                             host_pivots.end(), pivots.begin());
        }
        if (!host_properties.empty()) {
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, host_properties.begin(),
                             host_properties.end(), properties.begin());
        }

        std::vector<PopulationRemapView> views;
        views.reserve(layout.populations().size());
        for (std::size_t p = 0; p < layout.populations().size(); ++p) {
            auto view = population_remap_view(
                layout, layout.populations()[p].population_id);
            view.edges = edges.data() + edge_offsets[p];
            view.pivots = layout.populations()[p].grid.pivots().empty()
                ? nullptr : pivots.data() + pivot_offsets[p];
            view.property_component_offsets =
                layout.populations()[p].property_component_offsets.empty()
                ? nullptr : properties.data() + property_offsets[p];
            views.push_back(view);
        }
        populations.resize(views.size());
        if (!views.empty()) {
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, views.begin(),
                             views.end(), populations.begin());
        }
    }

    SBMLayout layout;
    int ncomp{0};
    std::vector<ConstraintGroup> groups;
    FlattenedConstraintSet flattened;
    amrex::Gpu::DeviceVector<ConstraintDescriptor> constraints;
    amrex::Gpu::DeviceVector<ConstraintTerm> terms;
    amrex::Gpu::DeviceVector<PopulationRemapView> populations;
    amrex::Gpu::DeviceVector<amrex::Real> edges;
    amrex::Gpu::DeviceVector<amrex::Real> pivots;
    amrex::Gpu::DeviceVector<int> properties;
    std::vector<std::unique_ptr<amrex::MultiFab>> invalid;
    mutable Counters counts;
};

SBMStateAdmission::SBMStateAdmission (const SBMLayout& layout)
    : m_impl(std::make_unique<Impl>(layout)) {}

SBMStateAdmission::~SBMStateAdmission () = default;

void SBMStateAdmission::define_level (
    const int level, const amrex::BoxArray& grids,
    const amrex::DistributionMapping& mapping)
{
    if (level < 0) throw std::invalid_argument("SBM admission level is negative");
    if (m_impl->invalid.size() <= static_cast<std::size_t>(level)) {
        m_impl->invalid.resize(static_cast<std::size_t>(level) + 1U);
    }
    if (m_impl->invalid[static_cast<std::size_t>(level)]) {
        throw std::logic_error("SBM admission level is already defined");
    }
    m_impl->invalid[static_cast<std::size_t>(level)] =
        std::make_unique<amrex::MultiFab>(grids, mapping, 1, 0);
}

void SBMStateAdmission::destroy_level (const int level)
{
    if (level >= 0 && static_cast<std::size_t>(level) < m_impl->invalid.size()) {
        m_impl->invalid[static_cast<std::size_t>(level)].reset();
    }
}

bool SBMStateAdmission::admissible (const amrex::MultiFab& spectrum,
                                    const int level,
                                    std::string& diagnostic) const
{
    diagnostic.clear();
    const bool level_defined = level >= 0 &&
        static_cast<std::size_t>(level) < m_impl->invalid.size() &&
        m_impl->invalid[static_cast<std::size_t>(level)] != nullptr;
    const bool local_layout_ok = level_defined &&
        spectrum.nComp() == m_impl->ncomp &&
        spectrum.boxArray() == m_impl->invalid[static_cast<std::size_t>(level)]->boxArray() &&
        spectrum.DistributionMap() ==
            m_impl->invalid[static_cast<std::size_t>(level)]->DistributionMap();
    int any_layout_bad = local_layout_ok ? 0 : 1;
    ++m_impl->counts.mpi_reductions;
    amrex::ParallelDescriptor::ReduceIntMax(any_layout_bad);
    if (any_layout_bad != 0) {
        diagnostic = "SBM device admission candidate does not match a defined level layout";
        return false;
    }
    ++m_impl->counts.device_admissions;
    auto& invalid = *m_impl->invalid[static_cast<std::size_t>(level)];
    invalid.setVal(amrex::Real(0.0));
    const auto* descriptors = m_impl->constraints.data();
    const auto* terms = m_impl->terms.data();
    const auto* populations = m_impl->populations.data();
    const int descriptor_count = static_cast<int>(m_impl->flattened.constraints.size());
    const int population_count = static_cast<int>(m_impl->populations.size());
    for (amrex::MFIter mfi(spectrum, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.tilebox();
        const auto state = spectrum.const_array(mfi);
        const auto bad = invalid.array(mfi);
        const int ncomp = m_impl->ncomp;
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int rejected = 0;
            for (int component = 0; component < ncomp; ++component) {
                if (!amrex::Math::isfinite(state(i, j, k, component))) rejected = 1;
            }
            const amrex::Real roundoff = amrex::Real(128.0) *
                std::numeric_limits<amrex::Real>::epsilon();
            for (int ci = 0; ci < descriptor_count; ++ci) {
                const auto& descriptor = descriptors[ci];
                amrex::Real margin = amrex::Real(0.0);
                amrex::Real scale = amrex::Real(0.0);
                for (int t = 0; t < descriptor.term_count; ++t) {
                    const auto& term = terms[descriptor.term_offset + t];
                    const amrex::Real value = state(i, j, k, term.component);
                    margin = std::fma(term.coefficient, value, margin);
                    scale += amrex::Math::abs(term.coefficient * value);
                }
                if (!amrex::Math::isfinite(margin) || margin < -roundoff * scale) {
                    rejected = 1;
                }
            }
            const remap_detail::RemapCellStateView cell_state{state, i, j, k};
            for (int p = 0; p < population_count; ++p) {
                const auto& population = populations[p];
                for (int bin = 0; bin < population.nbins; ++bin) {
                    if (!remap_detail::canonical_persisted_bin_state(
                            population, bin, cell_state, ncomp)) rejected = 1;
                }
            }
            bad(i, j, k, 0) = static_cast<amrex::Real>(rejected);
        });
    }

    int any_device_bad = invalid.max(0, 0, true) > amrex::Real(0.0) ? 1 : 0;
    ++m_impl->counts.mpi_reductions;
    amrex::ParallelDescriptor::ReduceIntMax(any_device_bad);
    const bool device_admissible = any_device_bad == 0;
    if (device_admissible) {
        ++m_impl->counts.successful_admissions;
        return true;
    }

    ++m_impl->counts.host_scans;
    std::uint64_t fab_bytes = 0;
    std::uint64_t d2h_bytes = 0;
    for (amrex::MFIter mfi(spectrum); mfi.isValid(); ++mfi) {
        const auto bytes = static_cast<std::uint64_t>(spectrum[mfi].size()) *
                           sizeof(amrex::Real);
        fab_bytes += bytes;
#ifdef AMREX_USE_GPU
        if (spectrum[mfi].arena()->isManaged() ||
            spectrum[mfi].arena()->isDevice()) d2h_bytes += bytes;
#endif
    }
    m_impl->counts.host_scan_fab_bytes += fab_bytes;
    m_impl->counts.host_d2h_bytes += d2h_bytes;
    std::string host_diagnostic;
    const bool local_host_admissible = authoritative_state_admissible(
        spectrum, m_impl->layout, level, &host_diagnostic);
    int any_host_bad = local_host_admissible ? 0 : 1;
    ++m_impl->counts.mpi_reductions;
    amrex::ParallelDescriptor::ReduceIntMax(any_host_bad);
    const bool host_admissible = any_host_bad == 0;
    if (device_admissible != host_admissible) {
        diagnostic = "SBM device and detailed host candidate-admission decisions disagree";
        return false;
    }
    diagnostic = host_diagnostic.empty()
        ? "SBM full candidate failed device canonical/group admission on another rank"
        : std::move(host_diagnostic);
    return false;
}

SBMStateAdmission::Counters SBMStateAdmission::counters () const noexcept
{
    return m_impl->counts;
}

std::string restart_schema (const SBMLayout& layout)
{
    std::ostringstream schema;
    schema << "ERF-SBM-RESTART-v1\n"
           << "layout=" << layout.schema_identity() << '\n'
           << "constraint-policy=linear-groups-plus-canonical-persisted-v1\n"
           << "projection=liquid-mass-sum-to-qc-qr-v1\n"
           << "transport=sbm-m3-mapped-group-fct-v1\n"
           << "high-order=endpoint-number-wenoz3-v1\n";
    return schema.str();
}

bool restart_schema_matches (const SBMLayout& layout, const std::string& persisted)
{
    return restart_schema(layout) == persisted;
}

bool authoritative_state_admissible (const amrex::MultiFab& spectrum,
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
    std::vector<PopulationRemapView> remap_views;
    remap_views.reserve(layout.populations().size());
    for (const auto& population : layout.populations()) {
        remap_views.push_back(population_remap_view(layout, population.population_id));
    }
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

            for (std::size_t population_index = 0;
                 population_index < remap_views.size(); ++population_index) {
                const auto& view = remap_views[population_index];
                if (view.moment_mode != MomentMode::TwoMoment) continue;
                for (int bin = 0; bin < view.nbins; ++bin) {
                    if (remap_detail::canonical_two_moment_bin_state(
                            view, bin, state.data(), static_cast<int>(state.size()))) continue;
                    std::ostringstream message;
                    message << "SBM authoritative restart state is inadmissible"
                            << ": level=" << level << ", cell=(";
                    for (int direction = 0; direction < AMREX_SPACEDIM; ++direction) {
                        if (direction != 0) message << ',';
                        message << cell[direction];
                    }
                    message << "), population=" << view.population_id
                            << ", bin=" << bin
                            << ", constraint=canonical-two-moment-bin-state";
                    return reject(message.str());
                }
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

            for (std::size_t population_index = 0;
                 population_index < remap_views.size(); ++population_index) {
                const auto& view = remap_views[population_index];
                if (view.moment_mode != MomentMode::OneMoment) continue;
                for (int bin = 0; bin < view.nbins; ++bin) {
                    if (remap_detail::canonical_persisted_bin_state(
                            view, bin, state.data(), static_cast<int>(state.size()))) continue;
                    std::ostringstream message;
                    message << "SBM authoritative restart state is inadmissible"
                            << ": level=" << level << ", cell=(";
                    for (int direction = 0; direction < AMREX_SPACEDIM; ++direction) {
                        if (direction != 0) message << ',';
                        message << cell[direction];
                    }
                    message << "), population=" << view.population_id
                            << ", bin=" << bin
                            << ", constraint=canonical-one-moment-bin-state";
                    return reject(message.str());
                }
            }
        }
    }
    return true;
}

bool restart_projection_matches (const amrex::MultiFab& spectrum,
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
