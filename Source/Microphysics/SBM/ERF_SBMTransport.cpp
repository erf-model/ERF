#include "ERF_SBMTransport.H"

#include "ERF_SBMRemapping.H"
#include "ERF_SBMRestart.H"
#include "AuxiliaryState/ERF_AuxiliaryMappedTransport.H"
#include "ERF_IndexDefines.H"
#include "Advection/ERF_AdvectionSrcForScalars.H"

#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFabUtil.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace erf_sbm {
namespace {

using amrex::Array4;
using amrex::GpuArray;
using amrex::MultiFab;
using amrex::Real;
using FaceArraySet = GpuArray<Array4<const Real>, AMREX_SPACEDIM>;

struct DeviceGroupInfo
{
    int constraints_begin{0};
    int linear_constraint_count{0};
    int supports_begin{0};
    int support_count{0};
    int chunk_index{0};
    int local_group_index{0};
    int ratio_count{0};
    int ratio_offset{0};
    int canonical_reserve{0};
    int canonical_mass_local{-1};
    int canonical_number_local{-1};
    Real canonical_mass_coefficient{Real(0.0)};
    Real canonical_number_coefficient{Real(0.0)};
};

struct LocalPropertySupport
{
    int property_global{-1};
    int carrier_global{-1};
    int property_local{-1};
    int carrier_local{-1};
    Real hard_min{Real(0.0)};
    Real hard_max{Real(0.0)};
    int has_hard_min{0};
    int has_hard_max{0};
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
clamp_unit (const Real value) noexcept
{
    return amrex::max(Real(0.0), amrex::min(Real(1.0), value));
}

void
build_intensive (MultiFab& intensive,
                 const MultiFab& state,
                 const MultiFab& conserved,
                 const amrex::Periodicity& periodicity,
                 const char* what)
{
    AMREX_ALWAYS_ASSERT(intensive.nComp() == state.nComp());
    AMREX_ALWAYS_ASSERT(conserved.nComp() > Rho_comp);
    AMREX_ALWAYS_ASSERT(erf_auxiliary::SameCellLayout(intensive, state));
    AMREX_ALWAYS_ASSERT(erf_auxiliary::SameCellLayout(intensive, conserved));
    std::string diagnostic;
    if (!erf_auxiliary::ValidatePositiveFiniteComponent(conserved, Rho_comp,
                                                        diagnostic)) {
        amrex::Abort(std::string("SBM ") + what + " density: " + diagnostic);
    }
    for (amrex::MFIter mfi(intensive, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto u = state.const_array(mfi);
        const auto rho = conserved.const_array(mfi);
        const auto out = intensive.array(mfi);
        const int ncomp = state.nComp();
        amrex::ParallelFor(
            bx, ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                out(i, j, k, n) = u(i, j, k, n) / rho(i, j, k, Rho_comp);
            });
    }
    intensive.FillBoundary(periodicity);
    if (!intensive.is_finite(0, intensive.nComp(), 0)) {
        amrex::Abort(std::string("SBM ") + what +
                     " intensive state is nonfinite");
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
face_correction (const erf_sbm::ConstraintDescriptor& descriptor,
                 const erf_sbm::ConstraintTerm* terms,
                 const FaceArraySet& low,
                 const FaceArraySet& high,
                 const int i,
                 const int j,
                 const int k,
                 const int dir,
                 const int side,
                 const Real tau,
                 const GpuArray<Real, AMREX_SPACEDIM>& dx_inv) noexcept
{
    const int fi = i + (dir == 0 && side > 0 ? 1 : 0);
    const int fj = j + (dir == 1 && side > 0 ? 1 : 0);
    const int fk = k + (dir == 2 && side > 0 ? 1 : 0);
    Real g_flux = Real(0.0);
    for (int t = 0; t < descriptor.term_count; ++t) {
        const auto& term = terms[descriptor.term_offset + t];
        g_flux = std::fma(term.coefficient,
                          high[dir](fi, fj, fk, term.component) -
                              low[dir](fi, fj, fk, term.component),
                          g_flux);
    }
    const Real outward_sign = side > 0 ? Real(1.0) : Real(-1.0);
    return -outward_sign * tau * dx_inv[dir] * g_flux;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
canonical_reserve_face_correction (const int mass_component,
                                   const int number_component,
                                   const Real mass_coefficient,
                                   const Real number_coefficient,
                                   const FaceArraySet& low,
                                   const FaceArraySet& high,
                                   const int i,
                                   const int j,
                                   const int k,
                                   const int dir,
                                   const int side,
                                   const Real tau,
                                   const GpuArray<Real, AMREX_SPACEDIM>& dx_inv)
    noexcept
{
    const int fi = i + (dir == 0 && side > 0 ? 1 : 0);
    const int fj = j + (dir == 1 && side > 0 ? 1 : 0);
    const int fk = k + (dir == 2 && side > 0 ? 1 : 0);
    const Real dmass = high[dir](fi, fj, fk, mass_component) -
                       low[dir](fi, fj, fk, mass_component);
    const Real dnumber = high[dir](fi, fj, fk, number_component) -
                         low[dir](fi, fj, fk, number_component);
    const Real g_flux = std::fma(mass_coefficient, dmass,
                                number_coefficient * dnumber);
    const Real outward_sign = side > 0 ? Real(1.0) : Real(-1.0);
    return -outward_sign * tau * dx_inv[dir] * g_flux;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
property_face_correction (const LocalPropertySupport& support,
                          const FaceArraySet& low,
                          const FaceArraySet& high,
                          const int i,
                          const int j,
                          const int k,
                          const int dir,
                          const int side,
                          const Real tau,
                          const GpuArray<Real, AMREX_SPACEDIM>& dx_inv,
                          const Real slope,
                          const bool upper) noexcept
{
    const int fi = i + (dir == 0 && side > 0 ? 1 : 0);
    const int fj = j + (dir == 1 && side > 0 ? 1 : 0);
    const int fk = k + (dir == 2 && side > 0 ? 1 : 0);
    const Real dproperty = high[dir](fi, fj, fk, support.property_local) -
                           low[dir](fi, fj, fk, support.property_local);
    const Real dcarrier = high[dir](fi, fj, fk, support.carrier_local) -
                          low[dir](fi, fj, fk, support.carrier_local);
    const Real g_flux =
        upper ? slope * dcarrier - dproperty : dproperty - slope * dcarrier;
    const Real outward_sign = side > 0 ? Real(1.0) : Real(-1.0);
    return -outward_sign * tau * dx_inv[dir] * g_flux;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
consider_donor_ratio (const Array4<const Real>& field,
                      const int i,
                      const int j,
                      const int k,
                      const LocalPropertySupport& support,
                      Real& lower,
                      Real& upper,
                      int& found) noexcept
{
    const Real carrier = field(i, j, k, support.carrier_global);
    const Real property = field(i, j, k, support.property_global);
    if (!amrex::Math::isfinite(carrier) || !amrex::Math::isfinite(property) ||
        carrier < Real(0.0) || property < Real(0.0)) {
        return;
    }
    if (carrier == Real(0.0)) {
        return;
    }
    const Real ratio = property / carrier;
    if (!amrex::Math::isfinite(ratio) || ratio < Real(0.0)) {
        return;
    }
    lower = found == 0 ? ratio : amrex::min(lower, ratio);
    upper = found == 0 ? ratio : amrex::max(upper, ratio);
    found = 1;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
include_neighbourhood (const Array4<const Real>& field,
                       const int i,
                       const int j,
                       const int k,
                       const LocalPropertySupport& support,
                       Real& lower,
                       Real& upper,
                       int& found) noexcept
{
    consider_donor_ratio(field, i, j, k, support, lower, upper, found);
    consider_donor_ratio(field, i - 1, j, k, support, lower, upper, found);
    consider_donor_ratio(field, i + 1, j, k, support, lower, upper, found);
    consider_donor_ratio(field, i, j - 1, k, support, lower, upper, found);
    consider_donor_ratio(field, i, j + 1, k, support, lower, upper, found);
    consider_donor_ratio(field, i, j, k - 1, support, lower, upper, found);
    consider_donor_ratio(field, i, j, k + 1, support, lower, upper, found);
}

} // namespace

struct SBMTransport::LevelStorage
{
    MultiFab anchor;
    MultiFab target;
    MultiFab input_intensive;
    MultiFab physical_input_intensive;
    MultiFab base_intensive;
    MultiFab anchor_intensive;
    MultiFab low_trial_h;
    MultiFab cell_ratios;
    MultiFab invalid;
    MultiFab outgoing_demand;
    MultiFab measure;
    erf_auxiliary::MappedFaceFluxRate low_rate;
    erf_auxiliary::MappedFaceFluxRate high_rate;
    erf_auxiliary::MappedFaceFluxRate accepted_rate;
    erf_auxiliary::MappedFaceFluxRate face_lambda;
    erf_auxiliary::MappedFaceFluxRate projected_rate;
    erf_auxiliary::CompletedStepFluxLedger projected_ledger;
    erf_auxiliary::CompletedStepFluxLedger spectral_ledger;
    bool measure_ready{false};
};

struct SBMTransport::DeviceMetadata
{
    std::vector<DeviceGroupInfo> host_groups;
    amrex::Gpu::DeviceVector<ConstraintDescriptor> local_constraints;
    amrex::Gpu::DeviceVector<ConstraintTerm> local_terms;
    amrex::Gpu::DeviceVector<LocalPropertySupport> local_supports;
    amrex::Gpu::DeviceVector<ConstraintDescriptor> global_constraints;
    amrex::Gpu::DeviceVector<ConstraintTerm> global_terms;
    amrex::Gpu::DeviceVector<PopulationRemapView> populations;
    amrex::Gpu::DeviceVector<Real> edges;
    amrex::Gpu::DeviceVector<Real> pivots;
    amrex::Gpu::DeviceVector<int> property_components;
};

SBMTransport::~SBMTransport () = default;

SBMTransport::SBMTransport (const SBMLayout& layout,
                            const int number_of_levels,
                            const int max_groups_per_chunk)
    : m_layout(layout),
      m_projection(layout),
      m_groups(make_constraint_groups(layout)),
      m_flat_constraints(make_constraint_descriptors(layout)),
      m_property_support(make_attached_property_support_descriptors(layout)),
      m_chunks(make_constraint_closure_chunks(layout, max_groups_per_chunk)),
      m_max_constraints_per_group(m_groups.size(), 0),
      m_levels(static_cast<std::size_t>(number_of_levels)),
      m_device_metadata(std::make_unique<DeviceMetadata>())
{
    AMREX_ALWAYS_ASSERT(number_of_levels > 0);
    if (AMREX_SPACEDIM != 3) {
        throw std::invalid_argument(
            "SBM M3 mapped transport currently requires three dimensions");
    }

    std::vector<ConstraintDescriptor> local_descriptors;
    std::vector<ConstraintTerm> local_terms;
    std::vector<LocalPropertySupport> local_supports;
    m_device_metadata->host_groups.resize(m_groups.size());
    std::vector<int> group_chunk(m_groups.size(), -1);
    std::vector<int> group_local(m_groups.size(), -1);

    for (std::size_t chunk_index = 0; chunk_index < m_chunks.size();
         ++chunk_index) {
        const auto& chunk = m_chunks[chunk_index];
        int chunk_ratio_components = 0;
        std::vector<int> global_to_local(
            static_cast<std::size_t>(m_layout.ncomp()), -1);
        for (std::size_t local = 0; local < chunk.components.size(); ++local) {
            global_to_local[static_cast<std::size_t>(chunk.components[local])] =
                static_cast<int>(local);
        }
        for (std::size_t local_group = 0;
             local_group < chunk.group_indices.size(); ++local_group) {
            const int group_index = chunk.group_indices[local_group];
            const auto& group =
                m_groups[static_cast<std::size_t>(group_index)];
            group_chunk[static_cast<std::size_t>(group_index)] =
                static_cast<int>(chunk_index);
            group_local[static_cast<std::size_t>(group_index)] =
                static_cast<int>(local_group);
            auto& device_group =
                m_device_metadata
                    ->host_groups[static_cast<std::size_t>(group_index)];
            device_group.constraints_begin =
                static_cast<int>(local_descriptors.size());
            device_group.supports_begin =
                static_cast<int>(local_supports.size());
            device_group.chunk_index = static_cast<int>(chunk_index);
            device_group.local_group_index = static_cast<int>(local_group);

            for (const auto& descriptor : m_flat_constraints.constraints) {
                if (descriptor.group_index != group_index) {
                    continue;
                }
                ConstraintDescriptor local_descriptor = descriptor;
                local_descriptor.term_offset =
                    static_cast<int>(local_terms.size());
                for (int term_index = 0; term_index < descriptor.term_count;
                     ++term_index) {
                    auto term =
                        m_flat_constraints.terms[static_cast<std::size_t>(
                            descriptor.term_offset + term_index)];
                    const int local_component =
                        global_to_local[static_cast<std::size_t>(
                            term.component)];
                    if (local_component < 0) {
                        throw std::logic_error("SBM closure chunk split an "
                                               "atomic constraint group");
                    }
                    term.component = local_component;
                    local_terms.push_back(term);
                }
                local_descriptor.group_index = static_cast<int>(local_group);
                local_descriptors.push_back(local_descriptor);
                ++device_group.linear_constraint_count;
            }

            for (const auto& support : m_property_support) {
                if (support.group_index != group_index) {
                    continue;
                }
                LocalPropertySupport local_support;
                local_support.property_global = support.property_component;
                local_support.carrier_global = support.carrier_component;
                local_support.property_local =
                    global_to_local[static_cast<std::size_t>(
                        support.property_component)];
                local_support.carrier_local =
                    global_to_local[static_cast<std::size_t>(
                        support.carrier_component)];
                local_support.hard_min = support.hard_min;
                local_support.hard_max = support.hard_max;
                local_support.has_hard_min = support.has_hard_min;
                local_support.has_hard_max = support.has_hard_max;
                if (local_support.property_local < 0 ||
                    local_support.carrier_local < 0) {
                    throw std::logic_error(
                        "SBM closure chunk split an attached-property group");
                }
                local_supports.push_back(local_support);
                ++device_group.support_count;
            }
            // Keep the strict persisted upper edge open during transport by
            // reserving eta of the endpoint_low margin for interior 2M bins.
            if (group.moment_mode == MomentMode::TwoMoment) {
                const auto population = std::find_if(
                    m_layout.populations().begin(), m_layout.populations().end(),
                    [&group](const PopulationLayout& candidate) {
                        return candidate.population_id == group.population_id;
                    });
                if (population == m_layout.populations().end()) {
                    throw std::logic_error(
                        "SBM transport group has no matching population layout");
                }
                if (group.bin < population->grid.nbins() - 1) {
                    const int mass_global = population->mass_offset + group.bin;
                    const int number_global =
                        population->number_offset + group.bin;
                    device_group.canonical_mass_local =
                        global_to_local[static_cast<std::size_t>(mass_global)];
                    device_group.canonical_number_local =
                        global_to_local[static_cast<std::size_t>(number_global)];
                    if (device_group.canonical_mass_local < 0 ||
                        device_group.canonical_number_local < 0) {
                        throw std::logic_error(
                            "SBM closure chunk split an interior two-moment "
                            "canonical-reserve group");
                    }
                    const Real lower = population->grid.edges()[
                        static_cast<std::size_t>(group.bin)];
                    const Real upper = population->grid.edges()[
                        static_cast<std::size_t>(group.bin + 1)];
                    const Real width = upper - lower;
                    const Real eta =
                        Real(128.0) * std::numeric_limits<Real>::epsilon();
                    device_group.canonical_mass_coefficient =
                        -(Real(1.0) + eta) / width;
                    device_group.canonical_number_coefficient =
                        (Real(1.0) - eta) * upper / width;
                    device_group.canonical_reserve = 1;
                }
            }
            device_group.ratio_count = device_group.linear_constraint_count +
                                       2 * device_group.support_count +
                                       device_group.canonical_reserve;
            device_group.ratio_offset = chunk_ratio_components;
            AMREX_ALWAYS_ASSERT(device_group.ratio_count >= 0);
            AMREX_ALWAYS_ASSERT(
                device_group.ratio_count <=
                std::numeric_limits<int>::max() - chunk_ratio_components);
            chunk_ratio_components += device_group.ratio_count;
            m_max_constraints_per_group[static_cast<std::size_t>(group_index)] =
                device_group.ratio_count;
        }
    }

    auto& device = *m_device_metadata;
    device.local_constraints.resize(local_descriptors.size());
    device.local_terms.resize(local_terms.size());
    device.local_supports.resize(local_supports.size());
    if (!local_descriptors.empty()) {
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, local_descriptors.begin(),
                         local_descriptors.end(),
                         device.local_constraints.begin());
    }
    if (!local_terms.empty()) {
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, local_terms.begin(),
                         local_terms.end(), device.local_terms.begin());
    }
    if (!local_supports.empty()) {
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, local_supports.begin(),
                         local_supports.end(), device.local_supports.begin());
    }
    device.global_constraints.resize(m_flat_constraints.constraints.size());
    device.global_terms.resize(m_flat_constraints.terms.size());
    if (!m_flat_constraints.constraints.empty()) {
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         m_flat_constraints.constraints.begin(),
                         m_flat_constraints.constraints.end(),
                         device.global_constraints.begin());
    }
    if (!m_flat_constraints.terms.empty()) {
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_flat_constraints.terms.begin(),
            m_flat_constraints.terms.end(), device.global_terms.begin());
    }

    std::vector<Real> host_edges;
    std::vector<Real> host_pivots;
    std::vector<int> host_property_components;
    std::vector<int> edge_offsets, pivot_offsets, property_offsets;
    for (const auto& population : m_layout.populations()) {
        edge_offsets.push_back(static_cast<int>(host_edges.size()));
        host_edges.insert(host_edges.end(), population.grid.edges().begin(),
                          population.grid.edges().end());
        pivot_offsets.push_back(static_cast<int>(host_pivots.size()));
        host_pivots.insert(host_pivots.end(), population.grid.pivots().begin(),
                           population.grid.pivots().end());
        property_offsets.push_back(
            static_cast<int>(host_property_components.size()));
        host_property_components.insert(
            host_property_components.end(),
            population.property_component_offsets.begin(),
            population.property_component_offsets.end());
    }
    device.edges.resize(host_edges.size());
    device.pivots.resize(host_pivots.size());
    device.property_components.resize(host_property_components.size());
    if (!host_edges.empty()) {
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, host_edges.begin(),
                         host_edges.end(), device.edges.begin());
    }
    if (!host_pivots.empty()) {
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, host_pivots.begin(),
                         host_pivots.end(), device.pivots.begin());
    }
    if (!host_property_components.empty()) {
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, host_property_components.begin(),
            host_property_components.end(), device.property_components.begin());
    }
    std::vector<PopulationRemapView> device_views;
    device_views.reserve(m_layout.populations().size());
    for (std::size_t p = 0; p < m_layout.populations().size(); ++p) {
        const auto& population = m_layout.populations()[p];
        auto view = population_remap_view(m_layout, population.population_id);
        view.edges = device.edges.data() + edge_offsets[p];
        view.pivots = population.grid.pivots().empty()
                          ? nullptr
                          : device.pivots.data() + pivot_offsets[p];
        view.property_component_offsets =
            population.property_component_offsets.empty()
                ? nullptr
                : device.property_components.data() + property_offsets[p];
        device_views.push_back(view);
    }
    device.populations.resize(device_views.size());
    if (!device_views.empty()) {
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, device_views.begin(),
                         device_views.end(), device.populations.begin());
    }
}

void
SBMTransport::define (const int level,
                      const amrex::BoxArray& cell_ba,
                      const amrex::DistributionMapping& dm)
{
    AMREX_ALWAYS_ASSERT(level >= 0 &&
                        level < static_cast<int>(m_levels.size()));
    AMREX_ALWAYS_ASSERT(m_levels[static_cast<std::size_t>(level)] == nullptr);
    int max_chunk_components = 0;
    int max_chunk_ratio_components = 0;
    for (const auto& chunk : m_chunks) {
        max_chunk_components = std::max(
            max_chunk_components, static_cast<int>(chunk.components.size()));
        int chunk_ratio_components = 0;
        for (const int group : chunk.group_indices) {
            chunk_ratio_components +=
                m_max_constraints_per_group[static_cast<std::size_t>(group)];
        }
        max_chunk_ratio_components =
            std::max(max_chunk_ratio_components, chunk_ratio_components);
    }
    auto data = std::make_unique<LevelStorage>();
    const int ncomp = m_layout.ncomp();
    data->anchor.define(cell_ba, dm, ncomp, 0);
    data->target.define(cell_ba, dm, ncomp, 0);
    data->input_intensive.define(cell_ba, dm, ncomp, 2);
    data->physical_input_intensive.define(cell_ba, dm, ncomp, 2);
    data->base_intensive.define(cell_ba, dm, ncomp, 2);
    data->anchor_intensive.define(cell_ba, dm, ncomp, 2);
    data->low_trial_h.define(cell_ba, dm, max_chunk_components, 0);
    data->cell_ratios.define(cell_ba, dm,
                             std::max(1, max_chunk_ratio_components), 1);
    data->invalid.define(cell_ba, dm, 1, 0);
    data->outgoing_demand.define(cell_ba, dm, 1, 0);
    data->measure.define(cell_ba, dm, 1, 0);
    data->low_rate.define(cell_ba, dm, max_chunk_components, 0);
    data->high_rate.define(cell_ba, dm, max_chunk_components, 0);
    data->accepted_rate.define(cell_ba, dm, max_chunk_components, 0);
    data->face_lambda.define(cell_ba, dm, 1, 0);
    data->projected_rate.define(cell_ba, dm, 2, 0);
    data->projected_ledger.define(cell_ba, dm, 2);
    if (m_levels.size() > 1) {
        data->spectral_ledger.define(cell_ba, dm, ncomp);
    }
    data->anchor.setVal(Real(0.0));
    data->target.setVal(Real(0.0));
    data->invalid.setVal(Real(0.0));
    data->projected_rate.setVal(Real(0.0));
    m_levels[static_cast<std::size_t>(level)] = std::move(data);
}

bool
SBMTransport::rebuild_static_measure (const int level,
                                      const MultiFab& detJ,
                                      const MultiFab& mx,
                                      const MultiFab& my,
                                      std::string& diagnostic)
{
    diagnostic.clear();
    if (!is_defined(level)) {
        diagnostic =
            "static mapped measure requires a defined SBM transport level";
        return false;
    }
    auto& data = *m_levels[static_cast<std::size_t>(level)];
    data.measure_ready = false;
    if (!erf_auxiliary::BuildMappedCellMeasure(data.measure, detJ, mx, my,
                                               diagnostic))
        return false;
    data.measure_ready = true;
    return true;
}

void
SBMTransport::advance_stage_from_host (
    const int level, const erf_auxiliary::HostIntegrator method,
    const int stage, const double step_old_time, const double input_time,
    const double target_time, const double host_stage_interval,
    SBMStateManager& state_manager, const MultiFab& S_old_cons,
    const MultiFab& S_new_cons, MultiFab& S_data_cons,
    const MultiFab& avg_xmom, const MultiFab& avg_ymom,
    const MultiFab& avg_zmom, const amrex::Geometry& geometry,
    const int qc_component, const int qr_component)
{
    // Preserve ERF's semantic state roles at the integration boundary:
    // S_old anchors the step, S_new is the stage predictor, and S_data is the
    // target. In particular, do not substitute S_data for the predictor
    // density used to form the high-order intensive spectrum.
    advance_stage(level, method, stage, step_old_time, input_time, target_time,
                  host_stage_interval, state_manager, S_old_cons, S_new_cons,
                  S_data_cons, avg_xmom, avg_ymom, avg_zmom, geometry,
                  qc_component, qr_component);
}

void
SBMTransport::advance_stage (const int level,
                             const erf_auxiliary::HostIntegrator method,
                             const int stage,
                             const double step_old_time,
                             const double input_time,
                             const double target_time,
                             const double host_stage_interval,
                             SBMStateManager& state_manager,
                             const MultiFab& conserved_anchor,
                             const MultiFab& conserved_input,
                             MultiFab& conserved_target,
                             const MultiFab& avg_xmom,
                             const MultiFab& avg_ymom,
                             const MultiFab& avg_zmom,
                             const amrex::Geometry& geometry,
                             const int qc_component,
                             const int qr_component)
{
    AMREX_ALWAYS_ASSERT(is_defined(level));
    auto& data = *m_levels[static_cast<std::size_t>(level)];
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        data.measure_ready,
        "SBM M3 cannot advance before its static mapped measure is built");
    AMREX_ALWAYS_ASSERT(state_manager.is_defined(level));
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        AMREX_ALWAYS_ASSERT(geometry.isPeriodic(dir));
    }

    erf_auxiliary::AuxiliaryStageRecipe recipe;
    std::string diagnostic;
    if (!erf_auxiliary::MakeAuxiliaryStageRecipe(
            method, stage, host_stage_interval, recipe, diagnostic)) {
        amrex::Abort("SBM M3 stage recipe: " + diagnostic);
    }
    if (!std::isfinite(step_old_time) || !std::isfinite(input_time) ||
        !std::isfinite(target_time)) {
        amrex::Abort("SBM M3 received a nonfinite semantic time");
    }
    if (!state_manager.step_active(level)) {
        amrex::Abort("SBM M3 stage requires a state-manager begin-step transition at ERF::Advance");
    }
    const MultiFab& spectrum = stage == 0
                                   ? state_manager.old_state(level)
                                   : state_manager.new_state(level);
    const double spectrum_time = stage == 0
                                     ? state_manager.old_time(level)
                                     : state_manager.new_time(level);
    const double expected_spectrum_time = stage == 0 ? step_old_time : input_time;
    if (spectrum_time != expected_spectrum_time || input_time != spectrum_time) {
        std::ostringstream message;
        message.precision(17);
        message << "SBM M3 spectral input time mismatch at level " << level
                << " stage " << stage << ": expected " << expected_spectrum_time
                << " and requested input time " << input_time
                << ", got " << spectrum_time;
        amrex::Abort(message.str());
    }
    AMREX_ALWAYS_ASSERT(spectrum.nComp() == m_layout.ncomp());
    AMREX_ALWAYS_ASSERT(erf_auxiliary::SameCellLayout(spectrum, data.target));
    AMREX_ALWAYS_ASSERT(
        erf_auxiliary::SameCellLayout(spectrum, conserved_anchor));
    AMREX_ALWAYS_ASSERT(
        erf_auxiliary::SameCellLayout(spectrum, conserved_input));
    AMREX_ALWAYS_ASSERT(
        erf_auxiliary::SameCellLayout(spectrum, conserved_target));
    AMREX_ALWAYS_ASSERT(erf_auxiliary::SameCellLayout(spectrum, data.measure));
    AMREX_ALWAYS_ASSERT(conserved_anchor.nComp() > Rho_comp &&
                        conserved_input.nComp() > Rho_comp &&
                        conserved_target.nComp() > Rho_comp);
    if (stage == 0) {
        if (data.projected_ledger.step_active()) {
            amrex::Abort("SBM M3 received stage 0 before the previous timestep "
                         "completed");
        }
        MultiFab::Copy(data.anchor, spectrum, 0, 0, m_layout.ncomp(), 0);
    } else if (!data.projected_ledger.step_active() ||
               data.projected_ledger.next_stage() != stage) {
        amrex::Abort("SBM M3 stage arrived before stage 0 or out of order");
    }
    if (data.spectral_ledger.is_defined() &&
        !data.spectral_ledger.begin_stage(method, stage, step_old_time,
                                          recipe, diagnostic)) {
        amrex::Abort("SBM M4a spectral face-ledger stage sequence: " + diagnostic);
    }

    const MultiFab& trial_state =
        recipe.limiter_trial_base == erf_auxiliary::LimiterTrialBase::Anchor
            ? data.anchor
            : spectrum;
    const MultiFab& trial_conserved =
        recipe.limiter_trial_base == erf_auxiliary::LimiterTrialBase::Anchor
            ? conserved_anchor
            : conserved_input;
    if (!erf_auxiliary::ValidatePositiveFiniteComponent(conserved_target,
                                                        Rho_comp, diagnostic)) {
        amrex::Abort("SBM M3 target density: " + diagnostic);
    }

    build_intensive(data.input_intensive, spectrum, conserved_input,
                    geometry.periodicity(), "input");
    MultiFab::Copy(data.physical_input_intensive, data.input_intensive, 0, 0,
                   m_layout.ncomp(), 2);
    build_intensive(data.base_intensive, trial_state, trial_conserved,
                    geometry.periodicity(), "limiter-trial base");
    build_intensive(data.anchor_intensive, data.anchor, conserved_anchor,
                    geometry.periodicity(), "old-step anchor");

    // Transform every 2M pair in place to endpoint-number coordinates.  The
    // canonical transform is shared with persisted-state reconstruction and
    // is intentionally performed before native scalar WENO-Z3.
    data.invalid.setVal(Real(0.0));
    const auto populations = m_device_metadata->populations.data();
    const int population_count =
        static_cast<int>(m_layout.populations().size());
    for (amrex::MFIter mfi(data.input_intensive, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto input = data.input_intensive.array(mfi);
        const auto invalid = data.invalid.array(mfi);
        const auto endpoint_populations = populations;
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                int cell_invalid = 0;
                for (int p = 0; p < population_count; ++p) {
                    const auto& pop = endpoint_populations[p];
                    if (pop.moment_mode != MomentMode::TwoMoment) {
                        continue;
                    }
                    for (int bin = 0; bin < pop.nbins; ++bin) {
                        const int mass_comp = pop.mass_offset + bin;
                        const int number_comp = pop.number_offset + bin;
                        erf_sbm::EndpointTransform endpoint;
                        if (!erf_sbm::try_two_moment_to_endpoints(
                                input(i, j, k, number_comp),
                                input(i, j, k, mass_comp), pop.edges[bin],
                                pop.edges[bin + 1], endpoint)) {
                            cell_invalid = 1;
                        } else {
                            input(i, j, k, mass_comp) = endpoint.L;
                            input(i, j, k, number_comp) = endpoint.H;
                        }
                    }
                }
                invalid(i, j, k, 0) = static_cast<Real>(cell_invalid);
            });
    }
    if (data.invalid.max(0) != Real(0.0)) {
        amrex::Abort("SBM M3 input spectrum could not be represented in "
                     "endpoint-number coordinates");
    }
    data.input_intensive.FillBoundary(geometry.periodicity());

    // Check the donor outgoing-demand condition against the semantic limiter
    // base before constructing any low-order face rates.
    const auto inv_dx = geometry.InvCellSizeArray();
    const Real tau = static_cast<Real>(recipe.limiter_trial_interval);
    data.outgoing_demand.setVal(Real(0.0));
    for (amrex::MFIter mfi(data.outgoing_demand, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto fx = avg_xmom.const_array(mfi);
        const auto fy = avg_ymom.const_array(mfi);
        const auto fz = avg_zmom.const_array(mfi);
        const auto rho = trial_conserved.const_array(mfi);
        const auto omega = data.measure.const_array(mfi);
        const auto out = data.outgoing_demand.array(mfi);
        const auto dx_inv = inv_dx;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j,
                                                    int k) noexcept {
            out(i, j, k, 0) = tau *
                erf_auxiliary::MappedOutgoingDemandRate(
                    fx(i + 1, j, k, 0), fx(i, j, k, 0),
                    fy(i, j + 1, k, 0), fy(i, j, k, 0),
                    fz(i, j, k + 1, 0), fz(i, j, k, 0),
                    omega(i, j, k, 0), rho(i, j, k, Rho_comp), dx_inv);
        });
    }
    const Real max_outgoing = data.outgoing_demand.max(0);
    const Real cfl_tolerance =
        Real(64.0) * std::numeric_limits<Real>::epsilon();
    if (!std::isfinite(max_outgoing) ||
        max_outgoing > Real(1.0) + cfl_tolerance) {
        std::ostringstream message;
        message << "SBM M3 donor outgoing-demand condition failed: max="
                << max_outgoing << " tau=" << tau
                << " (no SBM subcycling is performed)";
        amrex::Abort(message.str());
    }

    data.projected_rate.setVal(Real(0.0));
    const auto& device = *m_device_metadata;
    for (std::size_t chunk_index = 0; chunk_index < m_chunks.size();
         ++chunk_index) {
        const auto& chunk = m_chunks[chunk_index];
        std::vector<int> global_to_local(
            static_cast<std::size_t>(m_layout.ncomp()), -1);
        for (std::size_t local = 0; local < chunk.components.size(); ++local) {
            global_to_local[static_cast<std::size_t>(chunk.components[local])] =
                static_cast<int>(local);
        }
        data.low_rate.setVal(Real(0.0));
        data.high_rate.setVal(Real(0.0));
        data.accepted_rate.setVal(Real(0.0));

        // Low-order donor rates use the recipe-selected trial base.
        for (amrex::MFIter mfi(data.input_intensive, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.tilebox();
            const auto z = data.base_intensive.const_array(mfi);
            const auto fd_x = avg_xmom.const_array(mfi);
            const auto fd_y = avg_ymom.const_array(mfi);
            const auto fd_z = avg_zmom.const_array(mfi);
            auto low_x = data.low_rate.dir(0).array(mfi);
            auto low_y = data.low_rate.dir(1).array(mfi);
            auto low_z = data.low_rate.dir(2).array(mfi);
            for (std::size_t local = 0; local < chunk.components.size();
                 ++local) {
                const int component = chunk.components[local];
                const int local_component = static_cast<int>(local);
                const amrex::Box xbx = amrex::surroundingNodes(bx, 0);
                const amrex::Box ybx = amrex::surroundingNodes(bx, 1);
                const amrex::Box zbx = amrex::surroundingNodes(bx, 2);
                amrex::ParallelFor(
                    xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const Real carrier = fd_x(i, j, k, 0);
                        const Real donor = carrier >= Real(0.0)
                                               ? z(i - 1, j, k, component)
                                               : z(i, j, k, component);
                        low_x(i, j, k, local_component) = carrier * donor;
                    });
                amrex::ParallelFor(
                    ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const Real carrier = fd_y(i, j, k, 0);
                        const Real donor = carrier >= Real(0.0)
                                               ? z(i, j - 1, k, component)
                                               : z(i, j, k, component);
                        low_y(i, j, k, local_component) = carrier * donor;
                    });
                amrex::ParallelFor(
                    zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const Real carrier = fd_z(i, j, k, 0);
                        const Real donor = carrier >= Real(0.0)
                                               ? z(i, j, k - 1, component)
                                               : z(i, j, k, component);
                        low_z(i, j, k, local_component) = carrier * donor;
                    });
            }
        }

        // ERF-native WENO-Z3 is applied to all input intensives.  The 2M
        // mass/number slots currently contain B_lo and B_hi, so the same
        // builder also supplies the required endpoint-coordinate candidate.
        for (amrex::MFIter mfi(data.input_intensive, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.tilebox();
            const GpuArray<const Array4<Real>, AMREX_SPACEDIM> flux_views{
                {data.high_rate.dir(0).array(mfi),
                 data.high_rate.dir(1).array(mfi),
                 data.high_rate.dir(2).array(mfi)}};
            for (std::size_t local = 0; local < chunk.components.size();
                 ++local) {
                BuildScalarAdvectionFluxes(
                    bx, data.input_intensive.const_array(mfi),
                    chunk.components[local], flux_views,
                    static_cast<int>(local), avg_xmom.const_array(mfi),
                    avg_ymom.const_array(mfi), avg_zmom.const_array(mfi),
                    AdvType::Weno_3Z, AdvType::Weno_3Z, Real(0.0), Real(0.0));
            }
        }

        // Invert endpoint face rates into physical water-mass and number
        // rates before the one common group limiter is applied.
        for (const auto& population : m_layout.populations()) {
            if (population.moment_mode != MomentMode::TwoMoment) {
                continue;
            }
            for (int bin = 0; bin < population.grid.nbins(); ++bin) {
                const int mass_global = population.mass_offset + bin;
                const int number_global = population.number_offset + bin;
                const int mass_local =
                    global_to_local[static_cast<std::size_t>(mass_global)];
                const int number_local =
                    global_to_local[static_cast<std::size_t>(number_global)];
                if (mass_local < 0 || number_local < 0) {
                    continue;
                }
                const Real lower =
                    population.grid.edges()[static_cast<std::size_t>(bin)];
                const Real upper =
                    population.grid.edges()[static_cast<std::size_t>(bin + 1)];
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    auto& face_rate = data.high_rate.dir(dir);
                    for (amrex::MFIter mfi(face_rate, amrex::TilingIfNotGPU());
                         mfi.isValid(); ++mfi) {
                        const amrex::Box bx = mfi.tilebox();
                        const auto flux = face_rate.array(mfi);
                        amrex::ParallelFor(
                            bx,
                            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                                const Real blo = flux(i, j, k, mass_local);
                                const Real bhi = flux(i, j, k, number_local);
                                flux(i, j, k, mass_local) =
                                    lower * blo + upper * bhi;
                                flux(i, j, k, number_local) = blo + bhi;
                            });
                    }
                }
            }
        }

        data.invalid.setVal(Real(0.0));
        data.cell_ratios.setVal(Real(0.0));
        for (const int group_index : chunk.group_indices) {
            const auto& group = m_groups[static_cast<std::size_t>(group_index)];
            const auto& group_device =
                device.host_groups[static_cast<std::size_t>(group_index)];
            AMREX_ALWAYS_ASSERT(group_device.ratio_offset >= 0);
            AMREX_ALWAYS_ASSERT(group_device.ratio_count >= 0);
            AMREX_ALWAYS_ASSERT(
                group_device.ratio_offset <= data.cell_ratios.nComp());
            AMREX_ALWAYS_ASSERT(
                group_device.ratio_count <=
                data.cell_ratios.nComp() - group_device.ratio_offset);

            // Form the complete donor low-order trial H^L for this atomic
            // group from the semantic limiter-trial base.
            for (const int component : group.members) {
                const int local_component =
                    global_to_local[static_cast<std::size_t>(component)];
                AMREX_ALWAYS_ASSERT(local_component >= 0);
                for (amrex::MFIter mfi(data.low_trial_h,
                                       amrex::TilingIfNotGPU());
                     mfi.isValid(); ++mfi) {
                    const amrex::Box bx = mfi.tilebox();
                    const auto u = trial_state.const_array(mfi);
                    const auto omega = data.measure.const_array(mfi);
                    const auto fx = data.low_rate.dir(0).const_array(mfi);
                    const auto fy = data.low_rate.dir(1).const_array(mfi);
                    const auto fz = data.low_rate.dir(2).const_array(mfi);
                    const auto hlow = data.low_trial_h.array(mfi);
                    const Real dx = inv_dx[0], dy = inv_dx[1], dz = inv_dx[2];
                    amrex::ParallelFor(
                        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            const Real div =
                                erf_auxiliary::ComputationalMappedDivergence(
                                    fx(i + 1, j, k, local_component),
                                    fx(i, j, k, local_component),
                                    fy(i, j + 1, k, local_component),
                                    fy(i, j, k, local_component),
                                    fz(i, j, k + 1, local_component),
                                    fz(i, j, k, local_component), dx, dy, dz);
                            hlow(i, j, k, local_component) =
                                omega(i, j, k, 0) * u(i, j, k, component) -
                                tau * div;
                        });
                }
            }

            const ConstraintDescriptor* local_constraints =
                device.local_constraints.data();
            const ConstraintTerm* local_terms = device.local_terms.data();
            const LocalPropertySupport* local_supports =
                device.local_supports.data();
            for (amrex::MFIter mfi(data.cell_ratios, amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi) {
                const amrex::Box bx = mfi.tilebox();
                const auto hlow = data.low_trial_h.const_array(mfi);
                const FaceArraySet low_faces{
                    {data.low_rate.dir(0).const_array(mfi),
                     data.low_rate.dir(1).const_array(mfi),
                     data.low_rate.dir(2).const_array(mfi)}};
                const FaceArraySet high_faces{
                    {data.high_rate.dir(0).const_array(mfi),
                     data.high_rate.dir(1).const_array(mfi),
                     data.high_rate.dir(2).const_array(mfi)}};
                const auto base = data.base_intensive.const_array(mfi);
                // Keep donor-envelope ratios in physical (M,N) coordinates.
                // The high-order input scratch is endpoint-transformed below.
                const auto input =
                    data.physical_input_intensive.const_array(mfi);
                const auto anchor = data.anchor_intensive.const_array(mfi);
                const auto ratios = data.cell_ratios.array(mfi);
                const auto invalid = data.invalid.array(mfi);
                const auto dx = inv_dx;
                const int constraints_begin = group_device.constraints_begin;
                const int linear_count = group_device.linear_constraint_count;
                const int supports_begin = group_device.supports_begin;
                const int support_count = group_device.support_count;
                const int ratio_offset = group_device.ratio_offset;
                const int canonical_reserve =
                    group_device.canonical_reserve;
                const int canonical_mass_local =
                    group_device.canonical_mass_local;
                const int canonical_number_local =
                    group_device.canonical_number_local;
                const Real canonical_mass_coefficient =
                    group_device.canonical_mass_coefficient;
                const Real canonical_number_coefficient =
                    group_device.canonical_number_coefficient;
                const Real roundoff =
                    Real(128.0) * std::numeric_limits<Real>::epsilon();
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j,
                                                            int k) noexcept {
                    int bad = 0;
                    for (int ci = 0; ci < linear_count; ++ci) {
                        const auto& descriptor =
                            local_constraints[constraints_begin + ci];
                        Real margin = Real(0.0);
                        Real scale = Real(0.0);
                        for (int t = 0; t < descriptor.term_count; ++t) {
                            const auto& term =
                                local_terms[descriptor.term_offset + t];
                            const Real value = hlow(i, j, k, term.component);
                            margin = std::fma(term.coefficient, value, margin);
                            scale += amrex::Math::abs(term.coefficient * value);
                        }
                        Real adverse = Real(0.0);
                        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                            const Real left_delta = face_correction(
                                descriptor, local_terms, low_faces, high_faces,
                                i, j, k, dir, -1, tau, dx);
                            const Real right_delta = face_correction(
                                descriptor, local_terms, low_faces, high_faces,
                                i, j, k, dir, 1, tau, dx);
                            adverse += amrex::max(Real(0.0), -left_delta) +
                                       amrex::max(Real(0.0), -right_delta);
                        }
                        if (!amrex::Math::isfinite(margin) ||
                            !amrex::Math::isfinite(adverse) ||
                            margin < -roundoff * scale) {
                            bad = 1;
                        }
                        ratios(i, j, k, ratio_offset + ci) =
                            adverse > Real(0.0)
                                ? clamp_unit(margin / adverse)
                                : Real(1.0);
                    }

                    for (int si = 0; si < support_count; ++si) {
                        const auto& support =
                            local_supports[supports_begin + si];
                        Real donor_min = std::numeric_limits<Real>::infinity();
                        Real donor_max = -std::numeric_limits<Real>::infinity();
                        int found = 0;
                        include_neighbourhood(base, i, j, k, support, donor_min,
                                              donor_max, found);
                        include_neighbourhood(input, i, j, k, support,
                                              donor_min, donor_max, found);
                        include_neighbourhood(anchor, i, j, k, support,
                                              donor_min, donor_max, found);
                        Real support_min = found != 0 ? donor_min : Real(0.0);
                        Real support_max = found != 0 ? donor_max : Real(0.0);
                        if (found != 0 && support.has_hard_min != 0) {
                            support_min =
                                amrex::max(support_min, support.hard_min);
                        }
                        if (found != 0 && support.has_hard_max != 0) {
                            support_max =
                                amrex::min(support_max, support.hard_max);
                        }
                        if (!amrex::Math::isfinite(support_min) ||
                            !amrex::Math::isfinite(support_max) ||
                            support_min > support_max) {
                            bad = 1;
                        }

                        const Real property =
                            hlow(i, j, k, support.property_local);
                        const Real carrier =
                            hlow(i, j, k, support.carrier_local);
                        for (int side = 0; side < 2; ++side) {
                            const bool upper = side != 0;
                            const Real slope =
                                upper ? support_max : support_min;
                            const Real margin =
                                upper ? slope * carrier - property
                                      : property - slope * carrier;
                            const Real scale =
                                amrex::Math::abs(property) +
                                amrex::Math::abs(slope * carrier);
                            Real adverse = Real(0.0);
                            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                                const Real left_delta =
                                    property_face_correction(
                                        support, low_faces, high_faces, i, j, k,
                                        dir, -1, tau, dx, slope, upper);
                                const Real right_delta =
                                    property_face_correction(
                                        support, low_faces, high_faces, i, j, k,
                                        dir, 1, tau, dx, slope, upper);
                                adverse += amrex::max(Real(0.0), -left_delta) +
                                           amrex::max(Real(0.0), -right_delta);
                            }
                            if (!amrex::Math::isfinite(margin) ||
                                !amrex::Math::isfinite(adverse) ||
                                margin < -roundoff * scale) {
                                bad = 1;
                            }
                            const int ratio_component =
                                ratio_offset + linear_count + 2 * si + side;
                            ratios(i, j, k, ratio_component) =
                                adverse > Real(0.0)
                                    ? clamp_unit(margin / adverse)
                                    : Real(1.0);
                        }
                    }

                    if (canonical_reserve != 0) {
                        // This transport-only reserve joins the existing
                        // closed FCT ratios; final canonical admission stays
                        // the authoritative persisted-state check.
                        const Real reserve = std::fma(
                            canonical_mass_coefficient,
                            hlow(i, j, k, canonical_mass_local),
                            canonical_number_coefficient *
                                hlow(i, j, k, canonical_number_local));
                        Real adverse = Real(0.0);
                        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                            const Real left_delta =
                                canonical_reserve_face_correction(
                                    canonical_mass_local,
                                    canonical_number_local,
                                    canonical_mass_coefficient,
                                    canonical_number_coefficient, low_faces,
                                    high_faces, i, j, k, dir, -1, tau, dx);
                            const Real right_delta =
                                canonical_reserve_face_correction(
                                    canonical_mass_local,
                                    canonical_number_local,
                                    canonical_mass_coefficient,
                                    canonical_number_coefficient, low_faces,
                                    high_faces, i, j, k, dir, 1, tau, dx);
                            adverse += amrex::max(Real(0.0), -left_delta) +
                                       amrex::max(Real(0.0), -right_delta);
                        }
                        const int ratio_component =
                            ratio_offset + linear_count + 2 * support_count;
                        if (!amrex::Math::isfinite(reserve) ||
                            !amrex::Math::isfinite(adverse)) {
                            bad = 1;
                            ratios(i, j, k, ratio_component) = Real(0.0);
                        } else if (reserve <= Real(0.0)) {
                            ratios(i, j, k, ratio_component) = Real(0.0);
                        } else if (adverse <= Real(0.0)) {
                            ratios(i, j, k, ratio_component) = Real(1.0);
                        } else {
                            ratios(i, j, k, ratio_component) =
                                clamp_unit(reserve / adverse);
                        }
                    }
                    if (bad != 0) {
                        invalid(i, j, k, 0) = Real(1.0);
                    }
                });
            }
        }

        // Every complete group has now written its own ratio components.
        // Exchange the chunk's ratio halo once before any group reads it.
        data.cell_ratios.FillBoundary(geometry.periodicity());
        if (data.invalid.max(0) != Real(0.0)) {
            amrex::Abort("SBM M3 donor low-order trial violates a linear or "
                         "attached-property/canonical-reserve constraint");
        }

        for (const int group_index : chunk.group_indices) {
            const auto& group = m_groups[static_cast<std::size_t>(group_index)];
            const auto& group_device =
                device.host_groups[static_cast<std::size_t>(group_index)];
            AMREX_ALWAYS_ASSERT(group_device.ratio_offset >= 0);
            AMREX_ALWAYS_ASSERT(group_device.ratio_count >= 0);
            AMREX_ALWAYS_ASSERT(
                group_device.ratio_offset <= data.cell_ratios.nComp());
            AMREX_ALWAYS_ASSERT(
                group_device.ratio_count <=
                data.cell_ratios.nComp() - group_device.ratio_offset);

            // Each face/group lambda is the minimum budget from both
            // adjacent cells and every constraint in the complete group.
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                auto& lambda_field = data.face_lambda.dir(dir);
                for (amrex::MFIter mfi(lambda_field, amrex::TilingIfNotGPU());
                     mfi.isValid(); ++mfi) {
                    const amrex::Box bx = mfi.tilebox();
                    const auto ratios = data.cell_ratios.const_array(mfi);
                    const auto lambda = lambda_field.array(mfi);
                    const int ratio_count = group_device.ratio_count;
                    const int ratio_offset = group_device.ratio_offset;
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(
                                               int i, int j, int k) noexcept {
                        int li = i, lj = j, lk = k;
                        if (dir == 0) {
                            --li;
                        } else if (dir == 1) {
                            --lj;
                        } else {
                            --lk;
                        }
                        Real accepted = Real(1.0);
                        for (int c = 0; c < ratio_count; ++c) {
                            accepted =
                                amrex::min(accepted,
                                           ratios(li, lj, lk, ratio_offset + c));
                            accepted = amrex::min(
                                accepted, ratios(i, j, k, ratio_offset + c));
                        }
                        lambda(i, j, k, 0) = clamp_unit(accepted);
                    });
                }
            }

            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                auto& accepted_field = data.accepted_rate.dir(dir);
                const auto& low_field = data.low_rate.dir(dir);
                const auto& high_field = data.high_rate.dir(dir);
                const auto& lambda_field = data.face_lambda.dir(dir);
                for (amrex::MFIter mfi(accepted_field, amrex::TilingIfNotGPU());
                     mfi.isValid(); ++mfi) {
                    const amrex::Box bx = mfi.tilebox();
                    const auto low = low_field.const_array(mfi);
                    const auto high = high_field.const_array(mfi);
                    const auto lambda = lambda_field.const_array(mfi);
                    const auto accepted = accepted_field.array(mfi);
                    for (const int component : group.members) {
                        const int local_component =
                            global_to_local[static_cast<std::size_t>(
                                component)];
                        AMREX_ALWAYS_ASSERT(local_component >= 0);
                        amrex::ParallelFor(
                            bx,
                            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                                const Real alpha = lambda(i, j, k, 0);
                                accepted(i, j, k, local_component) =
                                    low(i, j, k, local_component) +
                                    alpha * (high(i, j, k, local_component) -
                                             low(i, j, k, local_component));
                            });
                    }
                }
            }

            // Accumulate only the accepted spectral mass rate into the
            // two-component fixed liquid projection ledger.
            const auto& projection_spec = m_layout.liquid_projection();
            if (group.population_id == projection_spec.population_id) {
                const int mass_component =
                    m_layout.mass_offset(group.population_id) + group.bin;
                const int local_mass =
                    global_to_local[static_cast<std::size_t>(mass_component)];
                const int projected_component =
                    group.bin < projection_spec.cloud_rain_split ? 0 : 1;
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    MultiFab::Saxpy(data.projected_rate.dir(dir), Real(1.0),
                                    data.accepted_rate.dir(dir), local_mass,
                                    projected_component, 1, 0);
                }
            }
        }

        // Apply the ERF stage recurrence component by component while using
        // this bounded complete-group chunk's accepted face rates.
        for (std::size_t local = 0; local < chunk.components.size(); ++local) {
            const int component = chunk.components[local];
            erf_auxiliary::AuxiliaryStageContext context;
            context.method = method;
            context.level = level;
            context.stage = stage;
            context.step_old_time = step_old_time;
            context.input_time = input_time;
            context.target_time = target_time;
            context.recurrence = recipe;
            context.state_anchor = {&data.anchor, component, step_old_time};
            context.state_input = {&spectrum, component, input_time};
            context.state_target = {&data.target, component, target_time};
            context.rho_anchor = {&conserved_anchor, Rho_comp, step_old_time};
            context.rho_input = {&conserved_input, Rho_comp, input_time};
            context.rho_target = {&conserved_target, Rho_comp, target_time};
            context.measure_anchor = {&data.measure, 0, step_old_time};
            context.measure_input = {&data.measure, 0, input_time};
            context.measure_target = {&data.measure, 0, target_time};
            context.carrier = {&avg_xmom, &avg_ymom, &avg_zmom};
            erf_auxiliary::ApplyAuxiliaryMappedStage(
                context, data.accepted_rate, inv_dx, static_cast<int>(local));
        }

        if (data.spectral_ledger.is_defined()) {
            for (std::size_t local = 0; local < chunk.components.size(); ++local) {
                if (!data.spectral_ledger.accumulate_stage_component(
                        data.accepted_rate, static_cast<int>(local),
                        chunk.components[local], diagnostic)) {
                    amrex::Abort("SBM M4a spectral chunk-to-ledger mapping: " +
                                 diagnostic);
                }
            }
        }
    }

    // Candidate-wide device admission combines finite-value checking, every
    // declared linear group inequality, and strict canonical persisted-state
    // validation.  Diagnostic host scanning occurs only after this fast gate
    // fails.
    const ConstraintDescriptor* global_constraints =
        device.global_constraints.data();
    const ConstraintTerm* global_terms = device.global_terms.data();
    const PopulationRemapView* remap_populations = device.populations.data();
    const int global_constraint_count =
        static_cast<int>(m_flat_constraints.constraints.size());
    const int remap_population_count =
        static_cast<int>(m_layout.populations().size());
    data.invalid.setVal(Real(0.0));
    for (amrex::MFIter mfi(data.target, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto candidate = data.target.const_array(mfi);
        const auto invalid = data.invalid.array(mfi);
        const int state_components = m_layout.ncomp();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j,
                                                    int k) noexcept {
            int bad = 0;
            for (int component = 0; component < state_components; ++component) {
                if (!amrex::Math::isfinite(candidate(i, j, k, component))) {
                    bad = 1;
                }
            }
            const Real roundoff =
                Real(128.0) * std::numeric_limits<Real>::epsilon();
            for (int ci = 0; ci < global_constraint_count; ++ci) {
                const auto& descriptor = global_constraints[ci];
                Real margin = Real(0.0);
                Real scale = Real(0.0);
                for (int t = 0; t < descriptor.term_count; ++t) {
                    const auto& term = global_terms[descriptor.term_offset + t];
                    const Real value = candidate(i, j, k, term.component);
                    margin = std::fma(term.coefficient, value, margin);
                    scale += amrex::Math::abs(term.coefficient * value);
                }
                if (!amrex::Math::isfinite(margin) ||
                    margin < -roundoff * scale)
                    bad = 1;
            }
            const erf_sbm::remap_detail::RemapCellStateView cell_state{
                candidate, i, j, k};
            for (int p = 0; p < remap_population_count; ++p) {
                const auto& population = remap_populations[p];
                for (int bin = 0; bin < population.nbins; ++bin) {
                    if (!erf_sbm::remap_detail::canonical_persisted_bin_state(
                            population, bin, cell_state, state_components)) {
                        bad = 1;
                    }
                }
            }
            invalid(i, j, k, 0) = static_cast<Real>(bad);
        });
    }
    if (data.invalid.max(0) != Real(0.0)) {
        std::string reason;
        const bool host_admissible = erf_sbm::authoritative_state_admissible(
            data.target, m_layout, level, &reason);
        amrex::ignore_unused(host_admissible);
        if (reason.empty()) {
            reason = "candidate failed device linear/canonical admission";
        }
        amrex::Abort("SBM M3 candidate rejected before commit: " + reason);
    }

    // Only fully admitted candidates may finish an authoritative spectral
    // ledger stage. Chunk calls above accumulate components without advancing
    // the host stage sequence.
    // Ledger admission is MPI-collective: every rank on this level must enter
    // despite local recoverable errors; do not branch around these calls using
    // rank-local predicates.
    if (data.spectral_ledger.is_defined() &&
        !data.spectral_ledger.finish_stage(diagnostic)) {
        amrex::Abort("SBM M4a spectral face-ledger stage finish: " + diagnostic);
    }
    // Retain the existing two-component projected ledger and its exact host
    // stage recipe independently of the full spectral integral.
    if (!data.projected_ledger.accept_stage(method, stage, step_old_time,
                                            recipe, data.projected_rate,
                                            diagnostic)) {
        amrex::Abort("SBM M3 projected face-ledger stage sequence: " +
                     diagnostic);
    }
    MultiFab::Copy(state_manager.new_target_storage(level), data.target, 0, 0,
                   m_layout.ncomp(), 0);
    const bool physical_step_complete =
        (method == erf_auxiliary::HostIntegrator::CompressibleRK3 && stage == 2) ||
        (method == erf_auxiliary::HostIntegrator::AnelasticHeun && stage == 1);
    if (!state_manager.accept_stage_target(level, target_time,
                                          physical_step_complete, diagnostic)) {
        amrex::Abort("SBM M4a accepted-state lifecycle commit: " + diagnostic);
    }
    state_manager.project_to_core(level, conserved_target, qc_component,
                                  qr_component);
}

void
SBMTransport::destroy (const int level)
{
    AMREX_ALWAYS_ASSERT(level >= 0 &&
                        level < static_cast<int>(m_levels.size()));
    m_levels[static_cast<std::size_t>(level)].reset();
}

bool
SBMTransport::is_defined (const int level) const
{
    if (level < 0 || level >= static_cast<int>(m_levels.size())) {
        throw std::out_of_range("SBM transport level is outside the workspace");
    }
    return static_cast<bool>(m_levels[static_cast<std::size_t>(level)]);
}

bool
SBMTransport::measure_is_ready (const int level) const
{
    if (!is_defined(level)) {
        return false;
    }
    return m_levels[static_cast<std::size_t>(level)]->measure_ready;
}

const MultiFab&
SBMTransport::static_measure (const int level) const
{
    if (!is_defined(level)) {
        throw std::logic_error("SBM transport level is not defined");
    }
    const auto& data = *m_levels[static_cast<std::size_t>(level)];
    if (!data.measure_ready) {
        throw std::logic_error("SBM static measure is not ready");
    }
    return data.measure;
}

const erf_auxiliary::CompletedStepFluxLedger&
SBMTransport::projected_ledger (const int level) const
{
    if (!is_defined(level)) {
        throw std::logic_error("SBM transport level is not defined");
    }
    return m_levels[static_cast<std::size_t>(level)]->projected_ledger;
}

const erf_auxiliary::CompletedStepFluxLedger&
SBMTransport::completed_spectral_ledger (const int level) const
{
    if (!is_defined(level)) {
        throw std::logic_error("SBM transport level is not defined");
    }
    const auto& ledger = m_levels[static_cast<std::size_t>(level)]->spectral_ledger;
    if (!ledger.is_defined()) {
        throw std::logic_error(
            "completed spectral transfer is available only on AMR-capable transport instances");
    }
    if (!ledger.step_complete()) {
        throw std::logic_error(
            "spectral transfer is available only after a complete host step");
    }
    return ledger;
}

} // namespace erf_sbm
