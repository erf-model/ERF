#include "ERF_SBMTransport.H"
#include "ERF_SBMAMRFluxRegister.H"
#ifdef ERF_SBM_AMR_QUALIFICATION
#include "ERF_SBMAMRQualification.H"
#endif

#include "ERF_SBMRemapping.H"
#include "ERF_SBMRestart.H"
#include "ERF_SBMAMRTime.H"
#include "AuxiliaryState/ERF_AuxiliaryMappedTransport.H"
#include "ERF_IndexDefines.H"
#include "Advection/ERF_AdvectionSrcForScalars.H"

#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
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

bool
coarse_stage_theta (const double requested_time,
                    const double coarse_old_time,
                    const double coarse_new_time,
                    const double dt_min,
                    double& theta,
                    std::string& diagnostic)
{
    diagnostic.clear();
    theta = 0.0;
    if (!std::isfinite(requested_time) ||
        !std::isfinite(coarse_old_time) ||
        !std::isfinite(coarse_new_time) ||
        !std::isfinite(dt_min) || !(dt_min > 0.0) ||
        !(coarse_new_time > coarse_old_time)) {
        diagnostic = "coarse SBM stage bracket has invalid times or timestep";
        return false;
    }
    if (semantic_time_matches(requested_time, coarse_old_time, dt_min)) {
        theta = 0.0;
        return true;
    }
    if (semantic_time_matches(requested_time, coarse_new_time, dt_min)) {
        theta = 1.0;
        return true;
    }
    if (!(requested_time > coarse_old_time &&
          requested_time < coarse_new_time)) {
        diagnostic = "fine SBM stage time is outside the accepted coarse bracket";
        return false;
    }
    theta = (requested_time - coarse_old_time) /
            (coarse_new_time - coarse_old_time);
    if (!std::isfinite(theta) || !(theta > 0.0 && theta < 1.0)) {
        diagnostic = "fine SBM stage interpolation fraction is invalid";
        return false;
    }
    return true;
}

// One collective: every staged endpoint value (valid and ghost) is finite and
// both staged densities are positive.
bool
coarse_stage_endpoints_valid (const MultiFab& coarse_stage,
                              const int component_count,
                              const int new_offset,
                              const int rho_old_local)
{
    amrex::ReduceOps<amrex::ReduceOpMin> reduce_op;
    amrex::ReduceData<int> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (amrex::MFIter mfi(coarse_stage, false); mfi.isValid(); ++mfi) {
        const auto values = coarse_stage.const_array(mfi);
        reduce_op.eval(mfi.fabbox(), reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple {
                const Real rho_old = values(i, j, k, rho_old_local);
                const Real rho_new = values(i, j, k, rho_old_local + 1);
                bool ok = amrex::Math::isfinite(rho_old) &&
                          amrex::Math::isfinite(rho_new) &&
                          rho_old > Real(0.0) && rho_new > Real(0.0);
                for (int n = 0; n < component_count && ok; ++n) {
                    ok = amrex::Math::isfinite(values(i, j, k, n)) &&
                         amrex::Math::isfinite(values(i, j, k, new_offset + n));
                }
                return {ok ? 1 : 0};
            });
    }
    int all_valid = amrex::get<0>(reduce_data.value(reduce_op));
    amrex::ParallelDescriptor::ReduceIntMin(all_valid);
    return all_valid != 0;
}

void
fill_coarse_parent_ghosts (MultiFab& destination,
                           MultiFab& missing_fine_ghosts,
                           const MultiFab& coarse_old_spectrum,
                           const MultiFab& coarse_new_spectrum,
                           const MultiFab& fine_source_state,
                           const MultiFab& fine_source_conserved,
                           const SBMAMRStageParent& parent,
                           const amrex::Geometry& fine_geometry,
                           const std::vector<ConstraintClosureChunk>& chunks,
                           const int max_chunk_components,
                           const double requested_time,
                           const double coarse_old_spectrum_time,
                           const double coarse_new_spectrum_time,
                           const int level,
                           const int stage,
                           const char* role,
                           MultiFab& coarse_pack,
                           MultiFab& coarse_stage,
                           const SBMAMRStageGhostProbe* probe = nullptr,
                           std::array<bool, 2>* split_values_checked = nullptr)
{
    using ProfileClock = std::chrono::steady_clock;
    const bool profile_parent_fill = parent.profile_parent_fill;
    const auto fill_start = profile_parent_fill ? ProfileClock::now()
                                                : ProfileClock::time_point{};
    double parallel_copy_seconds = 0.0;
    int parallel_copy_calls = 0;
    const int current_role = std::strcmp(role, "input") == 0 ? 0
        : (std::strcmp(role, "limiter-base") == 0 ? 1 : 2);
    const bool probe_this_role = probe == nullptr ||
        probe->requested_role == current_role;
    bool inject_rank_local_store_fault = false;
#ifdef ERF_SBM_AMR_QUALIFICATION
    // Rank 0 skips component 0's actual store at its first coarse-only ghost
    // while keeping the independent wrote_components tally.
    inject_rank_local_store_fault =
        parent.inject_parent_store_fault && level == 1 && stage == 0 &&
        current_role == 0 && amrex::ParallelDescriptor::MyProc() == 0;
#endif
    const double coarse_dt = parent.new_time - parent.old_time;
    const double dt_min = std::min(parent.fine_step_dt, coarse_dt);
    double theta = 0.0;
    std::string diagnostic;
    if (!coarse_stage_theta(requested_time, parent.old_time, parent.new_time,
                            dt_min, theta, diagnostic)) {
        amrex::Abort("SBM coarse/fine stage interpolation: " + diagnostic);
    }
    if (!semantic_time_matches(coarse_old_spectrum_time, parent.old_time,
                               dt_min) ||
        !semantic_time_matches(coarse_new_spectrum_time, parent.new_time,
                               dt_min)) {
        amrex::Abort("SBM coarse spectrum and host density endpoints do not "
                     "share semantic times");
    }
    if (parent.conserved_old == nullptr || parent.conserved_new == nullptr ||
        !parent.refinement_ratio.allGT(0) ||
        coarse_old_spectrum.nComp() != destination.nComp() ||
        coarse_new_spectrum.nComp() != destination.nComp() ||
        parent.conserved_old->nComp() <= Rho_comp ||
        parent.conserved_new->nComp() <= Rho_comp ||
        !erf_auxiliary::SameCellLayout(coarse_old_spectrum,
                                       coarse_new_spectrum) ||
        !erf_auxiliary::SameCellLayout(coarse_old_spectrum,
                                       *parent.conserved_old) ||
        !erf_auxiliary::SameCellLayout(coarse_old_spectrum,
                                       *parent.conserved_new) ||
        !erf_auxiliary::SameCellLayout(destination, missing_fine_ghosts) ||
        missing_fine_ghosts.nComp() < 3 ||
        !destination.nGrowVect().allGE(amrex::IntVect(2)) ||
        !missing_fine_ghosts.nGrowVect().allGE(amrex::IntVect(2))) {
        amrex::Abort("SBM coarse/fine stage interpolation received an "
                     "incompatible layout");
    }

    amrex::BoxArray coarse_on_fine_ba = destination.boxArray();
    coarse_on_fine_ba.coarsen(parent.refinement_ratio);
    const amrex::IntVect coarse_ghosts(2);
    if (coarse_stage.nComp() != 2 * max_chunk_components + 2 ||
        coarse_stage.boxArray() != coarse_on_fine_ba ||
        coarse_stage.DistributionMap() != destination.DistributionMap() ||
        coarse_stage.nGrowVect() != coarse_ghosts) {
        coarse_stage.define(coarse_on_fine_ba, destination.DistributionMap(),
                            2 * max_chunk_components + 2, coarse_ghosts);
    }
    coarse_stage.setVal(std::numeric_limits<Real>::quiet_NaN());
    const bool componentwise_parent_fill = parent.componentwise_parent_fill;
    if (!componentwise_parent_fill) {
        if (coarse_pack.nComp() != 2 * max_chunk_components + 2 ||
            coarse_pack.boxArray() != coarse_old_spectrum.boxArray() ||
            coarse_pack.DistributionMap() !=
                coarse_old_spectrum.DistributionMap()) {
            coarse_pack.define(coarse_old_spectrum.boxArray(),
                               coarse_old_spectrum.DistributionMap(),
                               2 * max_chunk_components + 2, 0);
            coarse_pack.setVal(std::numeric_limits<Real>::quiet_NaN());
        }
    }

    amrex::Box coarse_domain = fine_geometry.Domain();
    coarse_domain.coarsen(parent.refinement_ratio);
    const amrex::Periodicity coarse_periodicity(coarse_domain.length());
    const int rho_old_local = 2 * max_chunk_components;
    const int rho_new_local = rho_old_local + 1;
    auto timed_parallel_copy = [&] (auto&& operation) {
        if (profile_parent_fill) {
            const auto copy_start = ProfileClock::now();
            operation();
            parallel_copy_seconds +=
                std::chrono::duration<double>(ProfileClock::now() - copy_start).count();
        } else {
            operation();
        }
        ++parallel_copy_calls;
    };

    // Component 0 marks coarse-only ghosts, component 1 counts component
    // writes in this role, and component 2 marks an unrepresentable quotient.
    // Seed the ghost region, clear valid cells, then exchange the valid mask:
    // same-level and periodic donor ghosts become zero while coarse-only
    // ghosts stay marked.  Fine valid data and fine-fine/periodic donors are
    // never poisoned or replaced by parent values.
    missing_fine_ghosts.setVal(Real(0.0), 0, 3, 2);
    missing_fine_ghosts.setVal(Real(1.0), 0, 1, 2);
    missing_fine_ghosts.setVal(Real(0.0), 0, 1, 0);
    missing_fine_ghosts.FillBoundary(0, 1, fine_geometry.periodicity());
    amrex::IntVect fault_cell(std::numeric_limits<int>::min());
#ifdef ERF_SBM_AMR_QUALIFICATION
    if (inject_rank_local_store_fault) {
        fault_cell = qualification::first_coarse_only_ghost(missing_fine_ghosts);
    }
#endif

    if (probe != nullptr && probe_this_role &&
        probe->seed_finite_stale_value) {
        const auto probe_cell = probe->fine_cell;
        const int probe_component = probe->global_component;
        const Real sentinel = probe->stale_value;
        for (amrex::MFIter mfi(destination, false); mfi.isValid(); ++mfi) {
            const auto target = destination.array(mfi);
            const auto missing = missing_fine_ghosts.const_array(mfi);
            amrex::ParallelFor(
                mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    if (missing(i, j, k, 0) == Real(1.0) &&
                        i == probe_cell[0] && j == probe_cell[1] &&
                        k == probe_cell[2]) {
                        target(i, j, k, probe_component) = sentinel;
                    }
                });
        }
    }

    // Poison only the cells classified as coarse-only ghosts in this role.
    // This makes a missed destination store observable even when that FAB
    // held a finite value from an earlier role or stage.  The destination is
    // the only N-component storage touched; no all-component staging copy is
    // created.
    const Real poison = std::numeric_limits<Real>::quiet_NaN();
    for (amrex::MFIter mfi(destination, false); mfi.isValid(); ++mfi) {
        const amrex::Box fab_box = mfi.fabbox();
        const auto missing = missing_fine_ghosts.const_array(mfi);
        const auto target = destination.array(mfi);
        const int ncomp = destination.nComp();
        amrex::ParallelFor(
            fab_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                if (missing(i, j, k, 0) != Real(1.0)) return;
                for (int n = 0; n < ncomp; ++n) {
                    target(i, j, k, n) = poison;
                }
            });
    }

    std::vector<int> flattened_components;
    std::vector<int> chunk_offsets;
    flattened_components.reserve(static_cast<std::size_t>(destination.nComp()));
    chunk_offsets.reserve(chunks.size());
    for (const auto& chunk : chunks) {
        chunk_offsets.push_back(static_cast<int>(flattened_components.size()));
        flattened_components.insert(flattened_components.end(),
                                    chunk.components.begin(),
                                    chunk.components.end());
    }
    if (flattened_components.size() !=
        static_cast<std::size_t>(destination.nComp())) {
        amrex::Abort("SBM coarse/fine stage chunks do not cover every "
                     "spectrum component exactly once");
    }
    amrex::Gpu::DeviceVector<int> device_components(flattened_components.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, flattened_components.begin(),
                     flattened_components.end(), device_components.begin());

    const auto coarse_old_periodic = coarse_periodicity;
    for (std::size_t chunk_index = 0; chunk_index < chunks.size();
         ++chunk_index) {
        const auto& chunk = chunks[chunk_index];
        if (chunk.components.empty() ||
            static_cast<int>(chunk.components.size()) >
                max_chunk_components) {
            amrex::Abort("SBM coarse/fine stage chunk exceeds bounded scratch");
        }
        const int component_count =
            static_cast<int>(chunk.components.size());
        if (componentwise_parent_fill) {
            for (std::size_t local = 0; local < chunk.components.size(); ++local) {
                const int global = chunk.components[local];
                timed_parallel_copy([&] {
                    coarse_stage.ParallelCopy(
                        coarse_old_spectrum, global, static_cast<int>(local), 1,
                        amrex::IntVect(0), coarse_ghosts, coarse_old_periodic);
                });
                timed_parallel_copy([&] {
                    coarse_stage.ParallelCopy(
                        coarse_new_spectrum, global,
                        max_chunk_components + static_cast<int>(local), 1,
                        amrex::IntVect(0), coarse_ghosts, coarse_old_periodic);
                });
            }
            timed_parallel_copy([&] {
                coarse_stage.ParallelCopy(
                    *parent.conserved_old, Rho_comp, rho_old_local, 1,
                    amrex::IntVect(0), coarse_ghosts, coarse_old_periodic);
            });
            timed_parallel_copy([&] {
                coarse_stage.ParallelCopy(
                    *parent.conserved_new, Rho_comp, rho_new_local, 1,
                    amrex::IntVect(0), coarse_ghosts, coarse_old_periodic);
            });
        } else {
            const int component_offset = chunk_offsets[chunk_index];
            const int* component_map_local =
                device_components.data() + component_offset;
            for (amrex::MFIter mfi(coarse_pack, false); mfi.isValid(); ++mfi) {
                const amrex::Box bx = mfi.validbox();
                const auto old_values = coarse_old_spectrum.const_array(mfi);
                const auto new_values = coarse_new_spectrum.const_array(mfi);
                const auto rho_old = parent.conserved_old->const_array(mfi);
                const auto rho_new = parent.conserved_new->const_array(mfi);
                const auto pack = coarse_pack.array(mfi);
                amrex::ParallelFor(
                    bx, 2 * component_count + 2,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k,
                                         int local) noexcept {
                        if (local < component_count) {
                            const int global = component_map_local[local];
                            pack(i, j, k, local) =
                                old_values(i, j, k, global);
                        } else if (local < 2 * component_count) {
                            const int local_component =
                                local - component_count;
                            const int global =
                                component_map_local[local_component];
                            pack(i, j, k,
                                 max_chunk_components + local_component) =
                                new_values(i, j, k, global);
                        } else if (local == 2 * component_count) {
                            pack(i, j, k, rho_old_local) =
                                rho_old(i, j, k, Rho_comp);
                        } else {
                            pack(i, j, k, rho_new_local) =
                                rho_new(i, j, k, Rho_comp);
                        }
                    });
            }
            timed_parallel_copy([&] {
                coarse_stage.ParallelCopy(
                    coarse_pack, 0, 0, 2 * max_chunk_components + 2,
                    amrex::IntVect(0), coarse_ghosts, coarse_old_periodic);
            });
        }

        if (!coarse_stage_endpoints_valid(coarse_stage, component_count,
                                          max_chunk_components,
                                          rho_old_local)) {
            std::ostringstream message;
            message.precision(17);
            message << "SBM coarse/fine endpoint data are not finite with "
                    << "positive density (level=" << level << ", stage="
                    << stage << ", role=" << role << ", requested_time="
                    << requested_time << ", coarse_old=" << parent.old_time
                    << ", coarse_new=" << parent.new_time << ", rank="
                    << amrex::ParallelDescriptor::MyProc() << ')';
            amrex::Abort(message.str());
        }

        const int component_offset = chunk_offsets[chunk_index];
        const Real theta_real = static_cast<Real>(theta);
        const auto component_map = device_components.data();
        const amrex::IntVect refinement_ratio = parent.refinement_ratio;
        const amrex::IntVect probe_cell =
            probe == nullptr ? amrex::IntVect(-1) : probe->fine_cell;
        const int skipped_component =
            probe == nullptr ? -1 : probe->skipped_global_component;
        const bool skip_tally_only =
            probe != nullptr && probe->skip_write_tally_only;
        const int probe_stage = probe == nullptr ? -1 : probe->requested_stage;
        const bool probe_role_matches = probe_this_role;
        for (amrex::MFIter mfi(destination, false); mfi.isValid(); ++mfi) {
            const amrex::Box fab_box = mfi.fabbox();
            const auto coarse = coarse_stage.const_array(mfi);
            const auto missing = missing_fine_ghosts.array(mfi);
            const auto target = destination.array(mfi);
            const auto component_map_local = component_map + component_offset;
            amrex::ParallelFor(
                fab_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    if (missing(i, j, k, 0) != Real(1.0)) return;
                    amrex::IntVect parent_cell(i, j, k);
                    parent_cell.coarsen(refinement_ratio);
                    const Real rho0 = coarse(parent_cell[0], parent_cell[1],
                                             parent_cell[2], rho_old_local);
                    const Real rho1 = coarse(parent_cell[0], parent_cell[1],
                                             parent_cell[2], rho_new_local);
                    const Real rho_t = (Real(1.0) - theta_real) * rho0 +
                                       theta_real * rho1;
                    int wrote_components = 0;
                    for (int local = 0; local < component_count; ++local) {
                        const int global = component_map_local[local];
                        const bool probe_fault_cell =
                            global == skipped_component &&
                            probe_role_matches &&
                            (probe_stage < 0 || probe_stage == stage) &&
                            i == probe_cell[0] && j == probe_cell[1] &&
                            k == probe_cell[2];
                        const bool skip_assignment =
                            probe_fault_cell && !skip_tally_only;
                        const bool skip_rank_local_assignment =
                            inject_rank_local_store_fault && global == 0 &&
                            i == fault_cell[0] && j == fault_cell[1] &&
                            k == fault_cell[2];
                        const Real u0 = coarse(parent_cell[0], parent_cell[1],
                                               parent_cell[2], local);
                        const Real u1 = coarse(
                            parent_cell[0], parent_cell[1], parent_cell[2],
                            max_chunk_components + local);
                        const Real u_t = (Real(1.0) - theta_real) * u0 +
                                         theta_real * u1;
                        const Real intensive = u_t / rho_t;
                        const bool positive_support_expected =
                            ((theta_real < Real(1.0)) && u0 > Real(0.0)) ||
                            ((theta_real > Real(0.0)) && u1 > Real(0.0));
                        const bool quotient_ok =
                            amrex::Math::isfinite(rho_t) && rho_t > Real(0.0) &&
                            amrex::Math::isfinite(u_t) && u_t >= Real(0.0) &&
                            amrex::Math::isfinite(intensive) &&
                            intensive >= Real(0.0) &&
                            (!(positive_support_expected) ||
                             (u_t > Real(0.0) && intensive > Real(0.0)));
                        if (!skip_assignment && !skip_rank_local_assignment) {
                            if (quotient_ok) {
                                target(i, j, k, global) = intensive;
                            } else {
                                // Keep a generated quotient failure distinct
                                // from a missing destination store.
                                target(i, j, k, global) =
                                    std::numeric_limits<Real>::infinity();
                                missing(i, j, k, 2) = Real(1.0);
                            }
                        }
                        if (!(probe_fault_cell && skip_tally_only)) {
                            ++wrote_components;
                        }
                    }
                    missing(i, j, k, 1) += static_cast<Real>(wrote_components);
                });
        }
    }

    // Independently inspect actual destination values at all marked coarse-
    // only ghosts before any stencil can consume them.  The per-cell loop
    // avoids races when multiple components fail at one location; the mask
    // carries one scalar failure code for the collective reduction.
    for (amrex::MFIter mfi(destination, false); mfi.isValid(); ++mfi) {
        const amrex::Box fab_box = mfi.fabbox();
        const auto mask = missing_fine_ghosts.array(mfi);
        const int ncomp = destination.nComp();
        const auto target = destination.const_array(mfi);
        const int sample_component = probe == nullptr
            ? -1 : probe->global_component;
        const amrex::IntVect sample_cell = probe == nullptr
            ? amrex::IntVect(-1) : probe->fine_cell;
        const bool sample_role_matches = probe_this_role;
        const int sample_stage = probe == nullptr
            ? -1 : probe->requested_stage;
        auto* sampled_value = probe == nullptr
            ? nullptr : probe->sampled_intensive_device;
        auto* sampled_hit = probe == nullptr ? nullptr : probe->hit_device;
        amrex::ParallelFor(
            fab_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                if (mask(i, j, k, 0) != Real(1.0)) return;
                int destination_invalid = 0;
                const bool quotient_failed =
                    mask(i, j, k, 2) == Real(1.0);
                for (int component = 0; component < ncomp; ++component) {
                    const Real value = target(i, j, k, component);
                    if (!amrex::Math::isfinite(value) || value < Real(0.0)) {
                        // The quotient path writes +infinity deliberately as
                        // its invalid marker. Keep that distinct from the
                        // NaN poison that proves a destination store was
                        // skipped; any other bad component remains visible.
                        if (quotient_failed &&
                            value == std::numeric_limits<Real>::infinity()) {
                            continue;
                        }
                        destination_invalid = 1;
                        break;
                    }
                }
                if (destination_invalid != 0) {
                    // Encode both independent failure classes if necessary:
                    // 1 = quotient, 2 = invalid destination, 3 = both.
                    mask(i, j, k, 2) = quotient_failed
                        ? Real(3.0) : Real(2.0);
                }
                if (sampled_value != nullptr && sampled_hit != nullptr &&
                    sample_role_matches &&
                    (sample_stage < 0 || sample_stage == stage) &&
                    i == sample_cell[0] && j == sample_cell[1] &&
                    k == sample_cell[2]) {
                    *sampled_value = target(i, j, k, sample_component);
                    *sampled_hit = 1;
                }
            });
    }
    const Real destination_failure = missing_fine_ghosts.max(2, 2);
    if (destination_failure == Real(2.0) ||
        destination_failure == Real(3.0)) {
        std::ostringstream message;
        message.precision(17);
        message << "SBM coarse/fine stage left a coarse-only parent ghost "
                << "component unwritten or invalid in the destination"
                << (destination_failure == Real(3.0)
                        ? " (also found an invalid quotient)" : "")
                << " (level=" << level << ", stage=" << stage << ", role="
                << role << ", requested_time=" << requested_time
                << ", coarse_old=" << parent.old_time << ", coarse_new="
                << parent.new_time << ", rank="
                << amrex::ParallelDescriptor::MyProc()
                << (inject_rank_local_store_fault
                        ? ", assignment_only_fault=rank-local" : "")
                << ')';
        amrex::Abort(message.str());
    }
    if (destination_failure == Real(1.0)) {
        std::ostringstream message;
        message.precision(17);
        message << "SBM coarse/fine stage produced an unrepresentable "
                << "U(t)/rho_d(t) quotient (level=" << level << ", stage="
                << stage << ", role=" << role << ", requested_time="
                << requested_time << ", coarse_old=" << parent.old_time
                << ", coarse_new=" << parent.new_time << ", rank="
                << amrex::ParallelDescriptor::MyProc() << ')';
        amrex::Abort(message.str());
    }

    // Retain the independent progress token as a second check.  It cannot
    // stand in for the actual-destination scan above.
    for (amrex::MFIter mfi(destination, false); mfi.isValid(); ++mfi) {
        const amrex::Box fab_box = mfi.fabbox();
        const auto mask = missing_fine_ghosts.array(mfi);
        const int ncomp = destination.nComp();
        amrex::ParallelFor(
            fab_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                if (mask(i, j, k, 0) == Real(1.0)) {
                    mask(i, j, k, 0) =
                        (mask(i, j, k, 1) == static_cast<Real>(ncomp))
                            ? Real(0.0) : Real(1.0);
                }
            });
    }
    if (missing_fine_ghosts.max(0, 2) != Real(0.0)) {
        std::ostringstream message;
        message.precision(17);
        message << "SBM coarse/fine stage omitted a required parent ghost "
                << "component write (level=" << level << ", stage=" << stage
                << ", role=" << role << ", requested_time=" << requested_time
                << ", coarse_old=" << parent.old_time << ", coarse_new="
                << parent.new_time << ", rank="
                << amrex::ParallelDescriptor::MyProc() << ')';
        amrex::Abort(message.str());
    }
#ifdef ERF_SBM_AMR_QUALIFICATION
    if (split_values_checked != nullptr &&
        parent.validate_split_ghost_values && current_role == 0) {
        // Once at the coarse start, and once strictly inside the coarse
        // interval, where both accepted endpoints contribute.
        const bool at_start = stage == 0 &&
            semantic_time_matches(requested_time, parent.old_time, dt_min);
        const bool interior = theta > 0.0 && theta < 1.0;
        const std::size_t slot = at_start ? 0U : 1U;
        if ((at_start || interior) && !(*split_values_checked)[slot]) {
            amrex::Gpu::streamSynchronize();
            qualification::verify_split_ghost_values(
                destination, fine_source_state, fine_source_conserved,
                coarse_old_spectrum, coarse_new_spectrum, parent, fine_geometry,
                theta, stage, requested_time);
            (*split_values_checked)[slot] = true;
        }
    }
#else
    amrex::ignore_unused(fine_source_state, fine_source_conserved,
                         split_values_checked);
#endif
    // Device component maps and parent chunks remain alive through the queued
    // ghost kernels. Synchronize before their storage is released.
    amrex::Gpu::streamSynchronize();
    if (profile_parent_fill) {
        std::uint64_t local_coarse_stage_bytes = 0;
        std::uint64_t local_coarse_pack_bytes = 0;
        for (amrex::MFIter mfi(coarse_stage, false); mfi.isValid(); ++mfi) {
            local_coarse_stage_bytes += static_cast<std::uint64_t>(
                coarse_stage[mfi].nBytesOwned());
        }
        if (coarse_pack.isDefined()) {
            for (amrex::MFIter mfi(coarse_pack, false);
                 mfi.isValid(); ++mfi) {
                local_coarse_pack_bytes += static_cast<std::uint64_t>(
                    coarse_pack[mfi].nBytesOwned());
            }
        }
        std::cout << std::setprecision(9)
                  << "SBM_AMR_PROFILE_PARENT rank="
                  << amrex::ParallelDescriptor::MyProc()
                  << " level=" << level << " stage=" << stage
                  << " role=" << role
                  << " ncomp=" << destination.nComp()
                  << " max_chunk_components=" << max_chunk_components
                  << " chunks=" << chunks.size()
                  << " parent_fill_mode="
                  << (componentwise_parent_fill ? "componentwise" : "packed")
                  << " parallel_copy_calls=" << parallel_copy_calls
                  << " parallel_copy_wall_s=" << parallel_copy_seconds
                  << " parent_fill_wall_s="
                  << std::chrono::duration<double>(ProfileClock::now() -
                                                    fill_start).count()
                  << " coarse_stage_fab_bytes_local="
                  << local_coarse_stage_bytes
                  << " coarse_source_pack_fab_bytes_local="
                  << local_coarse_pack_bytes
                  << " rank_live_fab_bytes="
                  << amrex::TotalBytesAllocatedInFabs()
                  << std::endl;
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

// Unit-test observation probes; each returns immediately unless a test
// supplies device sample storage for the requested stage.
void
sample_probe_endpoint_value (const SBMAMRStageGhostProbe* probe,
                             const int stage,
                             const MultiFab& endpoint_input)
{
    if (probe == nullptr || probe->sampled_endpoint_device == nullptr ||
        probe->endpoint_hit_device == nullptr ||
        (probe->requested_stage >= 0 && probe->requested_stage != stage)) {
        return;
    }
    const amrex::IntVect cell = probe->fine_cell;
    const int component = probe->global_component;
    auto* value = probe->sampled_endpoint_device;
    auto* hit = probe->endpoint_hit_device;
    for (amrex::MFIter mfi(endpoint_input, false); mfi.isValid(); ++mfi) {
        if (!mfi.fabbox().contains(cell)) continue;
        const auto endpoint = endpoint_input.const_array(mfi);
        amrex::ParallelFor(amrex::Box(cell, cell),
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                *value = endpoint(i, j, k, component);
                *hit = 1;
            });
    }
}

void
sample_probe_stencil (const SBMAMRStageGhostProbe* probe,
                      const int stage,
                      const MultiFab& physical_input,
                      const MultiFab& endpoint_input)
{
    if (probe == nullptr || probe->sampled_face_direction != 0 ||
        probe->sampled_physical_stencil_device == nullptr ||
        probe->sampled_endpoint_stencil_device == nullptr ||
        (probe->sampled_stage >= 0 && probe->sampled_stage != stage)) {
        return;
    }
    // Positive-x WENO-Z3 stencil of face i: cells i-2, i-1, i.
    amrex::IntVect lo = probe->sampled_face;
    lo[0] -= 2;
    const amrex::Box stencil(lo, probe->sampled_face);
    const int first = lo[0];
    const int mass = probe->sampled_mass_component;
    const int number = probe->sampled_number_component;
    const int property = probe->sampled_property_component;
    auto* physical_out = probe->sampled_physical_stencil_device;
    auto* endpoint_out = probe->sampled_endpoint_stencil_device;
    for (amrex::MFIter mfi(endpoint_input, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.fabbox() & stencil;
        if (!bx.ok()) continue;
        const auto physical = physical_input.const_array(mfi);
        const auto endpoint = endpoint_input.const_array(mfi);
        amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const int offset = 3 * (i - first);
                physical_out[offset] = physical(i, j, k, mass);
                physical_out[offset + 1] = physical(i, j, k, number);
                physical_out[offset + 2] = physical(i, j, k, property);
                endpoint_out[offset] = endpoint(i, j, k, mass);
                endpoint_out[offset + 1] = endpoint(i, j, k, number);
                endpoint_out[offset + 2] = endpoint(i, j, k, property);
            });
    }
}

void
sample_probe_face_rates (const SBMAMRStageGhostProbe* probe,
                         const int stage,
                         const std::vector<int>& members,
                         const std::vector<int>& global_to_local,
                         const erf_auxiliary::MappedFaceFluxRate& accepted_rate,
                         const erf_auxiliary::MappedFaceFluxRate& high_rate,
                         const erf_auxiliary::MappedFaceFluxRate& low_rate,
                         const erf_auxiliary::MappedFaceFluxRate& face_lambda)
{
    if (probe == nullptr || probe->sampled_face_direction < 0 ||
        (probe->sampled_stage >= 0 && probe->sampled_stage != stage) ||
        probe->sampled_mass_rate_device == nullptr ||
        probe->sampled_number_rate_device == nullptr ||
        probe->sampled_property_rate_device == nullptr ||
        probe->sampled_high_mass_rate_device == nullptr ||
        probe->sampled_high_number_rate_device == nullptr ||
        probe->sampled_high_property_rate_device == nullptr ||
        probe->sampled_low_mass_rate_device == nullptr ||
        probe->sampled_low_number_rate_device == nullptr ||
        probe->sampled_low_property_rate_device == nullptr ||
        probe->sampled_lambda_device == nullptr ||
        probe->sampled_rate_hit_device == nullptr) {
        return;
    }
    const auto local_member = [&] (const int component) {
        const bool member = std::find(members.begin(), members.end(),
                                      component) != members.end();
        return member
            ? global_to_local[static_cast<std::size_t>(component)] : -1;
    };
    const int m = local_member(probe->sampled_mass_component);
    const int n = local_member(probe->sampled_number_component);
    const int p = local_member(probe->sampled_property_component);
    if (m < 0 || n < 0 || p < 0) return;

    const int dir = probe->sampled_face_direction;
    const amrex::IntVect face = probe->sampled_face;
    const amrex::Box face_box(face, face,
                              amrex::IntVect::TheDimensionVector(dir));
    const auto& accepted = accepted_rate.dir(dir);
    const auto& high = high_rate.dir(dir);
    const auto& low = low_rate.dir(dir);
    const auto& lambda_field = face_lambda.dir(dir);
    auto* am = probe->sampled_mass_rate_device;
    auto* an = probe->sampled_number_rate_device;
    auto* ap = probe->sampled_property_rate_device;
    auto* hm = probe->sampled_high_mass_rate_device;
    auto* hn = probe->sampled_high_number_rate_device;
    auto* hp = probe->sampled_high_property_rate_device;
    auto* lm = probe->sampled_low_mass_rate_device;
    auto* ln = probe->sampled_low_number_rate_device;
    auto* lp = probe->sampled_low_property_rate_device;
    auto* out_lambda = probe->sampled_lambda_device;
    auto* out_hit = probe->sampled_rate_hit_device;
    for (amrex::MFIter mfi(accepted, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        if (!mfi.tilebox().contains(face)) continue;
        const auto accepted_values = accepted.const_array(mfi);
        const auto high_values = high.const_array(mfi);
        const auto low_values = low.const_array(mfi);
        const auto limiter = lambda_field.const_array(mfi);
        amrex::ParallelFor(face_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                am[0] = accepted_values(i, j, k, m);
                an[0] = accepted_values(i, j, k, n);
                ap[0] = accepted_values(i, j, k, p);
                hm[0] = high_values(i, j, k, m);
                hn[0] = high_values(i, j, k, n);
                hp[0] = high_values(i, j, k, p);
                lm[0] = low_values(i, j, k, m);
                ln[0] = low_values(i, j, k, n);
                lp[0] = low_values(i, j, k, p);
                out_lambda[0] = limiter(i, j, k, 0);
                out_hit[0] = 1;
            });
    }
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
    MultiFab coarse_fill_mask;
    MultiFab coarse_source_pack;
    MultiFab coarse_stage;
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
    SBMAMRFluxRegister* interface_as_coarse{nullptr};
    SBMAMRFluxRegister* interface_as_fine{nullptr};
    bool measure_ready{false};
    std::array<double, 3> profile_stage_seconds{};
    std::chrono::steady_clock::time_point profile_fine_step_start{};
};

struct SBMTransport::DeviceMetadata
{
    std::vector<DeviceGroupInfo> host_groups;
    amrex::Gpu::DeviceVector<ConstraintDescriptor> local_constraints;
    amrex::Gpu::DeviceVector<ConstraintTerm> local_terms;
    amrex::Gpu::DeviceVector<LocalPropertySupport> local_supports;
    amrex::Gpu::DeviceVector<PopulationRemapView> populations;
    amrex::Gpu::DeviceVector<Real> edges;
    amrex::Gpu::DeviceVector<Real> pivots;
    amrex::Gpu::DeviceVector<int> property_components;
};

SBMTransport::~SBMTransport () = default;

SBMTransport::SBMTransport (const SBMLayout& layout,
                            const int number_of_levels,
                            const int max_groups_per_chunk,
                            const bool keep_spectral_ledger_reference)
    : m_layout(layout),
      m_projection(layout),
      m_groups(make_constraint_groups(layout)),
      m_flat_constraints(make_constraint_descriptors(layout)),
      m_property_support(make_attached_property_support_descriptors(layout)),
      m_chunks(make_constraint_closure_chunks(layout, max_groups_per_chunk)),
      m_max_groups_per_chunk(max_groups_per_chunk),
      m_max_constraints_per_group(m_groups.size(), 0),
      m_levels(static_cast<std::size_t>(number_of_levels)),
      m_device_metadata(std::make_unique<DeviceMetadata>()),
      m_admission(std::make_unique<SBMStateAdmission>(layout)),
      m_keep_spectral_ledger_reference(keep_spectral_ledger_reference)
{
    AMREX_ALWAYS_ASSERT(number_of_levels > 0);
    if (AMREX_SPACEDIM != 3) {
        throw std::invalid_argument(
            "SBM M3 mapped transport currently requires three dimensions");
    }
    // Local slots follow global order, so YA adds each moment family as one run.
    for (auto& chunk : m_chunks) {
        std::sort(chunk.components.begin(), chunk.components.end());
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
        m_max_chunk_components = std::max(
            m_max_chunk_components,
            static_cast<int>(chunk.components.size()));
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
    m_admission->define_level(level, cell_ba, dm);
    const int ncomp = m_layout.ncomp();
    data->anchor.define(cell_ba, dm, ncomp, 0);
    data->target.define(cell_ba, dm, ncomp, 0);
    data->input_intensive.define(cell_ba, dm, ncomp, 2);
    data->physical_input_intensive.define(cell_ba, dm, ncomp, 2);
    data->base_intensive.define(cell_ba, dm, ncomp, 2);
    data->anchor_intensive.define(cell_ba, dm, ncomp, 2);
    if (level > 0) {
        data->coarse_fill_mask.define(cell_ba, dm, 3, 2);
    }
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
    if (m_keep_spectral_ledger_reference && m_levels.size() > 1) {
        data->spectral_ledger.define(cell_ba, dm, ncomp);
    }
    data->anchor.setVal(Real(0.0));
    data->target.setVal(Real(0.0));
    data->invalid.setVal(Real(0.0));
    data->projected_rate.setVal(Real(0.0));
    m_levels[static_cast<std::size_t>(level)] = std::move(data);
}

void
SBMTransport::set_amr_flux_registers (
    const int level, SBMAMRFluxRegister* as_coarse,
    SBMAMRFluxRegister* as_fine)
{
    if (!is_defined(level)) {
        throw std::logic_error("SBM transport level must be defined before binding AMR registers");
    }
    auto& data = *m_levels[static_cast<std::size_t>(level)];
    data.interface_as_coarse = as_coarse;
    data.interface_as_fine = as_fine;
}

bool
SBMTransport::begin_physical_step (const int level, const double step_begin,
                                   const double step_end,
                                   std::string& diagnostic)
{
    diagnostic.clear();
    if (!is_defined(level)) {
        diagnostic = "SBM transport physical step requires a defined level";
        return false;
    }
    auto& data = *m_levels[static_cast<std::size_t>(level)];
#ifdef ERF_SBM_AMR_QUALIFICATION
    if (level == 0) {
        m_interval_limiter_active_groups = 0;
        m_interval_limiter_changed_rate = false;
        m_interval_carrier_linf = Real(0.0);
        m_interval_accepted_rate_linf = Real(0.0);
    }
#endif
    if (data.interface_as_coarse &&
        !data.interface_as_coarse->begin_level_step(
            true, step_begin, step_end, diagnostic)) {
        return false;
    }
    if (data.interface_as_fine &&
        !data.interface_as_fine->begin_level_step(
            false, step_begin, step_end, diagnostic)) {
        return false;
    }
    return true;
}

bool
SBMTransport::accept_physical_step (const int level,
                                    const double step_begin,
                                    const double step_end,
                                    std::string& diagnostic)
{
    diagnostic.clear();
    if (!is_defined(level)) {
        diagnostic = "SBM transport physical-step acceptance requires a defined level";
        return false;
    }
    auto& data = *m_levels[static_cast<std::size_t>(level)];
    if (data.interface_as_coarse &&
        !data.interface_as_coarse->accept_level_step(
            true, step_begin, step_end, diagnostic)) {
        return false;
    }
    if (data.interface_as_fine &&
        !data.interface_as_fine->accept_level_step(
            false, step_begin, step_end, diagnostic)) {
        return false;
    }
    return true;
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
    const int qc_component, const int qr_component,
    const SBMAMRStageParent* parent)
{
    // Preserve ERF's semantic state roles at the integration boundary:
    // S_old anchors the step, S_new is the stage predictor, and S_data is the
    // target. In particular, do not substitute S_data for the predictor
    // density used to form the high-order intensive spectrum.
    if (!m_host_integrator_reported) {
        amrex::Print() << "SBM transport host integrator: "
                       << erf_auxiliary::HostIntegratorName(method) << "\n";
        m_host_integrator_reported = true;
    }
    advance_stage(level, method, stage, step_old_time, input_time, target_time,
                  host_stage_interval, state_manager, S_old_cons, S_new_cons,
                  S_data_cons, avg_xmom, avg_ymom, avg_zmom, geometry,
                  qc_component, qr_component, parent);
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
                             const int qr_component,
                             const SBMAMRStageParent* parent)
{
    AMREX_ALWAYS_ASSERT(is_defined(level));
    auto& data = *m_levels[static_cast<std::size_t>(level)];
    const bool profile_parent_fill =
        parent != nullptr && parent->profile_parent_fill;
    using ProfileClock = std::chrono::steady_clock;
    const auto profile_stage_start = profile_parent_fill
        ? ProfileClock::now() : ProfileClock::time_point{};
    if (profile_parent_fill && stage == 0) {
        data.profile_stage_seconds.fill(0.0);
        data.profile_fine_step_start = profile_stage_start;
    }
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
    if (data.interface_as_coarse &&
        !data.interface_as_coarse->begin_stage(
            true, method, stage, step_old_time, diagnostic)) {
        amrex::Abort("SBM YA coarse-side stage sequence: " + diagnostic);
    }
    if (data.interface_as_fine &&
        !data.interface_as_fine->begin_stage(
            false, method, stage, step_old_time, diagnostic)) {
        amrex::Abort("SBM YA fine-side stage sequence: " + diagnostic);
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
        std::ostringstream message;
        message.precision(17);
        message << "SBM M3 target density at level " << level << " stage "
                << stage << " (input=" << input_time << ", target="
                << target_time << ") is invalid: " << diagnostic
                << "; valid-cell min=" << conserved_target.min(Rho_comp)
                << " max=" << conserved_target.max(Rho_comp);
        amrex::Abort(message.str());
    }

#ifdef ERF_SBM_AMR_QUALIFICATION
    if (m_zero_flow_diagnostics_enabled) {
        m_interval_carrier_linf = std::max(
            m_interval_carrier_linf,
            std::max({avg_xmom.norminf(0), avg_ymom.norminf(0),
                      avg_zmom.norminf(0)}));
    }
#endif

    build_intensive(data.input_intensive, spectrum, conserved_input,
                    geometry.periodicity(), "input");
    build_intensive(data.base_intensive, trial_state, trial_conserved,
                    geometry.periodicity(), "limiter-trial base");
    build_intensive(data.anchor_intensive, data.anchor, conserved_anchor,
                    geometry.periodicity(), "old-step anchor");

    if (level > 0) {
        if (parent == nullptr || parent->conserved_old == nullptr ||
            parent->conserved_new == nullptr ||
            !state_manager.old_valid(level - 1) ||
            !state_manager.new_valid(level - 1)) {
            amrex::Abort("fine SBM stage is missing its accepted coarse "
                         "spectral/density endpoint bracket");
        }
        const MultiFab& coarse_old_spectrum = state_manager.old_state(level - 1);
        const MultiFab& coarse_new_spectrum = state_manager.new_state(level - 1);
        fill_coarse_parent_ghosts(
            data.input_intensive, data.coarse_fill_mask,
            coarse_old_spectrum, coarse_new_spectrum, spectrum,
            conserved_input, *parent, geometry,
            m_chunks, m_max_chunk_components, input_time,
            state_manager.old_time(level - 1),
            state_manager.new_time(level - 1), level, stage, "input",
            data.coarse_source_pack, data.coarse_stage,
            parent->ghost_probe, &m_split_ghost_values_checked);
        const double base_time =
            recipe.limiter_trial_base ==
                    erf_auxiliary::LimiterTrialBase::Anchor
                ? step_old_time : input_time;
        fill_coarse_parent_ghosts(
            data.base_intensive, data.coarse_fill_mask,
            coarse_old_spectrum, coarse_new_spectrum, trial_state,
            trial_conserved, *parent, geometry,
            m_chunks, m_max_chunk_components, base_time,
            state_manager.old_time(level - 1),
            state_manager.new_time(level - 1), level, stage,
            "limiter-base", data.coarse_source_pack, data.coarse_stage,
            parent->ghost_probe);
        fill_coarse_parent_ghosts(
            data.anchor_intensive, data.coarse_fill_mask,
            coarse_old_spectrum, coarse_new_spectrum, data.anchor,
            conserved_anchor, *parent, geometry,
            m_chunks, m_max_chunk_components, step_old_time,
            state_manager.old_time(level - 1),
            state_manager.new_time(level - 1), level, stage, "anchor",
            data.coarse_source_pack, data.coarse_stage,
            parent->ghost_probe);
    } else if (parent != nullptr) {
        amrex::Abort("coarse SBM stage unexpectedly received a parent bracket");
    }
    MultiFab::Copy(data.physical_input_intensive, data.input_intensive, 0, 0,
                   m_layout.ncomp(), 2);

    // Transform every 2M pair in place to endpoint-number coordinates.  The
    // canonical transform is shared with persisted-state reconstruction and
    // is intentionally performed before native scalar WENO-Z3.
    data.invalid.setVal(Real(0.0));
    const auto populations = m_device_metadata->populations.data();
    const int population_count =
        static_cast<int>(m_layout.populations().size());
    // A fine level also transforms the parent-filled coarse-only ghosts that
    // FillBoundary below does not refill; a single level needs valid cells only.
    const bool transform_ghosts = level > 0;
    MultiFab& transform_invalid =
        transform_ghosts ? data.coarse_fill_mask : data.invalid;
    const int transform_ghost_cells = transform_ghosts ? 2 : 0;
    if (transform_ghosts) {
        data.coarse_fill_mask.setVal(Real(0.0), 0, 1, transform_ghost_cells);
    }
    for (amrex::MFIter mfi(data.input_intensive,
                           transform_ghosts ? false : amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        const amrex::Box bx = transform_ghosts ? mfi.fabbox() : mfi.tilebox();
        const auto input = data.input_intensive.array(mfi);
        const auto invalid = transform_invalid.array(mfi);
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
    if (parent != nullptr) {
        sample_probe_endpoint_value(parent->ghost_probe, stage,
                                    data.input_intensive);
        for (const auto* probe :
             {parent->ghost_probe, parent->secondary_ghost_probe}) {
            sample_probe_stencil(probe, stage, data.physical_input_intensive,
                                 data.input_intensive);
        }
    }
    if (transform_invalid.max(0, transform_ghost_cells) != Real(0.0)) {
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
            std::ostringstream message;
            message.precision(17);
            message << "SBM M3 donor low-order trial violates a linear or "
                    << "attached-property/canonical-reserve constraint "
                    << "(level=" << level << ", stage=" << stage
                    << ", input_time=" << input_time
                    << ", target_time=" << target_time
                    << ", chunk=" << chunk_index << ", rank="
                    << amrex::ParallelDescriptor::MyProc() << ')';
            amrex::Abort(message.str());
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

            // A small-grid scientific oracle may observe one group's native
            // post-FCT face values.
            if (parent != nullptr) {
                for (const auto* probe :
                     {parent->ghost_probe, parent->secondary_ghost_probe}) {
                    sample_probe_face_rates(
                        probe, stage, group.members, global_to_local,
                        data.accepted_rate, data.high_rate, data.low_rate,
                        data.face_lambda);
                }
            }

#ifdef ERF_SBM_AMR_QUALIFICATION
            // Qualification-only evidence: record a group when any accepted
            // face limiter is below one, and separately prove that at least
            // one accepted face rate differs from its unlimited high rate.
            data.outgoing_demand.setVal(Real(0.0));
            for (amrex::MFIter mfi(data.outgoing_demand,
                                   amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi) {
                const amrex::Box bx = mfi.tilebox();
                const auto marker = data.outgoing_demand.array(mfi);
                const auto lambda_x = data.face_lambda.dir(0).const_array(mfi);
                const auto lambda_y = data.face_lambda.dir(1).const_array(mfi);
                const auto lambda_z = data.face_lambda.dir(2).const_array(mfi);
                amrex::ParallelFor(bx,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        int active = 0;
                        active |= lambda_x(i, j, k, 0) < Real(1.0) ||
                                  lambda_x(i + 1, j, k, 0) < Real(1.0);
                        active |= lambda_y(i, j, k, 0) < Real(1.0) ||
                                  lambda_y(i, j + 1, k, 0) < Real(1.0);
                        active |= lambda_z(i, j, k, 0) < Real(1.0) ||
                                  lambda_z(i, j, k + 1, 0) < Real(1.0);
                        marker(i, j, k, 0) = active ? Real(1.0) : Real(0.0);
                    });
            }
            for (const int component : group.transport_members) {
                const int local_component =
                    global_to_local[static_cast<std::size_t>(component)];
                AMREX_ALWAYS_ASSERT(local_component >= 0);
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    const auto& high_field = data.high_rate.dir(dir);
                    const auto& accepted_field = data.accepted_rate.dir(dir);
                    const auto& lambda_field = data.face_lambda.dir(dir);
                    for (amrex::MFIter mfi(data.outgoing_demand,
                                           amrex::TilingIfNotGPU());
                         mfi.isValid(); ++mfi) {
                        const amrex::Box bx = mfi.tilebox();
                        const auto high = high_field.const_array(mfi);
                        const auto accepted = accepted_field.const_array(mfi);
                        const auto lambda = lambda_field.const_array(mfi);
                        const auto marker = data.outgoing_demand.array(mfi);
                        amrex::ParallelFor(bx,
                            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                                int changed = 0;
                                for (int side = 0; side < 2; ++side) {
                                    const int fi = i + (dir == 0 && side != 0);
                                    const int fj = j + (dir == 1 && side != 0);
                                    const int fk = k + (dir == 2 && side != 0);
                                    if (lambda(fi, fj, fk, 0) < Real(1.0) &&
                                        accepted(fi, fj, fk, local_component) !=
                                            high(fi, fj, fk, local_component)) {
                                        changed = 1;
                                    }
                                }
                                if (changed != 0) marker(i, j, k, 0) = Real(2.0);
                            });
                    }
                }
            }
            const Real limiter_marker = data.outgoing_demand.max(0, 0);
            if (limiter_marker >= Real(1.0)) {
                ++m_interval_limiter_active_groups;
            }
            if (limiter_marker >= Real(2.0)) {
                m_interval_limiter_changed_rate = true;
            }
#endif

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

        // Qualification-only zero-flow evidence records the accepted face
        // rates after complete-group FCT, including mass, number and attached
        // components.  The expensive collective norms are disabled in all
        // other qualification cases and production builds.
#ifdef ERF_SBM_AMR_QUALIFICATION
        if (m_zero_flow_diagnostics_enabled) {
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                const auto& accepted = data.accepted_rate.dir(dir);
                m_interval_accepted_rate_linf = std::max(
                    m_interval_accepted_rate_linf,
                    accepted.norm0(0, accepted.nComp(), amrex::IntVect(0)));
            }
        }
#endif

        // Hand the accepted, chunk-local rates directly to the YA interface
        // register.  The owner keeps these additions pending until the full
        // physical step is admitted by ERF::Advance.
        if (data.interface_as_coarse &&
            !data.interface_as_coarse->add_stage_rates(
                true, data.accepted_rate, chunk.components,
                recipe.completed_ledger_time, diagnostic)) {
            amrex::Abort("SBM YA coarse-side accepted-rate accumulation: " + diagnostic);
        }
        if (data.interface_as_fine &&
            !data.interface_as_fine->add_stage_rates(
                false, data.accepted_rate, chunk.components,
                recipe.completed_ledger_time, diagnostic)) {
            amrex::Abort("SBM YA fine-side accepted-rate accumulation: " + diagnostic);
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

    // Reuse the same persistent full-state gate used by AMR publication.
    // Successful candidates are decided on device without a host FAB scan.
    if (!m_admission->admissible(data.target, level, diagnostic)) {
        amrex::Abort("SBM M3 candidate rejected before commit: " + diagnostic);
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
    const int stage_count = erf_auxiliary::HostIntegratorStageCount(method);
    const bool physical_step_complete = stage + 1 == stage_count;
    if (!state_manager.accept_stage_target(level, target_time,
                                          physical_step_complete, diagnostic)) {
        amrex::Abort("SBM M4a accepted-state lifecycle commit: " + diagnostic);
    }
    if (data.interface_as_coarse &&
        !data.interface_as_coarse->finish_stage(true, diagnostic)) {
        amrex::Abort("SBM YA coarse-side stage admission: " + diagnostic);
    }
    if (data.interface_as_fine &&
        !data.interface_as_fine->finish_stage(false, diagnostic)) {
        amrex::Abort("SBM YA fine-side stage admission: " + diagnostic);
    }
    state_manager.project_to_core(level, conserved_target, qc_component,
                                  qr_component);
    if (profile_parent_fill) {
        AMREX_ALWAYS_ASSERT(stage >= 0 &&
                            stage < static_cast<int>(data.profile_stage_seconds.size()));
        data.profile_stage_seconds[static_cast<std::size_t>(stage)] =
            std::chrono::duration<double>(ProfileClock::now() -
                                          profile_stage_start).count();
        std::cout << std::setprecision(9)
                  << "SBM_AMR_PROFILE_STAGE rank="
                  << amrex::ParallelDescriptor::MyProc()
                  << " level=" << level << " stage=" << stage
                  << " ncomp=" << m_layout.ncomp()
                  << " stage_wall_s="
                  << data.profile_stage_seconds[static_cast<std::size_t>(stage)]
                  << std::endl;
        if (physical_step_complete) {
            double stage_sum = 0.0;
            for (int s = 0; s < stage_count; ++s) {
                stage_sum += data.profile_stage_seconds[static_cast<std::size_t>(s)];
            }
            std::cout << std::setprecision(9)
                      << "SBM_AMR_PROFILE_FINE_STEP rank="
                      << amrex::ParallelDescriptor::MyProc()
                      << " level=" << level
                      << " ncomp=" << m_layout.ncomp()
                      << " fine_step_wall_s="
                      << std::chrono::duration<double>(
                             ProfileClock::now() -
                             data.profile_fine_step_start).count()
                      << " summed_stage_wall_s=" << stage_sum
                      << std::endl;
        }
    }
}

void
SBMTransport::destroy (const int level)
{
    AMREX_ALWAYS_ASSERT(level >= 0 &&
                        level < static_cast<int>(m_levels.size()));
    m_admission->destroy_level(level);
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

SBMStateAdmission& SBMTransport::admission () noexcept
{
    return *m_admission;
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
