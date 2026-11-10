#include "ERF_SBMAMRQualification.H"

#include "ERF.H"
#include "ERF_IndexDefines.H"
#include "ERF_SBMLayout.H"
#include "ERF_SBMStateManager.H"
#include "ERF_SBMTransport.H"

#include <AMReX_BoxIterator.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace amrex;

namespace {

struct SplitGhostValueSample
{
    bool coarse_parent{false};
    int depth{0};
    int codimension{0};
    int destination_box{-1};
    int destination_rank{-1};
    IntVect ghost{0};
    int donor_box{-1};
    int donor_rank{-1};
    IntVect donor_cell{0};
    int comparison_box{-1};
    IntVect comparison_cell{0};
};

IntVect
wrap_periodic_cell (IntVect cell, const Geometry& geometry)
{
    const auto domain = geometry.Domain();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (!geometry.isPeriodic(dir)) continue;
        const int length = domain.length(dir);
        while (cell[dir] < domain.smallEnd(dir)) cell[dir] += length;
        while (cell[dir] > domain.bigEnd(dir)) cell[dir] -= length;
    }
    return cell;
}

int
containing_box (const BoxArray& boxes, const IntVect& cell)
{
    for (int index = 0; index < boxes.size(); ++index) {
        if (boxes[index].contains(cell)) return index;
    }
    return -1;
}

std::pair<int, IntVect>
adjacent_valid_cell (const IntVect& source, const BoxArray& boxes)
{
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        for (int sign : {1, -1}) {
            IntVect candidate = source;
            candidate[dir] += sign;
            const int owner_box = containing_box(boxes, candidate);
            if (owner_box >= 0) return {owner_box, candidate};
        }
    }
    return {-1, IntVect(0)};
}

std::vector<SplitGhostValueSample>
select_split_ghost_value_samples (const BoxArray& coarse_ba,
                                  const DistributionMapping& coarse_dm,
                                  const Geometry& fine_geometry,
                                  const BoxArray& fine_ba,
                                  const DistributionMapping& fine_dm,
                                  const IntVect& refinement_ratio)
{
    // Keep one remote sample for each halo depth and ghost codimension when
    // the fixed layout contains one. The split test requires face samples at
    // depths one and two; edges/corners are additionally emitted if present.
    std::array<std::array<std::array<bool, 4>, 2>, 2> found{};
    std::vector<SplitGhostValueSample> samples;
    for (int destination_box = 0; destination_box < fine_ba.size();
         ++destination_box) {
        const auto valid = fine_ba[destination_box];
        auto grown = valid;
        grown.grow(2);
        const int destination_rank = fine_dm[destination_box];
        for (BoxIterator bit(grown); bit.ok(); ++bit) {
            const auto ghost = bit();
            if (valid.contains(ghost)) continue;
            int depth = 0;
            int codimension = 0;
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                int outside = 0;
                if (ghost[dir] < valid.smallEnd(dir)) {
                    outside = valid.smallEnd(dir) - ghost[dir];
                } else if (ghost[dir] > valid.bigEnd(dir)) {
                    outside = ghost[dir] - valid.bigEnd(dir);
                }
                if (outside > 0) {
                    ++codimension;
                    depth = std::max(depth, outside);
                }
            }
            if (depth < 1 || depth > 2 || codimension < 1 ||
                codimension > AMREX_SPACEDIM) {
                continue;
            }
            const auto wrapped = wrap_periodic_cell(ghost, fine_geometry);
            const int fine_donor = containing_box(fine_ba, wrapped);
            const bool remote_fine = fine_donor >= 0 &&
                fine_donor != destination_box &&
                fine_dm[fine_donor] != destination_rank;
            IntVect donor_cell = wrapped;
            int donor_box = fine_donor;
            bool coarse_parent = false;
            int donor_rank = remote_fine ? fine_dm[fine_donor] : -1;
            if (fine_donor < 0) {
                auto parent_cell = wrapped;
                parent_cell.coarsen(refinement_ratio);
                const int parent_box = containing_box(coarse_ba, parent_cell);
                if (parent_box < 0 || coarse_dm[parent_box] == destination_rank) {
                    continue;
                }
                coarse_parent = true;
                donor_cell = parent_cell;
                donor_box = parent_box;
                donor_rank = coarse_dm[parent_box];
            } else if (!remote_fine) {
                continue;
            }
            const int type = coarse_parent ? 0 : 1;
            if (found[type][depth - 1][codimension]) continue;
            const auto comparison = adjacent_valid_cell(
                donor_cell, coarse_parent ? coarse_ba : fine_ba);
            if (comparison.first < 0) continue;
            SplitGhostValueSample sample;
            sample.coarse_parent = coarse_parent;
            sample.depth = depth;
            sample.codimension = codimension;
            sample.destination_box = destination_box;
            sample.destination_rank = destination_rank;
            sample.ghost = ghost;
            sample.donor_box = donor_box;
            sample.donor_rank = donor_rank;
            sample.donor_cell = donor_cell;
            sample.comparison_box = comparison.first;
            sample.comparison_cell = comparison.second;
            samples.push_back(sample);
            found[type][depth - 1][codimension] = true;
        }
    }
    for (int type = 0; type < 2; ++type) {
        for (int depth = 1; depth <= 2; ++depth) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                found[type][depth - 1][1],
                type == 0
                    ? "split SBM qualification lacks a remote coarse-parent face ghost at a required depth"
                    : "split SBM qualification lacks a remote fine-fine face ghost at a required depth");
        }
    }
    return samples;
}

Real
read_fab_value (const MultiFab& field, const int box_index,
                const IntVect& cell, const int component)
{
    for (MFIter mfi(field, false); mfi.isValid(); ++mfi) {
        if (mfi.index() != box_index) continue;
        const auto values = field.const_array(mfi);
        Real result = Real(0.0);
        Gpu::copy(Gpu::deviceToHost, values.ptr(cell, component),
                  values.ptr(cell, component) + 1, &result);
        Gpu::streamSynchronize();
        return result;
    }
    return Real(0.0);
}

void
validate_split_layout (const BoxArray& coarse_ba,
                       const DistributionMapping& coarse_dm,
                       const BoxArray& fine_ba,
                       const DistributionMapping& fine_dm,
                       const Geometry& fine_geom,
                       const IntVect& refinement_ratio)
{
    const int ranks = ParallelDescriptor::NProcs();
    std::vector<int> coarse_boxes(static_cast<std::size_t>(ranks), 0);
    std::vector<int> fine_boxes(static_cast<std::size_t>(ranks), 0);
    std::vector<Long> coarse_cells(static_cast<std::size_t>(ranks), 0);
    std::vector<Long> fine_cells(static_cast<std::size_t>(ranks), 0);
    for (int box = 0; box < coarse_ba.size(); ++box) {
        const int owner = coarse_dm[box];
        ++coarse_boxes[static_cast<std::size_t>(owner)];
        coarse_cells[static_cast<std::size_t>(owner)] += coarse_ba[box].numPts();
    }
    for (int box = 0; box < fine_ba.size(); ++box) {
        const int owner = fine_dm[box];
        ++fine_boxes[static_cast<std::size_t>(owner)];
        fine_cells[static_cast<std::size_t>(owner)] += fine_ba[box].numPts();
    }
    int coarse_owner_count = 0;
    int fine_owner_count = 0;
    for (int rank = 0; rank < ranks; ++rank) {
        coarse_owner_count += coarse_boxes[static_cast<std::size_t>(rank)] > 0;
        fine_owner_count += fine_boxes[static_cast<std::size_t>(rank)] > 0;
        std::ostringstream line;
        line << "SBM_AMR_QUAL_OWNER rank=" << rank
             << " coarse_boxes=" << coarse_boxes[static_cast<std::size_t>(rank)]
             << " coarse_cells=" << coarse_cells[static_cast<std::size_t>(rank)]
             << " fine_boxes=" << fine_boxes[static_cast<std::size_t>(rank)]
             << " fine_cells=" << fine_cells[static_cast<std::size_t>(rank)];
        AllPrint() << line.str() << '\n';
    }

    bool found_remote_parent = false;
    bool found_remote_fine_donor = false;
    std::ostringstream parent_sample;
    std::ostringstream fine_sample;
    for (int destination_box = 0;
         destination_box < fine_ba.size() &&
             (!found_remote_parent || !found_remote_fine_donor);
         ++destination_box) {
        const int destination_rank = fine_dm[destination_box];
        const Box valid = fine_ba[destination_box];
        Box grown = valid;
        grown.grow(2);
        for (BoxIterator bit(grown);
             bit.ok() && (!found_remote_parent || !found_remote_fine_donor);
             ++bit) {
            const IntVect raw_ghost = bit();
            if (valid.contains(raw_ghost)) continue;
            const IntVect wrapped = wrap_periodic_cell(raw_ghost, fine_geom);
            const int fine_donor = containing_box(fine_ba, wrapped);
            if (!found_remote_fine_donor && fine_donor >= 0 &&
                fine_donor != destination_box &&
                fine_dm[fine_donor] != destination_rank) {
                found_remote_fine_donor = true;
                fine_sample << "dest_fine_box=" << destination_box
                            << " dest_rank=" << destination_rank
                            << " ghost=" << raw_ghost
                            << " donor_fine_box=" << fine_donor
                            << " donor_rank=" << fine_dm[fine_donor]
                            << " donor_cell=" << wrapped;
            }
            if (!found_remote_parent && fine_donor < 0) {
                IntVect parent = wrapped;
                parent.coarsen(refinement_ratio);
                const int parent_box = containing_box(coarse_ba, parent);
                if (parent_box >= 0 && coarse_dm[parent_box] != destination_rank) {
                    found_remote_parent = true;
                    parent_sample << "dest_fine_box=" << destination_box
                                  << " dest_rank=" << destination_rank
                                  << " ghost=" << raw_ghost
                                  << " parent_cell=" << parent
                                  << " source_coarse_box=" << parent_box
                                  << " source_rank=" << coarse_dm[parent_box];
                }
            }
        }
    }

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        ranks >= 2 && coarse_owner_count >= 2 && fine_owner_count >= 2,
        "split SBM qualification must place coarse and fine valid cells on "
        "at least two ranks");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        found_remote_parent,
        "split SBM qualification found no remote coarse parent ghost donor");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        found_remote_fine_donor,
        "split SBM qualification found no remote same-level fine-fine ghost donor");
    Print() << "SBM_AMR_QUAL_SPLIT_LAYOUT coarse_ranks="
            << coarse_owner_count << " fine_ranks=" << fine_owner_count
            << " remote_parent=true remote_fine=true\n"
            << "SBM_AMR_QUAL_GHOST_SAMPLE type=coarse_parent "
            << parent_sample.str() << '\n'
            << "SBM_AMR_QUAL_GHOST_SAMPLE type=fine_fine "
            << fine_sample.str() << '\n';
}

void
validate_empty_fine_rank_layout (const BoxArray& fine_ba,
                                 const DistributionMapping& fine_dm)
{
    const int ranks = ParallelDescriptor::NProcs();
    std::vector<int> fine_boxes(static_cast<std::size_t>(ranks), 0);
    std::vector<Long> fine_cells(static_cast<std::size_t>(ranks), 0);
    for (int box = 0; box < fine_ba.size(); ++box) {
        const int owner = fine_dm[box];
        ++fine_boxes[static_cast<std::size_t>(owner)];
        fine_cells[static_cast<std::size_t>(owner)] += fine_ba[box].numPts();
    }
    int empty_ranks = 0;
    for (int rank = 0; rank < ranks; ++rank) {
        if (fine_boxes[static_cast<std::size_t>(rank)] == 0) ++empty_ranks;
        std::ostringstream line;
        line << "SBM_AMR_QUAL_EMPTY_OWNER rank=" << rank
             << " fine_boxes=" << fine_boxes[static_cast<std::size_t>(rank)]
             << " fine_cells=" << fine_cells[static_cast<std::size_t>(rank)];
        AllPrint() << line.str() << '\n';
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        ranks >= 2 && empty_ranks > 0,
        "empty-fine-rank SBM qualification did not leave a rank without fine cells");
    Print() << "SBM_AMR_QUAL_EMPTY_FINE_RANK empty_ranks=" << empty_ranks << '\n';
}

struct InventoryRoundoff
{
    std::uint64_t leaf_cell_terms{0};
    std::uint64_t operation_count{0};
    int reduction_depth{0};
    Real rounding_factor{0};
};

// Operation-count bound on the composite leaf inventory, shared by the
// per-step and regrid conservation checks.
InventoryRoundoff
composite_inventory_roundoff (const std::uint64_t coarse_cells,
                              const std::uint64_t fine_cells,
                              const IntVect& ratio,
                              const bool with_prolongation)
{
    std::uint64_t refinement_volume = 1U;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        refinement_volume *= static_cast<std::uint64_t>(ratio[dir]);
    }
    AMREX_ALWAYS_ASSERT(refinement_volume > 0U &&
                        fine_cells % refinement_volume == 0U &&
                        coarse_cells >= fine_cells / refinement_volume);
    InventoryRoundoff roundoff;
    roundoff.leaf_cell_terms = coarse_cells - fine_cells / refinement_volume + fine_cells;
    for (int ranks = ParallelDescriptor::NProcs(); ranks > 1;
         ranks = (ranks + 1) / 2) {
        ++roundoff.reduction_depth;
    }
    // Each positive leaf term uses up to three multiplies and one local
    // accumulation. Bound both initial and final inventories, the final
    // subtraction, and a pairwise MPI reduction with a conservative
    // eight operations per leaf plus two reduction-tree levels per depth;
    // a regrid adds the prolongation of a new fine cell from its parent.
    roundoff.operation_count =
        8U * roundoff.leaf_cell_terms +
        2U * static_cast<std::uint64_t>(roundoff.reduction_depth) + 1U +
        (with_prolongation ? refinement_volume + 6U : 0U);
    const Real q = static_cast<Real>(roundoff.operation_count) *
                   std::numeric_limits<Real>::epsilon();
    AMREX_ALWAYS_ASSERT(q < Real(1.0));
    roundoff.rounding_factor = q / (Real(1.0) - q);
    return roundoff;
}

} // namespace

namespace erf_sbm::qualification {

void
validate_configuration (const SolverChoice& choice, const int max_level)
{
    if (max_level == 0) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            !choice.sbm_amr_qualification_test &&
            !choice.sbm_amr_qualification_expect_zero_flow &&
            !choice.sbm_amr_qualification_expect_split &&
            !choice.sbm_amr_qualification_expect_empty_fine_rank &&
            !choice.sbm_amr_qualification_unique_profile &&
            !choice.sbm_amr_qualification_parent_store_fault &&
            !choice.sbm_amr_qualification_regrid &&
            !choice.sbm_amr_qualification_regrid_reject_fault &&
            !choice.sbm_amr_qualification_uniform_ratio,
            "the internal SBM AMR qualification flag requires amr.max_level=1");
        return;
    }
    // The production capability check (ERF_Constructors.cpp) has already
    // restricted the hierarchy, mesh, coupling, and integrators.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(max_level == 1 &&
                                     choice.sbm_amr_qualification_test,
        "the isolated SBM AMR qualification build permits only its "
        "explicit fixed two-level test deck");
    ParmParse pp_erf("erf");
    int regrid_interval = -1;
    pp_erf.query("regrid_int", regrid_interval);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        choice.sbm_amr_qualification_regrid ? regrid_interval > 0
                                            : regrid_interval == 0,
        "SBM AMR qualification requires erf.regrid_int>0 with "
        "erf.sbm_amr_qualification_regrid and erf.regrid_int=0 otherwise");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !choice.sbm_amr_qualification_regrid_reject_fault ||
        choice.sbm_amr_qualification_regrid,
        "the SBM regrid rejection fault requires erf.sbm_amr_qualification_regrid");
    const bool all_levels_anelastic_rk2 =
        !choice.compressibility_type.empty() &&
        std::all_of(choice.compressibility_type.begin(),
                    choice.compressibility_type.end(),
                    [] (CompressibilityType type) {
                        return type == CompressibilityType::Anelastic;
                    }) &&
        std::all_of(choice.timestepping_type.begin(),
                    choice.timestepping_type.end(),
                    [] (TimesteppingType type) {
                        return type == TimesteppingType::RK2;
                    });
    // The compressible host leaves a roundoff carrier at zero wind.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !choice.sbm_amr_qualification_expect_zero_flow ||
        all_levels_anelastic_rk2,
        "SBM zero-flow AMR qualification requires the exact-zero "
        "anelastic RK2 carrier on every level");
    // A uniform dry-air-relative ratio needs a common factor-free profile.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !choice.sbm_amr_qualification_uniform_ratio ||
        (!choice.sbm_amr_qualification_unique_profile &&
         !choice.sbm_amr_qualification_expect_zero_flow),
        "SBM uniform-ratio AMR qualification excludes the unique and "
        "zero-flow profiles");
}

void
add_layout_properties (SBMLayoutSpec& spec, const SolverChoice& choice)
{
    if (!choice.sbm_amr_qualification_test || choice.sbm_moment_mode != 2) {
        return;
    }
    AttachedPropertyDescriptor property;
    property.name = "qualification_number_tracer";
    property.semantic_id = "qualification_number_tracer";
    property.units = "m^-3";
    property.carrier_population = 0;
    property.kind = PropertyKind::NumberCarried;
    property.support = SupportRequirement::PositiveNumber;
    property.support_min = Real(0.0);
    spec.attached_properties.push_back(std::move(property));
}

void
verify_split_ghost_values (const MultiFab& destination,
                           const MultiFab& fine_source_state,
                           const MultiFab& fine_source_conserved,
                           const MultiFab& coarse_old_spectrum,
                           const MultiFab& coarse_new_spectrum,
                           const SBMAMRStageParent& parent,
                           const Geometry& fine_geometry,
                           const double theta,
                           const int stage,
                           const double requested_time)
{
    const auto samples = select_split_ghost_value_samples(
        coarse_old_spectrum.boxArray(), coarse_old_spectrum.DistributionMap(),
        fine_geometry, destination.boxArray(), destination.DistributionMap(),
        parent.refinement_ratio);
    const std::array<int, 2> components{0, destination.nComp() - 1};
    const Real tolerance = Real(512.0) * std::numeric_limits<Real>::epsilon();
    const int my_rank = ParallelDescriptor::MyProc();
    for (const auto& sample : samples) {
        for (const int component : components) {
            std::array<Real, 9> values{};
            Real actual = Real(0.0);
            if (my_rank == sample.destination_rank) {
                actual = read_fab_value(destination, sample.destination_box,
                                        sample.ghost, component);
            }
            if (sample.coarse_parent) {
                if (my_rank == sample.donor_rank) {
                    values[0] = read_fab_value(coarse_old_spectrum, sample.donor_box,
                                               sample.donor_cell, component);
                    values[1] = read_fab_value(coarse_new_spectrum, sample.donor_box,
                                               sample.donor_cell, component);
                    values[2] = read_fab_value(*parent.conserved_old, sample.donor_box,
                                               sample.donor_cell, Rho_comp);
                    values[3] = read_fab_value(*parent.conserved_new, sample.donor_box,
                                               sample.donor_cell, Rho_comp);
                }
                if (my_rank ==
                    coarse_old_spectrum.DistributionMap()[sample.comparison_box]) {
                    values[4] = read_fab_value(coarse_old_spectrum, sample.comparison_box,
                                               sample.comparison_cell, component);
                    values[5] = read_fab_value(coarse_new_spectrum, sample.comparison_box,
                                               sample.comparison_cell, component);
                    values[6] = read_fab_value(*parent.conserved_old, sample.comparison_box,
                                               sample.comparison_cell, Rho_comp);
                    values[7] = read_fab_value(*parent.conserved_new, sample.comparison_box,
                                               sample.comparison_cell, Rho_comp);
                }
                values[8] = actual;
                ParallelDescriptor::ReduceRealSum(values.data(),
                                                  static_cast<int>(values.size()));
                const Real theta_real = static_cast<Real>(theta);
                const Real rho = (Real(1.0) - theta_real) * values[2] +
                                 theta_real * values[3];
                const Real expected =
                    ((Real(1.0) - theta_real) * values[2] * values[0] +
                     theta_real * values[3] * values[1]) / rho;
                const Real comparison_rho =
                    (Real(1.0) - theta_real) * values[6] + theta_real * values[7];
                const Real comparison_expected =
                    ((Real(1.0) - theta_real) * values[6] * values[4] +
                     theta_real * values[7] * values[5]) / comparison_rho;
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    std::isfinite(values[8]) &&
                        std::abs(values[8] - expected) <=
                            tolerance * std::max(Real(1.0), std::abs(expected)),
                    "split SBM coarse-parent received value differs from independently reconstructed U(t)/rho(t)");
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    std::abs(expected - comparison_expected) >
                        Real(100.0) * tolerance *
                            std::max(Real(1.0), std::abs(expected)),
                    "split SBM coarse-parent spatial fixture did not distinguish adjacent donors");
                Print() << std::setprecision(17)
                        << "SBM_AMR_QUAL_GHOST_VALUE type=coarse_parent"
                        << " depth=" << sample.depth
                        << " codim=" << sample.codimension
                        << " component=" << component
                        << " dest_rank=" << sample.destination_rank
                        << " dest_box=" << sample.destination_box
                        << " ghost=" << sample.ghost
                        << " donor_rank=" << sample.donor_rank
                        << " donor_box=" << sample.donor_box
                        << " donor_cell=" << sample.donor_cell
                        << " requested_time=" << requested_time
                        << " received=" << values[8]
                        << " expected=" << expected
                        << " neighbor_expected=" << comparison_expected
                        << " contrast=" << std::abs(expected - comparison_expected)
                        << '\n';
            } else {
                if (my_rank == sample.donor_rank) {
                    values[0] = read_fab_value(fine_source_state, sample.donor_box,
                                               sample.donor_cell, component);
                    values[1] = read_fab_value(fine_source_conserved, sample.donor_box,
                                               sample.donor_cell, Rho_comp);
                }
                if (my_rank ==
                    fine_source_state.DistributionMap()[sample.comparison_box]) {
                    values[2] = read_fab_value(fine_source_state, sample.comparison_box,
                                               sample.comparison_cell, component);
                    values[3] = read_fab_value(fine_source_conserved, sample.comparison_box,
                                               sample.comparison_cell, Rho_comp);
                }
                values[4] = actual;
                ParallelDescriptor::ReduceRealSum(values.data(),
                                                  static_cast<int>(values.size()));
                const Real expected = values[0] / values[1];
                const Real comparison_expected = values[2] / values[3];
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    std::isfinite(values[4]) &&
                        std::abs(values[4] - expected) <=
                            tolerance * std::max(Real(1.0), std::abs(expected)),
                    "split SBM fine-fine received value differs from its remote valid-source intensive value");
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    std::abs(expected - comparison_expected) >
                        Real(100.0) * tolerance *
                            std::max(Real(1.0), std::abs(expected)),
                    "split SBM fine-fine spatial fixture did not distinguish adjacent donors");
                Print() << std::setprecision(17)
                        << "SBM_AMR_QUAL_GHOST_VALUE type=fine_fine"
                        << " depth=" << sample.depth
                        << " codim=" << sample.codimension
                        << " component=" << component
                        << " dest_rank=" << sample.destination_rank
                        << " dest_box=" << sample.destination_box
                        << " ghost=" << sample.ghost
                        << " donor_rank=" << sample.donor_rank
                        << " donor_box=" << sample.donor_box
                        << " donor_cell=" << sample.donor_cell
                        << " requested_time=" << requested_time
                        << " received=" << values[4]
                        << " expected=" << expected
                        << " neighbor_expected=" << comparison_expected
                        << " contrast=" << std::abs(expected - comparison_expected)
                        << '\n';
            }
        }
    }
    Print() << "SBM_AMR_QUAL_GHOST_VALUE_COMPLETE stage=" << stage
            << " requested_time=" << requested_time
            << " theta=" << theta
            << " samples=" << samples.size() << '\n';
}

IntVect
first_coarse_only_ghost (const MultiFab& mask)
{
    for (MFIter mfi(mask, false); mfi.isValid(); ++mfi) {
        const Box fab_box = mfi.fabbox();
        const auto npts = static_cast<std::size_t>(fab_box.numPts());
        std::vector<Real> host(npts);
        const Real* device = mask[mfi].dataPtr(0);
        Gpu::copy(Gpu::deviceToHost, device, device + npts, host.begin());
        Gpu::streamSynchronize();
        for (std::size_t offset = 0; offset < npts; ++offset) {
            if (host[offset] == Real(1.0)) {
                return fab_box.atOffset(static_cast<Long>(offset));
            }
        }
    }
    return IntVect(std::numeric_limits<int>::min());
}

} // namespace erf_sbm::qualification

void
ERF::sbm_qualification_apply_initial_profile (int lev)
{
    if (!solverChoice.sbm_amr_qualification_test) return;
    auto& spectrum = sbm_state_manager->new_state_for_initialization(lev);
    const int ncomp = spectrum.nComp();
    if (solverChoice.sbm_amr_qualification_unique_profile) {
        // Split-layout MPI oracle: a deterministic, spatially unique common
        // factor per cell keeps every 2M and attached-property support ratio.
        const Box domain = geom[lev].Domain();
        const Long nx = domain.length(0);
        const Long ny = domain.length(1);
        const Long nz = domain.length(2);
        const Long cells = nx * ny * nz;
        for (MFIter mfi(spectrum, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box bx = mfi.tilebox();
            const auto values = spectrum.array(mfi);
            ParallelFor(bx, ncomp,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    const Long ii = static_cast<Long>(i - domain.smallEnd(0));
                    const Long jj = static_cast<Long>(j - domain.smallEnd(1));
                    const Long kk = static_cast<Long>(k - domain.smallEnd(2));
                    const Long index = ii + nx * (jj + ny * kk);
                    const Real factor = Real(0.5) + Real(0.5) *
                        static_cast<Real>(index) / static_cast<Real>(cells);
                    values(i, j, k, n) *= factor;
                });
        }
    } else if (lev == 0 && !solverChoice.sbm_amr_qualification_expect_zero_flow &&
               !solverChoice.sbm_amr_qualification_uniform_ratio) {
        // A common-factor jump at the left edge of the fixed refined box gives
        // a sharp periodic translation profile with unchanged support ratios.
        const int jump_x = geom[lev].Domain().smallEnd(0) +
                           geom[lev].Domain().length(0) / 4;
        for (MFIter mfi(spectrum, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box bx = mfi.tilebox();
            const auto values = spectrum.array(mfi);
            ParallelFor(bx, ncomp,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    const Real factor = i < jump_x ? Real(0.1) : Real(1.0);
                    values(i, j, k, n) *= factor;
                });
        }
    }
}

std::vector<Real>
ERF::sbm_qualification_composite_inventory () const
{
    AMREX_ALWAYS_ASSERT(sbm_state_manager != nullptr && sbm_transport != nullptr);
    const int ncomp = sbm_state_manager->layout().ncomp();
    std::vector<Real> inventory(static_cast<std::size_t>(ncomp), Real(0.0));
    for (int component = 0; component < ncomp; ++component) {
        Real total = Real(0.0);
        for (int lev = 0; lev <= finest_level; ++lev) {
            // The SBM level layouts, not grids[], define the composite: during
            // a remake grids[lev] already names the new grids.
            const MultiFab& spectrum = sbm_state_manager->new_state(lev);
            const BoxArray& level_ba = spectrum.boxArray();
            const DistributionMapping& level_dm = spectrum.DistributionMap();
            MultiFab leaf_mask(level_ba, level_dm, 1, 0);
            leaf_mask.setVal(Real(1.0));
            if (lev < finest_level) {
                const MultiFab& fine_spectrum = sbm_state_manager->new_state(lev + 1);
                MultiFab fine_coverage(fine_spectrum.boxArray(),
                                       fine_spectrum.DistributionMap(), 1, 0);
                MultiFab covered(level_ba, level_dm, 1, 0);
                fine_coverage.setVal(Real(1.0));
                covered.setVal(Real(0.0));
                average_down(fine_coverage, covered, 0, 1, ref_ratio[lev]);
                for (MFIter mfi(leaf_mask, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box bx = mfi.tilebox();
                    const auto leaf = leaf_mask.array(mfi);
                    const auto mask = covered.const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        leaf(i, j, k, 0) = Real(1.0) - mask(i, j, k, 0);
                    });
                }
            }

            MultiFab weighted(level_ba, level_dm, 1, 0);
            const MultiFab& measure = sbm_transport->static_measure(lev);
            const Real cell_volume = geom[lev].CellSize(0) *
                                     geom[lev].CellSize(1) *
                                     geom[lev].CellSize(2);
            for (MFIter mfi(weighted, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box bx = mfi.tilebox();
                const auto u = spectrum.const_array(mfi);
                const auto omega = measure.const_array(mfi);
                const auto leaf = leaf_mask.const_array(mfi);
                const auto out = weighted.array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    out(i, j, k, 0) = cell_volume * omega(i, j, k, 0) *
                                      u(i, j, k, component) * leaf(i, j, k, 0);
                });
            }
            total += weighted.sum(0, false);
        }
        inventory[static_cast<std::size_t>(component)] = total;
    }
    return inventory;
}

void
ERF::sbm_qualification_record_initial_inventory ()
{
    if (!solverChoice.sbm_amr_qualification_test) return;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        finest_level == 1 ||
        (finest_level == 0 && solverChoice.sbm_amr_qualification_regrid),
        "SBM AMR qualification did not create its fixed fine level");
    sbm_qualification_initial_inventory = sbm_qualification_composite_inventory();
    Print() << "SBM_AMR_QUAL_INIT levels=" << finest_level + 1
            << " ncomp=" << sbm_state_manager->layout().ncomp() << '\n';
    for (std::size_t component = 0;
         component < sbm_qualification_initial_inventory.size(); ++component) {
        Print().SetPrecision(17)
            << "SBM_AMR_QUAL_INITIAL_INVENTORY component=" << component
            << " value=" << sbm_qualification_initial_inventory[component]
            << '\n';
    }
}

void
ERF::sbm_qualification_check_step (int nstep, Real correction_linf,
                                   std::uint64_t reflux_calls,
                                   std::uint64_t reflux_syncs)
{
    if (!solverChoice.sbm_amr_qualification_test) return;
    // A regrid deck may clear its fine level and later recreate it.  The
    // runner checks the fine-substep count each deck selects.
    const bool has_fine_level = finest_level == 1;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        has_fine_level ||
        (finest_level == 0 && solverChoice.sbm_amr_qualification_regrid),
        "SBM AMR qualification lost its fine level outside a regrid deck");
    if (solverChoice.sbm_amr_qualification_expect_split && nstep == 0) {
        validate_split_layout(grids[0], dmap[0], grids[1], dmap[1], geom[1],
                              ref_ratio[0]);
    } else if (solverChoice.sbm_amr_qualification_expect_empty_fine_rank &&
               nstep == 0) {
        validate_empty_fine_rank_layout(grids[1], dmap[1]);
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !has_fine_level || (reflux_calls > 0 && reflux_syncs > 0),
        "SBM AMR qualification did not exercise and synchronize its YA register");
    const bool limiter_changed_rate = sbm_transport->interval_limiter_changed_rate();
    const auto coarse_cells = static_cast<std::uint64_t>(grids[0].numPts());
    const auto fine_cells = has_fine_level
        ? static_cast<std::uint64_t>(grids[1].numPts()) : std::uint64_t(0);
    const auto roundoff = composite_inventory_roundoff(
        coarse_cells, fine_cells, ref_ratio[0], false);
    const std::uint64_t leaf_cell_terms = roundoff.leaf_cell_terms;
    const std::uint64_t inventory_operation_count = roundoff.operation_count;
    const int reduction_depth = roundoff.reduction_depth;
    const Real rounding_factor = roundoff.rounding_factor;
    const bool zero_flow = solverChoice.sbm_amr_qualification_expect_zero_flow;
    Real zero_flow_correction_tolerance = Real(0.0);
    if (zero_flow) {
        const Real max_carrier = sbm_transport->interval_carrier_linf();
        const Real max_accepted_rate = sbm_transport->interval_accepted_rate_linf();
        Real inventory_scale = std::numeric_limits<Real>::min();
        for (const Real value : sbm_qualification_initial_inventory) {
            inventory_scale = std::max(inventory_scale, std::abs(value));
        }
        zero_flow_correction_tolerance = rounding_factor * inventory_scale;
        Print() << "SBM_AMR_QUAL_ZERO_FLOW_CHECK correction_linf=" << correction_linf
                << " tolerance=" << zero_flow_correction_tolerance
                << " max_carrier=" << max_carrier
                << " max_accepted_rate=" << max_accepted_rate
                << " limiter_changed_rate=" << (limiter_changed_rate ? "true" : "false")
                << " operation_count=" << inventory_operation_count
                << " reduction_depth=" << reduction_depth << '\n';
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(max_carrier == Real(0.0),
            "zero-flow SBM AMR qualification has a nonzero mass carrier");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(max_accepted_rate == Real(0.0),
            "zero-flow SBM AMR qualification has a nonzero accepted face rate");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!limiter_changed_rate,
            "zero-flow SBM AMR qualification changed an accepted face rate");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            correction_linf <= zero_flow_correction_tolerance,
            "zero-flow SBM AMR qualification produced a correction beyond "
            "its operation-count roundoff bound");
    } else if (solverChoice.sbm_amr_qualification_uniform_ratio) {
        // A uniform dry-air-relative spectrum stays uniform only as far as
        // its fluxes and the dry-air fluxes agree on every level and
        // interface.  The host's own uniform scalar, carried and refluxed by
        // the same dry-air fluxes, measures how far that is.
        const auto relative_spread = [&] (auto&& field_at_level, int component) {
            Real ratio_min = std::numeric_limits<Real>::max();
            Real ratio_max = std::numeric_limits<Real>::lowest();
            for (int lev = 0; lev <= finest_level; ++lev) {
                const MultiFab& field = field_at_level(lev);
                MultiFab ratio(field.boxArray(), field.DistributionMap(), 1, 0);
                MultiFab::Copy(ratio, field, component, 0, 1, 0);
                MultiFab::Divide(ratio, vars_new[lev][Vars::cons], Rho_comp, 0, 1, 0);
                ratio_min = std::min(ratio_min, ratio.min(0));
                ratio_max = std::max(ratio_max, ratio.max(0));
            }
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ratio_min > Real(0.0),
                "SBM uniform-ratio AMR qualification needs a positive uniform "
                "spectrum and host scalar (prob.rad_0 covering the domain)");
            return (ratio_max - ratio_min) / ratio_max;
        };
        const Real host_spread = relative_spread(
            [&] (int lev) -> const MultiFab& { return vars_new[lev][Vars::cons]; },
            RhoScalar_comp);
        // Where the host scalar stays uniform to roundoff the carriers agree
        // and so must the spectrum.  Where it does not (the compressible host
        // does not reflux dry-air density), the two nonuniform ratios are then
        // carried by different advection schemes; the spectrum must track the
        // host scalar's spread to within one percent.
        const bool host_uniform = host_spread <= rounding_factor;
        const Real bound = host_uniform
            ? rounding_factor
            : Real(1.01) * host_spread + rounding_factor;
        Real max_spread = Real(0.0);
        for (int component = 0; component < sbm_state_manager->layout().ncomp();
             ++component) {
            const Real spread = relative_spread(
                [&] (int lev) -> const MultiFab& {
                    return sbm_state_manager->new_state(lev);
                }, component);
            max_spread = std::max(max_spread, spread);
            Print().SetPrecision(17)
                << "SBM_AMR_QUAL_UNIFORM_RATIO step=" << nstep + 1
                << " component=" << component << " relative_spread=" << spread
                << " host_scalar_relative_spread=" << host_spread
                << " host_uniform=" << (host_uniform ? "true" : "false")
                << " bound=" << bound << '\n';
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            max_spread <= bound,
            "SBM AMR qualification did not keep a uniform dry-air-relative "
            "spectrum as uniform as the host's own scalar");
    } else {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!has_fine_level || correction_linf > Real(0.0),
            "SBM AMR qualification saw no nonzero accepted YA correction");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            sbm_transport->interval_limiter_active_groups() > 0 &&
            limiter_changed_rate,
            "SBM AMR qualification did not observe an active limiter changing an accepted rate");
    }
    const auto inventory = sbm_qualification_composite_inventory();
    AMREX_ALWAYS_ASSERT(inventory.size() == sbm_qualification_initial_inventory.size());
    Print() << "SBM_AMR_QUAL step=" << nstep + 1
            << " levels=" << finest_level + 1
            << " fine_substeps=" << (has_fine_level ? nsubsteps[1] : 0)
            << " leaf_cell_terms=" << leaf_cell_terms
            << " operation_count=" << inventory_operation_count
            << " reduction_depth=" << reduction_depth
            << " reflux_calls=" << reflux_calls
            << " reflux_syncs=" << reflux_syncs
            << " nonzero_reflux_linf=" << correction_linf
            << " limiter_active_groups="
            << sbm_transport->interval_limiter_active_groups()
            << " limiter_changed_rate=" << (limiter_changed_rate ? "true" : "false")
            << " zero_flow=" << (zero_flow ? "true" : "false");
    if (zero_flow) {
        Print() << " zero_flow_tolerance=" << zero_flow_correction_tolerance
                << " max_carrier=" << sbm_transport->interval_carrier_linf()
                << " max_accepted_rate=" << sbm_transport->interval_accepted_rate_linf();
    }
    Print() << '\n';
    for (std::size_t component = 0; component < inventory.size(); ++component) {
        const Real initial = sbm_qualification_initial_inventory[component];
        const Real final = inventory[component];
        const Real tolerance = rounding_factor * std::max(
            std::abs(initial) + std::abs(final), std::numeric_limits<Real>::min());
        const Real error = std::abs(final - initial);
        const Real scale = std::abs(initial) + std::abs(final);
        // Full precision: the rank-parity comparator parses these values.
        Print().SetPrecision(17)
            << "SBM_AMR_QUAL_INVENTORY component=" << component
            << " initial=" << initial << " final=" << final
            << " abs_error=" << error
            << " magnitude_scale=" << scale
            << " relative_error=" << (scale > Real(0.0) ? error / scale : Real(0.0))
            << " tolerance=" << tolerance
            << " leaf_cell_terms=" << leaf_cell_terms
            << " operation_count=" << inventory_operation_count
            << " reduction_depth=" << reduction_depth << '\n';
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(error <= tolerance,
            "SBM AMR qualification composite leaf inventory changed beyond its roundoff bound");
    }
}

void
ERF::sbm_qualification_check_regrid (int lev, Long cells, Long retained,
                                     Long prolonged, Long released,
                                     const std::vector<Real>& inventory_before,
                                     const std::vector<std::uint64_t>& phase_bytes)
{
    if (!solverChoice.sbm_amr_qualification_test) return;
    AMREX_ALWAYS_ASSERT(phase_bytes.size() == 4);
    const auto inventory_after = sbm_qualification_composite_inventory();
    AMREX_ALWAYS_ASSERT(inventory_after.size() == inventory_before.size());

    const auto coarse_cells =
        static_cast<std::uint64_t>(sbm_state_manager->new_state(0).boxArray().numPts());
    const auto fine_cells = finest_level == 1
        ? static_cast<std::uint64_t>(sbm_state_manager->new_state(1).boxArray().numPts())
        : std::uint64_t(0);
    const auto roundoff = composite_inventory_roundoff(
        coarse_cells, fine_cells, ref_ratio[0], true);
    const std::uint64_t leaf_cell_terms = roundoff.leaf_cell_terms;
    const Real rounding_factor = roundoff.rounding_factor;
    Real max_error = Real(0.0);
    Real max_error_over_bound = Real(0.0);
    for (std::size_t component = 0; component < inventory_after.size(); ++component) {
        const Real before = inventory_before[component];
        const Real after = inventory_after[component];
        const Real bound = rounding_factor * std::max(
            std::abs(before) + std::abs(after), std::numeric_limits<Real>::min());
        const Real error = std::abs(after - before);
        max_error = std::max(max_error, error);
        max_error_over_bound = std::max(max_error_over_bound, error / bound);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(error <= bound,
            "SBM regrid changed the composite leaf inventory beyond its roundoff bound");
    }
    Print().SetPrecision(17)
        << "SBM_AMR_QUAL_REGRID level=" << lev << " step=" << istep[0]
        << " levels=" << finest_level + 1
        << " cells=" << cells << " retained=" << retained
        << " prolonged=" << prolonged << " released=" << released
        << " leaf_cell_terms=" << leaf_cell_terms
        << " max_inventory_abs_error=" << max_error
        << " max_inventory_error_over_bound=" << max_error_over_bound << '\n';
    AllPrint() << "SBM_AMR_QUAL_REGRID_MEMORY rank=" << ParallelDescriptor::MyProc()
               << " level=" << lev << " step=" << istep[0]
               << " baseline_fab_bytes=" << phase_bytes[0]
               << " validation_peak_fab_bytes=" << phase_bytes[1]
               << " before_commit_fab_bytes=" << phase_bytes[2]
               << " after_commit_fab_bytes=" << phase_bytes[3]
               << " fab_bytes_high_water=" << TotalBytesAllocatedInFabsHWM() << '\n';
}
