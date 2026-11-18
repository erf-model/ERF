#ifndef ERF_SBM_INVENTORY_H_
#define ERF_SBM_INVENTORY_H_

#include <AMReX.H>
#include <AMReX_IntVect.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_REAL.H>

#include <cstdint>
#include <limits>

namespace erf_sbm {

struct InventoryRoundoff
{
    std::uint64_t leaf_cell_terms{0};
    std::uint64_t operation_count{0};
    int reduction_depth{0};
    amrex::Real rounding_factor{0};
};

/**
 * Operation-count bound on the difference of two composite leaf inventories,
 * shared by the per-step, regrid and boundary-budget conservation checks.
 * Each positive leaf term uses up to three multiplies and one local
 * accumulation; bound both inventories, the final subtraction, and a pairwise
 * MPI reduction with eight operations per leaf plus two reduction-tree levels
 * per depth. A regrid adds the prolongation of a new fine cell from its
 * parent, every leaf update between the two inventories re-forms each leaf,
 * and each physical-face term of the boundary ledger takes a product, a mask
 * and two sums.
 */
inline InventoryRoundoff
composite_inventory_roundoff (const std::uint64_t coarse_cells,
                              const std::uint64_t fine_cells,
                              const amrex::IntVect& ratio,
                              const bool with_prolongation,
                              const std::uint64_t leaf_updates = 0,
                              const std::uint64_t boundary_face_terms = 0)
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
    for (int ranks = amrex::ParallelDescriptor::NProcs(); ranks > 1;
         ranks = (ranks + 1) / 2) {
        ++roundoff.reduction_depth;
    }
    roundoff.operation_count =
        8U * roundoff.leaf_cell_terms * (1U + leaf_updates) + 6U * boundary_face_terms +
        2U * static_cast<std::uint64_t>(roundoff.reduction_depth) + 1U +
        (with_prolongation ? refinement_volume + 6U : 0U);
    const amrex::Real q = static_cast<amrex::Real>(roundoff.operation_count) *
                          std::numeric_limits<amrex::Real>::epsilon();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(q < amrex::Real(1.0),
        "SBM composite inventory check: operation count times machine epsilon "
        "reaches one, so this grid cannot be checked at this precision");
    roundoff.rounding_factor = q / (amrex::Real(1.0) - q);
    return roundoff;
}

} // namespace erf_sbm

#endif
