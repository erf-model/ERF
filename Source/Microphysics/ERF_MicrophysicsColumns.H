/*! @file ERF_MicrophysicsColumns.H
 *  \brief Full-column layout for the column microphysics schemes (Morrison, WSM6, WDM6)
 */

#ifndef ERF_MICROPHYSICSCOLUMNS_H
#define ERF_MICROPHYSICSCOLUMNS_H

#include <array>
#include <memory>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

/**
 * Morrison, WSM6 and WDM6 sediment each column inside one per-column kernel over the
 * k range of the box they are handed: the box top is treated as the top of the column
 * (no inflow) and the box bottom as its ground (outflow booked as surface precipitation).
 * That is only right when no column is split between boxes.
 *
 * This class records, for one level's grids, whether any column is split between boxes
 * stacked in z, and if so builds a "column layout": a BoxArray covering exactly the same
 * cells in which every contiguous vertical run of covered cells lies inside one box. A
 * scheme then copies its fields onto that layout, runs its unchanged column kernels, and
 * copies the result back. The kernels only read and write their own column, so this gives
 * the same arithmetic as a run whose grids were never split in z.
 *
 * A run that does not reach the domain bottom (e.g. a refined patch aloft) has no ground,
 * so the schemes accumulate surface precipitation only where klo is the domain bottom.
 */
class MicrophysicsColumns
{
public:
    using FabPtr = std::shared_ptr<amrex::MultiFab>;

    //! Build from a level's grids; call from the scheme's Init (which also runs on regrid)
    void define (const amrex::BoxArray& ba, const amrex::Box& domain);

    //! True when some column is split between boxes, so the fields must be gathered
    [[nodiscard]] bool needs_gather () const { return m_needs_gather; }

    //! True when every box spans the domain in z (required by the Fortran bridges)
    [[nodiscard]] bool boxes_span_domain () const { return m_boxes_span_domain; }

    //! Swap each non-null field for a copy on the column layout; the originals go to saved
    template <std::size_t N>
    void gather (std::array<FabPtr,N>& vars, std::array<FabPtr,N>& saved) const
    {
        for (std::size_t n = 0; n < N; ++n) {
            saved[n] = vars[n];
            if (vars[n]) { vars[n] = gather_one(*vars[n]); }
        }
    }

    //! Copy the column-layout fields back into the originals and restore them
    template <std::size_t N>
    void scatter (std::array<FabPtr,N>& vars, std::array<FabPtr,N>& saved) const
    {
        for (std::size_t n = 0; n < N; ++n) {
            if (saved[n]) { saved[n]->ParallelCopy(*vars[n], 0, 0, vars[n]->nComp()); }
            vars[n] = saved[n];
            saved[n].reset();
        }
    }

    //! Copy of a field (cell-centered or z-face; with its ghost cells) on the column layout
    [[nodiscard]] FabPtr gather_one (const amrex::MultiFab& src) const;

    //! Copy of a nodal field (e.g. z_phys_nd) on the column layout, or null for null input
    [[nodiscard]] std::unique_ptr<amrex::MultiFab> gather_nodal (const amrex::MultiFab* src) const;

    //! Copy of a 2D (z-flattened) integer field (e.g. lmask) on the column layout, or null
    [[nodiscard]] std::unique_ptr<amrex::iMultiFab> gather_2d (const amrex::iMultiFab* src) const;

private:
    bool m_needs_gather = false;
    bool m_boxes_span_domain = true;
    amrex::BoxArray m_col_ba;
    amrex::DistributionMapping m_col_dm;
};

#endif
