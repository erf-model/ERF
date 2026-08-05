/**
 * \file ERF_EBRedistribute.H
 * \brief Declares embedded-boundary conservative redistribution helpers.
 */
#ifndef ERF_EBRedistribute_H_
#define ERF_EBRedistribute_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_BCRec.H>
#include <AMReX_MFIter.H>
#include <ERF_EB.H>

/**
 * \brief Apply EB state redistribution to a temporary tendency field.
 * \tparam EBFactType EB factory interface type, either AMReX native or ERF auxiliary.
 * \param ncomp Number of components to redistribute.
 * \param geom Geometry associated with the result field.
 * \param result Output redistributed field.
 * \param result_tmp Input temporary field before redistribution.
 * \param state State field used by redistribution slopes.
 * \param ebfact EB factory supplying flags and geometric moments.
 * \param bc Boundary conditions for the state components.
 * \param dt Time-step size used by the redistribution operator.
 * \param igrid Optional ERF integer-variable index for staggered grids; -1 means cell-centered.
 */
template <typename EBFactType>
void redistribute_term (int ncomp,
                        const amrex::Geometry& geom,
                        amrex::MultiFab& result,
                        amrex::MultiFab& result_tmp,
                        amrex::MultiFab const& state,
                        EBFactType const& ebfact,
                        amrex::BCRec const* bc,
                        double dt,
                        int const igrid = -1);

#endif
