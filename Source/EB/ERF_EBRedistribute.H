#ifndef ERF_EBRedistribute_H_
#define ERF_EBRedistribute_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_BCRec.H>
#include <AMReX_MFIter.H>
#include <ERF_EB.H>

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
