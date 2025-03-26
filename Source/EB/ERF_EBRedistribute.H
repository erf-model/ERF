#ifndef ERF_EBRedistribute_H_
#define ERF_EBRedistribute_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_BCRec.H>
#include <AMReX_MFIter.H>

void redistribute_term (int ncomp,
                        const amrex::Geometry& geom,
                        amrex::MultiFab& result,
                        amrex::MultiFab& result_tmp,
                        amrex::MultiFab const& state,
                        amrex::EBFArrayBoxFactory const& ebfact,
                        amrex::BCRec const* bc,
                        amrex::Real const dt);

void redistribute_term (amrex::MFIter const& mfi,
                        int ncomp,
                        const amrex::Geometry& geom,
                        amrex::MultiFab& result,
                        amrex::MultiFab& result_tmp,
                        amrex::MultiFab const& state,
                        amrex::EBFArrayBoxFactory const& ebfact,
                        amrex::BCRec const* bc,
                        amrex::Real const dt);

#endif