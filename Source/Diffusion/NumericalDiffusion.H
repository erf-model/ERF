#ifndef _NUMERICALDIFFUSION_H_
#define _NUMERICALDIFFUSION_H_

#include <AMReX.H>
#include <DataStruct.H>
#include <AMReX_MultiFab.H>

void NumericalDiffusion (const amrex::Box& bx,
                         const int n_start,
                         const int n_end,
                         const amrex::Real dt,
                         const SolverChoice& solverChoice,
                         const amrex::Array4<const amrex::Real>& data,
                         const amrex::Array4<      amrex::Real>& rhs,
                         const amrex::Array4<const amrex::Real>& mf_x,
                         const amrex::Array4<const amrex::Real>& mf_y,
                         const bool avg_mf_x_y,
                         const bool avg_mf_y_x);
#endif
