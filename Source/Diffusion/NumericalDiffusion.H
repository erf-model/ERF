#ifndef _NUMERICALDIFFUSION_H_
#define _NUMERICALDIFFUSION_H_

#include <AMReX.H>
#include <DataStruct.H>

void NumericalDiffusion (const amrex::Box& bx,
                         const int n_start,
                         const int n_end,
                         const amrex::Real dt,
                         const amrex::Array4<const amrex::Real>& data,
                         const amrex::Array4<      amrex::Real>& rhs,
                         const SolverChoice& solverChoice);
#endif
