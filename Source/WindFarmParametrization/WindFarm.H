#ifndef WINDFARM_H
#define WINDFARM_H

#include <DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

void advance_windfarm (int lev,
                       const amrex::Geometry& geom,
                       const amrex::Real& dt_advance,
                       amrex::MultiFab& cons_in,
                       amrex::MultiFab& U_old, amrex::MultiFab& V_old, amrex::MultiFab& W_old,
                       amrex::MultiFab& mf_vars_windfarm, const amrex::MultiFab& mf_Nturb,
                       SolverChoice& solver_choice);
#endif

