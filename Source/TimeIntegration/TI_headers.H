#ifndef _INTEGRATION_H_
#define _INTEGRATION_H_

#include <AMReX_MultiFab.H>
#include "DataStruct.H"

#ifdef ERF_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

/**
 * Function for computing the buoyancy term to be used in the evolution
 * equation for the z-component of momentum in the slow integrator.  There
 * are three options for how buoyancy is computed (two are the same in the absence of moisture).
 */
void make_buoyancy (amrex::Vector<  amrex::MultiFab>& S_data,
                    const           amrex::MultiFab & S_prim,
                    amrex::MultiFab& buoyancy,
                    const amrex::Geometry geom,
                    const SolverChoice& solverChoice,
                    const amrex::MultiFab* r0,
                    const int& qstate_size);
#endif
