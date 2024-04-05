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

void make_sources (int level, int nrk,
                   amrex::Real dt,
                   amrex::Vector<amrex::MultiFab>& S_data,
                   const amrex::MultiFab& S_prim,
                         amrex::MultiFab& source,
#if defined(ERF_USE_RRTMGP)
                   const amrex::MultiFab& qheating_rates,
#endif
                   const amrex::Geometry geom,
                   const SolverChoice& solverChoice,
                   std::unique_ptr<amrex::MultiFab>& mapfac_u,
                   std::unique_ptr<amrex::MultiFab>& mapfac_v,
                   const amrex::Real* dptr_rhotheta_src,
                   const amrex::Real* dptr_rhoqt_src,
                   const amrex::Real* dptr_wbar_sub,
                   const amrex::Vector<amrex::Real*> d_rayleigh_ptrs_at_lev);
#endif

#if defined(ERF_USE_NETCDF)
void
moist_set_rhs (
               const amrex::Box& tbx,
               const amrex::Array4<amrex::Real const>& old_cons,
               const amrex::Array4<amrex::Real const>& new_cons,
               const amrex::Array4<amrex::Real      >& cell_rhs,
               const amrex::Real& bdy_time_interval,
               const amrex::Real& start_bdy_time,
               const amrex::Real& new_stage_time,
               const amrex::Real& dt,
               int width, int set_width,
               const amrex::Box& domain,
               amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xlo,
               amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xhi,
               amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_ylo,
               amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_yhi);
#endif
