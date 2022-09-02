#ifndef _INTEGRATION_H_
#define _INTEGRATION_H_

#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_InterpFaceRegister.H>
#include "DataStruct.H"
#include "IndexDefines.H"
#include "ABLMost.H"

// This is the slow RHS when doing multi-rate, and the only RHS when doing RK3
void erf_slow_rhs_pre(int level,
                      amrex::Vector<amrex::MultiFab>& S_rhs,
                      const amrex::Vector<amrex::MultiFab>& S_data,
                      const amrex::MultiFab& S_prim,
                            amrex::Vector<amrex::MultiFab >& S_scratch,
                      const amrex::MultiFab& xvel,
                      const amrex::MultiFab& yvel,
                      const amrex::MultiFab& zvel,
                      std::unique_ptr<amrex::MultiFab>& z_t,
                      const amrex::MultiFab& source,
                      const amrex::MultiFab& eddyDiffs,
                      std::array< amrex::MultiFab, AMREX_SPACEDIM>& diffflux,
                      const amrex::Geometry geom,
                            amrex::InterpFaceRegister* ifr,
                      const SolverChoice& solverChoice,
                      std::unique_ptr<ABLMost>& most,
                      const amrex::Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d,
                      std::unique_ptr<amrex::MultiFab>& z0,
                      std::unique_ptr<amrex::MultiFab>& dJ,
                      const amrex::MultiFab* r0,
                      const amrex::MultiFab* p0,
                      const amrex::Real* dptr_rayleigh_tau,
                      const amrex::Real* dptr_rayleigh_ubar,
                      const amrex::Real* dptr_rayleigh_vbar,
                      const amrex::Real* dptr_rayleigh_thetabar);

void erf_slow_rhs_post(int level,
                       amrex::Vector<amrex::MultiFab>& S_rhs,
                       const amrex::Vector<amrex::MultiFab>& S_data,
                       const amrex::MultiFab& S_prim,
                             amrex::Vector<amrex::MultiFab >& S_scratch,
                       const amrex::MultiFab& xvel,
                       const amrex::MultiFab& yvel,
                       const amrex::MultiFab& zvel,
                       std::unique_ptr<amrex::MultiFab>& z_t,
                       const amrex::MultiFab& source,
                       const amrex::MultiFab& eddyDiffs,
                       std::array< amrex::MultiFab, AMREX_SPACEDIM>& diffflux,
                       const amrex::Geometry geom,
                       const SolverChoice& solverChoice,
                       std::unique_ptr<ABLMost>& most,
                       const amrex::Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d,
                       std::unique_ptr<amrex::MultiFab>& z0,
                       std::unique_ptr<amrex::MultiFab>& dJ);

void erf_fast_rhs_N (int level, const amrex::Real time,
                     amrex::Vector<amrex::MultiFab >& S_slow_rhs,
                     amrex::Vector<amrex::MultiFab >& S_stage_data,
                     const amrex::MultiFab& S_stage_prim,
                     amrex::Vector<amrex::MultiFab >& S_data,
                     amrex::Vector<amrex::MultiFab >& S_scratch,
                     const amrex::Geometry geom,
                     amrex::InterpFaceRegister* ifr,
                     const SolverChoice& solverChoice,
                     const amrex::MultiFab* z_t_pert,
                     std::unique_ptr<amrex::MultiFab>& z0,
                     std::unique_ptr<amrex::MultiFab>& dJ,
                     const amrex::MultiFab* r0,
                     const amrex::MultiFab* p0,
                     const amrex::MultiFab* pi0,
                     const amrex::Real fast_dt, const amrex::Real invfac);

void erf_fast_rhs_T (int level, const amrex::Real time,
                     amrex::Vector<amrex::MultiFab >& S_slow_rhs,
                     amrex::Vector<amrex::MultiFab >& S_stage_data,
                     const amrex::MultiFab& S_stage_prim,
                     amrex::Vector<amrex::MultiFab >& S_data,
                     amrex::Vector<amrex::MultiFab >& S_scratch,
                     const amrex::Geometry geom,
                     amrex::InterpFaceRegister* ifr,
                     const SolverChoice& solverChoice,
                     const amrex::MultiFab* z_t_pert,
                     std::unique_ptr<amrex::MultiFab>& z0,
                     std::unique_ptr<amrex::MultiFab>& dJ,
                     const amrex::MultiFab* r0,
                     const amrex::MultiFab* p0,
                     const amrex::MultiFab* pi0,
                     const amrex::Real fast_dt, const amrex::Real invfac);
#endif
