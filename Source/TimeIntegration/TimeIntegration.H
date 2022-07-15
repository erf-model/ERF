#ifndef _INTEGRATION_H_
#define _INTEGRATION_H_

#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_InterpFaceRegister.H>
#include "DataStruct.H"
#include "IndexDefines.H"
#include "ABLMost.H"

namespace IntVar {
    enum {
        cons = 0,
        xmom,
        ymom,
        zmom,
        xflux,
        yflux,
        zflux
    };
}

// This is the slow RHS when doing multi-rate, and the only RHS when doing RK3
void erf_slow_rhs(int level,
                  amrex::Vector<amrex::MultiFab>& S_rhs,
                  const amrex::Vector<amrex::MultiFab>& S_data,
                  const amrex::MultiFab& S_prim,
                        amrex::Vector<amrex::MultiFab >& S_scratch,
                  const amrex::MultiFab& xvel,
                  const amrex::MultiFab& yvel,
                  const amrex::MultiFab& zvel,
                  amrex::MultiFab& source,
                  std::array< amrex::MultiFab, AMREX_SPACEDIM>&  advflux,
                  std::array< amrex::MultiFab, AMREX_SPACEDIM>& diffflux,
                  const amrex::Geometry geom,
                        amrex::InterpFaceRegister* ifr,
                  const SolverChoice& solverChoice,
                  std::unique_ptr<ABLMost>& most,
                  const amrex::Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d,
                  const amrex::MultiFab* z0,
                  const amrex::MultiFab* dJ,
                  const amrex::MultiFab* r0,
                  const amrex::MultiFab* p0,
                  const amrex::Real* dptr_dens_hse,
                  const amrex::Real* dptr_pres_hse,
                  const amrex::Real* dptr_rayleigh_tau,
                  const amrex::Real* dptr_rayleigh_ubar,
                  const amrex::Real* dptr_rayleigh_vbar,
                  const amrex::Real* dptr_rayleigh_thetabar,
                  const int rhs_vars=RHSVar::all);

void erf_fast_rhs (int level,
                   amrex::Vector<amrex::MultiFab >& S_rhs,
                   amrex::Vector<amrex::MultiFab >& S_slow_rhs,
                   amrex::Vector<amrex::MultiFab >& S_stage_data,
                   const amrex::MultiFab& S_stage_prim,
                   const amrex::Vector<amrex::MultiFab >& S_data,
                         amrex::Vector<amrex::MultiFab >& S_scratch,
                   std::array< amrex::MultiFab, AMREX_SPACEDIM>&  advflux,
                   const amrex::Geometry geom,
                   amrex::InterpFaceRegister* ifr,
                   const SolverChoice& solverChoice,
                   const amrex::MultiFab* z0,
                   const amrex::MultiFab* dJ,
                   const amrex::MultiFab* r0,
                   const amrex::MultiFab* p0,
                   const amrex::Real* dptr_dens_hse,
                   const amrex::Real* dptr_pres_hse,
                   const amrex::Real fast_dt, const amrex::Real invfac);
#endif
