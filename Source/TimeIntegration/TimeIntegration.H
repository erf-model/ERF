#ifndef _INTEGRATION_H_
#define _INTEGRATION_H_

#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_InterpFaceRegister.H>
#include "DataStruct.H"
#include "IndexDefines.H"
#include "ABLMost.H"

using namespace amrex;

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

void erf_rhs(int level,
             Vector<MultiFab>& S_rhs,
             const Vector<MultiFab>& S_data,
             const MultiFab& S_prim,
             const MultiFab& xvel,
             const MultiFab& yvel,
             const MultiFab& zvel,
             MultiFab& source,
             std::array< MultiFab, AMREX_SPACEDIM>&  advflux,
             std::array< MultiFab, AMREX_SPACEDIM>& diffflux,
             const amrex::Geometry geom,
                   amrex::InterpFaceRegister* ifr,
             const SolverChoice& solverChoice,
             const ABLMost& most,
             const amrex::Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d,
             const amrex::Real* dptr_dens_hse,
             const amrex::Real* dptr_pres_hse,
             const amrex::Real* dptr_rayleigh_tau,
             const amrex::Real* dptr_rayleigh_ubar,
             const amrex::Real* dptr_rayleigh_vbar,
             const amrex::Real* dptr_rayleigh_thetabar);

void erf_fast_rhs(int level,
                  Vector<MultiFab >& S_rhs,
                  const Vector<MultiFab >& S_stage_data,
                  const Vector<MultiFab >& S_data,
                  std::array< MultiFab, AMREX_SPACEDIM>&  advflux,
                  const amrex::Geometry geom,
                        amrex::InterpFaceRegister* ifr,
                  const SolverChoice& solverChoice,
                  const amrex::Real* dptr_dens_hse, const amrex::Real* dptr_pres_hse);

#endif
