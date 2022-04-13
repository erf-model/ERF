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
#ifdef ERF_USE_TERRAIN
             const MultiFab& z_phys_nd,
             const MultiFab& detJ_cc,
             const MultiFab& r0,
             const MultiFab& p0,
#else
             const amrex::Real* dptr_dens_hse,
             const amrex::Real* dptr_pres_hse,
#endif
             const amrex::Real* dptr_rayleigh_tau,
             const amrex::Real* dptr_rayleigh_ubar,
             const amrex::Real* dptr_rayleigh_vbar,
             const amrex::Real* dptr_rayleigh_thetabar);

void erf_implicit_fast_rhs (int level,
                            Vector<MultiFab >& S_rhs,
                            Vector<MultiFab >& S_slow_rhs,
                            Vector<MultiFab >& S_stage_data,
                            const MultiFab& S_stage_prim,
                            const Vector<MultiFab >& S_data,
                            const Vector<MultiFab >& S_data_old,
                            std::array< MultiFab, AMREX_SPACEDIM>&  advflux,
                            const amrex::Geometry geom,
                                  amrex::InterpFaceRegister* ifr,
                            const SolverChoice& solverChoice,
#ifdef ERF_USE_TERRAIN
                            const MultiFab& r0,
                            const MultiFab& p0,
#else
                            const amrex::Real* dptr_dens_hse, const amrex::Real* dptr_pres_hse,
#endif
                            const amrex::Real time, const amrex::Real fast_dt);
#endif
