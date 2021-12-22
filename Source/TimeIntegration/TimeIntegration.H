#ifndef _INTEGRATION_H_
#define _INTEGRATION_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_InterpFaceRegister.H>
#include "DataStruct.H"
#include "IndexDefines.H"

using namespace amrex;

// cons: cell-centered state
// U,V,W: face-based velocity
// source: cell-centered source term

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

/**
 * Update the cell-centered state and face-based velocity using a time integrator
 *
 * \param[in]  S_old  state on cell centers
 * \param[in]  U_old  x-velocity on x-faces
 * \param[in]  V_old  y-velocity on y-faces
 * \param[in]  W_old  z-velocity on z-faces
 * \param[in]  source source term on cell centers
 *
 * \param[out] S_new  state on cell centers
 * \param[out] U_new  x-velocity on x-faces
 * \param[out] V_new  y-velocity on y-faces
 * \param[out] W_new  z-velocity on z-faces
 *
 */
void erf_advance(int level,
                 MultiFab& cons_old, MultiFab& cons_new,
                 MultiFab& U_old, MultiFab& V_old, MultiFab& W_old,
                 MultiFab& U_new, MultiFab& V_new, MultiFab& W_new,
                 MultiFab& xmom_crse, MultiFab& ymom_crse, MultiFab& zmom_crse,
                 MultiFab& source,
                 std::array<MultiFab, AMREX_SPACEDIM>& flux,
                 const amrex::Geometry crse_geom,
                 const amrex::Geometry fine_geom,
                 const amrex::IntVect ref_ratio,
                 const amrex::Real dt,
                 const amrex::Real time,
                       amrex::InterpFaceRegister* ifr,
                 const SolverChoice& solverChoice,
                 const amrex::Real* dptr_dens_hse,
                 const amrex::Real* dptr_pres_hse,
                 const amrex::Real* dptr_rayleigh_tau,
                 const amrex::Real* dptr_rayleigh_ubar,
                 const amrex::Real* dptr_rayleigh_vbar,
                 const amrex::Real* dptr_rayleigh_thetabar);

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
             const amrex::Geometry geom, const amrex::Real dt,
                   amrex::InterpFaceRegister* ifr,
             const SolverChoice& solverChoice,
             const bool lo_z_is_dirichlet,
             const bool hi_z_is_dirichlet,
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
