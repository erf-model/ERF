#ifndef _RK3_H_
#define _RK3_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <ERF.H>

using namespace amrex;

// cons: cell-centered state
// U,V,W: face-based velocity
// source: cell-centered source term

//void ConservedToPrimitive(MultiFab& prim_out, 
//                          MultiFab& xvel_out, MultiFab& yvel_out, MultiFab& zvel_out,
//                          MultiFab& cons_in,
//                          MultiFab& xmom_in, MultiFab& ymom_in, MultiFab& zmom_in);
                          
void MomentumToVelocity(MultiFab& xvel_out, MultiFab& yvel_out, MultiFab& zvel_out,
                        MultiFab& cons_in,
                        MultiFab& xmom_in, MultiFab& ymom_in, MultiFab& zmom_in);

void VelocityToMomentum(MultiFab& xvel_in , MultiFab& yvel_in, MultiFab& zvel_in,
                        MultiFab& cons_in,
                        MultiFab& xmom_out, MultiFab& ymom_out, MultiFab& zmom_out);

void RK3_advance(MultiFab& cons_old, MultiFab& cons_new, 
                 MultiFab& U_old, MultiFab& V_old, MultiFab& W_old, 
                 MultiFab& U_new, MultiFab& V_new, MultiFab& W_new, 
                 MultiFab& source,
                 MultiFab& eta, MultiFab& zeta, MultiFab& kappa,
                 std::array<MultiFab, AMREX_SPACEDIM>& faceflux,
                 std::array< MultiFab, 2 >& edgeflux_x,
                 std::array< MultiFab, 2 >& edgeflux_y,
                 std::array< MultiFab, 2 >& edgeflux_z,
                 std::array< MultiFab, AMREX_SPACEDIM>& cenflux,
                 const amrex::Geometry geom, const amrex::Real* dxp, const amrex::Real dt);

void RK3_stage  (MultiFab& cons_old, MultiFab& cons_new, 
                 MultiFab& xmom_old, MultiFab& ymom_old, MultiFab& zmom_old, 
                 MultiFab& xmom_new, MultiFab& ymom_new, MultiFab& zmom_new, 
                 MultiFab& xvel    , MultiFab& yvel    , MultiFab& zvel    , 
                 MultiFab& prim    , MultiFab& source,
                 MultiFab& eta, MultiFab& zeta, MultiFab& kappa,
                 std::array<MultiFab, AMREX_SPACEDIM>& faceflux,
                 std::array< MultiFab, 2 >& edgeflux_x,
                 std::array< MultiFab, 2 >& edgeflux_y,
                 std::array< MultiFab, 2 >& edgeflux_z,
                 std::array< MultiFab, AMREX_SPACEDIM>& cenflux,
                 const amrex::Geometry geom, const amrex::Real* dxp, const amrex::Real dt);

void CalcAdvFlux (const MultiFab& cons_in, 
                  const MultiFab& cons_x, const MultiFab& cons_y, const MultiFab& cons_z, 
                  const MultiFab& u_x, const MultiFab& v_y, const MultiFab& w_z, 
                  std::array<MultiFab, AMREX_SPACEDIM>& faceflux_in,
                  std::array< MultiFab, 2 >& edgeflux_x_in,
                  std::array< MultiFab, 2 >& edgeflux_y_in,
                  std::array< MultiFab, 2 >& edgeflux_z_in,
                  std::array< MultiFab, AMREX_SPACEDIM>& cenflux_in,
                  const amrex::Geometry geom,
                  const amrex::Real* dx, const amrex::Real dt);

void CalcDiffFlux(const MultiFab& cons_in, 
                  const MultiFab& cons_x, const MultiFab& cons_y, const MultiFab& cons_z, 
                  const MultiFab& u_x, const MultiFab& v_y, const MultiFab& w_z, 
                  const MultiFab& eta_in, const MultiFab& zeta_in, const MultiFab& kappa_in,
                  std::array<MultiFab, AMREX_SPACEDIM>& faceflux_in,
                  std::array< MultiFab, 2 >& edgeflux_x_in,
                  std::array< MultiFab, 2 >& edgeflux_y_in,
                  std::array< MultiFab, 2 >& edgeflux_z_in,
                  std::array< MultiFab, AMREX_SPACEDIM>& cenflux_in,
                  const amrex::Geometry geom,
                  const amrex::Real* dx, const amrex::Real dt);

//void conservedToPrimitive(MultiFab& prim_in, const MultiFab& cons_in);

//void primitiveToConserved(const MultiFab& prim, MultiFab& cons);

#endif
