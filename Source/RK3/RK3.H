#ifndef _RK3_H_
#define _RK3_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <ERF.H>

using namespace amrex;

void conservedToPrimitiveStag(MultiFab& prim_in, 
                              MultiFab& u_x, MultiFab& v_y, MultiFab& w_z,
                              MultiFab& cons_in,
                              MultiFab& cu_x, MultiFab& cu_y, MultiFab& cu_z);

void RK3_advance(MultiFab& cu, MultiFab& cu_new, 
                 MultiFab& cu_x, MultiFab& cu_y, MultiFab& cu_z, 
                 MultiFab& cunew_x, MultiFab& cunew_y, MultiFab& cunew_z, 
                 MultiFab& prim, 
                 MultiFab&  u_x, MultiFab&  v_y, MultiFab&  w_z, 
                 MultiFab&  unew_x, MultiFab&  vnew_y, MultiFab&  wnew_z, 
                 MultiFab& source,
                 MultiFab& eta, MultiFab& zeta, MultiFab& kappa,
                 std::array<MultiFab, AMREX_SPACEDIM>& faceflux,
                 std::array< MultiFab, 2 >& edgeflux_x,
                 std::array< MultiFab, 2 >& edgeflux_y,
                 std::array< MultiFab, 2 >& edgeflux_z,
                 std::array< MultiFab, AMREX_SPACEDIM>& cenflux,
                 const amrex::Geometry geom, const amrex::Real* dxp, const amrex::Real dt);

void calculateFluxStag(const MultiFab& cons_in, 
                       const MultiFab& cons_x, const MultiFab& cons_y, const MultiFab& cons_z, 
                       const MultiFab& prim_in, 
                       const MultiFab& u_x, const MultiFab& v_y, const MultiFab& w_z, 
                       const MultiFab& eta_in, const MultiFab& zeta_in, const MultiFab& kappa_in,
                       std::array<MultiFab, AMREX_SPACEDIM>& faceflux_in,
                       std::array< MultiFab, 2 >& edgeflux_x_in,
                       std::array< MultiFab, 2 >& edgeflux_y_in,
                       std::array< MultiFab, 2 >& edgeflux_z_in,
                       std::array< MultiFab, AMREX_SPACEDIM>& cenflux_in,
                       const amrex::Geometry geom,
                       const amrex::Real* dx, const amrex::Real dt);

void conservedToPrimitive(MultiFab& prim_in, const MultiFab& cons_in);

void primitiveToConserved(const MultiFab& prim, MultiFab& cons);

#endif
