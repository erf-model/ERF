#ifndef _PROBCOMMON_H_
#define _PROBCOMMON_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

void amrex_probinit (
  const amrex_real* problo,
  const amrex_real* probhi);

#ifdef ERF_USE_TERRAIN
void
erf_init_dens_hse(amrex::MultiFab& dens_hse,
                  const amrex::MultiFab& z_phys_nd,
                  const amrex::MultiFab& z_phys_cc,
                  amrex::Geometry const& geom);

void
init_custom_prob(
  const amrex::Box& bx,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& x_vel,
  amrex::Array4<amrex::Real> const& y_vel,
  amrex::Array4<amrex::Real> const& z_vel,
  amrex::Array4<amrex::Real> const& r_hse,
  amrex::Array4<amrex::Real> const& p_hse,
  amrex::Array4<amrex::Real const> const& z_nd,
  amrex::Array4<amrex::Real const> const& z_cc,
  amrex::GeometryData const& geomdata);

#else

void
erf_init_dens_hse(amrex::Real* dens_hse_ptr,
                  amrex::Geometry const& geom,
                  const int ng_dens_hse);
void
init_custom_prob(
  const amrex::Box& bx,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& x_vel,
  amrex::Array4<amrex::Real> const& y_vel,
  amrex::Array4<amrex::Real> const& z_vel,
  amrex::GeometryData const& geomdata);
#endif

void
erf_init_rayleigh(amrex::Vector<amrex::Real>& tau,
                  amrex::Vector<amrex::Real>& ubar,
                  amrex::Vector<amrex::Real>& vbar,
                  amrex::Vector<amrex::Real>& thetabar,
                  amrex::Geometry      const& geom);
#endif

