#ifndef _TRANSPORT_H_
#define _TRANSPORT_H_

#include <cmath>

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX.H>
#include <AMReX_Gpu.H>

#include "TransportParams.H"

void transport_init();
void transport_close();

AMREX_GPU_DEVICE
void get_transport_coeffs(
  amrex::Box const& bx,
  amrex::Array4<const amrex::Real> const& T_in,
  amrex::Array4<const amrex::Real> const& Rho_in,
  amrex::Array4<amrex::Real> const& mu_out,
  amrex::Array4<amrex::Real> const& xi_out,
  amrex::Array4<amrex::Real> const& lam_out,
  const transport_params::transport_param_values& trparms);

using namespace transport_params;

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
transport(
  bool wtr_get_xi,
  bool wtr_get_mu,
  bool wtr_get_lam,
  amrex::Real& Tloc,
  amrex::Real& Rloc,
  amrex::Real& mu,
  amrex::Real& xi,
  amrex::Real& lam,
  const transport_params::transport_param_values& trparms)

{
  if (wtr_get_lam) {
    lam = trparms.const_conductivity;
  }

  if (wtr_get_mu) {
    mu = trparms.const_viscosity;
  }

  if (wtr_get_xi) {
    xi = trparms.const_bulk_viscosity;
  }
}

#endif
