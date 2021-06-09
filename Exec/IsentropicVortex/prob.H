#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

#include "IndexDefines.H"
#include "Constants.H"
#include "EOS.H"
#include "Tagging.H"
#include "Transport.H"
#include "ProblemDerive.H"
#include "prob_parm.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
erf_vortex_Gaussian(
  amrex::Real x,
  amrex::Real y)
{
  // Evaluate Gaussian function
  const amrex::Real r2 = ((x-ProbParm::xc)*(x-ProbParm::xc) + (y-ProbParm::yc)*(y-ProbParm::yc))
                         / (ProbParm::R*ProbParm::R);
  return ProbParm::beta * std::exp(-r2/(2*ProbParm::sigma*ProbParm::sigma));
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
erf_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& x_vel,
  amrex::Array4<amrex::Real> const& y_vel,
  amrex::Array4<amrex::Real> const& z_vel,
  amrex::GeometryData const& geomdata)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
//  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0]; // cell center
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1]; // cell center

  // Calculate perturbation temperature
  const amrex::Real Omg = erf_vortex_Gaussian(x,y);
  const amrex::Real deltaT = -(ProbParm::gamma - 1.0)/(2.0*ProbParm::sigma*ProbParm::sigma) * Omg*Omg;
  // Set the density 
  const amrex::Real rho_norm = std::pow(1.0 + deltaT, ProbParm::inv_gm1);
  state(i, j, k, Density_comp) = rho_norm * ProbParm::rho_inf;
 
  // Initial _potential_ temperature
  const amrex::Real T = (1.0 + deltaT) * ProbParm::T_inf;
  const amrex::Real p = std::pow(rho_norm, Gamma) / Gamma  // isentropic relation
                        * ProbParm::rho_inf*ProbParm::a_inf*ProbParm::a_inf;
  state(i, j, k, RhoTheta_comp) = T * std::pow(p_0 / p, R_d/c_p); // T --> theta
  state(i, j, k, RhoTheta_comp) *= state(i, j, k, Density_comp);

  // Set scalar = 0 -- unused
  state(i, j, k, RhoScalar_comp) = 0.0;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
erf_init_xvel(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& x_vel,
  amrex::GeometryData const& geomdata)
{
  // Note that this is called on a box of x-faces
  const amrex::Real* prob_lo = geomdata.ProbLo();
//  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + i * dx[0]; // face center
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1]; // cell center
  const amrex::Real Omg = erf_vortex_Gaussian(x,y);

  // Set the x-velocity
  x_vel(i, j, k) = (ProbParm::M_inf * std::cos(ProbParm::alpha)
                 - (y - ProbParm::yc)/ProbParm::R * Omg) * ProbParm::a_inf;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
erf_init_yvel(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& y_vel,
  amrex::GeometryData const& geomdata)
{
  // Note that this is called on a box of y-faces
  const amrex::Real* prob_lo = geomdata.ProbLo();
//  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0]; // cell center
  const amrex::Real y = prob_lo[1] + j * dx[1]; // face center
  const amrex::Real Omg = erf_vortex_Gaussian(x,y);

  // Set the y-velocity
  y_vel(i, j, k) = (ProbParm::M_inf * std::sin(ProbParm::alpha)
                 + (x - ProbParm::xc)/ProbParm::R * Omg) * ProbParm::a_inf;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
erf_init_zvel(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& z_vel,
  amrex::GeometryData const& geomdata)
{
  // Note that this is called on a box of z-faces

  // Set the z-velocity
  z_vel(i, j, k) = 0.0;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata)
{
  for (int n = 0; n < NVAR; n++) {
    s_ext[n] = s_int[n];
  }
}

void erf_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
