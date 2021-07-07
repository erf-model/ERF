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
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  // Set the density 
  state(i, j, k, Rho_comp) = ProbParm::rho_0;
 
  // Initial potential temperature (Actually rho*theta)
  const amrex::Real p = ProbParm::rho_0 * ProbParm::V_0*ProbParm::V_0*
                        (
                           1.0 / (Gamma * ProbParm::M_0 * ProbParm::M_0)
                        + (1.0 / 16.0) * (cos(2 * x) + cos(2 * y)) * (cos(2 * z) + 2)
                        );
  state(i, j, k, RhoTheta_comp) = getRhoThetagivenP(p);

  // Set scalar = 0 everywhere
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
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.0) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  // Set the x-velocity
  x_vel(i, j, k) = ProbParm::V_0 * sin(x) * cos(y) * cos(z);
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
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.0) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  // Set the y-velocity
  y_vel(i, j, k) = - ProbParm::V_0 * cos(x) * sin(y) * cos(z);
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
