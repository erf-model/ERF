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
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];

  // Define a point (xc,yc,zc) at (0.25,0.5,0.5) of domain span
  const amrex::Real xc = 0.25 * (prob_lo[0] + prob_hi[0]);
  const amrex::Real yc = 0.5 * (prob_lo[1] + prob_hi[1]);

  // 2D problem
  const amrex::Real r = std::sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc));

  // Arbitrary radius of disc
  const amrex::Real r0 = 0.2 * (prob_hi[0] - prob_lo[0]);

  // Set the density
  state(i, j, k, Rho_comp) = ProbParm::rho_0;
 
  // Initial potential temperature
  state(i, j, k, RhoTheta_comp) = ProbParm::rho_0 * ProbParm::T_0;

  // Set scalar = non-zero in a disk of radius r0 and 0 elsewhere
  state(i, j, k, RhoScalar_comp) = ProbParm::A_0 * 0.25 * (1 + cos(PI * std::min(r, r0) / r0));
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
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];

  // Set the x-velocity
  x_vel(i, j, k) = -y + 0.5;
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

  // Set the y-velocity
  y_vel(i, j, k) = x - 0.5;
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
