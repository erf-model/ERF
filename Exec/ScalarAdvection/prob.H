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
pc_initdata(
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
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  // Define a point (xc,yc,zc) at the center of the domain 
  const amrex::Real xc = 0.5 * (prob_lo[0] + prob_hi[0]);
  const amrex::Real yc = 0.5 * (prob_lo[1] + prob_hi[1]);
  const amrex::Real zc = 0.5 * (prob_lo[2] + prob_hi[2]);

  const amrex::Real r  = (x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc);

  // Arbitrarily choose the radius of the bubble to be 0.05 times the length of the domain
  const amrex::Real r0 = 0.05 * (prob_hi[0] - prob_lo[0]);

  // Set the density 
  state(i, j, k, URHO) = ProbParm::rho_0;

  // Set the momentum based on the face-based velocities which we assume have already been set
  state(i, j, k, UMX) = 0.5 * (x_vel(i,j,k) + x_vel(i+1,j,k)) * state(i,j,k,URHO);
  state(i, j, k, UMY) = 0.5 * (y_vel(i,j,k) + y_vel(i,j+1,k)) * state(i,j,k,URHO);
  state(i, j, k, UMZ) = 0.5 * (z_vel(i,j,k) + z_vel(i,j,k+1)) * state(i,j,k,URHO);

  // (rho e)
  state(i, j, k, UEINT) = ProbParm::rhoe_0;
 
  // (rho E) = (rho e) + 1/2 rho vel^2
  state(i, j, k, UEDEN) = ProbParm::rhoe_0 + 0.5 * (1.0/state(i,j,k,URHO)) * (
                           state(i,j,k,UMX)*state(i,j,k,UMX) +
                           state(i,j,k,UMY)*state(i,j,k,UMY) +
                           state(i,j,k,UMZ)*state(i,j,k,UMZ) );

  // Initial temperature -- though we should make sure its consistent with (rho e) above
  state(i, j, k, UTEMP) = ProbParm::T_0;

  // Set scalar = A_0 in a ball of radius r0 and 0 elsewhere
  if (r < r0) 
     state(i, j, k, UFA) = ProbParm::A_0;
  else
     state(i, j, k, UFA) = 0.0;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_init_xvel(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& x_vel,
  amrex::Array4<amrex::Real> const& x_mom,
  amrex::GeometryData const& geomdata)
{
  // Note that this is called on a box of x-faces

  // Set the x-velocity
  x_vel(i, j, k) = 100.0;

  // Set the x-momentum
  x_mom(i, j, k) = 100.0;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_init_yvel(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& y_vel,
  amrex::Array4<amrex::Real> const& y_mom,
  amrex::GeometryData const& geomdata)
{
  // Note that this is called on a box of x-faces

  // Set the y-velocity
  y_vel(i, j, k) = 0.0;

  // Set the x-momentum
  y_mom(i, j, k) = 0.0;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_init_zvel(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& z_vel,
  amrex::Array4<amrex::Real> const& z_mom,
  amrex::GeometryData const& geomdata)
{
  // Note that this is called on a box of z-faces

  // Set the z-velocity
  z_vel(i, j, k) = 0.0;

  // Set the x-momentum
  z_mom(i, j, k) = 0.0;
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

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
