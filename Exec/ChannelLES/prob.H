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
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  // Define a point (xc,yc,zc) at the center of the domain 
  const amrex::Real xc = 0.5 * (prob_lo[0] + prob_hi[0]);
  const amrex::Real yc = 0.5 * (prob_lo[1] + prob_hi[1]);
  const amrex::Real zc = 0.5 * (prob_lo[2] + prob_hi[2]);

  const amrex::Real r  = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));

  // Arbitrarily choose the radius of the bubble to be 0.05 times the length of the domain

  // Set the density 
  state(i, j, k, Rho_comp) = ProbParm::rho_0;
 
  // Initial Rho0*Theta0
  state(i, j, k, RhoTheta_comp) = ProbParm::rho_0 * ProbParm::Theta_0;

  // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
  state(i, j, k, RhoScalar_comp) = ProbParm::A_0 * exp(-10.*r*r);
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

  // Set the x-velocity
  amrex::Real rand_double = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); // Between 0.0 and 1.0
  amrex::Real x_vel_prime = (rand_double*2.0 - 1.0)*ProbParm::U0_Pert_Mag;
//  if (i < 2 && j < 2 && k <2)
//      amrex::Print() << "random number in erf_init_xvel for (i, j, k) = (" << i << ", " << j << ", " << k <<"): "
//                     << rand_double << ", U' = " << x_vel_prime  << std::endl;
  x_vel(i, j, k) = ProbParm::U0 + x_vel_prime;
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

  // Set the y-velocity
  amrex::Real rand_double = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); // Between 0.0 and 1.0
  amrex::Real y_vel_prime = (rand_double*2.0 - 1.0)*ProbParm::V0_Pert_Mag;
//  if (i < 2 && j < 2 && k <2)
//      amrex::Print() << "random number in erf_init_yvel for (i, j, k) = (" << i << ", " << j << ", " << k <<"): "
//                     << rand_double << ", V' = " << y_vel_prime << std::endl;
  y_vel(i, j, k) = ProbParm::V0 + y_vel_prime;
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
  amrex::Real rand_double = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); // Between 0.0 and 1.0
  amrex::Real z_vel_prime = (rand_double*2.0 - 1.0)*ProbParm::W0_Pert_Mag;
//  if (i < 2 && j < 2 && k <2)
//      amrex::Print() << "random number in erf_init_zvel for (i, j, k) = (" << i << ", " << j << ", " << k <<"): "
//                     << rand_double << ", W' = " << z_vel_prime << std::endl;
  z_vel(i, j, k) = ProbParm::W0 + z_vel_prime;
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
