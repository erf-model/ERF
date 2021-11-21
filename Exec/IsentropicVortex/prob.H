#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

#include "IndexDefines.H"
#include "ERF_Constants.H"
#include "EOS.H"
#include "Tagging.H"
#include "ProblemDerive.H"

struct ProbParm {
  amrex::Real p_inf = p_0; //freestream pressure [Pa]
  amrex::Real T_inf = 300.0; //freestream temperature [K]
  amrex::Real M_inf = 0.2; //freestream Mach number [-]
  amrex::Real alpha = 0.0; //inflow angle, 0 --> x-aligned [rad]
  amrex::Real gamma = Gamma; //specific heat ratio [-]
  amrex::Real beta = 0.01; //non-dimensional max perturbation strength [-]
  amrex::Real sigma = 1.0; //Gaussian standard deviation, i.e., spreading parameter [-]
  amrex::Real R = 2.0; //characteristic length scale for grid [m]
  amrex::Real xc = 0.5; //normalized x-location of vortex center [-]
  amrex::Real yc = 0.5; //normalized y-location of vortex center [-]
  // calculated quantiites
  amrex::Real rho_inf; //characteristic density [m/s]
  amrex::Real a_inf; //speed of sound [m/s]
  amrex::Real inv_gm1; //1/(gamma - 1) [-]
}; // namespace ProbParm

extern ProbParm parms;

AMREX_GPU_DEVICE
amrex::Real
erf_vortex_Gaussian(
  amrex::Real x,
  amrex::Real y,
  const ProbParm& parms);

void
erf_init_dens_hse(amrex::Vector<amrex::Real>& dens_hse);

void
erf_init_prob(
  const amrex::Box& bx,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& x_vel,
  amrex::Array4<amrex::Real> const& y_vel,
  amrex::Array4<amrex::Real> const& z_vel,
  amrex::GeometryData const& geomdata);

AMREX_GPU_DEVICE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata);

void erf_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
