#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Random.H>

#include "IndexDefines.H"
#include "ERF_Constants.H"
#include "EOS.H"
#include "Tagging.H"
#include "ProblemDerive.H"

struct ProbParm {
  amrex::Real rho_0 = 0.0; // left density (g/cc)
  amrex::Real rhoe_0;
  amrex::Real Theta_0 = 1.0;
  amrex::Real A_0 = 1.0;

  amrex::Real U0 = 10.0;
  amrex::Real V0 = 0.0;
  amrex::Real W0 = 0.0;
  amrex::Real U0_Pert_Mag = 1.0;
  amrex::Real V0_Pert_Mag = 1.0;
  amrex::Real W0_Pert_Mag = 0.0;
}; // namespace ProbParm

extern ProbParm parms;

void
erf_init_hse(amrex::Vector<amrex::Real>& dens_hse,
             amrex::GeometryData const& geomdata);

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
