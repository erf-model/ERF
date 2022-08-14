#ifndef _PROB_H_
#define _PROB_H_

#include "ERF_Constants.H"

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

#endif
