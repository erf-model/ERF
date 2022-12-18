#ifndef _PROB_H_
#define _PROB_H_

#include "EOS.H"
#include "AMReX_REAL.H"

struct ProbParm {
  // background conditions
  // if init_type != "" then these are perturbations
  amrex::Real T_0 = 0.0;
  amrex::Real U_0 = 0.0;
  amrex::Real V_0 = 0.0;

  // center of thermal perturbation
  amrex::Real x_c = 0.0;
  amrex::Real y_c = 0.0;
  amrex::Real z_c = 0.0;

  // radial extent of thermal perturbation
  amrex::Real x_r = 0.0;
  amrex::Real y_r = 0.0;
  amrex::Real z_r = 0.0;

  // perturbation temperature
  amrex::Real T_pert = -15.0;

  // physical constants that may be overwritten
  amrex::Real c_p = 1004.5;
}; // namespace ProbParm

extern ProbParm parms;

/*
 * Calculate perturbation to potential temperature and density at a
 * given (x,y,z). Note that pressure is not perturbed.
 */
AMREX_FORCE_INLINE
void
perturb_rho_theta (const amrex::Real x,
                   const amrex::Real y,
                   const amrex::Real z,
                   const amrex::Real p_hse,
                   const amrex::Real r_hse,
                   const ProbParm pp,
                   amrex::Real& rho,
                   amrex::Real& rhotheta)
{
    // Temperature that satisfies the EOS given the hydrostatically balanced (r,p)
    const amrex::Real Tbar_hse = p_hse / (R_d * r_hse);

    // Perturbation air temperature
    // - The bubble is either cylindrical (for 2-D problems, if two
    //   radial extents are specified) or an ellipsoid (if all three
    //   radial extents are specified).
    amrex::Real L = 0.0;
    if (pp.x_r > 0) L += std::pow((x - pp.x_c)/pp.x_r, 2);
    if (pp.y_r > 0) L += std::pow((y - pp.y_c)/pp.y_r, 2);
    if (pp.z_r > 0) L += std::pow((z - pp.z_c)/pp.z_r, 2);
    L = std::sqrt(L);
    amrex::Real dT;
    if (L > 1.0) {
        dT = 0.0;
    }
    else {
        dT = pp.T_pert * (std::cos(PI*L) + 1.0)/2.0;
    }

    // Note: theta_perturbed is theta PLUS perturbation in theta
    amrex::Real theta_perturbed = (Tbar_hse+dT)*std::pow(p_0/p_hse, R_d/pp.c_p);

    // Update state -- this version perturbs rho but not p
    rhotheta = getRhoThetagivenP(p_hse);
    rho = rhotheta / theta_perturbed;
}

#endif

