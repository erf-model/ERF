#ifndef _PROB_H_
#define _PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "prob_common.H"
#include "EOS.H"

struct ProbParm : ProbParmDefaults {
    // background conditions
    // if init_type != "" then these are perturbations and should be 0
    amrex::Real T_0 = 300.0;
    amrex::Real U_0 = 0.0;
    amrex::Real V_0 = 0.0;
    amrex::Real W_0 = 0.0;

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
    bool T_pert_is_airtemp = true; // T_pert input is air temperature
    bool perturb_rho = true; // not rho*theta (i.e., p is constant); otherwise perturb rho*theta

    // rayleigh damping
    amrex::Real dampcoef = 0.0;    // inverse time scale [1/s]
    amrex::Real zdamp    = 5000.0; // damping depth [m] from model top

    // Moist bubble params
    bool do_moist_bubble = false;
    amrex::Real theta_pert  = 2.0;
    amrex::Real qt_init     = 0.02;
    amrex::Real eq_pot_temp = 320.0;
}; // namespace ProbParm

class Problem : public ProblemBase
{
public:
    Problem();

#include "Prob/init_density_hse_dry.H"

    void init_custom_pert (
        const amrex::Box&  bx,
        const amrex::Box& xbx,
        const amrex::Box& ybx,
        const amrex::Box& zbx,
        amrex::Array4<amrex::Real      > const& state,
        amrex::Array4<amrex::Real      > const& x_vel,
        amrex::Array4<amrex::Real      > const& y_vel,
        amrex::Array4<amrex::Real      > const& z_vel,
        amrex::Array4<amrex::Real      > const& r_hse,
        amrex::Array4<amrex::Real      > const& p_hse,
        amrex::Array4<amrex::Real const> const& z_nd,
        amrex::Array4<amrex::Real const> const& z_cc,
        amrex::GeometryData const& geomdata,
        amrex::Array4<amrex::Real const> const& mf_m,
        amrex::Array4<amrex::Real const> const& mf_u,
        amrex::Array4<amrex::Real const> const& mf_v,
        const SolverChoice& sc) override;

    void erf_init_dens_hse_moist (amrex::MultiFab& rho_hse,
                                  std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                                  amrex::Geometry const& geom) override;

    void init_custom_terrain (
        const amrex::Geometry& geom,
        amrex::MultiFab& z_phys_nd,
        const amrex::Real& time) override;

    amrex::Real compute_theta (const amrex::Real T_b, const amrex::Real p_b);

    amrex::Real compute_temperature (const amrex::Real p_b);

    amrex::Real compute_F_for_temp(const amrex::Real T_b, const amrex::Real p_b);

    amrex::Real compute_p_k (amrex::Real& p_k,
                             const amrex::Real p_k_minus_1,
                             amrex::Real& theta_k,
                             amrex::Real& rho_k,
                             amrex::Real& q_v_k,
                             amrex::Real& T_dp,
                             amrex::Real& T_b,
                             const amrex::Real dz,
                             const amrex::Real rho_k_minus_1);

    amrex::Real compute_F (const amrex::Real& p_k,
                           const amrex::Real& p_k_minus_1,
                           amrex::Real &theta_k,
                           amrex::Real& rho_k,
                           amrex::Real& q_v_k,
                           amrex::Real& T_dp,
                           amrex::Real& T_b,
                           const amrex::Real& dz,
                           const amrex::Real& rho_k_minus_1);

    void compute_rho (const amrex::Real& pressure,
                      amrex::Real &theta,
                      amrex::Real& rho,
                      amrex::Real& q_v,
                      amrex::Real& T_dp,
                      amrex::Real& T_b);

    void init_isentropic_hse_no_terrain(amrex::Real *theta,
                                        amrex::Real* r,
                                        amrex::Real* p,
                                        amrex::Real *q_v,
                                        const amrex::Real& dz,
                                        const int& khi);

#include "Prob/init_rayleigh_damping.H"

protected:
    std::string name() override { return "Bubble"; }

private:
    ProbParm parms;
};

/*
 * Calculate perturbation to potential temperature or density at a
 * given (x,y,z). If the perturb_rho flag is set, then rho is perturbed,
 * and pressure is not (rho*theta is constant); otherwise, rho is
 * constant and rho*theta is perturbed
 *
 * Note: the input p_hse and r_hse should be in HSE given an appropriate
 * background flow initialization.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
perturb_rho_theta (const amrex::Real x,
                   const amrex::Real y,
                   const amrex::Real z,
                   const amrex::Real p_hse,
                   const amrex::Real r_hse,
                   const ProbParm pp,
                   const amrex::Real rdOcp,
                   amrex::Real& rho,
                   amrex::Real& rhotheta)
{
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
        dT = pp.T_pert * std::pow(cos(PI*L/2.0),2);
    }

    // Temperature that satisfies the EOS given the hydrostatically balanced (r,p)
    const amrex::Real Tbar_hse = p_hse / (R_d * r_hse);

    // Note: theta_perturbed is theta PLUS perturbation in theta
    amrex::Real theta_perturbed;
    if (pp.T_pert_is_airtemp) {
        // dT is air temperature
        theta_perturbed = (Tbar_hse + dT)*std::pow(p_0/p_hse, rdOcp);
    } else {
        // dT is potential temperature
        theta_perturbed = Tbar_hse*std::pow(p_0/p_hse, rdOcp) + dT;
    }

    if (pp.perturb_rho)
    {
        // this version perturbs rho but not p (i.e., rho*theta)
        // - hydrostatic rebalance is needed (TODO: is this true?)
        // - this is the approach taken in the density current problem
        rhotheta = 0.0; // i.e., hydrostatically balanced pressure stays const
        rho = getRhoThetagivenP(p_hse) / theta_perturbed - r_hse;
    }
    else
    {
        // this version perturbs rho*theta (i.e., p) but not rho
        rho = 0.0; // i.e., hydrostatically balanced density stays const
        rhotheta = r_hse * theta_perturbed - getRhoThetagivenP(p_hse);
    }
}

#endif
