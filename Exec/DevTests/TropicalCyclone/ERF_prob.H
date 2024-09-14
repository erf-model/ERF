#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "ERF_prob_common.H"

struct ProbParm : ProbParmDefaults {
    amrex::Real QKE_0   = 0.1;

    amrex::Real Xc_0  = 0.0;
    amrex::Real Yc_0  = 0.0;
    amrex::Real VMAX  = 15.0;   // Maximum horizontal velocity in vortex
    amrex::Real RMAX  = 100.0;  // Radius of maximum winds
    amrex::Real RZERO = 800.0;  // Radius of zero wind speed
    amrex::Real ZZERO = 2000.0; // Height of zero wind speed

    // rayleigh damping
    amrex::Real U_0 = 0.0;
    amrex::Real V_0 = 0.0;
    amrex::Real W_0 = 0.0;
    amrex::Real dampcoef = 0.2; // inverse time scale [1/s]
    amrex::Real zdamp = 500.0; // damping depth [m] from model top

}; // namespace ProbParm

class Problem : public ProblemBase
{
public:
    Problem();

//#include "Prob/ERF_init_constant_density_hse.H" // not needed if init_type=="input_sounding"
#include "Prob/ERF_init_rayleigh_damping.H" // provides WRF-style damping layer, with U_0, V_0, W_0, and T_0 as ref values

    void init_custom_pert (
        const amrex::Box&  bx,
        const amrex::Box& xbx,
        const amrex::Box& ybx,
        const amrex::Box& zbx,
        amrex::Array4<amrex::Real const> const& state,
        amrex::Array4<amrex::Real      > const& state_pert,
        amrex::Array4<amrex::Real      > const& x_vel_pert,
        amrex::Array4<amrex::Real      > const& y_vel_pert,
        amrex::Array4<amrex::Real      > const& z_vel_pert,
        amrex::Array4<amrex::Real      > const& r_hse,
        amrex::Array4<amrex::Real      > const& p_hse,
        amrex::Array4<amrex::Real const> const& z_nd,
        amrex::Array4<amrex::Real const> const& z_cc,
        amrex::GeometryData const& geomdata,
        amrex::Array4<amrex::Real const> const& mf_m,
        amrex::Array4<amrex::Real const> const& mf_u,
        amrex::Array4<amrex::Real const> const& mf_v,
        const SolverChoice& sc) override;

protected:
    std::string name() override { return "Tropical Cyclone"; }

private:
    ProbParm parms;
};

#endif
