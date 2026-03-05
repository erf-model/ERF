#include "ERF_Prob.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_Constants.H"
#include "ERF_EOS.H"
#include "ERF_HSEUtils.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit(const amrex_real* problo, const amrex_real* probhi)
{
    return std::make_unique<Problem>(problo, probhi);
}

Problem::Problem(const Real* /*problo*/, const Real* /*probhi*/)
{
    ParmParse pp("prob");
    Real rho_0 =   1.0; int found_rho0 = pp.query("rho_0", rho_0);
    Real p_inf        ; int found_p0   = pp.query("p_inf", p_inf);

    Real   T_0 = 300.0; pp.query("T_0", T_0);

    if (!found_rho0 && found_p0) {
        rho_0 = p_inf / (R_d * T_0);
    }

    init_base_parms(rho_0, T_0);
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_saturation_pressure (const Real T_b)
{
    return erf_esatw(T_b)*100.0;
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_relative_humidity (bool /*use_empirical*/)
{
    return 1.0;
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real vapor_mixing_ratio (const Real p_b, const Real T_b, const Real RH)
{
    Real p_s = compute_saturation_pressure(T_b);
    Real p_v = compute_vapor_pressure(p_s, RH);
    Real q_v = Rd_on_Rv*p_v/(p_b - p_v);
    return q_v;
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_dewpoint_temperature (const Real T_b, const Real RH)
{
    Real T_dp, gamma, T;
    T = T_b - 273.15;

    Real b = 18.678, c = 257.14, d = 234.5;
    gamma = log(RH*exp((b - T/d)*T/(c + T)));

    T_dp = c*gamma/(b - gamma);

    return T_dp;
}

void
Problem::init_custom_pert (
    const Box& bx,
    Array4<Real const> const& /*state*/,
    Array4<Real      > const& state_pert,
    Array4<Real      > const& r_hse,
    Array4<Real      > const& p_hse,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& /*z_cc*/,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    const SolverChoice& sc,
    const int /*lev*/)
{
    ParmParse pp_erf("erf");
    std::string my_prob_name; pp_erf.get("prob_name",my_prob_name);

    if (my_prob_name == "Bubble") {
#include "Prob/ERF_InitCustomPert_Bubble.H"
    }

    Gpu::streamSynchronize();
}

void
Problem::init_custom_pert_vels (
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    Array4<Real      > const& x_vel_pert,
    Array4<Real      > const& y_vel_pert,
    Array4<Real      > const& z_vel_pert,
    Array4<Real const> const& /*z_nd*/,
    GeometryData const& /*geomdata*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& /*sc*/,
    const int /*lev*/)
{
    ParmParse pp_erf("erf");
    std::string my_prob_name; pp_erf.get("prob_name",my_prob_name);

    if (my_prob_name == "Bubble") {
#include "Prob/ERF_InitCustomPertVels_ConstantU.H"
    }

    Gpu::streamSynchronize();
}
