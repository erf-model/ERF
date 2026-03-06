#include "ERF_Prob.H"
#include "ERF_EOS.H"
#include "ERF_TerrainMetrics.H"

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

void
Problem::init_custom_pert (
    const Box& bx,
    Array4<Real const> const& /*state*/,
    Array4<Real      > const& state_pert,
    Array4<Real      > const& r_hse,
    Array4<Real      > const& p_hse,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& z_cc,
    GeometryData const& geomdata,
    Array4<Real const> const&   mf_m,
    const SolverChoice& sc,
    const int lev)
{
    ParmParse pp_erf("erf");
    std::string my_prob_name; pp_erf.get("prob_name",my_prob_name);
    std::string my_prob_name_ci = amrex::toLower(my_prob_name);

    if (my_prob_name_ci == "density current") {
#include "Prob/ERF_InitCustomPert_DensityCurrent.H"
    }
    else if ( (my_prob_name_ci == "advecting isentropic vortex") ||
              (my_prob_name_ci == "stationary isentropic vortex") ) {
#include "Prob/ERF_InitCustomPert_IsentropicVortex.H"
    }
    else if ( (my_prob_name_ci == "particles over flat ground") ||
              (my_prob_name_ci == "particles over witch of agnesi hill") ) {
#include "Prob/ERF_InitCustomPert_ParticleTests.H"
    }
    else if (my_prob_name_ci == "scalar advection/diffusion") {
#include "Prob/ERF_InitCustomPert_ScalarAdvDiff.H"
    }
    else if (my_prob_name_ci == "stokes second problem") {
#include "Prob/ERF_InitCustomPert_StokesSecondProblem.H"
    }
    else if (my_prob_name_ci == "taylor-green vortex") {
#include "Prob/ERF_InitCustomPert_TaylorGreenVortex.H"
    }
    else if (my_prob_name_ci == "turbulent inflow") {
#include "Prob/ERF_InitCustomPert_TurbulentInflow.H"
    }
    else if (my_prob_name_ci == "moving terrain") {
#include "Prob/ERF_InitCustomPert_MovingTerrain.H"
    }
    else if (my_prob_name_ci == "eb poiseuille") {
#include "Prob/ERF_InitCustomPert_EBPoiseuille.H"
    }
    else if (my_prob_name_ci == "flow in a box") {
#include "Prob/ERF_InitCustomPert_FlowInABox.H"
    }
    else if (my_prob_name_ci == "userdefined") {
#include "Prob/ERF_InitCustomPert_UserDefined.H"
    }

    amrex::Gpu::streamSynchronize();
}

void
Problem::init_custom_pert_vels (
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    Array4<Real      > const& x_vel_pert,
    Array4<Real      > const& y_vel_pert,
    Array4<Real      > const& z_vel_pert,
    Array4<Real const> const& z_nd,
    GeometryData const& geomdata,
    Array4<Real const> const& mf_u,
    Array4<Real const> const& mf_v,
    const SolverChoice& sc,
    const int /*lev*/)
{
    ParmParse pp("erf");
    std::string my_prob_name; pp.get("prob_name",my_prob_name);
    std::string my_prob_name_ci = amrex::toLower(my_prob_name);

    if ( (my_prob_name_ci == "couette flow") ||
         (my_prob_name_ci == "poiseuille flow") ) {
#include "Prob/ERF_InitCustomPertVels_CouettePoiseuille.H"
    }
    else if (my_prob_name_ci == "ekman spiral") {
#include "Prob/ERF_InitCustomPertVels_EkmanSpiral.H"
    }
    else if ( (my_prob_name_ci == "advecting isentropic vortex") ||
              (my_prob_name_ci == "stationary isentropic vortex") ) {
#include "Prob/ERF_InitCustomPertVels_IsentropicVortex.H"
    }
    else if ( (my_prob_name_ci == "particles over flat ground") ||
              (my_prob_name_ci == "particles over witch of agnesi hill") ) {
#include "Prob/ERF_InitCustomPertVels_ParticleTests.H"
    }
    else if (my_prob_name_ci == "scalar advection/diffusion") {
#include "Prob/ERF_InitCustomPertVels_ScalarAdvDiff.H"
    }
    else if (my_prob_name_ci == "taylor-green vortex") {
#include "Prob/ERF_InitCustomPertVels_TaylorGreenVortex.H"
    }
    else if ( (my_prob_name_ci == "terrain - 2d cylinder") ||
              (my_prob_name_ci == "eb square cylinder"   ) ||
              (my_prob_name_ci == "eb poiseuille"   ) ) {
#include "Prob/ERF_InitCustomPertVels_ConstantU.H"
    }
    else if (my_prob_name_ci == "terrain - 3d hemisphere") {
#include "Prob/ERF_InitCustomPertVels_Terrain3DHemisphere.H"
    }
    else if (my_prob_name_ci == "turbulent inflow") {
#include "Prob/ERF_InitCustomPertVels_TurbulentInflow.H"
    }
    else if ( (my_prob_name_ci == "flow over witch of agnesi hill") ||
              (my_prob_name_ci == "flow over schar mountain") ) {
#include "Prob/ERF_InitCustomPertVels_WitchOfAgnesi.H"
    }
    else if (my_prob_name_ci == "moving terrain") {
#include "Prob/ERF_InitCustomPertVels_MovingTerrain.H"
    }
    else if (my_prob_name_ci == "userdefined") {
#include "Prob/ERF_InitCustomPertVels_UserDefined.H"
    }

    amrex::Gpu::streamSynchronize();
}
