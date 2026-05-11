#include "ERF_Prob.H"
#include "ERF_EOS.H"
#include "ERF_TerrainMetrics.H"

using namespace amrex;

#ifdef ERF_REMORA_FORCE_PROBINIT_LINK
// Force archive extraction of this TU when ERF is linked as a static library
// inside a parent coupled executable and amrex_probinit is weak.
void erf_probinit_link_anchor_func () noexcept {}
#endif

std::unique_ptr<ProblemBase>
amrex_probinit (const amrex_real* problo, const amrex_real* probhi)
{
    return std::make_unique<Problem>(problo, probhi);
}

Problem::Problem (const Real* /*problo*/, const Real* /*probhi*/)
{
    ParmParse pp_prob("prob");
    Real rho_0 =   1.0; int found_rho0 = pp_prob.query("rho_0", rho_0);
    Real p_inf        ; int found_p0   = pp_prob.query("p_inf", p_inf);

    Real   T_0 = 300.0; pp_prob.query("T_0", T_0);

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
    Array4<Real const> const& state,
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

    if (my_prob_name_ci == "abl") {
#include "Prob/ERF_InitCustomPert_ABL.H"
    } else if (my_prob_name_ci == "density current") {
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
    else if (my_prob_name_ci == "bubble") {
#include "Prob/ERF_InitCustomPert_Bubble.H"
    }
    else if (my_prob_name_ci == "multispecies bubble") {
#include "Prob/ERF_InitCustomPert_MultiSpeciesBubble.H"
    }
    else if  (my_prob_name_ci == "bomex") {
#include "Prob/ERF_InitCustomPert_Bomex.H"
    }
    else if  (my_prob_name_ci == "rico") {
#include "Prob/ERF_InitCustomPert_RICO.H"
    }
    else if  (my_prob_name_ci == "sdm_congestus3d") {
#include "Prob/ERF_InitCustomPert_SDMCongestus3D.H"
    }
    else if  (my_prob_name_ci == "squallline") {
#include "Prob/ERF_InitCustomPert_SquallLine.H"
    }
    else if  (my_prob_name_ci == "supercell") {
#include "Prob/ERF_InitCustomPert_SuperCell.H"
    }
     else if  (my_prob_name_ci == "data_assimilation_isv") {
#include "Prob/ERF_InitCustomPert_DataAssimilation_ISV.H"
    }
    else if  (my_prob_name_ci == "sinusoidalmassflux") {
#include "Prob/ERF_InitCustomPert_Bomex.H"
    }
    else {
        Print() << "Problem name" << " \"" <<  my_prob_name_ci << "\" "
                << "is not known, no state perturbations added. \n";
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

    if (my_prob_name_ci == "abl") {
#include "Prob/ERF_InitCustomPertVels_ABL.H"
    }
    else if ( (my_prob_name_ci == "couette flow") ||
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
              (my_prob_name_ci == "eb poiseuille"        ) ||
              (my_prob_name_ci == "multispecies bubble"  ) ) {
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
    if (my_prob_name_ci == "bubble") {
#include "Prob/ERF_InitCustomPertVels_ConstantU.H"
    }
    else if  (my_prob_name_ci == "bomex") {
#include "Prob/ERF_InitCustomPertVels_Bomex.H"
    }
    else if  (my_prob_name_ci == "rico") {
#include "Prob/ERF_InitCustomPertVels_RICO.H"
    }
    else if  (my_prob_name_ci == "sdm_congestus3d") {
#include "Prob/ERF_InitCustomPertVels_SDMCongestus3D.H"
    }
    else if ( (my_prob_name_ci == "squallline") ||
              (my_prob_name_ci == "supercell") ) {
#include "Prob/ERF_InitCustomPertVels_SquallLine.H"
    }
    else if (my_prob_name_ci == "userdefined") {
#include "Prob/ERF_InitCustomPertVels_UserDefined.H"
    }
     else if  (my_prob_name_ci == "data_assimilation_isv") {
#include "Prob/ERF_InitCustomPertVels_DataAssimilation_ISV.H"
    }
    else if  (my_prob_name_ci == "sinusoidalmassflux") {
#include "Prob/ERF_InitCustomPertVels_Bomex.H"
    }
    else {
        Print() << "Problem name" << " \"" <<  my_prob_name_ci << "\" "
                << "is not known, no velocity perturbations added. \n";
    }

    amrex::Gpu::streamSynchronize();
}

void
Problem::update_rhotheta_sources (const Real& time,
                                  amrex::MultiFab* src,
                                  const Geometry& geom,
                                  std::unique_ptr<MultiFab>& z_phys_cc)
{
    if (src->empty()) return;

    const int khi       = geom.Domain().bigEnd()[2];

    // If the z coordinate varies in time and or space, then the the height
    // needs to be calculated at each time step. Here, we assume that only
    // grid stretching exists.

    Vector<Real> zlevels;
    zlevels.resize(khi+1);

    Gpu::DeviceVector<Real> d_zlevels;
    d_zlevels.resize(khi+1);

    reduce_to_max_per_height(zlevels, z_phys_cc);
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, zlevels.begin(), zlevels.end(), d_zlevels.begin());

    const Real* d_zlevels_arr = d_zlevels.dataPtr();

    ParmParse pp_erf("erf");
    std::string my_prob_name; pp_erf.get("prob_name",my_prob_name);
    std::string my_prob_name_ci = amrex::toLower(my_prob_name);

    if (my_prob_name_ci == "bomex") {
#include "Prob/ERF_UpdateRhoThetaSources_Bomex.H"
    } else if (my_prob_name_ci == "rico") {
#include "Prob/ERF_UpdateRhoThetaSources_RICO.H"
    } else if  (my_prob_name_ci == "sdm_congestus3d") {
#include "Prob/ERF_UpdateRhoThetaSources_SDMCongestus3D.H"
    } else if  (my_prob_name_ci == "sinusoidalmassflux") {
#include "Prob/ERF_UpdateRhoThetaSources_SineMassFlux.H"
    }
}

void
Problem::update_rhoqt_sources (const Real& time,
                               amrex::MultiFab* qsrc,
                               const Geometry& geom,
                               std::unique_ptr<MultiFab>& z_phys_cc)
{
    if (qsrc->empty()) return;

    const int khi       = geom.Domain().bigEnd()[2];

    // If the z coordinate varies in time and or space, then the the height
    // needs to be calculated at each time step. Here, we assume that only
    // grid stretching exists.

    Vector<Real> zlevels;
    zlevels.resize(khi+1);

    Gpu::DeviceVector<Real> d_zlevels;
    d_zlevels.resize(khi+1);

    reduce_to_max_per_height(zlevels, z_phys_cc);
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, zlevels.begin(), zlevels.end(), d_zlevels.begin());

    const Real* d_zlevels_arr = d_zlevels.dataPtr();

    ParmParse pp_erf("erf");
    std::string my_prob_name; pp_erf.get("prob_name",my_prob_name);
    std::string my_prob_name_ci = amrex::toLower(my_prob_name);

    if  (my_prob_name_ci == "bomex") {
#include "Prob/ERF_UpdateRhoQtSources_Bomex.H"
    } else if (my_prob_name_ci == "rico") {
#include "Prob/ERF_UpdateRhoQtSources_RICO.H"
    } else if  (my_prob_name_ci == "sdm_congestus3d") {
#include "Prob/ERF_UpdateRhoQtSources_SDMCongestus3D.H"
    } else if  (my_prob_name_ci == "sinusoidalmassflux") {
#include "Prob/ERF_UpdateRhoQtSources_SineMassFlux.H"
    }
}

//=============================================================================
// USER-DEFINED FUNCTION
//=============================================================================
void
Problem::update_w_subsidence (const Real& time,
                              Vector<Real>& wbar,
                              Gpu::DeviceVector<Real>& d_wbar,
                              const amrex::MultiFab& state,
                              const Geometry& geom,
                              std::unique_ptr<MultiFab>& z_phys_nd)
{
    if (wbar.empty()) return;

    const int khi       = geom.Domain().bigEnd()[2];

    // If the z coordinate varies in time and or space, then the the height
    // needs to be calculated at each time step. Here, we assume that only
    // grid stretching exists.
    Vector<Real> zlevels;
    zlevels.resize(khi+2);
    reduce_to_max_per_height(zlevels, z_phys_nd);

    ParmParse pp_erf("erf");
    std::string my_prob_name; pp_erf.get("prob_name",my_prob_name);
    std::string my_prob_name_ci = amrex::toLower(my_prob_name);

    if  (my_prob_name_ci == "bomex") {
#include "Prob/ERF_UpdateWSubsidence_Bomex.H"
    } else if  (my_prob_name_ci == "rico") {
#include "Prob/ERF_UpdateWSubsidence_RICO.H"
    } else if  (my_prob_name_ci == "sinusoidalmassflux") {
#include "Prob/ERF_UpdateWSubsidence_SineMassFlux.H"
    }
}

//=============================================================================
// USER-DEFINED FUNCTION
//=============================================================================
void
Problem::update_geostrophic_profile (const Real& /*time*/,
                                     Vector<Real>& u_geos,
                                     Gpu::DeviceVector<Real>& d_u_geos,
                                     Vector<Real>& v_geos,
                                     Gpu::DeviceVector<Real>& d_v_geos,
                                     const Geometry& geom,
                                     std::unique_ptr<MultiFab>& z_phys_cc)
{
    if (u_geos.empty()) return;

    const int khi       = geom.Domain().bigEnd()[2];

    // If the z coordinate varies in time and or space, then the the height
    // needs to be calculated at each time step. Here, we assume that only
    // grid stretching exists.
    Vector<Real> zlevels;
    zlevels.resize(khi+1);
    reduce_to_max_per_height(zlevels, z_phys_cc);

    ParmParse pp_erf("erf");
    std::string my_prob_name; pp_erf.get("prob_name",my_prob_name);
    std::string my_prob_name_ci = amrex::toLower(my_prob_name);

    if  (my_prob_name_ci == "bomex") {
#include "Prob/ERF_UpdateGeostrophicProfile_Bomex.H"
    } else if  (my_prob_name_ci == "rico") {
#include "Prob/ERF_UpdateGeostrophicProfile_RICO.H"
    } else if  (my_prob_name_ci == "sinusoidalmassflux") {
#include "Prob/ERF_UpdateGeostrophicProfile_SineMassFlux.H"
    }

    // Copy from host version to device version
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, u_geos.begin(), u_geos.end(), d_u_geos.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, v_geos.begin(), v_geos.end(), d_v_geos.begin());
}
