#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "ERF_ProbCommon.H"

struct ProbParm : ProbParmDefaults {
    amrex::Real rho_0 = 1.0;
    amrex::Real T_0   = 300.0;
    amrex::Real u_0   = 0.0;
    amrex::Real v_0   = 0.0;
    amrex::Real w_0   = 0.0;

    int SLM_num_ref_inputs = 1;
    std::string SLM_ref_sounding_file;
    std::string SLM_ref_flux_file;
    std::string SLM_ref_sst_file;
}; // namespace ProbParm

class Problem : public ProblemBase
{
public:
    Problem(const amrex::Real* /*problo*/, const amrex::Real* /*probhi*/);

#include "Prob/ERF_InitConstantDensityHSE.H"

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
        const SolverChoice& sc,
        const int lev) override;

    void read_SLM_inputs();
    void read_SLM_sounding_file(const std::string &fname);
    void read_SLM_flux_file(const std::string &fname);
    void read_SLM_sst_file(const std::string &fname);

protected:
    std::string name() override { return "LandSurfaceModel_SLM"; }

private:
    ProbParm parms;

    std::vector<std::vector<amrex::Real>> sounding;
    std::vector<std::vector<amrex::Real>> fluxes;
    std::vector<std::vector<amrex::Real>> sst;
};

#endif
