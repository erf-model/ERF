#include "ERF_Prob.H"
#include "ERF_EOS.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit(
    const amrex_real* /*problo*/,
    const amrex_real* /*probhi*/)
{
    return std::make_unique<Problem>();
}

Problem::Problem()
{
    // Parse params
    amrex::ParmParse pp("prob");
    pp.query("T_0", parms.T_0);
    pp.query("U_0", parms.U_0);
    pp.query("x_c", parms.x_c);
    pp.query("z_c", parms.z_c);
    pp.query("x_r", parms.x_r);
    pp.query("z_r", parms.z_r);
    pp.query("T_pert", parms.T_pert);

    init_base_parms(parms.rho_0, parms.T_0);
}
