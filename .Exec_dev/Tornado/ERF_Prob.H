#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include <AMReX_REAL.H>

#include "ERF_ProbCommon.H"

struct ProbParm : ProbParmDefaults {
    amrex::Real rho_0 = amrex::Real(1.0);
    amrex::Real T_0   = amrex::Real(300.0);
    amrex::Real V_0   = amrex::Real(1.0);
    amrex::Real KE_0  = amrex::Real(0.1);

    amrex::Real KE_decay_height = amrex::Real(-1);
    amrex::Real KE_decay_order  = amrex::Real( 1);
}; // namespace ProbParm

class Problem : public ProblemBase
{
public:
    Problem();

#include "Prob/ERF_InitCustomPert.H"

private:
    ProbParm parms;
};

#endif
