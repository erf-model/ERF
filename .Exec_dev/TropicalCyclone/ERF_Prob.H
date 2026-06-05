#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "ERF_ProbCommon.H"

struct ProbParm : ProbParmDefaults {
    amrex::Real KE_0   = amrex::Real(0.1);

    amrex::Real Xc_0  = amrex::Real(0.0);
    amrex::Real Yc_0  = amrex::Real(0.0);
    amrex::Real VMAX  = amrex::Real(15.0);   // Maximum horizontal velocity in vortex
    amrex::Real RMAX  = amrex::Real(100.0);  // Radius of maximum winds
    amrex::Real RZERO = amrex::Real(800.0);  // Radius of zero wind speed
    amrex::Real ZZERO = amrex::Real(2000.0); // Height of zero wind speed

    amrex::Real U_0 = amrex::Real(0.0);
    amrex::Real V_0 = amrex::Real(0.0);
    amrex::Real W_0 = amrex::Real(0.0);

}; // namespace ProbParm

class Problem : public ProblemBase
{
public:
    Problem();

#include "Prob/ERF_InitRayleighDamping.H" // provides WRF-style damping layer, with U_0, V_0, W_0, and T_0 as ref values
#include "Prob/ERF_InitCustomPert.H"

protected:

private:
    ProbParm parms;
};

#endif
