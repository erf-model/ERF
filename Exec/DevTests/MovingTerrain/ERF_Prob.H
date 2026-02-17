#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "ERF_ProbCommon.H"

struct ProbParm : ProbParmDefaults {
  amrex::Real Ampl = 0.16;
  amrex::Real wavelength = 100.0;
  amrex::Real T_0  = 300.;
}; // namespace ProbParm

class Problem : public ProblemBase
{
public:
    Problem();

#include "Prob/ERF_InitDensityHSE.H"
#include "Prob/ERF_InitCustomPert.H"

    void erf_init_rayleigh (
        amrex::Vector<amrex::Vector<amrex::Real> >& rayleigh_ptrs,
        amrex::Geometry const& geom,
        std::unique_ptr<amrex::MultiFab>& z_phys_nd,
        amrex::Real zdamp) override;

protected:

private:
    ProbParm parms;
};

#endif
