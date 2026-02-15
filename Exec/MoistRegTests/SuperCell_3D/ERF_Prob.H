#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"

#include "ERF_ProbCommon.H"
#include "ERF_EOS.H"
#include "ERF_IndexDefines.H"
#include "ERF_TileNoZ.H"

struct ProbParm : ProbParmDefaults {

  amrex::Real x_c = 0.0;
  amrex::Real y_c = 0.0;
  amrex::Real z_c = 1.5e3;

  amrex::Real x_r = 10.0e3;
  amrex::Real y_r = 10.0e3;
  amrex::Real z_r = 1.5e3;

  amrex::Real C_p = 1004.0;
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
