#include <AMReX_MultiFab.H>
#include <Src_headers.H>

using namespace amrex;

void CalcTurbPert (
  int level,
  const Geometry geom,
  const Box& bx,
  const Array4<Real>& cell_rhs,
  const Array4<const Real>& cell_data,
  TurbulentPerturbation& turbPert)
{
    // Domain cell size and real bounds
    //auto dx = geom.CellSizeArray();
    //auto ProbHiArr = geom.ProbHiArray();
    //auto ProbLoArr = geom.ProbLoArray();

    turbPert.calc_TurbPert_amplitude(level, bx, RhoTheta_comp, cell_rhs);
}
