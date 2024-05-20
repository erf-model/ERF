#include <AMReX_MultiFab.H>
#include <Src_headers.H>

using namespace amrex;

void CalcTurbPert (
  const Box& bx,
  const Array4<Real>& cell_rhs,
  const Array4<const Real>& cell_data)
{
    // Domain cell size and real bounds
    // auto dx = geom.CellSizeArray();
    // auto ProbHiArr = geom.ProbHiArray();
    // auto ProbLoArr = geom.ProbLoArray();

    int n = RhoTheta_comp;  

    Box mybox(IntVect(1,0,0), IntVect(3,2,2));

    if (bx.contains(mybox)) {

    ParallelFor(mybox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        cell_rhs(i, j, k, n) += 0.0;
    });
    } else {
        amrex::Abort("mybox not in bx!");
    }
}






