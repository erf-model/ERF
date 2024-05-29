#include <AMReX_MultiFab.H>
#include <Src_headers.H>

using namespace amrex;

void CalcTurbPert (
  const Geometry geom,
  const Box& bx,
  const Array4<Real>& cell_rhs,
  const Array4<const Real>& cell_data,
  const BoxArray& turb_ba)
{
#if 1
    // Domain cell size and real bounds
    //auto dx = geom.CellSizeArray();
    //auto ProbHiArr = geom.ProbHiArray();
    //auto ProbLoArr = geom.ProbLoArray();

    Print() << "\n";
    //Print() << "[Source/SourceTerms/ERF_TurbPert.cpp] turbPert_ba[0] contains "  << turb_ba << std::endl;
    Print() << "[Source/SourceTerms/ERF_TurbPert.cpp] Box being passed in reads: " << bx << "\n";

    for (int boxIdx = 0; boxIdx < turb_ba.size(); boxIdx++) // Is this the best way to iteratre through the box array?
    {
        if (bx.contains(turb_ba[boxIdx]))
        {
            Print() << "bx: " << bx << " -- contains perturbation box #" << boxIdx << ": " << turb_ba[boxIdx] << "\n";

            // Parallel iterate through the box array and assign values
            ParallelFor(turb_ba[boxIdx], [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                // Creating artificial values to populate the boxes
                //Real tmp_Src = (i + j*bx.size()[0] + k*bx.size()[0]*bx.size()[1]); // continuous count
                Real tmp_Src = (Real) (boxIdx+1.0)*1e5; // distinguish each box
     
                // Adding temperature source onto RHS
                cell_rhs(i, j, k, RhoTheta_comp) = tmp_Src;
            });
        }
    }
#endif
}
