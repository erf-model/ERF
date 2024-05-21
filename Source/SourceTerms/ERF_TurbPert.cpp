#include <AMReX_MultiFab.H>
#include <Src_headers.H>

using namespace amrex;

void CalcTurbPert (
  const Geometry geom,
  const Box& bx,
  const Array4<Real>& cell_rhs,
  const Array4<const Real>& cell_data)
{
    Print() << "Passed in Box: " << bx << "\n";

    // Domain cell size and real bounds
    auto dx = geom.CellSizeArray();
    auto ProbHiArr = geom.ProbHiArray();
    auto ProbLoArr = geom.ProbLoArray();

    // Box creation through index values (for now)
    Box pertRegionBox(IntVect(0,0,0), IntVect(15,31,31));
    BoxArray pertRegionBoxArray(pertRegionBox);  // creating a box array of this region
    pertRegionBoxArray.maxSize(16);
    Print() << "Dividing created box into box array of element: "  << pertRegionBoxArray.size() << "\n";

    // 1D array incremental values
    Real max_i = bx.size()[0];
    Real max_j = bx.size()[1];
    for (int boxIdx = 0; boxIdx < pertRegionBoxArray.size(); boxIdx++)
    {
        if (bx.contains(pertRegionBoxArray[boxIdx]))
        //if (bx.contains(pertRegionBox))
        {
            Print() << "My box: " << bx << " Contains sub box" << boxIdx << ": " << pertRegionBoxArray[boxIdx] << "\n";
            ParallelFor(pertRegionBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                // Creating incremental value pertRegionBoxArraysed on indexing
                //Real tmp_Src = (i + j*max_i + k*max_i*max_j);
                Real tmp_Src = (Real) boxIdx;
     
                // Adding temperature source onto RHS
                cell_rhs(i, j, k, RhoTheta_comp) = tmp_Src;
            });
        // Took these two lines out because it's okay if the region isn't in the right box
        // Depending on the if statement to do it's job 
        //} else {
        //    amrex::Abort("pertRegionBox not in bx!");
        }
    }
}






