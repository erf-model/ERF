#include <AMReX_MultiFab.H>
#include <Src_headers.H>

using namespace amrex;

//#define DEBUG_PERT_MSG

void CalcTurbPert (
  const Geometry geom,
  const Box& bx,
  const Array4<Real>& cell_rhs,
  const Array4<const Real>& cell_data)
{
    // Domain cell size and real bounds
    //auto dx = geom.CellSizeArray();
    //auto ProbHiArr = geom.ProbHiArray();
    //auto ProbLoArr = geom.ProbLoArray();

    // Box creation through index values (for now)
    int offset_pertBox = 0;
    Box pertRegionBox(IntVect(0+offset_pertBox,0,0), IntVect(7+offset_pertBox,31,31)); // for visualization purpose
    BoxArray pertRegionBoxArray(pertRegionBox);                                        // creating a box array of this region
    pertRegionBoxArray.maxSize(8);                                                     // Further dividing box into box array

    #if DEBUG_PERT_MSG
    Print() << "Dividing created box into box array of element: "  << pertRegionBoxArray.size() << "\n";
    #endif

    for (int boxIdx = 0; boxIdx < pertRegionBoxArray.size(); boxIdx++) // Is this the best way to iteratre through the box array?
    {
        if (bx.contains(pertRegionBoxArray[boxIdx]))
        {
            #if DEBUG_PERT_MSG
            Print() << "bx: " << bx << " -- contains perturbation box #" << boxIdx << ": " << pertRegionBoxArray[boxIdx] << "\n";
            #endif

            // Parallel iterate through the box array and assign values
            ParallelFor(pertRegionBoxArray[boxIdx], [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                // Creating artificial values to populate the boxes
                //Real tmp_Src = (i + j*bx.size()[0] + k*bx.size()[0]*bx.size()[1]); // continuous count
                Real tmp_Src = (Real) boxIdx*1e5; // distinguish each box
     
                // Adding temperature source onto RHS
                cell_rhs(i, j, k, RhoTheta_comp) = tmp_Src;
            });
        }
    }
}
