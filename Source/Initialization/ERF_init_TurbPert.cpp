#include <ERF.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

void
ERF::calcTurbulentPerturation (int lev)
{
#if 0
    for (int boxIdx = 0; boxIdx < turb_ba.size(); boxIdx++) // Is this the best way to iteratre through the box array?
    {
        if (bx.contains(turb_ba[boxIdx]))
        {
            Print() << "[Initialization/ERF_init_TurbPert.cpp] bx: " << bx <<
                       " -- contains perturbation box #" << boxIdx << ": " << turb_ba[boxIdx] << "\n";

            // Parallel iterate through the box array and assign values
            ParallelFor(turb_ba[boxIdx], [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                // Creating artificial values to populate the boxes
                //Real tmp_Src = (i + j*bx.size()[0] + k*bx.size()[0]*bx.size()[1]); // continuous count
                Real tmp_Src = (Real) boxIdx*1e5; // distinguish each box

                // Adding temperature source onto RHS
                cell_rhs(i, j, k, RhoTheta_comp) = tmp_Src;
            });
        }
    }
#endif
}
