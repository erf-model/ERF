#include <vector>
#include <stdexcept>
#include <cassert>

#include <ERF.H>
#include <ERF_DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include "ERF_DA_EnKFSRF.H"

using namespace amrex;

void 
read_in_observations(MultiFab& y_obs)
{


}

void
Apply_H(const MultiFab& x_mf, MultiFab& y_mf)
{
    // Define y_mf with same BoxArray and DistributionMapping,
    // but only 2 components and same number of ghost cells
    y_mf.define(x_mf.boxArray(),
                x_mf.DistributionMap(),
                2,                  // number of components
                x_mf.nGrow());      // match ghost cells

    // Copy components 2 and 3 from x_mf into y_mf (as 0 and 1)
    MultiFab::Copy(y_mf, x_mf,
                   2,  // src component (3rd)
                   0,  // dest component
                   2,  // number of components
                   x_mf.nGrow());
}
