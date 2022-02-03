//
// Created by pkjha on 2/3/22.
//

#ifndef ERF_WPS_INTERFACE_H
#define ERF_WPS_INTERFACE_H

#include <AMReX_MultiFab.H>
//#include "IO/NCInterface.H" // Uncomment when implementing NetCDF
// Can use NetCDF to Multifab conversion from AMR-wind
// Fore each MPI_rank, look for only it's (i, j, k) and fill them instead of filling entire Multifab by each MPI rank
// Parallel NetCDF to have each MPI rank read a portion ...particle to mesh redistribution
using namespace amrex;

class Metgrid_Output {
public:
    void readMetgridOutputFiles(const std::string &metgrid_files_loc);
    void convertNetCDFToPlotfiles(); // Output the read netCDF to plotfiles
    Metgrid_Output(); // We can probably pass teh path of the metgrid files here itself
private:
//std::vector<NCFile> metgrid_files; // Uncomment when implementing NetCDF

};

struct BoundaryData {
    MultiFab u, v, w, theta;
};

// Could have different multifab than BoundaryData
struct InitialData {
    MultiFab u, v, w, theta;
};

class ERF_IC_BC {
public:
    void readMetgridNetCDFFiles();
    void convertNetCDFBOundaryToMultifab();
    void convertNetCDFICToMultifab();
private:
    std::string metgrid_files_loc;
    Metgrid_Output metgridOutput;

    // Each element of the vector corresponds to a time stamp for mesa data, every 3hr or so
    std::vector<BoundaryData> x_lo;
    std::vector<BoundaryData> x_hi;
    std::vector<BoundaryData> y_lo;
    std::vector<BoundaryData> y_hi;

    // Initial Data (3-D flow-field)
    InitialData initialData;
};

#endif //ERF_WPS_INTERFACE_H