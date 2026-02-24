#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFab.H>

#include "ReadPlotFile.H"

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        amrex::MultiFab mf;

        PlotFileData pf("plt000000");

        ReadPlotFile("vars.txt", pf, mf);

        amrex::Print() << "MultiFab has "
                       << mf.nComp() << " components\n";

        Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};

        Geometry geom(pf.probDomain(0),
              RealBox(pf.probLo(), pf.probHi()),
              pf.coordSys(),
              is_periodic);

        // Read variable names
        Vector<std::string> varnames = ReadVarNames("vars.txt");

        // Write plotfile
        WriteSingleLevelPlotfile("plt_new", mf, varnames, geom, 0.0, 0);

    }
    amrex::Finalize();
}

