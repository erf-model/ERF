#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_PlotFileUtil.H>

#include "ReadPlotFile.H"

using namespace amrex;

int main (int argc, char* argv[])
{
        std::cout << "Before reading data 0" << std::endl;
    Initialize(argc, argv);

    {
        if (argc != 4)
        {
            Print() << "Usage: " << argv[0]
                    << " <fine-mesh-plot-file> <coarse-mesh-plot-file>\n";
            Finalize();
            return 1;
        }

        std::string fine_plotfile   = argv[2];
        std::string coarse_plotfile = argv[3];

        MultiFab mf_cc_coarse;
        MultiFab mf_cc_fine;

        std::cout << "Before reading data" << std::endl;

        // Read fine plotfile
        PlotFileData pf_fine(fine_plotfile);
        ReadPlotFile("vars.txt", pf_fine, mf_cc_fine);

        std::cout << "Fine data read complete" << std::endl;

        // Read coarse plotfile
        PlotFileData pf_coarse(coarse_plotfile);
        ReadPlotFile("vars.txt", pf_coarse, mf_cc_coarse);

        std::cout << "Coarse data read complete" << std::endl;

        Array<int,AMREX_SPACEDIM> is_periodic
        {
            AMREX_D_DECL(0,0,0)
        };

        Geometry geom_fine(
            pf_fine.probDomain(0),
            RealBox(pf_fine.probLo(),
                    pf_fine.probHi()),
            pf_fine.coordSys(),
            is_periodic);

        Geometry geom_coarse(
            pf_coarse.probDomain(0),
            RealBox(pf_coarse.probLo(),
                    pf_coarse.probHi()),
            pf_coarse.coordSys(),
            is_periodic);

        MultiFab mf_cc_coarse_avg(mf_cc_coarse.boxArray(),
                                  mf_cc_coarse.DistributionMap(),
                                  mf_cc_fine.nComp(),
                                  0);

        AverageFineToCoarse(mf_cc_fine,
                            mf_cc_coarse_avg,
                            geom_fine,
                            geom_coarse);

        std::cout << "Interpolation from fine to coarse is complete" << std::endl;

        Vector<std::string> varnames = ReadVarNames("vars.txt");

        WriteSingleLevelPlotfile("plt_coarse_avg",
                                  mf_cc_coarse_avg,
                                  varnames,
                                  geom_coarse,
                                  0.0,
                                  0);

        std::string filename_custom = "coarse_data.bin";

        ApplyNeumannBCs(geom_coarse, mf_cc_coarse);

        WriteCustomDataFile(
            geom_coarse,
            mf_cc_coarse_avg,
            filename_custom);

        // Example:
        // WriteSingleLevelPlotfile("plt_final",
        //                          mf_cc_fine_from_coarse,
        //                          varnames,
        //                          geom_fine,
        //                          0.0,
        //                          0);

    }

    Finalize();

    return 0;
}

