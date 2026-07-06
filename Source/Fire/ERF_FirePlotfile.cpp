#include "ERF_FirePlotfile.H"
#include "ERF_FirePlotfileCatalog.H"
#include "ERF_FireLayer.H"

#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

#include <iomanip>
#include <fstream>
#include <sstream>

using namespace amrex;

static void
write_fire_metadata_json(const std::string& plotfilename,
                         Real time, int step, int grid_ratio, int n_vars)
{
    if (!ParallelDescriptor::IOProcessor()) { return; }
    const std::string filename = plotfilename + "/FireMetadata.json";
    std::ofstream outfile(filename, std::ios::out | std::ios::trunc);
    if (!outfile.good()) { FileOpenFailed(filename); }
    outfile << "{\n";
    outfile << "  \"format_version\": 1,\n";
    outfile << "  \"time\": " << std::fixed << std::setprecision(15) << time << ",\n";
    outfile << "  \"step\": " << step << ",\n";
    outfile << "  \"grid_ratio\": " << grid_ratio << ",\n";
    outfile << "  \"n_variables\": " << n_vars << "\n";
    outfile << "}\n";
    if (!outfile.good()) { FileOpenFailed(filename); }
}

void
WriteFirePlotfile(const std::string& plotfile_prefix,
                  const FireLayer&   fire_layer,
                  Real               time,
                  int                step)
{
    const FireGrid& fg = fire_layer.get_fire_grid();

    Vector<std::string> varnames = fire_plotfile_var_names();
    int ncomp = fire_plotfile_ncomp();

    MultiFab mf(fg.ba, fg.dm, ncomp, 0);
    MultiFab::Copy(mf, *fire_layer.get_levelset(), 0, 0, 1, 0);
    MultiFab::Copy(mf, *fire_layer.get_ros(),      0, 1, 1, 0);
    MultiFab::Copy(mf, *fire_layer.get_wind_eff(), 0, 2, 2, 0);
    MultiFab::Copy(mf, *fire_layer.get_wind_ref(), 0, 4, 2, 0);
    MultiFab::Copy(mf, *fire_layer.get_slopes(),   0, 6, 2, 0);
    MultiFab::Copy(mf, *fire_layer.get_fuel_mc(),  0, 8, 3, 0);

    std::string plotfilename = Concatenate(plotfile_prefix, step, 5);

    // -----------------------------------------------------------------------
    // Step 1: IOProcessor creates BOTH the top-level plotfile directory AND
    //         the Level_0 subdirectory, then ALL ranks barrier before writing.
    // Using WriteSingleLevelPlotfile avoids the internal PreBuildDirectorHierarchy
    // race inside WriteMultiLevelPlotfile.
    // -----------------------------------------------------------------------
    if (ParallelDescriptor::IOProcessor()) {
        // Create top-level directory
        if (!amrex::UtilCreateDirectory(plotfilename, 0755)) {
            amrex::CreateDirectoryFailed(plotfilename);
        }
        // Create Level_0 subdirectory
        const std::string level_dir = plotfilename + "/Level_0";
        if (!amrex::UtilCreateDirectory(level_dir, 0755)) {
            amrex::CreateDirectoryFailed(level_dir);
        }
    }
    // All ranks wait until directories exist
    ParallelDescriptor::Barrier();

    // -----------------------------------------------------------------------
    // Step 2: Write using WriteSingleLevelPlotfile which does not call
    //         PreBuildDirectorHierarchy internally (directories already exist).
    // -----------------------------------------------------------------------
    WriteSingleLevelPlotfile(plotfilename, mf, varnames, fg.geom, time, step);

    if (ParallelDescriptor::IOProcessor()) {
        Print() << "[FIRE] Writing fire plotfile " << plotfilename << "\n";
    }

    // Barrier before JSON write so all ranks finish VisMF::Write first
    ParallelDescriptor::Barrier();
    write_fire_metadata_json(plotfilename, time, step, fg.C, ncomp);
}