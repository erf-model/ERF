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
    // Phase 5: new plotfile fields
    MultiFab::Copy(mf, *fire_layer.get_heat_flux(), 0, 11, 1, 0);
    MultiFab::Copy(mf, *fire_layer.get_fuel_load(), 0, 12, 1, 0);
    MultiFab::Copy(mf, *fire_layer.get_fireline_intensity(), 0, 13, 1, 0);
    MultiFab::Copy(mf, *fire_layer.get_flame_length(), 0, 14, 1, 0);

    std::string plotfilename = Concatenate(plotfile_prefix, step, 5);

    // Step 1: IOProcessor creates directories; all ranks barrier before writing.
    // Bypasses WriteMultiLevelPlotfile which races on >1 rank and renames
    // existing directories to .old.<pid>.
    if (ParallelDescriptor::IOProcessor()) {
        if (!amrex::UtilCreateDirectory(plotfilename, 0755)) {
            amrex::CreateDirectoryFailed(plotfilename);
        }
        const std::string level_dir = plotfilename + "/Level_0";
        if (!amrex::UtilCreateDirectory(level_dir, 0755)) {
            amrex::CreateDirectoryFailed(level_dir);
        }
    }
    ParallelDescriptor::Barrier();

    // Step 2: All ranks write their owned fabs (collective call)
    VisMF::Write(mf, plotfilename + "/Level_0/Cell");

    // Step 3: IOProcessor writes the standard AMReX Header using the same
    // WriteGenericPlotfileHeader function used by the rest of ERF, so the
    // format is identical to any other AMReX plotfile and ParaView reads it.
    if (ParallelDescriptor::IOProcessor()) {
        const std::string header_path = plotfilename + "/Header";
        std::ofstream hfile;
        hfile.open(header_path.c_str(),
                   std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
        if (!hfile.good()) { FileOpenFailed(header_path); }

        // Use AMReX's canonical header writer — guarantees ParaView compatibility
        Vector<BoxArray>    ba_vec   = {fg.ba};
        Vector<std::string> var_vec  = varnames;
        Vector<Geometry>    geom_vec = {fg.geom};
        Vector<int>         steps    = {step};
        Vector<IntVect>     rr       = {};   // no refinement — single level

        WriteGenericPlotfileHeader(hfile, 1, ba_vec, var_vec,
                                   geom_vec, time, steps, rr);
        hfile << "Level_0/Cell\n";  // MultiFab path at Level 0
        if (!hfile.good()) { FileOpenFailed(header_path); }

        Print() << "[FIRE] Writing fire plotfile " << plotfilename << "\n";
    }

    ParallelDescriptor::Barrier();
    write_fire_metadata_json(plotfilename, time, step, fg.C, ncomp);
}