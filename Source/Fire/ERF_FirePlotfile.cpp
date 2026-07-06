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

    // IOProcessor creates directories; all ranks barrier before writing.
    // We bypass WriteMultiLevelPlotfile entirely because it calls
    // PreBuildDirectorHierarchy with callBarrier=false (races on >1 rank)
    // and renames existing directories to .old.<pid>.
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

    // Write MultiFab data (all ranks participate)
    VisMF::Write(mf, plotfilename + "/Level_0/Cell");

    // Write AMReX Header (IOProcessor only)
    if (ParallelDescriptor::IOProcessor()) {
        const std::string header_path = plotfilename + "/Header";
        std::ofstream hfile(header_path, std::ios::out | std::ios::trunc);
        if (!hfile.good()) { FileOpenFailed(header_path); }

        hfile << "HyperCLaw-V1.1\n";
        hfile << ncomp << "\n";
        for (const auto& vn : varnames) { hfile << vn << "\n"; }
        hfile << "2\n";
        hfile << std::setprecision(17) << time << "\n";
        hfile << "0\n";

        const auto* plo  = fg.geom.ProbLo();
        const auto* phi_h = fg.geom.ProbHi();
        hfile << plo[0]  << " " << plo[1]  << "\n";
        hfile << phi_h[0] << " " << phi_h[1] << "\n";
        hfile << "\n";  // ref ratios — none for single level

        const Box& dom = fg.geom.Domain();
        hfile << "(" << dom.smallEnd(0) << "," << dom.smallEnd(1) << ",0) "
              << "(" << dom.bigEnd(0)   << "," << dom.bigEnd(1)   << ",0) "
              << "(0,0,0)\n";
        hfile << step << "\n";

        const auto* dx = fg.geom.CellSize();
        hfile << dx[0] << " " << dx[1] << " 1.0\n";
        hfile << "0\n";  // Cartesian
        hfile << "0\n";

        // Level 0 entry
        hfile << "0 " << fg.ba.size() << " "
              << std::setprecision(17) << time << "\n";
        hfile << step << "\n";
        for (int i = 0; i < fg.ba.size(); ++i) {
            const Box& b = fg.ba[i];
            RealBox rb(b, fg.geom.CellSize(), fg.geom.ProbLo());
            hfile << rb.lo(0) << " " << rb.hi(0) << "\n";
            hfile << rb.lo(1) << " " << rb.hi(1) << "\n";
        }
        hfile << "Level_0/Cell\n";

        if (!hfile.good()) { FileOpenFailed(header_path); }
        Print() << "[FIRE] Writing fire plotfile " << plotfilename << "\n";
    }

    ParallelDescriptor::Barrier();
    write_fire_metadata_json(plotfilename, time, step, fg.C, ncomp);
}