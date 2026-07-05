#include "ERF_FirePlotfile.H"
#include "ERF_FirePlotfileCatalog.H"
#include "ERF_FireLayer.H"

#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

#include <iomanip>
#include <fstream>
#include <sstream>

using namespace amrex;

/**
 * @brief Write fire metadata JSON sidecar.
 *
 * @param[in] plotfilename  Path to the fire plotfile directory
 * @param[in] time          Simulation time [s]
 * @param[in] step          Timestep number
 * @param[in] grid_ratio    Fire grid refinement factor C
 * @param[in] n_vars        Number of output variables
 */
static void
write_fire_metadata_json(const std::string& plotfilename,
                         Real time,
                         int step,
                         int grid_ratio,
                         int n_vars)
{
    // Write JSON sidecar on I/O processor only
    if (!ParallelDescriptor::IOProcessor()) {
        return;
    }

    const std::string filename = plotfilename + "/FireMetadata.json";
    std::ofstream outfile(filename, std::ios::out | std::ios::trunc);
    if (!outfile.good()) {
        FileOpenFailed(filename);
    }

    outfile << "{\n";
    outfile << "  \"format_version\": 1,\n";
    outfile << "  \"time\": " << std::fixed << std::setprecision(15) << time << ",\n";
    outfile << "  \"step\": " << step << ",\n";
    outfile << "  \"grid_ratio\": " << grid_ratio << ",\n";
    outfile << "  \"n_variables\": " << n_vars << "\n";
    outfile << "}\n";

    if (!outfile.good()) {
        FileOpenFailed(filename);
    }
}

void
WriteFirePlotfile(const std::string& plotfile_prefix,
                  const FireLayer&   fire_layer,
                  Real               time,
                  int                step)
{
    // Get fire grid
    const FireGrid& fg = fire_layer.get_fire_grid();

    // Build variable names
    Vector<std::string> varnames = fire_plotfile_var_names();
    int ncomp = fire_plotfile_ncomp();

    // Allocate output MultiFab on fire grid
    MultiFab mf(fg.ba, fg.dm, ncomp, 0);

    // Copy components in canonical order:
    // 1. fire_phi (1 comp)
    MultiFab::Copy(mf, *fire_layer.get_levelset(), 0, 0, 1, 0);

    // 2. fire_ros (1 comp)
    MultiFab::Copy(mf, *fire_layer.get_ros(), 0, 1, 1, 0);

    // 3-4. fire_wind_eff (2 comps)
    MultiFab::Copy(mf, *fire_layer.get_wind_eff(), 0, 2, 2, 0);

    // 5-6. fire_wind_ref (2 comps)
    MultiFab::Copy(mf, *fire_layer.get_wind_ref(), 0, 4, 2, 0);

    // 7-8. fire_slopes (2 comps)
    MultiFab::Copy(mf, *fire_layer.get_slopes(), 0, 6, 2, 0);

    // 9-11. fire_fuel_mc (3 comps)
    MultiFab::Copy(mf, *fire_layer.get_fuel_mc(), 0, 8, 3, 0);

    // Generate plotfile name
    std::string plotfilename = Concatenate(plotfile_prefix, step, 5);

    // Write plotfile
    Vector<MultiFab> mf_vec = {mf};
    Vector<Geometry> geom_vec = {fg.geom};
    WriteMultiLevelPlotfile(plotfilename, 1, mf_vec, varnames, geom_vec, time, {step}, {});

    // Print message on I/O processor
    if (ParallelDescriptor::IOProcessor()) {
        Print() << "[FIRE] Writing fire plotfile " << plotfilename << "\n";
    }

    // Write JSON metadata sidecar on I/O processor
    write_fire_metadata_json(plotfilename, time, step, fg.C, ncomp);
}
