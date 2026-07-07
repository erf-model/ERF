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

    bool has_spotting = (fire_layer.get_params().spotting.enable &&
                         fire_layer.get_albini_data() != nullptr);

    Vector<std::string> varnames = fire_plotfile_var_names(has_spotting);
    int ncomp = fire_plotfile_ncomp(has_spotting);

    MultiFab mf(fg.ba, fg.dm, ncomp, 0);
    
    // Component 0: fire_phi
    MultiFab::Copy(mf, *fire_layer.get_levelset(), 0, 0, 1, 0);
    // Component 1: fire_ros
    MultiFab::Copy(mf, *fire_layer.get_ros(),      0, 1, 1, 0);
    // Components 2-3: fire_wind_eff (u, v)
    MultiFab::Copy(mf, *fire_layer.get_wind_eff(), 0, 2, 2, 0);
    // Components 4-5: fire_wind_ref (u, v)
    MultiFab::Copy(mf, *fire_layer.get_wind_ref(), 0, 4, 2, 0);
    // Component 6: fire_wind_extract_z
    MultiFab::Copy(mf, *fire_layer.get_wind_extract_z(), 0, 6, 1, 0);
    // Components 7-8: fire_slopes (dz/dx, dz/dy)
    MultiFab::Copy(mf, *fire_layer.get_slopes(),   0, 7, 2, 0);
    // Components 9-11: fire_fuel_mc (M_1hr, M_10hr, M_100hr)
    MultiFab::Copy(mf, *fire_layer.get_fuel_mc(),  0, 9, 3, 0);
    // Component 12: fire_surface_temp_K
    MultiFab::Copy(mf, *fire_layer.get_surface_temp(), 0, 12, 1, 0);
    // Component 13: fire_surface_rh
    MultiFab::Copy(mf, *fire_layer.get_surface_rh(), 0, 13, 1, 0);
    // Component 14: fire_heat_flux
    MultiFab::Copy(mf, *fire_layer.get_heat_flux(), 0, 14, 1, 0);
    // Component 15: fire_fuel_load
    MultiFab::Copy(mf, *fire_layer.get_fuel_load(), 0, 15, 1, 0);
    // Component 16: fire_fireline_intensity
    MultiFab::Copy(mf, *fire_layer.get_fireline_intensity(), 0, 16, 1, 0);
    // Component 17: fire_flame_length
    MultiFab::Copy(mf, *fire_layer.get_flame_length(), 0, 17, 1, 0);
    // Component 18: fire_arrival_time
    MultiFab::Copy(mf, *fire_layer.get_arrival_time(), 0, 18, 1, 0);

    // Phase 8: Albini spotting diagnostics (components 19–22, when enabled)
    if (has_spotting) {
        MultiFab::Copy(mf, *fire_layer.get_albini_data(), 0, 19, 4, 0);
    }

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