/**
 * @file ERF_DustPlotfile.cpp
 * @brief Dust 2D plotfile writer on native dust grid.
 *
 * Follows ERF_FirePlotfile.cpp pattern:
 *   https://github.com/hgopalan/ERF/blob/ERF-Fire/Source/Fire/ERF_FirePlotfile.cpp
 *
 * Steps:
 *   1. IOProcessor creates plotfile/ and plotfile/Level_0/ directories.
 *   2. All ranks call VisMF::Write (collective).
 *   3. IOProcessor writes Header via WriteGenericPlotfileHeader.
 *   4. IOProcessor writes DustMetadata.json sidecar.
 *
 * The dust grid (m_dg.ba, m_dg.geom) is a 2D BoxArray independent of
 * the ERF AMR hierarchy, so WriteMultiLevelPlotfile is not used.
 */
#ifdef ERF_USE_DUST
#include "ERF_DustPlotfile.H"
#include "ERF_DustPlotfileCatalog.H"
#include "ERF_DustLayer.H"

#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <iomanip>
#include <fstream>
#include <string>

using namespace amrex;

static void write_dust_metadata_json(const std::string& plotfilename,
                                     Real time, int step, int grid_ratio, int n_vars)
{
    if (!ParallelDescriptor::IOProcessor()) return;
    const std::string filename = plotfilename + "/DustMetadata.json";
    std::ofstream out(filename, std::ios::out | std::ios::trunc);
    if (!out.good()) { FileOpenFailed(filename); }
    out << "{\n"
        << "  \"format_version\": 1,\n"
        << "  \"time\": " << std::fixed << std::setprecision(15) << time << ",\n"
        << "  \"step\": " << step << ",\n"
        << "  \"grid_ratio\": " << grid_ratio << ",\n"
        << "  \"n_variables\": " << n_vars << "\n"
        << "}\n";
    if (!out.good()) { FileOpenFailed(filename); }
}

void WriteDustPlotfile(const std::string& plotfile_prefix,
                       const DustLayer&   dust_layer,
                       Real               time,
                       int                step)
{
    const DustGrid& dg = dust_layer.get_dust_grid();

    Vector<std::string> varnames = dust_plotfile_var_names();
    int ncomp = dust_plotfile_ncomp();

    // Assemble all fields into one MultiFab.
    // Null fields (e.g. dust_surf_moist when no moisture scheme) are zeroed.
    MultiFab mf(dg.ba, dg.dm, ncomp, 0);
    mf.setVal(0.0);

    auto copy_if = [&](const MultiFab* src, int dst_comp) {
        if (src) MultiFab::Copy(mf, *src, 0, dst_comp, 1, 0);
    };

    copy_if(dust_layer.get_emission_flux(),    0);
    copy_if(dust_layer.get_ustar_in(),         1);
    copy_if(dust_layer.get_ustar_t(),          2);
    copy_if(dust_layer.get_deposition_rate(),  3);
    copy_if(dust_layer.get_conc_sfc(),         4);
    copy_if(dust_layer.get_surf_moist(),       5);
    copy_if(dust_layer.get_suppression(),      6);
    copy_if(dust_layer.get_retreat_flag(),     7);
    // Phase 17: EPA NAAQS fields
    copy_if(dust_layer.get_pm25(),             8);
    copy_if(dust_layer.get_pm10(),             9);
    copy_if(dust_layer.get_pm25_24h(),        10);
    copy_if(dust_layer.get_pm10_24h(),        11);
    copy_if(dust_layer.get_pm25_exceed(),     12);
    copy_if(dust_layer.get_pm10_exceed(),     13);
    // Phase 18: MSHA worker exposure fields
    copy_if(dust_layer.get_msha_dose(),       14);
    copy_if(dust_layer.get_msha_twa(),        15);
    copy_if(dust_layer.get_msha_exceed(),     16);
    copy_if(dust_layer.get_msha_shift_twa(),  17);
    // Phase 19: Lagrangian super-particle source-receptor attribution
#if defined(ERF_USE_PARTICLES)
    copy_if(dust_layer.get_source_map(),      18);
#else
    { MultiFab z(dg.ba,dg.dm,1,0); z.setVal(0.0);
     MultiFab::Copy(mf,z,0,18,1,0); }
#endif


    std::string plotfilename = Concatenate(plotfile_prefix, step, 5);

    // Step 1: IOProcessor creates directories; barrier before collective write.
    if (ParallelDescriptor::IOProcessor()) {
        if (!amrex::UtilCreateDirectory(plotfilename, 0755))
            amrex::CreateDirectoryFailed(plotfilename);
        const std::string level_dir = plotfilename + "/Level_0";
        if (!amrex::UtilCreateDirectory(level_dir, 0755))
            amrex::CreateDirectoryFailed(level_dir);
    }
    ParallelDescriptor::Barrier();

    // Step 2: All ranks write owned fabs (collective).
    VisMF::Write(mf, plotfilename + "/Level_0/Cell");

    // Step 3: IOProcessor writes standard AMReX Header.
    if (ParallelDescriptor::IOProcessor()) {
        const std::string header_path = plotfilename + "/Header";
        std::ofstream hfile(header_path.c_str(),
                            std::ofstream::out | std::ofstream::trunc |
                            std::ofstream::binary);
        if (!hfile.good()) { FileOpenFailed(header_path); }

        Vector<BoxArray>    ba_vec   = {dg.ba};
        Vector<std::string> var_vec  = varnames;
        Vector<Geometry>    geom_vec = {dg.geom};
        Vector<int>         steps    = {step};
        Vector<IntVect>     rr       = {};  // single level, no refinement

        WriteGenericPlotfileHeader(hfile, 1, ba_vec, var_vec,
                                   geom_vec, time, steps, rr);
        hfile << "Level_0/Cell\n";
        if (!hfile.good()) { FileOpenFailed(header_path); }

        Print() << "[DUST] Writing dust plotfile " << plotfilename << "\n";
    }
    ParallelDescriptor::Barrier();

    write_dust_metadata_json(plotfilename, time, step, dg.grid_ratio, ncomp);
}
#endif
