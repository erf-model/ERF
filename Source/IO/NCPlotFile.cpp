#include <iomanip>
#include <iostream>
#include <string>
#include <ctime>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_Utility.H>
#include <AMReX_buildInfo.H>
#include <AMReX_ParmParse.H>

#include "ERF.H"
#include "NCInterface.H"
#include "NCPlotFile.H"
#include "IndexDefines.H"

namespace {
 std::string nc_state_filename = "Plot_State_MF.nc";
} // namespace

using namespace amrex;

void
ERF::writeNCPlotFile(const std::string& dir, const Vector<const MultiFab*> &plotMF,
                     const Vector<std::string> &plot_var_names,
                     const Vector<int> level_steps, const Real time) const
{
     //
     // TODO: Right now this appears to be hard-wired for single-level so we'll leave it that way
     //
     int lev = 0;

     // get the processor number
     int iproc = amrex::ParallelContext::MyProcAll();

     // number of grid points at this block
     int num_pts;

     // number of cells in each block
     std::vector<int> n_cells;

     // set the full IO path for NetCDF output
     static const std::string BaseName = "/NetCDF";
     char buf[64];
     sprintf(buf, "Level_%d", lev);
     std::string LevelStr = buf;

     std::cout << "Writing NetCDF plot file, iproc= " << iproc << '\n';

     std::string FullPath = dir;
     if (!FullPath.empty() && FullPath[FullPath.size() - 1] != '/') FullPath += '/';
     FullPath += LevelStr;
     FullPath += BaseName;

     if (!amrex::UtilCreateDirectory(FullPath, 0755))
         amrex::CreateDirectoryFailed(FullPath);

     // open netcdf file to write data
     auto ncf = ncutils::NCFile::create_par(FullPath+'/'+nc_state_filename, NC_NETCDF4 | NC_MPIIO,
                                            amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);


     int nblocks = grids[lev].size();
     int flev    = finest_level;

     n_cells.push_back(amrex::length(grids[lev][0]).x);
     n_cells.push_back(amrex::length(grids[lev][0]).y);
     n_cells.push_back(amrex::length(grids[lev][0]).z);
     num_pts = (amrex::length(grids[lev][0]).x * amrex::length(grids[lev][0]).y * amrex::length(grids[lev][0]).z);

     int n_data_items = plotMF[lev]->nComp();

     const std::string nt_name   = "num_time_steps";
     const std::string ndim_name = "num_geo_dimensions";
     const std::string np_name   = "num_points_per_block";
     const std::string nb_name   = "num_blocks";
     const std::string nx_name   = "NX";
     const std::string ny_name   = "NY";
     const std::string nz_name   = "NZ";
     const std::string flev_name = "FINEST_LEVEL";

     ncf.enter_def_mode();
     ncf.put_attr("title", "ERF NetCDF Plot data output");
     ncf.def_dim(nt_name,   NC_UNLIMITED);
     ncf.def_dim(ndim_name, AMREX_SPACEDIM);
     ncf.def_dim(np_name,   num_pts);
     ncf.def_dim(nb_name,   nblocks);
     ncf.def_dim(flev_name, flev);

     ncf.def_dim(nx_name,   n_cells[0]);
     ncf.def_dim(ny_name,   n_cells[1]);
     ncf.def_dim(nz_name,   n_cells[2]);

     ncf.def_var("probLo"  ,   NC_FLOAT,  {ndim_name});
     ncf.def_var("probHi"  ,   NC_FLOAT,  {ndim_name});
     ncf.def_var("refRatio",   NC_INT,    {flev_name});
     ncf.def_var("levelSteps", NC_INT,    {flev_name});

     ncf.def_var("Geom.smallend", NC_INT, {flev_name, ndim_name});
     ncf.def_var("Geom.bigend"  , NC_INT, {flev_name, ndim_name});
     ncf.def_var("CellSize"     , NC_FLOAT, {flev_name, ndim_name});

     ncf.def_var("x_grid", NC_FLOAT, {nb_name, nx_name});
     ncf.def_var("y_grid", NC_FLOAT, {nb_name, ny_name});
     ncf.def_var("z_grid", NC_FLOAT, {nb_name, nz_name});

     for (int i = 0; i < plot_var_names.size(); i++) {
         ncf.def_var(plot_var_names[i], NC_FLOAT, {nb_name, np_name});
     }

     ncf.exit_def_mode();


    if (lev == 0) {
      //
      // Write out the netcdf plotfile head information.
      //
      if (n_data_items == 0)
        amrex::Error("Must specify at least one valid data item to plot");

      ncf.put_attr("number_variables", std::vector<int>{n_data_items});
      ncf.put_attr("space_dimension", std::vector<int>{AMREX_SPACEDIM});
      ncf.put_attr("current_time", std::vector<double>{time});
      ncf.put_attr("FinestLevel", std::vector<int>{finest_level});
      ncf.put_attr("CurrentLevel", std::vector<int>{lev});

      amrex::Vector<Real> probLo;
      amrex::Vector<Real> probHi;
      for (int i = 0; i < AMREX_SPACEDIM; i++) {
        probLo.push_back(geom[lev].ProbLo(i));
        probHi.push_back(geom[lev].ProbHi(i));
      }

      auto nc_probLo = ncf.var("probLo");
      nc_probLo.par_access(NC_COLLECTIVE);
      nc_probLo.put(probLo.data(), {0}, {AMREX_SPACEDIM});

      auto nc_probHi = ncf.var("probHi");
      nc_probHi.par_access(NC_COLLECTIVE);
      nc_probHi.put(probHi.data(), {0}, {AMREX_SPACEDIM});


      amrex::Vector<int> refRatio;
      for (int i = 0; i < flev; i++)
        refRatio.push_back(ref_ratio[i][0]);
      auto nc_refRatio = ncf.var("refRatio");
      nc_refRatio.par_access(NC_COLLECTIVE);
      nc_refRatio.put(refRatio.data(), {0}, {static_cast<long unsigned int>(flev)});

      amrex::Vector<int> levelSteps;
      for (int i = 0; i <= flev; i++)
        levelSteps.push_back(level_steps[i]);
      auto nc_levelSteps = ncf.var("levelSteps");
      nc_levelSteps.par_access(NC_COLLECTIVE);
      nc_levelSteps.put(levelSteps.data(), {0}, {static_cast<long unsigned int>(flev)});

      amrex::Vector<int> smallend;
      amrex::Vector<int> bigend;
      for (int i = 0; i <= flev; i++) {
        smallend.clear(); bigend.clear();
        for (int j = 0; j < AMREX_SPACEDIM; j++) {
           smallend.push_back(geom[i].Domain().smallEnd(j));
             bigend.push_back(geom[i].Domain().bigEnd(j));
        }
        auto nc_Geom_smallend = ncf.var("Geom.smallend");
        nc_Geom_smallend.par_access(NC_COLLECTIVE);
        nc_Geom_smallend.put(smallend.data(), {static_cast<long unsigned int>(i), 0}, {1, AMREX_SPACEDIM});

        auto nc_Geom_bigend = ncf.var("Geom.bigend");
        nc_Geom_bigend.par_access(NC_COLLECTIVE);
        nc_Geom_bigend.put(bigend.data(), {static_cast<long unsigned int>(i), 0}, {1, AMREX_SPACEDIM});
      }

      amrex::Vector<Real> CellSize;
      for (int i = 0; i <= flev; i++) {
        CellSize.clear();
        for (int j = 0; j < AMREX_SPACEDIM; j++) {
          CellSize.push_back(geom[lev].CellSize()[j]);
        }
        auto nc_CellSize = ncf.var("CellSize");
        nc_CellSize.par_access(NC_COLLECTIVE);
        nc_CellSize.put(CellSize.data(), {static_cast<long unsigned int>(i), 0}, {1, AMREX_SPACEDIM});
      }

      ncf.put_attr("DefaultGeometry", std::vector<int>{amrex::DefaultGeometry().Coord()});

    } // level == 0

    std::vector<Real> x_grid;
    std::vector<Real> y_grid;
    std::vector<Real> z_grid;
    for (int i = 0; i < grids[lev].size(); ++i) {
        RealBox gridloc = RealBox(grids[lev][i], geom[lev].CellSize(), geom[lev].ProbLo());

        x_grid.clear(); y_grid.clear(); z_grid.clear();
        for (auto k1 = 0; k1 < grids[lev][i].length(0); ++k1) {
          x_grid.push_back(gridloc.lo(0)+geom[lev].CellSize(0)*static_cast<Real>(k1));
        }
        for (auto k2 = 0; k2 < grids[lev][i].length(1); ++k2) {
          y_grid.push_back(gridloc.lo(1)+geom[lev].CellSize(1)*static_cast<Real>(k2));
        }
        for (auto k3 = 0; k3 < grids[lev][i].length(2); ++k3) {
          z_grid.push_back(gridloc.lo(2)+geom[lev].CellSize(2)*static_cast<Real>(k3));
        }

        auto xlen = static_cast<long unsigned int>(grids[lev][i].length(0));
        auto ylen = static_cast<long unsigned int>(grids[lev][i].length(1));
        auto zlen = static_cast<long unsigned int>(grids[lev][i].length(2));
        auto nc_x_grid = ncf.var("x_grid");
        nc_x_grid.par_access(NC_COLLECTIVE);
        nc_x_grid.put(x_grid.data(), {static_cast<long unsigned int>(i), 0}, {1, xlen});
        auto nc_y_grid = ncf.var("y_grid");
        nc_y_grid.par_access(NC_COLLECTIVE);
        nc_y_grid.put(y_grid.data(), {static_cast<long unsigned int>(i), 0}, {1, ylen});
        auto nc_z_grid = ncf.var("z_grid");
        nc_z_grid.par_access(NC_COLLECTIVE);
        nc_z_grid.put(z_grid.data(), {static_cast<long unsigned int>(i), 0}, {1, zlen});
   }

   size_t nfai = 0;
   const int ncomp = plotMF[lev]->nComp();

   for (amrex::MFIter fai(*plotMF[lev]); fai.isValid(); ++fai) {
       auto box             = fai.validbox();
       long unsigned numpts = box.numPts();
       for (int k(0); k < ncomp; ++k) {
          auto data = plotMF[lev]->get(fai).dataPtr(k);
          auto nc_plot_var = ncf.var(plot_var_names[k]);
          nc_plot_var.par_access(NC_INDEPENDENT);
          nc_plot_var.put(data, {nfai+iproc, 0}, {1, numpts});
       }
       nfai++;
   }
   ncf.close();
}
