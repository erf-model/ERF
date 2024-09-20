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
#include "ERF_NCInterface.H"
#include "ERF_NCPlotFile.H"
#include "ERF_IndexDefines.H"

using namespace amrex;

void
ERF::writeNCPlotFile (int lev, int which_subdomain, const std::string& dir,
                      const Vector<const MultiFab*> &plotMF,
                      const Vector<std::string> &plot_var_names,
                      const Vector<int>& /*level_steps*/, const Real time) const
{
     // get the processor number
     int iproc = amrex::ParallelContext::MyProcAll();
     int nproc = amrex::ParallelDescriptor::NProcs();

     // total number of cells in this "domain" at this level
     std::vector<int> n_cells;

     // number of points in each block at this level
     std::vector<int> offset;

     // set the full IO path for NetCDF output
     std::string FullPath = dir;
     if (lev == 0) {
         const std::string& extension = amrex::Concatenate("_d",lev+1,2);
         FullPath += extension + ".nc";
     } else {
         const std::string& extension = amrex::Concatenate("_d",lev+1+which_subdomain,2);
         FullPath += extension + ".nc";
     }

     Print() << "Writing level " << lev << " NetCDF plot file " << FullPath << std::endl;

     // open netcdf file to write data
     auto ncf = ncutils::NCFile::create_par(FullPath, NC_NETCDF4 | NC_MPIIO,
                                            amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

     int nblocks = grids[lev].size();
     auto dm = plotMF[lev]->DistributionMap();
     offset.reserve(nproc);
     for(auto n = 0; n < nproc; n++) {
         offset[n] = 0;
     }
     for(auto ib=0; ib<nblocks; ib++) {
        auto npts_per_block = grids[lev][ib].length(0)*grids[lev][ib].length(1)*grids[lev][ib].length(2);
        offset[dm[ib]] += npts_per_block;
     }

     // We only do single-level writes when using NetCDF format
     int flev = lev;

     Box subdomain;
     if (lev == 0) {
         subdomain = geom[lev].Domain();
     } else {
         subdomain = boxes_at_level[lev][which_subdomain];
     }

     int nx = subdomain.length(0);
     int ny = subdomain.length(1);
     int nz = subdomain.length(2);

     n_cells.push_back(nx);
     n_cells.push_back(ny);
     n_cells.push_back(nz);

     int num_pts = nx*ny*nz;

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

     ncf.def_var("Geom.smallend", NC_INT, {flev_name, ndim_name});
     ncf.def_var("Geom.bigend"  , NC_INT, {flev_name, ndim_name});
     ncf.def_var("CellSize"     , NC_FLOAT, {flev_name, ndim_name});

     ncf.def_var("x_grid", NC_FLOAT, {np_name});
     ncf.def_var("y_grid", NC_FLOAT, {np_name});
     ncf.def_var("z_grid", NC_FLOAT, {np_name});

     for (int i = 0; i < plot_var_names.size(); i++) {
         ncf.def_var(plot_var_names[i], NC_FLOAT, {np_name});
     }

     ncf.exit_def_mode();

     {
      // We are doing single-level writes but it doesn't have to be level 0
      //
      // Write out the netcdf plotfile head information.
      //
      if (n_data_items == 0)
        amrex::Error("Must specify at least one valid data item to plot");

      ncf.put_attr("number_variables", std::vector<int>{n_data_items});
      ncf.put_attr("space_dimension", std::vector<int>{AMREX_SPACEDIM});
      ncf.put_attr("current_time", std::vector<double>{time});
      ncf.put_attr("start_time", std::vector<double>{start_bdy_time});
      ncf.put_attr("CurrentLevel", std::vector<int>{flev});

      Real dx[AMREX_SPACEDIM];
      for (int i = 0; i < AMREX_SPACEDIM; i++)
         dx[i] = geom[lev].CellSize()[i];
      const auto *base = geom[lev].ProbLo();
      RealBox rb(subdomain,dx,base);

      amrex::Vector<Real> probLo;
      amrex::Vector<Real> probHi;
      for (int i = 0; i < AMREX_SPACEDIM; i++) {
        probLo.push_back(rb.lo(i));
        probHi.push_back(rb.hi(i));
      }

      auto nc_probLo = ncf.var("probLo");
      nc_probLo.par_access(NC_COLLECTIVE);
      nc_probLo.put(probLo.data(), {0}, {AMREX_SPACEDIM});

      auto nc_probHi = ncf.var("probHi");
      nc_probHi.par_access(NC_COLLECTIVE);
      nc_probHi.put(probHi.data(), {0}, {AMREX_SPACEDIM});

      amrex::Vector<int> smallend;
      amrex::Vector<int> bigend;
      for (int i = lev; i <= flev; i++) {
        smallend.clear(); bigend.clear();
        for (int j = 0; j < AMREX_SPACEDIM; j++) {
           smallend.push_back(subdomain.smallEnd(j));
             bigend.push_back(subdomain.bigEnd(j));
        }
        auto nc_Geom_smallend = ncf.var("Geom.smallend");
        nc_Geom_smallend.par_access(NC_COLLECTIVE);
        nc_Geom_smallend.put(smallend.data(), {static_cast<long unsigned int>(i-lev), 0}, {1, AMREX_SPACEDIM});

        auto nc_Geom_bigend = ncf.var("Geom.bigend");
        nc_Geom_bigend.par_access(NC_COLLECTIVE);
        nc_Geom_bigend.put(bigend.data(), {static_cast<long unsigned int>(i-lev), 0}, {1, AMREX_SPACEDIM});
      }

      amrex::Vector<Real> CellSize;
      for (int i = lev; i <= flev; i++) {
        CellSize.clear();
        for (double & j : dx) {
          CellSize.push_back(j);
        }
        auto nc_CellSize = ncf.var("CellSize");
        nc_CellSize.par_access(NC_COLLECTIVE);
        nc_CellSize.put(CellSize.data(), {static_cast<long unsigned int>(i-lev), 0}, {1, AMREX_SPACEDIM});
      }

      ncf.put_attr("DefaultGeometry", std::vector<int>{amrex::DefaultGeometry().Coord()});
    }

    std::vector<Real> x_grid;
    std::vector<Real> y_grid;
    std::vector<Real> z_grid;
    long unsigned goffset = 0;
    long unsigned glen    = 0;
    for (int i = 0; i < grids[lev].size(); ++i) {
        auto box = grids[lev][i];
        if (subdomain.contains(box)) {
            RealBox gridloc = RealBox(grids[lev][i], geom[lev].CellSize(), geom[lev].ProbLo());

            x_grid.clear(); y_grid.clear(); z_grid.clear();
            for (auto k1 = 0; k1 < grids[lev][i].length(0); ++k1) {
              for (auto k2 = 0; k2 < grids[lev][i].length(1); ++k2) {
                 for (auto k3 = 0; k3 < grids[lev][i].length(2); ++k3) {
                    x_grid.push_back(gridloc.lo(0)+geom[lev].CellSize(0)*static_cast<Real>(k1));
                    y_grid.push_back(gridloc.lo(1)+geom[lev].CellSize(1)*static_cast<Real>(k2));
                    z_grid.push_back(gridloc.lo(2)+geom[lev].CellSize(2)*static_cast<Real>(k3));
                 }
              }
            }

            goffset += glen;
            glen = grids[lev][i].length(0)*grids[lev][i].length(1)*grids[lev][i].length(2);

            auto nc_x_grid = ncf.var("x_grid");
            auto nc_y_grid = ncf.var("y_grid");
            auto nc_z_grid = ncf.var("z_grid");

            nc_x_grid.par_access(NC_INDEPENDENT);
            nc_y_grid.par_access(NC_INDEPENDENT);
            nc_z_grid.par_access(NC_INDEPENDENT);

            nc_x_grid.put(x_grid.data(), {goffset}, {glen});
            nc_y_grid.put(y_grid.data(), {goffset}, {glen});
            nc_z_grid.put(z_grid.data(), {goffset}, {glen});
       }
   }

   size_t nfai = 0;
   long unsigned numpts = 0;
   const int ncomp = plotMF[lev]->nComp();

   for (MFIter fai(*plotMF[lev]); fai.isValid(); ++fai) {
       auto box = fai.validbox();
       if (subdomain.contains(box)) {
           numpts = box.numPts();
           long unsigned diff = nfai*numpts;
           for(auto ip = 1; ip <= iproc; ++ip) diff += offset[ip-1];

           for (int k(0); k < ncomp; ++k) {
              const auto *data = plotMF[lev]->get(fai).dataPtr(k);
              auto nc_plot_var = ncf.var(plot_var_names[k]);
              nc_plot_var.par_access(NC_INDEPENDENT);
              nc_plot_var.put(data, {diff}, {numpts});
           }
           nfai++;
       }
   }
   ncf.close();
}
