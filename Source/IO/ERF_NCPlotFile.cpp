#include <iomanip>
#include <iostream>
#include <string>
#include <ctime>

#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>

#include "ERF_NCInterface.H"

using namespace amrex;

void
writeNCPlotFile (int lev, int which_subdomain, const std::string& dir,
                 const Vector<const MultiFab*> &plotMF,
                 const Vector<std::string> &plot_var_names,
                 const Vector<int>& /*level_steps*/,
                 Array<Real,AMREX_SPACEDIM> prob_lo,
                 Array<Real,AMREX_SPACEDIM> prob_hi,
                 Array<Real,AMREX_SPACEDIM> dx_in,
                 const Box& subdomain,
                 const Real time, const Real start_bdy_time)
{
    //
    // Set the full IO path for NetCDF output
    //
    std::string FullPath = dir;
    if (lev == 0) {
        const std::string& extension = amrex::Concatenate("_d",lev+1,2);
        FullPath += extension + ".nc";
    } else {
        const std::string& extension = amrex::Concatenate("_d",lev+1+which_subdomain,2);
        FullPath += extension + ".nc";
    }

    Print() << "Writing level " << lev << " NetCDF plot file " << FullPath << std::endl;

    //
    // Open netcdf file to write data
    //
    auto ncf = ncutils::NCFile::create_par(FullPath, NC_NETCDF4 | NC_MPIIO,
                                           amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

    auto ba = plotMF[lev]->boxArray();
    auto dm = plotMF[lev]->DistributionMap();

    int nblocks = ba.size();

    int nx = subdomain.length(0);
    int ny = subdomain.length(1);
    int nz = subdomain.length(2);

    int num_pts = nx*ny*nz;

    int n_data_items = plotMF[lev]->nComp();

    //
    // Start Define stuff
    //
    ncf.enter_def_mode();

    const std::string nt_name   = "num_time_steps";
    const std::string nb_name   = "num_blocks";
    const std::string np_name   = "num_pts";
    const std::string nx_name   = "nx";
    const std::string ny_name   = "ny";
    const std::string nz_name   = "nz";

    const std::string ndim_name = "num_geo_dims";

    ncf.put_attr("title", "ERF NetCDF Plot data output");

    ncf.def_dim(ndim_name, AMREX_SPACEDIM);
    ncf.def_dim(np_name  , num_pts);
    ncf.def_dim(nb_name  , nblocks);

    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim(nx_name, nx);
    ncf.def_dim(ny_name, ny);
    ncf.def_dim(nz_name, nz);

    ncf.def_var("probLo"  ,   NC_FLOAT,  {ndim_name});
    ncf.def_var("probHi"  ,   NC_FLOAT,  {ndim_name});

    ncf.def_var("Geom.smallend", NC_INT, {ndim_name});
    ncf.def_var("Geom.bigend"  , NC_INT, {ndim_name});
    ncf.def_var("CellSize"     , NC_FLOAT, {ndim_name});

    ncf.def_var("x_grid", NC_DOUBLE, {np_name});
    ncf.def_var("y_grid", NC_DOUBLE, {np_name});
    ncf.def_var("z_grid", NC_DOUBLE, {np_name});

    for (int i = 0; i < plot_var_names.size(); i++) {
        ncf.def_var(plot_var_names[i], NC_DOUBLE, {nz_name, ny_name, nx_name});
    }

    ncf.exit_def_mode();
    //
    // End Define stuff
    //

    // We are doing single-level writes but it doesn't have to be level 0
    //
    // Write out the netcdf plotfile head information.
    //
    if (n_data_items == 0) {
      amrex::Error("Must specify at least one valid data item to plot");
    }

    ncf.put_attr("number_variables", std::vector<int>{n_data_items});
    ncf.put_attr("space_dimension", std::vector<int>{AMREX_SPACEDIM});
    ncf.put_attr("current_time", std::vector<double>{time});
    ncf.put_attr("start_time", std::vector<double>{start_bdy_time});
    ncf.put_attr("CurrentLevel", std::vector<int>{lev});

    Real dx[AMREX_SPACEDIM];
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
       dx[i] = dx_in[i];
    }

    amrex::Vector<Real> probLo;
    amrex::Vector<Real> probHi;
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
      probLo.push_back(prob_lo[i]);
      probHi.push_back(prob_hi[i]);
    }

    auto nc_probLo = ncf.var("probLo");
    nc_probLo.par_access(NC_COLLECTIVE);
    nc_probLo.put(probLo.data(), {0}, {AMREX_SPACEDIM});

    auto nc_probHi = ncf.var("probHi");
    nc_probHi.par_access(NC_COLLECTIVE);
    nc_probHi.put(probHi.data(), {0}, {AMREX_SPACEDIM});

    amrex::Vector<int> smallend;
    amrex::Vector<int> bigend;
    smallend.clear(); bigend.clear();
    for (int j = 0; j < AMREX_SPACEDIM; j++) {
       smallend.push_back(subdomain.smallEnd(j));
         bigend.push_back(subdomain.bigEnd(j));
    }

    auto nc_Geom_smallend = ncf.var("Geom.smallend");
    nc_Geom_smallend.par_access(NC_COLLECTIVE);
    nc_Geom_smallend.put(smallend.data(), {0}, {AMREX_SPACEDIM});

    auto nc_Geom_bigend = ncf.var("Geom.bigend");
    nc_Geom_bigend.par_access(NC_COLLECTIVE);
    nc_Geom_bigend.put(bigend.data(), {0}, {AMREX_SPACEDIM});

    amrex::Vector<Real> CellSize;
    CellSize.clear();
    for (double & j : dx) {
        CellSize.push_back(j);
    }
    auto nc_CellSize = ncf.var("CellSize");
    nc_CellSize.par_access(NC_COLLECTIVE);
    nc_CellSize.put(CellSize.data(), {0}, {AMREX_SPACEDIM});

    ncf.put_attr("DefaultGeometry", std::vector<int>{amrex::DefaultGeometry().Coord()});

    std::vector<Real> x_grid;
    std::vector<Real> y_grid;
    std::vector<Real> z_grid;
    long unsigned goffset = 0;
    long unsigned glen    = 0;

    // *******************************************************************************
    // NOTE: the (x,y,z) output here are for a mesh withOUT terrain-fitted coordinates
    // *******************************************************************************
    for (int i = 0; i < ba.size(); ++i) {
        auto bx = ba[i];
        if (subdomain.contains(bx)) {
            x_grid.clear(); y_grid.clear(); z_grid.clear();
            for (auto k3 = 0; k3 < bx.length(2); ++k3) {
                for (auto k2 = 0; k2 < bx.length(1); ++k2) {
                    for (auto k1 = 0; k1 < bx.length(0); ++k1) {
                        x_grid.push_back(prob_lo[0]+dx[0]*(static_cast<Real>(k1)+0.5));
                        y_grid.push_back(prob_lo[1]+dx[1]*(static_cast<Real>(k2)+0.5));
                        z_grid.push_back(prob_lo[2]+dx[2]*(static_cast<Real>(k3)+0.5));
                     }
                }
            }

            goffset += glen;
            glen = bx.numPts();

            auto nc_x_grid = ncf.var("x_grid");
            auto nc_y_grid = ncf.var("y_grid");
            auto nc_z_grid = ncf.var("z_grid");

            nc_x_grid.par_access(NC_COLLECTIVE);
            nc_y_grid.par_access(NC_COLLECTIVE);
            nc_z_grid.par_access(NC_COLLECTIVE);

            nc_x_grid.put(x_grid.data(), {goffset}, {glen});
            nc_y_grid.put(y_grid.data(), {goffset}, {glen});
            nc_z_grid.put(z_grid.data(), {goffset}, {glen});
       }
   }

   const int ncomp = plotMF[lev]->nComp();

   for (MFIter mfi(*plotMF[lev]); mfi.isValid(); ++mfi)
   {
       auto bx = mfi.validbox();

       if (subdomain.contains(bx))

       {
           //
           // These are the dimensions of the data we write for only this box
           //
           long unsigned local_nx = bx.length()[0];
           long unsigned local_ny = bx.length()[1];
           long unsigned local_nz = bx.length()[2];

           long unsigned local_start_x  = static_cast<long unsigned>(bx.smallEnd()[0]-subdomain.smallEnd()[0]);
           long unsigned local_start_y  = static_cast<long unsigned>(bx.smallEnd()[1]-subdomain.smallEnd()[1]);
           long unsigned local_start_z  = static_cast<long unsigned>(bx.smallEnd()[2]-subdomain.smallEnd()[2]);

           for (int k(0); k < ncomp; ++k) {
               FArrayBox tmp;
               tmp.resize(bx, 1, amrex::The_Pinned_Arena());
               tmp.template copy<RunOn::Device>((*plotMF[lev])[mfi.index()], k, 0, 1);
               Gpu::streamSynchronize();

               auto nc_plot_var = ncf.var(plot_var_names[k]);
               nc_plot_var.par_access(NC_COLLECTIVE);
               nc_plot_var.put(tmp.dataPtr(), {local_start_z,local_start_y,local_start_x},
                                              {local_nz, local_ny, local_nx});
           }
       }
   }
   ncf.close();
}
