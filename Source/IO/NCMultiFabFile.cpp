#include <iomanip>
#include <iostream>
#include <sstream>
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

using namespace amrex;

void
ERF::ReadNCMultiFab (FabArray<FArrayBox> &mf,
                     const std::string  &mf_name,
                     int /*coordinatorProc*/,
                     int /*allow_empty_mf*/) {

    const std::string& FullPath = amrex::Concatenate(check_file,istep[0],5);

    Print() << "Reading MultiFab NetCDF checkpoint file to path: " << FullPath << "\n";

    if (amrex::ParallelDescriptor::IOProcessor())
    {
      static const std::string Suffix(mf_name);
      auto ncf = ncutils::NCFile::create(FullPath+Suffix, NC_CLOBBER | NC_NETCDF4);

      // int num_pts = static_cast<int>(ncf.dim("num_points").len());
      int ncomp   = static_cast<int>(ncf.dim("num_componts").len());
      int nbox    = static_cast<int>(ncf.dim("num_box").len());
      int ngrow   = 0;

      BoxList boxlist;
      amrex::Vector<std::string> plt_var_names;

      for (int nb(0); nb < nbox; ++nb) {
          amrex::IntVect smallend(AMREX_SPACEDIM);
          amrex::IntVect bigend(AMREX_SPACEDIM);
          amrex::IntVect btype(AMREX_SPACEDIM);

          auto nbb = static_cast<unsigned long int>(nb);
          ncf.var("SmallEnd").get(smallend.begin(), {0}, {nbb, AMREX_SPACEDIM});
          ncf.var("BigEnd"  ).get(bigend.begin()  , {0}, {nbb, AMREX_SPACEDIM});
          ncf.var("BoxType" ).get(btype.begin()   , {0}, {nbb, AMREX_SPACEDIM});

          std::stringstream ss;
          ss << btype;
          std::string typ_name = ss.str();
          plt_var_names.push_back(typ_name+"_"+std::to_string(nb));

          boxlist.push_back(amrex::Box(smallend, bigend, btype));
      }

      BoxArray boxarray(boxlist);
      DistributionMapping dm(boxarray);
      mf.define(boxarray, dm, ncomp, ngrow, MFInfo(), FArrayBoxFactory());

      for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
           auto ncomp_mf   = mf.nComp();
           auto num_pts_mf  = static_cast<long unsigned int>(mf.get(mfi).numPts());

           for (int k(0); k < ncomp_mf; ++k) {
               auto *dataPtr = mf.get(mfi).dataPtr(k);
               ncf.var(plt_var_names[mfi.index()*ncomp_mf+k]).get(dataPtr, {0}, {num_pts_mf});
            }
      }


   }
}


void
ERF::WriteNCMultiFab (const FabArray<FArrayBox> &fab,
                      const std::string& name,
                      bool /*set_ghost*/) {

    if (amrex::ParallelDescriptor::IOProcessor())
    {
      static const std::string Suffix{"_Data.nc"};
      auto ncf = ncutils::NCFile::create(name+Suffix, NC_CLOBBER | NC_NETCDF4);

      //
      // use separate name scope for data output
      //
      {
        // get the number points and size of different variable type, and geometric blocks
        const std::string nb_name     = "num_blocks";
        const std::string ndim_name   = "num_dimension";
        const std::string nvar_name   = "num_variables";

        // setup the plot variable names
        amrex::Vector<std::string> plt_var_names;
        amrex::Vector<std::string> npts_names;
        amrex::Vector<int> num_pts;
        for (MFIter mfi(fab); mfi.isValid(); ++mfi) {
           auto ncomp            = fab.nComp();
           std::string comp_name = std::to_string(mfi.index());
           int num_points        = fab.get(mfi).numPts();
           for (int k(0); k < ncomp; ++k) {
              plt_var_names.push_back("var_"+comp_name+"_"+std::to_string(k));
              npts_names.push_back("numpts_"+comp_name+"_"+std::to_string(k));
              num_pts.push_back(num_points);
           }
        }

        auto nvar      = plt_var_names.size();
        auto nbox      = fab.local_size();

        amrex::Vector<std::string> lo_names;
        amrex::Vector<std::string> hi_names;
        amrex::Vector<std::string> typ_names;
        for (int nb(0); nb < nbox; ++nb) {
           lo_names.push_back("SmallEnd_"+std::to_string(nb));
           hi_names.push_back("BigEnd_"+std::to_string(nb));
           typ_names.push_back("BoxType_"+std::to_string(nb));
        }

        ncf.enter_def_mode();
        ncf.put_attr("title", "ERF NetCDF MultiFab Data");

        ncf.def_dim(ndim_name, AMREX_SPACEDIM);
        ncf.def_dim(nb_name, nbox);
        ncf.def_dim(nvar_name, nvar);

        for (int k = 0; k < nvar; ++k) {
          ncf.def_dim(npts_names[k], num_pts[k]);
          ncf.def_var(plt_var_names[k], ncutils::NCDType::Real, {npts_names[k]});
        }

        for (int nb = 0; nb < nbox; ++nb) {
           ncf.def_var(lo_names[nb],  ncutils::NCDType::Int, {nb_name, ndim_name});
           ncf.def_var(hi_names[nb],  ncutils::NCDType::Int, {nb_name, ndim_name});
           ncf.def_var(typ_names[nb], ncutils::NCDType::Int, {nb_name, ndim_name});
        }

        ncf.exit_def_mode();

        for (MFIter mfi(fab); mfi.isValid(); ++mfi) {
            auto ncomp_mf   = fab.nComp();
            auto box        = fab.get(mfi).box();
            auto num_pts_mf = fab.get(mfi).numPts();

            amrex::IntVect smallend = box.smallEnd();
            amrex::IntVect bigend   = box.bigEnd();
            amrex::IntVect itype    = box.type();

            auto index = static_cast<long unsigned int>(mfi.index());
            ncf.var(lo_names[mfi.index()] ).put(smallend.begin(), {index, 0}, {1, AMREX_SPACEDIM});
            ncf.var(hi_names[mfi.index()] ).put(bigend.begin()  , {index, 0}, {1, AMREX_SPACEDIM});
            ncf.var(typ_names[mfi.index()]).put(itype.begin()   , {index, 0}, {1, AMREX_SPACEDIM});

            for (int k(0); k < ncomp_mf; ++k) {
                const auto *dataPtr = fab.get(mfi).dataPtr(k);
                ncf.var(plt_var_names[mfi.index()*ncomp_mf+k]).put(dataPtr, {0}, {static_cast<long unsigned int>(num_pts_mf)});
             }
        }
        ncf.close();
      }
   }
}
