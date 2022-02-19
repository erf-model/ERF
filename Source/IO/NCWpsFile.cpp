#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <ctime>

#include "NCInterface.H"
#include <ERF.H>

void ReadNCWpsFile(const std::string  &fname) {

    amrex::Print() << "Reading NetCDF WPS file: " << fname << "\n";
   
    if (amrex::ParallelDescriptor::IOProcessor())
    {
      auto ncf = ncutils::NCFile::create(fname, NC_CLOBBER | NC_NETCDF4);
      
      // get the dimension information
      int DateStrLen         = static_cast<int>(ncf.dim("DateStrLen").len());
      int Time               = static_cast<int>(ncf.dim("Time").len());
      int west_east          = static_cast<int>(ncf.dim("west_east").len());
      int south_north        = static_cast<int>(ncf.dim("south_north").len());
      int num_metgrid_levels = static_cast<int>(ncf.dim("num_metgrid_levels").len());
      int num_st_layers      = static_cast<int>(ncf.dim("num_st_layers").len());
      int num_sm_layers      = static_cast<int>(ncf.dim("num_sm_layers").len());
      int south_north_stag   = static_cast<int>(ncf.dim("south_north_stag").len());
      int west_east_stag     = static_cast<int>(ncf.dim("west_east_stag").len());
      int z_dimension0012    = static_cast<int>(ncf.dim("z-dimension0012").len());
      int z_dimension0016    = static_cast<int>(ncf.dim("z-dimension0016").len());
      int z_dimension0024    = static_cast<int>(ncf.dim("z-dimension0024").len());

      // get the pressure
      std::vector<size_t> pshape = ncf.var("PRES").shape();
      auto pressure = new typename ncutils::NCDType::RType[pshape.size()];
      ncf.var("PRES").get(pressure, {0}, {pshape.size()});

      // get uu
      std::vector<size_t> ushape = ncf.var("UU").shape();
      auto uu = new typename ncutils::NCDType::RType[ushape.size()];
      ncf.var("UU").get(uu, {0}, {ushape.size()});

      // get vv
      std::vector<size_t> vshape = ncf.var("VV").shape();
      auto vv = new typename ncutils::NCDType::RType[vshape.size()];
      ncf.var("VV").get(vv, {0}, {vshape.size()});

      // get tt
      std::vector<size_t> tshape = ncf.var("TT").shape();
      auto tt = new typename ncutils::NCDType::RType[tshape.size()];
      ncf.var("TT").get(tt, {0}, {tshape.size()});

      // get psfc
      std::vector<size_t> psshape = ncf.var("PSFC").shape();
      auto psfc = new typename ncutils::NCDType::RType[psshape.size()];
      ncf.var("PSFC").get(psfc, {0}, {psshape.size()});
   }
}
