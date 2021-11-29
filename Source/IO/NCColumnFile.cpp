#include "ERF.H"
#include "IOManager.H"
#include "NCInterface.H"
#include "IndexDefines.H"

#ifndef ERF_USE_NETCDF
// If not compiling with NetCDF, fail gracefully
void IOManager::createNCColumnFile(const std::string& colfile_name,
				   const amrex::Real xloc,
				   const amrex::Real yloc)
{
  amrex::Error("To save output data in NetCDF format, you must enable NetCDF at compile time");
}

void IOManager::writeToNCColumnFile(const std::string& colfile_name, const amrex::Real xloc, const amrex::Real yloc)
{
  amrex::Error("To save output data in NetCDF format, you must enable NetCDF at compile time");
}

#else
void IOManager::createNCColumnFile(const std::string& colfile_name,
				   const amrex::Real xloc,
				   const amrex::Real yloc)
{
  // TODO: Relax assumption of only a single level
  if (erf.parent->finestLevel() > 0) {
    amrex::Error("Cannot createNCColumnFile: present implementation assumes only a single level is present");
  }

  // Create file to which column data will be written every timestep
  if (amrex::ParallelDescriptor::IOProcessor()) {
    std::cout << "IO Processor" << std::endl;
    auto ncf = ncutils::NCFile::create(colfile_name, NC_CLOBBER | NC_NETCDF4);
    const std::string nt_name = "ntime";
    const std::string nh_name = "nheight";
    // Use one grow cell (on either side) to allow interpolation to boundaries
    const int nheights = erf.geom.Domain().length(2) + 2;
    ncf.enter_def_mode();
    ncf.put_attr("title", "ERF NetCDF Vertical Column Output");
    ncf.put_attr("units", "mks");
    amrex::Vector<amrex::Real> loc = {xloc, yloc};
    ncf.put_attr("location", loc);
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim(nh_name, nheights);
    ncf.def_var("times", NC_FLOAT, {nt_name});
    ncf.def_var("heights", NC_FLOAT, {nh_name});
    ncf.def_var("wrf_momentum_u", NC_FLOAT, {nt_name, nh_name});
    ncf.def_var("wrf_momentum_v", NC_FLOAT, {nt_name, nh_name});
    ncf.def_var("wrf_temperature", NC_FLOAT, {nt_name, nh_name});
    ncf.def_var("wrf_tflux", NC_FLOAT, {nt_name});
    ncf.exit_def_mode();

    // Put in the Z grid, but not any actual data yet
    amrex::Real zmin = erf.geom.ProbLo(2);
    amrex::Real dz = erf.geom.CellSize(2);
    amrex::Vector<amrex::Real> zvalues(nheights, zmin-0.5*dz);
    for (int ii = 0; ii < nheights; ++ii) {
      zvalues[ii] += ii * dz;
    }
    ncf.var("heights").put(zvalues.data());
    ncf.close();
  }
}

void IOManager::writeToNCColumnFile(const std::string& colfile_name, const amrex::Real xloc, const amrex::Real yloc)
{
  // TODO: Relax assumption of only a single level
  if (erf.parent->finestLevel() > 0) {
    amrex::Error("Cannot createNCColumnFile: present implementation assumes only a single level is present");
  }

  // All processors: look for the requested column and get data if it's there
  amrex::Box probBox = erf.geom.Domain();
  const size_t nheights = probBox.length(2) + 2;
  amrex::Vector<amrex::Real> Ucol(nheights, 0.0);
  amrex::Vector<amrex::Real> Vcol(nheights, 0.0);
  amrex::Vector<amrex::Real> Thetacol(nheights, 0.0);

  const amrex::Real problo_x = erf.geom.ProbLo(0);
  const amrex::Real problo_y = erf.geom.ProbLo(1);
  if (xloc < erf.geom.ProbLo(0) || xloc > erf.geom.ProbHi(0) ||
      yloc < erf.geom.ProbLo(1) || yloc > erf.geom.ProbHi(1)) {
    amrex::Error("Invalid xy location to save column data - outside of domain");
  }

  const amrex::Real x_from_lo = (xloc - erf.geom.ProbLo(0));
  const amrex::Real y_from_lo = (yloc - erf.geom.ProbLo(1));
  int iloc = probBox.smallEnd(0) + x_from_lo / probBox.length(0);
  int jloc = probBox.smallEnd(1) + y_from_lo / probBox.length(1);

  //std::cout << "Box " << probBox << std::endl;
  //std::cout << "geom " << erf.geom << std::endl;

  // These are valid in one grow cell?
  amrex::MultiFab& S_new = erf.get_new_data(State_Type);
  amrex::MultiFab& U_new = erf.get_new_data(X_Vel_Type);
  amrex::MultiFab& V_new = erf.get_new_data(Y_Vel_Type);
  //std::cout << std::endl << "GROW cells " << S_new.nGrow() << std::endl << std::endl;

  // No tiling - we're just interested in one location
  for ( MFIter mfi(S_new); mfi.isValid(); ++mfi){
    const amrex::Box &vbox = mfi.validbox();
    const amrex::Box &abox = mfi.fabbox();
    const amrex::Box &tbox = mfi.tilebox();
    const amrex::Box &ntxbox = mfi.nodaltilebox(0);
    const amrex::Box &ntybox = mfi.nodaltilebox(1);
    const amrex::Box &gntxbox = mfi.grownnodaltilebox(0);
    const amrex::Box &gntybox = mfi.grownnodaltilebox(1);
    //const amrex::Box &vboxu = U_new[mfi];//.validbox();
    //const amrex::Box &vboxv = V_new[mfi];//.validbox();
    std::cout << "MFI VBOX "  << vbox << std::endl
	      << "MFI ABOX " << abox << std::endl
	      << "MFI TBOX " << tbox << std::endl
	      << "MFI NTXBOX " << ntxbox << std::endl
	      << "MFI NTYBOX " << ntybox << std::endl
	      << "MFI GNTXBOX " << gntxbox << std::endl
	      << "MFI GNTYBOX " << gntybox << std::endl;
  }

  // IO processor only: write the relevant data to file
  if (amrex::ParallelDescriptor::IOProcessor()) {
    std::cout << "IO Processor - next step" << std::endl;
    auto ncf = ncutils::NCFile::open(colfile_name, NC_WRITE | NC_NETCDF4);
    size_t putloc = ncf.dim("ntime").len();

    // Time
    amrex::Real cumtime = erf.parent->cumTime();
    amrex::Vector<size_t> start_t {putloc};
    amrex::Vector<size_t> count_t {1};
    ncf.var("times").put(&cumtime, start_t, count_t);

    // T flux
    // TODO: Make this the actual flux rather than just a placeholder
    amrex::Real Tflux = 0.0;
    ncf.var("wrf_tflux").put(&Tflux, start_t, count_t);

    // U, V, Theta
    std::vector<size_t> start = {putloc, 0};
    std::vector<size_t> count = {1, nheights};
    ncf.var("wrf_momentum_u").put(Ucol.data(), start, count);
    ncf.var("wrf_momentum_v").put(Vcol.data(), start, count);
    ncf.var("wrf_temperature").put(Thetacol.data(), start, count);
    ncf.close();
  }
}

#endif
