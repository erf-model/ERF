#include "ERF.H"
#include "IOManager.H"
#include "NCInterface.H"
#include "IndexDefines.H"

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
    amrex::Error("Cannot writeToNCColumnFile: present implementation assumes only a single level is present");
  }

  // All processors: look for the requested column and get data if it's there
  amrex::Box probBox = erf.geom.Domain();
  const size_t nheights = probBox.length(2) + 2;
  amrex::Gpu::DeviceVector<amrex::Real> d_column_data(nheights*3, 0.0);
  amrex::Vector<amrex::Real> h_column_data(nheights*3, 0.0);
  amrex::Real* ucol = &d_column_data[0];
  amrex::Real* vcol = &d_column_data[nheights];
  amrex::Real* thetacol = &d_column_data[nheights*2];

  // Requested point must be inside problem domain
  if (xloc < erf.geom.ProbLo(0) || xloc > erf.geom.ProbHi(0) ||
      yloc < erf.geom.ProbLo(1) || yloc > erf.geom.ProbHi(1)) {
    amrex::Error("Invalid xy location to save column data - outside of domain");
  }

  // get indices and interpolation coefficients
  const amrex::Real x_cell_loc = probBox.smallEnd(0) + (xloc - erf.geom.ProbLo(0))* erf.geom.InvCellSize(0);
  const amrex::Real y_cell_loc = probBox.smallEnd(1) + (yloc - erf.geom.ProbLo(1))* erf.geom.InvCellSize(1);
  const int iloc = floor(x_cell_loc - 0.5);
  const int jloc = floor(y_cell_loc - 0.5);
  const amrex::Real alpha_x = x_cell_loc - 0.5 - iloc;
  const amrex::Real alpha_y = y_cell_loc - 0.5 - jloc;
  amrex::Array2D<amrex::Real, 0, 1, 0, 1> alpha_theta;
  alpha_theta(0,0) = (1.0 - alpha_x) * (1.0-alpha_y);
  alpha_theta(1,0) = (alpha_x) * (1.0-alpha_y);
  alpha_theta(0,1) = (1.0 - alpha_x) * (alpha_y);
  alpha_theta(1,1) = (alpha_x) * (alpha_y);
  // may need different indices for u,v due to not being collocated
  const int iloc_shift = floor(x_cell_loc) - iloc;
  const int jloc_shift = floor(y_cell_loc) - jloc;
  const amrex::Real alpha_x_u = x_cell_loc - iloc - iloc_shift;
  const amrex::Real alpha_y_v = y_cell_loc - jloc - jloc_shift;
  amrex::Array2D<amrex::Real, 0, 1, 0, 1> alpha_u;
  alpha_u(0,0) = (1.0 - alpha_x_u) * (1.0-alpha_y);
  alpha_u(1,0) = (alpha_x_u) * (1.0-alpha_y);
  alpha_u(0,1) = (1.0 - alpha_x_u) * (alpha_y);
  alpha_u(1,1) = (alpha_x_u) * (alpha_y);
  amrex::Array2D<amrex::Real, 0, 1, 0, 1> alpha_v;
  alpha_v(0,0) = (1.0 - alpha_x) * (1.0-alpha_y_v);
  alpha_v(1,0) = (alpha_x) * (1.0-alpha_y_v);
  alpha_v(0,1) = (1.0 - alpha_x) * (alpha_y_v);
  alpha_v(1,1) = (alpha_x) * (alpha_y_v);
  const int kstart = probBox.smallEnd(2)-1;
  const int kend = probBox.bigEnd(2)+1;
  const amrex::Box target_box(IntVect{iloc, jloc, kstart},
                              IntVect{iloc+1, jloc+1, kend});

  //  Need to get these valid in one grow cell for interpolation
  amrex::MultiFab& S_new = erf.get_new_data(State_Type);
  amrex::MultiFab& U_new = erf.get_new_data(X_Vel_Type);
  amrex::MultiFab& V_new = erf.get_new_data(Y_Vel_Type);
  S_new.FillBoundary(erf.geom.periodicity());
  U_new.FillBoundary(erf.geom.periodicity());
  V_new.FillBoundary(erf.geom.periodicity());
  amrex::Vector<MultiFab*> all_vars{&S_new, &U_new, &V_new};
  ERF::applyBCs(erf.geom,all_vars);

  // No tiling - we're just interested in one location
  for ( MFIter mfi(S_new); mfi.isValid(); ++mfi){
    const amrex::Array4<Real const> & state = S_new.array(mfi);
    const amrex::Array4<Real const> & velx = U_new.array(mfi);
    const amrex::Array4<Real const> & vely = V_new.array(mfi);

    // we want to include data at physical boundary ghost cells
    // for interpolation (i,j) or saving (k)
    amrex::Box vbox = mfi.validbox();
    for (int idir = 0; idir <=2; idir++) {
      if (vbox.smallEnd(idir) == probBox.smallEnd(idir)) {
        vbox.growLo(idir,1);
      }
      if (vbox.bigEnd(idir) == probBox.bigEnd(idir)) {
        vbox.growHi(idir,1);
      }
    }

    // Loop over indices where our box contains data needed for interpolation
    const amrex::Box overlap_box = vbox & target_box;
    if (overlap_box.ok()) {
      amrex::ParallelFor(overlap_box,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        const int idx_vec = k - kstart;
        const int ialpha = i - iloc;
        const int jalpha = j - jloc;
        ucol[idx_vec] += velx(i+iloc_shift,j,k) * alpha_u(ialpha, jalpha);
        vcol[idx_vec] += vely(i,j+jloc_shift,k) * alpha_v(ialpha, jalpha);
        thetacol[idx_vec] += state(i,j,k,RhoTheta_comp) / state(i,j,k,Rho_comp)
          * alpha_theta(ialpha, jalpha);
      });
    }
  }

  // Communicate values to CPU then to the IO processor
  amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_column_data.begin(),
    d_column_data.end(), h_column_data.begin());
  amrex::ParallelDescriptor::ReduceRealSum(h_column_data.data(),
    h_column_data.size(), amrex::ParallelDescriptor::IOProcessorNumber());

  // IO processor only: write the relevant data to file
  if (amrex::ParallelDescriptor::IOProcessor()) {
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
    ncf.var("wrf_momentum_u").put(&h_column_data[0], start, count);
    ncf.var("wrf_momentum_v").put(&h_column_data[nheights], start, count);
    ncf.var("wrf_temperature").put(&h_column_data[2*nheights], start, count);
    ncf.close();
  }
}
