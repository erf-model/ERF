#include "ERF.H"
#include "NCInterface.H"
#include "IndexDefines.H"

using namespace amrex;

/**
 * Creates a NetCDF file containing column data we write to during runtime
 *
 * @param lev Current level
 * @param colfile_name Name of the NetCDF file containing column data
 * @param xloc Location of the column in the x-dimension
 * @param yloc Location of the column in the y-dimension
 */
void
ERF::createNCColumnFile (int lev,
                         const std::string& colfile_name,
                         const Real xloc,
                         const Real yloc)
{
  // Create file to which column data will be written every timestep
  if (amrex::ParallelDescriptor::IOProcessor()) {
    auto ncf = ncutils::NCFile::create(colfile_name, NC_CLOBBER | NC_NETCDF4);
    const std::string nt_name = "ntime";
    const std::string nh_name = "nheight";
    // Use one grow cell (on either side) to allow interpolation to boundaries
    const int nheights = geom[lev].Domain().length(2) + 2;
    ncf.enter_def_mode();
    ncf.put_attr("title", "ERF NetCDF Vertical Column Output");
    ncf.put_attr("units", "mks");
    amrex::Vector<Real> loc = {xloc, yloc};
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
    Real zmin = geom[lev].ProbLo(2);
    Real dz = geom[lev].CellSize(2);
    amrex::Vector<Real> zvalues(nheights, zmin-0.5*dz);
    for (int ii = 0; ii < nheights; ++ii) {
      zvalues[ii] += ii * dz;
    }
    ncf.var("heights").put(zvalues.data());
    ncf.close();
  }
}

/**
 * Writes column data from a given level to the NetCDF column data file
 *
 * @param lev Current level
 * @param colfile_name Name of the NetCDF file containing column data
 * @param xloc Location of the column in the x-dimension
 * @param yloc Location of the column in the y-dimension
 * @param cumtime Current time
 */
void
ERF::writeToNCColumnFile (const int lev,
                          const std::string& colfile_name, const Real xloc, const Real yloc,
                          const Real cumtime)
{
  //
  // This routine assumes that we can grab the whole column of data from the MultiFabs at
  //     a single level, "lev".  This assumption is true as long as we don't refine only
  //     partway up a column, which is the plan.
  //

  // All processors: look for the requested column and get data if it's there
  amrex::Box probBox = geom[lev].Domain();
  const size_t nheights = probBox.length(2) + 2;
  amrex::Gpu::DeviceVector<Real> d_column_data(nheights*3, 0.0);
  amrex::Vector<Real> h_column_data(nheights*3, 0.0);
  Real* ucol = &d_column_data[0];
  Real* vcol = &d_column_data[nheights];
  Real* thetacol = &d_column_data[nheights*2];

  // Requested point must be inside problem domain
  if (xloc < geom[0].ProbLo(0) || xloc > geom[0].ProbHi(0) ||
      yloc < geom[0].ProbLo(1) || yloc > geom[0].ProbHi(1)) {
    amrex::Error("Invalid xy location to save column data - outside of domain");
  }

  // get indices and interpolation coefficients
  const Real x_cell_loc = probBox.smallEnd(0) + (xloc - geom[lev].ProbLo(0))* geom[lev].InvCellSize(0);
  const Real y_cell_loc = probBox.smallEnd(1) + (yloc - geom[lev].ProbLo(1))* geom[lev].InvCellSize(1);
  const int iloc = static_cast<int>(floor(x_cell_loc - 0.5));
  const int jloc = static_cast<int>(floor(y_cell_loc - 0.5));
  const Real alpha_x = x_cell_loc - 0.5 - iloc;
  const Real alpha_y = y_cell_loc - 0.5 - jloc;
  amrex::Array2D<Real, 0, 1, 0, 1> alpha_theta;
  alpha_theta(0,0) = (1.0 - alpha_x) * (1.0-alpha_y);
  alpha_theta(1,0) = (alpha_x) * (1.0-alpha_y);
  alpha_theta(0,1) = (1.0 - alpha_x) * (alpha_y);
  alpha_theta(1,1) = (alpha_x) * (alpha_y);
  // may need different indices for u,v due to not being collocated
  const int iloc_shift = static_cast<int>(floor(x_cell_loc)) - iloc;
  const int jloc_shift = static_cast<int>(floor(y_cell_loc)) - jloc;
  const Real alpha_x_u = x_cell_loc - iloc - iloc_shift;
  const Real alpha_y_v = y_cell_loc - jloc - jloc_shift;
  amrex::Array2D<Real, 0, 1, 0, 1> alpha_u;
  alpha_u(0,0) = (1.0 - alpha_x_u) * (1.0-alpha_y);
  alpha_u(1,0) = (alpha_x_u) * (1.0-alpha_y);
  alpha_u(0,1) = (1.0 - alpha_x_u) * (alpha_y);
  alpha_u(1,1) = (alpha_x_u) * (alpha_y);
  amrex::Array2D<Real, 0, 1, 0, 1> alpha_v;
  alpha_v(0,0) = (1.0 - alpha_x) * (1.0-alpha_y_v);
  alpha_v(1,0) = (alpha_x) * (1.0-alpha_y_v);
  alpha_v(0,1) = (1.0 - alpha_x) * (alpha_y_v);
  alpha_v(1,1) = (alpha_x) * (alpha_y_v);
  const int kstart = probBox.smallEnd(2)-1;
  const int kend = probBox.bigEnd(2)+1;
  const amrex::Box target_box(IntVect{iloc, jloc, kstart},
                              IntVect{iloc+1, jloc+1, kend});

  //  Need data in one grow cell for interpolation
  //  Note that vars_new is what's filled here; rU_new/rV_new/rW_new are just used as scratch space
  FillPatch(lev, t_new[lev], {&vars_new[lev][Vars::cons], &vars_new[lev][Vars::xvel],
                              &vars_new[lev][Vars::yvel], &vars_new[lev][Vars::zvel]},
                             {&vars_new[lev][Vars::cons], &rU_new[lev],
                              &rV_new[lev], &rW_new[lev]});

  MultiFab& S_new = vars_new[lev][Vars::cons];
  MultiFab& U_new = vars_new[lev][Vars::xvel];
  MultiFab& V_new = vars_new[lev][Vars::yvel];

  // No tiling - we're just interested in one location
  for ( MFIter mfi(S_new); mfi.isValid(); ++mfi){
    const amrex::Array4<Real const> & state = S_new.array(mfi);
    const amrex::Array4<Real const> & velx  = U_new.array(mfi);
    const amrex::Array4<Real const> & vely  = V_new.array(mfi);

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
      ParallelFor(overlap_box,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        const int idx_vec = k - kstart;
        const int ialpha = i - iloc;
        const int jalpha = j - jloc;

        auto tmpx = velx(i+iloc_shift,j,k) * alpha_u(ialpha, jalpha);
        Gpu::Atomic::Add(&(ucol[idx_vec]), tmpx);

        auto tmpy = vely(i,j+jloc_shift,k) * alpha_v(ialpha, jalpha);
        Gpu::Atomic::Add(&(vcol[idx_vec]), tmpy);

        auto tmpt = state(i,j,k,RhoTheta_comp) / state(i,j,k,Rho_comp)
                  * alpha_theta(ialpha, jalpha);
        Gpu::Atomic::Add(&(thetacol[idx_vec]), tmpt);
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
    amrex::Vector<size_t> start_t {putloc};
    amrex::Vector<size_t> count_t {1};
    ncf.var("times").put(&cumtime, start_t, count_t);

    // T flux
    // TODO: Make this the actual flux rather than just a placeholder
    Real Tflux = 0.0;
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
