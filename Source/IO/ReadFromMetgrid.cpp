#include <NCWpsFile.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <Metgrid_utils.H>

using namespace amrex;

#ifdef ERF_USE_NETCDF

void
read_from_metgrid (int lev, const Box& domain, const std::string& fname,
                   std::string& NC_dateTime, Real& NC_epochTime,
                   int& flag_psfc, int& flag_msfu, int& flag_msfv,  int& flag_msfm,
                   int& flag_hgt,  int& flag_sst,  int& flag_lmask,
                   int& NC_nx,     int& NC_ny,
                   Real& NC_dx,    Real& NC_dy,
                   FArrayBox& NC_xvel_fab, FArrayBox& NC_yvel_fab,
                   FArrayBox& NC_temp_fab, FArrayBox& NC_rhum_fab,
                   FArrayBox& NC_pres_fab, FArrayBox& NC_ght_fab,
                   FArrayBox& NC_hgt_fab,  FArrayBox& NC_psfc_fab,
                   FArrayBox& NC_msfu_fab, FArrayBox& NC_msfv_fab,
                   FArrayBox& NC_msfm_fab, FArrayBox& NC_sst_fab,
                   FArrayBox& NC_LAT_fab,  FArrayBox& NC_LON_fab,
                   IArrayBox& NC_lmask_iab,
                   Real& Latitude,
                   Real& Longitude,
                   Geometry& geom)
{
    Print() << "Loading header data from NetCDF file at level " << lev << std::endl;

    if (ParallelDescriptor::IOProcessor()) {
        auto ncf = ncutils::NCFile::open(fname, NC_CLOBBER | NC_NETCDF4);
        { // Global Attributes (int)
            std::vector<int> attr;
            ncf.get_attr("FLAG_PSFC", attr);     flag_psfc  = attr[0];

            /* UNCOMMENT FOR FLOWMAS
            flag_msfu  = 0; //ncf.get_attr("FLAG_MAPFAC_U", attr); flag_msfu  = attr[0];
            flag_msfv  = 0; //ncf.get_attr("FLAG_MAPFAC_V", attr); flag_msfv  = attr[0];
            flag_msfm = 0; //ncf.get_attr("FLAG_MAPFAC_M", attr); flag_msfm  = attr[0];
            flag_hgt = 1; //ncf.get_attr("FLAG_HGT_M", attr);    flag_hgt   = attr[0];
            flag_sst = 0; //ncf.get_attr("FLAG_SST", attr);      flag_sst   = attr[0];
            flag_lmask = 0; //ncf.get_attr("FLAG_LANDMASK", attr); flag_lmask = attr[0];
            */

            ncf.get_attr("FLAG_MAPFAC_U", attr); flag_msfu  = attr[0];
            ncf.get_attr("FLAG_MAPFAC_V", attr); flag_msfv  = attr[0];
            ncf.get_attr("FLAG_MAPFAC_M", attr); flag_msfm  = attr[0];
            ncf.get_attr("FLAG_HGT_M", attr);    flag_hgt   = attr[0];
            ncf.get_attr("FLAG_SST", attr);      flag_sst   = attr[0];
            ncf.get_attr("FLAG_LANDMASK", attr); flag_lmask = attr[0];


            ncf.get_attr("WEST-EAST_GRID_DIMENSION", attr);   NC_nx = attr[0];
            ncf.get_attr("SOUTH-NORTH_GRID_DIMENSION", attr); NC_ny = attr[0];
        }
        { // Global Attributes (string)
            NC_dateTime = ncf.get_attr("SIMULATION_START_DATE")+"UTC";
            const std::string dateTimeFormat = "%Y-%m-%d_%H:%M:%S%Z";
            NC_epochTime = getEpochTime(NC_dateTime, dateTimeFormat);
        }
        { // Global Attributes (Real)
            std::vector<Real> attr;
            ncf.get_attr("DX", attr); NC_dx = attr[0];
            ncf.get_attr("DY", attr); NC_dy = attr[0];
        }
        ncf.close();

        // Verify the inputs geometry matches what the NETCDF file has
        Real tol   = 1.0e-3;
        Real Len_x = NC_dx * Real(NC_nx-1);
        Real Len_y = NC_dy * Real(NC_ny-1);
        if (std::fabs(Len_x - (geom.ProbHi(0) - geom.ProbLo(0))) > tol) {
            Print() << "X problem extent " << (geom.ProbHi(0) - geom.ProbLo(0)) << " does not match NETCDF file "
                    << Len_x << "!\n";
            Print() << "dx: " << NC_dx << ' ' << "Nx: " << NC_nx-1 << "\n";
            Abort("Domain specification error");
        }
        if (std::fabs(Len_y - (geom.ProbHi(1) - geom.ProbLo(1))) > tol) {
            Print() << "Y problem extent " << (geom.ProbHi(1) - geom.ProbLo(1)) << " does not match NETCDF file "
                    << Len_y << "!\n";
            Print() << "dy: " << NC_dy << ' ' << "Ny: " << NC_ny-1 << "\n";
            Abort("Domain specification error");
        }
    }
    int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank
    ParallelDescriptor::Bcast(&flag_psfc,    1, ioproc);
    ParallelDescriptor::Bcast(&flag_msfu,    1, ioproc);
    ParallelDescriptor::Bcast(&flag_msfv,    1, ioproc);
    ParallelDescriptor::Bcast(&flag_msfm,    1, ioproc);
    ParallelDescriptor::Bcast(&flag_hgt,     1, ioproc);
    ParallelDescriptor::Bcast(&flag_sst,     1, ioproc);
    ParallelDescriptor::Bcast(&flag_lmask,   1, ioproc);
    ParallelDescriptor::Bcast(&NC_nx,        1, ioproc);
    ParallelDescriptor::Bcast(&NC_ny,        1, ioproc);
    ParallelDescriptor::Bcast(&NC_epochTime, 1, ioproc);
    ParallelDescriptor::Bcast(&NC_dx,        1, ioproc);
    ParallelDescriptor::Bcast(&NC_dy,        1, ioproc);

    Print() << "Loading initial data from NetCDF file at level " << lev << std::endl;

    Vector<FArrayBox*>  NC_fabs;
    Vector<IArrayBox*>  NC_iabs;
    Vector<std::string> NC_fnames;
    Vector<std::string> NC_inames;
    Vector<enum NC_Data_Dims_Type> NC_fdim_types;
    Vector<enum NC_Data_Dims_Type> NC_idim_types;

    NC_fabs.push_back(&NC_xvel_fab);      NC_fnames.push_back("UU");        NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    NC_fabs.push_back(&NC_yvel_fab);      NC_fnames.push_back("VV");        NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    NC_fabs.push_back(&NC_temp_fab);      NC_fnames.push_back("TT");        NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    NC_fabs.push_back(&NC_rhum_fab);      NC_fnames.push_back("RH");        NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    NC_fabs.push_back(&NC_pres_fab);      NC_fnames.push_back("PRES");      NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    NC_fabs.push_back(&NC_ght_fab);       NC_fnames.push_back("GHT");       NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    NC_fabs.push_back(&NC_LAT_fab);       NC_fnames.push_back("XLAT_V");    NC_fdim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);
    NC_fabs.push_back(&NC_LON_fab);       NC_fnames.push_back("XLONG_U");   NC_fdim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);

    if (flag_psfc)  { NC_fabs.push_back(&NC_psfc_fab);      NC_fnames.push_back("PSFC");      NC_fdim_types.push_back(NC_Data_Dims_Type::Time_SN_WE); }
    if (flag_msfu)  { NC_fabs.push_back(&NC_msfu_fab);      NC_fnames.push_back("MAPFAC_U");  NC_fdim_types.push_back(NC_Data_Dims_Type::Time_SN_WE); }
    if (flag_msfv)  { NC_fabs.push_back(&NC_msfv_fab);      NC_fnames.push_back("MAPFAC_V");  NC_fdim_types.push_back(NC_Data_Dims_Type::Time_SN_WE); }
    if (flag_msfm)  { NC_fabs.push_back(&NC_msfm_fab);      NC_fnames.push_back("MAPFAC_M");  NC_fdim_types.push_back(NC_Data_Dims_Type::Time_SN_WE); }
    if (flag_hgt)   { NC_fabs.push_back(&NC_hgt_fab);       NC_fnames.push_back("HGT_M");     NC_fdim_types.push_back(NC_Data_Dims_Type::Time_SN_WE); }
    if (flag_sst)   { NC_fabs.push_back(&NC_sst_fab);       NC_fnames.push_back("SST");       NC_fdim_types.push_back(NC_Data_Dims_Type::Time_SN_WE); }

    if (flag_lmask) { NC_iabs.push_back(&NC_lmask_iab);     NC_inames.push_back("LANDMASK");   NC_idim_types.push_back(NC_Data_Dims_Type::Time_SN_WE); }

    // Read the netcdf file and fill these FABs
    std::string Lat_var_name = "XLAT_V";
    std::string Lon_var_name = "XLONG_U";
    Print() << "Building initial FABS from file " << fname << std::endl;
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, Latitude, Longitude,
                                            Lat_var_name, Lon_var_name,
                                            fname, NC_fnames, NC_fdim_types, NC_fabs);

    // Read the netcdf file and fill these IABs
    Print() << "Building initial IABS from file " << fname << std::endl;
    BuildFABsFromNetCDFFile<IArrayBox,int>(domain, Latitude, Longitude,
                                           Lat_var_name, Lon_var_name,
                                           fname, NC_inames, NC_idim_types, NC_iabs);

    // TODO: FIND OUT IF WE NEED TO DIVIDE VELS BY MAPFAC
    //
    // Convert the velocities using the map factors
    //
    const Box& uubx = NC_xvel_fab.box();
    // const Array4<Real>    u_arr = NC_xvel_fab.array();
    // const Array4<Real> msfu_arr = NC_msfu_fab.array();
    ParallelFor(uubx, [=] AMREX_GPU_DEVICE (int , int , int )
    {
        // u_arr(i,j,k) /= msfu_arr(i,j,0);
    });

    const Box& vvbx = NC_yvel_fab.box();
    // const Array4<Real>    v_arr = NC_yvel_fab.array();
    // const Array4<Real> msfv_arr = NC_msfv_fab.array();
    ParallelFor(vvbx, [=] AMREX_GPU_DEVICE (int , int , int )
    {
        // v_arr(i,j,k) /= msfv_arr(i,j,0);
    });
}
#endif // ERF_USE_NETCDF
