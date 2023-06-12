#include "NCWpsFile.H"
#include "AMReX_FArrayBox.H"

using namespace amrex;

#ifdef ERF_USE_NETCDF
void
read_from_metgrid(int lev, const std::string& fname,
                  std::string& NC_datetime,
                  int& flag_psfc, int& flag_msfu, int& flag_msfv, int& flag_msfm,
                  int& flag_hgt,  int& NC_nx,     int& NC_ny,
                  Real& NC_dx,    Real& NC_dy,
                  FArrayBox& NC_xvel_fab, FArrayBox& NC_yvel_fab,
                  FArrayBox& NC_temp_fab, FArrayBox& NC_rhum_fab,
                  FArrayBox& NC_pres_fab, FArrayBox& NC_ght_fab,
                  FArrayBox& NC_hgt_fab,  FArrayBox& NC_psfc_fab,
                  FArrayBox& NC_msfu_fab, FArrayBox& NC_msfv_fab,
                  FArrayBox& NC_msfm_fab)
{
    amrex::Print() << "Loading initial data from NetCDF file at level " << lev << std::endl;

    int ncomp = 1;

    if (amrex::ParallelDescriptor::IOProcessor()) {
        auto ncf = ncutils::NCFile::open(fname, NC_CLOBBER | NC_NETCDF4);
        { // Global Attributes (int)
            std::vector<int> attr;
            ncf.get_attr("FLAG_PSFC", attr); flag_psfc = attr[0];
            ncf.get_attr("FLAG_MAPFAC_U", attr); flag_msfu = attr[0];
            ncf.get_attr("FLAG_MAPFAC_V", attr); flag_msfv = attr[0];
            ncf.get_attr("FLAG_MAPFAC_M", attr); flag_msfm = attr[0];
            ncf.get_attr("FLAG_HGT_M", attr); flag_hgt  = attr[0];
            ncf.get_attr("WEST-EAST_GRID_DIMENSION", attr); NC_nx = attr[0];
            ncf.get_attr("SOUTH-NORTH_GRID_DIMENSION", attr); NC_ny = attr[0];
        }
        { // Global Attributes (string)
            NC_datetime = ncf.get_attr("SIMULATION_START_DATE");
        }
        { // Global Attributes (Real)
            std::vector<Real> attr;
            ncf.get_attr("DX", attr); NC_dx = attr[0];
            ncf.get_attr("DY", attr); NC_dy = attr[0];
        }
        ncf.close();
    }

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_xvel_fab);      NC_names.push_back("UU");        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    NC_fabs.push_back(&NC_yvel_fab);      NC_names.push_back("VV");        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    NC_fabs.push_back(&NC_temp_fab);      NC_names.push_back("TT");        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    NC_fabs.push_back(&NC_rhum_fab);      NC_names.push_back("RH");        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    NC_fabs.push_back(&NC_pres_fab);      NC_names.push_back("PRES");      NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    NC_fabs.push_back(&NC_ght_fab);       NC_names.push_back("GHT");       NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);

    NC_fabs.push_back(&NC_hgt_fab);       NC_names.push_back("HGT_M");     NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);
    NC_fabs.push_back(&NC_psfc_fab);      NC_names.push_back("PSFC");      NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);
    NC_fabs.push_back(&NC_msfu_fab);      NC_names.push_back("MAPFAC_U");  NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);
    NC_fabs.push_back(&NC_msfv_fab);      NC_names.push_back("MAPFAC_V");  NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);
    NC_fabs.push_back(&NC_msfm_fab);      NC_names.push_back("MAPFAC_M");  NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile(fname, NC_names, NC_dim_types, NC_fabs);

    // TODO: Read DX and DY from met_em file.

    // TODO: FIND OUT IF WE NEED TO DIVIDE VELS BY MAPFAC
    //
    // Convert the velocities using the map factors
    //
    const Box& uubx = NC_xvel_fab.box();
    const Array4<Real>    u_arr = NC_xvel_fab.array();
    const Array4<Real> msfu_arr = NC_msfu_fab.array();
    ParallelFor(uubx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        // u_arr(i,j,k) /= msfu_arr(i,j,0);
    });

    const Box& vvbx = NC_yvel_fab.box();
    const Array4<Real>    v_arr = NC_yvel_fab.array();
    const Array4<Real> msfv_arr = NC_msfv_fab.array();
    ParallelFor(vvbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        // v_arr(i,j,k) /= msfv_arr(i,j,0);
    });
}
#endif // ERF_USE_NETCDF
