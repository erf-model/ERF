#include <ERF_NCWpsFile.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <ERF_MetgridUtils.H>

using namespace amrex;

#ifdef ERF_USE_NETCDF

Box
read_subdomain_from_metgrid(int /*lev*/, const std::string& fname, int& ratio, int& klo, int& khi)
{
    int is, js;
    int nx, ny, nz;
    if (ParallelDescriptor::IOProcessor()) {
        auto ncf = ncutils::NCFile::open(fname, NC_CLOBBER | NC_NETCDF4);
        { // Global Attributes (int)
            std::vector<int> attr;
            ncf.get_attr("WEST-EAST_GRID_DIMENSION"  , attr); nx = attr[0]-1;
            ncf.get_attr("SOUTH-NORTH_GRID_DIMENSION", attr); ny = attr[0]-1;
            ncf.get_attr("BOTTOM-TOP_GRID_DIMENSION" , attr); nz = attr[0]-1;
#ifndef AMREX_USE_GPU
            Print() << "Have read (nx,ny,nz) = " << nx << " " << ny << " " << nz << std::endl;
#endif
            ncf.get_attr("i_parent_start", attr)   ; is    = attr[0]-1;
            ncf.get_attr("j_parent_start", attr)   ; js    = attr[0]-1;
            ncf.get_attr("parent_grid_ratio", attr); ratio = attr[0];
        }
        ncf.close();

#ifndef AMREX_USE_GPU
        amrex::Print() << "Have read (parent_ilo,parent_jlo) = " << is << " " << js << std::endl;
        amrex::Print() << "Have read refinement ratio        = " << ratio << std::endl;
#endif

    }
    ParallelDescriptor::Bcast(&is   ,1,amrex::ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&js   ,1,amrex::ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&nx   ,1,amrex::ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&ny   ,1,amrex::ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&nz   ,1,amrex::ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&ratio,1,amrex::ParallelDescriptor::IOProcessorNumber());

    return Box( IntVect(ratio*is, ratio*js, klo), IntVect(ratio*is+nx-1, ratio*js+ny-1, khi) );
}

void
read_from_metgrid (int lev, int itime,
                   const Box& domain, const std::string& fname,
                   std::string& NC_dateTime, Real& NC_epochTime,
                   int& flag_psfc, int& flag_msf,
                   int& flag_sst,  int& flag_tsk, int& flag_lmask,
                   int& NC_nx,     int& NC_ny,
                   Real& NC_dx,    Real& NC_dy,
                   FArrayBox& NC_xvel_fab, FArrayBox& NC_yvel_fab,
                   FArrayBox& NC_temp_fab, FArrayBox& NC_rhum_fab,
                   FArrayBox& NC_pres_fab, FArrayBox& NC_ght_fab,
                   FArrayBox& NC_hgt_fab,  FArrayBox& NC_psfc_fab,
                   FArrayBox& NC_msfu_fab, FArrayBox& NC_msfv_fab,
                   FArrayBox& NC_msfm_fab,
                   FArrayBox& NC_sst_fab,  FArrayBox& NC_tsk_fab,
                   FArrayBox& NC_LAT_fab,  FArrayBox& NC_LON_fab,
                   IArrayBox& NC_lmask_iab,
                   Geometry& geom)
{
#ifndef AMREX_USE_GPU
    Print() << "Loading header data from NetCDF file at level " << lev << std::endl;
#endif

    if (ParallelDescriptor::IOProcessor()) {
        auto ncf = ncutils::NCFile::open(fname, NC_CLOBBER | NC_NETCDF4);

        { // Global Attributes (int)
            std::vector<int> attr;
            ncf.get_attr("WEST-EAST_GRID_DIMENSION", attr);   NC_nx = attr[0];
            ncf.get_attr("SOUTH-NORTH_GRID_DIMENSION", attr); NC_ny = attr[0];
        }

        { // Global Attributes (string)
            NC_dateTime = ncf.get_attr("SIMULATION_START_DATE");
            const std::string dateTimeFormat = "%Y-%m-%d_%H:%M:%S";
            NC_epochTime = getEpochTime(NC_dateTime, dateTimeFormat);
        }

        { // Global Attributes (Real)
            std::vector<Real> attr;
            ncf.get_attr("DX", attr); NC_dx = attr[0];
            ncf.get_attr("DY", attr); NC_dy = attr[0];
        }
        ncf.close();

        // Verify the inputs geometry matches what the NETCDF file has
        if (lev == 0) {
            Real tol   = 1.0e-3;
            Real Len_x = NC_dx * Real(NC_nx-1);
            Real Len_y = NC_dy * Real(NC_ny-1);
            if (std::fabs(Len_x - (geom.ProbHi(0) - geom.ProbLo(0))) > tol) {
#ifndef AMREX_USE_GPU
                Print() << "X problem extent " << (geom.ProbHi(0) - geom.ProbLo(0)) << " does not match NETCDF file "
                        << Len_x << "!\n";
                Print() << "dx: " << NC_dx << ' ' << "Nx: " << NC_nx-1 << "\n";
#endif
                Abort("Domain specification error");
            }
            if (std::fabs(Len_y - (geom.ProbHi(1) - geom.ProbLo(1))) > tol) {
#ifndef AMREX_USE_GPU
                Print() << "Y problem extent " << (geom.ProbHi(1) - geom.ProbLo(1)) << " does not match NETCDF file "
                        << Len_y << "!\n";
                Print() << "dy: " << NC_dy << ' ' << "Ny: " << NC_ny-1 << "\n";
#endif
                Abort("Domain specification error");
            }
        } // lev == 0
    } // IOProc

    int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank
    ParallelDescriptor::Bcast(&NC_nx,        1, ioproc);
    ParallelDescriptor::Bcast(&NC_ny,        1, ioproc);
    ParallelDescriptor::Bcast(&NC_epochTime, 1, ioproc);
    ParallelDescriptor::Bcast(&NC_dx,        1, ioproc);
    ParallelDescriptor::Bcast(&NC_dy,        1, ioproc);

#ifndef AMREX_USE_GPU
    Print() << "Loading initial data from NetCDF file at level " << lev << std::endl;
#endif

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
    NC_fabs.push_back(&NC_LAT_fab);       NC_fnames.push_back("XLAT_M");    NC_fdim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);
    NC_fabs.push_back(&NC_LON_fab);       NC_fnames.push_back("XLONG_M");   NC_fdim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);
    NC_fabs.push_back(&NC_hgt_fab);       NC_fnames.push_back("HGT_M");     NC_fdim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);

    if (itime == 0) {
        NC_fabs.push_back(&NC_psfc_fab);  NC_fnames.push_back("PSFC");      NC_fdim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);
    }
    NC_fabs.push_back(&NC_msfu_fab);      NC_fnames.push_back("MAPFAC_U");  NC_fdim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);
    NC_fabs.push_back(&NC_msfv_fab);      NC_fnames.push_back("MAPFAC_V");  NC_fdim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);
    NC_fabs.push_back(&NC_msfm_fab);      NC_fnames.push_back("MAPFAC_M");  NC_fdim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);
    NC_fabs.push_back(&NC_sst_fab);       NC_fnames.push_back("SST");       NC_fdim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);
    NC_fabs.push_back(&NC_tsk_fab);       NC_fnames.push_back("SKINTEMP");  NC_fdim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);

    NC_iabs.push_back(&NC_lmask_iab);     NC_inames.push_back("LANDMASK");   NC_idim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);

    // Read the netcdf file and fill these FABs
#ifndef AMREX_USE_GPU
    Print() << "Building initial FABS from file " << fname << std::endl;
#endif

    Vector<int> success; success.resize(NC_fabs.size());
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_fnames, NC_fdim_types, NC_fabs, success);
    for (int i = 0; i < success.size(); i++) {
        flag_psfc = (NC_fnames[i] == "PSFC"     && success[i] == 1) ? 1 : 0;
        flag_sst  = (NC_fnames[i] == "SST"      && success[i] == 1) ? 1 : 0;
        flag_tsk  = (NC_fnames[i] == "SKINTEMP" && success[i] == 1) ? 1 : 0;
        flag_msf  = (NC_fnames[i] == "MAPFAC_M" && success[i] == 1) ? 1 : 0;
    }

    // Read the netcdf file and fill these IABs
#ifndef AMREX_USE_GPU
    Print() << "Building initial IABS from file " << fname << std::endl;
#endif
    Vector<int> success_i; success_i.resize(NC_iabs.size());
    BuildFABsFromNetCDFFile<IArrayBox,int>(domain, fname, NC_inames, NC_idim_types, NC_iabs, success_i);
    for (int i = 0; i < success_i.size(); i++) {
        flag_lmask = (NC_inames[i] == "LANDMASK" && success_i[i] == 1) ? 1 : 0;
    }

#ifndef AMREX_USE_GPU
    Print() << " flag_psfc:  " << flag_psfc << std::endl;
    Print() << " flag_sst:   " << flag_sst << std::endl;
    Print() << " flag_tsk:   " << flag_tsk << std::endl;
    Print() << " flag_msf:   " << flag_msf << std::endl;
    Print() << " flag_lmask: " << flag_lmask << std::endl;
#endif

    // TODO: FIND OUT IF WE NEED TO DIVIDE VELS BY MAPFAC
    //
    // Convert the velocities using the map factors
    //
    //const Box& uubx = NC_xvel_fab.box();
    // const Array4<Real>    u_arr = NC_xvel_fab.array();
    // const Array4<Real> msfu_arr = NC_msfu_fab.array();
    //ParallelFor(uubx, [=] AMREX_GPU_DEVICE (int , int , int )
    //{
    //    u_arr(i,j,k) /= msfu_arr(i,j,0);
    //});

    //const Box& vvbx = NC_yvel_fab.box();
    // const Array4<Real>    v_arr = NC_yvel_fab.array();
    // const Array4<Real> msfv_arr = NC_msfv_fab.array();
    //ParallelFor(vvbx, [=] AMREX_GPU_DEVICE (int , int , int )
    //{
    //    v_arr(i,j,k) /= msfv_arr(i,j,0);
    //});
}
#endif // ERF_USE_NETCDF
