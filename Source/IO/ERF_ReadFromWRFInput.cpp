#include "ERF_NCWpsFile.H"
#include "AMReX_FArrayBox.H"
#include "ERF_DataStruct.H"

using namespace amrex;

#ifdef ERF_USE_NETCDF
void
read_from_wrfinput (int lev,
                    const Box& domain,
                    const std::string& fname,
                    FArrayBox& NC_fab,
                    const std::string& NC_name,
                    Geometry& geom)
{
    Print() << "Loading header data from NetCDF file at level " << lev << std::endl;

    if (ParallelDescriptor::IOProcessor()) {
        int  NC_nx, NC_ny;
        Real NC_dx, NC_dy;
        Real NC_epochTime;
        std::string NC_dateTime;
        auto ncf = ncutils::NCFile::open(fname, NC_CLOBBER | NC_NETCDF4);
        { // Global Attributes (int)
            std::vector<int> attr;
            ncf.get_attr("WEST-EAST_GRID_DIMENSION", attr);   NC_nx = attr[0];
            ncf.get_attr("SOUTH-NORTH_GRID_DIMENSION", attr); NC_ny = attr[0];
        }
        { // Global Attributes (Real)
            std::vector<Real> attr;
            ncf.get_attr("DX", attr); NC_dx = attr[0];
            ncf.get_attr("DY", attr); NC_dy = attr[0];
        }
        { // Global Attributes (string)
            NC_dateTime = ncf.get_attr("SIMULATION_START_DATE")+"UTC";
            const std::string dateTimeFormat = "%Y-%m-%d_%H:%M:%S%Z";
            NC_epochTime = getEpochTime(NC_dateTime, dateTimeFormat);
        }
        ncf.close();

        amrex::ignore_unused(NC_dateTime); amrex::ignore_unused(NC_epochTime);

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

    Print() << "Loading initial data from NetCDF file at level " << lev << std::endl;

    Vector<FArrayBox*> NC_fabs; NC_fabs.push_back(&NC_fab);
    Vector<std::string> NC_names; NC_names.push_back(NC_name);
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    if (NC_name == "ALB")       { NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); } // 0
    if (NC_name == "AL")        { NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); } // 1
    if (NC_name == "U")         { NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); } // 2
    if (NC_name == "V")         { NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); } // 3
    if (NC_name == "W")         { NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); } // 4
    if (NC_name == "T")         { NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); } // 5
    if (NC_name == "PH")        { NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); } // 6
    if (NC_name == "PHB")       { NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); } // 7
    if (NC_name == "PB")        { NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); } // 8
    if (NC_name == "P")         { NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); } // 9
    if (NC_name == "MUB")       { NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);    } // 10
    if (NC_name == "MAPFAC_UY") { NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);    } // 11
    if (NC_name == "MAPFAC_VY") { NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);    } // 12
    if (NC_name == "MAPFAC_MY") { NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);    } // 13
    if (NC_name == "SST")       { NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);    } // 14
    if (NC_name == "LANDMASK")  { NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);    } // 15
    if (NC_name == "C1H")       { NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT);       } // 16
    if (NC_name == "C2H")       { NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT);       } // 17
    if (NC_name == "XLAT_V")    { NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);    } // 18
    if (NC_name == "XLONG_U")   { NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);    } // 19
    if (NC_name == "QVAPOR")    { NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); } // 20
    if (NC_name == "QCLOUD")    { NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); } // 21
    if (NC_name == "QRAIN")     { NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); } // 22

    amrex::Print() <<" NC_name " << NC_name << std::endl;
    amrex::Print() <<" NC_fabs " << NC_fabs[0]->box() << std::endl;

    amrex::Print() <<" SIZE OF NAMES " << NC_names.size() << std::endl;
    amrex::Print() <<" SIZE OF  DIMS " << NC_dim_types.size() << std::endl;

    // Read the netcdf file and fill these FABs
    Print() << "Building initial FABS from file " << fname << std::endl;
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
}
#endif // ERF_USE_NETCDF
