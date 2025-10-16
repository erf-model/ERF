#include "ERF_NCWpsFile.H"
#include "AMReX_FArrayBox.H"
#include "ERF_DataStruct.H"

using namespace amrex;

#ifdef ERF_USE_NETCDF
Box
read_subdomain_from_wrfinput(int /*lev*/, const std::string& fname, int& ratio)
{
    int  is, js;
    int  nx, ny, nz;
    if (ParallelDescriptor::IOProcessor()) {
        auto ncf = ncutils::NCFile::open(fname, NC_CLOBBER | NC_NETCDF4);
        { // Global Attributes (int)
            std::vector<int> attr;
            ncf.get_attr("WEST-EAST_GRID_DIMENSION"  , attr); nx = attr[0]-1;
            ncf.get_attr("SOUTH-NORTH_GRID_DIMENSION", attr); ny = attr[0]-1;
            ncf.get_attr("BOTTOM-TOP_GRID_DIMENSION" , attr); nz = attr[0]-1;
            amrex::Print() << "Have read (nx,ny,nz) = " << nx << " " << ny << " " << nz << std::endl;
            ncf.get_attr("I_PARENT_START", attr)   ; is    = attr[0]-1;
            ncf.get_attr("J_PARENT_START", attr)   ; js    = attr[0]-1;
            ncf.get_attr("PARENT_GRID_RATIO", attr); ratio = attr[0];
        }
        ncf.close();
        amrex::Print() << "Have read (parent_ilo,parent_jlo) = " << is << " " << js << std::endl;
        amrex::Print() << "Have read refinement ratio        = " << ratio << std::endl;
    }

    amrex::ParallelDescriptor::Bcast(&is   ,1,amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::Bcast(&js   ,1,amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::Bcast(&nx   ,1,amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::Bcast(&ny   ,1,amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::Bcast(&nz   ,1,amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::Bcast(&ratio,1,amrex::ParallelDescriptor::IOProcessorNumber());

    return Box( IntVect(ratio*is,ratio*js,0), IntVect(ratio*is+nx-1, ratio*js+ny-1, nz-1) );
}

void
read_from_wrfinput (int lev,
                    const Box& subdomain,
                    const std::string& fname,
                    FArrayBox& NC_fab,
                    const std::string& NC_name,
                    Geometry& geom,
                    int& use_theta_m,
                    int& success)
{
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
            ncf.get_attr("USE_THETA_M", attr);                use_theta_m = attr[0];
        }
        { // Global Attributes (Real)
            std::vector<Real> attr;
            ncf.get_attr("DX", attr); NC_dx = attr[0];
            ncf.get_attr("DY", attr); NC_dy = attr[0];
        }
        { // Global Attributes (string)
            NC_dateTime = ncf.get_attr("SIMULATION_START_DATE");
            const std::string dateTimeFormat = "%Y-%m-%d_%H:%M:%S";
            NC_epochTime = getEpochTime(NC_dateTime, dateTimeFormat);
        }
        int west_east_stag   = static_cast<int>(ncf.dim("west_east_stag").len());
        int south_north_stag = static_cast<int>(ncf.dim("south_north_stag").len());
        int bottom_top       = static_cast<int>(ncf.dim("bottom_top").len());
        ncf.close();

        amrex::ignore_unused(NC_dateTime); amrex::ignore_unused(NC_epochTime);

        // Verify the inputs geometry matches what the NETCDF file has
        if (lev == 0) {
            AMREX_ASSERT(west_east_stag == NC_nx);
            AMREX_ASSERT(south_north_stag == NC_ny);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(west_east_stag == geom.Domain().length(0)+1,
                    "amr.n_cell[0] does not match wrfinput");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(south_north_stag == geom.Domain().length(1)+1,
                    "amr.n_cell[1] does not match wrfinput");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bottom_top >= geom.Domain().length(2),
                    "amr.n_cell[2] does not match wrfinput");
            if (geom.Domain().length(2) < bottom_top) {
                Warning("Fewer vertical levels in ERF domain than the wrfinput");
            }

            Real rtol  = 1.0e-7;
            Real Len_x = NC_dx * Real(NC_nx-1);
            Real Len_y = NC_dy * Real(NC_ny-1);
            if (std::fabs((Len_x - (geom.ProbHi(0) - geom.ProbLo(0))) / Len_x) > rtol) {
                Print() << "X problem extent " << (geom.ProbHi(0) - geom.ProbLo(0)) << " does not match NETCDF file "
                            << Len_x << "!\n";
                Print() << "dx: " << NC_dx << ' ' << "Nx: " << NC_nx-1 << "\n";
                Abort("Domain specification error");
            }
            if (std::fabs((Len_y - (geom.ProbHi(1) - geom.ProbLo(1))) / Len_y) > rtol) {
                Print() << "Y problem extent " << (geom.ProbHi(1) - geom.ProbLo(1)) << " does not match NETCDF file "
                        << Len_y << "!\n";
                Print() << "dy: " << NC_dy << ' ' << "Ny: " << NC_ny-1 << "\n";
                Abort("Domain specification error");
            }
        } // lev == 0
    } // IOProc

    Vector<FArrayBox*> NC_fabs; NC_fabs.push_back(&NC_fab);
    Vector<std::string> NC_names; NC_names.push_back(NC_name);
    Vector<enum NC_Data_Dims_Type> NC_dim_types;
    Vector<int> successes; successes.resize(NC_names.size());

    if (NC_name == "ALB" || NC_name == "AL" || NC_name == "U" ||  NC_name == "V" ||  NC_name == "W" ||
        NC_name == "THM" || NC_name == "PH" || NC_name == "PHB" || NC_name == "PB" || NC_name == "P" ||
        NC_name == "QVAPOR"   || NC_name == "QCLOUD" || NC_name == "QRAIN")
    {
        // Note: staggering is handled in `fill_fab_from_arrays`
        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    }
    else if (NC_name == "MAPFAC_U" || NC_name == "MAPFAC_V" || NC_name == "MAPFAC_M" ||
             NC_name == "MUB"      || NC_name == "SST"      || NC_name == "LANDMASK"  ||
             NC_name == "XLAT_V"   || NC_name == "XLONG_U"  || NC_name == "TSK" ||
             NC_name == "PSFC"     || NC_name == "IVGTYP"   || NC_name == "ISLTYP" ||
             NC_name == "LAI"      || NC_name == "VEGFRA"   || NC_name == "TMN" ||
             NC_name == "SHDMIN"   || NC_name == "SHDMAX")
    {
        // Note: staggering is handled in `fill_fab_from_arrays`
        NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);
    }
    else if (NC_name == "C1H" || NC_name == "C2H")
    {
        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT);
    }
    else if (NC_name == "TSLB" || NC_name == "SMOIS" || NC_name == "SH2O")
    {
        NC_dim_types.push_back(NC_Data_Dims_Type::Time_SL_SN_WE);
    }
    else if (NC_name == "ZS" || NC_name == "DZS") {
        NC_dim_types.push_back(NC_Data_Dims_Type::Time_SL);
    }
    else {
        amrex::Print() << " ERROR: no NC dim type for NC_name = '" << NC_name << "'" << std::endl;
    }

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(subdomain, fname, NC_names, NC_dim_types, NC_fabs, successes);

    // Success was already broadcast in ERF_NCWpsFile.H
    success = successes[0];

    // Broadcast use_theta_m
    ParallelDescriptor::Bcast(&use_theta_m, 1, ParallelDescriptor::IOProcessorNumber());
}
#endif // ERF_USE_NETCDF
