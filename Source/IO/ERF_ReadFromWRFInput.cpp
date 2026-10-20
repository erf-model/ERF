/**
 * \file ERF_ReadFromWRFInput.cpp
 */
#include "ERF_NCWpsFile.H"
#include "AMReX_FArrayBox.H"
#include "ERF_DataStruct.H"

using namespace amrex;

#ifdef ERF_USE_NETCDF
Box
read_subdomain_from_wrfinput (int /*lev*/, const std::string& fname, int& ratio)
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

/**
 * The layout in which a wrfinput variable is stored in the file.
 *
 * @param NC_name Variable name
 */
static NC_Data_Dims_Type
wrfinput_dim_type (const std::string& NC_name)
{
    if (NC_name == "ALB" || NC_name == "AL" || NC_name == "ALT" ||
        NC_name == "U" ||  NC_name == "V" ||  NC_name == "W" ||
        NC_name == "THM" || NC_name == "PH" || NC_name == "PHB" || NC_name == "PB" || NC_name == "P" ||
        NC_name == "QVAPOR"   || NC_name == "QCLOUD"  || NC_name == "QICE" ||
        NC_name == "QRAIN"    || NC_name == "QSNOW"   || NC_name == "QGRAUP")
    {
        // Note: staggering is handled in `fill_fab_from_arrays`
        return NC_Data_Dims_Type::Time_BT_SN_WE;
    }
    else if (NC_name == "MAPFAC_U" || NC_name == "MAPFAC_V" || NC_name == "MAPFAC_M" ||
             NC_name == "MUB"      || NC_name == "SST"      || NC_name == "LANDMASK"  ||
             NC_name == "XLAT"     || NC_name == "XLONG"    || NC_name == "TSK" ||
             NC_name == "PSFC"     || NC_name == "IVGTYP"   || NC_name == "ISLTYP" ||
             NC_name == "LAI"      || NC_name == "VEGFRA"   || NC_name == "TMN" ||
             NC_name == "SHDMIN"   || NC_name == "SHDMAX")
    {
        // Note: staggering is handled in `fill_fab_from_arrays`
        return NC_Data_Dims_Type::Time_SN_WE;
    }
    else if (NC_name == "C1H" || NC_name == "C2H" || NC_name == "RDNW")
    {
        return NC_Data_Dims_Type::Time_BT;
    }
    else if (NC_name == "TSLB" || NC_name == "SMOIS" || NC_name == "SH2O")
    {
        return NC_Data_Dims_Type::Time_SL_SN_WE;
    }
    else if (NC_name == "ZS" || NC_name == "DZS") {
        return NC_Data_Dims_Type::Time_SL;
    }
    amrex::Abort("wrfinput_dim_type: no NC dim type for NC_name = '" + NC_name + "'");
    return NC_Data_Dims_Type::Time;
}

/**
 * Check that the wrfinput file matches the problem geometry (at level 0) and read
 * USE_THETA_M from it.
 *
 * @param lev Level of the file
 * @param fname Name of the wrfinput file
 * @param geom Geometry of the level
 * @param use_theta_m Set on all ranks from the file's USE_THETA_M attribute
 */
void
check_wrfinput_geometry (int lev,
                         const std::string& fname,
                         const Geometry& geom,
                         int& use_theta_m)
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

            Real rtol  = Real(1.0e-7);
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

    ParallelDescriptor::Bcast(&use_theta_m, 1, ParallelDescriptor::IOProcessorNumber());
}

void
read_from_wrfinput (int lev,
                    const Box& subdomain_to_read,
                    const std::string& fname,
                    FArrayBox& NC_fab_to_read,
                    const std::string& NC_name,
                    Geometry& geom,
                    int& use_theta_m,
                    int& success)
{
    check_wrfinput_geometry(lev, fname, geom, use_theta_m);

    Vector<FArrayBox*> NC_fabs_to_read; NC_fabs_to_read.push_back(&NC_fab_to_read);
    Vector<std::string> NC_names; NC_names.push_back(NC_name);
    Vector<enum NC_Data_Dims_Type> NC_dim_types; NC_dim_types.push_back(wrfinput_dim_type(NC_name));
    Vector<int> successes; successes.resize(NC_names.size());

    // Read the netcdf file and fill these FABs.  The fab spans the whole file domain on
    // every rank and every caller only copies out of it, so keep it in (pinned) host
    // memory rather than allocating a full-domain copy on the device.
    BuildFABsFromNetCDFFile<FArrayBox,Real>(subdomain_to_read, fname, NC_names, NC_dim_types, NC_fabs_to_read, successes,
                                            /*keep_on_host=*/ true);

    // Success was already broadcast in ERF_NCWpsFile.H
    success = successes[0];
}

/**
 * Find the box that read_from_wrfinput would fill for a variable, from the file's
 * metadata alone; no field data is read.  The box is in the same index space, and has
 * the same staggering, as the fab that read_from_wrfinput returns.
 *
 * @param subdomain_to_read Box specifying the region covered by the file (from read_subdomain_from_wrfinput)
 * @param fname Name of the wrfinput file
 * @param NC_name Variable name
 * @param var_box Box of the variable, set on all ranks when the variable is present
 * @return 1 if the file holds the variable, else 0, on all ranks
 */
int
wrfinput_var_box (const Box& subdomain_to_read,
                  const std::string& fname,
                  const std::string& NC_name,
                  Box& var_box)
{
    const int ioproc = ParallelDescriptor::IOProcessorNumber();
    int success = 0;
    if (ParallelDescriptor::IOProcessor()) {
        auto ncf = ncutils::NCFile::open(fname, NC_NOWRITE);
        success = ncf.has_var(NC_name);
        if (success) {
            std::vector<size_t> shape = ncf.var(NC_name).shape();

            // This mirrors the box built by fill_fab_from_arrays
            NC_Data_Dims_Type dim_type = wrfinput_dim_type(NC_name);
            int ns1 = 1, ns2 = 1, ns3 = 1;
            if (dim_type == NC_Data_Dims_Type::Time_BT || dim_type == NC_Data_Dims_Type::Time_SL) {
                ns1 = static_cast<int>(shape[1]);
            } else if (dim_type == NC_Data_Dims_Type::Time_SN_WE) {
                ns2 = static_cast<int>(shape[1]);
                ns3 = static_cast<int>(shape[2]);
            } else {
                ns1 = static_cast<int>(shape[1]);
                ns2 = static_cast<int>(shape[2]);
                ns3 = static_cast<int>(shape[3]);
            }
            int khi = std::min(subdomain_to_read.bigEnd(2), ns1-1);

            var_box = Box(IntVect(0,0,0), IntVect(ns3-1,ns2-1,khi));
            if (NC_name == "PH" || NC_name == "PHB" || NC_name == "W") {
                var_box.setType(IndexType(IntVect(0,0,1)));
                var_box.setBig(2, khi+1);
            } else if (NC_name == "U" || NC_name == "MAPFAC_U") {
                var_box.setType(IndexType(IntVect(1,0,0)));
            } else if (NC_name == "V" || NC_name == "MAPFAC_V") {
                var_box.setType(IndexType(IntVect(0,1,0)));
            } else if (NC_name == "TSLB" || NC_name == "SMOIS" || NC_name == "SH2O" ||
                       NC_name == "ZS"   || NC_name == "DZS") {
                // soil variables are staggered in z
                var_box.setType(IndexType(IntVect(0,0,1)));
            }
            var_box.shift(subdomain_to_read.smallEnd());
        }
        ncf.close();
    }
    ParallelDescriptor::Bcast(&success, 1, ioproc);
    if (success) {
        ParallelDescriptor::Bcast(&var_box, 1, ioproc);
    }
    return success;
}

/**
 * Read one region of a wrfinput variable on every rank.  Each rank reads only the part
 * of the file that its own region covers, so no rank ever holds the whole field.
 *
 * Every rank opens the file itself.  So that a large job does not have every rank
 * hitting the file at once, the ranks read in turns, at most max_readers at a time.
 *
 * @param fname Name of the wrfinput file
 * @param NC_name Variable name
 * @param var_box Box of the whole variable, as returned by wrfinput_var_box
 * @param region Part of var_box (same index type) this rank wants; may be empty
 * @param fab Fab we are to fill; on GPU builds it lives in pinned host memory
 * @param max_readers Most ranks reading at once (<= 0 for no limit)
 */
void
read_region_from_wrfinput (const std::string& fname,
                           const std::string& NC_name,
                           const Box& var_box,
                           const Box& region,
                           FArrayBox& fab,
                           int max_readers)
{
    AMREX_ALWAYS_ASSERT(region.ixType() == var_box.ixType());
    const Box bx = region & var_box;

    const int nprocs = ParallelDescriptor::NProcs();
    const int myproc = ParallelDescriptor::MyProc();
    const int nwaves = (max_readers > 0) ? (nprocs + max_readers - 1) / max_readers : 1;

    for (int wave = 0; wave < nwaves; ++wave) {
        if ((myproc % nwaves) == wave && bx.ok()) {
            // Offsets into the file, which indexes from zero at var_box's low corner
            const int ilo = bx.smallEnd(0) - var_box.smallEnd(0);
            const int jlo = bx.smallEnd(1) - var_box.smallEnd(1);
            const int klo = bx.smallEnd(2) - var_box.smallEnd(2);
            const int ni  = bx.length(0);
            const int nj  = bx.length(1);
            const int nk  = bx.length(2);

            // Only the first time is ever used
            std::vector<size_t> start, count;
            NC_Data_Dims_Type dim_type = wrfinput_dim_type(NC_name);
            if (dim_type == NC_Data_Dims_Type::Time_BT || dim_type == NC_Data_Dims_Type::Time_SL) {
                start = {0, size_t(klo)};
                count = {1, size_t(nk)};
            } else if (dim_type == NC_Data_Dims_Type::Time_SN_WE) {
                start = {0, size_t(jlo), size_t(ilo)};
                count = {1, size_t(nj), size_t(ni)};
            } else {
                start = {0, size_t(klo), size_t(jlo), size_t(ilo)};
                count = {1, size_t(nk), size_t(nj), size_t(ni)};
            }

            std::vector<float> buf(bx.numPts());
            {
                auto prev_fpe_except = amrex::disableFPExcept(amrex::FPExcept::all);
                auto ncf = ncutils::NCFile::open(fname, NC_NOWRITE);
                ncf.var(NC_name).get(buf.data(), start, count);
                ncf.close();
                std::feclearexcept(FE_ALL_EXCEPT);
                amrex::setFPExcept(prev_fpe_except);
            }

#ifdef AMREX_USE_GPU
            fab.resize(bx, 1, The_Pinned_Arena());
#else
            fab.resize(bx, 1);
#endif
            const Array4<Real>& arr = fab.array();
            const Dim3 lo = lbound(bx);
            for (int k = 0; k < nk; ++k) {
                for (int j = 0; j < nj; ++j) {
                    for (int i = 0; i < ni; ++i) {
                        arr(lo.x+i,lo.y+j,lo.z+k) = static_cast<Real>(buf[(size_t(k)*nj + j)*ni + i]);
                    }
                }
            }
        }
        if (nwaves > 1) { ParallelDescriptor::Barrier(); }
    }
}
#endif // ERF_USE_NETCDF
