#include "ERF_SLM.H"

#include <AMReX_PlotFileUtil.H>
#include "ERF_NCInterface.H"

using namespace amrex;

#ifdef ERF_USE_NETCDF
void
SLM::writeSLM_NetCDF(const MultiFab& mf, const Vector<std::string>& varnames, const amrex::Real time,
                     const std::string plot_prefix, const int level_step)
{
    std::string plotfilename = plot_prefix + ".nc";

    ncutils::NCFile ncf =
        (first_step)
        ? ncutils::NCFile::create_par(plotfilename, NC_CLOBBER | NC_NETCDF4)
        : ncutils::NCFile::open_par(plotfilename, NC_WRITE);

    if (first_step)
    {
        writeNCHeader(ncf, m_geom);

        // Write out 3D variables that are constant in time
        for (const int var : const_vars)
        {
            amrex::Print() << " Writing CONSTANT 3D SLM MF '" << LsmVarName_Full[var] << "'" << std::endl;
            writeMFtoNC(ncf, lsm_fab_vars[var].get(), LsmVarName_Full[var], -1.0);
        }

        // Write out 2D variables that are constant in time
        writeMFtoNC(ncf, &LAI, "LAI", -1.0);
        writeMFtoNC(ncf, &BAI, "BAI", -1.0);
        writeMFtoNC(ncf, &ztop, "ztop", -1.0);
    }

    // Write out SLM 3D fields
    for (int var = 0; var < LsmVar_SLM::NumVars; ++var) {
        if (const_vars.find(var) != const_vars.end()) continue;
        //writeMFtoNC(ncf, lsm_fab_vars[var].get(), LsmVarName_Full[var], time, true); // write ghost cells
        writeMFtoNC(ncf, lsm_fab_vars[var].get(), LsmVarName_Full[var], time, false); // write ghost cells
    }

    // Write out SLM 2D fields
    // TODO: fix this
    //  2D fields are combined into the single mf with varnames components
    AMREX_ALWAYS_ASSERT(mf.nComp() == varnames.size());
    IntVect ng(0, 0, 0);
    for (int var = 0; var < mf.nComp(); ++var)
    {
        MultiFab fab(mf, make_alias, var, 1);
        writeMFtoNC(ncf, &fab, varnames[var], time);
    }

    writeMFtoNC(ncf, &cp_vege, "cp_vege", time);
    writeMFtoNC(ncf, &z0_sfc, "z0_sfc", time);

    //writeMFtoNC(ncf, &Khai_L, "Khai_L", -1.0);

    //writeMFtoNC(ncf, &ustar, "ustar", time); // ustar and tstar already saved in 2D fab array above

    // additional slm outputs
    for (int var = 0; var < slm_diag.nComp(); ++var)
    {
        MultiFab fab(slm_diag, make_alias, var, 1);
        writeMFtoNC(ncf, &fab, diag_names[var], time);
    }

    // Save soil temperature and moisture diagnostic variables
    for (int var = 0; var < soilt_vars.nComp(); ++var)
    {
        MultiFab fab (soilt_vars, make_alias, var, 1);
        writeMFtoNC(ncf, &fab, soilt_var_names[var], time);
    }

    for (int var = 0; var < soilw_vars.nComp() - 1; ++var)
    {
        MultiFab fab (soilw_vars, make_alias, var, 1);
        writeMFtoNC(ncf, &fab, soilw_var_names[var], time);
    }

    ncf.close();
}

void
SLM::writeNCHeader(ncutils::NCFile &nc_file, const amrex::Geometry &geom)
{
    // define data dimensions: time, x, y, z
    nc_file.enter_def_mode();

    // create time as unlimited dimension
    nc_file.def_dim("time", NC_UNLIMITED);
    nc_file.def_var("time", ncutils::NCDType::Real, {"time"});

    auto domain = geom.Domain();
    const int nx = domain.length(0);
    const int ny = domain.length(1);
    const int nz = m_nz_lsm; // Z dim is defined using SLM geometry instead of ERF's

    // define dimensions for data variables
    nc_file.def_dim("x", nx);
    nc_file.def_dim("y", ny);
    nc_file.def_dim("z", nz);

    // define coordinates for each dimension
    nc_file.def_var("x", ncutils::NCDType::Real, {"x"});
    nc_file.def_var("y", ncutils::NCDType::Real, {"y"});
    nc_file.def_var("z", ncutils::NCDType::Real, {"z"});

    // enable collective mode for parallel NetCDF
    nc_file.var("time").par_access(NC_COLLECTIVE);
    nc_file.var("x").par_access(NC_COLLECTIVE);
    nc_file.var("y").par_access(NC_COLLECTIVE);
    nc_file.var("z").par_access(NC_COLLECTIVE);

    // write other metadata

    nc_file.exit_def_mode();

    // each rank writes its portion of x,y,z coordinates
    // TODO: is there a better way to do this?
    amrex::Arena* Arena_Used = amrex::The_Arena();
#ifdef AMREX_USE_GPU
    Arena_Used = amrex::The_Pinned_Arena();
#endif

    //   create table data in pinned space for device
    amrex::TableData<Real, 1> x({0}, {nx}, Arena_Used);
    amrex::TableData<Real, 1> y({0}, {ny}, Arena_Used);
    amrex::TableData<Real, 1> z({0}, {nz}, Arena_Used);

    auto x_arr = x.table();
    auto y_arr = y.table();
    auto z_arr = z.table();

    const auto prob_lo = geom.ProbLoArray();
    const auto prob_hi = geom.ProbHiArray();
    const auto dx = geom.CellSizeArray();
    const int d_khi_lsm = khi_lsm;
    for(MFIter mfi(*lsm_fab_vars[0]); mfi.isValid(); ++mfi)
    {
        const auto &box = mfi.validbox();
        const auto node_z_arr = lsm_fab_vars[LsmVar_SLM::node_z]->const_array(mfi);
        const int bx_offset = box.smallEnd(0);
        const int by_offset = box.smallEnd(1);
        ParallelFor(box.length(0), [=] AMREX_GPU_DEVICE (int i)
        {
            x_arr(i) = prob_lo[0] + (bx_offset+i+myhalf)*dx[0];
        });

        ParallelFor(box.length(1), [=] AMREX_GPU_DEVICE (int j)
        {
            y_arr(j) = prob_lo[1] + (by_offset+j+myhalf)*dx[1];
        });

        ParallelFor(nz, [=] AMREX_GPU_DEVICE (int k)
        {
            // For SLM, z is node_z
            z_arr(k) = node_z_arr(box.smallEnd(0), box.smallEnd(1), (k*-1)+d_khi_lsm);
        });

        amrex::Gpu::Device::synchronize();

        nc_file.var("x").put(x_arr.p, {static_cast<unsigned long>(box.smallEnd()[0])}, {static_cast<unsigned long>(box.length()[0])});
        nc_file.var("y").put(y_arr.p, {static_cast<unsigned long>(box.smallEnd()[1])}, {static_cast<unsigned long>(box.length()[1])});
        nc_file.var("z").put(z_arr.p, {0}, {static_cast<unsigned long>(box.length()[2])});
    }

    amrex::ParallelDescriptor::Barrier();
}

void SLM::writeMFtoNC(ncutils::NCFile &nc_file, const MultiFab* mf,
                      const std::string name, Real time, bool write_ghost)
{
    IntVect ngrow = mf->nGrowVect();

    bool var_2d = false;
    if (mf->boxArray().minimalBox().length(2) == 1) var_2d = true;

    int khi = 0;
    if (var_2d) {
        khi = mf->boxArray().minimalBox().bigEnd(2) - ngrow[2];
    }

    if (!nc_file.has_var(name))
    {
        nc_file.enter_def_mode();

        // create time dimension in file
        // TODO: fix this - move time dim to NCInterface
        if (time > -1.0 && !nc_file.has_dim("time"))
        {
            nc_file.def_dim("time", NC_UNLIMITED);
            nc_file.def_var("time", ncutils::NCDType::Real, {"time"});
        }

        IntVect box_size = mf->boxArray().minimalBox().size();
        if (write_ghost)
        {
            box_size += 2 * ngrow;
        }

        bool add_dims = true;
        // Check if file has x and y dimensions
        if (nc_file.has_dim("x") && nc_file.has_dim("y"))
        {
            // Check that the dimension bounds match
            size_t x_size = nc_file.dim("x").len();
            size_t y_size = nc_file.dim("y").len();

            if (x_size == box_size[0] && y_size == box_size[1])
            {
                // reuse these dimensions for this dataset
                add_dims = false;
            }
        } else {
            // if the x and y dimensions don't exist, create them
            nc_file.def_dim("x", box_size[0]);
            nc_file.def_dim("y", box_size[1]);
            if (!var_2d) nc_file.def_dim("z", box_size[2]);
            add_dims = false;
        }

        std::vector<std::string> dim_names;
        if (add_dims) {
            dim_names = {name + "_z", name + "_y", name + "_x"}; // col major
            nc_file.def_dim(name + "_x", box_size[0]);
            nc_file.def_dim(name + "_y", box_size[1]);
            if (!var_2d) nc_file.def_dim(name + "_z", box_size[2]);
        } else {
            dim_names = {"z", "y", "x"}; // col major
        }

        if (var_2d) {
            dim_names.erase(dim_names.begin());
        }

        // add a dimension for MF components
        if (mf->nComp() > 1)
        {
            nc_file.def_dim(name + "_nComp", mf->nComp());
            dim_names.emplace(dim_names.begin(), name + "_nComp");
        }

        if (time > -1.0)
        {
            dim_names.emplace(dim_names.begin(), "time");
        }

        nc_file.def_var(name, ncutils::NCDType::Real, dim_names);
        nc_file.exit_def_mode();

        if (write_ghost)
        {
            nc_file.var(name).put_attr("growVect", std::vector<int>(ngrow.begin(), ngrow.end()));
        }
    }

    // Find the current time index in the file
    // This is needed to avoid the time dimension from growing each time a MF is written to file
    // TODO: fix this - move time dim to NCInterface
    int time_index = -1;
    if (time > -1.0)
    {
        size_t time_dim_len = nc_file.dim("time").len();
        size_t time_var_len = nc_file.var("time").shape()[0];

        Real last_time = -1.0;
        if (time_dim_len > 0)
        {
            std::vector<size_t> start = {time_dim_len - 1};
            std::vector<size_t> count = {1};
            nc_file.var("time").get(&last_time, start, count);
        }

        AMREX_ALWAYS_ASSERT(time_dim_len == time_var_len);

        if (last_time != time)
        {
            nc_file.var("time").par_access(NC_COLLECTIVE);
            nc_file.var("time").put(&time, {static_cast<size_t>(std::max(0, static_cast<int>(time_var_len)))}, {1});

            time_dim_len = nc_file.dim("time").len();
            time_var_len = nc_file.var("time").shape()[0];
        }
        time_index = std::max(0, static_cast<int>(time_dim_len - 1));
    }

    const MultiFab *data = mf;
    std::unique_ptr<MultiFab> mf_tmp;
    if (!write_ghost && ngrow != 0)
    {
        mf_tmp = std::make_unique<MultiFab>(mf->boxArray(), mf->DistributionMap(), mf->nComp(), 0, MFInfo(), mf->Factory());
        MultiFab::Copy(*mf_tmp, *mf, 0, 0, mf->nComp(), 0);
        data = mf_tmp.get();
    }

    // TODO: fix - this is required since file opened each time step, so the variable access prop is reset
    nc_file.var(name).par_access(NC_COLLECTIVE);

    amrex::Arena* Arena_Used = amrex::The_Arena();
#ifdef AMREX_USE_GPU
    Arena_Used = amrex::The_Pinned_Arena();
#endif

    for(MFIter mfi(*data); mfi.isValid(); ++mfi)
    {
        const auto &box = mfi.fabbox();
        FArrayBox tmp(box,  data->nComp(), Arena_Used);
        for(int comp = 0; comp < data->nComp(); ++comp)
        {
            //const auto *dataPtr = data->get(mfi).dataPtr(comp);
            const auto *dataPtr = tmp.dataPtr(comp);
            AMREX_ALWAYS_ASSERT(dataPtr != nullptr);

            // TODO: this assumes z always start at 0.. handle box.smallEnd()[2] better for 2D/3D MFs
            std::vector<size_t> starts = {0,
                                          static_cast<unsigned long>(box.smallEnd()[1]),
                                          static_cast<unsigned long>(box.smallEnd()[0])};
            std::vector<size_t> counts = {static_cast<unsigned long>(box.length()[2]),
                                          static_cast<unsigned long>(box.length()[1]),
                                          static_cast<unsigned long>(box.length()[0])};

            if (var_2d) {
                starts.erase(starts.begin());
                counts.erase(counts.begin());
            }

            if (data->nComp() > 1)
            {
                starts.emplace(starts.begin(), comp);
                counts.emplace(counts.begin(), 1);
            }

            if (time > -1.0)
            {
                starts.emplace(starts.begin(), time_index);
                counts.emplace(counts.begin(), 1);
            }

            auto mf_arr = data->array(mfi, comp);
            auto tmp_arr = tmp.array(comp);
            if (var_2d)
            {
                ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    tmp_arr(i, j, k) = mf_arr(i, j, khi);
                });
            } else {
                ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    tmp_arr(i, j, k) = mf_arr(i, j, k);
                });
            }

            nc_file.var(name).put(dataPtr, starts, counts);
        }
    }
}

void
SLM::writeMFVecToNC(ncutils::NCFile &nc_file,
                    const Vector<const MultiFab*> &mf_vec,
                    const Vector<std::string> &mf_names,
                    const Geometry& geom) const
{
    // Helper function to write a vector of MF to a NetCDF file
    // it is assumed that all MFs in mf_vec share geom
}
#endif


// utility to skip to next line in Header
//  -- taken from ERF_Checkpoint
void
SLM::GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}

void SLM::WriteCheckpoint(const int &lev, const std::string &checkpointname) const
{
    auto check_start = amrex::second();

    // write SLM header
    if (ParallelDescriptor::IOProcessor()) {

        amrex::Print() << " Writing SLM checkpoint at level " << lev << std::endl;

        std::string HeaderFileName(checkpointname + "/Level_" + std::to_string(lev) + "/SLM_Header");
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                std::ofstream::trunc |
                                                std::ofstream::binary);
        if(! HeaderFile.good()) {
            FileOpenFailed(HeaderFileName);
        }

        HeaderFile.precision(17);

        // write out title line
        HeaderFile << "Checkpoint file for SLM\n";

        // write out number of soil layers
        HeaderFile << m_nz_lsm << "\n";

        // write out array of dz
        for (int i = 0; i < static_cast<int>(m_dz_lsm.size()); ++i) {
            HeaderFile << m_dz_lsm[i] << " ";
        }
        HeaderFile << "\n";
        HeaderFile << "\n";

        // write out timestep and cal day
        HeaderFile << m_dt << "\n";
        HeaderFile << time << "\n";
        HeaderFile << start_time << "\n";
        HeaderFile << m_calday << "\n";
        HeaderFile << m_orbital_year << "\n";
        HeaderFile << m_orbital_mon << "\n";
        HeaderFile << m_orbital_day << "\n";
        HeaderFile << m_orbital_sec << "\n";
        HeaderFile << "\n";

        // write out constants
        //  - tabs_s
        //  - t00
        //  - LAI0
        //  - tausoil
        //  - zref
        //  - z0_soil
        //  - mws_mx0
        HeaderFile << tabs_s << "\n";
        HeaderFile << t00 << "\n";
        HeaderFile << LAI0 << "\n";
        HeaderFile << tausoil << "\n";
        HeaderFile << zref << "\n";
        HeaderFile << z0_soil << "\n";
        HeaderFile << mws_mx0 << "\n";
        HeaderFile << Rc_max << "\n";
        HeaderFile << T_opt << "\n";
        HeaderFile << "\n";

        // write out flags
        //  - dosoiltnudge
        //  - dosoilwnudge
        //  - reserved compatibility slot for the removed set_from_file flag
        //  - use_param_file
        //  - interpolate_lai
        //  - use_wrf_lai
        //  - use_wrfinput
        HeaderFile << dosoiltnudging << "\n";
        HeaderFile << dosoilwnudging << "\n";
        HeaderFile << 0 << "\n";
        HeaderFile << use_param_file << "\n";
        HeaderFile << interpolate_lai << "\n";
        HeaderFile << use_wrf_lai << "\n";
        HeaderFile << use_wrfinput << "\n";
        HeaderFile << "\n";

        // Write box arrays for SLM
        ba_lsm_2d.writeOn(HeaderFile);
        HeaderFile << '\n';
    }

    amrex::ParallelDescriptor::Barrier();

    const std::string prefix = "SLM_";

    // Use distribution map from data vars, since those are already created and loaded by ERF
    const DistributionMapping& dm = lsm_fab_vars[0]->DistributionMap();
    IntVect ng = IntVect(0,0,0);
    {
        MultiFab imf_lm = amrex::ToMultiFab(landmask);
        VisMF::Write(imf_lm, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "landmask"));

        MultiFab imf_lt = amrex::ToMultiFab(landtype);
        VisMF::Write(imf_lt, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "landtype"));

        MultiFab imf_vt = amrex::ToMultiFab(vegetype);
        VisMF::Write(imf_vt, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "vegetype"));
    }

    {
        // DST
        MultiFab mf(soilt_vars.boxArray(),dm,SLM_DST::NumVars,0);
        MultiFab::Copy(mf,soilt_vars,0,0,SLM_DST::NumVars,ng);
        VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "DST"));
    }

    {
        // DSW
        MultiFab mf(soilw_vars.boxArray(),dm,SLM_DSW::NumVars,0);
        MultiFab::Copy(mf,soilw_vars,0,0,SLM_DSW::NumVars,ng);
        VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "DSW"));
    }

    {
        // net_rad
        MultiFab mf(net_rad.boxArray(),dm,SLM_NetRad::NumVars,0);
        MultiFab::Copy(mf,net_rad,0,0,SLM_NetRad::NumVars,ng);
        VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "net_rad"));
    }

    {
        // slm_diag
        MultiFab mf(slm_diag.boxArray(),dm,SLM_Diag::NumVars,0);
        MultiFab::Copy(mf,slm_diag,0,0,SLM_Diag::NumVars,ng);
        VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "slm_diag"));
    }

    MultiFab mf(ba_lsm_2d,dm,1,ng);

    MultiFab::Copy(mf,LAI,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "LAI"));

    MultiFab::Copy(mf,SAI,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "SAI"));

    MultiFab::Copy(mf,sstxy,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "sstxy"));

    MultiFab::Copy(mf,t_canop,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "t_canop"));

    MultiFab::Copy(mf,mw,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "mw"));

    MultiFab::Copy(mf,mws,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "mws"));

    MultiFab::Copy(mf,t_skin,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "t_skin"));

    MultiFab::Copy(mf,t_ground_skin,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "t_ground_skin"));

    MultiFab::Copy(mf,t_cas,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "t_cas"));

    MultiFab::Copy(mf,q_cas,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "q_cas"));

    MultiFab::Copy(mf,vege_YES,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "vege_YES"));

    MultiFab::Copy(mf,cp_vege,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "cp_vege"));

    MultiFab::Copy(mf,z0_sfc,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "z0_sfc"));

    MultiFab::Copy(mf,Khai_L,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "Khai_L"));

    MultiFab::Copy(mf,phi_1,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "phi_1"));

    MultiFab::Copy(mf,phi_2,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "phi_2"));

    MultiFab::Copy(mf,IR_emis_vege,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "IR_emis_vege"));

    MultiFab::Copy(mf,IR_emis_soil,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "IR_emis_soil"));

    MultiFab::Copy(mf,ztop,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "ztop"));

    MultiFab::Copy(mf,disp_hgt,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "disp_hgt"));

    MultiFab::Copy(mf,Rgl,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "Rgl"));

    MultiFab::Copy(mf,Rc_min,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "Rc_min"));

    MultiFab::Copy(mf,hs_rc,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "hs_rc"));

    MultiFab::Copy(mf,precip_extinc,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "precip_extinc"));

    MultiFab::Copy(mf,mw_mx,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "mw_mx"));

    MultiFab::Copy(mf,mws_mx,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "mws_mx"));

    MultiFab::Copy(mf,BAI,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "BAI"));

    MultiFab::Copy(mf,mw_inc,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "mw_inc"));

    MultiFab::Copy(mf,evapo_dry,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "evapo_dry"));

    MultiFab::Copy(mf,shf_canop,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "shf_canop"));

    MultiFab::Copy(mf,shf_soil,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "shf_soil"));

    MultiFab::Copy(mf,shf_air,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "shf_air"));

    MultiFab::Copy(mf,lhf_canop,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "lhf_canop"));

    MultiFab::Copy(mf,lhf_soil,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "lhf_soil"));

    MultiFab::Copy(mf,lhf_air,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "lhf_air"));

    MultiFab::Copy(mf,albedovis_v,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "albedovis_v"));

    MultiFab::Copy(mf,albedonir_v,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "albedonir_v"));

    MultiFab::Copy(mf,albedovis_s,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "albedovis_s"));

    MultiFab::Copy(mf,albedonir_s,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "albedonir_s"));

    MultiFab::Copy(mf,r_a,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "r_a"));

    MultiFab::Copy(mf,r_b,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "r_b"));

    MultiFab::Copy(mf,r_c,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "r_c"));

    MultiFab::Copy(mf,r_d,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "r_d"));

    MultiFab::Copy(mf,r_soil,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "r_soil"));

    MultiFab::Copy(mf,wet_canop,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "wet_canop"));

    MultiFab::Copy(mf,zrefxy,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "zrefxy"));

    MultiFab::Copy(mf,albold_noahmp,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "albold_noahmp"));

    MultiFab::Copy(mf,tauss_noahmp,0,0,1,ng);
    VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "tauss_noahmp"));

    for (int i = 0; i < static_cast<int>(unmapped_fields.size()); i++) {
        MultiFab mf(lsm_fab_vars[unmapped_fields[i]]->boxArray(),dm,1,IntVect(1,1,1));
        MultiFab::Copy(mf,*(lsm_fab_vars[unmapped_fields[i]]),0,0,1,IntVect(1,1,1));
        VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "Data" + std::to_string(unmapped_fields[i])));
    }

    auto check_end = amrex::second() - check_start;
    ParallelDescriptor::ReduceRealMax(check_end,ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << "    SLM Checkpoint write time = " << check_end << " seconds." << '\n';
}

void SLM::ReadCheckpoint(const int &lev, const std::string &checkpointname)
{
    auto check_start = amrex::second();

    amrex::Print() << " Reading SLM checkpoint at level " << lev << std::endl;

    // Header
    std::string File(checkpointname + "/Level_" + std::to_string(lev) + "/SLM_Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // read in title line
    std::getline(is, line);

    // read in number of soil layers
    is >> m_nz_lsm;

    // read in array of dz
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            m_dz_lsm[i++] = std::stod(word);
        }
    }
    GotoNextLine(is);

    // Read in time
    is >> m_dt;
    is >> time;
    is >> start_time;
    is >> m_calday;
    is >> m_orbital_year;
    is >> m_orbital_mon;
    is >> m_orbital_day;
    is >> m_orbital_sec;
    GotoNextLine(is);

    // Read in constants
    is >> tabs_s;
    is >> t00;
    is >> LAI0;
    is >> tausoil;
    is >> zref;
    is >> z0_soil;
    is >> mws_mx0;
    is >> Rc_max;
    is >> T_opt;
    GotoNextLine(is);

    // Read in flags
    int reserved_slm_flag = 0;
    is >> dosoiltnudging;
    is >> dosoilwnudging;
    is >> reserved_slm_flag;
    is >> use_param_file;
    is >> interpolate_lai;
    is >> use_wrf_lai;
    is >> use_wrfinput;
    amrex::ignore_unused(reserved_slm_flag);
    GotoNextLine(is);

    // read in level 'lev' BoxArray from Header
    BoxArray ba;
    ba.readFrom(is);
    GotoNextLine(is);

    AMREX_ALWAYS_ASSERT(ba == ba_lsm_2d);

    const std::string prefix = "SLM_";

    // Use distribution map from data vars, since those are already created and loaded by ERF
    const DistributionMapping& dm = lsm_fab_vars[0]->DistributionMap();
    IntVect ng = IntVect(0,0,0);
    {
        MultiFab mf(ba_lsm_2d,dm,1,ng);
        VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "landmask"));
        landmask = amrex::cast<iMultiFab>(mf);

        VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "landtype"));
        landtype = amrex::cast<iMultiFab>(mf);

        VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "vegetype"));
        vegetype = amrex::cast<iMultiFab>(mf);
    }

    {
        // DST
        MultiFab mf(soilt_vars.boxArray(),dm,SLM_DST::NumVars,0);
        VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "DST"));
        MultiFab::Copy(soilt_vars,mf,0,0,SLM_DST::NumVars,ng);
    }

    {
        // DSW
        MultiFab mf(soilw_vars.boxArray(),dm,SLM_DSW::NumVars,0);
        VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "DSW"));
        MultiFab::Copy(soilw_vars,mf,0,0,SLM_DSW::NumVars,ng);
    }

    {
        // net_rad
        MultiFab mf(net_rad.boxArray(),dm,SLM_NetRad::NumVars,0);
        VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "net_rad"));
        MultiFab::Copy(net_rad,mf,0,0,SLM_NetRad::NumVars,ng);
    }

    {
        // slm_diag
        MultiFab mf(slm_diag.boxArray(),dm,SLM_Diag::NumVars,0);
        VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "slm_diag"));
        MultiFab::Copy(slm_diag,mf,0,0,SLM_Diag::NumVars,ng);
    }

    MultiFab mf(ba_lsm_2d,dm,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "LAI"));
    MultiFab::Copy(LAI,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "SAI"));
    MultiFab::Copy(SAI,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "sstxy"));
    MultiFab::Copy(sstxy,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "t_canop"));
    MultiFab::Copy(t_canop,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "mw"));
    MultiFab::Copy(mw,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "mws"));
    MultiFab::Copy(mws,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "t_skin"));
    MultiFab::Copy(t_skin,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "t_ground_skin"));
    MultiFab::Copy(t_ground_skin,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "t_cas"));
    MultiFab::Copy(t_cas,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "q_cas"));
    MultiFab::Copy(q_cas,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "vege_YES"));
    MultiFab::Copy(vege_YES,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "cp_vege"));
    MultiFab::Copy(cp_vege,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "z0_sfc"));
    MultiFab::Copy(z0_sfc,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "Khai_L"));
    MultiFab::Copy(Khai_L,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "phi_1"));
    MultiFab::Copy(phi_1,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "phi_2"));
    MultiFab::Copy(phi_2,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "IR_emis_vege"));
    MultiFab::Copy(IR_emis_vege,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "IR_emis_soil"));
    MultiFab::Copy(IR_emis_soil,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "ztop"));
    MultiFab::Copy(ztop,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "disp_hgt"));
    MultiFab::Copy(disp_hgt,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "Rgl"));
    MultiFab::Copy(Rgl,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "Rc_min"));
    MultiFab::Copy(Rc_min,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "hs_rc"));
    MultiFab::Copy(hs_rc,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "precip_extinc"));
    MultiFab::Copy(precip_extinc,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "mw_mx"));
    MultiFab::Copy(mw_mx,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "mws_mx"));
    MultiFab::Copy(mws_mx,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "BAI"));
    MultiFab::Copy(BAI,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "mw_inc"));
    MultiFab::Copy(mw_inc,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "evapo_dry"));
    MultiFab::Copy(evapo_dry,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "shf_canop"));
    MultiFab::Copy(shf_canop,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "shf_soil"));
    MultiFab::Copy(shf_soil,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "shf_air"));
    MultiFab::Copy(shf_air,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "lhf_canop"));
    MultiFab::Copy(lhf_canop,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "lhf_soil"));
    MultiFab::Copy(lhf_soil,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "lhf_air"));
    MultiFab::Copy(lhf_air,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "albedovis_v"));
    MultiFab::Copy(albedovis_v,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "albedonir_v"));
    MultiFab::Copy(albedonir_v,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "albedovis_s"));
    MultiFab::Copy(albedovis_s,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "albedonir_s"));
    MultiFab::Copy(albedonir_s,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "r_a"));
    MultiFab::Copy(r_a,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "r_b"));
    MultiFab::Copy(r_b,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "r_c"));
    MultiFab::Copy(r_c,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "r_d"));
    MultiFab::Copy(r_d,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "r_soil"));
    MultiFab::Copy(r_soil,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "wet_canop"));
    MultiFab::Copy(wet_canop,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "zrefxy"));
    MultiFab::Copy(zrefxy,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "albold_noahmp"));
    MultiFab::Copy(albold_noahmp,mf,0,0,1,ng);

    VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "tauss_noahmp"));
    MultiFab::Copy(tauss_noahmp,mf,0,0,1,ng);

    for (int i = 0; i < static_cast<int>(unmapped_fields.size()); i++) {
        MultiFab mf(lsm_fab_vars[unmapped_fields[i]]->boxArray(),dm,1,IntVect(1,1,1));
        VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "Data" + std::to_string(unmapped_fields[i])));
        MultiFab::Copy(*(lsm_fab_vars[unmapped_fields[i]]),mf,0,0,1,IntVect(1,1,1));
    }

    first_step = false;
    wrfinput_initialized = use_wrfinput;

    // Initialize common parameters
    init_layer_depths();
    init_soil_vars();
    rebuild_restart_fields();

    auto check_end = amrex::second() - check_start;
    ParallelDescriptor::ReduceRealMax(check_end,ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << "    SLM Checkpoint load time = " << check_end << " seconds." << '\n';
}
