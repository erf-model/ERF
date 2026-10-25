#include "ERF_SLM.H"

#include <AMReX_PlotFileUtil.H>
#include "ERF_NCInterface.H"
#include "ERF_ReadCheckpointMF.H"

using namespace amrex;

#ifdef ERF_USE_NETCDF
void
SLM::writeSLM_NetCDF (const MultiFab& mf, const Vector<std::string>& varnames, const amrex::Real cur_time,
                      const std::string plot_prefix, const int /*level_step*/)
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
        //writeMFtoNC(ncf, lsm_fab_vars[var].get(), LsmVarName_Full[var], cur_time, true); // write ghost cells
        writeMFtoNC(ncf, lsm_fab_vars[var].get(), LsmVarName_Full[var], cur_time, false); // write ghost cells
    }

    // Write out SLM 2D fields
    // TODO: fix this
    //  2D fields are combined into the single mf with varnames components
    AMREX_ALWAYS_ASSERT(mf.nComp() == varnames.size());
    for (int var = 0; var < mf.nComp(); ++var)
    {
        MultiFab fab(mf, make_alias, var, 1);
        writeMFtoNC(ncf, &fab, varnames[var], cur_time);
    }

    writeMFtoNC(ncf, &cp_vege, "cp_vege", cur_time);
    writeMFtoNC(ncf, &z0_sfc, "z0_sfc", cur_time);

    //writeMFtoNC(ncf, &Khai_L, "Khai_L", -1.0);

    //writeMFtoNC(ncf, &ustar, "ustar", cur_time); // ustar and tstar already saved in 2D fab array above

    // additional slm outputs
    for (int var = 0; var < slm_diag.nComp(); ++var)
    {
        MultiFab fab(slm_diag, make_alias, var, 1);
        writeMFtoNC(ncf, &fab, diag_names[var], cur_time);
    }

    // Save soil temperature and moisture diagnostic variables
    for (int var = 0; var < soilt_vars.nComp(); ++var)
    {
        MultiFab fab (soilt_vars, make_alias, var, 1);
        writeMFtoNC(ncf, &fab, soilt_var_names[var], cur_time);
    }

    for (int var = 0; var < soilw_vars.nComp() - 1; ++var)
    {
        MultiFab fab (soilw_vars, make_alias, var, 1);
        writeMFtoNC(ncf, &fab, soilw_var_names[var], cur_time);
    }

    ncf.close();
}

void
SLM::writeNCHeader (ncutils::NCFile &nc_file, const amrex::Geometry &geom)
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

    // The x, y and z coordinates are written below from inside an MFIter loop, so
    // each rank makes one put call per box that it owns.  Ranks can own different
    // numbers of boxes, so we must use independent rather than collective access,
    // otherwise the write hangs waiting for ranks that have no box left.
    nc_file.var("time").par_access(NC_COLLECTIVE);
    nc_file.var("x").par_access(NC_INDEPENDENT);
    nc_file.var("y").par_access(NC_INDEPENDENT);
    nc_file.var("z").par_access(NC_INDEPENDENT);

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

void SLM::writeMFtoNC (ncutils::NCFile &nc_file, const MultiFab* mf,
                       const std::string name, Real cur_time, bool write_ghost)
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
        if (cur_time > -1.0 && !nc_file.has_dim("time"))
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
            int x_size = static_cast<int>(nc_file.dim("x").len());
            int y_size = static_cast<int>(nc_file.dim("y").len());

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

        if (cur_time > -1.0)
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
    if (cur_time > -1.0)
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

        if (last_time != cur_time)
        {
            nc_file.var("time").par_access(NC_COLLECTIVE);
            nc_file.var("time").put(&cur_time, {static_cast<size_t>(std::max(0, static_cast<int>(time_var_len)))}, {1});

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

    // The data below is written from inside an MFIter loop, so each rank makes one put
    // call per box that it owns.  Ranks can own different numbers of boxes, so we must
    // use independent rather than collective access.  The access mode is reset each time
    // the file is opened, so it is set here on every call.
    nc_file.var(name).par_access(NC_INDEPENDENT);

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

            if (cur_time > -1.0)
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
            Gpu::streamSynchronize();

            nc_file.var(name).put(dataPtr, starts, counts);
        }
    }
}

void
SLM::writeMFVecToNC (ncutils::NCFile& /*nc_file*/,
                     const Vector<const MultiFab*>& /*mf_vec*/,
                     const Vector<std::string>& /*mf_names*/,
                     const Geometry& /*geom*/) const
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

void SLM::WriteCheckpoint (const int &lev, const std::string &checkpointname) const
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

        //
        // Where SLM is in its lifecycle.  With a wrfinput run the surface fields and
        // the canopy-adjusted reference height are not built in Init -- they need the
        // file, and are deferred to the first Advance -- so a checkpoint written
        // before any step has been taken holds zrefxy at the scalar slm.zref
        // placeholder rather than ztop + dz0.  The reader used to declare the model
        // initialized unconditionally, which is right for every later checkpoint and
        // wrong for that one.  Written last so a checkpoint from before this simply
        // ends after the BoxArray and the reader keeps its old default.
        //
        HeaderFile << first_step << " " << wrfinput_initialized << "\n";
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
        MultiFab mf_unmapped(lsm_fab_vars[unmapped_fields[i]]->boxArray(),dm,1,IntVect(1,1,1));
        MultiFab::Copy(mf_unmapped,*(lsm_fab_vars[unmapped_fields[i]]),0,0,1,IntVect(1,1,1));
        VisMF::Write(mf_unmapped, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "Data" + std::to_string(unmapped_fields[i])));
    }

    auto check_end = amrex::second() - check_start;
    ParallelDescriptor::ReduceRealMax(check_end,ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << "    SLM Checkpoint write time = " << check_end << " seconds." << '\n';
}

void SLM::ReadCheckpoint (const int &lev, const std::string &checkpointname)
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

    // read in number of soil layers and the array of dz. Init has already built the
    // soil geometry from slm.nsoil / slm.soil_dz, so these are checked against it
    // rather than adopted.
    int nz_chk = 0;
    is >> nz_chk;
    GotoNextLine(is);

    Vector<Real> dz_chk;
    std::getline(is, line);
    {
        std::istringstream lis(line);
        while (lis >> word) {
            dz_chk.push_back(static_cast<Real>(std::stod(word)));
        }
    }

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        nz_chk == m_nz_lsm && static_cast<int>(dz_chk.size()) == m_nz_lsm,
        "SLM::ReadCheckpoint: the checkpoint's number of soil layers differs from slm.nsoil");
    for (int k = 0; k < m_nz_lsm; ++k) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::abs(dz_chk[k] - m_dz_lsm[k]) <= Real(1.e-6) * m_dz_lsm[k],
            "SLM::ReadCheckpoint: the checkpoint's soil layer thicknesses differ from slm.soil_dz");
    }

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

    //
    // The grids the checkpoint was written on. They need not be this run's: a restart that
    // re-makes the level-0 grids (issue 4225) builds the model on the new ones, and every
    // field below is read on the grids in the file and redistributed onto those.
    //
    amrex::ignore_unused(ba);

    //
    // Lifecycle flags, if this checkpoint carries them (see WriteCheckpoint).  The
    // defaults are what the reader always did: treat the model as initialized.
    //
    // NOTE: the BoxArray read above is already followed by its own GotoNextLine, so
    //       the stream sits at the start of this line -- do not skip another.
    //
    int chk_first_step = 0;
    int chk_wrfinput_initialized = (use_wrfinput) ? 1 : 0;
    {
        int a, b;
        if (is >> a >> b) { chk_first_step = a; chk_wrfinput_initialized = b; }
    }

    const std::string prefix = "SLM_";

    // Every field is planar or a column of soil layers: periodic laterally as the domain is,
    // never in z
    const Periodicity period = erf_checkpoint::lateral_periodicity(m_geom);
    auto read_field = [&] (MultiFab& dst, const std::string& name, int ncomp,
                           const IntVect& ng_read, bool src_ghosts)
    {
        erf_checkpoint::read_checkpoint_multifab(
            dst, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + name),
            0, 0, ncomp, ng_read, src_ghosts, period);
    };

    // Use distribution map from data vars, since those are already created and loaded by ERF
    const DistributionMapping& dm = lsm_fab_vars[0]->DistributionMap();
    erf_checkpoint::read_checkpoint_imultifab(landmask, ba_lsm_2d, dm,
        MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "landmask"), period);
    erf_checkpoint::read_checkpoint_imultifab(landtype, ba_lsm_2d, dm,
        MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "landtype"), period);
    erf_checkpoint::read_checkpoint_imultifab(vegetype, ba_lsm_2d, dm,
        MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "vegetype"), period);

    // DST
    read_field(soilt_vars, "DST", SLM_DST::NumVars, IntVect(0), false);

    // DSW
    read_field(soilw_vars, "DSW", SLM_DSW::NumVars, IntVect(0), false);

    // net_rad
    read_field(net_rad, "net_rad", SLM_NetRad::NumVars, IntVect(0), false);

    // slm_diag
    read_field(slm_diag, "slm_diag", SLM_Diag::NumVars, IntVect(0), false);

    read_field(LAI, "LAI", 1, IntVect(0), false);

    read_field(SAI, "SAI", 1, IntVect(0), false);

    read_field(sstxy, "sstxy", 1, IntVect(0), false);

    read_field(t_canop, "t_canop", 1, IntVect(0), false);

    read_field(mw, "mw", 1, IntVect(0), false);

    read_field(mws, "mws", 1, IntVect(0), false);

    read_field(t_skin, "t_skin", 1, IntVect(0), false);

    read_field(t_ground_skin, "t_ground_skin", 1, IntVect(0), false);

    read_field(t_cas, "t_cas", 1, IntVect(0), false);

    read_field(q_cas, "q_cas", 1, IntVect(0), false);

    read_field(vege_YES, "vege_YES", 1, IntVect(0), false);

    read_field(cp_vege, "cp_vege", 1, IntVect(0), false);

    read_field(z0_sfc, "z0_sfc", 1, IntVect(0), false);

    read_field(Khai_L, "Khai_L", 1, IntVect(0), false);

    read_field(phi_1, "phi_1", 1, IntVect(0), false);

    read_field(phi_2, "phi_2", 1, IntVect(0), false);

    read_field(IR_emis_vege, "IR_emis_vege", 1, IntVect(0), false);

    read_field(IR_emis_soil, "IR_emis_soil", 1, IntVect(0), false);

    read_field(ztop, "ztop", 1, IntVect(0), false);

    read_field(disp_hgt, "disp_hgt", 1, IntVect(0), false);

    read_field(Rgl, "Rgl", 1, IntVect(0), false);

    read_field(Rc_min, "Rc_min", 1, IntVect(0), false);

    read_field(hs_rc, "hs_rc", 1, IntVect(0), false);

    read_field(precip_extinc, "precip_extinc", 1, IntVect(0), false);

    read_field(mw_mx, "mw_mx", 1, IntVect(0), false);

    read_field(mws_mx, "mws_mx", 1, IntVect(0), false);

    read_field(BAI, "BAI", 1, IntVect(0), false);

    read_field(mw_inc, "mw_inc", 1, IntVect(0), false);

    read_field(evapo_dry, "evapo_dry", 1, IntVect(0), false);

    read_field(shf_canop, "shf_canop", 1, IntVect(0), false);

    read_field(shf_soil, "shf_soil", 1, IntVect(0), false);

    read_field(shf_air, "shf_air", 1, IntVect(0), false);

    read_field(lhf_canop, "lhf_canop", 1, IntVect(0), false);

    read_field(lhf_soil, "lhf_soil", 1, IntVect(0), false);

    read_field(lhf_air, "lhf_air", 1, IntVect(0), false);

    read_field(albedovis_v, "albedovis_v", 1, IntVect(0), false);

    read_field(albedonir_v, "albedonir_v", 1, IntVect(0), false);

    read_field(albedovis_s, "albedovis_s", 1, IntVect(0), false);

    read_field(albedonir_s, "albedonir_s", 1, IntVect(0), false);

    read_field(r_a, "r_a", 1, IntVect(0), false);

    read_field(r_b, "r_b", 1, IntVect(0), false);

    read_field(r_c, "r_c", 1, IntVect(0), false);

    read_field(r_d, "r_d", 1, IntVect(0), false);

    read_field(r_soil, "r_soil", 1, IntVect(0), false);

    read_field(wet_canop, "wet_canop", 1, IntVect(0), false);

    read_field(zrefxy, "zrefxy", 1, IntVect(0), false);

    read_field(albold_noahmp, "albold_noahmp", 1, IntVect(0), false);

    read_field(tauss_noahmp, "tauss_noahmp", 1, IntVect(0), false);

    for (int i = 0; i < static_cast<int>(unmapped_fields.size()); i++) {
        read_field(*(lsm_fab_vars[unmapped_fields[i]]), "Data" + std::to_string(unmapped_fields[i]),
                   1, IntVect(1,1,1), true);
    }

    //
    // Restore where the model was, rather than asserting it was past its first step.
    // A checkpoint taken before any step has been advanced leaves first_step set, so
    // the first Advance after the restart does the deferred wrfinput read and rebuilds
    // zrefxy -- exactly what the run being continued did at that point.
    //
    first_step           = (chk_first_step != 0);
    wrfinput_initialized = (chk_wrfinput_initialized != 0);

    // A checkpoint is written after Update_Lsm_Vars, so tref, qref and pref hold the state
    has_atmos_state         = true;
    init_from_atmos_pending = false;

    // Initialize common parameters
    init_layer_depths();
    init_soil_vars();
    rebuild_restart_fields();

    auto check_end = amrex::second() - check_start;
    ParallelDescriptor::ReduceRealMax(check_end,ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << "    SLM Checkpoint load time = " << check_end << " seconds." << '\n';
}
