
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "ERF.H"
#include "AMReX_PlotFileUtil.H"

using namespace amrex;

/**
 * Utility to skip to next line in Header file input stream.
 */
void
ERF::GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}

/**
 * ERF function for writing a checkpoint file.
 */
void
ERF::WriteCheckpointFile () const
{
    // chk00010            write a checkpoint file with this root directory
    // chk00010/Header     this contains information you need to save (e.g., finest_level, t_new, etc.) and also
    //                     the BoxArrays at each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                these subdirectories will hold the MultiFab data at each level of refinement

    // checkpoint file name, e.g., chk00010
    const std::string& checkpointname = Concatenate(check_file,istep[0],5);

    Print() << "Writing native checkpoint " << checkpointname << "\n";

    const int nlevels = finest_level+1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then build
    // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    int ncomp_cons = vars_new[0][Vars::cons].nComp();

    // write Header file
    if (ParallelDescriptor::IOProcessor()) {

       std::string HeaderFileName(checkpointname + "/Header");
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
       HeaderFile << "Checkpoint file for ERF\n";

       // write out finest_level
       HeaderFile << finest_level << "\n";

       // write the number of components
       // for each variable we store

       // conservative, cell-centered vars
       HeaderFile << ncomp_cons << "\n";

       // x-velocity on faces
       HeaderFile << 1 << "\n";

       // y-velocity on faces
       HeaderFile << 1 << "\n";

       // z-velocity on faces
       HeaderFile << 1 << "\n";

       // write out array of istep
       for (int i = 0; i < istep.size(); ++i) {
           HeaderFile << istep[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of dt
       for (int i = 0; i < dt.size(); ++i) {
           HeaderFile << dt[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of t_new
       for (int i = 0; i < t_new.size(); ++i) {
           HeaderFile << t_new[i] << " ";
       }
       HeaderFile << "\n";

       // write the BoxArray at each level
       for (int lev = 0; lev <= finest_level; ++lev) {
           boxArray(lev).writeOn(HeaderFile);
           HeaderFile << '\n';
       }

       // Write separate file that tells how many components we have of the base state
       std::string BaseStateFileName(checkpointname + "/num_base_state_comps");
       std::ofstream BaseStateFile;
       BaseStateFile.open(BaseStateFileName.c_str(), std::ofstream::out   |
                                                     std::ofstream::trunc |
                                                     std::ofstream::binary);
       if(! BaseStateFile.good()) {
           FileOpenFailed(BaseStateFileName);
       } else {
           // write out number of components in base state
           BaseStateFile << BaseState::num_comps      << "\n";
           BaseStateFile << base_state[0].nGrowVect() << "\n";
       }
   }

    // write the MultiFab data to, e.g., chk00010/Level_0/
    // Here we make copies of the MultiFab with no ghost cells
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        MultiFab cons(grids[lev],dmap[lev],ncomp_cons,0);
        MultiFab::Copy(cons,vars_new[lev][Vars::cons],0,0,ncomp_cons,0);
        VisMF::Write(cons, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Cell"));

        MultiFab xvel(convert(grids[lev],IntVect(1,0,0)),dmap[lev],1,0);
        MultiFab::Copy(xvel,vars_new[lev][Vars::xvel],0,0,1,0);
        VisMF::Write(xvel, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "XFace"));

        MultiFab yvel(convert(grids[lev],IntVect(0,1,0)),dmap[lev],1,0);
        MultiFab::Copy(yvel,vars_new[lev][Vars::yvel],0,0,1,0);
        VisMF::Write(yvel, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "YFace"));

        MultiFab zvel(convert(grids[lev],IntVect(0,0,1)),dmap[lev],1,0);
        MultiFab::Copy(zvel,vars_new[lev][Vars::zvel],0,0,1,0);
        VisMF::Write(zvel, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "ZFace"));

        // Note that we write the ghost cells of the base state (unlike above)
        IntVect ng_base = base_state[lev].nGrowVect();
        int  ncomp_base = base_state[lev].nComp();
        MultiFab base(grids[lev],dmap[lev],ncomp_base,ng_base);
        MultiFab::Copy(base,base_state[lev],0,0,ncomp_base,ng_base);
        VisMF::Write(base, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "BaseState"));

        if (SolverChoice::mesh_type != MeshType::ConstantDz)  {
            // Note that we also write the ghost cells of z_phys_nd
            IntVect ng = z_phys_nd[lev]->nGrowVect();
            MultiFab z_height(convert(grids[lev],IntVect(1,1,1)),dmap[lev],1,ng);
            MultiFab::Copy(z_height,*z_phys_nd[lev],0,0,1,ng);
            VisMF::Write(z_height, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Z_Phys_nd"));
        }

        // We must read and write qmoist with ghost cells because we don't directly impose BCs on these vars
        // Write the moisture model restart variables
        std::vector<int> qmoist_indices;
        std::vector<std::string> qmoist_names;
        micro->Get_Qmoist_Restart_Vars(lev, solverChoice, qmoist_indices, qmoist_names);
        int qmoist_nvar = qmoist_indices.size();
        for (int var = 0; var < qmoist_nvar; var++) {
            const int ncomp  = 1;
            IntVect ng_moist = qmoist[lev][qmoist_indices[var]]->nGrowVect();
            MultiFab moist_vars(grids[lev],dmap[lev],ncomp,ng_moist);
            MultiFab::Copy(moist_vars,*(qmoist[lev][qmoist_indices[var]]),0,0,ncomp,ng_moist);
            VisMF::Write(moist_vars, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", qmoist_names[var]));
        }

#if defined(ERF_USE_WINDFARM)
        if(solverChoice.windfarm_type == WindFarmType::Fitch or
           solverChoice.windfarm_type == WindFarmType::EWP or
           solverChoice.windfarm_type == WindFarmType::SimpleAD){
            IntVect ng_turb = Nturb[lev].nGrowVect();
            MultiFab mf_Nturb(grids[lev],dmap[lev],1,ng_turb);
            MultiFab::Copy(mf_Nturb,Nturb[lev],0,0,1,ng_turb);
            VisMF::Write(mf_Nturb, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "NumTurb"));
        }
#endif

        if (solverChoice.lsm_type != LandSurfaceType::None) {
            for (int mvar(0); mvar<lsm_data[lev].size(); ++mvar) {
                BoxArray ba = lsm_data[lev][mvar]->boxArray();
                DistributionMapping dm = lsm_data[lev][mvar]->DistributionMap();
                IntVect ng = lsm_data[lev][mvar]->nGrowVect();
                int nvar = lsm_data[lev][mvar]->nComp();
                MultiFab lsm_vars(ba,dm,nvar,ng);
                MultiFab::Copy(lsm_vars,*(lsm_data[lev][mvar]),0,0,nvar,ng);
                VisMF::Write(lsm_vars, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "LsmVars"));
            }
        }

        IntVect ng = mapfac[lev][MapFacType::mx]->nGrowVect();
        MultiFab mf_m(ba2d[lev],dmap[lev],1,ng);
        MultiFab::Copy(mf_m,*mapfac[lev][MapFacType::mx],0,0,1,ng);
        VisMF::Write(mf_m, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "MapFactor_mx"));

        if (MapFacType::mx != MapFacType::my) {
            MultiFab::Copy(mf_m,*mapfac[lev][MapFacType::my],0,0,1,ng);
            VisMF::Write(mf_m, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "MapFactor_my"));
        }

        ng = mapfac[lev][MapFacType::ux]->nGrowVect();
        MultiFab mf_u(convert(ba2d[lev],IntVect(1,0,0)),dmap[lev],1,ng);
        MultiFab::Copy(mf_u,*mapfac[lev][MapFacType::ux],0,0,1,ng);
        VisMF::Write(mf_u, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "MapFactor_ux"));

        if (MapFacType::mx != MapFacType::my) {
            MultiFab::Copy(mf_u,*mapfac[lev][MapFacType::uy],0,0,1,ng);
            VisMF::Write(mf_u, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "MapFactor_uy"));
        }

        ng = mapfac[lev][MapFacType::vx]->nGrowVect();
        MultiFab mf_v(convert(ba2d[lev],IntVect(0,1,0)),dmap[lev],1,ng);
        MultiFab::Copy(mf_v,*mapfac[lev][MapFacType::vx],0,0,1,ng);
        VisMF::Write(mf_v, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "MapFactor_vx"));

        if (MapFacType::mx != MapFacType::my) {
            MultiFab::Copy(mf_v,*mapfac[lev][MapFacType::vy],0,0,1,ng);
            VisMF::Write(mf_v, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "MapFactor_vy"));
        }

        if (m_SurfaceLayer)  {
            amrex::Print() << "Writing SurfaceLayer variables" << std::endl;
            ng = vars_new[lev][Vars::cons].nGrowVect(); ng[2]=0;
            MultiFab   m_var(ba2d[lev],dmap[lev],1,ng);
            MultiFab* src = nullptr;

            // U*
            src = m_SurfaceLayer->get_u_star(lev);
            MultiFab::Copy(m_var,*src,0,0,1,ng);
            VisMF::Write(m_var, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Ustar"));

            // W*
            src = m_SurfaceLayer->get_w_star(lev);
            MultiFab::Copy(m_var,*src,0,0,1,ng);
            VisMF::Write(m_var, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Wstar"));

            // T*
            src = m_SurfaceLayer->get_t_star(lev);
            MultiFab::Copy(m_var,*src,0,0,1,ng);
            VisMF::Write(m_var, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Tstar"));

            // Q*
            src = m_SurfaceLayer->get_q_star(lev);
            MultiFab::Copy(m_var,*src,0,0,1,ng);
            VisMF::Write(m_var, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Qstar"));

            // Olen
            src = m_SurfaceLayer->get_olen(lev);
            MultiFab::Copy(m_var,*src,0,0,1,ng);
            VisMF::Write(m_var, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Olen"));

            // PBLH
            src = m_SurfaceLayer->get_pblh(lev);
            MultiFab::Copy(m_var,*src,0,0,1,ng);
            VisMF::Write(m_var, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "PBLH"));

            // Z0
            for (MFIter mfi(m_var); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.growntilebox();
                Array4<const Real> const& fab_arr = m_SurfaceLayer->get_z0(lev)->const_array();
                Array4<      Real> const&  mv_arr = m_var.array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    mv_arr(i,j,k) = fab_arr(i,j,k);
                });
            }
            VisMF::Write(m_var, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Z0"));
        }

        if (sst_lev[lev][0]) {
            int ntimes = 1;
            ng = vars_new[lev][Vars::cons].nGrowVect(); ng[2]=0;
            MultiFab sst_at_t(ba2d[lev],dmap[lev],1,ng);
            for (int nt(0); nt<ntimes; ++nt) {
                MultiFab::Copy(sst_at_t,*sst_lev[lev][nt],0,0,1,ng);
                VisMF::Write(sst_at_t, MultiFabFileFullPrefix(lev, checkpointname, "Level_",
                                                             "SST_" + std::to_string(nt)));
            }
        }

        if (tsk_lev[lev][0]) {
            int ntimes = 1;
            ng = vars_new[lev][Vars::cons].nGrowVect(); ng[2]=0;
            MultiFab tsk_at_t(ba2d[lev],dmap[lev],1,ng);
            for (int nt(0); nt<ntimes; ++nt) {
                MultiFab::Copy(tsk_at_t,*tsk_lev[lev][nt],0,0,1,ng);
                VisMF::Write(tsk_at_t, MultiFabFileFullPrefix(lev, checkpointname, "Level_",
                                                             "TSK_" + std::to_string(nt)));
            }
        }

        {
            int ntimes = 1;
            ng = vars_new[lev][Vars::cons].nGrowVect(); ng[2]=0;
            MultiFab lmask_at_t(ba2d[lev],dmap[lev],1,ng);
            for (int nt(0); nt<ntimes; ++nt) {
                for (MFIter mfi(lmask_at_t); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.growntilebox();
                    Array4<int>  const& src_arr = lmask_lev[lev][nt]->array(mfi);
                    Array4<Real> const& dst_arr = lmask_at_t.array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        dst_arr(i,j,k) = Real(src_arr(i,j,k));
                    });
                }
                VisMF::Write(lmask_at_t, MultiFabFileFullPrefix(lev, checkpointname, "Level_",
                                                              "LMASK_" + std::to_string(nt)));
            }
        }

#ifdef ERF_USE_NETCDF
        IntVect ngv = ng; ngv[2] = 0;

        // Write lat/lon if it exists
        if (lat_m[lev] && lon_m[lev] && solverChoice.has_lat_lon) {
            amrex::Print() << "Writing Lat/Lon variables" << std::endl;
            MultiFab lat(ba2d[lev],dmap[lev],1,ngv);
            MultiFab lon(ba2d[lev],dmap[lev],1,ngv);
            MultiFab::Copy(lat,*lat_m[lev],0,0,1,ngv);
            MultiFab::Copy(lon,*lon_m[lev],0,0,1,ngv);
            VisMF::Write(lat, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "LAT"));
            VisMF::Write(lon, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "LON"));
        }

        // Write sinPhi and cosPhi if it exists
        if (cosPhi_m[lev] && sinPhi_m[lev] && solverChoice.variable_coriolis) {
            amrex::Print() << "Writing Coriolis factors" << std::endl;
            MultiFab sphi(ba2d[lev],dmap[lev],1,ngv);
            MultiFab cphi(ba2d[lev],dmap[lev],1,ngv);
            MultiFab::Copy(sphi,*sinPhi_m[lev],0,0,1,ngv);
            MultiFab::Copy(cphi,*cosPhi_m[lev],0,0,1,ngv);
            VisMF::Write(sphi, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "SinPhi"));
            VisMF::Write(cphi, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "CosPhi"));
        }

        if (solverChoice.use_real_bcs && solverChoice.init_type == InitType::WRFInput) {
            MultiFab tmp1d(ba1d[lev],dmap[lev],1,ngv);
            MultiFab tmp2d(ba2d[lev],dmap[lev],1,ngv);

            MultiFab::Copy(tmp1d,*mf_C1H[lev],0,0,1,ngv);
            VisMF::Write(tmp1d, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "C1H"));

            MultiFab::Copy(tmp1d,*mf_C2H[lev],0,0,1,ngv);
            VisMF::Write(tmp1d, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "C2H"));

            MultiFab::Copy(tmp2d,*mf_MUB[lev],0,0,1,ngv);
            VisMF::Write(tmp2d, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "MUB"));
        }
#endif

    } // for lev

#ifdef ERF_USE_PARTICLES
   particleData.Checkpoint(checkpointname);
#endif

#if 0
#ifdef ERF_USE_NETCDF
   // Write bdy_data files
   if ( ParallelDescriptor::IOProcessor() &&
        ((solverChoice.init_type==InitType::WRFInput) || (solverChoice.init_type==InitType::Metgrid)) &&
         solverChoice.use_real_bcs )
   {
       // Vector dimensions
       int num_time = bdy_data_xlo.size();
       int num_var  = bdy_data_xlo[0].size();

       // Open header file and write to it
       std::ofstream bdy_h_file(MultiFabFileFullPrefix(0, checkpointname, "Level_", "bdy_H"));
       bdy_h_file << std::setprecision(1) << std::fixed;
       bdy_h_file << num_time << "\n";
       bdy_h_file << num_var  << "\n";
       bdy_h_file << start_bdy_time << "\n";
       bdy_h_file << bdy_time_interval << "\n";
       bdy_h_file << real_width << "\n";
       for (int ivar(0); ivar<num_var; ++ivar) {
           bdy_h_file << bdy_data_xlo[0][ivar].box() << "\n";
           bdy_h_file << bdy_data_xhi[0][ivar].box() << "\n";
           bdy_h_file << bdy_data_ylo[0][ivar].box() << "\n";
           bdy_h_file << bdy_data_yhi[0][ivar].box() << "\n";
       }

       // Open data file and write to it
       std::ofstream bdy_d_file(MultiFabFileFullPrefix(0, checkpointname, "Level_", "bdy_D"));
       for (int itime(0); itime<num_time; ++itime) {
           if (bdy_data_xlo[itime].size() > 0) {
               for (int ivar(0); ivar<num_var; ++ivar) {
                   bdy_data_xlo[itime][ivar].writeOn(bdy_d_file,0,1);
                   bdy_data_xhi[itime][ivar].writeOn(bdy_d_file,0,1);
                   bdy_data_ylo[itime][ivar].writeOn(bdy_d_file,0,1);
                   bdy_data_yhi[itime][ivar].writeOn(bdy_d_file,0,1);
               }
           }
       }
   }
#endif
#endif

}

/**
 * ERF function for reading data from a checkpoint file during restart.
 */
void
ERF::ReadCheckpointFile ()
{
    Print() << "Restart from native checkpoint " << restart_chkfile << "\n";

    // Header
    std::string File(restart_chkfile + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    int chk_ncomp_cons, chk_ncomp;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level;
    GotoNextLine(is);

    // read the number of components
    // for each variable we store

    // conservative, cell-centered vars
    is >> chk_ncomp_cons;
    GotoNextLine(is);

    // x-velocity on faces
    is >> chk_ncomp;
    GotoNextLine(is);
    AMREX_ASSERT(chk_ncomp == 1);

    // y-velocity on faces
    is >> chk_ncomp;
    GotoNextLine(is);
    AMREX_ASSERT(chk_ncomp == 1);

    // z-velocity on faces
    is >> chk_ncomp;
    GotoNextLine(is);
    AMREX_ASSERT(chk_ncomp == 1);

    // read in array of istep
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            istep[i++] = std::stoi(word);
        }
    }

    // read in array of dt
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            dt[i++] = std::stod(word);
        }
    }

    // read in array of t_new
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            t_new[i++] = std::stod(word);
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        // read in level 'lev' BoxArray from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // create a distribution mapping
        DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

        MakeNewLevelFromScratch (lev, t_new[lev], ba, dm);
    }

    // ncomp is only valid after we MakeNewLevelFromScratch (asks micro how many vars)
    // NOTE: Data is written over ncomp, so check that we match the header file
    int ncomp_cons = vars_new[0][Vars::cons].nComp();

    // NOTE: QKE was removed so this is for backward compatibility
    AMREX_ASSERT((chk_ncomp_cons==ncomp_cons) || ((chk_ncomp_cons-1)==ncomp_cons));
    //
    // See if we have a written separate file that tells how many components and how many ghost cells
    // we have of the base state
    //
    // If we can't find the file, then set the number of components to the original number = 3
    //
    int ncomp_base_to_read = 3;
    IntVect ng_base = IntVect{1};
    {
         std::string BaseStateFile(restart_chkfile + "/num_base_state_comps");

         if (amrex::FileExists(BaseStateFile))
         {
             Vector<char> BaseStatefileCharPtr;
             ParallelDescriptor::ReadAndBcastFile(BaseStateFile, BaseStatefileCharPtr);
             std::string BaseStatefileCharPtrString(BaseStatefileCharPtr.dataPtr());

             // We set this to the default value of 3 but allow it be larger if th0 and qv0 were written
             std::istringstream isb(BaseStatefileCharPtrString, std::istringstream::in);
             isb >> ncomp_base_to_read;
             isb >> ng_base;
         }
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // NOTE: For backward compatibility (chk file has QKE)
        if ((chk_ncomp_cons-1)==ncomp_cons) {
            MultiFab cons(grids[lev],dmap[lev],chk_ncomp_cons,0);
            VisMF::Read(cons, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Cell"));

            // Copy up to RhoKE_comp
            MultiFab::Copy(vars_new[lev][Vars::cons],cons,0,0,(RhoKE_comp+1),0);

            // Only if we have a PBL model do we need to copy QKE is src to KE in dst
            if ( (solverChoice.turbChoice[lev].pbl_type == PBLType::MYNN25) ||
                 (solverChoice.turbChoice[lev].pbl_type == PBLType::MYNNEDMF) ) {
                MultiFab::Copy(vars_new[lev][Vars::cons],cons,(RhoKE_comp+1),RhoKE_comp,1,0);
                vars_new[lev][Vars::cons].mult(0.5,RhoKE_comp,1,0);
            }

            // Copy other components
            int ncomp_remainder = ncomp_cons - (RhoKE_comp + 1);
            MultiFab::Copy(vars_new[lev][Vars::cons],cons,(RhoKE_comp+2),(RhoKE_comp+1),ncomp_remainder,0);

            vars_new[lev][Vars::cons].setBndry(1.0e34);
        } else {
            MultiFab cons(grids[lev],dmap[lev],ncomp_cons,0);
            VisMF::Read(cons, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Cell"));
            MultiFab::Copy(vars_new[lev][Vars::cons],cons,0,0,ncomp_cons,0);
            vars_new[lev][Vars::cons].setBndry(1.0e34);
        }

        MultiFab xvel(convert(grids[lev],IntVect(1,0,0)),dmap[lev],1,0);
        VisMF::Read(xvel, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "XFace"));
        MultiFab::Copy(vars_new[lev][Vars::xvel],xvel,0,0,1,0);
        vars_new[lev][Vars::xvel].setBndry(1.0e34);

        MultiFab yvel(convert(grids[lev],IntVect(0,1,0)),dmap[lev],1,0);
        VisMF::Read(yvel, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "YFace"));
        MultiFab::Copy(vars_new[lev][Vars::yvel],yvel,0,0,1,0);
        vars_new[lev][Vars::yvel].setBndry(1.0e34);

        MultiFab zvel(convert(grids[lev],IntVect(0,0,1)),dmap[lev],1,0);
        VisMF::Read(zvel, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "ZFace"));
        MultiFab::Copy(vars_new[lev][Vars::zvel],zvel,0,0,1,0);
        vars_new[lev][Vars::zvel].setBndry(1.0e34);

        // Note that we read the ghost cells of the base state (unlike above)

        // The original base state only had 3 components and 1 ghost cell -- we read this
        // here to be consistent with the old style
        MultiFab base(grids[lev],dmap[lev],ncomp_base_to_read,ng_base);
        VisMF::Read(base, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "BaseState"));

        MultiFab::Copy(base_state[lev],base,0,0,ncomp_base_to_read,ng_base);

        // Create theta0 from p0, rh0
        if (ncomp_base_to_read < 4) {
            for (MFIter mfi(base_state[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                // We only compute theta_0 on valid cells since we will impose domain BC's after restart
                const Box& bx = mfi.tilebox();
                Array4<Real> const& fab = base_state[lev].array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    fab(i,j,k,BaseState::th0_comp) = getRhoThetagivenP(fab(i,j,k,BaseState::p0_comp))
                                                     / fab(i,j,k,BaseState::r0_comp);
                });
            }
        }
        // Default theta0 to 0
        if (ncomp_base_to_read < 5) {
            for (MFIter mfi(base_state[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                // We only compute theta_0 on valid cells since we will impose domain BC's after restart
                const Box& bx = mfi.tilebox();
                Array4<Real> const& fab = base_state[lev].array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    fab(i,j,k,BaseState::qv0_comp) = 0.0;
                });
            }
        }
        base_state[lev].FillBoundary(geom[lev].periodicity());

        if (SolverChoice::mesh_type != MeshType::ConstantDz)  {
           // Note that we also read the ghost cells of z_phys_nd
           IntVect ng = z_phys_nd[lev]->nGrowVect();
           MultiFab z_height(convert(grids[lev],IntVect(1,1,1)),dmap[lev],1,ng);
           VisMF::Read(z_height, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Z_Phys_nd"));
           MultiFab::Copy(*z_phys_nd[lev],z_height,0,0,1,ng);
           update_terrain_arrays(lev);
        }

        // Read in the moisture model restart variables
        std::vector<int> qmoist_indices;
        std::vector<std::string> qmoist_names;
        micro->Get_Qmoist_Restart_Vars(lev, solverChoice, qmoist_indices, qmoist_names);
        int qmoist_nvar = qmoist_indices.size();
        for (int var = 0; var < qmoist_nvar; var++) {
            const int ncomp  = 1;
            IntVect ng_moist = qmoist[lev][qmoist_indices[var]]->nGrowVect();
            MultiFab moist_vars(grids[lev],dmap[lev],ncomp,ng_moist);
            VisMF::Read(moist_vars, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", qmoist_names[var]));
            MultiFab::Copy(*(qmoist[lev][qmoist_indices[var]]),moist_vars,0,0,ncomp,ng_moist);
        }

#if defined(ERF_USE_WINDFARM)
        if(solverChoice.windfarm_type == WindFarmType::Fitch or
           solverChoice.windfarm_type == WindFarmType::EWP or
           solverChoice.windfarm_type == WindFarmType::SimpleAD){
            IntVect ng = Nturb[lev].nGrowVect();
            MultiFab mf_Nturb(grids[lev],dmap[lev],1,ng);
            VisMF::Read(mf_Nturb, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "NumTurb"));
            MultiFab::Copy(Nturb[lev],mf_Nturb,0,0,1,ng);
        }
#endif

        if (solverChoice.lsm_type != LandSurfaceType::None) {
            for (int mvar(0); mvar<lsm_data[lev].size(); ++mvar) {
                BoxArray ba = lsm_data[lev][mvar]->boxArray();
                DistributionMapping dm = lsm_data[lev][mvar]->DistributionMap();
                IntVect ng = lsm_data[lev][mvar]->nGrowVect();
                int nvar = lsm_data[lev][mvar]->nComp();
                MultiFab lsm_vars(ba,dm,nvar,ng);
                VisMF::Read(lsm_vars, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "LsmVars"));
                MultiFab::Copy(*(lsm_data[lev][mvar]),lsm_vars,0,0,nvar,ng);
            }
        }


        IntVect ng = mapfac[lev][MapFacType::mx]->nGrowVect();
        MultiFab mf_m(ba2d[lev],dmap[lev],1,ng);
        VisMF::Read(mf_m, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "MapFactor_mx"));
        MultiFab::Copy(*mapfac[lev][MapFacType::mx],mf_m,0,0,1,ng);

        if (MapFacType::mx != MapFacType::my) {
            VisMF::Read(mf_m, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "MapFactor_my"));
            MultiFab::Copy(*mapfac[lev][MapFacType::my],mf_m,0,0,1,ng);
        }

        ng = mapfac[lev][MapFacType::ux]->nGrowVect();
        MultiFab mf_u(convert(ba2d[lev],IntVect(1,0,0)),dmap[lev],1,ng);
        VisMF::Read(mf_u, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "MapFactor_ux"));
        MultiFab::Copy(*mapfac[lev][MapFacType::ux],mf_u,0,0,1,ng);

        if (MapFacType::ux != MapFacType::uy) {
            VisMF::Read(mf_u, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "MapFactor_uy"));
            MultiFab::Copy(*mapfac[lev][MapFacType::uy],mf_u,0,0,1,ng);
        }

        ng = mapfac[lev][MapFacType::vx]->nGrowVect();
        MultiFab mf_v(convert(ba2d[lev],IntVect(0,1,0)),dmap[lev],1,ng);
        VisMF::Read(mf_v, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "MapFactor_vx"));
        MultiFab::Copy(*mapfac[lev][MapFacType::vx],mf_v,0,0,1,ng);

        if (MapFacType::vx != MapFacType::vy) {
            VisMF::Read(mf_v, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "MapFactor_vy"));
            MultiFab::Copy(*mapfac[lev][MapFacType::vy],mf_v,0,0,1,ng);
        }


        // NOTE: We read MOST data in ReadCheckpointFileMOST (see below)!

        // See if we wrote out SST data
        std::string FirstSSTFileName(restart_chkfile + "/Level_0/SST_0_H");
        if (amrex::FileExists(FirstSSTFileName))
        {
            amrex::Print() << "Reading SST data" << std::endl;
            int ntimes = 1;
            ng = vars_new[lev][Vars::cons].nGrowVect(); ng[2]=0;
            MultiFab sst_at_t(ba2d[lev],dmap[lev],1,ng);
            sst_lev[lev][0] = std::make_unique<MultiFab>(ba2d[lev],dmap[lev],1,ng);
            for (int nt(0); nt<ntimes; ++nt) {
                VisMF::Read(sst_at_t, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_",
                                                             "SST_" + std::to_string(nt)));
                MultiFab::Copy(*sst_lev[lev][nt],sst_at_t,0,0,1,ng);
            }
        }

        // See if we wrote out TSK data
        std::string FirstTSKFileName(restart_chkfile + "/Level_0/TSK_0_H");
        if (amrex::FileExists(FirstTSKFileName))
        {
            amrex::Print() << "Reading TSK data" << std::endl;
            int ntimes = 1;
            ng = vars_new[lev][Vars::cons].nGrowVect(); ng[2]=0;
            MultiFab tsk_at_t(ba2d[lev],dmap[lev],1,ng);
            tsk_lev[lev][0] = std::make_unique<MultiFab>(ba2d[lev],dmap[lev],1,ng);
            for (int nt(0); nt<ntimes; ++nt) {
                VisMF::Read(tsk_at_t, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_",
                                                             "TSK_" + std::to_string(nt)));
                MultiFab::Copy(*tsk_lev[lev][nt],tsk_at_t,0,0,1,ng);
            }
        }

        {
            amrex::Print() << "Reading LMASK data" << std::endl;
            int ntimes = 1;
            ng = vars_new[lev][Vars::cons].nGrowVect(); ng[2]=0;
            MultiFab lmask_at_t(ba2d[lev],dmap[lev],1,ng);
            for (int nt(0); nt<ntimes; ++nt) {
                VisMF::Read(lmask_at_t, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_",
                                                              "LMASK_" + std::to_string(nt)));
                for (MFIter mfi(lmask_at_t); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.growntilebox();
                    Array4<int>  const& dst_arr = lmask_lev[lev][nt]->array(mfi);
                    Array4<Real> const& src_arr = lmask_at_t.array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        dst_arr(i,j,k) = int(src_arr(i,j,k));
                    });
                }
            }
        }

#ifdef ERF_USE_NETCDF
        IntVect ngv = ng; ngv[2] = 0;

        // Read lat/lon if it exists
        if (solverChoice.has_lat_lon) {
            amrex::Print() << "Reading Lat/Lon variables" << std::endl;
            MultiFab lat(ba2d[lev],dmap[lev],1,ngv);
            MultiFab lon(ba2d[lev],dmap[lev],1,ngv);
            VisMF::Read(lat, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "LAT"));
            VisMF::Read(lon, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "LON"));
            lat_m[lev] = std::make_unique<MultiFab>(ba2d[lev],dmap[lev],1,ngv);
            lon_m[lev] = std::make_unique<MultiFab>(ba2d[lev],dmap[lev],1,ngv);
            MultiFab::Copy(*lat_m[lev],lat,0,0,1,ngv);
            MultiFab::Copy(*lon_m[lev],lon,0,0,1,ngv);
        }

        // Read sinPhi and cosPhi if it exists
        if (solverChoice.variable_coriolis) {
            amrex::Print() << "Reading Coriolis factors" << std::endl;
            MultiFab sphi(ba2d[lev],dmap[lev],1,ngv);
            MultiFab cphi(ba2d[lev],dmap[lev],1,ngv);
            VisMF::Read(sphi, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "SinPhi"));
            VisMF::Read(cphi, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "CosPhi"));
            sinPhi_m[lev] = std::make_unique<MultiFab>(ba2d[lev],dmap[lev],1,ngv);
            cosPhi_m[lev] = std::make_unique<MultiFab>(ba2d[lev],dmap[lev],1,ngv);
            MultiFab::Copy(*sinPhi_m[lev],sphi,0,0,1,ngv);
            MultiFab::Copy(*cosPhi_m[lev],cphi,0,0,1,ngv);
        }

        if (solverChoice.use_real_bcs && solverChoice.init_type == InitType::WRFInput) {
            MultiFab tmp1d(ba1d[lev],dmap[lev],1,ngv);
            MultiFab tmp2d(ba2d[lev],dmap[lev],1,ngv);

            VisMF::Read(tmp1d, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "C1H"));
            MultiFab::Copy(*mf_C1H[lev],tmp1d,0,0,1,ngv);

            VisMF::Read(tmp1d, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "C2H"));
            MultiFab::Copy(*mf_C2H[lev],tmp1d,0,0,1,ngv);

            VisMF::Read(tmp2d, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "MUB"));
            MultiFab::Copy(*mf_MUB[lev],tmp2d,0,0,1,ngv);
        }
#endif

    } // for lev

#ifdef ERF_USE_PARTICLES
    restartTracers((ParGDBBase*)GetParGDB(),restart_chkfile);
    if (Microphysics::modelType(solverChoice.moisture_type) == MoistureModelType::Lagrangian) {
        dynamic_cast<LagrangianMicrophysics&>(*micro).restartParticles((ParGDBBase*)GetParGDB(),restart_chkfile);
    }
#endif

#if 0
#ifdef ERF_USE_NETCDF
    // Read bdy_data files
    if ( ((solverChoice.init_type==InitType::WRFInput) || (solverChoice.init_type==InitType::Metgrid)) &&
         solverChoice.use_real_bcs )
    {
        int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank
        int num_time;
        int num_var;
        Vector<Box> bx_v;
        if (ParallelDescriptor::IOProcessor()) {
            // Open header file and read from it
            std::ifstream bdy_h_file(MultiFabFileFullPrefix(0, restart_chkfile, "Level_", "bdy_H"));
            bdy_h_file >> num_time;
            bdy_h_file >> num_var;
            bdy_h_file >> start_bdy_time;
            bdy_h_file >> bdy_time_interval;
            bdy_h_file >> real_width;
            bx_v.resize(4*num_var);
            for (int ivar(0); ivar<num_var; ++ivar) {
                bdy_h_file >> bx_v[4*ivar  ];
                bdy_h_file >> bx_v[4*ivar+1];
                bdy_h_file >> bx_v[4*ivar+2];
                bdy_h_file >> bx_v[4*ivar+3];
            }

            // IO size the FABs
            bdy_data_xlo.resize(num_time);
            bdy_data_xhi.resize(num_time);
            bdy_data_ylo.resize(num_time);
            bdy_data_yhi.resize(num_time);
            for (int itime(0); itime<num_time; ++itime) {
                bdy_data_xlo[itime].resize(num_var);
                bdy_data_xhi[itime].resize(num_var);
                bdy_data_ylo[itime].resize(num_var);
                bdy_data_yhi[itime].resize(num_var);
                for (int ivar(0); ivar<num_var; ++ivar) {
                    bdy_data_xlo[itime][ivar].resize(bx_v[4*ivar  ]);
                    bdy_data_xhi[itime][ivar].resize(bx_v[4*ivar+1]);
                    bdy_data_ylo[itime][ivar].resize(bx_v[4*ivar+2]);
                    bdy_data_yhi[itime][ivar].resize(bx_v[4*ivar+3]);
                }
            }

            // Open data file and read from it
            std::ifstream bdy_d_file(MultiFabFileFullPrefix(0, restart_chkfile, "Level_", "bdy_D"));
            for (int itime(0); itime<num_time; ++itime) {
                for (int ivar(0); ivar<num_var; ++ivar) {
                    bdy_data_xlo[itime][ivar].readFrom(bdy_d_file);
                    bdy_data_xhi[itime][ivar].readFrom(bdy_d_file);
                    bdy_data_ylo[itime][ivar].readFrom(bdy_d_file);
                    bdy_data_yhi[itime][ivar].readFrom(bdy_d_file);
                }
            }
        } // IO

        // Broadcast the data
        ParallelDescriptor::Barrier();
        ParallelDescriptor::Bcast(&start_bdy_time,1,ioproc);
        ParallelDescriptor::Bcast(&bdy_time_interval,1,ioproc);
        ParallelDescriptor::Bcast(&real_width,1,ioproc);
        ParallelDescriptor::Bcast(&num_time,1,ioproc);
        ParallelDescriptor::Bcast(&num_var,1,ioproc);

        // Everyone size their boxes
        bx_v.resize(4*num_var);

        ParallelDescriptor::Bcast(bx_v.dataPtr(),bx_v.size(),ioproc);

        // Everyone but IO size their FABs
        if (!ParallelDescriptor::IOProcessor()) {
          bdy_data_xlo.resize(num_time);
          bdy_data_xhi.resize(num_time);
          bdy_data_ylo.resize(num_time);
          bdy_data_yhi.resize(num_time);
          for (int itime(0); itime<num_time; ++itime) {
            bdy_data_xlo[itime].resize(num_var);
            bdy_data_xhi[itime].resize(num_var);
            bdy_data_ylo[itime].resize(num_var);
            bdy_data_yhi[itime].resize(num_var);
            for (int ivar(0); ivar<num_var; ++ivar) {
              bdy_data_xlo[itime][ivar].resize(bx_v[4*ivar  ]);
              bdy_data_xhi[itime][ivar].resize(bx_v[4*ivar+1]);
              bdy_data_ylo[itime][ivar].resize(bx_v[4*ivar+2]);
              bdy_data_yhi[itime][ivar].resize(bx_v[4*ivar+3]);
            }
          }
        }

        for (int itime(0); itime<num_time; ++itime) {
            for (int ivar(0); ivar<num_var; ++ivar) {
                ParallelDescriptor::Bcast(bdy_data_xlo[itime][ivar].dataPtr(),bdy_data_xlo[itime][ivar].box().numPts(),ioproc);
                ParallelDescriptor::Bcast(bdy_data_xhi[itime][ivar].dataPtr(),bdy_data_xhi[itime][ivar].box().numPts(),ioproc);
                ParallelDescriptor::Bcast(bdy_data_ylo[itime][ivar].dataPtr(),bdy_data_ylo[itime][ivar].box().numPts(),ioproc);
                ParallelDescriptor::Bcast(bdy_data_yhi[itime][ivar].dataPtr(),bdy_data_yhi[itime][ivar].box().numPts(),ioproc);
            }
        }
    } // init_type == WRFInput or Metgrid
#endif
#endif
}

/**
 * ERF function for reading additional data for MOST from a checkpoint file during restart.
 *
 * This is called after the ABLMost object is instantiated.
 */
void
ERF::ReadCheckpointFileSurfaceLayer ()
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        amrex::Print() << "Reading MOST variables" << std::endl;

        IntVect ng = vars_new[lev][Vars::cons].nGrowVect(); ng[2]=0;
        MultiFab  m_var(ba2d[lev],dmap[lev],1,ng);
        MultiFab* dst = nullptr;

        // U*
        dst = m_SurfaceLayer->get_u_star(lev);
        VisMF::Read(m_var, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Ustar"));
        MultiFab::Copy(*dst,m_var,0,0,1,ng);

        // W*
        dst = m_SurfaceLayer->get_w_star(lev);
        VisMF::Read(m_var, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Wstar"));
        MultiFab::Copy(*dst,m_var,0,0,1,ng);

        // T*
        dst = m_SurfaceLayer->get_t_star(lev);
        VisMF::Read(m_var, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Tstar"));
        MultiFab::Copy(*dst,m_var,0,0,1,ng);

        // Q*
        dst = m_SurfaceLayer->get_q_star(lev);
        VisMF::Read(m_var, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Qstar"));
        MultiFab::Copy(*dst,m_var,0,0,1,ng);

        // Olen
        dst = m_SurfaceLayer->get_olen(lev);
        VisMF::Read(m_var, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Olen"));
        MultiFab::Copy(*dst,m_var,0,0,1,ng);

        // PBLH
        dst = m_SurfaceLayer->get_pblh(lev);
        VisMF::Read(m_var, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "PBLH"));
        MultiFab::Copy(*dst,m_var,0,0,1,ng);

        // Z0
        VisMF::Read(m_var, MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Z0"));
        for (amrex::MFIter mfi(m_var); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.growntilebox();
            FArrayBox* most_z0 = (m_SurfaceLayer->get_z0(lev));
            most_z0->copy<RunOn::Device>(m_var[mfi], bx);
        }
    }
}
