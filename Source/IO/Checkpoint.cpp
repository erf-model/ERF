#include <ERF.H>
#include "AMReX_PlotFileUtil.H"

#include <iostream>
#include <fstream>

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
    const std::string& checkpointname = amrex::Concatenate(check_file,istep[0],5);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    const int nlevels = finest_level+1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then build
    // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    // write Header file
   if (ParallelDescriptor::IOProcessor()) {

       std::string HeaderFileName(checkpointname + "/Header");
       VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
       std::ofstream HeaderFile;
       HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
       HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                               std::ofstream::trunc |
                                               std::ofstream::binary);
       if( ! HeaderFile.good()) {
           amrex::FileOpenFailed(HeaderFileName);
       }

       HeaderFile.precision(17);

       // write out title line
       HeaderFile << "Checkpoint file for ERF\n";

       // write out finest_level
       HeaderFile << finest_level << "\n";

       // write the number of components
       // for each variable we store

       // conservative, cell-centered vars
       HeaderFile << Cons::NumVars << "\n";

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
   }

   // write the MultiFab data to, e.g., chk00010/Level_0/
   // Here we make copies of the MultiFab with no ghost cells
   for (int lev = 0; lev <= finest_level; ++lev)
   {
       MultiFab cons(grids[lev],dmap[lev],Cons::NumVars,0);
       MultiFab::Copy(cons,vars_new[lev][Vars::cons],0,0,NVAR,0);
       VisMF::Write(cons, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Cell"));

       MultiFab xvel(convert(grids[lev],IntVect(1,0,0)),dmap[lev],1,0);
       MultiFab::Copy(xvel,vars_new[lev][Vars::xvel],0,0,1,0);
       VisMF::Write(xvel, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "XFace"));

       MultiFab yvel(convert(grids[lev],IntVect(0,1,0)),dmap[lev],1,0);
       MultiFab::Copy(yvel,vars_new[lev][Vars::yvel],0,0,1,0);
       VisMF::Write(yvel, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "YFace"));

       MultiFab zvel(convert(grids[lev],IntVect(0,0,1)),dmap[lev],1,0);
       MultiFab::Copy(zvel,vars_new[lev][Vars::zvel],0,0,1,0);
       VisMF::Write(zvel, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "ZFace"));

       IntVect ng;
#ifdef ERF_USE_MOISTURE
       // We must read and write qmoist with ghost cells because we don't directly impose BCs on these vars
       ng = qmoist[lev].nGrowVect();
       MultiFab moist_vars(grids[lev],dmap[lev],qmoist[lev].nComp(),ng);
       MultiFab::Copy(moist_vars,qmoist[lev],0,0,qmoist[lev].nComp(),ng);
       VisMF::Write(moist_vars, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "MoistVars"));
#endif

       // Note that we write the ghost cells of the base state (unlike above)
       ng = base_state[lev].nGrowVect();
       MultiFab base(grids[lev],dmap[lev],base_state[lev].nComp(),ng);
       MultiFab::Copy(base,base_state[lev],0,0,base.nComp(),ng);
       VisMF::Write(base, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "BaseState"));

       if (solverChoice.use_terrain)  {
           // Note that we also write the ghost cells of z_phys_nd
           ng = z_phys_nd[lev]->nGrowVect();
           MultiFab z_height(convert(grids[lev],IntVect(1,1,1)),dmap[lev],1,ng);
           MultiFab::Copy(z_height,*z_phys_nd[lev],0,0,1,ng);
           VisMF::Write(z_height, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Z_Phys_nd"));
       }

       // Note that we also write the ghost cells of the mapfactors (2D)
       BoxList bl2d = grids[lev].boxList();
       for (auto& b : bl2d) {
           b.setRange(2,0);
       }
       BoxArray ba2d(std::move(bl2d));

       ng = mapfac_m[lev]->nGrowVect();
       MultiFab mf_m(ba2d,dmap[lev],1,ng);
       MultiFab::Copy(mf_m,*mapfac_m[lev],0,0,1,ng);
       VisMF::Write(mf_m, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "MapFactor_m"));

       ng = mapfac_u[lev]->nGrowVect();
       MultiFab mf_u(convert(ba2d,IntVect(1,0,0)),dmap[lev],1,ng);
       MultiFab::Copy(mf_u,*mapfac_u[lev],0,0,1,ng);
       VisMF::Write(mf_u, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "MapFactor_u"));

       ng = mapfac_v[lev]->nGrowVect();
       MultiFab mf_v(convert(ba2d,IntVect(0,1,0)),dmap[lev],1,ng);
       MultiFab::Copy(mf_v,*mapfac_v[lev],0,0,1,ng);
       VisMF::Write(mf_v, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "MapFactor_v"));
   }

#ifdef ERF_USE_PARTICLES
   if (use_tracer_particles) {
       tracer_particles->Checkpoint(checkpointname, "tracers", true, tracer_particle_varnames);
   }
#endif

#ifdef ERF_USE_NETCDF
   // Write bdy_data files
   if (ParallelDescriptor::IOProcessor() && (init_type == "real")) {

     // Vector dimensions
     int num_time = bdy_data_xlo.size();
     int num_var  = bdy_data_xlo[0].size();

     // Open header file and write to it
     std::ofstream bdy_h_file(amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "bdy_H"));
     bdy_h_file << std::setprecision(1) << std::fixed;
     bdy_h_file << num_time << "\n";
     bdy_h_file << num_var  << "\n";
     bdy_h_file << start_bdy_time << "\n";
     bdy_h_file << bdy_time_interval << "\n";
     bdy_h_file << wrfbdy_width << "\n";
     for (int ivar(0); ivar<num_var; ++ivar) {
       bdy_h_file << bdy_data_xlo[0][ivar].box() << "\n";
       bdy_h_file << bdy_data_xhi[0][ivar].box() << "\n";
       bdy_h_file << bdy_data_ylo[0][ivar].box() << "\n";
       bdy_h_file << bdy_data_yhi[0][ivar].box() << "\n";
     }

     // Open data file and write to it
     std::ofstream bdy_d_file(amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "bdy_D"));
     for (int itime(0); itime<num_time; ++itime) {
       for (int ivar(0); ivar<num_var; ++ivar) {
         bdy_data_xlo[itime][ivar].writeOn(bdy_d_file,0,1);
         bdy_data_xhi[itime][ivar].writeOn(bdy_d_file,0,1);
         bdy_data_ylo[itime][ivar].writeOn(bdy_d_file,0,1);
         bdy_data_yhi[itime][ivar].writeOn(bdy_d_file,0,1);
       }
     }
   }
#endif

}

/**
 * ERF function for reading data from a checkpoint file during restart.
 */
void
ERF::ReadCheckpointFile ()
{
    amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

    // Header
    std::string File(restart_chkfile + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    int chk_ncomp;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level;
    GotoNextLine(is);

    // read the number of components
    // for each variable we store

    // conservative, cell-centered vars
    is >> chk_ncomp;
    GotoNextLine(is);
    AMREX_ASSERT(chk_ncomp == Cons::NumVars);

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

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        MultiFab cons(grids[lev],dmap[lev],Cons::NumVars,0);
        VisMF::Read(cons, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Cell"));
        MultiFab::Copy(vars_new[lev][Vars::cons],cons,0,0,Cons::NumVars,0);

        MultiFab xvel(convert(grids[lev],IntVect(1,0,0)),dmap[lev],1,0);
        VisMF::Read(xvel, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "XFace"));
        MultiFab::Copy(vars_new[lev][Vars::xvel],xvel,0,0,1,0);

        MultiFab yvel(convert(grids[lev],IntVect(0,1,0)),dmap[lev],1,0);
        VisMF::Read(yvel, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "YFace"));
        MultiFab::Copy(vars_new[lev][Vars::yvel],yvel,0,0,1,0);

        MultiFab zvel(convert(grids[lev],IntVect(0,0,1)),dmap[lev],1,0);
        VisMF::Read(zvel, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "ZFace"));
        MultiFab::Copy(vars_new[lev][Vars::zvel],zvel,0,0,1,0);

        IntVect ng;
#ifdef ERF_USE_MOISTURE
        ng = qmoist[lev].nGrowVect();
        MultiFab moist_vars(grids[lev],dmap[lev],qmoist[lev].nComp(),ng);
        VisMF::Read(moist_vars, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "MoistVars"));
        MultiFab::Copy(qmoist[lev],moist_vars,0,0,qmoist[lev].nComp(),ng);
#endif

        // Note that we read the ghost cells of the base state (unlike above)
        ng = base_state[lev].nGrowVect();
        MultiFab base(grids[lev],dmap[lev],base_state[lev].nComp(),ng);
        VisMF::Read(base, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "BaseState"));
        MultiFab::Copy(base_state[lev],base,0,0,base.nComp(),ng);
        base_state[lev].FillBoundary(geom[lev].periodicity());

        if (solverChoice.use_terrain)  {
           // Note that we also read the ghost cells of z_phys_nd
           ng = z_phys_nd[lev]->nGrowVect();
           MultiFab z_height(convert(grids[lev],IntVect(1,1,1)),dmap[lev],1,ng);
           VisMF::Read(z_height, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Z_Phys_nd"));
           MultiFab::Copy(*z_phys_nd[lev],z_height,0,0,1,ng);
           update_terrain_arrays(lev, t_new[lev]);
        }

        // Note that we read the ghost cells of the mapfactors
        BoxList bl2d = grids[lev].boxList();
        for (auto& b : bl2d) {
            b.setRange(2,0);
        }
        BoxArray ba2d(std::move(bl2d));

        ng = mapfac_m[lev]->nGrowVect();
        MultiFab mf_m(ba2d,dmap[lev],1,ng);
        VisMF::Read(mf_m, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "MapFactor_m"));
        MultiFab::Copy(*mapfac_m[lev],mf_m,0,0,1,ng);

        ng = mapfac_u[lev]->nGrowVect();
        MultiFab mf_u(convert(ba2d,IntVect(1,0,0)),dmap[lev],1,ng);
        VisMF::Read(mf_u, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "MapFactor_u"));
        MultiFab::Copy(*mapfac_u[lev],mf_u,0,0,1,ng);

        ng = mapfac_v[lev]->nGrowVect();
        MultiFab mf_v(convert(ba2d,IntVect(0,1,0)),dmap[lev],1,ng);
        VisMF::Read(mf_v, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "MapFactor_v"));
        MultiFab::Copy(*mapfac_v[lev],mf_v,0,0,1,ng);
    }

#ifdef ERF_USE_PARTICLES
   if (use_tracer_particles) {
       tracer_particles = std::make_unique<TerrainFittedPC>(Geom(0), dmap[0], grids[0]);
       std::string tracer_file("tracers");
       tracer_particles->Restart(restart_chkfile, tracer_file);
   }
#endif

#ifdef ERF_USE_NETCDF
    // Read bdy_data files
    if (init_type == "real") {
        int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank
        int num_time;
        int num_var;
        Vector<Box> bx_v;
        if (ParallelDescriptor::IOProcessor()) {
            // Open header file and read from it
            std::ifstream bdy_h_file(amrex::MultiFabFileFullPrefix(0, restart_chkfile, "Level_", "bdy_H"));
            bdy_h_file >> num_time;
            bdy_h_file >> num_var;
            bdy_h_file >> start_bdy_time;
            bdy_h_file >> bdy_time_interval;
            bdy_h_file >> wrfbdy_width;
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
            std::ifstream bdy_d_file(amrex::MultiFabFileFullPrefix(0, restart_chkfile, "Level_", "bdy_D"));
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
        ParallelDescriptor::Bcast(&wrfbdy_width,1,ioproc);
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
    } // init real
#endif
}
