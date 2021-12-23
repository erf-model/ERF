#include <ERF.H>

// utility to skip to next line in Header
void
ERF::GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}

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
   for (int lev = 0; lev <= finest_level; ++lev) {

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
   }
}

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

    int ngrow_state = ComputeGhostCells(solverChoice.spatial_order)+1;
    int ngrow_vels  = ComputeGhostCells(solverChoice.spatial_order);

    for (int lev = 0; lev <= finest_level; ++lev) {

        // read in level 'lev' BoxArray from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // create a distribution mapping
        DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

        // set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        // build MultiFab data
        int ncomp = Cons::NumVars;

        auto& lev_old = vars_old[lev];
        auto& lev_new = vars_new[lev];

        lev_new[Vars::cons].define(grids[lev], dmap[lev], ncomp, ngrow_state);
        lev_old[Vars::cons].define(grids[lev], dmap[lev], ncomp, ngrow_state);

        //!don: get the ghost cells right here
        lev_new[Vars::xvel].define(convert(grids[lev], IntVect(1,0,0)), dmap[lev], 1, ngrow_vels);
        lev_old[Vars::xvel].define(convert(grids[lev], IntVect(1,0,0)), dmap[lev], 1, ngrow_vels);

        lev_new[Vars::yvel].define(convert(grids[lev], IntVect(0,1,0)), dmap[lev], 1, ngrow_vels);
        lev_old[Vars::yvel].define(convert(grids[lev], IntVect(0,1,0)), dmap[lev], 1, ngrow_vels);

        lev_new[Vars::zvel].define(convert(grids[lev], IntVect(0,0,1)), dmap[lev], 1, ngrow_vels);
        lev_old[Vars::zvel].define(convert(grids[lev], IntVect(0,0,1)), dmap[lev], 1, ngrow_vels);
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

        // Copy from new into old just in case
        MultiFab::Copy(vars_old[lev][Vars::cons],vars_new[lev][Vars::cons],0,0,NVAR,0);
        MultiFab::Copy(vars_old[lev][Vars::xvel],vars_new[lev][Vars::xvel],0,0,1,0);
        MultiFab::Copy(vars_old[lev][Vars::yvel],vars_new[lev][Vars::yvel],0,0,1,0);
        MultiFab::Copy(vars_old[lev][Vars::zvel],vars_new[lev][Vars::zvel],0,0,1,0);
    }
}
