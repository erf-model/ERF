#include <ERF.H>
#include <NCInterface.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

/**
 * Writes a checkpoint file in NetCDF format
 */
void
ERF::WriteNCCheckpointFile () const
{
    // checkpoint file name, e.g., chk00010
    const std::string& checkpointname = amrex::Concatenate(check_file,istep[0],5);

    amrex::Print() << "Writing NetCDF checkpoint " << checkpointname << "\n";

    const int nlevels = finest_level+1;

    // ---- ParallelDescriptor::IOProcessor() creates the directories
    PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    // write Header file
   if (ParallelDescriptor::IOProcessor()) {

       std::string HeaderFileName(checkpointname + "/Header.nc");

       auto ncf = ncutils::NCFile::create(HeaderFileName, NC_CLOBBER | NC_NETCDF4);

       const std::string ndim_name  = "num_dimension";
       const std::string nl_name    = "finest_levels";
       const std::string ng_name    = "num_grids";
       const std::string nvar_name  = "num_vars";
       const std::string ndt_name   = "num_dt";
       const std::string nstep_name = "num_istep";
       const std::string ntime_name = "num_newtime";

       const int ndt   = dt.size();
       const int nstep = istep.size();
       const int ntime = t_new.size();

       amrex::Vector<int> nbox(nlevels);
       amrex::Vector<std::string> nbox_name(nlevels);
       amrex::Vector<amrex::Vector<std::string> > lo_names(nlevels);
       amrex::Vector<amrex::Vector<std::string> > hi_names(nlevels);
       amrex::Vector<amrex::Vector<std::string> > typ_names(nlevels);
       for (auto lev{0}; lev <= finest_level; ++lev) {
           nbox[lev]      = boxArray(lev).size();
           nbox_name[lev] = "NBox_"+std::to_string(lev);
           for (int nb(0); nb < boxArray(lev).size(); ++nb) {
              lo_names[lev] .push_back("SmallEnd_"+std::to_string(lev)+"_"+std::to_string(nb));
              hi_names[lev] .push_back("BigEnd_"+std::to_string(lev)+"_"+std::to_string(nb));
              typ_names[lev].push_back("BoxType_"+std::to_string(lev)+"_"+std::to_string(nb));
           }
       }

       ncf.enter_def_mode();
       ncf.put_attr("title", "ERF NetCDF CheckPoint Header");

       ncf.def_dim(ndim_name,  AMREX_SPACEDIM);
       ncf.def_dim(nl_name,    nlevels);
       ncf.def_dim(ndt_name,   ndt);
       ncf.def_dim(nstep_name, nstep);
       ncf.def_dim(ntime_name, ntime);

       for (auto lev{0}; lev <= finest_level; ++lev) {
           ncf.def_dim(nbox_name[lev], nbox[lev]);
           for (int nb(0); nb < boxArray(lev).size(); ++nb) {
              ncf.def_var(lo_names[lev][nb],  ncutils::NCDType::Int, {nbox_name[lev], ndim_name});
              ncf.def_var(hi_names[lev][nb],  ncutils::NCDType::Int, {nbox_name[lev], ndim_name});
              ncf.def_var(typ_names[lev][nb], ncutils::NCDType::Int, {nbox_name[lev], ndim_name});
           }
       }

       ncf.def_var("istep", ncutils::NCDType::Int,  {nstep_name});
       ncf.def_var("dt"   , ncutils::NCDType::Real, {ndt_name}  );
       ncf.def_var("tnew" , ncutils::NCDType::Real, {ntime_name});

       ncf.exit_def_mode();

       // output headfile in NetCDF format
       ncf.var("istep").put(istep.data(), {0}, {static_cast<long unsigned int>(nstep)});
       ncf.var("dt")   .put(dt.data(),    {0}, {static_cast<long unsigned int>(ndt)});
       ncf.var("tnew") .put(t_new.data(), {0}, {static_cast<long unsigned int>(ntime)});
//       ncf.var("nbox") .put(nbox.begin(),  {0}, {nstep});

       for (auto lev{0}; lev <= finest_level; ++lev) {
           auto box_array = boxArray(lev);
           for (int nb(0); nb < box_array.size(); ++nb) {
              auto nbb = static_cast<long unsigned int>(nb);
              auto box = box_array[nb];
              ncf.var(lo_names[lev][nb] ).put(box.smallEnd().begin(), {nbb, 0}, {1, AMREX_SPACEDIM});
              ncf.var(hi_names[lev][nb] ).put(box.bigEnd().begin()  , {nbb, 0}, {1, AMREX_SPACEDIM});
              ncf.var(typ_names[lev][nb]).put(box.type().begin()    , {nbb, 0}, {1, AMREX_SPACEDIM});
           }
       }
   }

   // write the MultiFab data to, e.g., chk00010/Level_0/
   // Here we make copies of the MultiFab with no ghost cells
   for (int lev = 0; lev <= finest_level; ++lev) {

       MultiFab cons(grids[lev],dmap[lev],Cons::NumVars,0);
       MultiFab::Copy(cons,vars_new[lev][Vars::cons],0,0,NVAR,0);
       WriteNCMultiFab(cons, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Cell"));

       MultiFab xvel(convert(grids[lev],IntVect(1,0,0)),dmap[lev],1,0);
       MultiFab::Copy(xvel,vars_new[lev][Vars::xvel],0,0,1,0);
       WriteNCMultiFab(xvel, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "XFace"));

       MultiFab yvel(convert(grids[lev],IntVect(0,1,0)),dmap[lev],1,0);
       MultiFab::Copy(yvel,vars_new[lev][Vars::yvel],0,0,1,0);
       WriteNCMultiFab(yvel, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "YFace"));

       MultiFab zvel(convert(grids[lev],IntVect(0,0,1)),dmap[lev],1,0);
       MultiFab::Copy(zvel,vars_new[lev][Vars::zvel],0,0,1,0);
       WriteNCMultiFab(zvel, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "ZFace"));
   }
}

/**
 * Read NetCDF checkpoint to restart ERF
 */
void
ERF::ReadNCCheckpointFile ()
{
    amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

    // Header
    std::string HeaderFileName(restart_chkfile + "/Header.nc");

    auto ncf = ncutils::NCFile::open(HeaderFileName, NC_CLOBBER | NC_NETCDF4);

    const std::string nl_name    = "finest_levels";
    const std::string ng_name    = "num_grids";
    const std::string nvar_name  = "num_vars";
    const std::string ndt_name   = "num_dt";
    const std::string nstep_name = "num_istep";
    const std::string ntime_name = "num_newtime";

    const int nvar         = static_cast<int>(ncf.dim(nvar_name).len());
    AMREX_ALWAYS_ASSERT(nvar == Cons::NumVars);

    const int ndt          = static_cast<int>(ncf.dim(ndt_name).len());
    const int nstep        = static_cast<int>(ncf.dim(nstep_name).len());
    const int ntime        = static_cast<int>(ncf.dim(ntime_name).len());

    // Assert we are reading in data with the same finest_level as we have
    AMREX_ALWAYS_ASSERT(finest_level == static_cast<int>(ncf.dim(nl_name).len()));

    // output headfile in NetCDF format
    ncf.var("istep").get(istep.data(), {0}, {static_cast<long unsigned int>(nstep)});
    ncf.var("dt")   .get(dt.data(),    {0}, {static_cast<long unsigned int>(ndt)});
    ncf.var("t_new").get(t_new.data(), {0}, {static_cast<long unsigned int>(ntime)});

    int ngrow_state = ComputeGhostCells(solverChoice) + 1;
    int ngrow_vels  = ComputeGhostCells(solverChoice);

    for (int lev = 0; lev <= finest_level; ++lev) {

        int num_box = static_cast<int>(ncf.dim("Nbox_"+std::to_string(lev)).len());

        // read in level 'lev' BoxArray from Header
        BoxArray ba;

        for (int nb(0); nb < num_box; ++nb) {
           amrex::IntVect lo(AMREX_SPACEDIM);
           amrex::IntVect hi(AMREX_SPACEDIM);
           amrex::IntVect typ(AMREX_SPACEDIM);

           auto lo_name  = "SmallEnd_"+std::to_string(lev)+"_"+std::to_string(nb);
           auto hi_name  = "BigEnd_"+std::to_string(lev)+"_"+std::to_string(nb);
           auto typ_name = "BoxType_"+std::to_string(lev)+"_"+std::to_string(nb);

           ncf.var(lo_name) .get(lo.begin(), {0}, {AMREX_SPACEDIM});
           ncf.var(hi_name) .get(hi.begin(), {0}, {AMREX_SPACEDIM});
           ncf.var(typ_name).get(typ.begin(),{0}, {AMREX_SPACEDIM});

           amrex::Box box = amrex::Box(lo, hi, typ);
           ba.set(nb, box);
        }

//        ba.readFrom(is);
//        GotoNextLine(is);

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
        WriteNCMultiFab(cons, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Cell"));
        MultiFab::Copy(vars_new[lev][Vars::cons],cons,0,0,Cons::NumVars,0);

        MultiFab xvel(convert(grids[lev],IntVect(1,0,0)),dmap[lev],1,0);
        WriteNCMultiFab(xvel, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Cell"));
        MultiFab::Copy(vars_new[lev][Vars::xvel],xvel,0,0,1,0);

        MultiFab yvel(convert(grids[lev],IntVect(0,1,0)),dmap[lev],1,0);
        WriteNCMultiFab(yvel, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Cell"));
        MultiFab::Copy(vars_new[lev][Vars::yvel],yvel,0,0,1,0);

        MultiFab zvel(convert(grids[lev],IntVect(0,0,1)),dmap[lev],1,0);
        WriteNCMultiFab(zvel, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Cell"));
        MultiFab::Copy(vars_new[lev][Vars::zvel],zvel,0,0,1,0);

        // Copy from new into old just in case
        MultiFab::Copy(vars_old[lev][Vars::cons],vars_new[lev][Vars::cons],0,0,NVAR,0);
        MultiFab::Copy(vars_old[lev][Vars::xvel],vars_new[lev][Vars::xvel],0,0,1,0);
        MultiFab::Copy(vars_old[lev][Vars::yvel],vars_new[lev][Vars::yvel],0,0,1,0);
        MultiFab::Copy(vars_old[lev][Vars::zvel],vars_new[lev][Vars::zvel],0,0,1,0);
    }
}
