/**
 * \file ERF.cpp
 */

#include <prob_common.H>
#include <EOS.H>
#include <ERF.H>

#include <AMReX_buildInfo.H>

#include <TerrainMetrics.H>

#ifdef ERF_USE_MULTIBLOCK
#include <MultiBlockContainer.H>
#endif

using namespace amrex;

amrex::Real ERF::startCPUTime        = 0.0;
amrex::Real ERF::previousCPUTimeUsed = 0.0;

Vector<AMRErrorTag> ERF::ref_tags;

SolverChoice ERF::solverChoice;

// Time step control
amrex::Real ERF::cfl           =  0.8;
amrex::Real ERF::fixed_dt      = -1.0;
amrex::Real ERF::fixed_fast_dt = -1.0;
amrex::Real ERF::init_shrink   =  1.0;
amrex::Real ERF::change_max    =  1.1;
int         ERF::fixed_mri_dt_ratio = 0;
bool        ERF::use_lowM_dt = false;

// Type of mesh refinement algorithm
std::string ERF::coupling_type = "OneWay";
int         ERF::do_reflux     = 0;
int         ERF::do_avg_down   = 0;

// Dictate verbosity in screen output
int         ERF::verbose       = 0;

// Use the native ERF MRI integrator
int         ERF::use_native_mri = 1;

// Frequency of diagnostic output
int         ERF::sum_interval  = -1;
amrex::Real ERF::sum_per       = -1.0;

// Native AMReX vs NetCDF
std::string ERF::plotfile_type    = "amrex";

// init_type:  "custom", "ideal", "real", "input_sounding"
std::string ERF::init_type        = "custom";

// NetCDF wrfinput (initialization) file(s)
amrex::Vector<amrex::Vector<std::string>> ERF::nc_init_file = {{""}}; // Must provide via input

// NetCDF wrfbdy (lateral boundary) file
std::string ERF::nc_bdy_file = ""; // Must provide via input

// Text input_sounding file
std::string ERF::input_sounding_file = ""; // Must provide via input

// 1D NetCDF output (for ingestion by AMR-Wind)
int         ERF::output_1d_column = 0;
int         ERF::column_interval  = -1;
amrex::Real ERF::column_per       = -1.0;
amrex::Real ERF::column_loc_x     = 0.0;
amrex::Real ERF::column_loc_y     = 0.0;
std::string ERF::column_file_name = "column_data.nc";

// 2D BndryRegister output (for ingestion by AMR-Wind)
int         ERF::output_bndry_planes            = 0;
int         ERF::bndry_output_planes_interval   = -1;
amrex::Real ERF::bndry_output_planes_per        = -1.0;
amrex::Real ERF::bndry_output_planes_start_time =  0.0;

// 2D BndryRegister input
int         ERF::input_bndry_planes             = 0;

amrex::Vector<std::string> BCNames = {"xlo", "ylo", "zlo", "xhi", "yhi", "zhi"};

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
ERF::ERF ()
{
    if (amrex::ParallelDescriptor::IOProcessor()) {
        const char* erf_hash = amrex::buildInfoGetGitHash(1);
        const char* amrex_hash = amrex::buildInfoGetGitHash(2);
        const char* buildgithash = amrex::buildInfoGetBuildGitHash();
        const char* buildgitname = amrex::buildInfoGetBuildGitName();

        if (strlen(erf_hash) > 0) {
          amrex::Print() << "\n"
                         << "ERF git hash: " << erf_hash << "\n";
        }
        if (strlen(amrex_hash) > 0) {
          amrex::Print() << "AMReX git hash: " << amrex_hash << "\n";
        }
        if (strlen(buildgithash) > 0) {
          amrex::Print() << buildgitname << " git hash: " << buildgithash << "\n";
        }

        amrex::Print() << "\n";
    }

    ReadParameters();
    const std::string& pv1 = "plot_vars_1"; setPlotVariables(pv1,plot_var_names_1);
    const std::string& pv2 = "plot_vars_2"; setPlotVariables(pv2,plot_var_names_2);

    amrex_probinit(geom[0].ProbLo(),geom[0].ProbHi());

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= max_level; ++lev) {
        nsubsteps[lev] = MaxRefRatio(lev-1);
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);
    dt_mri_ratio.resize(nlevs_max, 1);

    vars_new.resize(nlevs_max);
    vars_old.resize(nlevs_max);

    for (int lev = 0; lev < nlevs_max; ++lev) {
        vars_new[lev].resize(Vars::NumTypes);
        vars_old[lev].resize(Vars::NumTypes);
    }

    mri_integrator_mem.resize(nlevs_max);
    sri_integrator_mem.resize(nlevs_max);

    flux_registers.resize(nlevs_max);

    // Initialize tagging criteria for mesh refinement
    refinement_criteria_setup();

    // We have already read in the ref_Ratio (via amr.ref_ratio =) but we need to enforce
    //     that there is no refinement in the vertical so we test on that here.
    for (int lev = 0; lev < max_level; ++lev)
    {
       amrex::Print() << "Refinement ratio at level " << lev << " set to be " <<
          ref_ratio[lev][0]  << " " << ref_ratio[lev][1]  <<  " " << ref_ratio[lev][2] << std::endl;

       if (ref_ratio[lev][2] != 1)
       {
           amrex::Error("We don't allow refinement in the vertical -- make sure to set ref_ratio = 1 in z");
       }
    }
}

ERF::~ERF ()
{
}

// advance solution to final time
void
ERF::Evolve ()
{
    Real cur_time = t_new[0];

    // Take one coarse timestep by calling timeStep -- which recursively calls timeStep
    //      for finer levels (with or without subcycling)
    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        amrex::Print() << "\nCoarse STEP " << step+1 << " starts ..." << std::endl;

        ComputeDt();

        // Make sure we have read enough of the boundary plane data to make it through this timestep
        if (input_bndry_planes)
        {
            m_r2d->read_input_files(cur_time,dt[0],m_bc_extdir_vals);
        }

        int lev = 0;
        int iteration = 1;
        timeStep(lev, cur_time, iteration);

        cur_time  += dt[0];

        amrex::Print() << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time
                       << " DT = " << dt[0]  << std::endl;

        if (plot_int_1 > 0 && (step+1) % plot_int_1 == 0) {
            last_plot_file_step_1 = step+1;
            WritePlotFile(1,plot_var_names_1);
        }
        if (plot_int_2 > 0 && (step+1) % plot_int_2 == 0) {
            last_plot_file_step_2 = step+1;
            WritePlotFile(2,plot_var_names_2);
        }

        if (check_int > 0 && (step+1) % check_int == 0) {
            last_check_file_step = step+1;
#ifdef ERF_USE_NETCDF
            if (check_type == "netcdf") {
               WriteNCCheckpointFile();
            }
#endif
            if (check_type == "native") {
               WriteCheckpointFile();
            }
        }

        post_timestep(step, cur_time, dt[0]);

#ifdef AMREX_MEM_PROFILING
        {
            std::ostringstream ss;
            ss << "[STEP " << step+1 << "]";
            MemProfiler::report(ss.str());
        }
#endif

        if (cur_time >= stop_time - 1.e-6*dt[0]) break;
    }

    if (plot_int_1 > 0 && istep[0] > last_plot_file_step_1) {
        WritePlotFile(1,plot_var_names_1);
    }
    if (plot_int_2 > 0 && istep[0] > last_plot_file_step_2) {
        WritePlotFile(2,plot_var_names_2);
    }

    if (check_int > 0 && istep[0] > last_check_file_step) {
#ifdef ERF_USE_NETCDF
        if (check_type == "netcdf") {
           WriteNCCheckpointFile();
        }
#endif
        if (check_type == "native") {
           WriteCheckpointFile();
        }
    }

}

// Called after every coarse timestep
void
ERF::post_timestep (int nstep, Real time, Real dt_lev0)
{
    BL_PROFILE("ERF::post_timestep()");

    if (do_reflux)
    {
        for (int lev = finest_level-1; lev >= 0; lev--)
        {
            // This call refluxes from the lev/lev+1 interface onto lev
            get_flux_reg(lev+1).Reflux(vars_new[lev][Vars::cons],1.0, 0, 0, NVAR, geom[lev]);

            // We need to do this before anything else because refluxing changes the
            // values of coarse cells underneath fine grids with the assumption they'll
            // be over-written by averaging down
            AverageDownTo(lev);
        }
    }

    if (is_it_time_for_action(nstep, time, dt_lev0, sum_interval, sum_per)) {
        sum_integrated_quantities(time);
    }

    if (output_1d_column) {
#ifdef ERF_USE_NETCDF
      if (is_it_time_for_action(nstep, time, dt_lev0, column_interval, column_per))
      {
         int lev_column = 0;
         for (int lev = finest_level; lev >= 0; lev--)
         {
            Real dx_lev = geom[lev].CellSize(0);
            Real dy_lev = geom[lev].CellSize(1);
            int i_lev = static_cast<int>(std::floor(column_loc_x / dx_lev));
            int j_lev = static_cast<int>(std::floor(column_loc_y / dy_lev));
            if (grids[lev].contains(IntVect(i_lev,j_lev,0))) lev_column = lev;
         }
         writeToNCColumnFile(lev_column, column_file_name, column_loc_x, column_loc_y, time);
      }
#else
      amrex::Abort("To output 1D column files ERF must be compiled with NetCDF");
#endif
    }

    if (output_bndry_planes)
    {
      if (is_it_time_for_action(istep[0], time, dt_lev0, bndry_output_planes_interval, bndry_output_planes_per) &&
          time >= bndry_output_planes_start_time)
      {
         m_w2d->write_planes(istep[0], time, vars_new);
      }
    }
}

// This is called from main.cpp and handles all initialization, whether from start or restart
void
ERF::InitData ()
{
    // Initialize the start time for our CPU-time tracker
    startCPUTime = amrex::ParallelDescriptor::second();

    // Map the words in the inputs file to BC types, then translate
    //     those types into what they mean for each variable
    init_bcs();

    // Verify BCs are compatible sith solver choice
    if (solverChoice.pbl_type == PBLType::MYNN25 &&
        phys_bc_type[Orientation(Direction::z,Orientation::low)] != ERF_BC::MOST) {
        amrex::Abort("MYNN2.5 PBL Model requires MOST at lower boundary");
    }

    last_plot_file_step_1 = -1;
    last_plot_file_step_2 = -1;
    last_check_file_step = -1;

    if (restart_chkfile == "") {
        // start simulation from the beginning

        const Real time = 0.0;
        InitFromScratch(time);

#ifdef ERF_USE_MULTIBLOCK
        // Multiblock: hook to set BL & comms once ba/dm are known
        if(domain_p[0].bigEnd(0) < 500 ) {
            m_mbc->SetBoxLists();
            m_mbc->SetBlockCommMetaData();
        }
#endif

        // For now we initialize rho_KE to 0
        for (int lev = finest_level-1; lev >= 0; --lev)
            vars_new[lev][Vars::cons].setVal(0.0,RhoKE_comp,1,0);

        if (init_type == "ideal" && solverChoice.use_terrain) {
            amrex::Abort("We do not currently support init_type = ideal with terrain");
        } else if (init_type == "input_sounding" && solverChoice.use_terrain) {
            amrex::Abort("We do not currently support init_type = input_sounding with terrain");
        }

        if (solverChoice.use_terrain) {
            if (init_type != "real") {
                for (int lev = 0; lev <= finest_level; lev++)
                {
                    init_custom_terrain(geom[lev],*z_phys_nd[lev]);
                    init_terrain_grid(geom[lev],*z_phys_nd[lev]);
                    make_metrics(geom[lev],*z_phys_nd[lev],*z_phys_cc[lev],*detJ_cc[lev]);
                }
            }
        }

        for (int lev = 0; lev <= finest_level; lev++)
            init_only(lev, time);

        AverageDown();

        if(solverChoice.use_terrain) {
            if (init_type == "real") {
                for (int lev = 0; lev <= finest_level; lev++)
                {
                    make_metrics(geom[lev],*z_phys_nd[lev],*z_phys_cc[lev],*detJ_cc[lev]);
                }
            }
        }

    } else { // Restart from a checkpoint

        restart();

        if(solverChoice.use_terrain) {
            // This must come after the call to restart because that
            //      is where we read in the mesh data
            for (int lev = finest_level-1; lev >= 0; --lev)
                make_metrics(geom[lev],*z_phys_nd[lev],*z_phys_cc[lev],*detJ_cc[lev]);
        }
    }

    if (input_bndry_planes) {
        // Create the ReadBndryPlanes object so we can handle reading of boundary plane data
        amrex::Print() << "Defining r2d for the first time " << std::endl;
        m_r2d = std::make_unique< ReadBndryPlanes>(geom[0]);

        // Read the "time.dat" file to know what data is available
        m_r2d->read_time_file();

        amrex::Real dt_dummy = 1.e200;
        m_r2d->read_input_files(t_new[0],dt_dummy,m_bc_extdir_vals);
    }

    // Initialize flux registers (whether we start from scratch or restart)
    if (do_reflux) {
        flux_registers[0] = 0;
        for (int lev = 1; lev <= finest_level; lev++)
        {
            flux_registers[lev] = new FluxRegister(grids[lev], dmap[lev], ref_ratio[lev-1], lev, NVAR);
        }
    }

    initHSE();

    // Configure ABLMost params if used MostWall boundary condition
    // NOTE: we must set up the MOST routine before calling WritePlotFile because
    //       WritePlotFile calls FillPatch in order to compute gradients
    if (phys_bc_type[Orientation(Direction::z,Orientation::low)] == ERF_BC::MOST)
    {
        m_most = std::make_unique<ABLMost>(geom);
        for (int lev = 0; lev <= finest_level; lev++)
            setupABLMost(lev);
    }

    if (restart_chkfile == "" && check_int > 0)
    {
#ifdef ERF_USE_NETCDF
        if (check_type == "netcdf") {
           WriteNCCheckpointFile();
        }
#endif
        if (check_type == "native") {
           WriteCheckpointFile();
        }
        last_check_file_step = 0;
    }

    if ( (restart_chkfile == "") ||
         (restart_chkfile != "" && plot_file_on_restart) )
    {
        if (plot_int_1 > 0)
        {
            WritePlotFile(1,plot_var_names_1);
            last_plot_file_step_1 = istep[0];
        }
        if (plot_int_2 > 0)
        {
            WritePlotFile(2,plot_var_names_2);
            last_plot_file_step_2 = istep[0];
        }
    }

    if (solverChoice.use_rayleigh_damping)
    {
        initRayleigh();
    }

    if (is_it_time_for_action(istep[0], t_new[0], dt[0], sum_interval, sum_per)) {
        sum_integrated_quantities(t_new[0]);
    }

    // We only write the file at level 0 for now
    if (output_bndry_planes)
    {
        // Create the WriteBndryPlanes object so we can handle writing of boundary plane data
        m_w2d = std::make_unique<WriteBndryPlanes>(grids,geom);

        amrex::Real time = 0.;
        if (time >= bndry_output_planes_start_time) {
            m_w2d->write_planes(0, time, vars_new);
        }
    }

    // Fill ghost cells/faces
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        auto& lev_new = vars_new[lev];
        auto& lev_old = vars_old[lev];

        FillPatch(lev, t_new[lev], lev_new);

        // Copy from new into old just in case
        int ngs   = lev_new[Vars::cons].nGrow();
        int ngvel = lev_new[Vars::xvel].nGrow();
        MultiFab::Copy(lev_old[Vars::cons],lev_new[Vars::cons],0,0,NVAR,ngs);
        MultiFab::Copy(lev_old[Vars::xvel],lev_new[Vars::xvel],0,0,1,ngvel);
        MultiFab::Copy(lev_old[Vars::yvel],lev_new[Vars::yvel],0,0,1,ngvel);
        MultiFab::Copy(lev_old[Vars::zvel],lev_new[Vars::zvel],0,0,1,IntVect(ngvel,ngvel,0));
    }
}

void
ERF::restart()
{
#ifdef ERF_USE_NETCDF
    if (restart_type == "netcdf") {
       ReadNCCheckpointFile();
    }
#endif
    if (restart_type == "native") {
       ReadCheckpointFile();
    }

    // We set this here so that we don't over-write the checkpoint file we just started from
    last_check_file_step = istep[0];
}

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void
ERF::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
                                    const DistributionMapping& dm)
{
    const auto& crse_new = vars_new[lev-1];
    auto& lev_new = vars_new[lev];
    auto& lev_old = vars_old[lev];

    lev_new[Vars::cons].define(ba, dm, crse_new[Vars::cons].nComp(), crse_new[Vars::cons].nGrowVect());
    lev_old[Vars::cons].define(ba, dm, crse_new[Vars::cons].nComp(), crse_new[Vars::cons].nGrowVect());

    lev_new[Vars::xvel].define(convert(ba, IntVect(1,0,0)), dm, 1, crse_new[Vars::xvel].nGrowVect());
    lev_old[Vars::xvel].define(convert(ba, IntVect(1,0,0)), dm, 1, crse_new[Vars::xvel].nGrowVect());

    lev_new[Vars::yvel].define(convert(ba, IntVect(0,1,0)), dm, 1, crse_new[Vars::yvel].nGrowVect());
    lev_old[Vars::yvel].define(convert(ba, IntVect(0,1,0)), dm, 1, crse_new[Vars::yvel].nGrowVect());

    lev_new[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, crse_new[Vars::zvel].nGrowVect());
    lev_old[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, crse_new[Vars::zvel].nGrowVect());

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    FillCoarsePatchAllVars(lev, time, vars_new[lev]);

    // Initialize the integrator memory
    amrex::Vector<amrex::MultiFab> int_state; // integration state data structure example
    int_state.push_back(MultiFab(lev_new[Vars::cons], amrex::make_alias, 0, Cons::NumVars)); // cons
    int_state.push_back(MultiFab(convert(ba,IntVect(1,0,0)), dm, 1, lev_new[Vars::xvel].nGrow())); // xmom
    int_state.push_back(MultiFab(convert(ba,IntVect(0,1,0)), dm, 1, lev_new[Vars::yvel].nGrow())); // ymom
    int_state.push_back(MultiFab(convert(ba,IntVect(0,0,1)), dm, 1, lev_new[Vars::zvel].nGrow())); // zmom
    int_state.push_back(MultiFab(convert(ba,IntVect(1,0,0)), dm, Cons::NumVars, 1)); // x-fluxes
    int_state.push_back(MultiFab(convert(ba,IntVect(0,1,0)), dm, Cons::NumVars, 1)); // y-fluxes
    int_state.push_back(MultiFab(convert(ba,IntVect(0,0,1)), dm, Cons::NumVars, 1)); // z-fluxes
    if (use_native_mri) {
        mri_integrator_mem[lev] = std::make_unique<MRISplitIntegrator<amrex::Vector<amrex::MultiFab> > >(int_state);
    } else {
        sri_integrator_mem[lev] = std::make_unique<SRIIntegrator<amrex::Vector<amrex::MultiFab> > >(int_state);
    }
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
ERF::RemakeLevel (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm)
{
    Vector<MultiFab> temp_lev_new(Vars::NumTypes);
    Vector<MultiFab> temp_lev_old(Vars::NumTypes);

    int ngrow_state = ComputeGhostCells(solverChoice.spatial_order)+1;
    int ngrow_vels  = ComputeGhostCells(solverChoice.spatial_order);

    temp_lev_new[Vars::cons].define(ba, dm, Cons::NumVars, ngrow_state);
    temp_lev_old[Vars::cons].define(ba, dm, Cons::NumVars, ngrow_state);

    temp_lev_new[Vars::xvel].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);
    temp_lev_old[Vars::xvel].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);

    temp_lev_new[Vars::yvel].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);
    temp_lev_old[Vars::yvel].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);

    temp_lev_new[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, IntVect(ngrow_vels,ngrow_vels,0));
    temp_lev_old[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, IntVect(ngrow_vels,ngrow_vels,0));

    // This will fill the temporary MultiFabs with data from vars_new
    FillPatch(lev, time, temp_lev_new);
    FillPatch(lev, time, temp_lev_old);

    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx) {
        std::swap(temp_lev_new[var_idx], vars_new[lev][var_idx]);
        std::swap(temp_lev_old[var_idx], vars_old[lev][var_idx]);
    }

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    // Initialize the integrator memory
    amrex::Vector<amrex::MultiFab> int_state; // integration state data structure example
    int_state.push_back(MultiFab(temp_lev_new[Vars::cons], amrex::make_alias, 0, Cons::NumVars)); // cons
    int_state.push_back(MultiFab(convert(ba,IntVect(1,0,0)), dm, 1, temp_lev_new[Vars::xvel].nGrow())); // xmom
    int_state.push_back(MultiFab(convert(ba,IntVect(0,1,0)), dm, 1, temp_lev_new[Vars::yvel].nGrow())); // ymom
    int_state.push_back(MultiFab(convert(ba,IntVect(0,0,1)), dm, 1, temp_lev_new[Vars::zvel].nGrow())); // zmom
    int_state.push_back(MultiFab(convert(ba,IntVect(1,0,0)), dm, Cons::NumVars, 1)); // x-fluxes
    int_state.push_back(MultiFab(convert(ba,IntVect(0,1,0)), dm, Cons::NumVars, 1)); // y-fluxes
    int_state.push_back(MultiFab(convert(ba,IntVect(0,0,1)), dm, Cons::NumVars, 1)); // z-fluxes
    if (use_native_mri) {
        mri_integrator_mem[lev] = std::make_unique<MRISplitIntegrator<amrex::Vector<amrex::MultiFab> > >(int_state);
    } else {
        sri_integrator_mem[lev] = std::make_unique<SRIIntegrator<amrex::Vector<amrex::MultiFab> > >(int_state);
    }
}

// Delete level data
// overrides the pure virtual function in AmrCore
void
ERF::ClearLevel (int lev)
{
    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx) {
        vars_new[lev][var_idx].clear();
        vars_old[lev][var_idx].clear();
    }

    // Clears the integrator memory
    if (use_native_mri) {
        mri_integrator_mem[lev].reset();
    } else {
        sri_integrator_mem[lev].reset();
    }
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// This is called both for initialization and for restart
// (overrides the pure virtual function in AmrCore)
// main.cpp --> ERF::InitData --> InitFromScratch --> MakeNewGrids --> MakeNewLevelFromScratch
//                                       restart  --> MakeNewGrids --> MakeNewLevelFromScratch
void ERF::MakeNewLevelFromScratch (int lev, Real /*time*/, const BoxArray& ba,
                                   const DistributionMapping& dm)
{
    // Set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);

    // The number of ghost cells for density must be 1 greater than that for velocity
    //     so that we can go back in forth betwen velocity and momentum on all faces
    int ngrow_state = ComputeGhostCells(solverChoice.spatial_order)+1;
    int ngrow_vels  = ComputeGhostCells(solverChoice.spatial_order);

    auto& lev_new = vars_new[lev];
    auto& lev_old = vars_old[lev];

    lev_new[Vars::cons].define(ba, dm, Cons::NumVars, ngrow_state);
    lev_old[Vars::cons].define(ba, dm, Cons::NumVars, ngrow_state);

    lev_new[Vars::xvel].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);
    lev_old[Vars::xvel].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);

    lev_new[Vars::yvel].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);
    lev_old[Vars::yvel].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);

    lev_new[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, IntVect(ngrow_vels,ngrow_vels,0));
    lev_old[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, IntVect(ngrow_vels,ngrow_vels,0));

    z_phys_nd.resize(lev+1);
    z_phys_cc.resize(lev+1);
    detJ_cc.resize(lev+1);

    BoxList bl2d = ba.boxList();
    for (auto& b : bl2d) {
        b.setRange(2,0);
    }
    BoxArray ba2d(std::move(bl2d));

    mapfac_m.resize(lev+1);
    mapfac_u.resize(lev+1);
    mapfac_v.resize(lev+1);
    mapfac_m[lev].reset(new MultiFab(ba2d,dm,1,0));
    mapfac_u[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,0));
    mapfac_v[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,0));
    mapfac_m[lev]->setVal(1.);
    mapfac_u[lev]->setVal(1.);
    mapfac_v[lev]->setVal(1.);

    pres_hse.resize(lev+1);
    dens_hse.resize(lev+1);
    pres_hse[lev].define(ba,dm,1,1);
    dens_hse[lev].define(ba,dm,1,1);
    pres_hse[lev].setVal(0.);
    dens_hse[lev].setVal(0.);

    if (solverChoice.use_terrain) {
        z_phys_cc[lev].reset(new MultiFab(ba,dm,1,1));
          detJ_cc[lev].reset(new MultiFab(ba,dm,1,1));

        BoxArray ba_nd(ba);
        ba_nd.surroundingNodes();

        // We need this to be one greater than the ghost cells to handle levels > 0
        int ngrow = ComputeGhostCells(solverChoice.spatial_order)+2;
        z_phys_nd[lev].reset(new MultiFab(ba_nd,dm,1,IntVect(ngrow,ngrow,1)));
    } else {
        z_phys_nd[lev] = nullptr;
        z_phys_cc[lev] = nullptr;
          detJ_cc[lev] = nullptr;
    }

    // Initialize the integrator memory
    amrex::Vector<amrex::MultiFab> int_state; // integration state data structure example
    int_state.push_back(MultiFab(lev_new[Vars::cons], amrex::make_alias, 0, Cons::NumVars)); // cons
    int_state.push_back(MultiFab(convert(ba,IntVect(1,0,0)), dm, 1, lev_new[Vars::xvel].nGrow())); // xmom
    int_state.push_back(MultiFab(convert(ba,IntVect(0,1,0)), dm, 1, lev_new[Vars::yvel].nGrow())); // ymom
    int_state.push_back(MultiFab(convert(ba,IntVect(0,0,1)), dm, 1, lev_new[Vars::zvel].nGrow())); // zmom
    int_state.push_back(MultiFab(convert(ba,IntVect(1,0,0)), dm, Cons::NumVars, 1)); // x-fluxes
    int_state.push_back(MultiFab(convert(ba,IntVect(0,1,0)), dm, Cons::NumVars, 1)); // y-fluxes
    int_state.push_back(MultiFab(convert(ba,IntVect(0,0,1)), dm, Cons::NumVars, 1)); // z-fluxes
    if (use_native_mri) {
        mri_integrator_mem[lev] = std::make_unique<MRISplitIntegrator<amrex::Vector<amrex::MultiFab> > >(int_state);
    } else {
        sri_integrator_mem[lev] = std::make_unique<SRIIntegrator<amrex::Vector<amrex::MultiFab> > >(int_state);
    }
}

void
ERF::init_only(int lev, Real time)
{
    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    auto& lev_new = vars_new[lev];

    // Loop over grids at this level to initialize our grid data
    lev_new[Vars::cons].setVal(0.0);
    lev_new[Vars::xvel].setVal(0.0);
    lev_new[Vars::yvel].setVal(0.0);
    lev_new[Vars::zvel].setVal(0.0);

    if (init_type == "input_sounding") {
        init_from_input_sounding(lev);
    } else if (init_type == "custom") {
        init_custom(lev);
#ifdef ERF_USE_NETCDF
    } else if (init_type == "ideal" || init_type == "real") {
        init_from_wrfinput(lev);
#endif
    }

    // Ensure that the face-based data are the same on both sides of a periodic domain.
    // The data associated with the lower grid ID is considered the correct value.
    lev_new[Vars::xvel].OverrideSync(geom[lev].periodicity());
    lev_new[Vars::yvel].OverrideSync(geom[lev].periodicity());
    lev_new[Vars::zvel].OverrideSync(geom[lev].periodicity());
}

// read in some parameters from inputs file
void
ERF::ReadParameters ()
{
    {
        ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);
    }

    ParmParse pp(pp_prefix);
    {
        // The type of the file we restart from
        pp.query("restart_type", restart_type);

        pp.query("regrid_int", regrid_int);
        pp.query("check_file", check_file);
        pp.query("check_type", check_type);
        pp.query("check_int", check_int);

        pp.query("restart", restart_chkfile);

        if (pp.contains("data_log"))
        {
            int num_datalogs = pp.countval("data_log");
            datalog.resize(num_datalogs);
            datalogname.resize(num_datalogs);
            pp.queryarr("data_log",datalogname,0,num_datalogs);
            for (int i = 0; i < num_datalogs; i++)
                setRecordDataInfo(i,datalogname[i]);
        }

        // Verbosity
        pp.query("v", verbose);

        // Use the native ERF MRI integrator
        pp.query("use_native_mri", use_native_mri);

        // Frequency of diagnostic output
        pp.query("sum_interval", sum_interval);
        pp.query("sum_period"  , sum_per);

        // Time step controls
        pp.query("cfl", cfl);
        pp.query("init_shrink", init_shrink);
        pp.query("change_max", change_max);

        pp.query("fixed_dt", fixed_dt);
        pp.query("fixed_fast_dt", fixed_fast_dt);
        pp.query("fixed_mri_dt_ratio", fixed_mri_dt_ratio);
        pp.query("use_lowM_dt", use_lowM_dt);

        if (fixed_dt > 0. && fixed_fast_dt > 0.) {
            if (fixed_mri_dt_ratio > 0 &&
               fixed_dt / fixed_fast_dt != fixed_mri_dt_ratio)
            {
                amrex::Abort("Dt is over-specfied");
            } else {
                fixed_mri_dt_ratio = static_cast<int>(fixed_dt / fixed_fast_dt);
                if (fixed_mri_dt_ratio%2 != 0)
                    fixed_mri_dt_ratio += 1; // This ratio must be even
            }
        }

        AMREX_ASSERT(cfl > 0. || fixed_dt > 0.);

        // Mesh refinement
        pp.query("coupling_type",coupling_type);
        if (coupling_type == "OneWay")
        {
            do_reflux = 0;
            do_avg_down = 0;
        }
        else if (coupling_type == "TwoWay")
        {
            do_reflux = 1;
            do_avg_down = 1;
        } else {
            amrex::Error("Unknown coupling type");
        }

        // How to initialize
        pp.query("init_type",init_type);
        if (init_type != "custom" &&
            init_type != "ideal" &&
            init_type != "real" &&
            init_type != "input_sounding")
        {
            amrex::Error("init_type must be custom, ideal, real, or input_sounding");
        }

        // We use this to keep track of how many boxes we read in from WRF initialization
        num_files_at_level.resize(max_level+1,0);

        // We use this to keep track of how many boxes are specified thru the refinement indicators
        num_boxes_at_level.resize(max_level+1,0);
            boxes_at_level.resize(max_level+1);

        // We always have exactly one file at level 0
        num_boxes_at_level[0] = 1;

#ifdef ERF_USE_NETCDF
        nc_init_file.resize(max_level+1);

        // NetCDF wrfinput initialization files -- possibly multiple files at each of multiple levels
        //        but we always have exactly one file at level 0
        for (int lev = 0; lev <= max_level; lev++)
        {
            const std::string nc_file_names = amrex::Concatenate("nc_init_file_",lev,1);
            if (pp.contains(nc_file_names.c_str()))
            {
                int num_files = pp.countval(nc_file_names.c_str());
                num_files_at_level[lev] = num_files;
                nc_init_file[lev].resize(num_files);
                pp.queryarr(nc_file_names.c_str(), nc_init_file[lev],0,num_files);
                for (int j = 0; j < num_files; j++)
                    amrex::Print() << "Reading NC init file names at level " << lev << " and index " << j << " : " << nc_init_file[lev][j] << std::endl;
            }
        }

        // NetCDF wrfbdy lateral boundary file
        pp.query("nc_bdy_file", nc_bdy_file);
#endif

        // Text input_sounding file
        pp.query("input_sounding_file", input_sounding_file);

        // Output format
        pp.query("plotfile_type", plotfile_type);
        if (plotfile_type != "amrex" &&
            plotfile_type != "netcdf" && plotfile_type != "NetCDF" &&
            plotfile_type != "hdf5"   && plotfile_type != "HDF5" )
        {
            amrex::Print() << "User selected plotfile_type = " << plotfile_type << std::endl;
            amrex::Abort("Dont know this plotfile_type");
        }
        pp.query("plot_file_1", plot_file_1);
        pp.query("plot_file_2", plot_file_2);
        pp.query("plot_int_1", plot_int_1);
        pp.query("plot_int_2", plot_int_2);

        pp.query("output_1d_column", output_1d_column);
        pp.query("column_per", column_per);
        pp.query("column_interval", column_interval);
        pp.query("column_loc_x", column_loc_x);
        pp.query("column_loc_y", column_loc_y);
        pp.query("column_file_name", column_file_name);

        // Specify information about outputting planes of data
        pp.query("output_bndry_planes", output_bndry_planes);
        pp.query("bndry_output_planes_interval", bndry_output_planes_interval);
        pp.query("bndry_output_planes_per", bndry_output_planes_per);
        pp.query("bndry_output_start_time", bndry_output_planes_start_time);

        // Specify whether ingest boundary planes of data
        pp.query("input_bndry_planes", input_bndry_planes);
    }

#ifdef ERF_USE_MULTIBLOCK
    solverChoice.pp_prefix = pp_prefix;
#endif

    solverChoice.init_params();
}

// Create horizontal average quantities
void
ERF::MakeHorizontalAverages ()
{
    // We need to create horizontal averages for 5 variables:
    // density, temperature, pressure, qc, qv (if present)

    // First, average down all levels
    AverageDown();

    // get the number of cells in z at level 0
    int dir_z = AMREX_SPACEDIM-1;
    auto domain = geom[0].Domain();
    int size_z = domain.length(dir_z);
    int start_z = domain.smallEnd()[dir_z];
    Real area_z = static_cast<Real>(domain.length(0));
    area_z *= domain.length(1);

    // resize the level 0 horizontal average vectors
    h_havg_density.resize(size_z, 0.0_rt);
    h_havg_temperature.resize(size_z, 0.0_rt);
    h_havg_pressure.resize(size_z, 0.0_rt);
#ifdef ERF_USE_MOISTURE
    h_havg_qv.resize(size_z, 0.0_rt);
    h_havg_qc.resize(size_z, 0.0_rt);
#endif

    // get the cell centered data and construct sums
    auto& mf_cons = vars_new[0][Vars::cons];
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf_cons); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const IntVect& se = box.smallEnd();
        const IntVect& be = box.bigEnd();
        auto  arr_cons = mf_cons[mfi].array();

#ifdef ERF_USE_MOISTURE
        FArrayBox fab_reduce(box, 5);
#else
        FArrayBox fab_reduce(box, 3);
#endif
        Elixir elx_reduce = fab_reduce.elixir();
        auto arr_reduce = fab_reduce.array();

        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Real dens = arr_cons(i, j, k, Cons::Rho);
            arr_reduce(i, j, k, 0) = dens;
            arr_reduce(i, j, k, 1) = arr_cons(i, j, k, Cons::RhoTheta) / dens;
            arr_reduce(i, j, k, 2) = getPgivenRTh(arr_cons(i, j, k, Cons::RhoTheta));
#ifdef ERF_USE_MOISTURE
            arr_reduce(i, j, k, 3) = arr_cons(i, j, k, Cons::RhoQv) / dens;
            arr_reduce(i, j, k, 4) = arr_cons(i, j, k, Cons::RhoQc) / dens;
#endif
        });

        for (int k=se[dir_z]; k <= be[dir_z]; ++k) {
            Box kbox(box); kbox.setSmall(dir_z,k); kbox.setBig(dir_z,k);
            h_havg_density     [k-start_z] += fab_reduce.sum<RunOn::Device>(kbox,0);
            h_havg_temperature [k-start_z] += fab_reduce.sum<RunOn::Device>(kbox,1);
            h_havg_pressure    [k-start_z] += fab_reduce.sum<RunOn::Device>(kbox,2);
#ifdef ERF_USE_MOISTURE
            h_havg_qv          [k-start_z] += fab_reduce.sum<RunOn::Device>(kbox,3);
            h_havg_qc          [k-start_z] += fab_reduce.sum<RunOn::Device>(kbox,4);
#endif
        }
    }

    // combine sums from different MPI ranks
    ParallelDescriptor::ReduceRealSum(h_havg_density.dataPtr(), h_havg_density.size());
    ParallelDescriptor::ReduceRealSum(h_havg_temperature.dataPtr(), h_havg_temperature.size());
    ParallelDescriptor::ReduceRealSum(h_havg_pressure.dataPtr(), h_havg_pressure.size());
#ifdef ERF_USE_MOISTURE
    ParallelDescriptor::ReduceRealSum(h_havg_qv.dataPtr(), h_havg_qv.size());
    ParallelDescriptor::ReduceRealSum(h_havg_qc.dataPtr(), h_havg_qc.size());
#endif

    // divide by the total number of cells we are averaging over
    for (int k = 0; k < size_z; ++k) {
        h_havg_density[k]     /= area_z;
        h_havg_temperature[k] /= area_z;
        h_havg_pressure[k]    /= area_z;
#ifdef ERF_USE_MOISTURE
        h_havg_qv[k]          /= area_z;
        h_havg_qc[k]          /= area_z;
#endif
    }

    // resize device vectors
    d_havg_density.resize(size_z, 0.0_rt);
    d_havg_temperature.resize(size_z, 0.0_rt);
    d_havg_pressure.resize(size_z, 0.0_rt);
#ifdef ERF_USE_MOISTURE
    d_havg_qv.resize(size_z, 0.0_rt);
    d_havg_qc.resize(size_z, 0.0_rt);
#endif

    // copy host vectors to device vectors
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_density.begin(), h_havg_density.end(), d_havg_density.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_temperature.begin(), h_havg_temperature.end(), d_havg_temperature.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_pressure.begin(), h_havg_pressure.end(), d_havg_pressure.begin());
#ifdef ERF_USE_MOISTURE
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_qv.begin(), h_havg_qv.end(), d_havg_qv.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_qc.begin(), h_havg_qc.end(), d_havg_qc.begin());
#endif
}

// Set covered coarse cells to be the average of overlying fine cells for all levels
void
ERF::AverageDown ()
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
        AverageDownTo(lev);
    }
}

// Set covered coarse cells to be the average of overlying fine cells at level crse_lev
void
ERF::AverageDownTo (int crse_lev)
{
    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx) {
        const BoxArray& ba(vars_new[crse_lev][var_idx].boxArray());
        if (ba[0].type() == IntVect::TheZeroVector())
            amrex::average_down(vars_new[crse_lev+1][var_idx], vars_new[crse_lev][var_idx],
                                0, vars_new[crse_lev][var_idx].nComp(), refRatio(crse_lev));
        else // We assume the arrays are face-centered if not cell-centered
            amrex::average_down_faces(vars_new[crse_lev+1][var_idx], vars_new[crse_lev][var_idx],
                                      refRatio(crse_lev),geom[crse_lev]);
    }
}

void
ERF::setupABLMost (int lev)
{
    MultiFab& S_new = vars_new[lev][Vars::cons];
    MultiFab& U_new = vars_new[lev][Vars::xvel];
    MultiFab& V_new = vars_new[lev][Vars::yvel];
    MultiFab& W_new = vars_new[lev][Vars::zvel];

    // Multifab to store primitive Theta, which is what we want to average
    MultiFab Theta_prim(S_new.boxArray(), S_new.DistributionMap(), 1, 0);
    MultiFab::Copy(Theta_prim, S_new, Cons::RhoTheta, 0, 1, 0);
    MultiFab::Divide(Theta_prim, S_new, Cons::Rho, 0, 1, 0);

    m_most->update_fluxes(lev,Theta_prim,U_new,V_new,W_new);
}

#ifdef ERF_USE_MULTIBLOCK
// constructor used when ERF is created by a multiblock driver
ERF::ERF (const amrex::RealBox& rb, int max_level_in,
          const amrex::Vector<int>& n_cell_in, int coord,
          const amrex::Vector<amrex::IntVect>& ref_ratios,
          const amrex::Array<int,AMREX_SPACEDIM>& is_per,
          std::string prefix)
    : AmrCore(rb, max_level_in, n_cell_in, coord, ref_ratios, is_per)
{
    SetParmParsePrefix(prefix);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        const char* erf_hash = amrex::buildInfoGetGitHash(1);
        const char* amrex_hash = amrex::buildInfoGetGitHash(2);
        const char* buildgithash = amrex::buildInfoGetBuildGitHash();
        const char* buildgitname = amrex::buildInfoGetBuildGitName();

        if (strlen(erf_hash) > 0) {
          amrex::Print() << "\n"
                         << "ERF git hash: " << erf_hash << "\n";
        }
        if (strlen(amrex_hash) > 0) {
          amrex::Print() << "AMReX git hash: " << amrex_hash << "\n";
        }
        if (strlen(buildgithash) > 0) {
          amrex::Print() << buildgitname << " git hash: " << buildgithash << "\n";
        }

        amrex::Print() << "\n";
    }

    ReadParameters();
    const std::string& pv1 = "plot_vars_1"; setPlotVariables(pv1,plot_var_names_1);
    const std::string& pv2 = "plot_vars_2"; setPlotVariables(pv2,plot_var_names_2);

    amrex_probinit(geom[0].ProbLo(),geom[0].ProbHi());

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= max_level; ++lev) {
        nsubsteps[lev] = MaxRefRatio(lev-1);
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);
    dt_mri_ratio.resize(nlevs_max, 1);

    vars_new.resize(nlevs_max);
    vars_old.resize(nlevs_max);

    for (int lev = 0; lev < nlevs_max; ++lev) {
        vars_new[lev].resize(Vars::NumTypes);
        vars_old[lev].resize(Vars::NumTypes);
    }

    mri_integrator_mem.resize(nlevs_max);
    sri_integrator_mem.resize(nlevs_max);

    // Multiblock: public domain sizes (need to know which vars are nodal)
    Box nbx;
    domain_p.push_back(geom[0].Domain());
    nbx = convert(domain_p[0],IntVect(1,0,0));
    domain_p.push_back(nbx);
    nbx = convert(domain_p[0],IntVect(0,1,0));
    domain_p.push_back(nbx);
    nbx = convert(domain_p[0],IntVect(0,0,1));
    domain_p.push_back(nbx);

    flux_registers.resize(nlevs_max);

    // Initialize tagging criteria for mesh refinement
    refinement_criteria_setup();

    // We have already read in the ref_Ratio (via amr.ref_ratio =) but we need to enforce
    //     that there is no refinement in the vertical so we test on that here.
    for (int lev = 0; lev < max_level; ++lev)
    {
       amrex::Print() << "Refinement ratio at level " << lev << " set to be " <<
          ref_ratio[lev][0]  << " " << ref_ratio[lev][1]  <<  " " << ref_ratio[lev][2] << std::endl;

       if (ref_ratio[lev][2] != 1)
       {
           amrex::Error("We don't allow refinement in the vertical -- make sure to set ref_ratio = 1 in z");
       }
    }
}

// advance solution over specified block steps
void
ERF::Evolve_MB (int MBstep, int max_block_step)
{
    Real cur_time = t_new[0];

    int step;

    // Take one coarse timestep by calling timeStep -- which recursively calls timeStep
    // for finer levels (with or without subcycling)
    for (int Bstep(0); Bstep < max_block_step && cur_time < stop_time; ++Bstep)
    {
        step = Bstep + MBstep - 1;

        amrex::Print() << "\nCoarse STEP " << step+1 << " starts ..." << std::endl;

        ComputeDt();

        // Make sure we have read enough of the boundary plane data to make it through this timestep
        if (input_bndry_planes)
        {
            m_r2d->read_input_files(cur_time,dt[0],m_bc_extdir_vals);
        }

        int lev = 0;
        int iteration = 1;
        timeStep(lev, cur_time, iteration);


        // DEBUG
        // Multiblock: hook for erf2 to fill from erf1
        if(domain_p[0].bigEnd(0) < 500) {
            for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx)
                m_mbc->FillPatchBlocks(var_idx,var_idx);
        }


        cur_time  += dt[0];

         amrex::Print() << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time
                       << " DT = " << dt[0]  << std::endl;

        if (plot_int_1 > 0 && (step+1) % plot_int_1 == 0) {
            last_plot_file_step_1 = step+1;
            WritePlotFile(1,plot_var_names_1);
        }
        if (plot_int_2 > 0 && (step+1) % plot_int_2 == 0) {
            last_plot_file_step_2 = step+1;
            WritePlotFile(2,plot_var_names_2);
        }

        if (check_int > 0 && (step+1) % check_int == 0) {
            last_check_file_step = step+1;
#ifdef ERF_USE_NETCDF
            if (check_type == "netcdf") {
               WriteNCCheckpointFile();
            }
#endif
            if (check_type == "native") {
               WriteCheckpointFile();
            }
        }

        post_timestep(step, cur_time, dt[0]);

#ifdef AMREX_MEM_PROFILING
        {
            std::ostringstream ss;
            ss << "[STEP " << step+1 << "]";
            MemProfiler::report(ss.str());
        }
#endif

        if (cur_time >= stop_time - 1.e-6*dt[0]) break;
    }

    if (plot_int_1 > 0 && istep[0] > last_plot_file_step_1) {
        WritePlotFile(1,plot_var_names_1);
    }
    if (plot_int_2 > 0 && istep[0] > last_plot_file_step_2) {
        WritePlotFile(2,plot_var_names_2);
    }

    if (check_int > 0 && istep[0] > last_check_file_step) {
#ifdef ERF_USE_NETCDF
        if (check_type == "netcdf") {
           WriteNCCheckpointFile();
        }
#endif
        if (check_type == "native") {
           WriteCheckpointFile();
        }
    }

}
#endif
