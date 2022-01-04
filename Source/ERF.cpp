/**
 * \file ERF.cpp
 */

#include <ERF.H>
#include <prob_common.H>
#include <AMReX_buildInfo.H>

using namespace amrex;

bool ERF::lo_z_is_dirichlet = false;
bool ERF::hi_z_is_dirichlet = false;

amrex::Real ERF::startCPUTime        = 0.0;
amrex::Real ERF::previousCPUTimeUsed = 0.0;

Vector<AMRErrorTag> ERF::ref_tags;

SolverChoice ERF::solverChoice;
ABLFieldInit ERF::ablinit;
bool         ERF::init_abl = false;

// Create dens_hse and pres_hse with one ghost cell
int ERF::ng_dens_hse = 1;
int ERF::ng_pres_hse = 1;

// Time step control
amrex::Real ERF::cfl         =  0.8;
amrex::Real ERF::fixed_dt    = -1.0;
amrex::Real ERF::init_shrink =  1.0;
amrex::Real ERF::change_max  =  1.1;

// Type of mesh refinement algorithm
std::string ERF::coupling_type = "OneWay";
int         ERF::do_reflux     = 0;
int         ERF::do_avg_down   = 0;

// Dictate verbosity in screen output
int         ERF::verbose       = 0;

// Frequency of diagnostic output
int         ERF::sum_interval  = -1;
amrex::Real ERF::sum_per       = -1.0;

// Native AMReX vs NetCDF
std::string ERF::plotfile_type   = "amrex";

// 1D NetCDF output (for ingestion by AMR-Wind)
int         ERF::output_1d_column = 0;
int         ERF::column_interval  = -1;
amrex::Real ERF::column_per       = -1.0;
amrex::Real ERF::column_loc_x     = 0.0;
amrex::Real ERF::column_loc_y     = 0.0;
std::string ERF::column_file_name = "column_data.nc";

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
    setPlotVariables();

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

    vars_new.resize(nlevs_max);
    vars_old.resize(nlevs_max);

    for (int lev = 0; lev < nlevs_max; ++lev) {
        vars_new[lev].resize(Vars::NumTypes);
        vars_old[lev].resize(Vars::NumTypes);
    }

    integrator.resize(nlevs_max);

    flux_registers.resize(nlevs_max);

    // Map the words in the inputs file to BC types, then translate
    //     those types into what they mean for each variable
    init_bcs();

    // For now we hard-wire interpolation type between levels
    mapper = &cell_cons_interp;

    // Initialize tagging criteria for mesh refinement
    refinement_criteria_setup();
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

        int lev = 0;
        int iteration = 1;
        timeStep(lev, cur_time, iteration);

        cur_time += dt[0];

        amrex::Print() << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time
                       << " DT = " << dt[0]  << std::endl;

        // sync up time
        for (lev = 0; lev <= finest_level; ++lev) {
            t_new[lev] = cur_time;
        }

        if (plot_int > 0 && (step+1) % plot_int == 0) {
            last_plot_file_step = step+1;
            WritePlotFile();
        }

        if (check_int > 0 && (step+1) % check_int == 0) {
            last_check_file_step = step+1;
            WriteCheckpointFile();
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

    if (plot_int > 0 && istep[0] > last_plot_file_step) {
        WritePlotFile();
    }

    if (check_int > 0 && istep[0] > last_check_file_step) {
        WriteCheckpointFile();
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
}

// initializes multilevel data
void
ERF::InitData ()
{
    // Initialize the start time for our CPU-time tracker
    startCPUTime = amrex::ParallelDescriptor::second();

    // For now we initialize rho_KE to 0
    for (int lev = finest_level-1; lev >= 0; --lev)
        vars_new[lev][Vars::cons].setVal(0.0,RhoKE_comp,1,0);

    // This defaults to false; is only true if we want to call ABLFieldInit::init_params
    ParmParse pp("erf");
    pp.query("init_abl", init_abl);

    if (init_abl)
    {
        ablinit.init_params();
    }

    last_plot_file_step = -1;
    last_check_file_step = -1;

    if (restart_chkfile == "") {
        // start simulation from the beginning
        const Real time = 0.0;
        InitFromScratch(time);
        AverageDown();

        for (int lev = finest_level-1; lev >= 0; --lev)
            vars_new[lev][Vars::cons].setVal(0.0,RhoKE_comp,1,0);

        if (check_int > 0)
        {
            WriteCheckpointFile();
            last_check_file_step = 0;
        }

    } else { // Restart from a checkpoint

        ReadCheckpointFile();

        // We set this here so that we don't over-write the checkpoint file we just started from
        last_check_file_step = istep[0];

        // Create the time integrator for this level
        for (int lev = 0; lev <= finest_level; lev++) {
            integrator[lev] = std::make_unique<TimeIntegrator<Vector<MultiFab> > >(vars_old[lev]);
        }
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

    if (plot_int > 0)
    {
        if ( (restart_chkfile == "") ||
             (restart_chkfile != "" && plot_file_on_restart) )
        {
            WritePlotFile();
            last_plot_file_step = istep[0];
        }
    }

    if (solverChoice.use_rayleigh_damping)
    {
        initRayleigh();
    }

    if (is_it_time_for_action(istep[0], t_new[0], dt[0], sum_interval, sum_per)) {
        sum_integrated_quantities(t_new[0]);
    }
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

    // also create the time integrator for this level
    integrator[lev] = std::make_unique<TimeIntegrator<Vector<MultiFab> > >(vars_old[lev]);
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
ERF::RemakeLevel (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm)
{
    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx) {
        const     int ncomp  = vars_new[lev][var_idx].nComp();
        const IntVect nghost = vars_new[lev][var_idx].nGrowVect();

        MultiFab new_v(ba, dm, ncomp, nghost);
        MultiFab old_v(ba, dm, ncomp, nghost);

        FillPatch(lev, time, new_v, 0, ncomp, var_idx);

        std::swap(new_v, vars_new[lev][var_idx]);
        std::swap(old_v, vars_old[lev][var_idx]);
    }

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    // also recreate the time integrator for this level
    integrator[lev] = std::make_unique<TimeIntegrator<Vector<MultiFab> > >(vars_old[lev]);
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
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void ERF::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
                                          const DistributionMapping& dm)
{
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

    lev_new[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, ngrow_vels);
    lev_old[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, ngrow_vels);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    // Loop over grids at this level to initialize our grid data
#ifdef _OPENMP
#pragma omp parallel
#endif
    lev_new[Vars::cons].setVal(1.e20);
    lev_new[Vars::xvel].setVal(1.e20);
    lev_new[Vars::yvel].setVal(1.e20);
    lev_new[Vars::zvel].setVal(1.e20);

    for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      const Box& bx = mfi.tilebox();
      const auto& cons_arr = lev_new[Vars::cons].array(mfi);
      const auto& xvel_arr = lev_new[Vars::xvel].array(mfi);
      const auto& yvel_arr = lev_new[Vars::yvel].array(mfi);
      const auto& zvel_arr = lev_new[Vars::zvel].array(mfi);

      erf_init_prob(bx, cons_arr, xvel_arr, yvel_arr, zvel_arr, geom[lev].data());
    }

    // Ensure that the face-based data are the same on both sides of a periodic domain.
    // The data associated with the lower grid ID is considered the correct value.
    lev_new[Vars::xvel].OverrideSync(geom[lev].periodicity());
    lev_new[Vars::yvel].OverrideSync(geom[lev].periodicity());
    lev_new[Vars::zvel].OverrideSync(geom[lev].periodicity());

    FillIntermediatePatch(lev, time, lev_new[Vars::cons], 0, Cons::NumVars, Vars::cons);
    FillIntermediatePatch(lev, time, lev_new[Vars::xvel], 0, 1, Vars::xvel);
    FillIntermediatePatch(lev, time, lev_new[Vars::yvel], 0, 1, Vars::yvel);
    FillIntermediatePatch(lev, time, lev_new[Vars::zvel], 0, 1, Vars::zvel);

    // Copy from new into old just in case
    int ngs   = lev_new[Vars::cons].nGrow();
    int ngvel = lev_new[Vars::xvel].nGrow();
    MultiFab::Copy(lev_old[Vars::cons],lev_new[Vars::cons],0,0,NVAR,ngs);
    MultiFab::Copy(lev_old[Vars::xvel],lev_new[Vars::xvel],0,0,1,ngvel);
    MultiFab::Copy(lev_old[Vars::yvel],lev_new[Vars::yvel],0,0,1,ngvel);
    MultiFab::Copy(lev_old[Vars::zvel],lev_new[Vars::zvel],0,0,1,ngvel);

    // also create the time integrator for this level
    integrator[lev] = std::make_unique<TimeIntegrator<Vector<MultiFab> > >(vars_old[lev]);
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

    {
        ParmParse pp("amr"); // Traditionally, these have prefix, amr.

        pp.query("interpolation_type", interpolation_type);

        pp.query("regrid_int", regrid_int);
        pp.query("plot_file", plot_file);
        pp.query("plot_int", plot_int);
        pp.query("check_file", check_file);
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
    }

    {
        ParmParse pp("erf");

        // Verbosity
        pp.query("v", verbose);

        // Frequency of diagnostic output
        pp.query("sum_interval", sum_interval);
        pp.query("sum_period"  , sum_per);

        // Time step controls
        pp.query("cfl", cfl);
        pp.query("init_shrink", init_shrink);
        pp.query("change_max", change_max);

        pp.query("fixed_dt", fixed_dt);

        AMREX_ASSERT(cfl > 0. || fixed_dt > 0.);
    }

    {  // Mesh refinement
        ParmParse pp("erf");

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
    }
    {  // Output format
        ParmParse pp("erf");
        pp.query("plotfile_type", plotfile_type);

        pp.query("output_1d_column", output_1d_column);
        pp.query("column_per", column_per);
        pp.query("column_interval", column_interval);
        pp.query("column_loc_x", column_loc_x);
        pp.query("column_loc_y", column_loc_y);
        pp.query("column_file_name", column_file_name);
    }


    //!don: set these properly
    lo_z_is_dirichlet = false;
    hi_z_is_dirichlet = false;

    solverChoice.init_params();
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
