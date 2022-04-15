/**
 * \file ERF.cpp
 */

#include <ERF.H>
#include <prob_common.H>
#include <AMReX_buildInfo.H>
#include <EOS.H>

using namespace amrex;

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

// init_type:  "ideal" vs "real"
std::string ERF::init_type        = "ideal";

// NetCDF initialization file
std::string ERF::nc_init_file = ""; // Must provide via input

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
    dt_mri_ratio.resize(nlevs_max, 1);

    vars_new.resize(nlevs_max);
    vars_old.resize(nlevs_max);

    for (int lev = 0; lev < nlevs_max; ++lev) {
        vars_new[lev].resize(Vars::NumTypes);
        vars_old[lev].resize(Vars::NumTypes);
    }

#ifdef ERF_USE_TERRAIN
    z_phys_nd.resize(nlevs_max);
    z_phys_cc.resize(nlevs_max);
    detJ_cc.resize(nlevs_max);
    dens_hse.resize(nlevs_max);
    pres_hse.resize(nlevs_max);
#endif

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

        if (plot_int > 0 && (step+1) % plot_int == 0) {
            last_plot_file_step = step+1;
            WritePlotFile();
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

    if (plot_int > 0 && istep[0] > last_plot_file_step) {
        WritePlotFile();
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

    if (input_bndry_planes) {
        // Create the ReadBndryPlanes object so we can handle reading of boundary plane data
        amrex::Print() << "Defining r2d for the first time " << std::endl;
        m_r2d = std::make_unique< ReadBndryPlanes>(geom[0]);

        // Read the "time.dat" file to know what data is available
        m_r2d->read_time_file();
    }

    // Map the words in the inputs file to BC types, then translate
    //     those types into what they mean for each variable
    init_bcs();

    // Verify BCs are compatible sith solver choice
    if (solverChoice.pbl_type == PBLType::MYNN25 &&
        phys_bc_type[Orientation(Direction::z,Orientation::low)] != BC::MOST) {
        amrex::Abort("MYNN2.5 PBL Model requires MOST at lower boundary");
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

    } else { // Restart from a checkpoint
#ifdef ERF_USE_NETCDF
        if (plot_type == "netcdf") {
           ReadNCCheckpointFile();
        }
#endif
        if (plot_type == "native") {
           ReadCheckpointFile();
        }
        // We set this here so that we don't over-write the checkpoint file we just started from
        last_check_file_step = istep[0];
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

#ifdef ERF_USE_TERRAIN
    pres_hse[lev].define(ba,dm,1,1);
    dens_hse[lev].define(ba,dm,1,1);
    z_phys_cc[lev].define(ba,dm,1,0);
    detJ_cc[lev].define(ba,dm,1,0);

    BoxArray ba_nd(ba);
    ba_nd.surroundingNodes();
    z_phys_nd[lev].define(ba_nd,dm,1,0);
#endif

    // Read enough of the boundary plane data to make it through the FillPatch at this time
    if (input_bndry_planes)
    {
        amrex::Real dt_dummy = 1.e-200;
        m_r2d->read_input_files(time,dt_dummy,m_bc_extdir_vals);
    }

#ifdef ERF_USE_NETCDF
    // We only want to read the file once -- here we fill one FArrayBox (per variable) that spans the domain
    if (lev == 0) {
        if (init_type == "real") {
            if (nc_init_file.empty())
                amrex::Error("NetCDF initialization file name must be provided via input");
            read_from_netcdf();
        }
    }
#endif

    // Loop over grids at this level to initialize our grid data
    lev_new[Vars::cons].setVal(0.0);
    lev_new[Vars::xvel].setVal(0.0);
    lev_new[Vars::yvel].setVal(0.0);
    lev_new[Vars::zvel].setVal(0.0);

    if (init_type == "ideal") {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
          const Box& bx = mfi.tilebox();
          const auto& cons_arr = lev_new[Vars::cons].array(mfi);
          const auto& xvel_arr = lev_new[Vars::xvel].array(mfi);
          const auto& yvel_arr = lev_new[Vars::yvel].array(mfi);
          const auto& zvel_arr = lev_new[Vars::zvel].array(mfi);

          erf_init_prob(bx, cons_arr, xvel_arr, yvel_arr, zvel_arr, geom[lev].data());
        }
#ifdef ERF_USE_TERRAIN
        init_ideal_terrain(lev);
        make_metrics(lev);
#endif
#ifdef ERF_USE_NETCDF
    } else if (init_type == "real") {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
          const Box& bx = mfi.tilebox();
          FArrayBox& cons_fab = lev_new[Vars::cons][mfi];
          FArrayBox& xvel_fab = lev_new[Vars::xvel][mfi];
          FArrayBox& yvel_fab = lev_new[Vars::yvel][mfi];
          FArrayBox& zvel_fab = lev_new[Vars::zvel][mfi];

#ifdef ERF_USE_TERRAIN
          FArrayBox& z_phys_nd_fab = z_phys_nd[mfi];
          init_from_netcdf(bx, cons_fab, xvel_fab, yvel_fab, zvel_fab,
                           z_phys_nd_fab);
#else
          init_from_netcdf(bx, cons_fab, xvel_fab, yvel_fab, zvel_fab);
#endif
        }
#ifdef ERF_USE_TERRAIN
        make_metrics(lev);
#endif
#endif
    }

    // Ensure that the face-based data are the same on both sides of a periodic domain.
    // The data associated with the lower grid ID is considered the correct value.
    lev_new[Vars::xvel].OverrideSync(geom[lev].periodicity());
    lev_new[Vars::yvel].OverrideSync(geom[lev].periodicity());
    lev_new[Vars::zvel].OverrideSync(geom[lev].periodicity());

    // configure ABLMost params if used MostWall boundary condition
    for (OrientationIter oitr; oitr; ++oitr) {
        const Orientation face = oitr();
        if (phys_bc_type[face] == BC::MOST && lev == 0) setupABLMost(lev);
    }

    // Fill ghost cells/faces
    FillPatch(lev, time, lev_new[Vars::cons], 0, Cons::NumVars, Vars::cons);
    FillPatch(lev, time, lev_new[Vars::xvel], 0, 1, Vars::xvel);
    FillPatch(lev, time, lev_new[Vars::yvel], 0, 1, Vars::yvel);
    FillPatch(lev, time, lev_new[Vars::zvel], 0, 1, Vars::zvel);

    // Copy from new into old just in case
    int ngs   = lev_new[Vars::cons].nGrow();
    int ngvel = lev_new[Vars::xvel].nGrow();
    MultiFab::Copy(lev_old[Vars::cons],lev_new[Vars::cons],0,0,NVAR,ngs);
    MultiFab::Copy(lev_old[Vars::xvel],lev_new[Vars::xvel],0,0,1,ngvel);
    MultiFab::Copy(lev_old[Vars::yvel],lev_new[Vars::yvel],0,0,1,ngvel);
    MultiFab::Copy(lev_old[Vars::zvel],lev_new[Vars::zvel],0,0,1,ngvel);
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

        pp.query("regrid_int", regrid_int);
        pp.query("plot_file", plot_file);
        pp.query("plot_type", plot_type);
        pp.query("plot_int", plot_int);
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
    }

    {
        ParmParse pp("erf");

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
    }

    {  // How to initialize
        ParmParse pp("erf");
        pp.query("init_type",init_type);
        if (init_type != "ideal" && init_type != "real")
        {
            amrex::Error("init_type must be ideal or real");
        }

        // NetCDF initialization file
        pp.query("nc_init_file", nc_init_file);
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

        // Specify information about outputting planes of data
        pp.query("output_bndry_planes", output_bndry_planes);
        pp.query("bndry_output_planes_interval", bndry_output_planes_interval);
        pp.query("bndry_output_planes_per", bndry_output_planes_per);
        pp.query("bndry_output_start_time", bndry_output_planes_start_time);

        // Specify whether ingest boundary planes of data
        pp.query("input_bndry_planes", input_bndry_planes);
    }

    solverChoice.init_params();
}

// Create horizontal average quantities
void
ERF::MakeHorizontalAverages ()
{
    // We need to create horizontal averages for 5 variables:
    // density, temperature, pressure, qc, qv

    // First, average down all levels
    AverageDown();

    // get the number of cells in z at level 0
    int dir_z = AMREX_SPACEDIM-1;
    auto domain = geom[0].Domain();
    int size_z = domain.length(dir_z);
    int start_z = domain.smallEnd()[dir_z];

    Real area_z = static_cast<Real>(domain.length(0));
#if (AMREX_SPACEDIM == 3)
    area_z *= domain.length(1);
#endif

    // resize the level 0 horizontal average vectors
    h_havg_density.resize(size_z, 0.0_rt);
    h_havg_temperature.resize(size_z, 0.0_rt);
    h_havg_pressure.resize(size_z, 0.0_rt);
    h_havg_qv.resize(size_z, 0.0_rt);
    h_havg_qc.resize(size_z, 0.0_rt);

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

        FArrayBox fab_reduce(box, 5);
        Elixir elx_reduce = fab_reduce.elixir();
        auto arr_reduce = fab_reduce.array();

        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Real dens = arr_cons(i, j, k, Cons::Rho);
            arr_reduce(i, j, k, 0) = dens;
            arr_reduce(i, j, k, 1) = arr_cons(i, j, k, Cons::RhoTheta) / dens;
            arr_reduce(i, j, k, 2) = getPgivenRTh(arr_cons(i, j, k, Cons::RhoTheta));
            arr_reduce(i, j, k, 3) = arr_cons(i, j, k, Cons::RhoQv) / dens;
            arr_reduce(i, j, k, 4) = arr_cons(i, j, k, Cons::RhoQc) / dens;
        });

        for (int k=se[dir_z]; k <= be[dir_z]; ++k) {
            Box kbox(box); kbox.setSmall(dir_z,k); kbox.setBig(dir_z,k);
            h_havg_density     [k-start_z] += fab_reduce.sum<RunOn::Device>(kbox,0);
            h_havg_temperature [k-start_z] += fab_reduce.sum<RunOn::Device>(kbox,1);
            h_havg_pressure    [k-start_z] += fab_reduce.sum<RunOn::Device>(kbox,2);
            h_havg_qv          [k-start_z] += fab_reduce.sum<RunOn::Device>(kbox,3);
            h_havg_qc          [k-start_z] += fab_reduce.sum<RunOn::Device>(kbox,4);
        }
    }

    // combine sums from different MPI ranks
    ParallelDescriptor::ReduceRealSum(h_havg_density.dataPtr(), h_havg_density.size());
    ParallelDescriptor::ReduceRealSum(h_havg_temperature.dataPtr(), h_havg_temperature.size());
    ParallelDescriptor::ReduceRealSum(h_havg_pressure.dataPtr(), h_havg_pressure.size());
    ParallelDescriptor::ReduceRealSum(h_havg_qv.dataPtr(), h_havg_qv.size());
    ParallelDescriptor::ReduceRealSum(h_havg_qc.dataPtr(), h_havg_qc.size());

    // divide by the total number of cells we are averaging over
    for (int k = 0; k < size_z; ++k) {
        h_havg_density[k]     /= area_z;
        h_havg_temperature[k] /= area_z;
        h_havg_pressure[k]    /= area_z;
        h_havg_qv[k]          /= area_z;
        h_havg_qc[k]          /= area_z;
    }

    // resize device vectors
    d_havg_density.resize(size_z, 0.0_rt);
    d_havg_temperature.resize(size_z, 0.0_rt);
    d_havg_pressure.resize(size_z, 0.0_rt);
    d_havg_qv.resize(size_z, 0.0_rt);
    d_havg_qc.resize(size_z, 0.0_rt);

    // copy host vectors to device vectors
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_density.begin(), h_havg_density.end(), d_havg_density.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_temperature.begin(), h_havg_temperature.end(), d_havg_temperature.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_pressure.begin(), h_havg_pressure.end(), d_havg_pressure.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_qv.begin(), h_havg_qv.end(), d_havg_qv.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_qc.begin(), h_havg_qc.end(), d_havg_qc.begin());
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
    amrex::ParmParse pp("erf");

    pp.query("most.surf_temp", most.surf_temp);
    pp.query("most.zref", most.zref);
    pp.query("most.z0", most.z0);

    MultiFab& S_new = vars_new[lev][Vars::cons];
    MultiFab& U_new = vars_new[lev][Vars::xvel];
    MultiFab& V_new = vars_new[lev][Vars::yvel];
    MultiFab& W_new = vars_new[lev][Vars::zvel];

    PlaneAverage save (&S_new, geom[0], 2, true);
    PlaneAverage vxave(&U_new, geom[0], 2, true);
    PlaneAverage vyave(&V_new, geom[0], 2, true);
    PlaneAverage vzave(&W_new, geom[0], 2, true);
    VelPlaneAverage vmagave({&U_new,&V_new,&W_new}, geom[0], 2, true);

    save. compute_averages(ZDir(), save.field());
    vxave.compute_averages(ZDir(), vxave.field());
    vyave.compute_averages(ZDir(), vyave.field());
    vzave.compute_averages(ZDir(), vzave.field());
    vmagave.compute_hvelmag_averages(ZDir(), 0, 1, vmagave.field());

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom[0].CellSizeArray();

    most.vel_mean[0] = vxave.line_average_interpolated(most.zref, 0);
    most.vel_mean[1] = vyave.line_average_interpolated(most.zref, 0);
    most.vel_mean[2] = vzave.line_average_interpolated(most.zref, 0);
    most.vmag_mean   = vmagave.line_hvelmag_average_interpolated(most.zref);
    most.theta_mean  = save.line_average_interpolated(most.zref, Cons::RhoTheta);

    printf("vmag_mean=%13.6e,theta_mean=%13.6e, zref=%13.6e, dx=%13.6e\n",most.vmag_mean,most.theta_mean,most.zref,dx[2]);

    most.update_fluxes();
}

#ifdef ERF_USE_NETCDF
void
ERF::read_from_netcdf()
{
    auto domain = geom[0].Domain();

    // We allocate these here so they exist on all ranks
    Box ubx(domain); ubx.surroundingNodes(0); NC_xvel_fab.resize(ubx,1);
    Box vbx(domain); vbx.surroundingNodes(1); NC_yvel_fab.resize(vbx,1);
    Box wbx(domain); wbx.surroundingNodes(2); NC_zvel_fab.resize(wbx,1);
    NC_rho_fab.resize(domain,1);
    NC_rhotheta_fab.resize(domain,1);

#ifdef ERF_USE_TERRAIN
    NC_PH_fab.resize(wbx,1);
    NC_PHB_fab.resize(wbx,1);
#endif

    if (ParallelDescriptor::IOProcessor())
    {
        Vector<FArrayBox*> NC_fabs;
        Vector<std::string> NC_names;

        NC_fabs.push_back(&NC_xvel_fab    ); NC_names.push_back("U");
        NC_fabs.push_back(&NC_yvel_fab)    ; NC_names.push_back("V");
        NC_fabs.push_back(&NC_zvel_fab)    ; NC_names.push_back("W");
        NC_fabs.push_back(&NC_rho_fab)     ; NC_names.push_back("ALB");
        NC_fabs.push_back(&NC_rhotheta_fab); NC_names.push_back("T_INIT");
#ifdef ERF_USE_TERRAIN
        NC_fabs.push_back(&NC_PH_fab); NC_names.push_back("PH");
        NC_fabs.push_back(&NC_PHB_fab); NC_names.push_back("PHB");
#endif

        // Read the netcdf file and fill these FABs
        BuildFABsFromIdealOutputFile(nc_init_file, NC_names, NC_fabs, NC_Data_Dims_Type::Time_BT_SN_WE);

        // Convert to rho by inverting
        NC_rho_fab.invert(1.0);

        // The ideal.exe NetCDF file has this ref value subtracted from theta or T_INIT. Need to add in ERF.
        const Real theta_ref = 300.0;
        NC_rhotheta_fab.plus(theta_ref);

        // Now multiply by rho to get (rho theta) instead of theta
        NC_rhotheta_fab.mult(NC_rho_fab,0,0,1);

#if 0
        // Read the perturbation density inverse
        BuildFabFromIdealOutputFile(nc_init_file, "AL", NC_rho_inv_pert_fab,
                                    NC_Data_Dims_Type::Time_BT_SN_WE);

        // Read Base Pressure
        BuildFabFromIdealOutputFile(nc_init_file, "PB", NC_p_base_fab,
                                    NC_Data_Dims_Type::Time_BT_SN_WE);

        // Read Perturbation Pressure
        BuildFabFromIdealOutputFile(nc_init_file, "P", NC_p_pert_fab,
                                    NC_Data_Dims_Type::Time_BT_SN_WE);

        // Read Surface Pressure
        BuildFabFromIdealOutputFile(nc_init_file, "PSFC", NC_p_surf_fab,
                                    NC_Data_Dims_Type::Time_SN_WE);

        // Read eta at cell-center levels
        // Likewise, we can read eta at z-stagger levels, but since pressure is defined at cell-centers we skip that
        // eta = (p(z) - p_top)/(p_surf - p_top) => 1 at surface and 0 at top
        BuildFabFromIdealOutputFile(nc_init_file, "ZNU", NC_eta_fab,
                                    NC_Data_Dims_Type::Time_BT);

        // Read base-state geopotential
        BuildFabFromIdealOutputFile(nc_init_file, "PHB", NC_phb_fab,
                                    NC_Data_Dims_Type::Time_BT_SN_WE);
#endif

#if 0
    //Construct physical z from base-state geopotential. We can possibly include pertubation geopotential
    BoxArray ba(phb_nc);
    DistributionMapping dm { ba };
    z_phy_mf.define(ba, dm, 1, 0);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(z_phy_mf); mfi.isValid(); ++mfi )
    {
        const Box& z_tilebx = mfi.tilebox();
        z_phy_arr = z_phy_mf.array(mfi);
        amrex::ParallelFor(z_tilebx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            z_phy_arr(i, j, k) = phb_arr(i, j, k) / CONST_GRAV;
        });
        // Alternatively, we can compute z = f(eta, pb+p, p_surf, theta)
    }
#endif

    } // if ParalleDescriptor::IOProcessor()

    // We put a barrier here so the rest of the processors wait to do anything until they have the data
    amrex::ParallelDescriptor::Barrier();

    // When an FArrayBox is built, space is allocated on every rank.  However, we only
    //    filled the data in these FABs on the IOProcessor.  So here we broadcast
    //    the data to every rank.

    int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank
    ParallelDescriptor::Bcast(NC_xvel_fab.dataPtr(),NC_xvel_fab.box().numPts(),ioproc);
    ParallelDescriptor::Bcast(NC_yvel_fab.dataPtr(),NC_yvel_fab.box().numPts(),ioproc);
    ParallelDescriptor::Bcast(NC_zvel_fab.dataPtr(),NC_zvel_fab.box().numPts(),ioproc);
    ParallelDescriptor::Bcast(NC_rho_fab.dataPtr(),NC_rho_fab.box().numPts(),ioproc);
    ParallelDescriptor::Bcast(NC_rhotheta_fab.dataPtr(),NC_rhotheta_fab.box().numPts(),ioproc);
#ifdef ERF_USE_TERRAIN
    ParallelDescriptor::Bcast(NC_PHB_fab.dataPtr(),NC_PHB_fab.box().numPts(),ioproc);
    ParallelDescriptor::Bcast(NC_PH_fab.dataPtr() ,NC_PH_fab.box().numPts() ,ioproc);
#endif

#if 0
    ParallelDescriptor::Bcast(NC_rho_inv_pert_fab.dataPtr(),NC_rho_inv_pert_fab.box().numPts(),ioproc);
    ParallelDescriptor::Bcast(NC_p_base_fab.dataPtr(),NC_p_base_fab.box().numPts(),ioproc);
    ParallelDescriptor::Bcast(NC_p_pert_fab.dataPtr(),NC_p_pert_fab.box().numPts(),ioproc);
    ParallelDescriptor::Bcast(NC_p_surf_fab.dataPtr(),NC_p_surf_fab.box().numPts(),ioproc);
    ParallelDescriptor::Bcast(NC_p_surf_fab.dataPtr(),NC_p_surf_fab.box().numPts(),ioproc);
    ParallelDescriptor::Bcast(NC_eta_fab.dataPtr(),NC_eta_fab.box().numPts(),ioproc);
#endif
    amrex::Print() << "Successfully loaded data from the 'ideal.exe' / 'real.exe' output NetCDF file" << std::endl << std::endl;
}

void
ERF::init_from_netcdf(const amrex::Box& bx, FArrayBox& state,
                      FArrayBox& x_vel, FArrayBox& y_vel, FArrayBox& z_vel
#ifdef ERF_USE_TERRAIN
                     ,FArrayBox& z_phys
#endif
                                       )
{
    //
    // FArrayBox to FArrayBox copy does "copy on intersection"
    // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
    //
    // This copies x-vel
    x_vel.copy(NC_xvel_fab);

    // This copies y-vel
    y_vel.copy(NC_yvel_fab);

    // This copies z-vel
    z_vel.copy(NC_zvel_fab);

    // We first initialize all state variables to zero
    state.setVal(0.);

    // This copies the density
    state.copy(NC_rho_fab,0,Rho_comp,1);

    // This copies (rho*theta)
    state.copy(NC_rhotheta_fab,0,RhoTheta_comp,1);

#ifdef ERF_USE_TERRAIN
    // This copies from NC_zphys on z-faces to z_phys_nd on nodes
    Array4<Real>    z_arr = z_phys.array();
    Array4<Real> nc_phb_arr = NC_PHB.array();
    Array4<Real> nc_ph_arr  = NC_PH.array();
    const Box& bx(z_phys.box());
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        z_phys_arr(i, j, k) = 0.25 * ( nc_ph_arr(i,j,k) +  nc_ph_arr(i-1,j,k) +  nc_ph_arr(i,j-1,k) +  nc_ph_arr(i-1,j-1,k) +
                                      nc_phb_arr(i,j,k) + nc_phb_arr(i-1,j,k) + nc_phb_arr(i,j-1,k) + nc_phb_arr(i-1,j-1,k) );
    });
#endif
}
#endif
