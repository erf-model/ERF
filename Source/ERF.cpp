/**
 * \file ERF.cpp
 */

#include "prob_common.H"
#include <EOS.H>
#include <ERF.H>

#include <AMReX_buildInfo.H>

#include <Utils.H>
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

// Type of mesh refinement algorithm
std::string ERF::coupling_type = "OneWay";
int         ERF::do_reflux     = 0;
int         ERF::do_avg_down   = 0;

// Dictate verbosity in screen output
int         ERF::verbose       = 0;

// Use the native ERF MRI integrator
int         ERF::use_native_mri = 1;
int         ERF::no_substepping = 0;
int         ERF::force_stage1_single_substep = 1;

// Frequency of diagnostic output
int         ERF::sum_interval  = -1;
amrex::Real ERF::sum_per       = -1.0;

// Native AMReX vs NetCDF
std::string ERF::plotfile_type    = "amrex";

// init_type:  "ideal", "real", "input_sounding", "metgrid" or ""
std::string ERF::init_type        = "";

// NetCDF wrfinput (initialization) file(s)
amrex::Vector<amrex::Vector<std::string>> ERF::nc_init_file = {{""}}; // Must provide via input

// NetCDF wrfbdy (lateral boundary) file
std::string ERF::nc_bdy_file = ""; // Must provide via input

// Text input_sounding file
std::string ERF::input_sounding_file = "input_sounding";

// Flag to trigger initialization from input_sounding like WRF's ideal.exe
bool ERF::init_sounding_ideal = false;

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

    grids_to_evolve.resize(nlevs_max);

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);
    dt_mri_ratio.resize(nlevs_max, 1);

    vars_new.resize(nlevs_max);
    vars_old.resize(nlevs_max);

    rU_new.resize(nlevs_max);
    rV_new.resize(nlevs_max);
    rW_new.resize(nlevs_max);

    rU_old.resize(nlevs_max);
    rV_old.resize(nlevs_max);
    rW_old.resize(nlevs_max);

    for (int lev = 0; lev < nlevs_max; ++lev) {
        vars_new[lev].resize(Vars::NumTypes);
        vars_old[lev].resize(Vars::NumTypes);
    }

#if defined(ERF_USE_MOISTURE)
    qv.resize(nlevs_max);
    qc.resize(nlevs_max);
    qi.resize(nlevs_max);
    qrain.resize(nlevs_max);
    qsnow.resize(nlevs_max);
    qgraup.resize(nlevs_max);
#endif

    mri_integrator_mem.resize(nlevs_max);
    physbcs.resize(nlevs_max);

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

        post_timestep(step, cur_time, dt[0]);

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
        write_1D_profiles(time);
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

    // Moving terrain
    if ( solverChoice.use_terrain &&  (solverChoice.terrain_type == 1) )
    {
      for (int lev = finest_level; lev >= 0; lev--)
      {
        // Copy z_phs_nd and detJ_cc at end of timestep
        MultiFab::Copy(*z_phys_nd[lev], *z_phys_nd_new[lev], 0, 0, 1, z_phys_nd[lev]->nGrowVect());
        MultiFab::Copy(  *detJ_cc[lev],   *detJ_cc_new[lev], 0, 0, 1,   detJ_cc[lev]->nGrowVect());
        MultiFab::Copy(base_state[lev],base_state_new[lev],0,0,3,1);

        make_zcc(geom[lev],*z_phys_nd[lev],*z_phys_cc[lev]);
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

        if ( (init_type == "ideal" || init_type == "input_sounding") && solverChoice.use_terrain) {
            amrex::Abort("We do not currently support init_type = ideal or input_sounding with terrain");
        }

        if (!solverChoice.use_terrain && solverChoice.terrain_type != 0) {
            amrex::Abort("We do not allow terrain_type != 0 with use_terrain = false");
        }

        if (solverChoice.use_terrain) {
            if (init_type != "real" && init_type != "metgrid") {
                for (int lev = 0; lev <= finest_level; lev++)
                {
                    init_custom_terrain(geom[lev],*z_phys_nd[lev],time);
                    init_terrain_grid(geom[lev],*z_phys_nd[lev]);
                    make_J(geom[lev],*z_phys_nd[lev],*detJ_cc[lev]);
                    make_zcc(geom[lev],*z_phys_nd[lev],*z_phys_cc[lev]);
                }
            }
        }

        // Note that make_J and make_zcc area now called inside init_from_wrfinput
        for (int lev = 0; lev <= finest_level; lev++)
            init_only(lev, time);

        // For now we initialize rho_KE to 0
        Real RhoKE_0 = 0.0;
        ParmParse pp(pp_prefix);
        pp.query("RhoKE_0", RhoKE_0);
        int lb = std::max(finest_level-1,0);
        for (int lev(lb); lev >= 0; --lev)
            vars_new[lev][Vars::cons].setVal(RhoKE_0,RhoKE_comp,1,0);

        AverageDown();

    } else { // Restart from a checkpoint

        restart();

        if (solverChoice.use_terrain) {
            // This must come after the call to restart because that
            //      is where we read in the mesh data
            for (int lev = finest_level; lev >= 0; --lev) {
                make_J  (geom[lev],*z_phys_nd[lev],*detJ_cc[lev]);
                make_zcc(geom[lev],*z_phys_nd[lev],*z_phys_cc[lev]);
            }
        }
    }

    if (input_bndry_planes) {
        // Create the ReadBndryPlanes object so we can handle reading of boundary plane data
        amrex::Print() << "Defining r2d for the first time " << std::endl;
        m_r2d = std::make_unique< ReadBndryPlanes>(geom[0], solverChoice.rdOcp);

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

    // If we are reading initial data from wrfinput, the base state is defined there.
    // If we are reading initial data from metgrid output, the base state is defined there and we will rebalance
    //    after interpolation.
    // If we are reading initial data from an input_sounding, then the base state is calculated by
    //   InputSoundingData.calc_rho_p().
    if ( (init_type != "real") && (init_type != "metgrid") && (!init_sounding_ideal)) {
        initHSE();
    }

#ifdef ERF_USE_MOISTURE
    // Initialize microphysics here
    micro.define(solverChoice);

    // Call Init which will call Diagnose to fill qv, qc, qi
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        micro.Init(vars_new[lev][Vars::cons],
                   qc[lev],
                   qv[lev],
                   qi[lev],
                   Geom(lev),
                   0.0 // dummy value, not needed just to diagnose
                  );
        micro.Update(vars_new[lev][Vars::cons],
                     qv[lev],
                     qc[lev],
                     qi[lev],
                     qrain[lev],
                     qsnow[lev],
                     qgraup[lev]);

#if 0
        for (MFIter mfi(qv[lev]); mfi.isValid(); ++mfi) {
            const Box& box = mfi.growntilebox(IntVect(1,1,0));
            // const Box& box = mfi.tilebox();
            auto qv_arr = qv[lev].array(mfi);
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
               amrex::Print() << "QV " << IntVect(i,j,k) << " " << std::endl;
               amrex::Print() << "QV " << qv_arr(i,j,k) << " " << std::endl;
            });
         }
#endif
    }
#endif

    // Configure ABLMost params if used MostWall boundary condition
    // NOTE: we must set up the MOST routine before calling WritePlotFile because
    //       WritePlotFile calls FillPatch in order to compute gradients
    if (phys_bc_type[Orientation(Direction::z,Orientation::low)] == ERF_BC::MOST)
    {
        m_most = std::make_unique<ABLMost>(geom,vars_old,Theta_prim,z_phys_nd);

        // We now configure ABLMost params here so that we can print the averages at t=0
        // Note we don't fill ghost cells here because this is just for diagnostics
        int lev = 0; amrex::IntVect ng = IntVect(0,0,0);
        MultiFab S(vars_new[lev][Vars::cons],make_alias,0,2);
        MultiFab::Copy(  *Theta_prim[lev], S, Cons::RhoTheta, 0, 1, ng);
        MultiFab::Divide(*Theta_prim[lev], S, Cons::Rho     , 0, 1, ng);
        m_most->update_mac_ptrs(lev, vars_new, Theta_prim);
        m_most->update_fluxes(lev);
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
        if (init_type == "input_sounding")
        {
            // overwrite Ubar, Tbar, and thetabar with input profiles
            bool restarting = (restart_chkfile != "");
            setRayleighRefFromSounding(restarting);
        }

    }

    if (is_it_time_for_action(istep[0], t_new[0], dt[0], sum_interval, sum_per)) {
        sum_integrated_quantities(t_new[0]);
        write_1D_profiles(t_new[0]);
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

    ComputeDt();

    // Fill ghost cells/faces
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        auto& lev_new = vars_new[lev];

        FillPatch(lev, t_new[lev],
                  {&lev_new[Vars::cons],&lev_new[Vars::xvel],&lev_new[Vars::yvel],&lev_new[Vars::zvel]});

        // We need to fill the ghost cell values of the base state in case it wasn't
        //    done in the initialization
        base_state[lev].FillBoundary(geom[lev].periodicity());

        // For moving terrain only
        if (solverChoice.terrain_type > 0) {
            MultiFab::Copy(base_state_new[lev],base_state[lev],0,0,3,1);
            base_state_new[lev].FillBoundary(geom[lev].periodicity());
        }
    }

#ifdef ERF_USE_POISSON_SOLVE
    if (restart_chkfile == "")
    {
        // Note -- this projection is only defined for no terrain
        if (solverChoice.project_initial_velocity) {
            AMREX_ALWAYS_ASSERT(solverChoice.use_terrain == 0);
            project_initial_velocities();
        }
    }
#endif
    // Copy from new into old just in case
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        auto& lev_new = vars_new[lev];
        auto& lev_old = vars_old[lev];

        int ngs   = lev_new[Vars::cons].nGrow();
        int ngvel = lev_new[Vars::xvel].nGrow();

        MultiFab::Copy(lev_old[Vars::cons],lev_new[Vars::cons],0,0,NVAR,ngs);
        MultiFab::Copy(lev_old[Vars::xvel],lev_new[Vars::xvel],0,0,1,ngvel);
        MultiFab::Copy(lev_old[Vars::yvel],lev_new[Vars::yvel],0,0,1,ngvel);
        MultiFab::Copy(lev_old[Vars::zvel],lev_new[Vars::zvel],0,0,1,IntVect(ngvel,ngvel,0));
    }

    // Set these up here because we need to know which MPI rank "cell" is on...
    ParmParse pp("erf");
    if (pp.contains("data_log"))
    {
        int num_datalogs = pp.countval("data_log");
        datalog.resize(num_datalogs);
        datalogname.resize(num_datalogs);
        pp.queryarr("data_log",datalogname,0,num_datalogs);
        for (int i = 0; i < num_datalogs; i++)
            setRecordDataInfo(i,datalogname[i]);
    }

    if (pp.contains("sample_log") && pp.contains("sample_point"))
    {
        int lev = 0;

        int num_samplepts = pp.countval("sample_point") / AMREX_SPACEDIM;
        if (num_samplepts > 0) {
            Vector<int> index; index.resize(num_samplepts*AMREX_SPACEDIM);
            samplepoint.resize(num_samplepts);

            pp.queryarr("sample_point",index,0,num_samplepts*AMREX_SPACEDIM);
            for (int i = 0; i < num_samplepts; i++) {
                IntVect iv(index[AMREX_SPACEDIM*i+0],index[AMREX_SPACEDIM*i+1],index[AMREX_SPACEDIM*i+2]);
                samplepoint[i] = iv;
            }
        }

        int num_samplelogs = pp.countval("sample_log");
        AMREX_ALWAYS_ASSERT(num_samplelogs == num_samplepts);
        if (num_samplelogs > 0) {
            samplelog.resize(num_samplelogs);
            samplelogname.resize(num_samplelogs);
            pp.queryarr("sample_log",samplelogname,0,num_samplelogs);

            for (int i = 0; i < num_samplelogs; i++) {
                setRecordSampleInfo(i,lev,samplepoint[i],samplelogname[i]);
            }
        }

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

    // ********************************************************************************************
    // These are the persistent containers for the old and new data
    // ********************************************************************************************
    lev_new[Vars::cons].define(ba, dm, crse_new[Vars::cons].nComp(), crse_new[Vars::cons].nGrowVect());
    lev_old[Vars::cons].define(ba, dm, crse_new[Vars::cons].nComp(), crse_new[Vars::cons].nGrowVect());

    lev_new[Vars::xvel].define(convert(ba, IntVect(1,0,0)), dm, 1, crse_new[Vars::xvel].nGrowVect());
    lev_old[Vars::xvel].define(convert(ba, IntVect(1,0,0)), dm, 1, crse_new[Vars::xvel].nGrowVect());

    lev_new[Vars::yvel].define(convert(ba, IntVect(0,1,0)), dm, 1, crse_new[Vars::yvel].nGrowVect());
    lev_old[Vars::yvel].define(convert(ba, IntVect(0,1,0)), dm, 1, crse_new[Vars::yvel].nGrowVect());

    lev_new[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, crse_new[Vars::zvel].nGrowVect());
    lev_old[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, crse_new[Vars::zvel].nGrowVect());

    // ********************************************************************************************
    // These are just used for scratch in the time integrator but we might as well define them here
    // ********************************************************************************************
    rU_old[lev].define(convert(ba, IntVect(1,0,0)), dm, 1, crse_new[Vars::xvel].nGrowVect());
    rU_new[lev].define(convert(ba, IntVect(1,0,0)), dm, 1, crse_new[Vars::xvel].nGrowVect());

    rV_old[lev].define(convert(ba, IntVect(0,1,0)), dm, 1, crse_new[Vars::yvel].nGrowVect());
    rV_new[lev].define(convert(ba, IntVect(0,1,0)), dm, 1, crse_new[Vars::yvel].nGrowVect());

    rW_old[lev].define(convert(ba, IntVect(0,0,1)), dm, 1, crse_new[Vars::zvel].nGrowVect());
    rW_new[lev].define(convert(ba, IntVect(0,0,1)), dm, 1, crse_new[Vars::zvel].nGrowVect());

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    define_grids_to_evolve(lev);

    FillCoarsePatch(lev, time, {&lev_new[Vars::cons],&lev_new[Vars::xvel],
                                &lev_new[Vars::yvel],&lev_new[Vars::zvel]});

    initialize_integrator(lev, lev_new[Vars::cons], lev_new[Vars::xvel]);
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
ERF::RemakeLevel (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm)
{
    define_grids_to_evolve(lev);

    Vector<MultiFab> temp_lev_new(Vars::NumTypes);
    Vector<MultiFab> temp_lev_old(Vars::NumTypes);

    int ngrow_state = ComputeGhostCells(solverChoice.horiz_spatial_order, solverChoice.vert_spatial_order)+1;
    int ngrow_vels  = ComputeGhostCells(solverChoice.horiz_spatial_order, solverChoice.vert_spatial_order);

    temp_lev_new[Vars::cons].define(ba, dm, Cons::NumVars, ngrow_state);
    temp_lev_old[Vars::cons].define(ba, dm, Cons::NumVars, ngrow_state);

    temp_lev_new[Vars::xvel].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);
    temp_lev_old[Vars::xvel].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);

    temp_lev_new[Vars::yvel].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);
    temp_lev_old[Vars::yvel].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);

    temp_lev_new[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, IntVect(ngrow_vels,ngrow_vels,0));
    temp_lev_old[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, IntVect(ngrow_vels,ngrow_vels,0));

    // ********************************************************************************************
    // These are just used for scratch in the time integrator but we might as well define them here
    // ********************************************************************************************
    rU_old[lev].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);
    rU_new[lev].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);

    rV_old[lev].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);
    rV_new[lev].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);

    rW_old[lev].define(convert(ba, IntVect(0,0,1)), dm, 1, ngrow_vels);
    rW_new[lev].define(convert(ba, IntVect(0,0,1)), dm, 1, ngrow_vels);

    // ********************************************************************************************
    // This will fill the temporary MultiFabs with data from vars_new
    // ********************************************************************************************
    FillPatch(lev, time, {&temp_lev_new[Vars::cons],&temp_lev_new[Vars::xvel],
                          &temp_lev_new[Vars::yvel],&temp_lev_new[Vars::zvel]});

    // ********************************************************************************************
    // Copy from new into old just in case
    // ********************************************************************************************
    MultiFab::Copy(temp_lev_old[Vars::cons],temp_lev_new[Vars::cons],0,0,NVAR,ngrow_state);
    MultiFab::Copy(temp_lev_old[Vars::xvel],temp_lev_new[Vars::xvel],0,0,   1,ngrow_vels);
    MultiFab::Copy(temp_lev_old[Vars::yvel],temp_lev_new[Vars::yvel],0,0,   1,ngrow_vels);
    MultiFab::Copy(temp_lev_old[Vars::zvel],temp_lev_new[Vars::zvel],0,0,   1,IntVect(ngrow_vels,ngrow_vels,0));

    // ********************************************************************************************
    // Now swap the pointers
    // ********************************************************************************************
    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx) {
        std::swap(temp_lev_new[var_idx], vars_new[lev][var_idx]);
        std::swap(temp_lev_old[var_idx], vars_old[lev][var_idx]);
    }

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    initialize_integrator(lev, temp_lev_new[Vars::cons],temp_lev_new[Vars::xvel]);
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

    rU_new[lev].clear();
    rU_old[lev].clear();
    rV_new[lev].clear();
    rV_old[lev].clear();
    rW_new[lev].clear();
    rW_old[lev].clear();

    // Clears the integrator memory
    mri_integrator_mem[lev].reset();
    physbcs[lev].reset();

    grids_to_evolve[lev].clear();
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

    define_grids_to_evolve(lev);

    // The number of ghost cells for density must be 1 greater than that for velocity
    //     so that we can go back in forth betwen velocity and momentum on all faces
    int ngrow_state = ComputeGhostCells(solverChoice.horiz_spatial_order, solverChoice.vert_spatial_order)+1;
    int ngrow_vels  = ComputeGhostCells(solverChoice.horiz_spatial_order, solverChoice.vert_spatial_order);

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

    // ********************************************************************************************
    // These are just used for scratch in the time integrator but we might as well define them here
    // ********************************************************************************************
    rU_old[lev].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);
    rU_new[lev].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);

    rV_old[lev].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);
    rV_new[lev].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);

    rW_old[lev].define(convert(ba, IntVect(0,0,1)), dm, 1, ngrow_vels);
    rW_new[lev].define(convert(ba, IntVect(0,0,1)), dm, 1, ngrow_vels);

    //********************************************************************************************
    // Microphysics
    // *******************************************************************************************
#if defined(ERF_USE_MOISTURE)
    qv[lev].define(ba, dm, 1, ngrow_state);
    qc[lev].define(ba, dm, 1, ngrow_state);
    qi[lev].define(ba, dm, 1, ngrow_state);
    qrain[lev].define(ba, dm, 1, ngrow_state);
    qsnow[lev].define(ba, dm, 1, ngrow_state);
    qgraup[lev].define(ba, dm, 1, ngrow_state);
#endif

    // ********************************************************************************************
    // Diffusive terms
    // ********************************************************************************************
    bool l_use_terrain = solverChoice.use_terrain;
    bool l_use_diff    = ( (solverChoice.molec_diff_type != MolecDiffType::None) ||
                           (solverChoice.les_type        !=       LESType::None) ||
                           (solverChoice.pbl_type        !=       PBLType::None) );
    bool l_use_kturb   = ( (solverChoice.les_type != LESType::None)   ||
                           (solverChoice.pbl_type != PBLType::None) );
    bool l_use_ddorf   = (solverChoice.les_type == LESType::Deardorff);

    BoxArray ba12 = convert(ba, IntVect(1,1,0));
    BoxArray ba13 = convert(ba, IntVect(1,0,1));
    BoxArray ba23 = convert(ba, IntVect(0,1,1));

    Tau11_lev.resize(lev+1); Tau22_lev.resize(lev+1); Tau33_lev.resize(lev+1);
    Tau12_lev.resize(lev+1); Tau21_lev.resize(lev+1);
    Tau13_lev.resize(lev+1); Tau31_lev.resize(lev+1);
    Tau23_lev.resize(lev+1); Tau32_lev.resize(lev+1);

    eddyDiffs_lev.resize(lev+1);
    SmnSmn_lev.resize(lev+1);

    if (l_use_diff) {
        Tau11_lev[lev].reset( new MultiFab(ba  , dm, 1, IntVect(1,1,0)) );
        Tau22_lev[lev].reset( new MultiFab(ba  , dm, 1, IntVect(1,1,0)) );
        Tau33_lev[lev].reset( new MultiFab(ba  , dm, 1, IntVect(1,1,0)) );
        Tau12_lev[lev].reset( new MultiFab(ba12, dm, 1, IntVect(1,1,0)) );
        Tau13_lev[lev].reset( new MultiFab(ba13, dm, 1, IntVect(1,1,0)) );
        Tau23_lev[lev].reset( new MultiFab(ba23, dm, 1, IntVect(1,1,0)) );
        if (l_use_terrain) {
            Tau21_lev[lev].reset( new MultiFab(ba12, dm, 1, IntVect(1,1,0)) );
            Tau31_lev[lev].reset( new MultiFab(ba13, dm, 1, IntVect(1,1,0)) );
            Tau32_lev[lev].reset( new MultiFab(ba23, dm, 1, IntVect(1,1,0)) );
        } else {
            Tau21_lev[lev] = nullptr;
            Tau31_lev[lev] = nullptr;
            Tau32_lev[lev] = nullptr;
        }
    } else {
      Tau11_lev[lev] = nullptr; Tau22_lev[lev] = nullptr; Tau33_lev[lev] = nullptr;
      Tau12_lev[lev] = nullptr; Tau21_lev[lev] = nullptr;
      Tau13_lev[lev] = nullptr; Tau31_lev[lev] = nullptr;
      Tau23_lev[lev] = nullptr; Tau32_lev[lev] = nullptr;
    }

    if (l_use_kturb) {
      eddyDiffs_lev[lev].reset( new MultiFab(ba, dm, EddyDiff::NumDiffs, 1) );
      if(l_use_ddorf) {
          SmnSmn_lev[lev].reset( new MultiFab(ba, dm, 1, 0) );
      } else {
          SmnSmn_lev[lev] = nullptr;
      }
    } else {
      eddyDiffs_lev[lev] = nullptr;
      SmnSmn_lev[lev]    = nullptr;
    }

    // ********************************************************************************************
    // Metric terms
    // ********************************************************************************************
    z_phys_nd.resize(lev+1);
    z_phys_cc.resize(lev+1);
    detJ_cc.resize(lev+1);

    z_phys_nd_new.resize(lev+1);
    detJ_cc_new.resize(lev+1);

    z_phys_nd_src.resize(lev+1);
    detJ_cc_src.resize(lev+1);

    z_t_rk.resize(lev+1);

    // ********************************************************************************************
    // Map factors
    // ********************************************************************************************
    BoxList bl2d = ba.boxList();
    for (auto& b : bl2d) {
        b.setRange(2,0);
    }
    BoxArray ba2d(std::move(bl2d));

    mapfac_m.resize(lev+1);
    mapfac_u.resize(lev+1);
    mapfac_v.resize(lev+1);
    mapfac_m[lev].reset(new MultiFab(ba2d,dm,1,3));
    mapfac_u[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,3));
    mapfac_v[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,3));
    if(solverChoice.test_mapfactor) {
        mapfac_m[lev]->setVal(0.5);
        mapfac_u[lev]->setVal(0.5);
        mapfac_v[lev]->setVal(0.5);
    }
    else {
        mapfac_m[lev]->setVal(1.);
        mapfac_u[lev]->setVal(1.);
        mapfac_v[lev]->setVal(1.);
    }

    // ********************************************************************************************
    // Base state holds r_0, pres_0, pi_0 (in that order)
    // ********************************************************************************************
    base_state.resize(lev+1);
    base_state[lev].define(ba,dm,3,1);
    base_state[lev].setVal(0.);

    if (solverChoice.use_terrain && solverChoice.terrain_type > 0) {
        base_state_new.resize(lev+1);
        base_state_new[lev].define(ba,dm,3,1);
        base_state_new[lev].setVal(0.);
    }

    if (solverChoice.use_terrain) {
        z_phys_cc[lev].reset(new MultiFab(ba,dm,1,1));
          detJ_cc[lev].reset(new MultiFab(ba,dm,1,1));

        if (solverChoice.terrain_type > 0) {
            detJ_cc_new[lev].reset(new MultiFab(ba,dm,1,1));
            detJ_cc_src[lev].reset(new MultiFab(ba,dm,1,1));
            z_t_rk[lev].reset(new MultiFab( convert(ba, IntVect(0,0,1)), dm, 1, 1 ));
        }

        BoxArray ba_nd(ba);
        ba_nd.surroundingNodes();

        // We need this to be one greater than the ghost cells to handle levels > 0
        int ngrow = ComputeGhostCells(solverChoice.horiz_spatial_order, solverChoice.vert_spatial_order)+2;
        z_phys_nd[lev].reset(new MultiFab(ba_nd,dm,1,IntVect(ngrow,ngrow,1)));
        if (solverChoice.terrain_type > 0) {
            z_phys_nd_new[lev].reset(new MultiFab(ba_nd,dm,1,IntVect(ngrow,ngrow,1)));
            z_phys_nd_src[lev].reset(new MultiFab(ba_nd,dm,1,IntVect(ngrow,ngrow,1)));
        }
    } else {
            z_phys_nd[lev] = nullptr;
            z_phys_cc[lev] = nullptr;
              detJ_cc[lev] = nullptr;

        z_phys_nd_new[lev] = nullptr;
          detJ_cc_new[lev] = nullptr;

        z_phys_nd_src[lev] = nullptr;
          detJ_cc_src[lev] = nullptr;

               z_t_rk[lev] = nullptr;
    }

    // ********************************************************************************************
    // Define Theta_prim storage if using MOST BC
    // ********************************************************************************************
    Theta_prim.resize(lev+1);
    if (phys_bc_type[Orientation(Direction::z,Orientation::low)] == ERF_BC::MOST) {
      Theta_prim[lev].reset(new MultiFab(ba,dm,1,{ngrow_state,ngrow_state,0}));
    } else {
      Theta_prim[lev] = nullptr;
    }

    // ********************************************************************************************
    // Initialize the integrator class
    // ********************************************************************************************
    initialize_integrator(lev, lev_new[Vars::cons],lev_new[Vars::xvel]);
}

void
ERF::initialize_integrator(int lev, MultiFab& cons_mf, MultiFab& vel_mf)
{
    const BoxArray& ba(cons_mf.boxArray());
    const DistributionMapping& dm(cons_mf.DistributionMap());

    // Initialize the integrator memory
    int use_fluxes = (finest_level > 0);
    amrex::Vector<amrex::MultiFab> int_state; // integration state data structure example
    int_state.push_back(MultiFab(cons_mf, amrex::make_alias, 0, Cons::NumVars)); // cons
    int_state.push_back(MultiFab(convert(ba,IntVect(1,0,0)), dm, 1, vel_mf.nGrow())); // xmom
    int_state.push_back(MultiFab(convert(ba,IntVect(0,1,0)), dm, 1, vel_mf.nGrow())); // ymom
    int_state.push_back(MultiFab(convert(ba,IntVect(0,0,1)), dm, 1, vel_mf.nGrow())); // zmom
    if (use_fluxes) {
        int_state.push_back(MultiFab(convert(ba,IntVect(1,0,0)), dm, Cons::NumVars, 1)); // x-fluxes
        int_state.push_back(MultiFab(convert(ba,IntVect(0,1,0)), dm, Cons::NumVars, 1)); // y-fluxes
        int_state.push_back(MultiFab(convert(ba,IntVect(0,0,1)), dm, Cons::NumVars, 1)); // z-fluxes
    }

    mri_integrator_mem[lev] = std::make_unique<MRISplitIntegrator<amrex::Vector<amrex::MultiFab> > >(int_state);
    mri_integrator_mem[lev]->setNoSubstepping(no_substepping);
    mri_integrator_mem[lev]->setForceFirstStageSingleSubstep(force_stage1_single_substep);

    physbcs[lev] = std::make_unique<ERFPhysBCFunct> (lev, geom[lev], domain_bcs_type, domain_bcs_type_d,
                                                     solverChoice.terrain_type, m_bc_extdir_vals, m_bc_neumann_vals,
                                                     z_phys_nd[lev], detJ_cc[lev]);
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

#if defined(ERF_USE_MOISTURE)
    qv[lev].setVal(0.0);
    qc[lev].setVal(0.0);
    qi[lev].setVal(0.0);
    qrain[lev].setVal(0.0);
    qsnow[lev].setVal(0.0);
    qgraup[lev].setVal(0.0);
#endif

    // Initialize background flow (optional)
    if (init_type == "input_sounding") {
        init_from_input_sounding(lev);
#ifdef ERF_USE_NETCDF
    } else if (init_type == "ideal" || init_type == "real") {
        init_from_wrfinput(lev);
    } else if (init_type == "metgrid") {
        init_from_metgrid(lev);
#endif
    }

    // Add problem-specific flow features
    // If init_type is specified, then this is a perturbation to the background
    // flow from either input_sounding data or WRF WPS outputs (wrfinput_d0*)
    init_custom(lev);

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
    ParmParse pp_amr("amr");
    {
        // The type of the file we restart from
        pp.query("restart_type", restart_type);

        pp.query("regrid_int", regrid_int);
        pp.query("check_file", check_file);
        pp.query("check_type", check_type);

        // The regression tests use "amr.restart" and "amr.check_int" so we allow
        //    for those or "erf.restart" / "erf.check_int" with the former taking
        //    precedenceif both are specified
        pp.query("check_int", check_int);
        pp_amr.query("check_int", check_int);

        pp.query("restart", restart_chkfile);
        pp_amr.query("restart", restart_chkfile);

        // Verbosity
        pp.query("v", verbose);

        // Use the native ERF MRI integrator
        pp.query("use_native_mri", use_native_mri);
        pp.query("no_substepping", no_substepping);
        pp.query("force_stage1_single_substep", force_stage1_single_substep);

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

        // If this is set, it must be even
        if (fixed_mri_dt_ratio > 0 && (fixed_mri_dt_ratio%2 != 0) )
        {
            amrex::Abort("If you specify fixed_mri_dt_ratio, it must be even");
        }

        // If both fixed_dt and fast_dt are specified, their ratio must be an even integer
        if (fixed_dt > 0. && fixed_fast_dt > 0. && fixed_mri_dt_ratio <= 0)
        {
            Real eps = 1.e-12;
            int ratio = static_cast<int>( ( (1.0+eps) * fixed_dt ) / fixed_fast_dt );
            if (fixed_dt / fixed_fast_dt != ratio)
            {
                amrex::Abort("Ratio of fixed_dt to fixed_fast_dt must be an even integer");
            }
        }

        // If all three are specified, they must be consistent
        if (fixed_dt > 0. && fixed_fast_dt > 0. &&  fixed_mri_dt_ratio > 0)
        {
            if (fixed_dt / fixed_fast_dt != fixed_mri_dt_ratio)
            {
                amrex::Abort("Dt is over-specfied");
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
        if (init_type != "" &&
            init_type != "ideal" &&
            init_type != "real" &&
            init_type != "metgrid" &&
            init_type != "input_sounding")
        {
            amrex::Error("if specified, init_type must be ideal, real, metgrid or input_sounding");
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

        // Flag to trigger initialization from input_sounding like WRF's ideal.exe
        pp.query("init_sounding_ideal", init_sounding_ideal);

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

// Create horizontal average quantities for 5 variables:
//        density, temperature, pressure, qc, qv (if present)
void
ERF::MakeHorizontalAverages ()
{
    int lev = 0;

    // First, average down all levels
    AverageDown();

    MultiFab mf(grids[lev], dmap[lev], 5, 0);

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto  fab_arr = mf.array(mfi);
        auto cons_arr = vars_new[lev][Vars::cons].array(mfi);
#if defined(ERF_USE_MOISTURE)
        auto   qv_arr = qv[lev].array(mfi);
#endif
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Real dens = cons_arr(i, j, k, Cons::Rho);
            fab_arr(i, j, k, 0) = dens;
            fab_arr(i, j, k, 1) = cons_arr(i, j, k, Cons::RhoTheta) / dens;
#if defined(ERF_USE_MOISTURE)
            fab_arr(i, j, k, 2) = getPgivenRTh(cons_arr(i, j, k, Cons::RhoTheta), qv_arr(i,j,k));
#else
            fab_arr(i, j, k, 2) = getPgivenRTh(cons_arr(i, j, k, Cons::RhoTheta));
#endif
#if defined(ERF_USE_MOISTURE)
            fab_arr(i, j, k, 3) = cons_arr(i, j, k, Cons::RhoQt) / dens;
            fab_arr(i, j, k, 4) = cons_arr(i, j, k, Cons::RhoQp) / dens;
#endif
        });
    }

    int zdir = 2;
    auto domain = geom[0].Domain();

    // Sum in the horizontal plane
    Gpu::HostVector<Real> h_avg_density     = sumToLine(mf,0,1,domain,zdir);
    Gpu::HostVector<Real> h_avg_temperature = sumToLine(mf,1,1,domain,zdir);
    Gpu::HostVector<Real> h_avg_pressure    = sumToLine(mf,2,1,domain,zdir);
#ifdef ERF_USE_MOISTURE
    Gpu::HostVector<Real> h_avg_qv          = sumToLine(mf,3,1,domain,zdir);
    Gpu::HostVector<Real> h_avg_qc          = sumToLine(mf,4,1,domain,zdir);
#endif

    // Divide by the total number of cells we are averaging over
     int size_z = domain.length(zdir);
    Real area_z = static_cast<Real>(domain.length(0)*domain.length(1));
    int klen = static_cast<int>(h_avg_density.size());
    for (int k = 0; k < klen; ++k) {
        h_havg_density[k]     /= area_z;
        h_havg_temperature[k] /= area_z;
        h_havg_pressure[k]    /= area_z;
#if defined(ERF_USE_MOISTURE)
        h_havg_qc[k]          /= area_z;
        h_havg_qv[k]          /= area_z;
#endif
    }

    // resize device vectors
    d_havg_density.resize(size_z, 0.0_rt);
    d_havg_temperature.resize(size_z, 0.0_rt);
    d_havg_pressure.resize(size_z, 0.0_rt);
#if defined(ERF_USE_MOISTURE)
    d_havg_qv.resize(size_z, 0.0_rt);
    d_havg_qc.resize(size_z, 0.0_rt);
#endif

    // copy host vectors to device vectors
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_density.begin(), h_havg_density.end(), d_havg_density.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_temperature.begin(), h_havg_temperature.end(), d_havg_temperature.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_pressure.begin(), h_havg_pressure.end(), d_havg_pressure.begin());
#if defined(ERF_USE_MOISTURE)
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_qv.begin(), h_havg_qv.end(), d_havg_qv.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_qc.begin(), h_havg_qc.end(), d_havg_qc.begin());
#endif
}

// Create horizontal average quantities for the MultiFab passed in
// NOTE: this does not create device versions of the 1d arrays
void
ERF::MakeDiagnosticAverage (Vector<Real>& h_havg, MultiFab& S, int n)
{
    // Get the number of cells in z at level 0
    int dir_z = AMREX_SPACEDIM-1;
    auto domain = geom[0].Domain();
    int size_z = domain.length(dir_z);
    int start_z = domain.smallEnd()[dir_z];
    Real area_z = static_cast<Real>(domain.length(0)*domain.length(1));

    // resize the level 0 horizontal average vectors
    h_havg.resize(size_z, 0.0_rt);

    // Get the cell centered data and construct sums
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const IntVect& se = box.smallEnd();
        const IntVect& be = box.bigEnd();

        auto      fab_arr = S[mfi].array();

        FArrayBox fab_reduce(box, 1);
        Elixir elx_reduce = fab_reduce.elixir();
        auto arr_reduce   = fab_reduce.array();

        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            arr_reduce(i, j, k, 0) = fab_arr(i,j,k,n);
        });

        for (int k=se[dir_z]; k <= be[dir_z]; ++k) {
            Box kbox(box); kbox.setSmall(dir_z,k); kbox.setBig(dir_z,k);
            h_havg[k-start_z] += fab_reduce.sum<RunOn::Device>(kbox,0);
        }
    }

    // combine sums from different MPI ranks
    ParallelDescriptor::ReduceRealSum(h_havg.dataPtr(), h_havg.size());

    // divide by the total number of cells we are averaging over
    for (int k = 0; k < size_z; ++k) {
        h_havg[k]     /= area_z;
    }
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
ERF::define_grids_to_evolve (int lev)
{
   Box domain(geom[lev].Domain());
   if (lev == 0 && ( init_type == "real" || init_type == "metgrid" ) )
   {
      Box shrunk_domain(domain);
      shrunk_domain.grow(0,-1);
      shrunk_domain.grow(1,-1);
      grids_to_evolve[lev] = amrex::intersect(grids[lev],shrunk_domain);
   } else if (lev == 1) {
      Box shrunk_domain(boxes_at_level[lev][0]);
      shrunk_domain.grow(0,-1);
      shrunk_domain.grow(1,-1);
      grids_to_evolve[lev] = amrex::intersect(grids[lev],shrunk_domain);
#if 0
      if (num_boxes_at_level[lev] > 1) {
          for (int i = 1; i < num_boxes_at_level[lev]; i++) {
              Box shrunk_domain(boxes_at_level[lev][i]);
              shrunk_domain.grow(0,-1);
              shrunk_domain.grow(1,-1);
              grids_to_evolve[lev] = amrex::intersect(grids_to_evolve[lev],shrunk_domain);
          }
      }
#endif
   } else {
      // Just copy grids...
      grids_to_evolve[lev] = grids[lev];
   }
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

#if defined(ERF_USE_MOISTURE)
    qv.resize(nlevs_max);
    qc.resize(nlevs_max);
    qi.resize(nlevs_max);
    qrain.resize(nlevs_max);
    qsnow.resize(nlevs_max);
    qgraup.resize(nlevs_max);
#endif

    mri_integrator_mem.resize(nlevs_max);
    physbcs.resize(nlevs_max);

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

        post_timestep(step, cur_time, dt[0]);

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
