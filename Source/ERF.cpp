/**
 * \file ERF.cpp
 */

/**
 * Main class in ERF code, instantiated from main.cpp
*/


#include <ERF_EOS.H>
#include <ERF.H>
#include <AMReX_buildInfo.H>
#include <ERF_Utils.H>
#include <ERF_TerrainMetrics.H>
#include <memory>

using namespace amrex;

Real ERF::startCPUTime        = 0.0;
Real ERF::previousCPUTimeUsed = 0.0;

Vector<AMRErrorTag> ERF::ref_tags;

SolverChoice ERF::solverChoice;

// Time step control
Real ERF::cfl           =  0.8;
Real ERF::init_shrink   =  1.0;
Real ERF::change_max    =  1.1;
int  ERF::fixed_mri_dt_ratio = 0;

// Dictate verbosity in screen output
int ERF::verbose        = 0;
int  ERF::mg_verbose    = 0;
bool ERF::use_fft       = false;

// Frequency of diagnostic output
int  ERF::sum_interval  = -1;
Real ERF::sum_per       = -1.0;

int  ERF::pert_interval = -1;

// Native AMReX vs NetCDF
std::string ERF::plotfile_type    = "amrex";

InitType ERF::init_type;

// use_real_bcs: only true if 1) ( (init_type == InitType::Real) or (init_type == InitGrid::Metgrid) )
//                        AND 2) we want to use the bc's from the WRF bdy file
bool ERF::use_real_bcs;

// NetCDF wrfinput (initialization) file(s)
Vector<Vector<std::string>> ERF::nc_init_file = {{""}}; // Must provide via input

// NetCDF wrfbdy (lateral boundary) file
std::string ERF::nc_bdy_file; // Must provide via input

// Flag to trigger initialization from input_sounding like WRF's ideal.exe
bool ERF::init_sounding_ideal = false;

// 1D NetCDF output (for ingestion by AMR-Wind)
int  ERF::output_1d_column = 0;
int  ERF::column_interval  = -1;
Real ERF::column_per       = -1.0;
Real ERF::column_loc_x     = 0.0;
Real ERF::column_loc_y     = 0.0;
std::string ERF::column_file_name = "column_data.nc";

// 2D BndryRegister output (for ingestion by AMR-Wind)
int  ERF::output_bndry_planes            = 0;
int  ERF::bndry_output_planes_interval   = -1;
Real ERF::bndry_output_planes_per        = -1.0;
Real ERF::bndry_output_planes_start_time =  0.0;

// 2D BndryRegister input
int  ERF::input_bndry_planes             = 0;

Vector<std::string> BCNames = {"xlo", "ylo", "zlo", "xhi", "yhi", "zhi"};

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
ERF::ERF ()
{
    ERF_shared();
}

void
ERF::ERF_shared ()
{
    if (ParallelDescriptor::IOProcessor()) {
        const char* erf_hash = buildInfoGetGitHash(1);
        const char* amrex_hash = buildInfoGetGitHash(2);
        const char* buildgithash = buildInfoGetBuildGitHash();
        const char* buildgitname = buildInfoGetBuildGitName();

        if (strlen(erf_hash) > 0) {
          Print() << "\n"
                         << "ERF git hash: " << erf_hash << "\n";
        }
        if (strlen(amrex_hash) > 0) {
          Print() << "AMReX git hash: " << amrex_hash << "\n";
        }
        if (strlen(buildgithash) > 0) {
          Print() << buildgitname << " git hash: " << buildgithash << "\n";
        }

        Print() << "\n";
    }

    int nlevs_max = max_level + 1;

#ifdef ERF_USE_WINDFARM
    Nturb.resize(nlevs_max);
    vars_windfarm.resize(nlevs_max);
    SMark.resize(nlevs_max);
#endif

#if defined(ERF_USE_RRTMGP)
    qheating_rates.resize(nlevs_max);
    sw_lw_fluxes.resize(nlevs_max);
    solar_zenith.resize(nlevs_max);
#endif

    // NOTE: size lsm before readparams (chooses the model at all levels)
    lsm.ReSize(nlevs_max);
    lsm_data.resize(nlevs_max);
    lsm_flux.resize(nlevs_max);

    ReadParameters();
    initializeMicrophysics(nlevs_max);

#ifdef ERF_USE_WINDFARM
    initializeWindFarm(nlevs_max);
#endif

    const std::string& pv1 = "plot_vars_1"; setPlotVariables(pv1,plot_var_names_1);
    const std::string& pv2 = "plot_vars_2"; setPlotVariables(pv2,plot_var_names_2);

    // Initialize staggered vertical levels for grid stretching or terrain, and
    // to simplify Rayleigh damping layer calculations.
    zlevels_stag.resize(max_level+1);
    init_zlevels(zlevels_stag,
                 geom,
                 refRatio(),
                 solverChoice.grid_stretching_ratio,
                 solverChoice.zsurf,
                 solverChoice.dz0);

    if (solverChoice.use_terrain) {
        int nz = geom[0].Domain().length(2) + 1; // staggered
        if (std::fabs(zlevels_stag[0][nz-1]-geom[0].ProbHi(2)) > 1.0e-4) {
            Print() << "Note: prob_hi[2]=" << geom[0].ProbHi(2)
                << " does not match highest requested z level " << zlevels_stag[0][nz-1]
                << std::endl;
        }
        if (std::fabs(zlevels_stag[0][0]-geom[0].ProbLo(2)) > 1.0e-4) {
            Print() << "Note: prob_lo[2]=" << geom[0].ProbLo(2)
                << " does not match lowest requested level " << zlevels_stag[0][0]
                << std::endl;
        }

        // Redefine the problem domain here?
    }

    prob = amrex_probinit(geom[0].ProbLo(),geom[0].ProbHi());

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

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

    // We resize this regardless in order to pass it without error
    pp_inc.resize(nlevs_max);

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

    // Time integrator
    mri_integrator_mem.resize(nlevs_max);

    // Physical boundary conditions
    physbcs_cons.resize(nlevs_max);
    physbcs_u.resize(nlevs_max);
    physbcs_v.resize(nlevs_max);
    physbcs_w.resize(nlevs_max);
    physbcs_base.resize(nlevs_max);

    // Planes to hold Dirichlet values at boundaries
    xvel_bc_data.resize(nlevs_max);
    yvel_bc_data.resize(nlevs_max);
    zvel_bc_data.resize(nlevs_max);

    advflux_reg.resize(nlevs_max);

    // Stresses
    Tau11_lev.resize(nlevs_max); Tau22_lev.resize(nlevs_max); Tau33_lev.resize(nlevs_max);
    Tau12_lev.resize(nlevs_max); Tau21_lev.resize(nlevs_max);
    Tau13_lev.resize(nlevs_max); Tau31_lev.resize(nlevs_max);
    Tau23_lev.resize(nlevs_max); Tau32_lev.resize(nlevs_max);
    SFS_hfx1_lev.resize(nlevs_max); SFS_hfx2_lev.resize(nlevs_max); SFS_hfx3_lev.resize(nlevs_max);
    SFS_diss_lev.resize(nlevs_max);
    SFS_q1fx1_lev.resize(nlevs_max); SFS_q1fx2_lev.resize(nlevs_max); SFS_q1fx3_lev.resize(nlevs_max);
    SFS_q2fx3_lev.resize(nlevs_max);
    eddyDiffs_lev.resize(nlevs_max);
    SmnSmn_lev.resize(nlevs_max);

    // Sea surface temps
    sst_lev.resize(nlevs_max);
    lmask_lev.resize(nlevs_max);

    // Metric terms
    z_phys_nd.resize(nlevs_max);
    z_phys_cc.resize(nlevs_max);
    detJ_cc.resize(nlevs_max);
    ax.resize(nlevs_max);
    ay.resize(nlevs_max);
    az.resize(nlevs_max);

    z_phys_nd_new.resize(nlevs_max);
    detJ_cc_new.resize(nlevs_max);
    ax_new.resize(nlevs_max);
    ay_new.resize(nlevs_max);
    az_new.resize(nlevs_max);

    z_phys_nd_src.resize(nlevs_max);
    detJ_cc_src.resize(nlevs_max);
    ax_src.resize(nlevs_max);
    ay_src.resize(nlevs_max);
    az_src.resize(nlevs_max);

    z_t_rk.resize(nlevs_max);

    // Mapfactors
    mapfac_m.resize(nlevs_max);
    mapfac_u.resize(nlevs_max);
    mapfac_v.resize(nlevs_max);

    // Thin immersed body
    xflux_imask.resize(nlevs_max);
    yflux_imask.resize(nlevs_max);
    zflux_imask.resize(nlevs_max);
    //overset_imask.resize(nlevs_max);
    thin_xforce.resize(nlevs_max);
    thin_yforce.resize(nlevs_max);
    thin_zforce.resize(nlevs_max);

    // Base state
    base_state.resize(nlevs_max);
    base_state_new.resize(nlevs_max);

    // Wave coupling data
    Hwave.resize(nlevs_max);
    Lwave.resize(nlevs_max);
    for (int lev = 0; lev < max_level; ++lev)
    {
        Hwave[lev] = nullptr;
        Lwave[lev] = nullptr;
    }
    Hwave_onegrid.resize(nlevs_max);
    Lwave_onegrid.resize(nlevs_max);
    for (int lev = 0; lev < max_level; ++lev)
    {
        Hwave_onegrid[lev] = nullptr;
        Lwave_onegrid[lev] = nullptr;
    }
    // Theta prim for MOST
    Theta_prim.resize(nlevs_max);

    // Qv prim for MOST
    Qv_prim.resize(nlevs_max);

    // Qr prim for MOST
    Qr_prim.resize(nlevs_max);

    // Time averaged velocity field
    vel_t_avg.resize(nlevs_max);
    t_avg_cnt.resize(nlevs_max);

#ifdef ERF_USE_NETCDF
    // Size lat long arrays if using netcdf
    lat_m.resize(nlevs_max);
    lon_m.resize(nlevs_max);
    for (int lev = 0; lev < max_level; ++lev)
    {
        lat_m[lev] = nullptr;
        lon_m[lev] = nullptr;
    }
#endif

    // Initialize tagging criteria for mesh refinement
    refinement_criteria_setup();

    for (int lev = 0; lev < max_level; ++lev)
    {
       Print() << "Refinement ratio at level " << lev+1 << " set to be " <<
          ref_ratio[lev][0]  << " " << ref_ratio[lev][1]  <<  " " << ref_ratio[lev][2] << std::endl;
    }

    // We define m_factory even with no EB
    m_factory.resize(max_level+1);

#ifdef AMREX_USE_EB
    // We will create each of these in MakeNewLevel.../RemakeLevel

    // This is needed before initializing level MultiFabs
    MakeEBGeometry();
#endif
}

ERF::~ERF () = default;

// advance solution to final time
void
ERF::Evolve ()
{
    BL_PROFILE_VAR("ERF::Evolve()", evolve);

    Real cur_time = t_new[0];

    // Take one coarse timestep by calling timeStep -- which recursively calls timeStep
    //      for finer levels (with or without subcycling)
    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        Print() << "\nCoarse STEP " << step+1 << " starts ..." << std::endl;

        ComputeDt(step);

        // Make sure we have read enough of the boundary plane data to make it through this timestep
        if (input_bndry_planes)
        {
            m_r2d->read_input_files(cur_time,dt[0],m_bc_extdir_vals);
        }

        int lev = 0;
        int iteration = 1;
        timeStep(lev, cur_time, iteration);

        cur_time  += dt[0];

        Print() << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time
                << " DT = " << dt[0]  << std::endl;

        post_timestep(step, cur_time, dt[0]);

        if (writeNow(cur_time, dt[0], step+1, m_plot_int_1, m_plot_per_1)) {
            last_plot_file_step_1 = step+1;
            WritePlotFile(1,plot_var_names_1);
        }
        if (writeNow(cur_time, dt[0], step+1, m_plot_int_2, m_plot_per_2)) {
            last_plot_file_step_2 = step+1;
            WritePlotFile(2,plot_var_names_2);
        }

        if (writeNow(cur_time, dt[0], step+1, m_check_int, m_check_per)) {
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

    // Write plotfiles at final time
    if ( (m_plot_int_1 > 0 || m_plot_per_1 > 0.) && istep[0] > last_plot_file_step_1 ) {
        WritePlotFile(1,plot_var_names_1);
    }
    if ( (m_plot_int_2 > 0 || m_plot_per_2 > 0.) && istep[0] > last_plot_file_step_2) {
        WritePlotFile(2,plot_var_names_2);
    }

    if ( (m_check_int > 0 || m_check_per > 0.) && istep[0] > last_check_file_step) {
#ifdef ERF_USE_NETCDF
        if (check_type == "netcdf") {
           WriteNCCheckpointFile();
        }
#endif
        if (check_type == "native") {
           WriteCheckpointFile();
        }
    }

    BL_PROFILE_VAR_STOP(evolve);
}

// Called after every coarse timestep
void
ERF::post_timestep (int nstep, Real time, Real dt_lev0)
{
    BL_PROFILE("ERF::post_timestep()");

#ifdef ERF_USE_PARTICLES
    particleData.Redistribute();
#endif

    if (solverChoice.coupling_type == CouplingType::TwoWay)
    {
        bool use_terrain = solverChoice.use_terrain;
        int ncomp = vars_new[0][Vars::cons].nComp();
        for (int lev = finest_level-1; lev >= 0; lev--)
        {
            // The quantity that is conserved is not (rho S), but rather (rho S / m^2) where
            // m is the map scale factor at cell centers
            // Here we pre-divide (rho S) by m^2 before refluxing
            for (MFIter mfi(vars_new[lev][Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                const Array4<      Real>   cons_arr = vars_new[lev][Vars::cons].array(mfi);
                const Array4<const Real> mapfac_arr = mapfac_m[lev]->const_array(mfi);
                if (use_terrain) {
                    const Array4<const Real>   detJ_arr = detJ_cc[lev]->const_array(mfi);
                    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        cons_arr(i,j,k,n) *= detJ_arr(i,j,k) / (mapfac_arr(i,j,0)*mapfac_arr(i,j,0));
                    });
                } else {
                    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        cons_arr(i,j,k,n) /= (mapfac_arr(i,j,0)*mapfac_arr(i,j,0));
                    });
                }
            } // mfi

            // This call refluxes from the lev/lev+1 interface onto lev
            getAdvFluxReg(lev+1)->Reflux(vars_new[lev][Vars::cons], 0, 0, ncomp);

            // Here we multiply (rho S) by m^2 after refluxing
            for (MFIter mfi(vars_new[lev][Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                const Array4<      Real>   cons_arr = vars_new[lev][Vars::cons].array(mfi);
                const Array4<const Real> mapfac_arr = mapfac_m[lev]->const_array(mfi);
                if (use_terrain) {
                    const Array4<const Real>   detJ_arr = detJ_cc[lev]->const_array(mfi);
                    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        cons_arr(i,j,k,n) *= (mapfac_arr(i,j,0)*mapfac_arr(i,j,0)) / detJ_arr(i,j,k);
                    });
                } else {
                    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        cons_arr(i,j,k,n) *= (mapfac_arr(i,j,0)*mapfac_arr(i,j,0));
                    });
                }
            } // mfi

            // We need to do this before anything else because refluxing changes the
            // values of coarse cells underneath fine grids with the assumption they'll
            // be over-written by averaging down
            AverageDownTo(lev,0,ncomp);
        }
    }

    if (is_it_time_for_action(nstep, time, dt_lev0, sum_interval, sum_per)) {
        sum_integrated_quantities(time);
    }

    if (solverChoice.pert_type == PerturbationType::Source ||
        solverChoice.pert_type == PerturbationType::Direct) {
        if (is_it_time_for_action(nstep, time, dt_lev0, pert_interval, -1.)) {
            turbPert.debug(time);
        }
    }

    if (profile_int > 0 && (nstep+1) % profile_int == 0) {
        if (destag_profiles) {
            // all variables cell-centered
            write_1D_profiles(time);
        } else {
            // some variables staggered
            write_1D_profiles_stag(time);
        }
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
      Abort("To output 1D column files ERF must be compiled with NetCDF");
#endif
    }

    if (output_bndry_planes)
    {
      if (is_it_time_for_action(istep[0], time, dt_lev0, bndry_output_planes_interval, bndry_output_planes_per) &&
          time >= bndry_output_planes_start_time)
      {
         bool is_moist = (micro->Get_Qstate_Size() > 0);
         m_w2d->write_planes(istep[0], time, vars_new, is_moist);
      }
    }

    // Write plane/line sampler data
    if (is_it_time_for_action(nstep, time, dt_lev0, sampler_interval, sampler_per) && (data_sampler) ) {
        data_sampler->get_sample_data(geom, vars_new);
        data_sampler->write_sample_data(t_new, istep, ref_ratio, geom);
    }

    // Moving terrain
    if ( solverChoice.use_terrain &&  (solverChoice.terrain_type == TerrainType::Moving) )
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
} // post_timestep

// This is called from main.cpp and handles all initialization, whether from start or restart
void
ERF::InitData ()
{
    BL_PROFILE_VAR("ERF::InitData()", InitData);
    InitData_pre();
    InitData_post();
    BL_PROFILE_VAR_STOP(InitData);
}
// This is called from main.cpp and handles all initialization, whether from start or restart
void
ERF::InitData_pre ()
{
    // Initialize the start time for our CPU-time tracker
    startCPUTime = ParallelDescriptor::second();

    // Create the ReadBndryPlanes object so we can read boundary plane data
    // m_r2d is used by init_bcs so we must instantiate this class before
    if (input_bndry_planes) {
        Print() << "Defining r2d for the first time " << std::endl;
        m_r2d = std::make_unique<ReadBndryPlanes>(geom[0], solverChoice.rdOcp);
    }

    if (!solverChoice.use_terrain && solverChoice.terrain_type != TerrainType::None) {
        Abort("We do not allow terrain_type to be moving or static with use_terrain = false");
    }

    last_plot_file_step_1 = -1;
    last_plot_file_step_2 = -1;
    last_check_file_step  = -1;

    if (restart_chkfile.empty()) {
        // start simulation from the beginning

        const Real time = start_time;
        InitFromScratch(time);
    } else {
        // For initialization this is done in init_only; it is done here for restart
        init_bcs();
    }

    // Verify BCs are compatible with solver choice
    for (int lev(0); lev <= max_level; ++lev) {
        if ( ( (solverChoice.turbChoice[lev].pbl_type == PBLType::MYNN25) ||
               (solverChoice.turbChoice[lev].pbl_type == PBLType::YSU)       ) &&
            phys_bc_type[Orientation(Direction::z,Orientation::low)] != ERF_BC::MOST ) {
            Abort("MYNN2.5/YSU PBL Model requires MOST at lower boundary");
        }
    }
}

void
ERF::InitData_post ()
{
    if (restart_chkfile.empty()) {
        if (solverChoice.use_terrain) {
            if (init_type == InitType::Ideal) {
                Abort("We do not currently support init_type = ideal with terrain");
            }
        }

        //
        // Make sure that detJ and z_phys_cc are the average of the data on a finer level if there is one
        //
        if (solverChoice.use_terrain != 0) {
            for (int crse_lev = finest_level-1; crse_lev >= 0; crse_lev--) {
                average_down(  *detJ_cc[crse_lev+1],   *detJ_cc[crse_lev], 0, 1, refRatio(crse_lev));
                average_down(*z_phys_cc[crse_lev+1], *z_phys_cc[crse_lev], 0, 1, refRatio(crse_lev));
            }
        }

        // If using the Deardoff LES model,
        // we initialize rho_KE to be nonzero (and positive) so that we end up
        // with reasonable values for the eddy diffusivity and the MOST fluxes
        // (~ 1/diffusivity) do not blow up
        Real RhoKE_0;
        ParmParse pp(pp_prefix);
        if (pp.query("RhoKE_0", RhoKE_0)) {
            // Uniform initial rho*e field
            for (int lev = 0; lev <= finest_level; lev++) {
                if (solverChoice.turbChoice[lev].les_type == LESType::Deardorff) {
                    Print() << "Initializing uniform rhoKE=" << RhoKE_0
                        << " on level " << lev
                        << std::endl;
                    vars_new[lev][Vars::cons].setVal(RhoKE_0,RhoKE_comp,1,0);
                } else {
                    vars_new[lev][Vars::cons].setVal(0.0,RhoKE_comp,1,0);
                }
            }
        }

        Real KE_0;
        if (pp.query("KE_0", KE_0)) {
            // Uniform initial e field
            for (int lev = 0; lev <= finest_level; lev++) {
                auto& lev_new = vars_new[lev];
                if (solverChoice.turbChoice[lev].les_type == LESType::Deardorff) {
                    Print() << "Initializing uniform KE=" << KE_0
                        << " on level " << lev
                        << std::endl;
                    for (MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                        const Box &bx = mfi.tilebox();
                        const auto &cons_arr = lev_new[Vars::cons].array(mfi);
                        // We want to set the lateral BC values, too
                        Box gbx = bx; // Copy constructor
                        gbx.grow(0,1); gbx.grow(1,1); // Grow by one in the lateral directions
                        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            cons_arr(i,j,k,RhoKE_comp) = cons_arr(i,j,k,Rho_comp) * KE_0;
                        });
                    } // mfi
                } else {
                    lev_new[Vars::cons].setVal(0.0,RhoKE_comp,1,0);
                }
            } // lev
        }

        Real QKE_0;
        if (pp.query("QKE_0", QKE_0)) {
            Print() << "Initializing uniform QKE=" << QKE_0 << std::endl;
            for (int lev = 0; lev <= finest_level; lev++) {
                auto& lev_new = vars_new[lev];
                for (MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box &bx = mfi.tilebox();
                    const auto &cons_arr = lev_new[Vars::cons].array(mfi);
                    // We want to set the lateral BC values, too
                    Box gbx = bx; // Copy constructor
                    gbx.grow(0,1); gbx.grow(1,1); // Grow by one in the lateral directions
                    ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        cons_arr(i,j,k,RhoQKE_comp) = cons_arr(i,j,k,Rho_comp) * QKE_0;
                    });
                } // mfi
            }
        }

        if (solverChoice.coupling_type == CouplingType::TwoWay) {
            AverageDown();
        }

        if ((solverChoice.advChoice.zero_xflux.size() > 0) ||
            (solverChoice.advChoice.zero_yflux.size() > 0) ||
            (solverChoice.advChoice.zero_zflux.size() > 0))
        {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(finest_level == 0,
                "Thin immersed body with refinement not currently supported.");
            if (solverChoice.use_terrain == 1) {
                amrex::Print() << "NOTE: Thin immersed body with terrain has not been tested." << std::endl;
            }
        }

#ifdef ERF_USE_PARTICLES
        if (Microphysics::modelType(solverChoice.moisture_type) == MoistureModelType::Lagrangian) {
            for (int lev = 0; lev <= finest_level; lev++) {
                dynamic_cast<LagrangianMicrophysics&>(*micro).initParticles(z_phys_nd[lev]);
            }
        }
#endif

    } else { // Restart from a checkpoint

        restart();

        // Create the physbc objects for {cons, u, v, w, base state}
        for (int lev(0); lev <= max_level; ++lev) {
            make_physbcs(lev);
        }
    }

#ifdef ERF_USE_PARTICLES
    /* If using a Lagrangian microphysics model, its particle container has now been
       constructed and initialized (calls to micro->Init). So, add its pointer to
       ERF::particleData and remove its name from list of unallocated particle containers. */
    if (Microphysics::modelType(solverChoice.moisture_type) == MoistureModelType::Lagrangian) {
        const auto& pc_name( dynamic_cast<LagrangianMicrophysics&>(*micro).getName() );
        const auto& pc_ptr( dynamic_cast<LagrangianMicrophysics&>(*micro).getParticleContainer() );
        particleData.pushBack(pc_name, pc_ptr);
        particleData.getNamesUnalloc().remove(pc_name);
    }
#endif

    if (input_bndry_planes) {
        // Read the "time.dat" file to know what data is available
        m_r2d->read_time_file();

        // We haven't populated dt yet, set to 0 to ensure assert doesn't crash
        Real dt_dummy = 0.0;
        m_r2d->read_input_files(t_new[0],dt_dummy,m_bc_extdir_vals);
    }

    if (solverChoice.custom_rhotheta_forcing)
    {
        h_rhotheta_src.resize(max_level+1, Vector<Real>(0));
        d_rhotheta_src.resize(max_level+1, Gpu::DeviceVector<Real>(0));
        for (int lev = 0; lev <= finest_level; lev++) {
            const int domlen = geom[lev].Domain().length(2);
            h_rhotheta_src[lev].resize(domlen, 0.0_rt);
            d_rhotheta_src[lev].resize(domlen, 0.0_rt);
            prob->update_rhotheta_sources(t_new[0],
                                          h_rhotheta_src[lev], d_rhotheta_src[lev],
                                          geom[lev], z_phys_cc[lev]);
        }
    }

    if (solverChoice.have_geo_wind_profile)
    {
        h_u_geos.resize(max_level+1, Vector<Real>(0));
        d_u_geos.resize(max_level+1, Gpu::DeviceVector<Real>(0));
        h_v_geos.resize(max_level+1, Vector<Real>(0));
        d_v_geos.resize(max_level+1, Gpu::DeviceVector<Real>(0));
        for (int lev = 0; lev <= finest_level; lev++) {
            const int domlen = geom[lev].Domain().length(2);
            h_u_geos[lev].resize(domlen, 0.0_rt);
            d_u_geos[lev].resize(domlen, 0.0_rt);
            h_v_geos[lev].resize(domlen, 0.0_rt);
            d_v_geos[lev].resize(domlen, 0.0_rt);
            if (solverChoice.custom_geostrophic_profile) {
                prob->update_geostrophic_profile(t_new[0],
                                                 h_u_geos[lev], d_u_geos[lev],
                                                 h_v_geos[lev], d_v_geos[lev],
                                                 geom[lev], z_phys_cc[lev]);
            } else {
                if (solverChoice.use_terrain > 0) {
                    amrex::Print() << "Note: 1-D geostrophic wind profile input is only defined for grid stretching, not real terrain" << std::endl;
                }
                init_geo_wind_profile(solverChoice.abl_geo_wind_table,
                                      h_u_geos[lev], d_u_geos[lev],
                                      h_v_geos[lev], d_v_geos[lev],
                                      geom[lev],
                                      zlevels_stag[0]);
            }
        }
    }

    if (solverChoice.custom_moisture_forcing)
    {
        h_rhoqt_src.resize(max_level+1, Vector<Real>(0));
        d_rhoqt_src.resize(max_level+1, Gpu::DeviceVector<Real>(0));
        for (int lev = 0; lev <= finest_level; lev++) {
            const int domlen = geom[lev].Domain().length(2);
            h_rhoqt_src[lev].resize(domlen, 0.0_rt);
            d_rhoqt_src[lev].resize(domlen, 0.0_rt);
            prob->update_rhoqt_sources(t_new[0],
                                       h_rhoqt_src[lev], d_rhoqt_src[lev],
                                       geom[lev], z_phys_cc[lev]);
        }
    }

    if (solverChoice.custom_w_subsidence)
    {
        h_w_subsid.resize(max_level+1, Vector<Real>(0));
        d_w_subsid.resize(max_level+1, Gpu::DeviceVector<Real>(0));
        for (int lev = 0; lev <= finest_level; lev++) {
            const int domlen = geom[lev].Domain().length(2) + 1; // lives on z-faces
            h_w_subsid[lev].resize(domlen, 0.0_rt);
            d_w_subsid[lev].resize(domlen, 0.0_rt);
            prob->update_w_subsidence(t_new[0],
                                      h_w_subsid[lev], d_w_subsid[lev],
                                      geom[lev], z_phys_cc[lev]);
        }
    }

    if (solverChoice.rayleigh_damp_U ||solverChoice.rayleigh_damp_V ||
        solverChoice.rayleigh_damp_W ||solverChoice.rayleigh_damp_T)
    {
        initRayleigh();
        if (init_type == InitType::Input_Sounding)
        {
            // Overwrite ubar, vbar, and thetabar with input profiles;
            // wbar is assumed to be 0. Note: the tau coefficient set by
            // prob->erf_init_rayleigh() is still used
            bool restarting = (!restart_chkfile.empty());
            setRayleighRefFromSounding(restarting);
        }
    }

    // Read in sponge data from input file
    if(solverChoice.spongeChoice.sponge_type == "input_sponge")
    {
        initSponge();
        bool restarting = (!restart_chkfile.empty());
        setSpongeRefFromSounding(restarting);
    }

    if (is_it_time_for_action(istep[0], t_new[0], dt[0], sum_interval, sum_per)) {
        sum_integrated_quantities(t_new[0]);
    }

    if (solverChoice.pert_type == PerturbationType::Source ||
        solverChoice.pert_type == PerturbationType::Direct) {
        if (is_it_time_for_action(istep[0], t_new[0], dt[0], pert_interval, -1.)) {
            turbPert.debug(t_new[0]);
        }
    }

    // We only write the file at level 0 for now
    if (output_bndry_planes)
    {
        // Create the WriteBndryPlanes object so we can handle writing of boundary plane data
        m_w2d = std::make_unique<WriteBndryPlanes>(grids,geom);

        Real time = 0.;
        if (time >= bndry_output_planes_start_time) {
            bool is_moist = (micro->Get_Qstate_Size() > 0);
            m_w2d->write_planes(0, time, vars_new, is_moist);
        }
    }

#ifdef ERF_USE_POISSON_SOLVE
    if (restart_chkfile == "")
    {
        // Note -- this projection is only defined for no terrain
        if (solverChoice.project_initial_velocity) {
            AMREX_ALWAYS_ASSERT(solverChoice.use_terrain == 0);
            Real dummy_dt = 1.0;
            for (int lev = 0; lev <= finest_level; ++lev)
            {
                project_velocities(lev, dummy_dt, vars_new[lev], pp_inc[lev]);
                pp_inc[lev].setVal(0.);
            }
        }
    }
#endif

    // Copy from new into old just in case
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        auto& lev_new = vars_new[lev];
        auto& lev_old = vars_old[lev];

        int ncomp = lev_new[Vars::cons].nComp();

        MultiFab::Copy(lev_old[Vars::cons],lev_new[Vars::cons],0,0,ncomp,lev_new[Vars::cons].nGrowVect());
        MultiFab::Copy(lev_old[Vars::xvel],lev_new[Vars::xvel],0,0,    1,lev_new[Vars::xvel].nGrowVect());
        MultiFab::Copy(lev_old[Vars::yvel],lev_new[Vars::yvel],0,0,    1,lev_new[Vars::yvel].nGrowVect());
        MultiFab::Copy(lev_old[Vars::zvel],lev_new[Vars::zvel],0,0,    1,lev_new[Vars::zvel].nGrowVect());
    }

    // Compute the minimum dz in the domain at each level (to be used for setting the timestep)
    dz_min.resize(max_level+1);
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        dz_min[lev] = geom[lev].CellSize(2);
        if ( solverChoice.use_terrain ) {
            dz_min[lev] *= (*detJ_cc[lev]).min(0);
        }
    }

    ComputeDt();

    // Fill ghost cells/faces
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        if (lev > 0 && cf_width >= 0) {
            Construct_ERFFillPatchers(lev);
        }

        auto& lev_new = vars_new[lev];

        //
        // Fill boundary conditions -- not sure why we need this here
        //
        bool fillset = false;
        FillPatch(lev, t_new[lev],
                  {&lev_new[Vars::cons],&lev_new[Vars::xvel],&lev_new[Vars::yvel],&lev_new[Vars::zvel]},
                  {&lev_new[Vars::cons],&rU_new[lev],&rV_new[lev],&rW_new[lev]},
                  fillset);

        //
        // We do this here to make sure level (lev-1) boundary conditions are filled
        // before we interpolate to level (lev) ghost cells
        //
        if (lev < finest_level) {
            auto& lev_old = vars_old[lev];
            MultiFab::Copy(lev_old[Vars::cons],lev_new[Vars::cons],0,0,lev_old[Vars::cons].nComp(),lev_old[Vars::cons].nGrowVect());
            MultiFab::Copy(lev_old[Vars::xvel],lev_new[Vars::xvel],0,0,lev_old[Vars::xvel].nComp(),lev_old[Vars::xvel].nGrowVect());
            MultiFab::Copy(lev_old[Vars::yvel],lev_new[Vars::yvel],0,0,lev_old[Vars::yvel].nComp(),lev_old[Vars::yvel].nGrowVect());
            MultiFab::Copy(lev_old[Vars::zvel],lev_new[Vars::zvel],0,0,lev_old[Vars::zvel].nComp(),lev_old[Vars::zvel].nGrowVect());
        }

        //
        // We fill the ghost cell values of the base state in case it wasn't done in the initialization
        //
        base_state[lev].FillBoundary(geom[lev].periodicity());

        // For moving terrain only
        if (solverChoice.terrain_type != TerrainType::Static) {
            MultiFab::Copy(base_state_new[lev],base_state[lev],0,0,3,1);
            base_state_new[lev].FillBoundary(geom[lev].periodicity());
        }

    }

    // Allow idealized cases over water, used to set lmask
    ParmParse pp("erf");
    int is_land;
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        if (pp.query("is_land", is_land, lev)) {
            if (is_land == 1) {
                amrex::Print() << "Level " << lev << " is land" << std::endl;
            } else if (is_land == 0) {
                amrex::Print() << "Level " << lev << " is water" << std::endl;
            } else {
                Error("is_land should be 0 or 1");
            }
            lmask_lev[lev][0]->setVal(is_land);
            lmask_lev[lev][0]->FillBoundary(geom[lev].periodicity());
        }
    }

#ifdef ERF_USE_WW3_COUPLING
    int lev = 0;
    amrex::Print() <<  " About to call send_to_ww3 from ERF.cpp" << std::endl;
    send_to_ww3(lev);
    amrex::Print() <<  " About to call read_waves from ERF.cpp"  << std::endl;
    read_waves(lev);
   // send_to_ww3(lev);
#endif

    // Configure ABLMost params if used MostWall boundary condition
    // NOTE: we must set up the MOST routine after calling FillPatch
    //       in order to have lateral ghost cells filled (MOST + terrain interp).
    //       FillPatch does not call MOST, FillIntermediatePatch does.
    if (phys_bc_type[Orientation(Direction::z,Orientation::low)] == ERF_BC::MOST)
    {
        bool use_exp_most = solverChoice.use_explicit_most;
        bool use_rot_most = solverChoice.use_rotate_most;
        if (use_exp_most) {
            Print() << "Using MOST with explicitly included surface stresses" << std::endl;
            if (use_rot_most) {
                Print() << "Using MOST with surface stress rotations" << std::endl;
            }
        }

        m_most = std::make_unique<ABLMost>(geom, use_exp_most, use_rot_most,
                                           vars_old, Theta_prim, Qv_prim, Qr_prim, z_phys_nd,
                                           sst_lev, lmask_lev, lsm_data, lsm_flux,
                                           Hwave, Lwave, eddyDiffs_lev
#ifdef ERF_USE_NETCDF
                                           ,start_bdy_time, bdy_time_interval
#endif
                                           );


        if (restart_chkfile != "") {
            // Update surface fields if needed
            ReadCheckpointFileMOST();
        }

        // We now configure ABLMost params here so that we can print the averages at t=0
        // Note we don't fill ghost cells here because this is just for diagnostics
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            Real time  = t_new[lev];
            IntVect ng = Theta_prim[lev]->nGrowVect();

            MultiFab::Copy(  *Theta_prim[lev], vars_new[lev][Vars::cons], RhoTheta_comp, 0, 1, ng);
            MultiFab::Divide(*Theta_prim[lev], vars_new[lev][Vars::cons],      Rho_comp, 0, 1, ng);

            if (solverChoice.moisture_type != MoistureType::None) {
                ng = Qv_prim[lev]->nGrowVect();

                MultiFab::Copy(  *Qv_prim[lev], vars_new[lev][Vars::cons], RhoQ1_comp, 0, 1, ng);
                MultiFab::Divide(*Qv_prim[lev], vars_new[lev][Vars::cons],   Rho_comp, 0, 1, ng);

                int rhoqr_comp = solverChoice.RhoQr_comp;
                if (rhoqr_comp > -1) {
                    MultiFab::Copy(  *Qr_prim[lev], vars_new[lev][Vars::cons], rhoqr_comp, 0, 1, ng);
                    MultiFab::Divide(*Qr_prim[lev], vars_new[lev][Vars::cons],   Rho_comp, 0, 1, ng);
                } else {
                    Qr_prim[lev]->setVal(0.0);
                }
            }
            m_most->update_mac_ptrs(lev, vars_new, Theta_prim, Qv_prim, Qr_prim);

            if (restart_chkfile == "") {
                // Only do this if starting from scratch; if restarting, then
                // we don't want to call update_fluxes multiple times because
                // it will change u* and theta* from their previous values
                m_most->update_pblh(lev, vars_new, z_phys_cc[lev].get(),
                                    solverChoice.RhoQv_comp, solverChoice.RhoQr_comp);
                m_most->update_fluxes(lev, time);
            }
        }
    }

    // Update micro vars before first plot file
    if (solverChoice.moisture_type != MoistureType::None) {
        for (int lev = 0; lev <= finest_level; ++lev) micro->Update_Micro_Vars_Lev(lev, vars_new[lev][Vars::cons]);
    }

    // Fill time averaged velocities before first plot file
    if (solverChoice.time_avg_vel) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            Time_Avg_Vel_atCC(dt[lev], t_avg_cnt[lev], vel_t_avg[lev].get(),
                              vars_new[lev][Vars::xvel],
                              vars_new[lev][Vars::yvel],
                              vars_new[lev][Vars::zvel]);
        }
    }

    // check for additional plotting variables that are available after particle containers
    // are setup.
    const std::string& pv1 = "plot_vars_1"; appendPlotVariables(pv1,plot_var_names_1);
    const std::string& pv2 = "plot_vars_2"; appendPlotVariables(pv2,plot_var_names_2);

    if ( restart_chkfile.empty() && (m_check_int > 0 || m_check_per > 0.) )
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

    if ( (restart_chkfile.empty()) ||
         (!restart_chkfile.empty() && plot_file_on_restart) )
    {
        if (m_plot_int_1 > 0 || m_plot_per_1 > 0.)
        {
            WritePlotFile(1,plot_var_names_1);
            last_plot_file_step_1 = istep[0];
        }
        if (m_plot_int_2 > 0 || m_plot_per_2 > 0.)
        {
            WritePlotFile(2,plot_var_names_2);
            last_plot_file_step_2 = istep[0];
        }
    }

    // Set these up here because we need to know which MPI rank "cell" is on...
    if (pp.contains("data_log"))
    {
        int num_datalogs = pp.countval("data_log");
        datalog.resize(num_datalogs);
        datalogname.resize(num_datalogs);
        pp.queryarr("data_log",datalogname,0,num_datalogs);
        for (int i = 0; i < num_datalogs; i++)
            setRecordDataInfo(i,datalogname[i]);
    }

    if (restart_chkfile.empty() && profile_int > 0) {
        if (destag_profiles) {
            // all variables cell-centered
            write_1D_profiles(t_new[0]);
        } else {
            // some variables staggered
            write_1D_profiles_stag(t_new[0]);
        }
    }

    if (pp.contains("sample_point_log") && pp.contains("sample_point"))
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

        int num_sampleptlogs = pp.countval("sample_point_log");
        AMREX_ALWAYS_ASSERT(num_sampleptlogs == num_samplepts);
        if (num_sampleptlogs > 0) {
            sampleptlog.resize(num_sampleptlogs);
            sampleptlogname.resize(num_sampleptlogs);
            pp.queryarr("sample_point_log",sampleptlogname,0,num_sampleptlogs);

            for (int i = 0; i < num_sampleptlogs; i++) {
                setRecordSamplePointInfo(i,lev,samplepoint[i],sampleptlogname[i]);
            }
        }

    }

    if (pp.contains("sample_line_log") && pp.contains("sample_line"))
    {
        int lev = 0;

        int num_samplelines = pp.countval("sample_line") / AMREX_SPACEDIM;
        if (num_samplelines > 0) {
            Vector<int> index; index.resize(num_samplelines*AMREX_SPACEDIM);
            sampleline.resize(num_samplelines);

            pp.queryarr("sample_line",index,0,num_samplelines*AMREX_SPACEDIM);
            for (int i = 0; i < num_samplelines; i++) {
                IntVect iv(index[AMREX_SPACEDIM*i+0],index[AMREX_SPACEDIM*i+1],index[AMREX_SPACEDIM*i+2]);
                sampleline[i] = iv;
            }
        }

        int num_samplelinelogs = pp.countval("sample_line_log");
        AMREX_ALWAYS_ASSERT(num_samplelinelogs == num_samplelines);
        if (num_samplelinelogs > 0) {
            samplelinelog.resize(num_samplelinelogs);
            samplelinelogname.resize(num_samplelinelogs);
            pp.queryarr("sample_line_log",samplelinelogname,0,num_samplelinelogs);

            for (int i = 0; i < num_samplelinelogs; i++) {
                setRecordSampleLineInfo(i,lev,sampleline[i],samplelinelogname[i]);
            }
        }

    }

    // Create object to do line and plane sampling if needed
    bool do_line = false; bool do_plane = false;
    pp.query("do_line_sampling",do_line); pp.query("do_plane_sampling",do_plane);
    if (do_line || do_plane) { data_sampler = std::make_unique<SampleData>(do_line, do_plane); }

#ifdef ERF_USE_EB
    bool write_eb_surface = false;
    pp.query("write_eb_surface", write_eb_surface);
    if (write_eb_surface) WriteMyEBSurface();
#endif

}

// Initialize microphysics object
void
ERF::initializeMicrophysics (const int& a_nlevsmax /*!< number of AMR levels */)
{
    if (Microphysics::modelType(solverChoice.moisture_type) == MoistureModelType::Eulerian) {

        micro = std::make_unique<EulerianMicrophysics>(a_nlevsmax, solverChoice.moisture_type);

    } else if (Microphysics::modelType(solverChoice.moisture_type) == MoistureModelType::Lagrangian) {
#ifdef ERF_USE_PARTICLES

        micro = std::make_unique<LagrangianMicrophysics>(a_nlevsmax, solverChoice.moisture_type);
        /* Lagrangian microphysics models will have a particle container; it needs to be added
           to ERF::particleData */
        const auto& pc_name( dynamic_cast<LagrangianMicrophysics&>(*micro).getName() );
        /* The particle container has not yet been constructed and initialized, so just add
           its name here for now (so that functions to set plotting variables can see it). */
        particleData.addName( pc_name );

#else
        Abort("Lagrangian microphysics can be used when compiled with ERF_USE_PARTICLES");
#endif
    }

    qmoist.resize(a_nlevsmax);
    return;
}


#ifdef ERF_USE_WINDFARM
void
ERF::initializeWindFarm(const int& a_nlevsmax/*!< number of AMR levels */ )
{
    windfarm = std::make_unique<WindFarm>(a_nlevsmax, solverChoice.windfarm_type);
}
#endif

void
ERF::restart ()
{
    // TODO: This could be deleted since ba/dm are not created yet?
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        auto& lev_new = vars_new[lev];
        auto& lev_old = vars_old[lev];
        lev_new[Vars::cons].setVal(0.); lev_old[Vars::cons].setVal(0.);
        lev_new[Vars::xvel].setVal(0.); lev_old[Vars::xvel].setVal(0.);
        lev_new[Vars::yvel].setVal(0.); lev_old[Vars::yvel].setVal(0.);
        lev_new[Vars::zvel].setVal(0.); lev_old[Vars::zvel].setVal(0.);
    }

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

// This is called only if starting from scratch (from ERF::MakeNewLevelFromScratch)
//
// If we are restarting, the base state is read from the restart file, including
// ghost cell data.
void
ERF::init_only (int lev, Real time)
{
    // Map the words in the inputs file to BC types, then translate
    //     those types into what they mean for each variable
    // This must be called before initHSE (where the base state is initialized)
    if (lev == 0 && init_type != InitType::Ideal) {
        init_bcs();
    }

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    auto& lev_new = vars_new[lev];
    auto& lev_old = vars_old[lev];

    // Loop over grids at this level to initialize our grid data
    lev_new[Vars::cons].setVal(0.0); lev_old[Vars::cons].setVal(0.0);
    lev_new[Vars::xvel].setVal(0.0); lev_old[Vars::xvel].setVal(0.0);
    lev_new[Vars::yvel].setVal(0.0); lev_old[Vars::yvel].setVal(0.0);
    lev_new[Vars::zvel].setVal(0.0); lev_old[Vars::zvel].setVal(0.0);

    // Initialize background flow (optional)
    if (init_type == InitType::Input_Sounding) {
        // The base state is initialized by integrating vertically through the
        // input sounding, if the init_sounding_ideal flag is set; otherwise
        // it is set by initHSE()

        // The physbc's need the terrain but are needed for initHSE
        // We have already made the terrain in the call to init_zphys
        //    in MakeNewLevelFromScratch
        make_physbcs(lev);

        // Now init the base state and the data itself
        init_from_input_sounding(lev);

        if (init_sounding_ideal) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(solverChoice.use_gravity,
                "Gravity should be on to be consistent with sounding initialization.");
        } else {
            initHSE();
        }

#ifdef ERF_USE_NETCDF
    } else if (init_type == InitType::Ideal || init_type == InitType::Real) {
        // The base state is initialized from WRF wrfinput data, output by
        // ideal.exe or real.exe
        init_from_wrfinput(lev);

        // The physbc's need the terrain but are needed for initHSE
        if (init_type == InitType::Ideal) {
            make_physbcs(lev);
            initHSE(lev);
        }

    } else if (init_type == InitType::Metgrid) {
        // The base state is initialized from data output by WPS metgrid;
        // we will rebalance after interpolation
        init_from_metgrid(lev);
#endif
    } else if (init_type == InitType::Uniform) {
        // Initialize a uniform background field and base state based on the
        // problem-specified reference density and temperature

        // The physbc's need the terrain but are needed for initHSE
        make_physbcs(lev);

        init_uniform(lev);
        initHSE(lev);
    } else {
        // No background flow initialization specified, initialize the
        // background field to be equal to the base state, calculated from the
        // problem-specific erf_init_dens_hse

        // The bc's need the terrain but are needed for initHSE
        make_physbcs(lev);

        // We will initialize the state from the background state so must set that first
        initHSE(lev);
        init_from_hse(lev);
    }

    // Add problem-specific flow features
    //
    // Notes:
    // - This calls init_custom_pert that is defined for each problem
    // - This may modify the base state
    // - The fields set by init_custom_pert are **perturbations** to the
    //   background flow set based on init_type
    init_custom(lev);

    // Ensure that the face-based data are the same on both sides of a periodic domain.
    // The data associated with the lower grid ID is considered the correct value.
    lev_new[Vars::xvel].OverrideSync(geom[lev].periodicity());
    lev_new[Vars::yvel].OverrideSync(geom[lev].periodicity());
    lev_new[Vars::zvel].OverrideSync(geom[lev].periodicity());

   if(solverChoice.spongeChoice.sponge_type == "input_sponge"){
        input_sponge(lev);
   }

    // Initialize turbulent perturbation
    if (solverChoice.pert_type == PerturbationType::Source ||
        solverChoice.pert_type == PerturbationType::Direct) {
        if (lev == 0) {
            turbPert_update(lev, 0.);
            turbPert_amplitude(lev);
        }
    }
}

// Read in some parameters from inputs file
void
ERF::ReadParameters ()
{
    {
        ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);

        pp.query("start_time", start_time); // This is optional, it defaults to 0
    }

    ParmParse pp(pp_prefix);
    ParmParse pp_amr("amr");
    {
        // The type of the file we restart from
        pp.query("restart_type", restart_type);

        pp.query("regrid_int", regrid_int);
        pp.query("check_file", check_file);
        pp.query("check_type", check_type);

        // The regression tests use "amr.restart" and "amr.m_check_int" so we allow
        //    for those or "erf.restart" / "erf.m_check_int" with the former taking
        //    precedence if both are specified
        pp.query("check_int", m_check_int);
        pp.query("check_per", m_check_per);
        pp_amr.query("check_int", m_check_int);
        pp_amr.query("check_per", m_check_per);

        pp.query("restart", restart_chkfile);
        pp_amr.query("restart", restart_chkfile);

        // Verbosity
        pp.query("v", verbose);
        pp.query("mg_v", mg_verbose);
        pp.query("use_fft", use_fft);

        // Frequency of diagnostic output
        pp.query("sum_interval", sum_interval);
        pp.query("sum_period"  , sum_per);

        pp.query("pert_interval", pert_interval);

        // Time step controls
        pp.query("cfl", cfl);
        pp.query("init_shrink", init_shrink);
        pp.query("change_max", change_max);

        fixed_dt.resize(max_level+1,-1.);
        fixed_fast_dt.resize(max_level+1,-1.);

        pp.query("fixed_dt", fixed_dt[0]);
        pp.query("fixed_fast_dt", fixed_fast_dt[0]);

        for (int lev = 1; lev <= max_level; lev++)
        {
            fixed_dt[lev]      = fixed_dt[lev-1]     / static_cast<Real>(MaxRefRatio(lev-1));
            fixed_fast_dt[lev] = fixed_fast_dt[lev-1] / static_cast<Real>(MaxRefRatio(lev-1));
        }

        pp.query("fixed_mri_dt_ratio", fixed_mri_dt_ratio);

        // How to initialize
        init_type = InitType::None;
        pp.query_enum_case_insensitive("init_type",init_type);

        // Should we use the bcs we've read in from wrfbdy or metgrid files?
        // We default to yes if we have them, but the user can override that option
        use_real_bcs = ( (init_type == InitType::Real) || (init_type == InitType::Metgrid) );
        pp.query("use_real_bcs",use_real_bcs);

        // We use this to keep track of how many boxes we read in from WRF initialization
        num_files_at_level.resize(max_level+1,0);

        // We use this to keep track of how many boxes are specified thru the refinement indicators
        num_boxes_at_level.resize(max_level+1,0);
            boxes_at_level.resize(max_level+1);

        // We always have exactly one file at level 0
        num_boxes_at_level[0] = 1;
        boxes_at_level[0].resize(1);
        boxes_at_level[0][0] = geom[0].Domain();

#ifdef ERF_USE_NETCDF
        nc_init_file.resize(max_level+1);

        // NetCDF wrfinput initialization files -- possibly multiple files at each of multiple levels
        //        but we always have exactly one file at level 0
        for (int lev = 0; lev <= max_level; lev++)
        {
            const std::string nc_file_names = Concatenate("nc_init_file_",lev,1);
            if (pp.contains(nc_file_names.c_str()))
            {
                int num_files = pp.countval(nc_file_names.c_str());
                num_files_at_level[lev] = num_files;
                nc_init_file[lev].resize(num_files);
                pp.queryarr(nc_file_names.c_str(), nc_init_file[lev],0,num_files);
                for (int j = 0; j < num_files; j++)
                    Print() << "Reading NC init file names at level " << lev << " and index " << j << " : " << nc_init_file[lev][j] << std::endl;
            }
        }

        // NetCDF wrfbdy lateral boundary file
        pp.query("nc_bdy_file", nc_bdy_file);
#endif

        // Flag to trigger initialization from input_sounding like WRF's ideal.exe
        pp.query("init_sounding_ideal", init_sounding_ideal);

        // Output format
        pp.query("plotfile_type", plotfile_type);
        pp.query("plot_file_1",   plot_file_1);
        pp.query("plot_file_2",   plot_file_2);
        pp.query("plot_int_1" , m_plot_int_1);
        pp.query("plot_int_2" , m_plot_int_2);
        pp.query("plot_per_1",  m_plot_per_1);
        pp.query("plot_per_2",  m_plot_per_2);

        pp.query("expand_plotvars_to_unif_rr",m_expand_plotvars_to_unif_rr);

        if ( (m_plot_int_1 > 0 && m_plot_per_1 > 0) ||
             (m_plot_int_2 > 0 && m_plot_per_2 > 0.) ) {
            Abort("Must choose only one of plot_int or plot_per");
        }

        pp.query("profile_int", profile_int);
        pp.query("destag_profiles", destag_profiles);

        pp.query("plot_lsm", plot_lsm);
#ifdef ERF_USE_RRTMGP
        pp.query("plot_rad", plot_rad);
#endif

        pp.query("output_1d_column", output_1d_column);
        pp.query("column_per", column_per);
        pp.query("column_interval", column_interval);
        pp.query("column_loc_x", column_loc_x);
        pp.query("column_loc_y", column_loc_y);
        pp.query("column_file_name", column_file_name);

        // Sampler output frequency
        pp.query("sampler_per", sampler_per);
        pp.query("sampler_interval", sampler_interval);

        // Specify information about outputting planes of data
        pp.query("output_bndry_planes", output_bndry_planes);
        pp.query("bndry_output_planes_interval", bndry_output_planes_interval);
        pp.query("bndry_output_planes_per", bndry_output_planes_per);
        pp.query("bndry_output_start_time", bndry_output_planes_start_time);

        // Specify whether ingest boundary planes of data
        pp.query("input_bndry_planes", input_bndry_planes);

        // Query the set and total widths for wrfbdy interior ghost cells
        pp.query("real_width", real_width);
        pp.query("real_set_width", real_set_width);

        // Query the set and total widths for crse-fine interior ghost cells
        pp.query("cf_width", cf_width);
        pp.query("cf_set_width", cf_set_width);

        // AmrMesh iterate on grids?
        bool iterate(true);
        pp_amr.query("iterate_grids",iterate);
        if (!iterate) SetIterateToFalse();
    }

#ifdef ERF_USE_PARTICLES
    readTracersParams();
#endif

#ifdef ERF_USE_MULTIBLOCK
    solverChoice.pp_prefix = pp_prefix;
#endif

    solverChoice.init_params(max_level);

    // No moving terrain with init real (we must do this after init_params
    //    because that is where we set terrain_type
    if (init_type == InitType::Real && solverChoice.terrain_type == TerrainType::Moving) {
        Abort("Moving terrain is not supported with init real");
    }

    // What type of land surface model to use
    // NOTE: Must be checked after init_params
    if (solverChoice.lsm_type == LandSurfaceType::SLM) {
        lsm.SetModel<SLM>();
        Print() << "SLM land surface model!\n";
    } else if (solverChoice.lsm_type == LandSurfaceType::MM5) {
        lsm.SetModel<MM5>();
        Print() << "MM5 land surface model!\n";
    } else if (solverChoice.lsm_type == LandSurfaceType::None) {
        lsm.SetModel<NullSurf>();
        Print() << "Null land surface model!\n";
    } else {
        Abort("Dont know this LandSurfaceType!") ;
    }

    if (verbose > 0) {
        solverChoice.display(max_level);
    }

    ParameterSanityChecks();
}

// Read in some parameters from inputs file
void
ERF::ParameterSanityChecks ()
{
    AMREX_ALWAYS_ASSERT(cfl > 0. || fixed_dt[0] > 0.);

    // We don't allow use_real_bcs to be true if init_type is not either InitType::Rreal or InitType::Metgrid
    AMREX_ALWAYS_ASSERT(!use_real_bcs || ((init_type == InitType::Real) || (init_type == InitType::Metgrid)) );

    AMREX_ALWAYS_ASSERT(real_width >= 0);
    AMREX_ALWAYS_ASSERT(real_set_width >= 0);
    AMREX_ALWAYS_ASSERT(real_width >= real_set_width);

    if (cf_width < 0 || cf_set_width < 0 || cf_width < cf_set_width) {
        Abort("You must set cf_width >= cf_set_width >= 0");
    }
    if (max_level > 0 && cf_set_width > 0) {
        for (int lev = 1; lev <= max_level; lev++) {
            if (cf_set_width%ref_ratio[lev-1][0] != 0 ||
                cf_set_width%ref_ratio[lev-1][1] != 0 ||
                cf_set_width%ref_ratio[lev-1][2] != 0 ) {
                Abort("You must set cf_width to be a multiple of ref_ratio");
            }
        }
    }

    if (plotfile_type != "amrex" &&
        plotfile_type != "netcdf" && plotfile_type != "NetCDF" &&
        plotfile_type != "hdf5"   && plotfile_type != "HDF5" )
    {
        Print() << "User selected plotfile_type = " << plotfile_type << std::endl;
        Abort("Dont know this plotfile_type");
    }

    // If fixed_mri_dt_ratio is set, it must be even
    if (fixed_mri_dt_ratio > 0 && (fixed_mri_dt_ratio%2 != 0) )
    {
        Abort("If you specify fixed_mri_dt_ratio, it must be even");
    }

    for (int lev = 0; lev <= max_level; lev++)
    {
        // We ignore fixed_fast_dt if not substepping
        if (solverChoice.substepping_type[lev] == SubsteppingType::None) {
            fixed_fast_dt[lev] = -1.0;
        }

        // If both fixed_dt and fast_dt are specified, their ratio must be an even integer
        if (fixed_dt[lev] > 0. && fixed_fast_dt[lev] > 0. && fixed_mri_dt_ratio <= 0)
        {
            Real eps = 1.e-12;
            int ratio = static_cast<int>( ( (1.0+eps) * fixed_dt[lev] ) / fixed_fast_dt[lev] );
            if (fixed_dt[lev] / fixed_fast_dt[lev] != ratio)
            {
                Abort("Ratio of fixed_dt to fixed_fast_dt must be an even integer");
            }
        }

        // If all three are specified, they must be consistent
        if (fixed_dt[lev] > 0. && fixed_fast_dt[lev] > 0. &&  fixed_mri_dt_ratio > 0)
        {
            if (fixed_dt[lev] / fixed_fast_dt[lev] != fixed_mri_dt_ratio)
            {
                Abort("Dt is over-specfied");
            }
        }
    } // lev

    if (solverChoice.coupling_type == CouplingType::TwoWay && cf_width > 0) {
        Abort("For two-way coupling you must set cf_width = 0");
    }
}

// Create horizontal average quantities for 5 variables:
//        density, temperature, pressure, qc, qv (if present)
void
ERF::MakeHorizontalAverages ()
{
    int lev = 0;

    // First, average down all levels (if doing two-way coupling)
    if (solverChoice.coupling_type == CouplingType::TwoWay) {
        AverageDown();
    }

    MultiFab mf(grids[lev], dmap[lev], 5, 0);

    int zdir = 2;
    auto domain = geom[0].Domain();

    bool use_moisture = (solverChoice.moisture_type != MoistureType::None);
    bool is_anelastic = (solverChoice.anelastic[lev] == 1);

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto  fab_arr = mf.array(mfi);
        auto const  hse_arr = base_state[lev].const_array(mfi);
        auto const cons_arr = vars_new[lev][Vars::cons].const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Real dens = cons_arr(i, j, k, Rho_comp);
            fab_arr(i, j, k, 0) = dens;
            fab_arr(i, j, k, 1) = cons_arr(i, j, k, RhoTheta_comp) / dens;
            if (!use_moisture) {
                if (is_anelastic) {
                    fab_arr(i,j,k,2) = hse_arr(i,j,k,1);
                } else {
                    fab_arr(i,j,k,2) = getPgivenRTh(cons_arr(i,j,k,RhoTheta_comp));
                }
            }
        });
    }

    if (use_moisture)
    {
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto  fab_arr = mf.array(mfi);
            auto const  hse_arr = base_state[lev].const_array(mfi);
            auto const cons_arr = vars_new[lev][Vars::cons].const_array(mfi);
            int ncomp = vars_new[lev][Vars::cons].nComp();

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real dens = cons_arr(i, j, k, Rho_comp);
                if (is_anelastic) {
                    fab_arr(i,j,k,2) = hse_arr(i,j,k,1);
                } else {
                    Real qv = cons_arr(i, j, k, RhoQ1_comp) / dens;
                    fab_arr(i, j, k, 2) = getPgivenRTh(cons_arr(i, j, k, RhoTheta_comp), qv);
                }
                fab_arr(i, j, k, 3) = (ncomp > RhoQ1_comp ? cons_arr(i, j, k, RhoQ1_comp) / dens : 0.0);
                fab_arr(i, j, k, 4) = (ncomp > RhoQ2_comp ? cons_arr(i, j, k, RhoQ2_comp) / dens : 0.0);
            });
        }

        Gpu::HostVector<Real> h_avg_qv          = sumToLine(mf,3,1,domain,zdir);
        Gpu::HostVector<Real> h_avg_qc          = sumToLine(mf,4,1,domain,zdir);
    }

    // Sum in the horizontal plane
    Gpu::HostVector<Real> h_avg_density     = sumToLine(mf,0,1,domain,zdir);
    Gpu::HostVector<Real> h_avg_temperature = sumToLine(mf,1,1,domain,zdir);
    Gpu::HostVector<Real> h_avg_pressure    = sumToLine(mf,2,1,domain,zdir);

    // Divide by the total number of cells we are averaging over
    int size_z = domain.length(zdir);
    Real area_z = static_cast<Real>(domain.length(0)*domain.length(1));
    int klen = static_cast<int>(h_avg_density.size());

    for (int k = 0; k < klen; ++k) {
        h_havg_density[k]     /= area_z;
        h_havg_temperature[k] /= area_z;
        h_havg_pressure[k]    /= area_z;
        if (solverChoice.moisture_type != MoistureType::None)
        {
            h_havg_qc[k]          /= area_z;
            h_havg_qv[k]          /= area_z;
        }
    } // k

    // resize device vectors
    d_havg_density.resize(size_z, 0.0_rt);
    d_havg_temperature.resize(size_z, 0.0_rt);
    d_havg_pressure.resize(size_z, 0.0_rt);

    // copy host vectors to device vectors
    Gpu::copy(Gpu::hostToDevice, h_havg_density.begin(), h_havg_density.end(), d_havg_density.begin());
    Gpu::copy(Gpu::hostToDevice, h_havg_temperature.begin(), h_havg_temperature.end(), d_havg_temperature.begin());
    Gpu::copy(Gpu::hostToDevice, h_havg_pressure.begin(), h_havg_pressure.end(), d_havg_pressure.begin());

    if (solverChoice.moisture_type != MoistureType::None)
    {
        d_havg_qv.resize(size_z, 0.0_rt);
        d_havg_qc.resize(size_z, 0.0_rt);
        Gpu::copy(Gpu::hostToDevice, h_havg_qv.begin(), h_havg_qv.end(), d_havg_qv.begin());
        Gpu::copy(Gpu::hostToDevice, h_havg_qc.begin(), h_havg_qc.end(), d_havg_qc.begin());
    }
}

// Create horizontal average quantities for the MultiFab passed in
// NOTE: this does not create device versions of the 1d arrays
// NOLINTNEXTLINE
void // NOLINTNEXTLINE
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
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const IntVect& se = box.smallEnd();
        const IntVect& be = box.bigEnd();

        auto      fab_arr = S[mfi].array();

        FArrayBox fab_reduce(box, 1, The_Async_Arena());
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
    AMREX_ALWAYS_ASSERT(solverChoice.coupling_type == CouplingType::TwoWay);
    int  src_comp = 0;
    int  num_comp = vars_new[0][Vars::cons].nComp();
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
        AverageDownTo(lev,src_comp,num_comp);
    }
}

// Set covered coarse cells to be the average of overlying fine cells at level crse_lev
void
ERF::AverageDownTo (int crse_lev, int scomp, int ncomp) // NOLINT
{
    AMREX_ALWAYS_ASSERT(solverChoice.coupling_type == CouplingType::TwoWay);

    // ******************************************************************************************
    // First do cell-centered quantities
    // The quantity that is conserved is not (rho S), but rather (rho S / m^2) where
    // m is the map scale factor at cell centers
    // Here we pre-divide (rho S) by m^2 before average down
    // ******************************************************************************************
    for (int lev = crse_lev; lev <= crse_lev+1; lev++) {
      for (MFIter mfi(vars_new[lev][Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const Array4<      Real>   cons_arr = vars_new[lev][Vars::cons].array(mfi);
        const Array4<const Real> mapfac_arr = mapfac_m[lev]->const_array(mfi);
        if (solverChoice.use_terrain) {
            const Array4<const Real>   detJ_arr = detJ_cc[lev]->const_array(mfi);
            ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                cons_arr(i,j,k,scomp+n) *= detJ_arr(i,j,k) / (mapfac_arr(i,j,0)*mapfac_arr(i,j,0));
            });
        } else {
            ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                cons_arr(i,j,k,scomp+n) /= (mapfac_arr(i,j,0)*mapfac_arr(i,j,0));
            });
        }
      } // mfi
    } // lev

    average_down(vars_new[crse_lev+1][Vars::cons],
                 vars_new[crse_lev  ][Vars::cons],
                 scomp, ncomp, refRatio(crse_lev));

    // Here we multiply (rho S) by m^2 after average down
    for (int lev = crse_lev; lev <= crse_lev+1; lev++) {
      for (MFIter mfi(vars_new[lev][Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const Array4<      Real>   cons_arr = vars_new[lev][Vars::cons].array(mfi);
        const Array4<const Real> mapfac_arr = mapfac_m[lev]->const_array(mfi);
        if (solverChoice.use_terrain) {
            const Array4<const Real>   detJ_arr = detJ_cc[lev]->const_array(mfi);
            ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                cons_arr(i,j,k,scomp+n) *= (mapfac_arr(i,j,0)*mapfac_arr(i,j,0)) / detJ_arr(i,j,k);
            });
        } else {
            ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                cons_arr(i,j,k,scomp+n) *= (mapfac_arr(i,j,0)*mapfac_arr(i,j,0));
            });
        }
      } // mfi
    } // lev

    // ******************************************************************************************
    // Now average down momenta.
    // Note that vars_new holds velocities not momenta, but we want to do conservative
    //    averaging so we first convert to momentum, then average down, then convert
    //    back to velocities -- only on the valid region
    // ******************************************************************************************
    for (int lev = crse_lev; lev <= crse_lev+1; lev++)
    {
        // FillBoundary for density so we can go back and forth between velocity and momentum
        vars_new[lev][Vars::cons].FillBoundary(geom[lev].periodicity());

        VelocityToMomentum(vars_new[lev][Vars::xvel], IntVect(0,0,0),
                           vars_new[lev][Vars::yvel], IntVect(0,0,0),
                           vars_new[lev][Vars::zvel], IntVect(0,0,0),
                           vars_new[lev][Vars::cons],
                             rU_new[lev],
                             rV_new[lev],
                             rW_new[lev],
                           Geom(lev).Domain(),
                           domain_bcs_type);
    }

    average_down_faces(rU_new[crse_lev+1], rU_new[crse_lev], refRatio(crse_lev), geom[crse_lev]);
    average_down_faces(rV_new[crse_lev+1], rV_new[crse_lev], refRatio(crse_lev), geom[crse_lev]);
    average_down_faces(rW_new[crse_lev+1], rW_new[crse_lev], refRatio(crse_lev), geom[crse_lev]);

    for (int lev = crse_lev; lev <= crse_lev+1; lev++) {
        MomentumToVelocity(vars_new[lev][Vars::xvel],
                           vars_new[lev][Vars::yvel],
                           vars_new[lev][Vars::zvel],
                           vars_new[lev][Vars::cons],
                             rU_new[lev],
                             rV_new[lev],
                             rW_new[lev],
                           Geom(lev).Domain(),
                           domain_bcs_type);
    }
}

void
ERF::Construct_ERFFillPatchers (int lev)
{
    auto& fine_new = vars_new[lev];
    auto& crse_new = vars_new[lev-1];
    auto& ba_fine  = fine_new[Vars::cons].boxArray();
    auto& ba_crse  = crse_new[Vars::cons].boxArray();
    auto& dm_fine  = fine_new[Vars::cons].DistributionMap();
    auto& dm_crse  = crse_new[Vars::cons].DistributionMap();

    int ncomp = vars_new[lev][Vars::cons].nComp();

    FPr_c.emplace_back(ba_fine, dm_fine, geom[lev]  ,
                       ba_crse, dm_crse, geom[lev-1],
                       -cf_width, -cf_set_width, ncomp, &cell_cons_interp);
    FPr_u.emplace_back(convert(ba_fine, IntVect(1,0,0)), dm_fine, geom[lev]  ,
                       convert(ba_crse, IntVect(1,0,0)), dm_crse, geom[lev-1],
                       -cf_width, -cf_set_width, 1, &face_cons_linear_interp);
    FPr_v.emplace_back(convert(ba_fine, IntVect(0,1,0)), dm_fine, geom[lev]  ,
                       convert(ba_crse, IntVect(0,1,0)), dm_crse, geom[lev-1],
                       -cf_width, -cf_set_width, 1, &face_cons_linear_interp);
    FPr_w.emplace_back(convert(ba_fine, IntVect(0,0,1)), dm_fine, geom[lev]  ,
                       convert(ba_crse, IntVect(0,0,1)), dm_crse, geom[lev-1],
                       -cf_width, -cf_set_width, 1, &face_cons_linear_interp);
}

void
ERF::Define_ERFFillPatchers (int lev)
{
    auto& fine_new = vars_new[lev];
    auto& crse_new = vars_new[lev-1];
    auto& ba_fine  = fine_new[Vars::cons].boxArray();
    auto& ba_crse  = crse_new[Vars::cons].boxArray();
    auto& dm_fine  = fine_new[Vars::cons].DistributionMap();
    auto& dm_crse  = crse_new[Vars::cons].DistributionMap();

    int ncomp = fine_new[Vars::cons].nComp();

    FPr_c[lev-1].Define(ba_fine, dm_fine, geom[lev]  ,
                        ba_crse, dm_crse, geom[lev-1],
                        -cf_width, -cf_set_width, ncomp, &cell_cons_interp);
    FPr_u[lev-1].Define(convert(ba_fine, IntVect(1,0,0)), dm_fine, geom[lev]  ,
                        convert(ba_crse, IntVect(1,0,0)), dm_crse, geom[lev-1],
                        -cf_width, -cf_set_width, 1, &face_cons_linear_interp);
    FPr_v[lev-1].Define(convert(ba_fine, IntVect(0,1,0)), dm_fine, geom[lev]  ,
                        convert(ba_crse, IntVect(0,1,0)), dm_crse, geom[lev-1],
                        -cf_width, -cf_set_width, 1, &face_cons_linear_interp);
    FPr_w[lev-1].Define(convert(ba_fine, IntVect(0,0,1)), dm_fine, geom[lev]  ,
                        convert(ba_crse, IntVect(0,0,1)), dm_crse, geom[lev-1],
                        -cf_width, -cf_set_width, 1, &face_cons_linear_interp);
}

#ifdef ERF_USE_MULTIBLOCK
// constructor used when ERF is created by a multiblock driver
ERF::ERF (const RealBox& rb, int max_level_in,
          const Vector<int>& n_cell_in, int coord,
          const Vector<IntVect>& ref_ratio,
          const Array<int,AMREX_SPACEDIM>& is_per,
          std::string prefix)
    : AmrCore(rb, max_level_in, n_cell_in, coord, ref_ratio, is_per)
{
    SetParmParsePrefix(prefix);

    // Multiblock: public domain sizes (need to know which vars are nodal)
    Box nbx;
    domain_p.push_back(geom[0].Domain());
    nbx = convert(domain_p[0],IntVect(1,0,0));
    domain_p.push_back(nbx);
    nbx = convert(domain_p[0],IntVect(0,1,0));
    domain_p.push_back(nbx);
    nbx = convert(domain_p[0],IntVect(0,0,1));
    domain_p.push_back(nbx);

    ERF_shared();
}
#endif

bool
ERF::writeNow(const Real cur_time, const Real dt_lev, const int nstep, const int plot_int, const Real plot_per)
{
    bool write_now = false;

    if ( plot_int > 0 && (nstep % plot_int == 0) ) {
        write_now = true;

    } else if (plot_per > 0.0) {

        // Check to see if we've crossed a plot_per interval by comparing
        // the number of intervals that have elapsed for both the current
        // time and the time at the beginning of this timestep.

        const Real eps = std::numeric_limits<Real>::epsilon() * Real(10.0) * std::abs(cur_time);

        int num_per_old = static_cast<int>(std::floor((cur_time-eps-dt_lev) / plot_per));
        int num_per_new = static_cast<int>(std::floor((cur_time-eps       ) / plot_per));

        // Before using these, however, we must test for the case where we're
        // within machine epsilon of the next interval. In that case, increment
        // the counter, because we have indeed reached the next plot_per interval
        // at this point.

        const Real next_plot_time = (num_per_old + 1) * plot_per;

        if ((num_per_new == num_per_old) && std::abs(cur_time - next_plot_time) <= eps)
        {
            num_per_new += 1;
        }

        // Similarly, we have to account for the case where the old time is within
        // machine epsilon of the beginning of this interval, so that we don't double
        // count that time threshold -- we already plotted at that time on the last timestep.

        if ((num_per_new != num_per_old) && std::abs((cur_time - dt_lev) - next_plot_time) <= eps)
            num_per_old += 1;

        if (num_per_old != num_per_new)
            write_now = true;
    }
    return write_now;
}
