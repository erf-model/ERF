/**
 * \file ERF.cpp
 */

/**
 * Main class in ERF code, instantiated from main.cpp
*/

#include <EOS.H>
#include <ERF.H>

#include <AMReX_buildInfo.H>

#include <Utils.H>
#include <TerrainMetrics.H>
#include <memory>

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

// Dictate verbosity in screen output
int         ERF::verbose       = 0;

// Frequency of diagnostic output
int         ERF::sum_interval  = -1;
amrex::Real ERF::sum_per       = -1.0;

// Native AMReX vs NetCDF
std::string ERF::plotfile_type    = "amrex";

// init_type:  "uniform", "ideal", "real", "input_sounding", "metgrid" or ""
std::string ERF::init_type;

// use_real_bcs: only true if 1) ( (init_type == real) or (init_type == metgrid) )
//                        AND 2) we want to use the bc's from the WRF bdy file
bool ERF::use_real_bcs;

// NetCDF wrfinput (initialization) file(s)
amrex::Vector<amrex::Vector<std::string>> ERF::nc_init_file = {{""}}; // Must provide via input

// NetCDF wrfbdy (lateral boundary) file
std::string ERF::nc_bdy_file; // Must provide via input

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

    int nlevs_max = max_level + 1;

    // NOTE: size micro before readparams (chooses the model at all levels)
    micro.ReSize(nlevs_max);
    qmoist.resize(nlevs_max);

#ifdef ERF_USE_WINDFARM
    if(solverChoice.windfarm_type == WindFarmType::Fitch){
        Nturb.resize(nlevs_max);
        vars_fitch.resize(nlevs_max);
    }
#endif

#if defined(ERF_USE_RRTMGP)
    qheating_rates.resize(nlevs_max);
#endif

    // NOTE: size lsm before readparams (chooses the model at all levels)
    lsm.ReSize(nlevs_max);
    lsm_data.resize(nlevs_max);
    lsm_flux.resize(nlevs_max);

    ReadParameters();
    const std::string& pv1 = "plot_vars_1"; setPlotVariables(pv1,plot_var_names_1);
    const std::string& pv2 = "plot_vars_2"; setPlotVariables(pv2,plot_var_names_2);

    // Initialize staggered vertical levels for grid stretching or terrain.

    if (solverChoice.use_terrain) {
        init_zlevels(zlevels_stag,
                     geom[0],
                     solverChoice.grid_stretching_ratio,
                     solverChoice.zsurf,
                     solverChoice.dz0);

        int nz = geom[0].Domain().length(2) + 1; // staggered
        if (zlevels_stag[nz-1] != geom[0].ProbHi(2)) {
            amrex::Print() << "WARNING: prob_hi[2]=" << geom[0].ProbHi(2)
                << " does not match highest requested z level " << zlevels_stag[nz-1]
                << std::endl;
        }
        if (zlevels_stag[0] != geom[0].ProbLo(2)) {
            amrex::Print() << "WARNING: prob_lo[2]=" << geom[0].ProbLo(2)
                << " does not match lowest requested level " << zlevels_stag[0]
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

    mri_integrator_mem.resize(nlevs_max);
    physbcs.resize(nlevs_max);

    advflux_reg.resize(nlevs_max);

    // Stresses
    Tau11_lev.resize(nlevs_max); Tau22_lev.resize(nlevs_max); Tau33_lev.resize(nlevs_max);
    Tau12_lev.resize(nlevs_max); Tau21_lev.resize(nlevs_max);
    Tau13_lev.resize(nlevs_max); Tau31_lev.resize(nlevs_max);
    Tau23_lev.resize(nlevs_max); Tau32_lev.resize(nlevs_max);
    SFS_hfx1_lev.resize(nlevs_max); SFS_hfx2_lev.resize(nlevs_max); SFS_hfx3_lev.resize(nlevs_max);
    SFS_diss_lev.resize(nlevs_max);
    eddyDiffs_lev.resize(nlevs_max);
    SmnSmn_lev.resize(nlevs_max);

    // Sea surface temps
    sst_lev.resize(nlevs_max);
    lmask_lev.resize(nlevs_max);

    // Metric terms
    z_phys_nd.resize(nlevs_max);
    z_phys_cc.resize(nlevs_max);
    detJ_cc.resize(nlevs_max);
    z_phys_nd_new.resize(nlevs_max);
    detJ_cc_new.resize(nlevs_max);
    z_phys_nd_src.resize(nlevs_max);
    detJ_cc_src.resize(nlevs_max);
    z_t_rk.resize(nlevs_max);

    // Mapfactors
    mapfac_m.resize(nlevs_max);
    mapfac_u.resize(nlevs_max);
    mapfac_v.resize(nlevs_max);

    // Base state
    base_state.resize(nlevs_max);
    base_state_new.resize(nlevs_max);

    // Theta prim for MOST
    Theta_prim.resize(nlevs_max);

    // Initialize tagging criteria for mesh refinement
    refinement_criteria_setup();

    // We have already read in the ref_Ratio (via amr.ref_ratio =) but we need to enforce
    //     that there is no refinement in the vertical so we test on that here.
    for (int lev = 0; lev < max_level; ++lev)
    {
       amrex::Print() << "Refinement ratio at level " << lev+1 << " set to be " <<
          ref_ratio[lev][0]  << " " << ref_ratio[lev][1]  <<  " " << ref_ratio[lev][2] << std::endl;

       if (ref_ratio[lev][2] != 1)
       {
           amrex::Error("We don't allow refinement in the vertical -- make sure to set ref_ratio = 1 in z");
       }
    }
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
        int  src_comp_reflux = 0;
        int  num_comp_reflux = vars_new[0][Vars::cons].nComp();
        bool use_terrain = solverChoice.use_terrain;
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
                    ParallelFor(bx, num_comp_reflux, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        cons_arr(i,j,k,src_comp_reflux+n) *= detJ_arr(i,j,k) / (mapfac_arr(i,j,0)*mapfac_arr(i,j,0));
                    });
                } else {
                    ParallelFor(bx, num_comp_reflux, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        cons_arr(i,j,k,src_comp_reflux+n) /= (mapfac_arr(i,j,0)*mapfac_arr(i,j,0));
                    });
                }
            } // mfi

            // This call refluxes from the lev/lev+1 interface onto lev
            getAdvFluxReg(lev+1)->Reflux(vars_new[lev][Vars::cons],
                                         src_comp_reflux, src_comp_reflux, num_comp_reflux);

            // Here we multiply (rho S) by m^2 after refluxing
            for (MFIter mfi(vars_new[lev][Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                const Array4<      Real>   cons_arr = vars_new[lev][Vars::cons].array(mfi);
                const Array4<const Real> mapfac_arr = mapfac_m[lev]->const_array(mfi);
                if (use_terrain) {
                    const Array4<const Real>   detJ_arr = detJ_cc[lev]->const_array(mfi);
                    ParallelFor(bx, num_comp_reflux, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        cons_arr(i,j,k,src_comp_reflux+n) *= (mapfac_arr(i,j,0)*mapfac_arr(i,j,0)) / detJ_arr(i,j,k);
                    });
                } else {
                    ParallelFor(bx, num_comp_reflux, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        cons_arr(i,j,k,src_comp_reflux+n) *= (mapfac_arr(i,j,0)*mapfac_arr(i,j,0));
                    });
                }
            } // mfi

            // We need to do this before anything else because refluxing changes the
            // values of coarse cells underneath fine grids with the assumption they'll
            // be over-written by averaging down
            AverageDownTo(lev, src_comp_reflux, num_comp_reflux);
        }
    }

    if (is_it_time_for_action(nstep, time, dt_lev0, sum_interval, sum_per)) {
        sum_integrated_quantities(time);
    }

    if (profile_int > 0 && (nstep+1) % profile_int == 0) {
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

    // Initialize the start time for our CPU-time tracker
    startCPUTime = amrex::ParallelDescriptor::second();

    // Map the words in the inputs file to BC types, then translate
    //     those types into what they mean for each variable
    init_bcs();

    // Verify BCs are compatible with solver choice
    for (int lev(0); lev <= max_level; ++lev) {
        if ( ( (solverChoice.turbChoice[lev].pbl_type == PBLType::MYNN25) ||
               (solverChoice.turbChoice[lev].pbl_type == PBLType::YSU)       ) &&
            phys_bc_type[Orientation(Direction::z,Orientation::low)] != ERF_BC::MOST ) {
            amrex::Abort("MYNN2.5/YSU PBL Model requires MOST at lower boundary");
        }
    }

    if (!solverChoice.use_terrain && solverChoice.terrain_type != TerrainType::Static) {
        amrex::Abort("We do not allow non-static terrain_type with use_terrain = false");
    }

    last_plot_file_step_1 = -1;
    last_plot_file_step_2 = -1;
    last_check_file_step  = -1;

    if (restart_chkfile.empty()) {
        // start simulation from the beginning

        const Real time = start_time;
        InitFromScratch(time);

#ifdef ERF_USE_MULTIBLOCK
        // Multiblock: hook to set BL & comms once ba/dm are known
        if(domain_p[0].bigEnd(0) < 500 ) {
            m_mbc->SetBoxLists();
            m_mbc->SetBlockCommMetaData();
        }
#endif

        if (solverChoice.use_terrain) {
            if (init_type == "ideal") {
                amrex::Abort("We do not currently support init_type = ideal with terrain");
            }
        }

        // If using the Deardoff LES model,
        // we initialize rho_KE to be nonzero (and positive) so that we end up
        // with reasonable values for the eddy diffusivity and the MOST fluxes
        // (~ 1/diffusivity) do not blow up
        Real RhoKE_0 = 0.1;
        ParmParse pp(pp_prefix);
        if (pp.query("RhoKE_0", RhoKE_0)) {
            // Uniform initial rho*e field
            for (int lev = 0; lev <= finest_level; lev++) {
                if (solverChoice.turbChoice[lev].les_type == LESType::Deardorff) {
                    vars_new[lev][Vars::cons].setVal(RhoKE_0,RhoKE_comp,1,0);
                } else {
                    vars_new[lev][Vars::cons].setVal(0.0,RhoKE_comp,1,0);
                }
            }
        } else {
            // default: uniform initial e field
            Real KE_0 = 0.1;
            pp.query("KE_0", KE_0);
            for (int lev = 0; lev <= finest_level; lev++) {
                auto& lev_new = vars_new[lev];
                if (solverChoice.turbChoice[lev].les_type == LESType::Deardorff) {
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
            int src_comp_reflux = 0;
            int num_comp_reflux = vars_new[0][Vars::cons].nComp();
            AverageDown(src_comp_reflux, num_comp_reflux);
        }

#ifdef ERF_USE_PARTICLES
        initializeTracers((amrex::ParGDBBase*)GetParGDB(),z_phys_nd);
#endif

    } else { // Restart from a checkpoint

        restart();


        // TODO: Check if this is needed. I don't think it is since we now
        //       advect all the scalars...

        // Need to fill ghost cells here since we will use this qmoist in advance
        if (solverChoice.moisture_type != MoistureType::None)
        {
            for (int lev = 0; lev <= finest_level; lev++) {
                FillPatchMoistVars(lev, *(qmoist[lev][0])); // qv component
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
    if (solverChoice.coupling_type == CouplingType::TwoWay) {
        advflux_reg[0] = nullptr;
        int ncomp_reflux = vars_new[0][Vars::cons].nComp();
        for (int lev = 1; lev <= finest_level; lev++)
        {
            advflux_reg[lev] = new YAFluxRegister(grids[lev], grids[lev-1],
                                                   dmap[lev],  dmap[lev-1],
                                                   geom[lev],  geom[lev-1],
                                              ref_ratio[lev-1], lev, ncomp_reflux);
        }
    }

    // Configure ABLMost params if used MostWall boundary condition
    // NOTE: we must set up the MOST routine before calling WritePlotFile because
    //       WritePlotFile calls FillPatch in order to compute gradients
    if (phys_bc_type[Orientation(Direction::z,Orientation::low)] == ERF_BC::MOST)
    {
        m_most = std::make_unique<ABLMost>(geom, vars_old, Theta_prim, z_phys_nd,
                                           sst_lev, lmask_lev, lsm_data, lsm_flux
#ifdef ERF_USE_NETCDF
                                           ,start_bdy_time, bdy_time_interval
#endif
                                           );

        // We now configure ABLMost params here so that we can print the averages at t=0
        // Note we don't fill ghost cells here because this is just for diagnostics
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            amrex::IntVect ng = IntVect(0,0,0);
            MultiFab S(vars_new[lev][Vars::cons],make_alias,0,2);
            MultiFab::Copy(  *Theta_prim[lev], S, RhoTheta_comp, 0, 1, ng);
            MultiFab::Divide(*Theta_prim[lev], S, Rho_comp     , 0, 1, ng);
            m_most->update_mac_ptrs(lev, vars_new, Theta_prim);
            m_most->update_fluxes(lev, t_new[lev]);
        }
    }

    if (solverChoice.custom_rhotheta_forcing)
    {
        h_rhotheta_src.resize(max_level+1, amrex::Vector<Real>(0));
        d_rhotheta_src.resize(max_level+1, amrex::Gpu::DeviceVector<Real>(0));
        for (int lev = 0; lev <= finest_level; lev++)
        {
            const int domlen = geom[lev].Domain().length(2);
            h_rhotheta_src[lev].resize(domlen, 0.0_rt);
            d_rhotheta_src[lev].resize(domlen, 0.0_rt);
            prob->update_rhotheta_sources(t_new[0],
                                          h_rhotheta_src[lev], d_rhotheta_src[lev],
                                          geom[lev], z_phys_cc[lev]);
        }
    }

    if (solverChoice.use_rayleigh_damping)
    {
        initRayleigh();
        if (init_type == "input_sounding")
        {
            // Overwrite ubar, vbar, and thetabar with input profiles; wbar is
            // assumed to be 0. Note: the tau coefficient set by
            // prob->erf_init_rayleigh() is still used
            bool restarting = (!restart_chkfile.empty());
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

#ifdef ERF_USE_POISSON_SOLVE
    if (restart_chkfile == "")
    {
        // Note -- this projection is only defined for no terrain
        if (solverChoice.project_initial_velocity) {
            AMREX_ALWAYS_ASSERT(solverChoice.use_terrain == 0);
            project_velocities(vars_new);
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

    // Compute the minimum dz in the domain (to be used for setting the timestep)
    dz_min = geom[0].CellSize(2);
    if ( solverChoice.use_terrain ) {
        dz_min *= (*detJ_cc[0]).min(0);
    }

    ComputeDt();

    // Fill ghost cells/faces
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        if (lev > 0 && cf_width > 0) {
            Construct_ERFFillPatchers(lev);
        }


        //
        // We don't use the FillPatcher in this call because
        //    we don't need to fill the interior data at this point.
        //
        bool fillset = false;
        auto& lev_new = vars_new[lev];
        FillPatch(lev, t_new[lev],
                  {&lev_new[Vars::cons],&lev_new[Vars::xvel],&lev_new[Vars::yvel],&lev_new[Vars::zvel]},
                  fillset);

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

    // Update micro vars before first plot file
    if (solverChoice.moisture_type != MoistureType::None) {
        for (int lev = 0; lev <= finest_level; ++lev) micro.Update_Micro_Vars_Lev(lev, vars_new[lev][Vars::cons]);
    }


    if (restart_chkfile.empty() && check_int > 0)
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

    BL_PROFILE_VAR_STOP(InitData);
}

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
    if (init_type == "input_sounding") {
        // The base state is initialized by integrating vertically through the
        // input sounding, if the init_sounding_ideal flag is set; otherwise
        // it is set by initHSE()
        init_from_input_sounding(lev);
        if (!init_sounding_ideal) initHSE();

#ifdef ERF_USE_NETCDF
    } else if (init_type == "ideal" || init_type == "real") {
        // The base state is initialized from WRF wrfinput data, output by
        // ideal.exe or real.exe
        init_from_wrfinput(lev);
        if (init_type == "ideal") initHSE();

    } else if (init_type == "metgrid") {
        // The base state is initialized from data output by WPS metgrid;
        // we will rebalance after interpolation
        init_from_metgrid(lev);
#endif
    } else if (init_type == "uniform") {
        // Initialize a uniform background field and base state based on the
        // problem-specified reference density and temperature
        init_uniform(lev);
        initHSE(lev);
    } else {
        // No background flow initialization specified, initialize the
        // background field to be equal to the base state, calculated from the
        // problem-specific erf_init_dens_hse
        initHSE(lev); // need to call this first
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
}

// read in some parameters from inputs file
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

        // The regression tests use "amr.restart" and "amr.check_int" so we allow
        //    for those or "erf.restart" / "erf.check_int" with the former taking
        //    precedenceif both are specified
        pp.query("check_int", check_int);
        pp_amr.query("check_int", check_int);

        pp.query("restart", restart_chkfile);
        pp_amr.query("restart", restart_chkfile);

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

        AMREX_ALWAYS_ASSERT(cfl > 0. || fixed_dt > 0.);

        // How to initialize
        pp.query("init_type",init_type);
        if (!init_type.empty() &&
            init_type != "uniform" &&
            init_type != "ideal" &&
            init_type != "real" &&
            init_type != "metgrid" &&
            init_type != "input_sounding")
        {
            amrex::Error("if specified, init_type must be uniform, ideal, real, metgrid or input_sounding");
        }

        // Should we use the bcs we've read in from wrfbdy or metgrid files?
        // We default to yes if we have them, but the user can override that option
        use_real_bcs = ( (init_type == "real") || (init_type == "metgrid") );
        pp.query("use_real_bcs",use_real_bcs);

        // We don't allow use_real_bcs to be true if init_type is not either real or metgrid
        AMREX_ALWAYS_ASSERT(!use_real_bcs || ((init_type == "real") || (init_type == "metgrid")) );

        // No moving terrain with init real
        if (init_type == "real" && solverChoice.terrain_type != TerrainType::Static) {
            amrex::Abort("Moving terrain is not supported with init real");
        }

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

        pp.query("profile_int", profile_int);

        pp.query("plot_lsm", plot_lsm);

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

        // Query the set and total widths for wrfbdy interior ghost cells
        pp.query("real_width", real_width);
        pp.query("real_set_width", real_set_width);
        AMREX_ALWAYS_ASSERT(real_width >= 0);
        AMREX_ALWAYS_ASSERT(real_set_width >= 0);
        AMREX_ALWAYS_ASSERT(real_width >= real_set_width);

        // Query the set and total widths for crse-fine interior ghost cells
        pp.query("cf_width", cf_width);
        pp.query("cf_set_width", cf_set_width);
        if (cf_width < 0 || cf_set_width < 0 || cf_width < cf_set_width) {
            amrex::Abort("You must set cf_width >= cf_set_width >= 0");
        }

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

    // What type of moisture model to use
    // NOTE: Must be checked after init_params
    if (solverChoice.moisture_type == MoistureType::SAM) {
        micro.SetModel<SAM>();
        amrex::Print() << "SAM moisture model!\n";
    } else if (solverChoice.moisture_type == MoistureType::Kessler) {
        micro.SetModel<Kessler>();
        amrex::Print() << "Kessler moisture model!\n";
    } else if (solverChoice.moisture_type == MoistureType::FastEddy) {
        micro.SetModel<FastEddy>();
        amrex::Print() << "FastEddy moisture model!\n";
    } else if (solverChoice.moisture_type == MoistureType::None) {
        micro.SetModel<NullMoist>();
        amrex::Print() << "WARNING: Compiled with moisture but using NullMoist model!\n";
    } else {
        amrex::Abort("Dont know this moisture_type!") ;
    }

    // What type of land surface model to use
    // NOTE: Must be checked after init_params
    if (solverChoice.lsm_type == LandSurfaceType::SLM) {
        lsm.SetModel<SLM>();
        amrex::Print() << "SLM land surface model!\n";
    } else if (solverChoice.lsm_type == LandSurfaceType::MM5) {
        lsm.SetModel<MM5>();
        amrex::Print() << "MM5 land surface model!\n";
    } else if (solverChoice.lsm_type == LandSurfaceType::None) {
        lsm.SetModel<NullSurf>();
        amrex::Print() << "Null land surface model!\n";
    } else {
        amrex::Abort("Dont know this moisture_type!") ;
    }

    if (verbose > 0) {
        solverChoice.display();
    }

    if (solverChoice.coupling_type == CouplingType::TwoWay && cf_width > 0) {
        amrex::Abort("For two-way coupling you must set cf_width = 0");
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
        int  src_comp_reflux = 0;
        int  num_comp_reflux = vars_new[lev][Vars::cons].nComp();
        AverageDown(src_comp_reflux, num_comp_reflux);
    }

    MultiFab mf(grids[lev], dmap[lev], 5, 0);

    int zdir = 2;
    auto domain = geom[0].Domain();

    bool use_moisture = (solverChoice.moisture_type != MoistureType::None);

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto  fab_arr = mf.array(mfi);
        auto const cons_arr = vars_new[lev][Vars::cons].const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Real dens = cons_arr(i, j, k, Rho_comp);
            fab_arr(i, j, k, 0) = dens;
            fab_arr(i, j, k, 1) = cons_arr(i, j, k, RhoTheta_comp) / dens;
            if (!use_moisture) {
                fab_arr(i, j, k, 2) = getPgivenRTh(cons_arr(i, j, k, RhoTheta_comp));
            }
        });
    }

    if (use_moisture)
    {
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto  fab_arr = mf.array(mfi);
            auto const cons_arr = vars_new[lev][Vars::cons].const_array(mfi);
            auto const   qv_arr = qmoist[lev][0]->const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real dens = cons_arr(i, j, k, Rho_comp);
                fab_arr(i, j, k, 2) = getPgivenRTh(cons_arr(i, j, k, RhoTheta_comp), qv_arr(i,j,k));
                fab_arr(i, j, k, 3) = cons_arr(i, j, k, RhoQ1_comp) / dens;
                fab_arr(i, j, k, 4) = cons_arr(i, j, k, RhoQ2_comp) / dens;
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
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_density.begin(), h_havg_density.end(), d_havg_density.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_temperature.begin(), h_havg_temperature.end(), d_havg_temperature.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_pressure.begin(), h_havg_pressure.end(), d_havg_pressure.begin());

    if (solverChoice.moisture_type != MoistureType::None)
    {
        d_havg_qv.resize(size_z, 0.0_rt);
        d_havg_qc.resize(size_z, 0.0_rt);
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_qv.begin(), h_havg_qv.end(), d_havg_qv.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_qc.begin(), h_havg_qc.end(), d_havg_qc.begin());
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
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
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
ERF::AverageDown (int scomp, int ncomp)
{
    AMREX_ALWAYS_ASSERT(solverChoice.coupling_type == CouplingType::TwoWay);
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
        AverageDownTo(lev, scomp, ncomp);
    }
}

// Set covered coarse cells to be the average of overlying fine cells at level crse_lev
void
ERF::AverageDownTo (int crse_lev, int scomp, int ncomp) // NOLINT
{
    AMREX_ALWAYS_ASSERT(solverChoice.coupling_type == CouplingType::TwoWay);

    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx) {
        const BoxArray& ba(vars_new[crse_lev][var_idx].boxArray());
        if (ba[0].type() == IntVect::TheZeroVector())
        {
            // The quantity that is conserved is not (rho S), but rather (rho S / m^2) where
            // m is the map scale factor at cell centers
            // Here we pre-divide (rho S) by m^2 before average down
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

            amrex::average_down(vars_new[crse_lev+1][var_idx],
                                vars_new[crse_lev  ][var_idx],
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

        } else { // We assume the arrays are face-centered if not cell-centered, and
                 //    only average down the momenta if using two-way coupling
            if (solverChoice.coupling_type == CouplingType::TwoWay) {
                amrex::average_down_faces(vars_new[crse_lev+1][var_idx], vars_new[crse_lev][var_idx],
                                          refRatio(crse_lev),geom[crse_lev]);
            }
        }
    }
}

void
ERF::Construct_ERFFillPatchers (int lev)
{
    AMREX_ALWAYS_ASSERT(solverChoice.coupling_type == CouplingType::OneWay);

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
                       -cf_width, -cf_set_width, 1, &face_linear_interp);
    FPr_v.emplace_back(convert(ba_fine, IntVect(0,1,0)), dm_fine, geom[lev]  ,
                       convert(ba_crse, IntVect(0,1,0)), dm_crse, geom[lev-1],
                       -cf_width, -cf_set_width, 1, &face_linear_interp);
    FPr_w.emplace_back(convert(ba_fine, IntVect(0,0,1)), dm_fine, geom[lev]  ,
                       convert(ba_crse, IntVect(0,0,1)), dm_crse, geom[lev-1],
                       -cf_width, -cf_set_width, 1, &face_linear_interp);
}

void
ERF::Define_ERFFillPatchers (int lev)
{
    AMREX_ALWAYS_ASSERT(solverChoice.coupling_type == CouplingType::OneWay);

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
                        -cf_width, -cf_set_width, 1, &face_linear_interp);
    FPr_v[lev-1].Define(convert(ba_fine, IntVect(0,1,0)), dm_fine, geom[lev]  ,
                        convert(ba_crse, IntVect(0,1,0)), dm_crse, geom[lev-1],
                        -cf_width, -cf_set_width, 1, &face_linear_interp);
    FPr_w[lev-1].Define(convert(ba_fine, IntVect(0,0,1)), dm_fine, geom[lev]  ,
                        convert(ba_crse, IntVect(0,0,1)), dm_crse, geom[lev-1],
                        -cf_width, -cf_set_width, 1, &face_linear_interp);
}

void
ERF::Register_ERFFillPatchers (int lev)
{
    auto& lev_new = vars_new[lev];
    auto& lev_old = vars_old[lev];
    FPr_c[lev].RegisterCoarseData({&lev_old[Vars::cons], &lev_new[Vars::cons]}, {t_old[lev], t_new[lev]});
    FPr_u[lev].RegisterCoarseData({&lev_old[Vars::xvel], &lev_new[Vars::xvel]}, {t_old[lev], t_new[lev]});
    FPr_v[lev].RegisterCoarseData({&lev_old[Vars::yvel], &lev_new[Vars::yvel]}, {t_old[lev], t_new[lev]});
    FPr_w[lev].RegisterCoarseData({&lev_old[Vars::zvel], &lev_new[Vars::zvel]}, {t_old[lev], t_new[lev]});
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

    int nlevs_max = max_level + 1;

    // NOTE: size micro before readparams (chooses the model at all levels)
    micro.ReSize(nlevs_max);
    qmoist.resize(nlevs_max);

#ifdef ERF_USE_WINDFARM
    if(solverChoice.windfarm_type == WindFarmType::Fitch){
        Nturb.resize(nlevs_max);
        vars_fitch.resize(nlevs_max);
    }
#endif

#if defined(ERF_USE_RRTMGP)
    qheating_rates.resize(nlevs_max);
#endif

    // NOTE: size micro before readparams (chooses the model at all levels)
    lsm.ReSize(nlevs_max);
    lsm_data.resize(nlevs_max);
    lsm_flux.resize(nlevs_max);

    ReadParameters();
    const std::string& pv1 = "plot_vars_1"; setPlotVariables(pv1,plot_var_names_1);
    const std::string& pv2 = "plot_vars_2"; setPlotVariables(pv2,plot_var_names_2);

    prob = amrex_probinit(geom[0].ProbLo(), geom[0].ProbHi());

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

    for (int lev = 0; lev < nlevs_max; ++lev) {
        vars_new[lev].resize(Vars::NumTypes);
        vars_old[lev].resize(Vars::NumTypes);
    }

    rU_new.resize(nlevs_max);
    rV_new.resize(nlevs_max);
    rW_new.resize(nlevs_max);

    rU_old.resize(nlevs_max);
    rV_old.resize(nlevs_max);
    rW_old.resize(nlevs_max);

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

    advflux_reg.resize(nlevs_max);

    // Stresses
    Tau11_lev.resize(nlevs_max); Tau22_lev.resize(nlevs_max); Tau33_lev.resize(nlevs_max);
    Tau12_lev.resize(nlevs_max); Tau21_lev.resize(nlevs_max);
    Tau13_lev.resize(nlevs_max); Tau31_lev.resize(nlevs_max);
    Tau23_lev.resize(nlevs_max); Tau32_lev.resize(nlevs_max);
    SFS_hfx1_lev.resize(nlevs_max); SFS_hfx2_lev.resize(nlevs_max); SFS_hfx3_lev.resize(nlevs_max);
    SFS_diss_lev.resize(nlevs_max);
    eddyDiffs_lev.resize(nlevs_max);
    SmnSmn_lev.resize(nlevs_max);

    // Sea surface temps
    sst_lev.resize(nlevs_max);
    lmask_lev.resize(nlevs_max);

    // Metric terms
    z_phys_nd.resize(nlevs_max);
    z_phys_cc.resize(nlevs_max);
    detJ_cc.resize(nlevs_max);
    z_phys_nd_new.resize(nlevs_max);
    detJ_cc_new.resize(nlevs_max);
    z_phys_nd_src.resize(nlevs_max);
    detJ_cc_src.resize(nlevs_max);
    z_t_rk.resize(nlevs_max);

    // Mapfactors
    mapfac_m.resize(nlevs_max);
    mapfac_u.resize(nlevs_max);
    mapfac_v.resize(nlevs_max);

    // Base state
    base_state.resize(nlevs_max);
    base_state_new.resize(nlevs_max);

    // Theta prim for MOST
    Theta_prim.resize(nlevs_max);

    // Initialize tagging criteria for mesh refinement
    refinement_criteria_setup();

    // We have already read in the ref_Ratio (via amr.ref_ratio =) but we need to enforce
    //     that there is no refinement in the vertical so we test on that here.
    for (int lev = 0; lev < max_level; ++lev)
    {
       amrex::Print() << "Refinement ratio at level " << lev+1 << " set to be " <<
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
