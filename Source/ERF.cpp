/**
 * \file ERF.cpp
 */

/**
 * Main class in ERF code, instantiated from main.cpp
*/

#include <memory>

#include "ERF_EOS.H"
#include "ERF.H"
#include "AMReX_buildInfo.H"
#include "AMReX_Random.H"
#include "AMReX_WriteEBSurface.H"
#include "AMReX_EB2_IF_Box.H"
#include "AMReX_EB2_IF_Sphere.H"

#include "ERF_EpochTime.H"
#include "ERF_Utils.H"
#include "ERF_TerrainMetrics.H"
#include "ERF_EBIFTerrain.H"
#include "ERF_HurricaneDiagnostics.H"

#ifdef ERF_USE_NETCDF
#include "ERF_ReadFromWRFInput.H"
#include "ERF_ReadFromWRFBdy.H"
#endif

using namespace amrex;

Real ERF::startCPUTime        = 0.0;
Real ERF::previousCPUTimeUsed = 0.0;

Vector<AMRErrorTag> ERF::ref_tags;

SolverChoice ERF::solverChoice;

Real ERF::start_time    = 0.0;
Real ERF::stop_time     = std::numeric_limits<amrex::Real>::max();

#ifdef ERF_USE_NETCDF
Real ERF::start_bdy_time     = 0.0;
Real ERF::start_low_time     = 0.0;

Real ERF::bdy_time_interval  = std::numeric_limits<amrex::Real>::max();
Real ERF::low_time_interval  = std::numeric_limits<amrex::Real>::max();
#endif

// Time step control
Real ERF::cfl            = 0.8;
Real ERF::sub_cfl        = 1.0;
Real ERF::init_shrink    = 1.0;
Real ERF::change_max     = 1.1;
Real ERF::dt_max_initial = 2.0e100;
Real ERF:: dt_max        = 1.0e9;

int  ERF::fixed_mri_dt_ratio = 0;

// Dictate verbosity in screen output
int  ERF::verbose       = 0;
int  ERF::mg_verbose    = 0;
bool ERF::use_fft       = false;

// Should we check the solution for NaNs every time step?
// 1: check state/vels after dycore, state after microphysics, and state/vels at end of full time step
// 2: add checks of state before dycore and of slow rhs
int ERF::check_for_nans = 0;

// Frequency of diagnostic output
int  ERF::sum_interval  = -1;
Real ERF::sum_per       = -1.0;

int  ERF::pert_interval = -1;

int ERF::last_plot3d_file_step_1 = -1;
int ERF::last_plot3d_file_step_2 = -1;
int ERF::last_plot2d_file_step_1 = -1;
int ERF::last_plot2d_file_step_2 = -1;
int ERF::last_check_file_step    = -1;

Real ERF::last_plot3d_file_time_1 = 0.0;
Real ERF::last_plot3d_file_time_2 = 0.0;
Real ERF::last_plot2d_file_time_1 = 0.0;
Real ERF::last_plot2d_file_time_2 = 0.0;
Real ERF::last_check_file_time    = 0.0;

bool ERF::plot_file_on_restart = true;

// Native AMReX vs NetCDF
PlotFileType ERF::plotfile3d_type_1  = PlotFileType::None;
PlotFileType ERF::plotfile3d_type_2  = PlotFileType::None;
PlotFileType ERF::plotfile2d_type_1  = PlotFileType::None;
PlotFileType ERF::plotfile2d_type_2  = PlotFileType::None;

StateInterpType ERF::interpolation_type;

// NetCDF wrfinput (initialization) file(s)
Vector<Vector<std::string>> ERF::nc_init_file = {{""}}; // Must provide via input
Vector<Vector<int>>         ERF::have_read_nc_init_file = {{0}};

// NetCDF wrfbdy (lateral boundary) file
std::string ERF::nc_bdy_file; // Must provide via input

// NetCDF wrflow (bottom boundary) file
std::string ERF::nc_low_file; // Must provide via input

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
    int fix_random_seed = 0;
    ParmParse pp("erf"); pp.query("fix_random_seed", fix_random_seed);
    // Note that the value of 1024UL is not significant -- the point here is just to set the
    // same seed for all MPI processes for the purpose of regression testing
    if (fix_random_seed) {
        Print() << "Fixing the random seed" << std::endl;
        InitRandom(1024UL);
    }

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

    qheating_rates.resize(nlevs_max);
    rad_fluxes.resize(nlevs_max);
    sw_lw_fluxes.resize(nlevs_max);
    solar_zenith.resize(nlevs_max);

    // NOTE: size lsm before readparams (chooses the model at all levels)
    lsm.ReSize(nlevs_max);
    lsm_data.resize(nlevs_max);
    lsm_flux.resize(nlevs_max);

    // NOTE: size canopy model before readparams (if file exists, we construct)
    m_forest_drag.resize(nlevs_max);
    for (int lev = 0; lev <= max_level; ++lev) { m_forest_drag[lev] = nullptr;}

    ReadParameters();
    initializeMicrophysics(nlevs_max);

#ifdef ERF_USE_WINDFARM
    initializeWindFarm(nlevs_max);
#endif

#ifdef ERF_USE_SHOC
    shoc_interface.resize(nlevs_max);
    if (solverChoice.use_shoc) {
        for (int lev = 0; lev <= max_level; ++lev) {
            shoc_interface[lev] = std::make_unique<SHOCInterface>(lev, solverChoice);
        }
    }
#endif

    rad.resize(nlevs_max);
    for (int lev = 0; lev <= max_level; ++lev) {
        if (solverChoice.rad_type == RadiationType::RRTMGP) {
#ifdef ERF_USE_RRTMGP
            rad[lev] = std::make_unique<Radiation>(lev, solverChoice);
            // pass radiation datalog frequency to model - RRTMGP needs to know when to save data for profiles
            rad[lev]->setDataLogFrequency(rad_datalog_int);
#endif
        } else if (solverChoice.rad_type != RadiationType::None) {
            Abort("Don't know this radiation model!");
        }
    }

    const std::string& pv3d_1 = "plot_vars_1"  ; setPlotVariables(pv3d_1,plot3d_var_names_1);
    const std::string& pv3d_2 = "plot_vars_2"  ; setPlotVariables(pv3d_2,plot3d_var_names_2);
    const std::string& pv2d_1 = "plot2d_vars_1"; setPlotVariables2D(pv2d_1,plot2d_var_names_1);
    const std::string& pv2d_2 = "plot2d_vars_2"; setPlotVariables2D(pv2d_2,plot2d_var_names_2);

    // This is only used when we have mesh_type == MeshType::StretchedDz
    stretched_dz_h.resize(nlevs_max);
    stretched_dz_d.resize(nlevs_max);

    // Initialize staggered vertical levels for grid stretching or terrain, and
    // to simplify Rayleigh damping layer calculations.
    zlevels_stag.resize(max_level+1);
    init_zlevels(zlevels_stag,
                 stretched_dz_h,
                 stretched_dz_d,
                 geom,
                 refRatio(),
                 solverChoice.grid_stretching_ratio,
                 solverChoice.zsurf,
                 solverChoice.dz0);

    if (SolverChoice::mesh_type == MeshType::StretchedDz ||
        SolverChoice::mesh_type == MeshType::VariableDz) {
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

    // Get lo/hi indices for massflux calc
    if ((solverChoice.const_massflux_u != 0) || (solverChoice.const_massflux_v != 0)) {
        if (solverChoice.mesh_type == MeshType::ConstantDz) {
            const Real massflux_zlo = solverChoice.const_massflux_layer_lo - geom[0].ProbLo(2);
            const Real massflux_zhi = solverChoice.const_massflux_layer_hi - geom[0].ProbLo(2);
            const Real dz = geom[0].CellSize(2);
            if (massflux_zlo == -1e34) {
                solverChoice.massflux_klo = geom[0].Domain().smallEnd(2);
            } else {
                solverChoice.massflux_klo = static_cast<int>(std::ceil(massflux_zlo / dz - 0.5));
            }
            if (massflux_zhi ==  1e34) {
                solverChoice.massflux_khi = geom[0].Domain().bigEnd(2);
            } else {
                solverChoice.massflux_khi = static_cast<int>(std::floor(massflux_zhi / dz - 0.5));
            }
        } else if (solverChoice.mesh_type == MeshType::StretchedDz) {
            const Real massflux_zlo = solverChoice.const_massflux_layer_lo;
            const Real massflux_zhi = solverChoice.const_massflux_layer_hi;
            solverChoice.massflux_klo = geom[0].Domain().smallEnd(2);
            solverChoice.massflux_khi = geom[0].Domain().bigEnd(2) + 1;
            for (int k=0; k <= geom[0].Domain().bigEnd(2)+1; ++k) {
                if (zlevels_stag[0][k] <= massflux_zlo) solverChoice.massflux_klo = k;
                if (zlevels_stag[0][k] <= massflux_zhi) solverChoice.massflux_khi = k;
            }
        } else { // solverChoice.mesh_type == MeshType::VariableDz
            Error("Const massflux with variable dz not supported -- planar averages are on k rather than constant-z planes");
        }

        Print() << "Constant mass flux based on k in ["
            << solverChoice.massflux_klo << ", " << solverChoice.massflux_khi << "]" << std::endl;
    }

    prob = amrex_probinit(geom[0].ProbLo(),geom[0].ProbHi());

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, std::min(1.e100,dt_max_initial));
    dt_mri_ratio.resize(nlevs_max, 1);

    vars_new.resize(nlevs_max);
    vars_old.resize(nlevs_max);
    gradp.resize(nlevs_max);

    // We resize this regardless in order to pass it without error
    pp_inc.resize(nlevs_max);

    // Used in the fast substepping only
    lagged_delta_rt.resize(nlevs_max);
    avg_xmom.resize(nlevs_max);
    avg_ymom.resize(nlevs_max);
    avg_zmom.resize(nlevs_max);

    rU_new.resize(nlevs_max);
    rV_new.resize(nlevs_max);
    rW_new.resize(nlevs_max);

    rU_old.resize(nlevs_max);
    rV_old.resize(nlevs_max);
    rW_old.resize(nlevs_max);

    // xmom_crse_rhs.resize(nlevs_max);
    // ymom_crse_rhs.resize(nlevs_max);
    zmom_crse_rhs.resize(nlevs_max);

    for (int lev = 0; lev < nlevs_max; ++lev) {
        vars_new[lev].resize(Vars::NumTypes);
        vars_old[lev].resize(Vars::NumTypes);
        gradp[lev].resize(AMREX_SPACEDIM);
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
    th_bc_data.resize(nlevs_max);

    advflux_reg.resize(nlevs_max);

    // Stresses
    Tau.resize(nlevs_max);
    Tau_corr.resize(nlevs_max);
    SFS_hfx1_lev.resize(nlevs_max); SFS_hfx2_lev.resize(nlevs_max); SFS_hfx3_lev.resize(nlevs_max);
    SFS_diss_lev.resize(nlevs_max);
    SFS_q1fx1_lev.resize(nlevs_max); SFS_q1fx2_lev.resize(nlevs_max); SFS_q1fx3_lev.resize(nlevs_max);
    SFS_q2fx3_lev.resize(nlevs_max);
    eddyDiffs_lev.resize(nlevs_max);
    SmnSmn_lev.resize(nlevs_max);

    // Sea surface temps
    sst_lev.resize(nlevs_max);
    tsk_lev.resize(nlevs_max);
    lmask_lev.resize(nlevs_max);

    // Land and soil grid type and urban fractions
    land_type_lev.resize(nlevs_max);
    soil_type_lev.resize(nlevs_max);
    urb_frac_lev.resize(nlevs_max);

    // Metric terms
    z_phys_nd.resize(nlevs_max);
    z_phys_cc.resize(nlevs_max);
    detJ_cc.resize(nlevs_max);
    ax.resize(nlevs_max);
    ay.resize(nlevs_max);
    az.resize(nlevs_max);

    z_phys_nd_new.resize(nlevs_max);
    detJ_cc_new.resize(nlevs_max);

    z_phys_nd_src.resize(nlevs_max);
    z_phys_cc_src.resize(nlevs_max);
    detJ_cc_src.resize(nlevs_max);
    ax_src.resize(nlevs_max);
    ay_src.resize(nlevs_max);
    az_src.resize(nlevs_max);

    z_t_rk.resize(nlevs_max);

    terrain_blanking.resize(nlevs_max);

    // Wall distance
    walldist.resize(nlevs_max);

    // BoxArrays to make MultiFabs needed to convert WRFBdy data
    ba1d.resize(nlevs_max);
    ba2d.resize(nlevs_max);

    // MultiFabs needed to convert WRFBdy data
    mf_PSFC.resize(nlevs_max);

    // Map factors
    mapfac.resize(nlevs_max);

    // Fine mask
    fine_mask.resize(nlevs_max);

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

    // Size lat long arrays and default to null pointers
    lat_m.resize(nlevs_max);
    lon_m.resize(nlevs_max);
    for (int lev = 0; lev < max_level; ++lev) {
        lat_m[lev] = nullptr;
        lon_m[lev] = nullptr;
    }

    // Variable coriolis
    sinPhi_m.resize(nlevs_max);
    cosPhi_m.resize(nlevs_max);
    for (int lev = 0; lev < max_level; ++lev) {
        sinPhi_m[lev] = nullptr;
        cosPhi_m[lev] = nullptr;
    }

    // Initialize tagging criteria for mesh refinement
    refinement_criteria_setup();

    for (int lev = 0; lev < max_level; ++lev)
    {
       Print() << "Refinement ratio at level " << lev+1 << " set to be " <<
          ref_ratio[lev][0]  << " " << ref_ratio[lev][1]  <<  " " << ref_ratio[lev][2] << std::endl;
    }

    // We will create each of these in MakeNewLevelFromScratch
    eb.resize(max_level+1);
    for (int lev = 0; lev < max_level + 1; lev++){
        eb[lev] = std::make_unique<eb_>();
    }

    //
    // Construct the EB data structures and store in a separate class
    //
    // This is needed before initializing level MultiFabs
    if ( solverChoice.terrain_type == TerrainType::EB ||
         solverChoice.terrain_type == TerrainType::ImmersedForcing)
    {
        std::string geometry ="terrain";
        ParmParse pp("eb2");
        pp.queryAdd("geometry", geometry);

        constexpr int ngrow_for_eb = 4;  // This is the default in amrex but we need to explicitly pass it here since
                               // we want to also pass the build_coarse_level_by_coarsening argument
        const bool build_eb_for_multigrid = (solverChoice.terrain_type == TerrainType::EB &&
                                            ((solverChoice.project_initial_velocity[0] == 1) ||
                                            solverChoice.anelastic[0] == 1));
        // Note this just needs to be an integer > number of V-cycles one might use
        const int max_coarsening_level = (build_eb_for_multigrid) ? 100 : 0;
        const bool build_coarse_level_by_coarsening(false);

        // Define GeometryShop using the implicit function
        if (geometry == "terrain") {
            Box terrain_bx(surroundingNodes(geom[max_level].Domain())); terrain_bx.grow(3);
            FArrayBox terrain_fab(makeSlab(terrain_bx,2,0),1);
            Real dummy_time = 0.0;
            prob->init_terrain_surface(geom[max_level], terrain_fab, dummy_time);
            TerrainIF implicit_fun(terrain_fab, geom[max_level], stretched_dz_d[max_level]);
            auto gshop = EB2::makeShop(implicit_fun);
            if (build_eb_for_multigrid) {
                EB2::Build(gshop, geom[max_level], max_level, max_coarsening_level,
                            ngrow_for_eb, build_coarse_level_by_coarsening);
            } else {
                EB2::Build(gshop, this->Geom(), ngrow_for_eb);
            }
        } else if (geometry == "box") {
            RealArray box_lo{0.0, 0.0, 0.0};
            RealArray box_hi{0.0, 0.0, 0.0};
            pp.query("box_lo", box_lo);
            pp.query("box_hi", box_hi);
            EB2::BoxIF implicit_fun(box_lo, box_hi, false);
            auto gshop = EB2::makeShop(implicit_fun);
            if (build_eb_for_multigrid) {
                EB2::Build(gshop, geom[max_level], max_level, max_coarsening_level,
                            ngrow_for_eb, build_coarse_level_by_coarsening);
            } else {
                EB2::Build(gshop, this->Geom(), ngrow_for_eb);
            }
        } else if (geometry == "sphere") {
            auto ProbLoArr = geom[max_level].ProbLoArray();
            auto ProbHiArr = geom[max_level].ProbHiArray();
            const Real xcen = 0.5 * (ProbLoArr[0] + ProbHiArr[0]);
            const Real ycen = 0.5 * (ProbLoArr[1] + ProbHiArr[1]);
            RealArray sphere_center = {xcen, ycen, 0.0};
            EB2::SphereIF implicit_fun(0.5, sphere_center, false);
            auto gshop = EB2::makeShop(implicit_fun);
            if (build_eb_for_multigrid) {
                EB2::Build(gshop, geom[max_level], max_level, max_coarsening_level,
                            ngrow_for_eb, build_coarse_level_by_coarsening);
            } else {
                EB2::Build(gshop, this->Geom(), ngrow_for_eb);
            }
        }
    }

    if ( solverChoice.buildings_type == BuildingsType::ImmersedForcing) {
        constexpr int ngrow_for_eb = 4;
        Box buildings_bx(surroundingNodes(geom[max_level].Domain())); buildings_bx.grow(3);
        FArrayBox buildings_fab(makeSlab(buildings_bx,2,0),1);
        Real dummy_time = 0.0;
        prob->init_buildings_surface(geom[max_level], buildings_fab, dummy_time);
        TerrainIF implicit_fun(buildings_fab, geom[max_level], stretched_dz_d[max_level]);
        auto gshop = EB2::makeShop(implicit_fun);
        EB2::Build(gshop, this->Geom(), ngrow_for_eb);
    }
    forecast_state_1.resize(nlevs_max);
    forecast_state_2.resize(nlevs_max);
    forecast_state_interp.resize(nlevs_max);
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
    for (int step = istep[0]; step < max_step && start_time+cur_time < stop_time; ++step)
    {
        if (use_datetime) {
            Print() << "\n" << getTimestamp(start_time+cur_time, datetime_format)
                    << " (" << cur_time << " s elapsed)" << std::endl;
        }
        Print() << "\nCoarse STEP " << step+1 << " starts ..." << std::endl;

        ComputeDt(step);

        // Make sure we have read enough of the boundary plane data to make it through this timestep
        if (input_bndry_planes)
        {
            m_r2d->read_input_files(cur_time,dt[0],m_bc_extdir_vals);
        }

#ifdef ERF_USE_PARTICLES
        // We call this every time step with the knowledge that the particles may be
        //    initialized at a later time than the simulation start time.
        // The ParticleContainer carries a "start time" so the initialization will happen
        //    only when a) time > start_time, and b) particles have not yet been initialized
        initializeTracers((ParGDBBase*)GetParGDB(),z_phys_nd,cur_time);
#endif

        if(solverChoice.init_type == InitType::HindCast and
          solverChoice.hindcast_lateral_forcing) {
            for(int lev=0;lev<finest_level+1;lev++){
                WeatherDataInterpolation(lev,cur_time,z_phys_nd,false);
            }
        }

        auto dEvolveTime0 = amrex::second();

        int iteration = 1;
        timeStep(0, cur_time, iteration);

        cur_time  += dt[0];

        Print() << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time
                << " DT = " << dt[0]  << std::endl;

        if (check_for_nans > 0) {
            amrex::Print() << "Testing new state and vels for NaNs at end of timestep" << std::endl;
            for (int lev = 0; lev <= finest_level; ++lev) {
                check_state_for_nans(vars_new[lev][IntVars::cons]);
                check_vels_for_nans(vars_new[lev][Vars::xvel],vars_new[lev][Vars::yvel],vars_new[lev][Vars::zvel]);
            }
        }

        if (verbose > 0)
        {
            auto dEvolveTime = amrex::second() - dEvolveTime0;
            ParallelDescriptor::ReduceRealMax(dEvolveTime,ParallelDescriptor::IOProcessorNumber());
            amrex::Print() << "Timestep time = " << dEvolveTime << " seconds." << '\n';
        }

        post_timestep(step, cur_time, dt[0]);

        if (writeNow(cur_time, step+1, m_plot3d_int_1, m_plot3d_per_1, dt[0], last_plot3d_file_time_1)) {
            last_plot3d_file_step_1 = step+1;
            Write3DPlotFile(1,plotfile3d_type_1,plot3d_var_names_1);
            for (int lev = 0; lev <= finest_level; ++lev) {lsm.Plot(lev, step+1);}
            if (m_plot3d_per_1 > 0.) {last_plot3d_file_time_1 += m_plot3d_per_1;}
        }
        if (writeNow(cur_time, step+1, m_plot3d_int_2, m_plot3d_per_2, dt[0], last_plot3d_file_time_2)) {
            last_plot3d_file_step_2 = step+1;
            Write3DPlotFile(2,plotfile3d_type_2,plot3d_var_names_2);
            for (int lev = 0; lev <= finest_level; ++lev) {lsm.Plot(lev, step+1);}
            if (m_plot3d_per_2 > 0.) {last_plot3d_file_time_2 += m_plot3d_per_2;}
        }

        if (writeNow(cur_time, step+1, m_plot2d_int_1, m_plot2d_per_1, dt[0], last_plot2d_file_time_1)) {
            last_plot2d_file_step_1 = step+1;
            Write2DPlotFile(1,plotfile2d_type_1,plot2d_var_names_1);
            if (m_plot2d_per_1 > 0.) {last_plot2d_file_time_1 += m_plot2d_per_1;}
        }

        if (writeNow(cur_time, step+1, m_plot2d_int_2, m_plot2d_per_2, dt[0], last_plot2d_file_time_2)) {
            last_plot2d_file_step_2 = step+1;
            Write2DPlotFile(2,plotfile2d_type_2,plot2d_var_names_2);
            if (m_plot2d_per_2 > 0.) {last_plot2d_file_time_2 += m_plot2d_per_2;}
        }

        for (int i = 0; i < m_subvol_int.size(); i++) {
            if (writeNow(cur_time, step+1, m_subvol_int[i], m_subvol_per[i], dt[0], last_subvol_time[i])) {
                last_subvol_step[i] = step+1;
                WriteSubvolume(i,subvol3d_var_names);
                if (m_subvol_per[i] > 0.) {last_subvol_time[i] += m_subvol_per[i];}
            }
        }

        if (writeNow(cur_time, step+1, m_check_int, m_check_per, dt[0], last_check_file_time)) {
            last_check_file_step = step+1;
            WriteCheckpointFile();
            if (m_check_per > 0.) {last_check_file_time += m_check_per;}
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
    if ( (m_plot3d_int_1 > 0 || m_plot3d_per_1 > 0.) && istep[0] > last_plot3d_file_step_1 ) {
        Write3DPlotFile(1,plotfile3d_type_1,plot3d_var_names_1);
        if (m_plot3d_per_1 > 0.) {last_plot3d_file_time_1 += m_plot3d_per_1;}
    }
    if ( (m_plot3d_int_2 > 0 || m_plot3d_per_2 > 0.) && istep[0] > last_plot3d_file_step_2) {
        Write3DPlotFile(2,plotfile3d_type_1,plot3d_var_names_2);
        if (m_plot3d_per_2 > 0.) {last_plot3d_file_time_2 += m_plot3d_per_2;}
    }
    if ( (m_plot2d_int_1 > 0 || m_plot2d_per_1 > 0.) && istep[0] > last_plot2d_file_step_1 ) {
        Write2DPlotFile(1,plotfile2d_type_1,plot2d_var_names_1);
        if (m_plot2d_per_1 > 0.) {last_plot2d_file_time_1 += m_plot2d_per_1;}
    }
    if ( (m_plot2d_int_2 > 0 || m_plot2d_per_2 > 0.) && istep[0] > last_plot2d_file_step_2) {
        Write2DPlotFile(2,plotfile2d_type_1,plot2d_var_names_2);
        if (m_plot2d_per_2 > 0.) {last_plot2d_file_time_2 += m_plot2d_per_2;}
    }

    for (int i = 0; i < m_subvol_int.size(); i++) {
        if ( (m_subvol_int[i] > 0 || m_subvol_per[i] > 0.) && istep[0] > last_subvol_step[i]) {
            WriteSubvolume(i,subvol3d_var_names);
            if (m_subvol_per[i] > 0.) {last_subvol_time[i] += m_subvol_per[i];}
        }
    }

    if ( (m_check_int > 0 || m_check_per > 0.) && istep[0] > last_check_file_step) {
        WriteCheckpointFile();
        if (m_check_per > 0.) {last_check_file_time += m_check_per;}
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
        int ncomp = vars_new[0][Vars::cons].nComp();
        for (int lev = finest_level-1; lev >= 0; lev--)
        {
            // The quantity that is conserved is not (rho S), but rather (rho S / m^2) where
            // m is the map scale factor at cell centers
            // Here we pre-divide (rho S) by m^2 before refluxing
            for (MFIter mfi(vars_new[lev][Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                const Array4<      Real> cons_arr = vars_new[lev][Vars::cons].array(mfi);
                const Array4<const Real>  mfx_arr = mapfac[lev][MapFacType::m_x]->const_array(mfi);
                const Array4<const Real>  mfy_arr = mapfac[lev][MapFacType::m_y]->const_array(mfi);
                if (SolverChoice::mesh_type == MeshType::ConstantDz) {
                    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        cons_arr(i,j,k,n) /= (mfx_arr(i,j,0)*mfy_arr(i,j,0));
                    });
                } else {
                    const Array4<const Real>   detJ_arr = detJ_cc[lev]->const_array(mfi);
                    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        cons_arr(i,j,k,n) *= detJ_arr(i,j,k) / (mfx_arr(i,j,0)*mfy_arr(i,j,0));
                    });
                }
            } // mfi

            // This call refluxes all "slow" cell-centered variables
            // (i.e. not density or (rho theta) or velocities) from the lev/lev+1 interface onto lev
            getAdvFluxReg(lev+1)->Reflux(vars_new[lev][Vars::cons], 2, 2, ncomp-2);

            // Here we multiply (rho S) by m^2 after refluxing
            for (MFIter mfi(vars_new[lev][Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                const Array4<      Real>   cons_arr = vars_new[lev][Vars::cons].array(mfi);
                const Array4<const Real>  mfx_arr = mapfac[lev][MapFacType::m_x]->const_array(mfi);
                const Array4<const Real>  mfy_arr = mapfac[lev][MapFacType::m_y]->const_array(mfi);
                if (SolverChoice::mesh_type == MeshType::ConstantDz) {
                    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        cons_arr(i,j,k,n) *= (mfx_arr(i,j,0)*mfy_arr(i,j,0));
                    });
                } else {
                    const Array4<const Real>   detJ_arr = detJ_cc[lev]->const_array(mfi);
                    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        cons_arr(i,j,k,n) *= (mfx_arr(i,j,0)*mfy_arr(i,j,0)) / detJ_arr(i,j,k);
                    });
                }
            } // mfi

            // We need to do this before anything else because refluxing changes the
            // values of coarse cells underneath fine grids with the assumption they'll
            // be over-written by averaging down
            int src_comp;
            if (solverChoice.anelastic[lev]) {
                src_comp = 1;
            } else {
                src_comp = 0;
            }
            int num_comp = ncomp - src_comp;
            AverageDownTo(lev,src_comp,num_comp);
        }
    }

    if (is_it_time_for_action(nstep, time, dt_lev0, sum_interval, sum_per)) {
        sum_integrated_quantities(time);
        sum_derived_quantities(time);
        sum_energy_quantities(time);
    }

    if (solverChoice.pert_type == PerturbationType::Source ||
        solverChoice.pert_type == PerturbationType::Direct ||
        solverChoice.pert_type == PerturbationType::CPM) {
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

    if (solverChoice.rad_type != RadiationType::None)
    {
        if ( rad_datalog_int > 0 &&
             (((nstep+1) % rad_datalog_int == 0) || (nstep==0)) ) {
            if (rad[0]->hasDatalog()) {
                rad[0]->WriteDataLog(time+start_time);
            }
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
         bool is_moist = (micro->Get_Qstate_Moist_Size() > 0);
         m_w2d->write_planes(istep[0], time, vars_new, is_moist);
      }
    }

    // Write plane/line sampler data
    if (line_sampler && is_it_time_for_action(nstep+1, time, dt_lev0, line_sampling_interval, line_sampling_per)) {
        line_sampler->get_sample_data(geom, vars_new);
        line_sampler->write_sample_data(t_new, istep, ref_ratio, geom);
    }
    if (plane_sampler && is_it_time_for_action(nstep+1, time, dt_lev0, plane_sampling_interval, plane_sampling_per)) {
        plane_sampler->get_sample_data(geom, vars_new);
        plane_sampler->write_sample_data(t_new, istep, ref_ratio, geom);
    }

    // Moving terrain
    if ( solverChoice.terrain_type == TerrainType::MovingFittedMesh )
    {
      for (int lev = finest_level; lev >= 0; lev--)
      {
        // Copy z_phs_nd and detJ_cc at end of timestep
        MultiFab::Copy(*z_phys_nd[lev], *z_phys_nd_new[lev], 0, 0, 1, z_phys_nd[lev]->nGrowVect());
        MultiFab::Copy(  *detJ_cc[lev],   *detJ_cc_new[lev], 0, 0, 1,   detJ_cc[lev]->nGrowVect());
        MultiFab::Copy(base_state[lev],base_state_new[lev],0,0,BaseState::num_comps,base_state[lev].nGrowVect());

        make_zcc(geom[lev],*z_phys_nd[lev],*z_phys_cc[lev]);
      }
    }

    if ( solverChoice.io_hurricane_eye_tracker and (nstep == 0 or (nstep+1)%m_plot3d_int_1 == 0) )
    {
        int levc=finest_level;

        HurricaneEyeTracker(geom[levc],
                            vars_new[levc],
                            solverChoice.moisture_type,
                            solverChoice.hindcast_lateral_forcing? &forecast_state_interp[levc] : nullptr,
                            solverChoice.hurricane_eye_latitude,
                            solverChoice.hurricane_eye_longitude,
                            hurricane_eye_track_xy,
                            hurricane_eye_track_latlon,
                            hurricane_tracker_circle);

        MultiFab& U_new = vars_new[levc][Vars::xvel];
        MultiFab& V_new = vars_new[levc][Vars::yvel];
        MultiFab& W_new = vars_new[levc][Vars::zvel];

        MultiFab mf_cc_vel(grids[levc], dmap[levc], AMREX_SPACEDIM, IntVect(0,0,0));
        average_face_to_cellcenter(mf_cc_vel,0,{AMREX_D_DECL(&U_new,&V_new,&W_new)},0);

        HurricaneMaxVelTracker(geom[levc],
                               mf_cc_vel,
                               t_new[0],
                               hurricane_eye_track_xy,
                               hurricane_maxvel_vs_time);

        std::string filename_tracker = MakeVTKFilename_TrackerCircle(nstep);
        std::string filename_xy      = MakeVTKFilename_EyeTracker_xy(nstep);
        std::string filename_latlon  = MakeFilename_EyeTracker_latlon(nstep);
        std::string filename_maxvel  = MakeFilename_EyeTracker_maxvel(nstep);
        if (ParallelDescriptor::IOProcessor()) {
            WriteVTKPolyline(filename_tracker, hurricane_tracker_circle);
            WriteVTKPolyline(filename_xy, hurricane_eye_track_xy);
            WriteLinePlot(filename_latlon, hurricane_eye_track_latlon);
            WriteLinePlot(filename_maxvel, hurricane_maxvel_vs_time);
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

    if (restart_chkfile.empty()) {
        // Start simulation from the beginning
        InitFromScratch(0.0);
    } else {
        // For initialization this is done in init_only; it is done here for restart
        init_bcs();
    }

    solverChoice.check_params(max_level,geom,phys_bc_type);
}

void
ERF::InitData_post ()
{
    if (solverChoice.advChoice.have_zero_flux_faces)
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(finest_level == 0,
            "Thin immersed body with refinement not currently supported.");
        if (SolverChoice::mesh_type != MeshType::ConstantDz) {
            amrex::Print() << "NOTE: Thin immersed body with non-constant dz has not been tested." << std::endl;
        }
    }

    if (!restart_chkfile.empty()) {
        restart();
    }
    //
    // Make sure that detJ and z_phys_cc are the average of the data on a finer level if there is one and if two way coupling
    //
    if (SolverChoice::mesh_type != MeshType::ConstantDz) {
        if (solverChoice.coupling_type == CouplingType::TwoWay) {
            for (int crse_lev = finest_level-1; crse_lev >= 0; crse_lev--) {
                average_down(  *detJ_cc[crse_lev+1],   *detJ_cc[crse_lev], 0, 1, refRatio(crse_lev));
                average_down(*z_phys_cc[crse_lev+1], *z_phys_cc[crse_lev], 0, 1, refRatio(crse_lev));
            }
        }
        for (int crse_lev = finest_level-1; crse_lev >= 0; crse_lev--) {
              detJ_cc[crse_lev]->FillBoundary(geom[crse_lev].periodicity());
            z_phys_cc[crse_lev]->FillBoundary(geom[crse_lev].periodicity());
        }
    }

#ifdef ERF_IMPLICIT_W
    if (SolverChoice::mesh_type == MeshType::VariableDz &&
        (solverChoice.vert_implicit_fac[0] > 0 ||
         solverChoice.vert_implicit_fac[1] > 0 ||
         solverChoice.vert_implicit_fac[2] > 0  )       &&
        solverChoice.implicit_momentum_diffusion)
    {
        Warning("Doing implicit solve for u, v, and w with terrain -- this has not been tested");
    }
#endif

    //
    // Copy vars_new into vars_old, then use vars_old to fill covered cells in vars_new during AverageDown
    //
    if (SolverChoice::terrain_type == TerrainType::EB) {
        for (int lev = 0; lev <= finest_level; lev++) {
            int ncomp_cons = vars_new[lev][Vars::cons].nComp();
            MultiFab::Copy(vars_old[lev][Vars::cons],vars_new[lev][Vars::cons],0,0,ncomp_cons,vars_new[lev][Vars::cons].nGrowVect());
        }
    }

    if (restart_chkfile.empty()) {
        if (solverChoice.coupling_type == CouplingType::TwoWay) {
            AverageDown();
        }
    }

#ifdef ERF_USE_PARTICLES
    if (restart_chkfile.empty()) {
        if (Microphysics::modelType(solverChoice.moisture_type) == MoistureModelType::Lagrangian) {
            if (solverChoice.moisture_tight_coupling) {
                Warning("Tight coupling has not been tested with Lagrangian microphysics");
            }

            for (int lev = 0; lev <= finest_level; lev++) {
                dynamic_cast<LagrangianMicrophysics&>(*micro).initParticles(z_phys_nd[lev]);
            }
        }
    }
#endif

    if (!restart_chkfile.empty()) { // Restart from a checkpoint

        // Create the physbc objects for {cons, u, v, w, base state}
        // We fill the additional base state ghost cells just in case we have read the old format
        for (int lev(0); lev <= finest_level; ++lev) {
            make_physbcs(lev);
            (*physbcs_base[lev])(base_state[lev],0,base_state[lev].nComp(),base_state[lev].nGrowVect());
        }

        if (solverChoice.do_forest_drag) {
            for (int lev(0); lev <= finest_level; ++lev) {
                m_forest_drag[lev]->define_drag_field(grids[lev], dmap[lev], geom[lev],
                                                      z_phys_cc[lev].get(), z_phys_nd[lev].get());
            }
        }

#ifdef ERF_USE_NETCDF
        //
        // Create the needed bdy_data_xlo etc ... since we don't read it in from checkpoint any more
        // This follows init_from_wrfinput()
        //
        bool use_moist = (solverChoice.moisture_type != MoistureType::None);
        if (solverChoice.use_real_bcs) {

            if ( geom[0].isPeriodic(0) || geom[0].isPeriodic(1) ) {
                 amrex::Error("Cannot set periodic lateral boundary conditions when reading in real boundary values");
            }

            bdy_time_interval = read_times_from_wrfbdy(nc_bdy_file,
                                                       bdy_data_xlo,bdy_data_xhi,bdy_data_ylo,bdy_data_yhi,
                                                       start_bdy_time);
            Real dT = bdy_time_interval;

            int n_time_old = static_cast<int>(t_new[0] /  dT);

            int lev = 0;

            int ntimes = std::min(n_time_old+3, static_cast<int>(bdy_data_xlo.size()));

            for (int itime = n_time_old; itime < ntimes; itime++)
            {
                read_from_wrfbdy(itime,nc_bdy_file,geom[0].Domain(),
                                 bdy_data_xlo,bdy_data_xhi,bdy_data_ylo,bdy_data_yhi,
                                 real_width);
                convert_all_wrfbdy_data(itime, geom[0].Domain(), bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi,
                                        *mf_MUB, *mf_C1H, *mf_C2H,
                                        vars_new[lev][Vars::xvel], vars_new[lev][Vars::yvel], vars_new[lev][Vars::cons],
                                        geom[lev], use_moist);
            } // itime
        } // use_real_bcs

        if (!nc_low_file.empty())
        {
            low_time_interval = read_times_from_wrflow(nc_low_file,
                                                       low_data_zlo,
                                                       start_low_time);
            Real dT = low_time_interval;

            int lev = 0;
            sst_lev[lev].resize(low_data_zlo.size());
            tsk_lev[lev].resize(low_data_zlo.size());

            int n_time_old = static_cast<int>(t_new[0] /  dT);

            int ntimes = std::min(n_time_old+2, static_cast<int>(low_data_zlo.size()));

            for (int itime = n_time_old; itime < ntimes; itime++)
            {
                read_from_wrflow(itime, nc_low_file, geom[lev].Domain(), low_data_zlo);

                // Need to read PSFC
                FArrayBox NC_fab_var_file;
                for (int idx = 0; idx < num_boxes_at_level[lev]; idx++) {
                    int success, use_theta_m;
                    read_from_wrfinput(lev, boxes_at_level[lev][idx], nc_init_file[lev][0],
                                       NC_fab_var_file, "PSFC", geom[lev],
                                       use_theta_m, success);
                    auto& var_fab = NC_fab_var_file;
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                    for ( MFIter mfi(*mf_PSFC[lev], false); mfi.isValid(); ++mfi )
                    {
                        FArrayBox &cur_fab = (*mf_PSFC[lev])[mfi];
                        cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
                    }
                    var_fab.clear();
                }

                update_sst_tsk(itime, geom[lev], ba2d[lev],
                               sst_lev[lev], tsk_lev[lev],
                               m_SurfaceLayer, low_data_zlo,
                               vars_new[lev][Vars::cons], *mf_PSFC[lev],
                               solverChoice.rdOcp, lmask_lev[lev][0], use_moist);
            } // itime
        }
#endif
    } // end restart

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
                if (SolverChoice::mesh_type == MeshType::VariableDz) {
                    amrex::Print() << "Note: 1-D geostrophic wind profile input is not defined for real terrain" << std::endl;
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
                                      geom[lev], z_phys_nd[lev]);
        }
    }

    if (solverChoice.dampingChoice.rayleigh_damp_U ||solverChoice.dampingChoice.rayleigh_damp_V ||
        solverChoice.dampingChoice.rayleigh_damp_W ||solverChoice.dampingChoice.rayleigh_damp_T)
    {
        initRayleigh();
        if (solverChoice.init_type == InitType::Input_Sounding)
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

    if (solverChoice.pert_type == PerturbationType::Source ||
        solverChoice.pert_type == PerturbationType::Direct ||
        solverChoice.pert_type == PerturbationType::CPM) {
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
            bool is_moist = (micro->Get_Qstate_Moist_Size() > 0);
            m_w2d->write_planes(0, time, vars_new, is_moist);
        }
    }

    // Fill boundary conditions in vars_new
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        auto& lev_new = vars_new[lev];

        // ***************************************************************************
        // Physical bc's at domain boundary
        // ***************************************************************************
        IntVect ngvect_cons = vars_new[lev][Vars::cons].nGrowVect();
        IntVect ngvect_vels = vars_new[lev][Vars::xvel].nGrowVect();

        int ncomp_cons = lev_new[Vars::cons].nComp();
        bool do_fb     = true;

#ifdef ERF_USE_NETCDF
        // We call this here because it is an ERF routine
        if (solverChoice.use_real_bcs && (lev==0)) {
            int icomp_cons = 0;
            bool cons_only = false;
            Vector<MultiFab*> mfs_vec = {&lev_new[Vars::cons],&lev_new[Vars::xvel],
                                         &lev_new[Vars::yvel],&lev_new[Vars::zvel]};
            if (solverChoice.upwind_real_bcs) {
                fill_from_realbdy_upwind(mfs_vec,t_new[lev],cons_only,icomp_cons,
                                         ncomp_cons,ngvect_cons,ngvect_vels);
            } else {
                fill_from_realbdy(mfs_vec,t_new[lev],cons_only,icomp_cons,
                                  ncomp_cons,ngvect_cons,ngvect_vels);
            }
            do_fb = false;
    }
#endif

        (*physbcs_cons[lev])(lev_new[Vars::cons],lev_new[Vars::xvel],lev_new[Vars::yvel],0,ncomp_cons,
                             ngvect_cons,t_new[lev],BCVars::cons_bc,do_fb);
        (   *physbcs_u[lev])(lev_new[Vars::xvel],lev_new[Vars::xvel],lev_new[Vars::yvel],
                             ngvect_vels,t_new[lev],BCVars::xvel_bc,do_fb);
        (   *physbcs_v[lev])(lev_new[Vars::yvel],lev_new[Vars::xvel],lev_new[Vars::yvel],
                             ngvect_vels,t_new[lev],BCVars::yvel_bc,do_fb);
        (   *physbcs_w[lev])(lev_new[Vars::zvel],lev_new[Vars::xvel],lev_new[Vars::yvel],
                             ngvect_vels,t_new[lev],BCVars::zvel_bc,do_fb);
    }

    //
    // If we are starting from scratch, we have the option to project the initial velocity field
    //    regardless of how we initialized.  Note that project_initial_velocity operates on vars_new.
    // pp_inc is used as scratch space here; we zero it out after the projection
    //
    if (restart_chkfile == "")
    {
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            if (solverChoice.project_initial_velocity[lev] == 1) {
                Real dummy_dt = 1.0;
                if (verbose > 0) {
                    amrex::Print() << "Projecting initial velocity field at level " << lev << std::endl;
                }

                project_initial_velocity(lev, t_new[lev], dummy_dt);

                pp_inc[lev].setVal(0.);
                gradp[lev][GpVars::gpx].setVal(0.);
                gradp[lev][GpVars::gpy].setVal(0.);
                gradp[lev][GpVars::gpz].setVal(0.);
            } // project
        } // lev
    }

    // Copy from new into old just in case (after filling boundary conditions and possibly projecting)
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        int nc = vars_new[lev][Vars::cons].nComp();

        MultiFab::Copy(vars_old[lev][Vars::cons],vars_new[lev][Vars::cons],0,0,nc,vars_new[lev][Vars::cons].nGrowVect());
        MultiFab::Copy(vars_old[lev][Vars::xvel],vars_new[lev][Vars::xvel],0,0, 1,vars_new[lev][Vars::xvel].nGrowVect());
        MultiFab::Copy(vars_old[lev][Vars::yvel],vars_new[lev][Vars::yvel],0,0, 1,vars_new[lev][Vars::yvel].nGrowVect());
        MultiFab::Copy(vars_old[lev][Vars::zvel],vars_new[lev][Vars::zvel],0,0, 1,vars_new[lev][Vars::zvel].nGrowVect());
    }

    // Compute the minimum dz in the domain at each level (to be used for setting the timestep)
    dz_min.resize(max_level+1);
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        dz_min[lev] = geom[lev].CellSize(2);
        if ( SolverChoice::mesh_type != MeshType::ConstantDz ) {
            dz_min[lev] *= (*detJ_cc[lev]).min(0);
        }
    }

    // We don't need to recompute dt[lev] on restart because we read it in from the checkpoint file.
    if (restart_chkfile.empty()) {
        ComputeDt();
    }

    // Check the viscous limit
    DiffChoice dc = solverChoice.diffChoice;
    if (dc.molec_diff_type == MolecDiffType::Constant ||
        dc.molec_diff_type == MolecDiffType::ConstantAlpha) {
        Real delta = std::min({geom[finest_level].CellSize(0),
                               geom[finest_level].CellSize(1),
                               dz_min[finest_level]});
        if (dc.dynamic_viscosity == 0) {
            Print() << "Note: Molecular diffusion specified but dynamic_viscosity has not been specified" << std::endl;
        } else {
            Real nu = dc.dynamic_viscosity / dc.rho0_trans;
            Real viscous_limit = 2.0 * delta*delta / nu;
            Print() << "smallest grid spacing at level " << finest_level << " = " << delta << std::endl;
            Print() << "dt at level " << finest_level << " = " << dt[finest_level] << std::endl;
            Print() << "Viscous CFL is " << dt[finest_level] / viscous_limit << std::endl;
            if (fixed_dt[finest_level] >= viscous_limit) {
                Warning("Specified fixed_dt is above the viscous limit");
            } else if (dt[finest_level] >= viscous_limit) {
                Warning("Adaptive dt based on convective CFL only is above the viscous limit");
            }
        }
    }

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
        if (lev == 0) {
            FillPatchCrseLevel(lev, t_new[lev],
                      {&lev_new[Vars::cons],&lev_new[Vars::xvel],&lev_new[Vars::yvel],&lev_new[Vars::zvel]});
        } else {
            FillPatchFineLevel(lev, t_new[lev],
                      {&lev_new[Vars::cons],&lev_new[Vars::xvel],&lev_new[Vars::yvel],&lev_new[Vars::zvel]},
                      {&lev_new[Vars::cons],&rU_new[lev],&rV_new[lev],&rW_new[lev]},
                      base_state[lev], base_state[lev],
                      fillset);
        }

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
        if (solverChoice.terrain_type == TerrainType::MovingFittedMesh) {
            MultiFab::Copy(base_state_new[lev],base_state[lev],0,0,BaseState::num_comps,base_state[lev].nGrowVect());
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

    // If lev > 0, we need to fill bc's by interpolation from coarser grid
    for (int lev = 1; lev <= finest_level; ++lev)
    {
        Interp2DArrays(lev,grids[lev],dmap[lev]);
    } // lev

#ifdef ERF_USE_WW3_COUPLING
    int my_lev = 0;
    amrex::Print() <<  " About to call send_to_ww3 from ERF.cpp" << std::endl;
    send_to_ww3(my_lev);
    amrex::Print() <<  " About to call read_waves from ERF.cpp"  << std::endl;
    read_waves(my_lev);
   // send_to_ww3(my_lev);
#endif

    // Configure SurfaceLayer params if used
    // NOTE: we must set up the MOST routine after calling FillPatch
    //       in order to have lateral ghost cells filled (MOST + terrain interp).
    if (phys_bc_type[Orientation(Direction::z,Orientation::low)] == ERF_BC::surface_layer)
    {
        bool has_diff = ( (solverChoice.diffChoice.molec_diff_type != MolecDiffType::None) ||
                          (solverChoice.turbChoice[0].les_type  != LESType::None)          ||
                          (solverChoice.turbChoice[0].rans_type != RANSType::None)         ||
                          (solverChoice.turbChoice[0].pbl_type  != PBLType::None) );
        AMREX_ALWAYS_ASSERT(has_diff);

        bool rotate = solverChoice.use_rotate_surface_flux;
        if (rotate) {
            Print() << "Using surface layer model with stress rotations" << std::endl;
        }

        //
        // This constructor will make the SurfaceLayer object but not allocate the arrays at each level.
        //
        m_SurfaceLayer = std::make_unique<SurfaceLayer>(geom, rotate, pp_prefix, Qv_prim,
                                                        z_phys_nd,
                                                        solverChoice.mesh_type,
                                                        solverChoice.terrain_type,
                                                        start_time, stop_time
#ifdef ERF_USE_NETCDF
                                                        , bdy_time_interval
#endif
                                                        );
        // This call will allocate the arrays at each level. If we regrid later, either changing
        // the number of levels or just the grids at each existing level, we will call an update routine
        // to redefine the internal arrays in m_SurfaceLayer.
        for (int lev = 0; lev <= finest_level; lev++)
        {
            Vector<MultiFab*> mfv_old = {&vars_old[lev][Vars::cons], &vars_old[lev][Vars::xvel],
                                         &vars_old[lev][Vars::yvel], &vars_old[lev][Vars::zvel]};
            m_SurfaceLayer->make_SurfaceLayer_at_level(lev,finest_level+1,
                                                       mfv_old, Theta_prim[lev], Qv_prim[lev],
                                                       Qr_prim[lev], z_phys_nd[lev],
                                                       Hwave[lev].get(),Lwave[lev].get(),eddyDiffs_lev[lev].get(),
                                                       lsm_data[lev], lsm_data_name, lsm_flux[lev], lsm_flux_name,
                                                       sst_lev[lev], tsk_lev[lev], lmask_lev[lev]);
        }

        // If initializing from an input_sounding, make sure the surface layer
        // is using the same surface conditions
        if (solverChoice.init_type == InitType::Input_Sounding) {
            const Real theta0 = input_sounding_data.theta_ref_inp_sound;
            const Real qv0    = input_sounding_data.qv_ref_inp_sound;
            for (int lev = 0; lev <= finest_level; lev++) {
                m_SurfaceLayer->set_t_surf(lev, theta0);
                m_SurfaceLayer->set_q_surf(lev, qv0);
            }
        }

        if (restart_chkfile != "") {
            // Update surface fields if needed (and available)
            ReadCheckpointFileSurfaceLayer();
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

                int rhoqr_comp = solverChoice.moisture_indices.qr;
                if (rhoqr_comp > -1) {
                    MultiFab::Copy(  *Qr_prim[lev], vars_new[lev][Vars::cons], rhoqr_comp, 0, 1, ng);
                    MultiFab::Divide(*Qr_prim[lev], vars_new[lev][Vars::cons],   Rho_comp, 0, 1, ng);
                } else {
                    Qr_prim[lev]->setVal(0.0);
                }
            }
            m_SurfaceLayer->update_mac_ptrs(lev, vars_new, Theta_prim, Qv_prim, Qr_prim);

            if (restart_chkfile == "") {
                // Only do this if starting from scratch; if restarting, then
                // we don't want to call update_fluxes multiple times because
                // it will change u* and theta* from their previous values
                m_SurfaceLayer->update_pblh(lev, vars_new, z_phys_cc[lev].get(),
                                            solverChoice.moisture_indices);
                m_SurfaceLayer->update_fluxes(lev, time, vars_new[lev][Vars::cons], z_phys_nd[lev]);

                // Initialize tke(x,y,z) as a function of u*(x,y)
                if (solverChoice.turbChoice[lev].init_tke_from_ustar) {
                    Real qkefac = 1.0;
                    if (solverChoice.turbChoice[lev].pbl_type == PBLType::MYNN25 ||
                        solverChoice.turbChoice[lev].pbl_type == PBLType::MYNNEDMF)
                    {
                        // https://github.com/NCAR/MYNN-EDMF/blob/90f36c25259ec1960b24325f5b29ac7c5adeac73/module_bl_mynnedmf.F90#L1325-L1333
                        const Real B1 = solverChoice.turbChoice[lev].pbl_mynn.B1;
                        qkefac = 1.5 * std::pow(B1, 2.0/3.0);
                    }
                    m_SurfaceLayer->init_tke_from_ustar(lev, vars_new[lev][Vars::cons], z_phys_nd[lev], qkefac);
                }
            }
        }
    } // end if (phys_bc_type[Orientation(Direction::z,Orientation::low)] == ERF_BC::surface_layer)

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
    const std::string& pv3d_1 = "plot_vars_1"  ; appendPlotVariables(pv3d_1,plot3d_var_names_1);
    const std::string& pv3d_2 = "plot_vars_2"  ; appendPlotVariables(pv3d_2,plot3d_var_names_2);
    const std::string& pv2d_1 = "plot2d_vars_1"; appendPlotVariables(pv2d_1,plot2d_var_names_1);
    const std::string& pv2d_2 = "plot2d_vars_2"; appendPlotVariables(pv2d_2,plot2d_var_names_2);

    if ( restart_chkfile.empty() && (m_check_int > 0 || m_check_per > 0.) )
    {
        WriteCheckpointFile();
        last_check_file_step = 0;
        if (m_check_per > 0.) {last_check_file_time += m_check_per;}
    }

    if ( (restart_chkfile.empty()) ||
         (!restart_chkfile.empty() && plot_file_on_restart) )
    {
        if (m_plot3d_int_1 > 0 || m_plot3d_per_1 > 0.)
        {
            Write3DPlotFile(1,plotfile3d_type_1,plot3d_var_names_1);
            if (m_plot3d_per_1 > 0.) {last_plot3d_file_time_1 += m_plot3d_per_1;}
            last_plot3d_file_step_1 = istep[0];
        }
        if (m_plot3d_int_2 > 0 || m_plot3d_per_2 > 0.)
        {
            Write3DPlotFile(2,plotfile3d_type_2,plot3d_var_names_2);
            if (m_plot3d_per_2 > 0.) {last_plot3d_file_time_2 += m_plot3d_per_2;}
            last_plot3d_file_step_2 = istep[0];
        }
        if (m_plot2d_int_1 > 0 || m_plot2d_per_1 > 0.)
        {
            Write2DPlotFile(1,plotfile2d_type_1,plot2d_var_names_1);
            if (m_plot2d_per_1 > 0.) {last_plot2d_file_time_1 += m_plot2d_per_1;}
            last_plot2d_file_step_1 = istep[0];
        }
        if (m_plot2d_int_2 > 0 || m_plot2d_per_2 > 0.)
        {
            Write2DPlotFile(2,plotfile2d_type_2,plot2d_var_names_2);
            if (m_plot2d_per_2 > 0.) {last_plot2d_file_time_2 += m_plot2d_per_2;}
            last_plot2d_file_step_2 = istep[0];
        }
        for (int i = 0; i < m_subvol_int.size(); i++) {
            if (m_subvol_int[i] > 0 || m_subvol_per[i] > 0.) {
                WriteSubvolume(i,subvol3d_var_names);
                last_subvol_step[i] = istep[0];
                if (m_subvol_per[i] > 0.) {last_subvol_time[i] += m_subvol_per[i];}
            }
        }
    }

    // Set these up here because we need to know which MPI rank "cell" is on...
    if (pp.contains("data_log"))
    {
        int num_datalogs = pp.countval("data_log");
        datalog.resize(num_datalogs);
        datalogname.resize(num_datalogs);
        pp.queryarr("data_log",datalogname,0,num_datalogs);
        for (int i = 0; i < num_datalogs; i++) {
            setRecordDataInfo(i,datalogname[i]);
        }
    }

    if (pp.contains("der_data_log"))
    {
        int num_der_datalogs = pp.countval("der_data_log");
        der_datalog.resize(num_der_datalogs);
        der_datalogname.resize(num_der_datalogs);
        pp.queryarr("der_data_log",der_datalogname,0,num_der_datalogs);
        for (int i = 0; i < num_der_datalogs; i++) {
            setRecordDerDataInfo(i,der_datalogname[i]);
        }
    }

    if (pp.contains("energy_data_log"))
    {
        int num_energy_datalogs = pp.countval("energy_data_log");
        tot_e_datalog.resize(num_energy_datalogs);
        tot_e_datalogname.resize(num_energy_datalogs);
        pp.queryarr("energy_data_log",tot_e_datalogname,0,num_energy_datalogs);
        for (int i = 0; i < num_energy_datalogs; i++) {
            setRecordEnergyDataInfo(i,tot_e_datalogname[i]);
        }
    }

    if (solverChoice.rad_type != RadiationType::None)
    {
        // Create data log for radiation model if requested
        rad[0]->setupDataLog();
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

    if (is_it_time_for_action(istep[0], t_new[0], dt[0], sum_interval, sum_per)) {
        sum_integrated_quantities(t_new[0]);
        sum_derived_quantities(t_new[0]);
        sum_energy_quantities(t_new[0]);
    }

    // Create object to do line and plane sampling if needed
    bool do_line = false; bool do_plane = false;
    pp.query("do_line_sampling",do_line); pp.query("do_plane_sampling",do_plane);
    if (do_line) {
        if (line_sampling_interval < 0 && line_sampling_per < 0) {
            Abort("Need to specify line_sampling_interval or line_sampling_per");
        }
        line_sampler = std::make_unique<LineSampler>();
        line_sampler->write_coords(z_phys_cc);
    }
    if (do_plane) {
        if (plane_sampling_interval < 0 && plane_sampling_per < 0) {
            Abort("Need to specify plane_sampling_interval or plane_sampling_per");
        }
        plane_sampler = std::make_unique<PlaneSampler>();
    }

    if ( solverChoice.terrain_type == TerrainType::EB ||
         solverChoice.terrain_type == TerrainType::ImmersedForcing  ||
         solverChoice.buildings_type == BuildingsType::ImmersedForcing )
    {
        bool write_eb_surface = false;
        pp.query("write_eb_surface", write_eb_surface);
        if (write_eb_surface) {
            if (verbose > 0) {
                amrex::Print() << "Writing the geometry to a vtp file.\n" << std::endl;
            }
            WriteEBSurface(grids[finest_level],dmap[finest_level],Geom(finest_level),&EBFactory(finest_level));
        }
    }

}

void
ERF::Interp2DArrays (int lev, const BoxArray& my_ba2d, const DistributionMapping& my_dm)
{
    if (lon_m[lev-1] && !lon_m[lev]) {
        auto ngv = lon_m[lev-1]->nGrowVect(); ngv[2] = 0;
        lon_m[lev] = std::make_unique<MultiFab>(my_ba2d,my_dm,1,ngv);
        InterpFromCoarseLevel(*lon_m[lev], ngv, IntVect(0,0,0), // do not fill ghost cells outside the domain
                              *lon_m[lev-1], 0, 0, 1,
                              geom[lev-1], geom[lev],
                              refRatio(lev-1), &cell_cons_interp,
                              domain_bcs_type, BCVars::cons_bc);
    }
    if (lat_m[lev-1] && !lat_m[lev]) {
        auto ngv = lat_m[lev-1]->nGrowVect(); ngv[2] = 0;
        lat_m[lev] = std::make_unique<MultiFab>(my_ba2d,my_dm,1,ngv);
        InterpFromCoarseLevel(*lat_m[lev], ngv, IntVect(0,0,0), // do not fill ghost cells outside the domain
                              *lat_m[lev-1], 0, 0, 1,
                              geom[lev-1], geom[lev],
                              refRatio(lev-1), &cell_cons_interp,
                              domain_bcs_type, BCVars::cons_bc);
    }
    if (sinPhi_m[lev-1] && !sinPhi_m[lev]) {
        auto ngv = sinPhi_m[lev-1]->nGrowVect(); ngv[2] = 0;
        sinPhi_m[lev] = std::make_unique<MultiFab>(my_ba2d,my_dm,1,ngv);
        InterpFromCoarseLevel(*sinPhi_m[lev], ngv, IntVect(0,0,0), // do not fill ghost cells outside the domain
                              *sinPhi_m[lev-1], 0, 0, 1,
                              geom[lev-1], geom[lev],
                              refRatio(lev-1), &cell_cons_interp,
                              domain_bcs_type, BCVars::cons_bc);
    }
    if (cosPhi_m[lev-1] && !cosPhi_m[lev]) {
        auto ngv = cosPhi_m[lev-1]->nGrowVect(); ngv[2] = 0;
        cosPhi_m[lev] = std::make_unique<MultiFab>(my_ba2d,my_dm,1,ngv);
        InterpFromCoarseLevel(*cosPhi_m[lev], ngv, IntVect(0,0,0), // do not fill ghost cells outside the domain
                              *cosPhi_m[lev-1], 0, 0, 1,
                              geom[lev-1], geom[lev],
                              refRatio(lev-1), &cell_cons_interp,
                              domain_bcs_type, BCVars::cons_bc);
    }
    if (sst_lev[lev-1][0] && !sst_lev[lev][0]) {
        int ntimes = sst_lev[lev-1].size();
        sst_lev[lev].resize(ntimes);
        auto ngv = sst_lev[lev-1][0]->nGrowVect(); ngv[2] = 0;
        for (int n = 0; n < ntimes; n++) {
            sst_lev[lev][n] = std::make_unique<MultiFab>(my_ba2d,my_dm,1,ngv);
            InterpFromCoarseLevel(*sst_lev[lev][n], ngv, IntVect(0,0,0), // do not fill ghost cells outside the domain
                                  *sst_lev[lev-1][n], 0, 0, 1,
                                  geom[lev-1], geom[lev],
                                  refRatio(lev-1), &cell_cons_interp,
                                  domain_bcs_type, BCVars::cons_bc);
        }
    }
    if (tsk_lev[lev-1][0] && !tsk_lev[lev][0]) {
        int ntimes = tsk_lev[lev-1].size();
        tsk_lev[lev].resize(ntimes);
        auto ngv = tsk_lev[lev-1][0]->nGrowVect(); ngv[2] = 0;
        for (int n = 0; n < ntimes; n++) {
            tsk_lev[lev][n] = std::make_unique<MultiFab>(my_ba2d,my_dm,1,ngv);
            InterpFromCoarseLevel(*tsk_lev[lev][n], ngv, IntVect(0,0,0), // do not fill ghost cells outside the domain
                                  *tsk_lev[lev-1][n], 0, 0, 1,
                                  geom[lev-1], geom[lev],
                                  refRatio(lev-1), &cell_cons_interp,
                                  domain_bcs_type, BCVars::cons_bc);
        }
    }

    Real time_for_fp = 0.; // This is not actually used
    Vector<Real> ftime    = {time_for_fp, time_for_fp};
    Vector<Real> ctime    = {time_for_fp, time_for_fp};
    if (lat_m[lev]) {
        // Call FillPatchTwoLevels which ASSUMES that all ghost cells at lev-1 have already been filled
        Vector<MultiFab*> fmf = {lat_m[lev  ].get(), lat_m[lev  ].get()};
        Vector<MultiFab*> cmf = {lat_m[lev-1].get(), lat_m[lev-1].get()};
        IntVect ngv = lat_m[lev]->nGrowVect(); ngv[2] = 0;
        Interpolater* mapper = &cell_cons_interp;
        FillPatchTwoLevels(*lat_m[lev].get(), ngv, IntVect(0,0,0),
                           time_for_fp, cmf, ctime, fmf, ftime,
                           0, 0, 1, geom[lev-1], geom[lev],
                           refRatio(lev-1), mapper, domain_bcs_type,
                           BCVars::cons_bc);
    }
    if (lon_m[lev]) {
        // Call FillPatchTwoLevels which ASSUMES that all ghost cells at lev-1 have already been filled
        Vector<MultiFab*> fmf = {lon_m[lev  ].get(), lon_m[lev  ].get()};
        Vector<MultiFab*> cmf = {lon_m[lev-1].get(), lon_m[lev-1].get()};
        IntVect ngv = lon_m[lev]->nGrowVect(); ngv[2] = 0;
        Interpolater* mapper = &cell_cons_interp;
        FillPatchTwoLevels(*lon_m[lev].get(), ngv, IntVect(0,0,0),
                           time_for_fp, cmf, ctime, fmf, ftime,
                           0, 0, 1, geom[lev-1], geom[lev],
                           refRatio(lev-1), mapper, domain_bcs_type,
                           BCVars::cons_bc);
    } // lon_m
    if (sinPhi_m[lev]) {
        // Call FillPatchTwoLevels which ASSUMES that all ghost cells at lev-1 have already been filled
        Vector<MultiFab*> fmf = {sinPhi_m[lev  ].get(), sinPhi_m[lev  ].get()};
        Vector<MultiFab*> cmf = {sinPhi_m[lev-1].get(), sinPhi_m[lev-1].get()};
        IntVect ngv = sinPhi_m[lev]->nGrowVect(); ngv[2] = 0;
        Interpolater* mapper = &cell_cons_interp;
        FillPatchTwoLevels(*sinPhi_m[lev].get(), ngv, IntVect(0,0,0),
                           time_for_fp, cmf, ctime, fmf, ftime,
                           0, 0, 1, geom[lev-1], geom[lev],
                           refRatio(lev-1), mapper, domain_bcs_type,
                           BCVars::cons_bc);
    } // sinPhi
    if (cosPhi_m[lev]) {
        // Call FillPatchTwoLevels which ASSUMES that all ghost cells at lev-1 have already been filled
        Vector<MultiFab*> fmf = {cosPhi_m[lev  ].get(), cosPhi_m[lev  ].get()};
        Vector<MultiFab*> cmf = {cosPhi_m[lev-1].get(), cosPhi_m[lev-1].get()};
        IntVect ngv = cosPhi_m[lev]->nGrowVect(); ngv[2] = 0;
        Interpolater* mapper = &cell_cons_interp;
        FillPatchTwoLevels(*cosPhi_m[lev].get(), ngv, IntVect(0,0,0),
                           time_for_fp, cmf, ctime, fmf, ftime,
                           0, 0, 1, geom[lev-1], geom[lev],
                           refRatio(lev-1), mapper, domain_bcs_type,
                           BCVars::cons_bc);
    } // cosPhi
    if (sst_lev[lev][0]) {
        // Call FillPatchTwoLevels which ASSUMES that all ghost cells at lev-1 have already been filled
    int ntimes = sst_lev[lev].size();
    for (int n = 0; n < ntimes; n++) {
            Vector<MultiFab*> fmf = {sst_lev[lev  ][n].get(), sst_lev[lev  ][n].get()};
            Vector<MultiFab*> cmf = {sst_lev[lev-1][n].get(), sst_lev[lev-1][n].get()};
            IntVect ngv = sst_lev[lev][n]->nGrowVect(); ngv[2] = 0;
            Interpolater* mapper = &cell_cons_interp;
            FillPatchTwoLevels(*sst_lev[lev][n].get(), ngv, IntVect(0,0,0),
                               time_for_fp, cmf, ctime, fmf, ftime,
                               0, 0, 1, geom[lev-1], geom[lev],
                               refRatio(lev-1), mapper, domain_bcs_type,
                               BCVars::cons_bc);
        } // ntimes
    } // sst_lev
    if (tsk_lev[lev][0]) {
        // Call FillPatchTwoLevels which ASSUMES that all ghost cells at lev-1 have already been filled
    int ntimes = tsk_lev[lev].size();
    for (int n = 0; n < ntimes; n++) {
            Vector<MultiFab*> fmf = {tsk_lev[lev  ][n].get(), tsk_lev[lev  ][n].get()};
            Vector<MultiFab*> cmf = {tsk_lev[lev-1][n].get(), tsk_lev[lev-1][n].get()};
            IntVect ngv = tsk_lev[lev][n]->nGrowVect(); ngv[2] = 0;
            Interpolater* mapper = &cell_cons_interp;
            FillPatchTwoLevels(*tsk_lev[lev][n].get(), ngv, IntVect(0,0,0),
                               time_for_fp, cmf, ctime, fmf, ftime,
                               0, 0, 1, geom[lev-1], geom[lev],
                               refRatio(lev-1), mapper, domain_bcs_type,
                               BCVars::cons_bc);
        } // ntimes
    } // tsk_lev
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
    auto dRestartTime0 = amrex::second();

    ReadCheckpointFile();

    if (regrid_level_0_on_restart) {
        //
        // Coarsening before we split the grids ensures that each resulting
        // grid will have an even number of cells in each direction.
        //
        BoxArray new_ba(amrex::coarsen(Geom(0).Domain(),2));
        //
        // Now split up into list of grids within max_grid_size[0] limit.
        //
        new_ba.maxSize(max_grid_size[0]/2);
        //
        // Now refine these boxes back to level 0.
        //
        new_ba.refine(2);

        if (refine_grid_layout) {
            ChopGrids(0, new_ba, ParallelDescriptor::NProcs());
        }

        if (new_ba != grids[0]) {
            DistributionMapping new_dm(new_ba);
            RemakeLevel(0,t_new[0],new_ba,new_dm);
        }
    }

#ifdef ERF_USE_PARTICLES
    // We call this here without knowing whether the particles have already been initialized or not
    initializeTracers((ParGDBBase*)GetParGDB(),z_phys_nd,t_new[0]);
#endif

    Real cur_time = t_new[0];
    if (m_check_per    > 0.) {last_check_file_time    = cur_time;}
    if (m_plot2d_per_1 > 0.) {last_plot2d_file_time_1 = std::floor(cur_time/m_plot2d_per_1) * m_plot2d_per_1;}
    if (m_plot2d_per_2 > 0.) {last_plot2d_file_time_2 = std::floor(cur_time/m_plot2d_per_2) * m_plot2d_per_2;}
    if (m_plot3d_per_1 > 0.) {last_plot3d_file_time_1 = std::floor(cur_time/m_plot3d_per_1) * m_plot3d_per_1;}
    if (m_plot3d_per_2 > 0.) {last_plot3d_file_time_2 = std::floor(cur_time/m_plot3d_per_2) * m_plot3d_per_2;}

    if (m_check_int    > 0.) {last_check_file_step    = istep[0];}
    if (m_plot2d_int_1 > 0.) {last_plot2d_file_step_1 = istep[0];}
    if (m_plot2d_int_2 > 0.) {last_plot2d_file_step_2 = istep[0];}
    if (m_plot3d_int_1 > 0.) {last_plot3d_file_step_1 = istep[0];}
    if (m_plot3d_int_2 > 0.) {last_plot3d_file_step_2 = istep[0];}

    if (verbose > 0)
    {
        auto dRestartTime = amrex::second() - dRestartTime0;
        ParallelDescriptor::ReduceRealMax(dRestartTime,ParallelDescriptor::IOProcessorNumber());
        amrex::Print() << "Restart time = " << dRestartTime << " seconds." << '\n';
    }
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
    if (solverChoice.init_type == InitType::Input_Sounding) {
        // The physbc's need the terrain but are needed for initHSE
        // We have already made the terrain in the call to init_zphys
        //    in MakeNewLevelFromScratch
        make_physbcs(lev);

        // Now init the base state and the data itself
        init_from_input_sounding(lev);

        // The base state has been initialized by integrating vertically
        // through the sounding for ideal (like WRF) or isentropic approaches
        if (solverChoice.sounding_type == SoundingType::Ideal ||
            solverChoice.sounding_type == SoundingType::Isentropic ||
            solverChoice.sounding_type == SoundingType::DryIsentropic) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(solverChoice.use_gravity,
                "Gravity should be on to be consistent with sounding initialization.");
        } else { // SoundingType::ConstantDensity
            AMREX_ASSERT_WITH_MESSAGE(!solverChoice.use_gravity,
                "Constant density probably doesn't make sense with gravity");
            initHSE();
        }

#ifdef ERF_USE_NETCDF
    }
    else if (solverChoice.init_type == InitType::WRFInput && !nc_init_file[lev].empty())
    {
        // The base state is initialized from WRF wrfinput data, output by
        // ideal.exe or real.exe

        init_from_wrfinput(lev, *mf_C1H, *mf_C2H, *mf_MUB, *mf_PSFC[lev]);

        if (lev==0) {
            if ((start_time > 0) && (start_time != start_bdy_time)) {
                Print() << "Ignoring specified start_time="
                        << std::setprecision(timeprecision) << start_time
                        << std::endl;
            }
        }

        start_time = start_bdy_time;

        use_datetime = true;

        // The physbc's need the terrain but are needed for initHSE
        if (!solverChoice.use_real_bcs) {
            make_physbcs(lev);
        }
    }
    else if (solverChoice.init_type == InitType::WRFInput && nc_init_file[lev].empty())
    {
        amrex::Abort("This pathway is not quite implemented yet");
    }
    else if (solverChoice.init_type == InitType::NCFile)
    {
        // The state is initialized by reading from a Netcdf file
        init_from_ncfile(lev);

        // The physbc's need the terrain but are needed for initHSE
        make_physbcs(lev);
    }
    else if (solverChoice.init_type == InitType::Metgrid)
    {
        // The base state is initialized from data output by WPS metgrid;
        // we will rebalance after interpolation
        init_from_metgrid(lev);
#endif
    } else if (solverChoice.init_type == InitType::Uniform) {
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
    if (solverChoice.init_type != InitType::NCFile) {
        init_custom(lev);
    }

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
        solverChoice.pert_type == PerturbationType::Direct ||
        solverChoice.pert_type == PerturbationType::CPM) {
        turbPert_update(lev, 0.);
        turbPert_amplitude(lev);
    }

    // Set initial velocity field for immersed cells to be close to 0
    if (solverChoice.terrain_type == TerrainType::ImmersedForcing ||
        solverChoice.buildings_type == BuildingsType::ImmersedForcing) {
        init_immersed_forcing(lev);
    }
}

// Read in some parameters from inputs file
void
ERF::ReadParameters ()
{
    {
        ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
        pp.query("max_step", max_step);
        if (max_step < 0) {
            max_step = std::numeric_limits<int>::max();
        }

        // TODO: more robust general datetime parsing
        std::string start_datetime, stop_datetime;
        if (pp.query("start_datetime", start_datetime)) {
            if (start_datetime.length() == 16) { // YYYY-MM-DD HH:MM
                start_datetime += ":00"; // add seconds
            }
            if (start_datetime.length() != 19) {
                Print() << "Got start_datetime = \"" << start_datetime
                    << "\", format should be " << datetime_format << std::endl;
                exit(0);
            }
            start_time = getEpochTime(start_datetime, datetime_format);
            Print() << "Start datetime : " << start_datetime << std::endl;

            if (pp.query("stop_datetime", stop_datetime)) {
                if (stop_datetime.length() == 16) { // YYYY-MM-DD HH:MM
                    stop_datetime += ":00"; // add seconds
                }
                if (stop_datetime.length() != 19) {
                    Print() << "Got stop_datetime = \"" << stop_datetime
                        << "\", format should be " << datetime_format << std::endl;
                    exit(0);
                }
                stop_time = getEpochTime(stop_datetime, datetime_format);
                Print() << "Stop  datetime : " << start_datetime << std::endl;
            } else if (pp.query("stop_time", stop_time)) {
                Print() << "Sim length     : " << stop_time << " s" << std::endl;
                stop_time += start_time;
            }

            use_datetime = true;

        } else {
            pp.query("stop_time", stop_time);
            pp.query("start_time", start_time); // This is optional, it defaults to 0
        }
    }

    ParmParse pp(pp_prefix);
    ParmParse pp_amr("amr");
    {
        pp.query("regrid_level_0_on_restart", regrid_level_0_on_restart);
        pp.query("regrid_int", regrid_int);
        pp.query("check_file", check_file);

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
#ifndef ERF_USE_FFT
        if (use_fft) {
            Abort("You must build with USE_FFT in order to set use_fft = true in your inputs file");
        }
#endif

        // Check for NaNs?
        pp.query("check_for_nans", check_for_nans);

        // Frequency of diagnostic output
        pp.query("sum_interval", sum_interval);
        pp.query("sum_period"  , sum_per);

        pp.query("pert_interval", pert_interval);

        // Time step controls
        pp.query("cfl", cfl);
        pp.query("substepping_cfl", sub_cfl);
        pp.query("init_shrink", init_shrink);
        pp.query("change_max", change_max);
        pp.query("dt_max_initial", dt_max_initial);
        pp.query("dt_max", dt_max);

        fixed_dt.resize(max_level+1,-1.);
        fixed_fast_dt.resize(max_level+1,-1.);

        pp.query("fixed_dt", fixed_dt[0]);
        pp.query("fixed_fast_dt", fixed_fast_dt[0]);

        int nlevs_max = max_level + 1;
        istep.resize(nlevs_max, 0);
        nsubsteps.resize(nlevs_max, 1);
        // This is the default
        for (int lev = 1; lev <= max_level; ++lev) {
            nsubsteps[lev] = MaxRefRatio(lev-1);
        }

        if (max_level > 0) {
            ParmParse pp_erf("erf");
            int count = pp_erf.countval("dt_ref_ratio");
            if (count > 0) {
                Vector<int> nsub;
                nsub.resize(nlevs_max, 0);
                if (count == 1) {
                    pp_erf.queryarr("dt_ref_ratio", nsub, 0, 1);
                    for (int lev = 1; lev <= max_level; ++lev) {
                        nsubsteps[lev] = nsub[0];
                    }
                } else {
                    pp_erf.queryarr("dt_ref_ratio", nsub, 0, max_level);
                    for (int lev = 1; lev <= max_level; ++lev) {
                        nsubsteps[lev] = nsub[lev-1];
                    }
                }
            }
        }

        // Make sure we do this after we have defined nsubsteps above
        for (int lev = 1; lev <= max_level; lev++)
        {
            fixed_dt[lev]      = fixed_dt[lev-1]      / static_cast<Real>(nsubsteps[lev]);
            fixed_fast_dt[lev] = fixed_fast_dt[lev-1] / static_cast<Real>(nsubsteps[lev]);
        }

        pp.query("fixed_mri_dt_ratio", fixed_mri_dt_ratio);

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
        have_read_nc_init_file.resize(max_level+1);

        // NetCDF wrfinput initialization files -- possibly multiple files at each of multiple levels
        //        but we always have exactly one file at level 0
        for (int lev = 0; lev <= max_level; lev++) {
            const std::string nc_file_names = Concatenate("nc_init_file_",lev,1);
            if (pp.contains(nc_file_names.c_str())) {
                int num_files = pp.countval(nc_file_names.c_str());
                num_files_at_level[lev] = num_files;
                nc_init_file[lev].resize(num_files);
                have_read_nc_init_file[lev].resize(num_files);
                pp.queryarr(nc_file_names.c_str(), nc_init_file[lev],0,num_files);
                for (int j = 0; j < num_files; j++) {
                    Print() << "Reading NC init file names at level " << lev << " and index " << j << " : " << nc_init_file[lev][j] << std::endl;
                    have_read_nc_init_file[lev][j] = 0;
                } // j
            } // if pp.contains
        } // lev

        // NetCDF wrfbdy lateral boundary file
        if (pp.query("nc_bdy_file", nc_bdy_file)) {
            Print() << "Reading NC bdy file name " << nc_bdy_file << std::endl;
        }

        // NetCDF wrflow lateral boundary file
        if (pp.query("nc_low_file", nc_low_file)) {
            Print() << "Reading NC low file name " << nc_low_file << std::endl;
        }

#endif

        // Options for vertical interpolation of met_em*.nc data.
        pp.query("metgrid_debug_quiescent",  metgrid_debug_quiescent);
        pp.query("metgrid_debug_isothermal", metgrid_debug_isothermal);
        pp.query("metgrid_debug_dry",        metgrid_debug_dry);
        pp.query("metgrid_debug_psfc",       metgrid_debug_psfc);
        pp.query("metgrid_debug_msf",        metgrid_debug_msf);
        pp.query("metgrid_interp_theta",     metgrid_interp_theta);
        pp.query("metgrid_basic_linear",     metgrid_basic_linear);
        pp.query("metgrid_use_below_sfc",    metgrid_use_below_sfc);
        pp.query("metgrid_use_sfc",          metgrid_use_sfc);
        pp.query("metgrid_retain_sfc",       metgrid_retain_sfc);
        pp.query("metgrid_proximity",        metgrid_proximity);
        pp.query("metgrid_order",            metgrid_order);
        pp.query("metgrid_force_sfc_k",      metgrid_force_sfc_k);

        // Set default to FullState for now ... later we will try Perturbation
        interpolation_type = StateInterpType::FullState;
        pp.query_enum_case_insensitive("interpolation_type"  ,interpolation_type);

        PlotFileType plotfile3d_type_temp = PlotFileType::None;
        pp.query_enum_case_insensitive("plotfile_type"  ,plotfile3d_type_temp);
        pp.query_enum_case_insensitive("plotfile_type_1",plotfile3d_type_1);
        pp.query_enum_case_insensitive("plotfile_type_2",plotfile3d_type_2);

        PlotFileType plotfile2d_type_temp = PlotFileType::None;
        pp.query_enum_case_insensitive("plotfile2d_type"  ,plotfile2d_type_temp);
        pp.query_enum_case_insensitive("plotfile2d_type_1",plotfile2d_type_1);
        pp.query_enum_case_insensitive("plotfile2d_type_2",plotfile2d_type_2);
        //
        // This option is for backward consistency -- if only plotfile_type is set,
        //     then it will be used for both 1 and 2 if and only if they are not set
        //
        // Default is native amrex if no type is specified
        //
        if (plotfile3d_type_temp == PlotFileType::None) {
            if (plotfile3d_type_1 == PlotFileType::None) {
                plotfile3d_type_1  = PlotFileType::Amrex;
            }
            if (plotfile3d_type_2 == PlotFileType::None) {
                plotfile3d_type_2  = PlotFileType::Amrex;
            }
        } else {
            if (plotfile3d_type_1 == PlotFileType::None) {
                plotfile3d_type_1  = plotfile3d_type_temp;
            } else {
                Abort("You must set either plotfile_type or plotfile_type_1, not both");
            }
            if (plotfile3d_type_2 == PlotFileType::None) {
                plotfile3d_type_2  = plotfile3d_type_temp;
            } else {
                Abort("You must set either plotfile_type or plotfile_type_2, not both");
            }
        }
        if (plotfile2d_type_temp == PlotFileType::None) {
            if (plotfile2d_type_1 == PlotFileType::None) {
                plotfile2d_type_1  = PlotFileType::Amrex;
            }
            if (plotfile2d_type_2 == PlotFileType::None) {
                plotfile2d_type_2  = PlotFileType::Amrex;
            }
        } else {
            if (plotfile2d_type_1 == PlotFileType::None) {
                plotfile2d_type_1  = plotfile2d_type_temp;
            } else {
                Abort("You must set either plotfile2d_type or plotfile2d_type_1, not both");
            }
            if (plotfile2d_type_2 == PlotFileType::None) {
                plotfile2d_type_2  = plotfile2d_type_temp;
            } else {
                Abort("You must set either plotfile2d_type or plotfile2d_type_2, not both");
            }
        }
#ifndef ERF_USE_NETCDF
        if (plotfile3d_type_1 == PlotFileType::Netcdf ||
            plotfile3d_type_2 == PlotFileType::Netcdf ||
            plotfile2d_type_1 == PlotFileType::Netcdf ||
            plotfile2d_type_2 == PlotFileType::Netcdf) {
            Abort("Plotfile type = Netcdf is not allowed without USE_NETCDF = TRUE");
        }
#endif

        pp.query("plot_file_1"  ,   plot3d_file_1);
        pp.query("plot_file_2"  ,   plot3d_file_2);
        pp.query("plot2d_file_1",   plot2d_file_1);
        pp.query("plot2d_file_2",   plot2d_file_2);

        pp.query("plot_int_1"   , m_plot3d_int_1);
        pp.query("plot_int_2"   , m_plot3d_int_2);
        pp.query("plot_per_1"   , m_plot3d_per_1);
        pp.query("plot_per_2"   , m_plot3d_per_2);

        pp.query("plot2d_int_1" , m_plot2d_int_1);
        pp.query("plot2d_int_2" , m_plot2d_int_2);
        pp.query("plot2d_per_1",  m_plot2d_per_1);
        pp.query("plot2d_per_2",  m_plot2d_per_2);

        pp.query("subvol_file",   subvol_file);

        // Should we use format like plt1970-01-01_00:00:00.000000 (if true) or plt00001 (if false)
        pp.query("use_real_time_in_pltname", use_real_time_in_pltname);

        // If use_real_time_in_pltname is false, how many digits should we use for the timestep?
        pp.query("file_name_digits", file_name_digits);

        // Default if subvol_int not specified
        m_subvol_int.resize(1); m_subvol_int[0] = -1;
        m_subvol_per.resize(1); m_subvol_per[0] = -1.0;
        last_subvol_step.resize(1);
        last_subvol_time.resize(1);

        int nsi = pp.countval("subvol_int");
        int nsr = pp.countval("subvol_per");

        // We must specify only subvol_int OR subvol_per
        AMREX_ALWAYS_ASSERT (!(nsi > 0 && nsr > 0));

        int nsub = -1;
        if (nsi > 0 || nsr > 0) {
            ParmParse pp_sv("erf.subvol");
            int n1 = pp_sv.countval("origin"); int n2 = pp_sv.countval("nxnynz"); int n3 = pp_sv.countval("dxdydz");
            if (n1 != n2 || n1 != n3 || n2 != n3) {
                Abort("WriteSubvolume: must have same number of entries in origin, nxnynz, and dxdydz.");
            }
            if ( n1%AMREX_SPACEDIM != 0) {
                Abort("WriteSubvolume: origin, nxnynz, and dxdydz must have multiples of AMReX_SPACEDIM");
            }
            nsub = n1/AMREX_SPACEDIM;
            m_subvol_int.resize(nsub);
            last_subvol_step.resize(nsub);
            last_subvol_time.resize(nsub);
            m_subvol_int.resize(nsub);
            m_subvol_per.resize(nsub);
        }

        if (nsi > 0) {
            for (int i = 1; i < nsub; i++) m_subvol_per[i] = -1.0;
            if ( nsi == 1) {
                m_subvol_int[0] = -1;
                pp.get("subvol_int" , m_subvol_int[0]);
            } else if ( nsi == nsub) {
                pp.getarr("subvol_int" , m_subvol_int);
            } else {
                Abort("There must either be a single value of subvol_int or one for every subdomain");
            }
        }

        if (nsr > 0) {
            for (int i = 1; i < nsub; i++) m_subvol_int[i] = -1.0;
            if ( nsr == 1) {
                m_subvol_per[0] = -1.0;
                pp.get("subvol_per" , m_subvol_per[0]);
            } else if ( nsr == nsub) {
                pp.getarr("subvol_per" , m_subvol_per);
            } else {
                Abort("There must either be a single value of subvol_per or one for every subdomain");
            }
        }

        setSubVolVariables("subvol_sampling_vars",subvol3d_var_names);

        pp.query("expand_plotvars_to_unif_rr",m_expand_plotvars_to_unif_rr);

        pp.query("plot_face_vels",m_plot_face_vels);

        if ( (m_plot3d_int_1 > 0 && m_plot3d_per_1 > 0) ||
             (m_plot3d_int_2 > 0 && m_plot3d_per_2 > 0.) ) {
            Abort("Must choose only one of plot_int or plot_per");
        }
        if ( (m_plot2d_int_1 > 0 && m_plot2d_per_1 > 0) ||
             (m_plot2d_int_2 > 0 && m_plot2d_per_2 > 0.) ) {
            Abort("Must choose only one of plot_int or plot_per");
        }

        pp.query("profile_int", profile_int);
        pp.query("destag_profiles", destag_profiles);

        pp.query("plot_lsm", plot_lsm);
#ifdef ERF_USE_RRTMGP
        pp.query("plot_rad", plot_rad);
#endif
        pp.query("profile_rad_int", rad_datalog_int);

        pp.query("output_1d_column", output_1d_column);
        pp.query("column_per", column_per);
        pp.query("column_interval", column_interval);
        pp.query("column_loc_x", column_loc_x);
        pp.query("column_loc_y", column_loc_y);
        pp.query("column_file_name", column_file_name);

        // Sampler output frequency
        pp.query("line_sampling_per", line_sampling_per);
        pp.query("line_sampling_interval", line_sampling_interval);
        pp.query("plane_sampling_per", plane_sampling_per);
        pp.query("plane_sampling_interval", plane_sampling_interval);

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

        // If using real boundaries, do we extrapolate w (or set to 0)
        pp.query("real_extrap_w", real_extrap_w);

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

    solverChoice.init_params(max_level,pp_prefix);

#ifndef ERF_USE_NETCDF
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(( (solverChoice.init_type != InitType::WRFInput) &&
                                       (solverChoice.init_type != InitType::Metgrid ) &&
                                       (solverChoice.init_type != InitType::NCFile  )  ),
                                     "init_type cannot be 'WRFInput', 'MetGrid' or 'NCFile' if we don't build with netcdf!");
#endif

    // Query the canopy model file name
    std::string forestfile;
    solverChoice.do_forest_drag = pp.query("forest_file", forestfile);
    if (solverChoice.do_forest_drag) {
        for (int lev = 0; lev <= max_level; ++lev) {
            m_forest_drag[lev] = std::make_unique<ForestDrag>(forestfile);
        }
    }

    // If init from WRFInput or Metgrid make sure a valid file name is present at level 0.
    // We allow for the possibility that finer levels may use native refinement rather than reading from a file
    if ((solverChoice.init_type == InitType::WRFInput) ||
        (solverChoice.init_type == InitType::Metgrid)  ||
        (solverChoice.init_type == InitType::NCFile) ) {
        int num_files = nc_init_file[0].size();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(num_files>0, "A file name must be present at level 0 for init type WRFInput, Metgrid or NCFile.");
        for (int j = 0; j < num_files; j++) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!nc_init_file[0][j].empty(), "Valid file name must be present at level 0 for init type WRFInput, Metgrid or NCFile.");
        } //j
    } // InitType

    // What type of land surface model to use
    // NOTE: Must be checked after init_params
    if (solverChoice.lsm_type == LandSurfaceType::SLM) {
        lsm.SetModel<SLM>();
        Print() << "SLM land surface model!\n";
    } else if (solverChoice.lsm_type == LandSurfaceType::MM5) {
        lsm.SetModel<MM5>();
        Print() << "MM5 land surface model!\n";
#ifdef ERF_USE_NOAHMP
    } else if (solverChoice.lsm_type == LandSurfaceType::NOAHMP) {
        lsm.SetModel<NOAHMP>();
        Print() << "Noah-MP land surface model!\n";
#endif
    } else if (solverChoice.lsm_type == LandSurfaceType::None) {
        lsm.SetModel<NullSurf>();
        Print() << "Null land surface model!\n";
    } else {
        Abort("Dont know this LandSurfaceType!") ;
    }

    if (verbose > 0) {
        solverChoice.display(max_level,pp_prefix);
    }

    ParameterSanityChecks();
}

// Read in some parameters from inputs file
void
ERF::ParameterSanityChecks ()
{
    AMREX_ALWAYS_ASSERT(cfl > 0. || fixed_dt[0] > 0.);

    // We don't allow use_real_bcs to be true if init_type is not either InitType::WRFInput or InitType::Metgrid
    AMREX_ALWAYS_ASSERT( !solverChoice.use_real_bcs ||
                        ((solverChoice.init_type == InitType::WRFInput) || (solverChoice.init_type == InitType::Metgrid)) );

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
                    fab_arr(i,j,k,2) = hse_arr(i,j,k,BaseState::p0_comp);
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
                    fab_arr(i,j,k,2) = hse_arr(i,j,k,BaseState::p0_comp);
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
ERF::writeNow(const Real cur_time, const int nstep, const int plot_int, const Real plot_per,
              const Real dt_0, Real& next_file_time)
{
    bool write_now = false;

    if ( plot_int > 0) {

        write_now = (nstep % plot_int == 0);

    } else if (plot_per > 0.0) {

        amrex::Print() << "CUR NEXT PER " << cur_time << " " << next_file_time << " " << plot_per << std::endl;

        // Only write now if nstep newly matches the number of elapsed periods
        write_now = (cur_time > (next_file_time - Real(0.1)*dt_0));
    }

    return write_now;
}

void
ERF::check_state_for_nans(MultiFab const& S)
{
    int ncomp = S.nComp();
    for (int lev = 0; lev <= finest_level; lev++)
    {
        //
        // Test at the end of every full timestep whether the solution data contains NaNs
        //
        bool any_have_nans = false;
        for (int i = 0; i < ncomp; i++) {
            if (S.contains_nan(i,1,0))
            {
                amrex::Print() << "Component " << i << " of conserved variables contains NaNs" << '\n';
                any_have_nans = true;
            }
        }
        if (any_have_nans) {
            exit(0);
        }
    }
}

void
ERF::check_vels_for_nans(MultiFab const& xvel, MultiFab const& yvel, MultiFab const& zvel)
{
    //
    // Test at the end of every full timestep whether the solution data contains NaNs
    //
    bool any_have_nans = false;
    if (xvel.contains_nan(0,1,0))
    {
        amrex::Print() << "x-velocity contains NaNs " << '\n';
        any_have_nans = true;
    }
    if (yvel.contains_nan(0,1,0))
    {
        amrex::Print() << "y-velocity contains NaNs" << '\n';
        any_have_nans = true;
    }
    if (zvel.contains_nan(0,1,0))
    {
        amrex::Print() << "z-velocity contains NaNs" << '\n';
        any_have_nans = true;
    }
    if (any_have_nans) {
        exit(0);
    }
}

void
ERF::check_for_low_temp(amrex::MultiFab& S)
{
    // *****************************************************************************
    // Test for low temp (low is defined as beyond the microphysics range of validity)
    // *****************************************************************************
    //
    // This value is defined in erf_dtesati in Source/Utils/ERF_MicrophysicsUtils.H
    Real t_low = 273.16 - 85.;
    //
    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox();
        const Array4<Real> &s_arr  = S.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real rho      = s_arr(i, j, k, Rho_comp);
            const Real rhotheta = s_arr(i, j, k, RhoTheta_comp);
            const Real qv       = s_arr(i, j, k, RhoQ1_comp);

            Real temp = getTgivenRandRTh(rho, rhotheta, qv);

            if (temp < t_low) {
#ifdef AMREX_USE_GPU
                AMREX_DEVICE_PRINTF("Temperature too low in cell: %d %d %d %e \n", i,j,k,temp);
#else
                printf("Temperature too low in cell: %d %d %d \n", i,j,k);
                printf("Based on temp / rhotheta / rho %e %e %e \n", temp,rhotheta,rho);
                Abort();
#endif
            }
        });
    }
}

void
ERF::check_for_negative_theta(amrex::MultiFab& S)
{
    // *****************************************************************************
    // Test for negative (rho theta)
    // *****************************************************************************
    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox();
        const Array4<Real> &s_arr  = S.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real rhotheta = s_arr(i, j, k, RhoTheta_comp);
            if (rhotheta <= 0.) {
#ifdef AMREX_USE_GPU
                AMREX_DEVICE_PRINTF("RhoTheta is negative at %d %d %d %e \n", i,j,k,rhotheta);
#else
                printf("RhoTheta is negative at %d %d %d %e \n", i,j,k,rhotheta);
                Abort("Bad theta in check_for_negative_theta");
#endif
            }
            });
    } // mfi
}
