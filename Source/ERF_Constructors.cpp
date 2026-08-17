/**
 * \file ERF_Constructors.cpp
 */

/**
 * Main class in ERF code, instantiated from main.cpp
*/

#include "ERF.H"
#include "AMReX_buildInfo.H"
#include "AMReX_EB2_IF_Box.H"
#include "AMReX_EB2_IF_Sphere.H"
#include "AMReX_EB2_IF_Plane.H"
#include "ERF_EBIFTerrain.H"

using namespace amrex;

// Constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
ERF::ERF ()
{
    ERF_shared();
}

#ifdef ERF_USE_MULTIBLOCK
// Constructor used when ERF is created by a multiblock driver
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

void
ERF::ERF_shared ()
{
    // Seeding lives here so every constructor gets it
    int fix_random_seed = 0;
    long random_seed = -1;
    {
        ParmParse pp("erf");
        pp.query("fix_random_seed", fix_random_seed);
        pp.query("random_seed", random_seed);
    }
    // Note that the value of 1024UL is not significant -- the point here is just to set the
    // same seed for all MPI processes for the purpose of regression testing
    if (fix_random_seed) {
        Print() << "Fixing the random seed" << std::endl;
        InitRandom(1024UL, ParallelDescriptor::NProcs(), 1024UL);
    } else if (random_seed >= 0) {
        // User-supplied seed: vary the random sampling per run (e.g. across ensemble
        // realizations), still offset by rank so the ranks draw independent streams.
        Print() << "Using user random seed " << random_seed << std::endl;
        auto s = static_cast<unsigned long>(random_seed);
        InitRandom(s + static_cast<unsigned long>(ParallelDescriptor::MyProc()) + 1UL,
                   ParallelDescriptor::NProcs(), s * 1234567UL + 12345UL);
    }

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

    // NOTE: size lsm before readparams (chooses the model at all levels)
    lsm.ReSize(nlevs_max);
    lsm_data.resize(nlevs_max);
    lsm_flux.resize(nlevs_max);

    nudge_data.resize(nlevs_max);
    lsf_data.resize(nlevs_max);

    rhotheta_src.resize(nlevs_max);
    rhoqt_src.resize(nlevs_max);

    // NOTE: size canopy model before readparams (if file exists, we construct)
    m_forest_drag.resize(nlevs_max);
    for (int lev = 0; lev <= max_level; ++lev) { m_forest_drag[lev] = nullptr;}

    ReadParameters();
    // Create one invocation identity after inputs are available and before
    // InitData can read restart metadata or write an output on restart.
    execution_provenance = erf_provenance::initialize_execution_provenance();
    initializeMicrophysics(nlevs_max);

#ifdef ERF_USE_WINDFARM
    initializeWindFarm(nlevs_max);
#endif

#ifdef ERF_USE_EAMXX_SHOC
    eamxx_shoc_interface.resize(nlevs_max);
    for (int lev = 0; lev <= max_level; ++lev) {
        if (solverChoice.turbChoice[lev].uses_eamxx_shoc()) {
            eamxx_shoc_interface[lev] = std::make_unique<SHOCInterface>(lev, solverChoice);
        }
    }
#endif

    native_shoc_driver.resize(nlevs_max);
    for (int lev = 0; lev <= max_level; ++lev) {
        if (solverChoice.turbChoice[lev].uses_native_shoc()) {
            native_shoc_driver[lev] = std::make_unique<ShocDriver>(lev, solverChoice);
        }
    }

    rad.resize(nlevs_max);
    for (int lev = 0; lev <= max_level; ++lev) {
        if (solverChoice.rad_type == RadiationType::RRTMGP) {
#ifdef ERF_USE_RRTMGP
            rad[lev] = std::make_unique<Radiation>(lev, solverChoice);
            // pass radiation datalog frequency to model - RRTMGP needs to know when to save data for profiles
            rad[lev]->setDataLogFrequency(rad_datalog_int);
#endif
        } else if (solverChoice.rad_type == RadiationType::Simple) {
            rad[lev] = std::make_unique<RadiationSimple>(lev, solverChoice);
            rad[lev]->setDataLogFrequency(rad_datalog_int);
        } else if (solverChoice.rad_type != RadiationType::None) {
            Abort("Don't know this radiation model!");
        }
    }
    const std::string& pv3d_1 = "plot_vars_1"  ; setPlotVariables(pv3d_1,plot3d_var_names_1);
    const std::string& pv3d_2 = "plot_vars_2"  ; setPlotVariables(pv3d_2,plot3d_var_names_2);

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
        if (std::fabs(zlevels_stag[0][nz-1]-geom[0].ProbHi(2)) > Real(1.0e-4)) {
            Print() << "Note: prob_hi[2]=" << geom[0].ProbHi(2)
                << " does not match highest requested z level " << zlevels_stag[0][nz-1]
                << std::endl;
        }
        if (std::fabs(zlevels_stag[0][0]-geom[0].ProbLo(2)) > Real(1.0e-4)) {
            Print() << "Note: prob_lo[2]=" << geom[0].ProbLo(2)
                << " does not match lowest requested level " << zlevels_stag[0][0]
                << std::endl;
        }

        // Redefine the problem domain here?
    }

    // Get lo/hi indices for massflux calc
    if ((solverChoice.const_massflux_u != 0) || (solverChoice.const_massflux_v != 0)) {
        if (solverChoice.mesh_type == MeshType::ConstantDz) {
            const bool zlo_unset = (solverChoice.const_massflux_layer_lo == amrex::Real(-bogus_large_value));
            const bool zhi_unset = (solverChoice.const_massflux_layer_hi == amrex::Real( bogus_large_value));
            const Real massflux_zlo = solverChoice.const_massflux_layer_lo - geom[0].ProbLo(2);
            const Real massflux_zhi = solverChoice.const_massflux_layer_hi - geom[0].ProbLo(2);
            const Real dz = geom[0].CellSize(2);
            if (zlo_unset) {
                solverChoice.massflux_klo = geom[0].Domain().smallEnd(2);
            } else {
                solverChoice.massflux_klo = static_cast<int>(std::ceil(massflux_zlo / dz - myhalf));
            }
            if (zhi_unset) {
                solverChoice.massflux_khi = geom[0].Domain().bigEnd(2);
            } else {
                solverChoice.massflux_khi = static_cast<int>(std::floor(massflux_zhi / dz - myhalf));
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

#ifdef ERF_REMORA_FORCE_PROBINIT_LINK
    extern void erf_probinit_link_anchor_func () noexcept;
    erf_probinit_link_anchor_func();
#endif
    prob = amrex_probinit(geom[0].ProbLo(),geom[0].ProbHi());

    ParmParse pp_erf("erf");
    std::string prob_name;
    pp_erf.query("prob_name", prob_name);
    const std::string prob_name_ci = amrex::toLower(prob_name);
    if (prob_name_ci == "cloud chamber" || prob_name_ci == "cloudchamber") {
        cloud_chamber_config = erf_cloud_chamber::parse_config(
            geom[0].ProbLo(), geom[0].ProbHi());
    }
    {
        int budget_interval = 0;
        pp_erf.query("cloud_chamber_budget_interval", budget_interval);
        if (budget_interval > 0) {
            if (!cloud_chamber_config.active ||
                !cloud_chamber_config.physical_initialization ||
                (solverChoice.moisture_type != MoistureType::None &&
                 solverChoice.moisture_type != MoistureType::SatAdj)) {
                Error("Cloud Chamber: cloud_chamber_budget_interval requires physical_temperature_rh with no moisture or SatAdj");
            }
            cloud_chamber_budget = std::make_unique<CloudChamberBudget>(budget_interval);
        }
    }

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    t_new.resize(nlevs_max, zero);
    t_old.resize(nlevs_max, -bogus_large_value);
    dt.resize(nlevs_max, std::min(static_cast<double>(bogus_large_value),dt_max_initial));
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
    SFS_hfx1_lev.resize(nlevs_max);  SFS_hfx2_lev.resize(nlevs_max);  SFS_hfx3_lev.resize(nlevs_max);
    SFS_diss_lev.resize(nlevs_max);
    SFS_q1fx1_lev.resize(nlevs_max); SFS_q1fx2_lev.resize(nlevs_max); SFS_q1fx3_lev.resize(nlevs_max);
    SFS_q2fx3_lev.resize(nlevs_max);
    eddyDiffs_lev.resize(nlevs_max);
    SmnSmn_lev.resize(nlevs_max);
    Tau_EB.resize(nlevs_max);
    hfx3_EB.resize(nlevs_max);

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
    terrain_blanking_xface.resize(nlevs_max);
    terrain_blanking_yface.resize(nlevs_max);
    terrain_blanking_zface.resize(nlevs_max);

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

    // Rayleigh damping
    h_rayleigh_ptrs.resize(nlevs_max);
    d_rayleigh_ptrs.resize(nlevs_max);
    h_sinesq_ptrs.resize(nlevs_max);
    d_sinesq_ptrs.resize(nlevs_max);
    h_sinesq_stag_ptrs.resize(nlevs_max);
    d_sinesq_stag_ptrs.resize(nlevs_max);

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
    std::string geometry ="terrain";
    ParmParse pp_eb2("eb2");
    pp_eb2.queryAdd("geometry", geometry);
    if ( solverChoice.terrain_type == TerrainType::EB ||
         solverChoice.terrain_type == TerrainType::ImmersedForcing)
    {
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
            double dummy_time = 0.0;
            prob->init_terrain_surface(geom[max_level], terrain_fab, dummy_time);
            TerrainIF implicit_fun(terrain_fab, geom[max_level], stretched_dz_d[max_level]);
            auto gshop = EB2::makeShop(implicit_fun);
            if (build_eb_for_multigrid) {
                EB2::Build(gshop, geom[max_level], max_level, max_coarsening_level,
                            ngrow_for_eb, build_coarse_level_by_coarsening);
            } else {
                EB2::Build(gshop, this->Geom(), ngrow_for_eb);
#if USE_FC_FACTORY
                EB2::BuildFC();
#endif
            }
        } else if (geometry == "plane") {
            RealArray plane_point{zero, zero, zero};
            RealArray plane_normal{zero, zero, -one}; // pointing into the solid region
            pp_eb2.query("plane_point", plane_point);
            pp_eb2.query("plane_normal", plane_normal);
            EB2::PlaneIF implicit_fun(plane_point, plane_normal, true);
            auto gshop = EB2::makeShop(implicit_fun);
            if (build_eb_for_multigrid) {
                EB2::Build(gshop, geom[max_level], max_level, max_coarsening_level,
                            ngrow_for_eb, build_coarse_level_by_coarsening);
            } else {
                EB2::Build(gshop, this->Geom(), ngrow_for_eb);
#if USE_FC_FACTORY
                EB2::BuildFC();
#endif
            }
        } else if (geometry == "box") {
            RealArray box_lo{zero, zero, zero};
            RealArray box_hi{zero, zero, zero};
            pp_eb2.query("box_lo", box_lo);
            pp_eb2.query("box_hi", box_hi);
            EB2::BoxIF implicit_fun(box_lo, box_hi, false);
            auto gshop = EB2::makeShop(implicit_fun);
            if (build_eb_for_multigrid) {
                EB2::Build(gshop, geom[max_level], max_level, max_coarsening_level,
                            ngrow_for_eb, build_coarse_level_by_coarsening);
            } else {
                EB2::Build(gshop, this->Geom(), ngrow_for_eb);
#if USE_FC_FACTORY
                EB2::BuildFC();
#endif
            }
        } else if (geometry == "sphere") {
            auto ProbLoArr = geom[max_level].ProbLoArray();
            auto ProbHiArr = geom[max_level].ProbHiArray();
            const Real xcen = myhalf * (ProbLoArr[0] + ProbHiArr[0]);
            const Real ycen = myhalf * (ProbLoArr[1] + ProbHiArr[1]);
            RealArray sphere_center = {xcen, ycen, zero};
            EB2::SphereIF implicit_fun(myhalf, sphere_center, false);
            auto gshop = EB2::makeShop(implicit_fun);
            if (build_eb_for_multigrid) {
                EB2::Build(gshop, geom[max_level], max_level, max_coarsening_level,
                            ngrow_for_eb, build_coarse_level_by_coarsening);
            } else {
                EB2::Build(gshop, this->Geom(), ngrow_for_eb);
#if USE_FC_FACTORY
                EB2::BuildFC();
#endif
            }
        }
    }

    if ( solverChoice.buildings_type == BuildingsType::ImmersedForcing) {
        constexpr int ngrow_for_eb = 4;
        if (geometry == "terrain") {
            Box buildings_bx(surroundingNodes(geom[max_level].Domain())); buildings_bx.grow(3);
            FArrayBox buildings_fab(makeSlab(buildings_bx,2,0),1);
            double dummy_time = 0.0;
            prob->init_buildings_surface(geom[max_level], buildings_fab, dummy_time);
            TerrainIF implicit_fun(buildings_fab, geom[max_level], stretched_dz_d[max_level]);
            auto gshop = EB2::makeShop(implicit_fun);
            EB2::Build(gshop, this->Geom(), ngrow_for_eb);
#if USE_FC_FACTORY
            EB2::BuildFC();
#endif
        } else if (geometry == "plane") {
            amrex::Abort("plane geometry is not supported with ImmersedForcing for buildings");
        } else if (geometry == "box") {
            RealArray box_lo{zero, zero, zero};
            RealArray box_hi{zero, zero, zero};
            pp_eb2.query("box_lo", box_lo);
            pp_eb2.query("box_hi", box_hi);
            EB2::BoxIF implicit_fun(box_lo, box_hi, false);
            auto gshop = EB2::makeShop(implicit_fun);
            EB2::Build(gshop, this->Geom(), ngrow_for_eb);
#if USE_FC_FACTORY
            EB2::BuildFC();
#endif
        } else if (geometry == "sphere") {
            amrex::Abort("sphere geometry is not supported with ImmersedForcing for buildings");
        }
    }

    forecast_state_1.resize(nlevs_max);
    forecast_state_2.resize(nlevs_max);
    forecast_state_interp.resize(nlevs_max);

    surface_state_1.resize(nlevs_max);
    surface_state_2.resize(nlevs_max);
    surface_state_interp.resize(nlevs_max);
}

ERF::~ERF () = default;
