/**
 * \file ERF_make_new_level.cpp
 */

/**
 * Routines that make and fill new levels after initialization, restart or regridding
*/

#include "prob_common.H"
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

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// This is called both for initialization and for restart
// (overrides the pure virtual function in AmrCore)
// main.cpp --> ERF::InitData --> InitFromScratch --> MakeNewGrids --> MakeNewLevelFromScratch
//                                       restart  --> MakeNewGrids --> MakeNewLevelFromScratch
void ERF::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
                                   const DistributionMapping& dm)
{
    // Set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);

    // The number of ghost cells for density must be 1 greater than that for velocity
    //     so that we can go back in forth betwen velocity and momentum on all faces
    int ngrow_state = ComputeGhostCells(solverChoice) + 1;
    int ngrow_vels  = ComputeGhostCells(solverChoice);

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
    qmoist[lev].define(ba, dm, 6, ngrow_state); // qv, qc, qi, qr, qs, qg
#endif

    // Build the data structures for calculating diffusive/turbulent terms
    update_arrays(lev, ba, dm);

    // ********************************************************************************************
    // Map factors
    // ********************************************************************************************
    BoxList bl2d = ba.boxList();
    for (auto& b : bl2d) {
        b.setRange(2,0);
    }
    BoxArray ba2d(std::move(bl2d));

    mapfac_m[lev] = std::make_unique<MultiFab>(ba2d,dm,1,3);
    mapfac_u[lev] = std::make_unique<MultiFab>(convert(ba2d,IntVect(1,0,0)),dm,1,3);
    mapfac_v[lev] = std::make_unique<MultiFab>(convert(ba2d,IntVect(0,1,0)),dm,1,3);
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
    base_state[lev].define(ba,dm,3,1);
    base_state[lev].setVal(0.);

    if (solverChoice.use_terrain && solverChoice.terrain_type > 0) {
        base_state_new[lev].define(ba,dm,3,1);
        base_state_new[lev].setVal(0.);
    }

    if (solverChoice.use_terrain) {
        z_phys_cc[lev] = std::make_unique<MultiFab>(ba,dm,1,1);
          detJ_cc[lev] = std::make_unique<MultiFab>(ba,dm,1,1);

        if (solverChoice.terrain_type > 0) {
            detJ_cc_new[lev] = std::make_unique<MultiFab>(ba,dm,1,1);
            detJ_cc_src[lev] = std::make_unique<MultiFab>(ba,dm,1,1);
            z_t_rk[lev] = std::make_unique<MultiFab>( convert(ba, IntVect(0,0,1)), dm, 1, 1 );
        }

        BoxArray ba_nd(ba);
        ba_nd.surroundingNodes();

        // We need this to be one greater than the ghost cells to handle levels > 0
        int ngrow = ComputeGhostCells(solverChoice) + 2;
        z_phys_nd[lev] = std::make_unique<MultiFab>(ba_nd,dm,1,IntVect(ngrow,ngrow,1));
        if (solverChoice.terrain_type > 0) {
            z_phys_nd_new[lev] = std::make_unique<MultiFab>(ba_nd,dm,1,IntVect(ngrow,ngrow,1));
            z_phys_nd_src[lev] = std::make_unique<MultiFab>(ba_nd,dm,1,IntVect(ngrow,ngrow,1));
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
    if (phys_bc_type[Orientation(Direction::z,Orientation::low)] == ERF_BC::MOST) {
      Theta_prim[lev] = std::make_unique<MultiFab>(ba,dm,1,IntVect(ngrow_state,ngrow_state,0));
    } else {
      Theta_prim[lev] = nullptr;
    }

    // ********************************************************************************************
    // Initialize the integrator class
    // ********************************************************************************************
    initialize_integrator(lev, lev_new[Vars::cons],lev_new[Vars::xvel]);

    // ********************************************************************************************
    // Initialize the data itself
    // If (init_type == "real") then we are initializing terrain and the initial data in
    //                          the same call so we must call init_only before update_terrain_arrays
    // If (init_type != "real") then we want to initialize the terrain before the initial data
    //                          since we may need to use the grid information before constructing
    //                          initial idealized data
    // ********************************************************************************************
    if ((init_type == "real") || (init_type == "metgrid")) {

        // If called from restart, the data structures for terrain-related quantities
        //    are built in the ReadCheckpoint routine.  Otherwise we build them here.
        if (restart_chkfile.empty()) {
            init_only(lev, start_time);
            update_terrain_arrays(lev, time);
        }

    } else {
        // Build the data structures for terrain-related quantities
        update_terrain_arrays(lev, time);

        // Initialize the solution data itself
        if (restart_chkfile.empty()) {
            init_only(lev, start_time);
        }
    }
}

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data (overrides the pure virtual function in AmrCore)
// regrid  --> RemakeLevel            (if level already existed)
// regrid  --> MakeNewLevelFromCoarse (if adding new level)
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

    // Build the data structures for terrain-related quantities
    update_terrain_arrays(lev, time);

    // Build the data structures for calculating diffusive/turbulent terms
    update_arrays(lev, ba, dm);

    define_grids_to_evolve(lev, grids[lev]);

    FillCoarsePatch(lev, time, {&lev_new[Vars::cons],&lev_new[Vars::xvel],
                                &lev_new[Vars::yvel],&lev_new[Vars::zvel]});

    initialize_integrator(lev, lev_new[Vars::cons], lev_new[Vars::xvel]);
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data (overrides the pure virtual function in AmrCore)
// regrid  --> RemakeLevel            (if level already existed)
// regrid  --> MakeNewLevelFromCoarse (if adding new level)
void
ERF::RemakeLevel (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm)
{
    define_grids_to_evolve(lev, ba);

    Vector<MultiFab> temp_lev_new(Vars::NumTypes);
    Vector<MultiFab> temp_lev_old(Vars::NumTypes);

    int ngrow_state = ComputeGhostCells(solverChoice) + 1;
    int ngrow_vels  = ComputeGhostCells(solverChoice);

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

    // Build the data structures for calculating diffusive/turbulent terms
    update_arrays(lev, ba, dm);

    // Build the data structures for terrain-related quantities
    update_terrain_arrays(lev, time);

    // ********************************************************************************************
    // Update the base state at this level
    // ********************************************************************************************
    base_state[lev].define(ba,dm,3,1);
    base_state[lev].setVal(0.);
    initHSE(lev);

    // ********************************************************************************************
    // Map factors
    // ********************************************************************************************
    BoxList bl2d = ba.boxList();
    for (auto& b : bl2d) {
        b.setRange(2,0);
    }
    BoxArray ba2d(std::move(bl2d));

    mapfac_m[lev] = std::make_unique<MultiFab>(ba2d,dm,1,3);
    mapfac_u[lev] = std::make_unique<MultiFab>(convert(ba2d,IntVect(1,0,0)),dm,1,3);
    mapfac_v[lev] = std::make_unique<MultiFab>(convert(ba2d,IntVect(0,1,0)),dm,1,3);
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

    initialize_integrator(lev, vars_new[lev][Vars::cons], vars_new[lev][Vars::xvel]);
}


void
ERF::update_arrays (int lev, const BoxArray& ba, const DistributionMapping& dm)
{
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

    if (l_use_diff) {
        Tau11_lev[lev] = std::make_unique<MultiFab>( ba  , dm, 1, IntVect(1,1,0) );
        Tau22_lev[lev] = std::make_unique<MultiFab>( ba  , dm, 1, IntVect(1,1,0) );
        Tau33_lev[lev] = std::make_unique<MultiFab>( ba  , dm, 1, IntVect(1,1,0) );
        Tau12_lev[lev] = std::make_unique<MultiFab>( ba12, dm, 1, IntVect(1,1,0) );
        Tau13_lev[lev] = std::make_unique<MultiFab>( ba13, dm, 1, IntVect(1,1,0) );
        Tau23_lev[lev] = std::make_unique<MultiFab>( ba23, dm, 1, IntVect(1,1,0) );
        if (l_use_terrain) {
            Tau21_lev[lev] = std::make_unique<MultiFab>( ba12, dm, 1, IntVect(1,1,0) );
            Tau31_lev[lev] = std::make_unique<MultiFab>( ba13, dm, 1, IntVect(1,1,0) );
            Tau32_lev[lev] = std::make_unique<MultiFab>( ba23, dm, 1, IntVect(1,1,0) );
        } else {
            Tau21_lev[lev] = nullptr;
            Tau31_lev[lev] = nullptr;
            Tau32_lev[lev] = nullptr;
        }
        SFS_hfx1_lev[lev] = std::make_unique<MultiFab>( ba  , dm, 1, IntVect(1,1,0) );
        SFS_hfx2_lev[lev] = std::make_unique<MultiFab>( ba  , dm, 1, IntVect(1,1,0) );
        SFS_hfx3_lev[lev] = std::make_unique<MultiFab>( ba  , dm, 1, IntVect(1,1,0) );
        SFS_diss_lev[lev] = std::make_unique<MultiFab>( ba  , dm, 1, IntVect(1,1,0) );
    } else {
      Tau11_lev[lev] = nullptr; Tau22_lev[lev] = nullptr; Tau33_lev[lev] = nullptr;
      Tau12_lev[lev] = nullptr; Tau21_lev[lev] = nullptr;
      Tau13_lev[lev] = nullptr; Tau31_lev[lev] = nullptr;
      Tau23_lev[lev] = nullptr; Tau32_lev[lev] = nullptr;
      SFS_hfx1_lev[lev] = nullptr; SFS_hfx2_lev[lev] = nullptr; SFS_hfx3_lev[lev] = nullptr;
      SFS_diss_lev[lev] = nullptr;
    }

    if (l_use_kturb) {
      eddyDiffs_lev[lev] = std::make_unique<MultiFab>( ba, dm, EddyDiff::NumDiffs, 1 );
      if(l_use_ddorf) {
          SmnSmn_lev[lev] = std::make_unique<MultiFab>( ba, dm, 1, 0 );
      } else {
          SmnSmn_lev[lev] = nullptr;
      }
    } else {
      eddyDiffs_lev[lev] = nullptr;
      SmnSmn_lev[lev]    = nullptr;
    }
}

void
ERF::update_terrain_arrays (int lev, Real time)
{
    if (solverChoice.use_terrain) {
        if (init_type != "real" && init_type != "metgrid") {
            init_custom_terrain(geom[lev],*z_phys_nd[lev],time);
            init_terrain_grid(geom[lev],*z_phys_nd[lev]);
        }
        if (lev>0 && (init_type == "real" || init_type == "metgrid")) {
            PhysBCFunctNoOp empty_bc;
            Vector<MultiFab*> fmf = {z_phys_nd[lev].get(), z_phys_nd[lev].get()};
            Vector<Real> ftime    = {t_old[lev], t_new[lev]};
            Vector<MultiFab*> cmf = {z_phys_nd[lev-1].get(), z_phys_nd[lev-1].get()};
            Vector<Real> ctime    = {t_old[lev-1], t_new[lev-1]};

            amrex::FillPatchTwoLevels(*z_phys_nd[lev], time, cmf, ctime, fmf, ftime,
                                      0, 0, 1, geom[lev-1], geom[lev],
                                      empty_bc, 0, empty_bc, 0, refRatio(lev-1),
                                      &node_bilinear_interp, domain_bcs_type, 0);
        }
        make_J(geom[lev],*z_phys_nd[lev],*detJ_cc[lev]);
        make_zcc(geom[lev],*z_phys_nd[lev],*z_phys_cc[lev]);
    }
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
    mri_integrator_mem[lev]->setNoSubstepping(solverChoice.no_substepping);
    mri_integrator_mem[lev]->setIncompressible(solverChoice.incompressible);
    mri_integrator_mem[lev]->setForceFirstStageSingleSubstep(solverChoice.force_stage1_single_substep);

    physbcs[lev] = std::make_unique<ERFPhysBCFunct> (lev, geom[lev], domain_bcs_type, domain_bcs_type_d,
                                                     solverChoice.terrain_type, m_bc_extdir_vals, m_bc_neumann_vals,
                                                     z_phys_nd[lev], detJ_cc[lev]);
}

//
// Delete level data (overrides the pure virtual function in AmrCore)
//
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
