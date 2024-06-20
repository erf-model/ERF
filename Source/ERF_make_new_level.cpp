/**
 * \file ERF_make_new_level.cpp
 */

/**
 * Routines that make and fill new levels after initialization, restart or regridding
 * The routines here call common routines in ERF_make_new_arrays.cpp
*/

#include "prob_common.H"
#include <ERF.H>
#include <AMReX_buildInfo.H>
#include <Utils.H>
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

    // amrex::Print() <<" BA FROM SCRATCH AT LEVEL " << lev << " " << ba << std::endl;

#ifdef AMREX_USE_EB
    m_factory[lev] = makeEBFabFactory(geom[lev], grids[lev], dmap[lev],
                                      {nghost_eb_basic(),
                                       nghost_eb_volume(),
                                       nghost_eb_full()},
                                       EBSupport::full);
#else
    m_factory[lev] = std::make_unique<FArrayBoxFactory>();
#endif

    // The number of ghost cells for density must be 1 greater than that for velocity
    //     so that we can go back in forth betwen velocity and momentum on all faces
    // int ngrow_state = ComputeGhostCells(solverChoice.advChoice, solverChoice.use_NumDiff) + 1;
    // int ngrow_vels  = ComputeGhostCells(solverChoice.advChoice, solverChoice.use_NumDiff);

    auto& lev_new = vars_new[lev];
    auto& lev_old = vars_old[lev];

    //********************************************************************************************
    // This allocates all kinds of things, including but not limited to: solution arrays,
    //      terrain arrays and metrics,a nd base state.
    // *******************************************************************************************
    init_stuff(lev, ba, dm, lev_new, lev_old);

    //********************************************************************************************
    // Land Surface Model
    // *******************************************************************************************
    int lsm_size  = lsm.Get_Data_Size();
    lsm_data[lev].resize(lsm_size);
    lsm_flux[lev].resize(lsm_size);
    lsm.Define(lev, solverChoice);
    if (solverChoice.lsm_type != LandSurfaceType::None)
    {
        lsm.Init(lev, vars_new[lev][Vars::cons], Geom(lev), 0.0); // dummy dt value
    }
    for (int mvar(0); mvar<lsm_data[lev].size(); ++mvar) {
        lsm_data[lev][mvar] = lsm.Get_Data_Ptr(lev,mvar);
        lsm_flux[lev][mvar] = lsm.Get_Flux_Ptr(lev,mvar);
    }

    // ********************************************************************************************
    // Build the data structures for calculating diffusive/turbulent terms
    // ********************************************************************************************
    update_diffusive_arrays(lev, ba, dm);

    // ********************************************************************************************
    // Build the data structures for holding sea surface temps
    // ********************************************************************************************
    sst_lev[lev].resize(1);     sst_lev[lev][0] = nullptr;

    //********************************************************************************************
    // Thin immersed body
    // *******************************************************************************************
#if 0
    if ((solverChoice.advChoice.zero_xflux.size() > 0) ||
        (solverChoice.advChoice.zero_yflux.size() > 0) ||
        (solverChoice.advChoice.zero_zflux.size() > 0))
    {
        overset_imask[lev] = std::make_unique<iMultiFab>(ba,dm,1,0);
        overset_imask[lev]->setVal(1); // == value is unknown (to be solved)
    }
#endif

    if (solverChoice.advChoice.zero_xflux.size() > 0) {
        amrex::Print() << "Setting up thin immersed body for "
            << solverChoice.advChoice.zero_xflux.size() << " xfaces" << std::endl;
        BoxArray ba_xf(ba);
        ba_xf.surroundingNodes(0);
        thin_xforce[lev] = std::make_unique<MultiFab>(ba_xf,dm,1,0);
        thin_xforce[lev]->setVal(0.0);
        xflux_imask[lev] = std::make_unique<iMultiFab>(ba_xf,dm,1,0);
        xflux_imask[lev]->setVal(1);
        for ( MFIter mfi(*xflux_imask[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            Array4<int> const& imask_arr = xflux_imask[lev]->array(mfi);
            //Array4<int> const& imask_cell_arr = overset_imask[lev]->array(mfi);
            Box xbx = mfi.nodaltilebox(0);
            for (int iv=0; iv < solverChoice.advChoice.zero_xflux.size(); ++iv) {
                const auto& faceidx = solverChoice.advChoice.zero_xflux[iv];
                if ((faceidx[0] >= xbx.smallEnd(0)) && (faceidx[0] <= xbx.bigEnd(0)) &&
                    (faceidx[1] >= xbx.smallEnd(1)) && (faceidx[1] <= xbx.bigEnd(1)) &&
                    (faceidx[2] >= xbx.smallEnd(2)) && (faceidx[2] <= xbx.bigEnd(2)))
                {
                    imask_arr(faceidx[0],faceidx[1],faceidx[2]) = 0;
                    //imask_cell_arr(faceidx[0],faceidx[1],faceidx[2]) = 0;
                    //imask_cell_arr(faceidx[0]-1,faceidx[1],faceidx[2]) = 0;
                    amrex::AllPrint() << "  mask xface at " << faceidx << std::endl;
                }
            }
        }
    } else {
        thin_xforce[lev] = nullptr;
        xflux_imask[lev] = nullptr;
    }

    if (solverChoice.advChoice.zero_yflux.size() > 0) {
        amrex::Print() << "Setting up thin interface boundary for "
            << solverChoice.advChoice.zero_yflux.size() << " yfaces" << std::endl;
        BoxArray ba_yf(ba);
        ba_yf.surroundingNodes(1);
        thin_yforce[lev] = std::make_unique<MultiFab>(ba_yf,dm,1,0);
        thin_yforce[lev]->setVal(0.0);
        yflux_imask[lev] = std::make_unique<iMultiFab>(ba_yf,dm,1,0);
        yflux_imask[lev]->setVal(1);
        for ( MFIter mfi(*yflux_imask[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            Array4<int> const& imask_arr = yflux_imask[lev]->array(mfi);
            //Array4<int> const& imask_cell_arr = overset_imask[lev]->array(mfi);
            Box ybx = mfi.nodaltilebox(1);
            for (int iv=0; iv < solverChoice.advChoice.zero_yflux.size(); ++iv) {
                const auto& faceidx = solverChoice.advChoice.zero_yflux[iv];
                if ((faceidx[0] >= ybx.smallEnd(0)) && (faceidx[0] <= ybx.bigEnd(0)) &&
                    (faceidx[1] >= ybx.smallEnd(1)) && (faceidx[1] <= ybx.bigEnd(1)) &&
                    (faceidx[2] >= ybx.smallEnd(2)) && (faceidx[2] <= ybx.bigEnd(2)))
                {
                    imask_arr(faceidx[0],faceidx[1],faceidx[2]) = 0;
                    //imask_cell_arr(faceidx[0],faceidx[1],faceidx[2]) = 0;
                    //imask_cell_arr(faceidx[0],faceidx[1]-1,faceidx[2]) = 0;
                    amrex::AllPrint() << "  mask yface at " << faceidx << std::endl;
                }
            }
        }
    } else {
        thin_yforce[lev] = nullptr;
        yflux_imask[lev] = nullptr;
    }

    if (solverChoice.advChoice.zero_zflux.size() > 0) {
        amrex::Print() << "Setting up thin interface boundary for "
            << solverChoice.advChoice.zero_zflux.size() << " zfaces" << std::endl;
        BoxArray ba_zf(ba);
        ba_zf.surroundingNodes(2);
        thin_zforce[lev] = std::make_unique<MultiFab>(ba_zf,dm,1,0);
        thin_zforce[lev]->setVal(0.0);
        zflux_imask[lev] = std::make_unique<iMultiFab>(ba_zf,dm,1,0);
        zflux_imask[lev]->setVal(1);
        for ( MFIter mfi(*zflux_imask[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            Array4<int> const& imask_arr = zflux_imask[lev]->array(mfi);
            //Array4<int> const& imask_cell_arr = overset_imask[lev]->array(mfi);
            Box zbx = mfi.nodaltilebox(2);
            for (int iv=0; iv < solverChoice.advChoice.zero_zflux.size(); ++iv) {
                const auto& faceidx = solverChoice.advChoice.zero_zflux[iv];
                if ((faceidx[0] >= zbx.smallEnd(0)) && (faceidx[0] <= zbx.bigEnd(0)) &&
                    (faceidx[1] >= zbx.smallEnd(1)) && (faceidx[1] <= zbx.bigEnd(1)) &&
                    (faceidx[2] >= zbx.smallEnd(2)) && (faceidx[2] <= zbx.bigEnd(2)))
                {
                    imask_arr(faceidx[0],faceidx[1],faceidx[2]) = 0;
                    //imask_cell_arr(faceidx[0],faceidx[1],faceidx[2]) = 0;
                    //imask_cell_arr(faceidx[0],faceidx[1],faceidx[2]-1) = 0;
                    amrex::AllPrint() << "  mask zface at " << faceidx << std::endl;
                }
            }
        }
    } else {
        thin_zforce[lev] = nullptr;
        zflux_imask[lev] = nullptr;
    }

    //********************************************************************************************
    // Microphysics
    // *******************************************************************************************
    int q_size  = micro->Get_Qmoist_Size(lev);
    qmoist[lev].resize(q_size);
    micro->Define(lev, solverChoice);
    if (solverChoice.moisture_type != MoistureType::None)
    {
        micro->Init(lev, vars_new[lev][Vars::cons],
                    grids[lev], Geom(lev), 0.0,
                    z_phys_nd[lev], detJ_cc[lev]); // dummy dt value
    }
    for (int mvar(0); mvar<qmoist[lev].size(); ++mvar) {
        qmoist[lev][mvar] = micro->Get_Qmoist_Ptr(lev,mvar);
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

    // ********************************************************************************************
    // If we are making a new level then the FillPatcher for this level hasn't been allocated yet
    // ********************************************************************************************
    if (lev > 0 && cf_width >= 0) {
        Construct_ERFFillPatchers(lev);
           Define_ERFFillPatchers(lev);
    }

    // ********************************************************************************************
    // Initialize the boundary conditions
    // ********************************************************************************************
    initialize_bcs(lev);

#ifdef ERF_USE_PARTICLES
    if (restart_chkfile.empty()) {
        if (lev == 0) {
            initializeTracers((ParGDBBase*)GetParGDB(),z_phys_nd);
        } else {
            particleData.Redistribute();
        }
    }
#endif
}

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data (overrides the pure virtual function in AmrCore)
// regrid  --> RemakeLevel            (if level already existed)
// regrid  --> MakeNewLevelFromCoarse (if adding new level)
void
ERF::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
                             const DistributionMapping& dm)
{
    AMREX_ALWAYS_ASSERT(lev > 0);

    // amrex::Print() <<" NEW BA FROM COARSE AT LEVEL " << lev << " " << ba << std::endl;

    //********************************************************************************************
    // This allocates all kinds of things, including but not limited to: solution arrays,
    //      terrain arrays and metrics,a nd base state.
    // *******************************************************************************************
    init_stuff(lev, ba, dm, vars_new[lev], vars_old[lev]);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    // ********************************************************************************************
    // Build the data structures for terrain-related quantities
    // ********************************************************************************************
    update_terrain_arrays(lev, time);

    //
    // Make sure that detJ and z_phys_cc are the average of the data on a finer level if there is one
    //
    if (solverChoice.use_terrain != 0) {
        for (int crse_lev = lev-1; crse_lev >= 0; crse_lev--) {
            average_down(  *detJ_cc[crse_lev+1],   *detJ_cc[crse_lev], 0, 1, refRatio(crse_lev));
            average_down(*z_phys_cc[crse_lev+1], *z_phys_cc[crse_lev], 0, 1, refRatio(crse_lev));
        }
    }

    //********************************************************************************************
    // Microphysics
    // *******************************************************************************************
    int q_size  = micro->Get_Qmoist_Size(lev);
    qmoist[lev].resize(q_size);
    micro->Define(lev, solverChoice);
    if (solverChoice.moisture_type != MoistureType::None)
    {
        micro->Init(lev, vars_new[lev][Vars::cons],
                    grids[lev], Geom(lev), 0.0,
                    z_phys_nd[lev], detJ_cc[lev]); // dummy dt value
    }
    for (int mvar(0); mvar<qmoist[lev].size(); ++mvar) {
        qmoist[lev][mvar] = micro->Get_Qmoist_Ptr(lev,mvar);
    }

    // ********************************************************************************************
    // Update the base state at this level
    // ********************************************************************************************
    // base_state[lev].define(ba,dm,3,1);
    // base_state[lev].setVal(0.);
    initHSE(lev);

    // ********************************************************************************************
    // Build the data structures for calculating diffusive/turbulent terms
    // ********************************************************************************************
    update_diffusive_arrays(lev, ba, dm);

    // *****************************************************************************************************
    // Initialize the boundary conditions (after initializing the terrain but before calling FillCoarsePatch
    // *****************************************************************************************************
    initialize_bcs(lev);

    // ********************************************************************************************
    // Fill data at the new level by interpolation from the coarser level
    // Note that internal to FillCoarsePatch we will convert velocity to momentum,
    //      then interpolate momentum, then convert momentum back to velocity
    // Also note that FillCoarsePatch is hard-wired to act only on lev_new at coarse and fine
    // ********************************************************************************************
    FillCoarsePatch(lev, time);

    // ********************************************************************************************
    // Initialize the integrator class
    // ********************************************************************************************
    dt_mri_ratio[lev] = dt_mri_ratio[lev-1];
    initialize_integrator(lev, vars_new[lev][Vars::cons], vars_new[lev][Vars::xvel]);

    // ********************************************************************************************
    // If we are making a new level then the FillPatcher for this level hasn't been allocated yet
    // ********************************************************************************************
    if (lev > 0 && cf_width >= 0) {
        Construct_ERFFillPatchers(lev);
           Define_ERFFillPatchers(lev);
    }

#ifdef ERF_USE_PARTICLES
    // particleData.Redistribute();
#endif
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data (overrides the pure virtual function in AmrCore)
// regrid  --> RemakeLevel            (if level already existed)
// regrid  --> MakeNewLevelFromCoarse (if adding new level)
void
ERF::RemakeLevel (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm)
{
    // amrex::Print() <<" REMAKING WITH NEW BA AT LEVEL " << lev << " " << ba << std::endl;

    BoxArray            ba_old(vars_new[lev][Vars::cons].boxArray());
    DistributionMapping dm_old(vars_new[lev][Vars::cons].DistributionMap());

    int     ncomp_cons  = vars_new[lev][Vars::cons].nComp();
    IntVect ngrow_state = vars_new[lev][Vars::cons].nGrowVect();

    // int ngrow_state = ComputeGhostCells(solverChoice.advChoice, solverChoice.use_NumDiff) + 1;
    int ngrow_vels  = ComputeGhostCells(solverChoice.advChoice, solverChoice.use_NumDiff);

    Vector<MultiFab> temp_lev_new(Vars::NumTypes);
    Vector<MultiFab> temp_lev_old(Vars::NumTypes);

    //********************************************************************************************
    // This allocates all kinds of things, including but not limited to: solution arrays,
    //      terrain arrays and metrics, and base state.
    // *******************************************************************************************
    init_stuff(lev, ba, dm, temp_lev_new, temp_lev_old);

    // *****************************************************************************************************
    // Initialize the boundary conditions (after initializing the terrain but before calling FillCoarsePatch
    // *****************************************************************************************************
    initialize_bcs(lev);

    // ********************************************************************************************
    // This will fill the temporary MultiFabs with data from vars_new
    // ********************************************************************************************
    FillPatch(lev, time, {&temp_lev_new[Vars::cons],&temp_lev_new[Vars::xvel],
                          &temp_lev_new[Vars::yvel],&temp_lev_new[Vars::zvel]},
                         {&temp_lev_new[Vars::cons],&rU_new[lev],&rV_new[lev],&rW_new[lev]},
                          false);

    // ********************************************************************************************
    // Copy from new into old just in case
    // ********************************************************************************************
    MultiFab::Copy(temp_lev_old[Vars::cons],temp_lev_new[Vars::cons],0,0,ncomp_cons,ngrow_state);
    MultiFab::Copy(temp_lev_old[Vars::xvel],temp_lev_new[Vars::xvel],0,0,    1,ngrow_vels);
    MultiFab::Copy(temp_lev_old[Vars::yvel],temp_lev_new[Vars::yvel],0,0,    1,ngrow_vels);
    MultiFab::Copy(temp_lev_old[Vars::zvel],temp_lev_new[Vars::zvel],0,0,    1,IntVect(ngrow_vels,ngrow_vels,0));

    // ********************************************************************************************
    // Now swap the pointers
    // ********************************************************************************************
    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx) {
        std::swap(temp_lev_new[var_idx], vars_new[lev][var_idx]);
        std::swap(temp_lev_old[var_idx], vars_old[lev][var_idx]);
    }

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    // ********************************************************************************************
    // Build the data structures for calculating diffusive/turbulent terms
    // ********************************************************************************************
    update_diffusive_arrays(lev, ba, dm);

    // ********************************************************************************************
    // Build the data structures for terrain-related quantities
    // ********************************************************************************************
    update_terrain_arrays(lev, time);

    //
    // Make sure that detJ and z_phys_cc are the average of the data on a finer level if there is one
    //
    if (solverChoice.use_terrain != 0) {
        for (int crse_lev = lev-1; crse_lev >= 0; crse_lev--) {
            average_down(  *detJ_cc[crse_lev+1],   *detJ_cc[crse_lev], 0, 1, refRatio(crse_lev));
            average_down(*z_phys_cc[crse_lev+1], *z_phys_cc[crse_lev], 0, 1, refRatio(crse_lev));
        }
    }

    //********************************************************************************************
    // Microphysics
    // *******************************************************************************************
    int q_size = micro->Get_Qmoist_Size(lev);
    qmoist[lev].resize(q_size);
    micro->Define(lev, solverChoice);
    if (solverChoice.moisture_type != MoistureType::None)
    {
        micro->Init(lev, vars_new[lev][Vars::cons],
                    grids[lev], Geom(lev), 0.0,
                    z_phys_nd[lev], detJ_cc[lev]); // dummy dt value
    }
    for (int mvar(0); mvar<qmoist[lev].size(); ++mvar) {
        qmoist[lev][mvar] = micro->Get_Qmoist_Ptr(lev,mvar);
    }

    // ********************************************************************************************
    // Update the base state at this level
    // ********************************************************************************************
    // base_state[lev].define(ba,dm,3,1);
    // base_state[lev].setVal(0.);
    initHSE(lev);

    // ********************************************************************************************
    // Initialize the integrator class
    // ********************************************************************************************
    initialize_integrator(lev, vars_new[lev][Vars::cons], vars_new[lev][Vars::xvel]);

    // We need to re-define the FillPatcher if the grids have changed
    if (lev > 0 && cf_width >= 0) {
        bool ba_changed = (ba != ba_old);
        bool dm_changed = (dm != dm_old);
        if (ba_changed || dm_changed) {
          Define_ERFFillPatchers(lev);
        }
    }

#ifdef ERF_USE_PARTICLES
    particleData.Redistribute();
#endif
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

    base_state[lev].clear();

    rU_new[lev].clear();
    rU_old[lev].clear();
    rV_new[lev].clear();
    rV_old[lev].clear();
    rW_new[lev].clear();
    rW_old[lev].clear();

#ifdef ERF_USE_POISSON_SOLVE
    pp_inc[lev].clear();
#endif

    // Clears the integrator memory
    mri_integrator_mem[lev].reset();

    // Clears the physical boundary condition routines
    physbcs_cons[lev].reset();
    physbcs_u[lev].reset();
    physbcs_v[lev].reset();
    physbcs_w[lev].reset();
    physbcs_w_no_terrain[lev].reset();

    // Clears the flux register array
    advflux_reg[lev]->reset();
}
