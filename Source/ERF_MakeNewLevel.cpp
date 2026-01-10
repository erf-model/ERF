/**
 * \file ERF_MakeNewLevel.cpp
 */

/**
 * Routines that make and fill new levels after initialization, restart or regridding
 * The routines here call common routines in ERF_MakeNewArrays.cpp
*/

#include <memory>

#include "AMReX_buildInfo.H"

#include "ERF.H"
#include "ERF_Utils.H"
#include "ERF_ProbCommon.H"

using namespace amrex;

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// This is called both for initialization and for restart
// (overrides the pure virtual function in AmrCore)
// main.cpp --> ERF::InitData --> InitFromScratch --> MakeNewGrids --> MakeNewLevelFromScratch
//                                       restart  --> MakeNewGrids --> MakeNewLevelFromScratch
void ERF::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba_in,
                                   const DistributionMapping& dm_in)
{
    BoxArray ba;
    DistributionMapping dm;
    Box domain(Geom(0).Domain());
    if (lev == 0 && restart_chkfile.empty() &&
        (max_grid_size[0][0] >= domain.length(0)) &&
        (max_grid_size[0][1] >= domain.length(1)) &&
        ba_in.size() != ParallelDescriptor::NProcs())
    {
        // We only decompose in z if max_grid_size_z indicates we should
        bool decompose_in_z = (max_grid_size[0][2] < domain.length(2));

        ba = ERFPostProcessBaseGrids(Geom(0).Domain(),decompose_in_z);
        dm = DistributionMapping(ba);
    } else {
        ba = ba_in;
        dm = dm_in;
    }

    // ********************************************************************************************
    // Define grids[lev] to be ba
    // ********************************************************************************************
    SetBoxArray(lev, ba);

    // ********************************************************************************************
    // Define dmap[lev] to be dm
    // ********************************************************************************************
    SetDistributionMap(lev, dm);

    if (verbose) {
        amrex::Print() <<            "BA FROM SCRATCH AT LEVEL " << lev << " " << ba << std::endl;
        // amrex::Print() <<" SIMPLIFIED BA FROM SCRATCH AT LEVEL " << lev << " " << ba.simplified_list() << std::endl;
    }

    subdomains.resize(lev+1);
    if ( (lev == 0) || (
         (solverChoice.anelastic[lev] == 0) && (solverChoice.project_initial_velocity[lev] == 0) &&
         (solverChoice.init_type != InitType::WRFInput) && (solverChoice.init_type != InitType::Metgrid) ) ) {
        BoxArray dom(geom[lev].Domain());
        subdomains[lev].push_back(dom);
    } else {
        //
        // Create subdomains at each level within the domain such that
        // 1) all boxes in a given subdomain are "connected"
        // 2) no boxes in a subdomain touch any boxes in any other subdomain
        //
        make_subdomains(ba.simplified_list(), subdomains[lev]);
    }

    if (lev == 0) init_bcs();

    if ( solverChoice.terrain_type == TerrainType::EB ||
         solverChoice.terrain_type == TerrainType::ImmersedForcing ||
         solverChoice.buildings_type == BuildingsType::ImmersedForcing)
    {
        const amrex::EB2::IndexSpace& ebis = amrex::EB2::IndexSpace::top();
        const EB2::Level& eb_level = ebis.getLevel(geom[lev]);
        if (solverChoice.terrain_type == TerrainType::EB) {
            eb[lev]->make_all_factories(lev, geom[lev], grids[lev], dmap[lev], eb_level);
        } else if (solverChoice.terrain_type == TerrainType::ImmersedForcing ||
                   solverChoice.buildings_type == BuildingsType::ImmersedForcing) {
            eb[lev]->make_cc_factory(lev, geom[lev], grids[lev], dmap[lev], eb_level);
        }
    }

    auto& lev_new = vars_new[lev];
    auto& lev_old = vars_old[lev];

    //********************************************************************************************
    // This allocates all kinds of things, including but not limited to: solution arrays,
    //      terrain arrays, metric terms and base state.
    // *******************************************************************************************
    init_stuff(lev, ba, dm, lev_new, lev_old, base_state[lev], z_phys_nd[lev]);

    //********************************************************************************************
    // Land Surface Model
    // *******************************************************************************************
    int lsm_data_size  = lsm.Get_Data_Size();
    int lsm_flux_size  = lsm.Get_Flux_Size();
    lsm_data[lev].resize(lsm_data_size);
    lsm_data_name.resize(lsm_data_size);
    lsm_flux[lev].resize(lsm_flux_size);
    lsm_flux_name.resize(lsm_flux_size);
    lsm.Define(lev, solverChoice);
    if (solverChoice.lsm_type != LandSurfaceType::None)
    {
        lsm.Init(lev, vars_new[lev][Vars::cons], Geom(lev), 0.0); // dummy dt value
    }
    for (int mvar(0); mvar<lsm_data[lev].size(); ++mvar) {
        lsm_data[lev][mvar] = lsm.Get_Data_Ptr(lev,mvar);
        lsm_data_name[mvar] = lsm.Get_DataName(mvar);
    }
    for (int mvar(0); mvar<lsm_flux[lev].size(); ++mvar) {
        lsm_flux[lev][mvar] = lsm.Get_Flux_Ptr(lev,mvar);
        lsm_flux_name[mvar] = lsm.Get_FluxName(mvar);
    }



    // ********************************************************************************************
    // Build the data structures for calculating diffusive/turbulent terms
    // ********************************************************************************************
    update_diffusive_arrays(lev, ba, dm);

    // ********************************************************************************************
    // Build the data structures for holding sea surface temps and skin temps
    // ********************************************************************************************
    sst_lev[lev].resize(1);     sst_lev[lev][0] = nullptr;
    tsk_lev[lev].resize(1);     tsk_lev[lev][0] = nullptr;

    // ********************************************************************************************
    // Thin immersed body
    // *******************************************************************************************
    init_thin_body(lev, ba, dm);

    // ********************************************************************************************
    // Initialize the integrator class
    // ********************************************************************************************
    initialize_integrator(lev, lev_new[Vars::cons],lev_new[Vars::xvel]);

    // ********************************************************************************************
    // Initialize the data itself
    // If (init_type == InitType::WRFInput) then we are initializing terrain and the initial data in
    //                                      the same call so we must call init_only before update_terrain_arrays
    // If (init_type != InitType::WRFInput) then we want to initialize the terrain before the initial data
    //                                      since we may need to use the grid information before constructing
    //                                      initial idealized data
    // ********************************************************************************************
    if (restart_chkfile.empty()) {
        if ( (solverChoice.init_type == InitType::WRFInput) || (solverChoice.init_type == InitType::Metgrid) )
        {
            AMREX_ALWAYS_ASSERT(solverChoice.terrain_type == TerrainType::StaticFittedMesh);
            init_only(lev, time);
            init_zphys(lev, time);
            update_terrain_arrays(lev);
            make_physbcs(lev);
        } else {
            init_zphys(lev, time);
            update_terrain_arrays(lev);
            // Note that for init_type != InitType::WRFInput and != InitType::Metgrid,
            // make_physbcs is called inside init_only
            init_only(lev, time);
        }
    } else {
        // if restarting and nudging from input sounding, load the input sounding files
        if (lev == 0 && solverChoice.init_type == InitType::Input_Sounding && solverChoice.nudging_from_input_sounding)
        {
            if (input_sounding_data.input_sounding_file.empty()) {
                Error("input_sounding file name must be provided via input");
            }

            input_sounding_data.resize_arrays();

            // this will interpolate the input profiles to the nominal height levels
            // (ranging from 0 to the domain top)
            for (int n = 0; n < input_sounding_data.n_sounding_files; n++) {
                input_sounding_data.read_from_file(geom[lev], zlevels_stag[lev], n);
            }

            // this will calculate the hydrostatically balanced density and pressure
            // profiles following WRF ideal.exe
            if (solverChoice.sounding_type == SoundingType::Ideal) {
                input_sounding_data.calc_rho_p(0);
            } else if (solverChoice.sounding_type == SoundingType::Isentropic ||
                       solverChoice.sounding_type == SoundingType::DryIsentropic) {
                input_sounding_data.assume_dry = (solverChoice.sounding_type == SoundingType::DryIsentropic);
                input_sounding_data.calc_rho_p_isentropic(0);
            }
        }

        // We re-create terrain_blanking on restart rather than storing it in the checkpoint
        if (solverChoice.terrain_type == TerrainType::ImmersedForcing ||
            solverChoice.buildings_type == BuildingsType::ImmersedForcing) {
            int ngrow = ComputeGhostCells(solverChoice) + 2;
            terrain_blanking[lev]->setVal(1.0);
            MultiFab::Subtract(*terrain_blanking[lev], EBFactory(lev).getVolFrac(), 0, 0, 1, ngrow);
            terrain_blanking[lev]->FillBoundary(geom[lev].periodicity());
        }
    }

     // Read in tables needed for windfarm simulations
    // fill in Nturb multifab - number of turbines in each mesh cell
    // write out the vtk files for wind turbine location and/or
    // actuator disks
    #ifdef ERF_USE_WINDFARM
        init_windfarm(lev);
    #endif

    // ********************************************************************************************
    // Build the data structures for canopy model (depends upon z_phys)
    // ********************************************************************************************
    if (restart_chkfile.empty()) {
        if (solverChoice.do_forest_drag) {
            m_forest_drag[lev]->define_drag_field(ba, dm, geom[lev], z_phys_cc[lev].get(), z_phys_nd[lev].get());
        }
    }

    //********************************************************************************************
    // Create wall distance field for RANS model (depends upon z_phys)
    // *******************************************************************************************
    if (solverChoice.turbChoice[lev].rans_type != RANSType::None) {
        // Handle bottom boundary
        poisson_wall_dist(lev);

        // Correct the wall distance for immersed bodies
        if (solverChoice.advChoice.have_zero_flux_faces) {
            thinbody_wall_dist(walldist[lev],
                               solverChoice.advChoice.zero_xflux,
                               solverChoice.advChoice.zero_yflux,
                               solverChoice.advChoice.zero_zflux,
                               geom[lev],
                               z_phys_cc[lev]);
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

    //********************************************************************************************
    // Radiation
    // *******************************************************************************************
    if (solverChoice.rad_type != RadiationType::None)
    {
        rad[lev]->Init(geom[lev], ba, &vars_new[lev][Vars::cons]);
    }

    // ********************************************************************************************
    // If we are making a new level then the FillPatcher for this level hasn't been allocated yet
    // ********************************************************************************************
    if (lev > 0 && cf_width >= 0) {
        Construct_ERFFillPatchers(lev);
           Define_ERFFillPatchers(lev);
    }

#ifdef ERF_USE_PARTICLES
    if (restart_chkfile.empty()) {
        if (lev == 0) {
            initializeTracers((ParGDBBase*)GetParGDB(),z_phys_nd,time);
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

    if (verbose) {
        amrex::Print() <<" NEW BA FROM COARSE AT LEVEL " << lev << " " << ba << std::endl;
    }

    //
    // Grow the subdomains vector and build the subdomains vector at this level
    //
    subdomains.resize(lev+1);
    //
    // Create subdomains at each level within the domain such that
    // 1) all boxes in a given subdomain are "connected"
    // 2) no boxes in a subdomain touch any boxes in any other subdomain
    //
    if ( (solverChoice.anelastic[lev] == 0) && (solverChoice.project_initial_velocity[lev] == 0) ) {
        BoxArray dom(geom[lev].Domain());
        subdomains[lev].push_back(dom);
    } else {
        make_subdomains(ba.simplified_list(), subdomains[lev]);
    }

    if (lev == 0) init_bcs();

    //********************************************************************************************
    // This allocates all kinds of things, including but not limited to: solution arrays,
    //      terrain arrays, ba2d, metric terms and base state.
    // *******************************************************************************************
    init_stuff(lev, ba, dm, vars_new[lev], vars_old[lev], base_state[lev], z_phys_nd[lev]);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    // ********************************************************************************************
    // Build the data structures for metric quantities used with terrain-fitted coordinates
    // ********************************************************************************************
    if ( solverChoice.terrain_type == TerrainType::EB ||
         solverChoice.terrain_type == TerrainType::ImmersedForcing ||
         solverChoice.buildings_type == BuildingsType::ImmersedForcing)
    {
        const amrex::EB2::IndexSpace& ebis = amrex::EB2::IndexSpace::top();
        const EB2::Level& eb_level = ebis.getLevel(geom[lev]);
        if (solverChoice.terrain_type == TerrainType::EB) {
            eb[lev]->make_all_factories(lev, geom[lev], ba, dm, eb_level);
        } else if (solverChoice.terrain_type == TerrainType::ImmersedForcing ||
                   solverChoice.buildings_type == BuildingsType::ImmersedForcing) {
            eb[lev]->make_cc_factory(lev, geom[lev], ba, dm, eb_level);
        }
    }
    init_zphys(lev, time);
    update_terrain_arrays(lev);

    //
    // Make sure that detJ and z_phys_cc are the average of the data on a finer level if there is one
    //     *and* if there is two-way coupling
    //
    if ( (SolverChoice::mesh_type != MeshType::ConstantDz) && (solverChoice.coupling_type == CouplingType::TwoWay) ) {
        for (int crse_lev = lev-1; crse_lev >= 0; crse_lev--) {
            average_down(  *detJ_cc[crse_lev+1],   *detJ_cc[crse_lev], 0, 1, refRatio(crse_lev));
            average_down(*z_phys_cc[crse_lev+1], *z_phys_cc[crse_lev], 0, 1, refRatio(crse_lev));
        }
    }

    // ********************************************************************************************
    // Build the data structures for canopy model (depends upon z_phys)
    // ********************************************************************************************
    if (solverChoice.do_forest_drag) {
        m_forest_drag[lev]->define_drag_field(ba, dm, geom[lev], z_phys_cc[lev].get(), z_phys_nd[lev].get());
    }

    //********************************************************************************************
    // Radiation
    // *******************************************************************************************
    if (solverChoice.rad_type != RadiationType::None)
    {
        rad[lev]->Init(geom[lev], ba, &vars_new[lev][Vars::cons]);
    }

    // *****************************************************************************************************
    // Initialize the boundary conditions (after initializing the terrain but before calling
    //     initHSE or FillCoarsePatch)
    // *****************************************************************************************************
    make_physbcs(lev);

    // ********************************************************************************************
    // Update the base state at this level by interpolation from coarser level
    // ********************************************************************************************
    InterpFromCoarseLevel(base_state[lev], base_state[lev].nGrowVect(),
                          IntVect(0,0,0), // do not fill ghost cells outside the domain
                          base_state[lev-1], 0, 0, base_state[lev].nComp(),
                          geom[lev-1], geom[lev],
                          refRatio(lev-1), &cell_cons_interp,
                          domain_bcs_type, BCVars::cons_bc);

    // Impose bc's outside the domain
    (*physbcs_base[lev])(base_state[lev],0,base_state[lev].nComp(),base_state[lev].nGrowVect());

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
    // Build the data structures for calculating diffusive/turbulent terms
    // ********************************************************************************************
    update_diffusive_arrays(lev, ba, dm);

    // ********************************************************************************************
    // Build the data structures for holding sea surface temps and skin temps
    // ********************************************************************************************
    sst_lev[lev].resize(1);     sst_lev[lev][0] = nullptr;
    tsk_lev[lev].resize(1);     tsk_lev[lev][0] = nullptr;

    // ********************************************************************************************
    // Fill data at the new level by interpolation from the coarser level
    // Note that internal to FillCoarsePatch we will convert velocity to momentum,
    //      then interpolate momentum, then convert momentum back to velocity
    // Also note that FillCoarsePatch is hard-wired to act only on lev_new at coarse and fine
    // ********************************************************************************************

#ifdef ERF_USE_NETCDF
    if ( ( (solverChoice.init_type == InitType::WRFInput) || (solverChoice.init_type == InitType::Metgrid) ) &&
         !nc_init_file[lev].empty() )
    {
        // Just making sure that ghost cells aren't uninitialized...
        vars_new[lev][Vars::cons].setVal(0.0); vars_old[lev][Vars::cons].setVal(0.0);
        vars_new[lev][Vars::xvel].setVal(0.0); vars_old[lev][Vars::xvel].setVal(0.0);
        vars_new[lev][Vars::yvel].setVal(0.0); vars_old[lev][Vars::yvel].setVal(0.0);
        vars_new[lev][Vars::zvel].setVal(0.0); vars_old[lev][Vars::zvel].setVal(0.0);

        AMREX_ALWAYS_ASSERT(solverChoice.terrain_type == TerrainType::StaticFittedMesh);
        if (solverChoice.init_type == InitType::Metgrid) {
            init_from_metgrid(lev);
        } else if (solverChoice.init_type == InitType::WRFInput) {
            init_from_wrfinput(lev, *mf_C1H, *mf_C2H, *mf_MUB, *mf_PSFC[lev]);
        }
        init_zphys(lev, time);
        update_terrain_arrays(lev);
        make_physbcs(lev);

        dz_min[lev] = (*detJ_cc[lev]).min(0) * geom[lev].CellSize(2);

    } else {
#endif
    //
    // Interpolate the solution data
    //
    FillCoarsePatch(lev, time);
    //
    // Interpolate the 2D arrays at the lower boundary
    // Note that ba2d is constructed already in init_stuff, but we have not yet defined dmap[lev]
    // so we must explicitly pass dm.
    Interp2DArrays(lev,ba2d[lev],dm);
#ifdef ERF_USE_NETCDF
    }
#endif

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

    //********************************************************************************************
    // Land Surface Model
    // *******************************************************************************************
    int lsm_data_size  = lsm.Get_Data_Size();
    int lsm_flux_size  = lsm.Get_Flux_Size();
    lsm_data[lev].resize(lsm_data_size);
    lsm_data_name.resize(lsm_data_size);
    lsm_flux[lev].resize(lsm_flux_size);
    lsm_flux_name.resize(lsm_flux_size);
    lsm.Define(lev, solverChoice);
    if (solverChoice.lsm_type != LandSurfaceType::None)
    {
        lsm.Init(lev, vars_new[lev][Vars::cons], Geom(lev), 0.0); // dummy dt value
    }
    for (int mvar(0); mvar<lsm_data[lev].size(); ++mvar) {
        lsm_data[lev][mvar] = lsm.Get_Data_Ptr(lev,mvar);
        lsm_data_name[mvar] = lsm.Get_DataName(mvar);
    }
    for (int mvar(0); mvar<lsm_flux[lev].size(); ++mvar) {
        lsm_flux[lev][mvar] = lsm.Get_Flux_Ptr(lev,mvar);
        lsm_flux_name[mvar] = lsm.Get_FluxName(mvar);
    }

    // ********************************************************************************************
    // Create the SurfaceLayer arrays at this (new) level
    // ********************************************************************************************
    if (phys_bc_type[Orientation(Direction::z,Orientation::low)] == ERF_BC::surface_layer) {
        Vector<MultiFab*> mfv_old = {&vars_old[lev][Vars::cons], &vars_old[lev][Vars::xvel],
                                     &vars_old[lev][Vars::yvel], &vars_old[lev][Vars::zvel]};
        m_SurfaceLayer->make_SurfaceLayer_at_level(lev,lev+1,
                                                   mfv_old, Theta_prim[lev], Qv_prim[lev],
                                                   Qr_prim[lev], z_phys_nd[lev],
                                                   Hwave[lev].get(), Lwave[lev].get(), eddyDiffs_lev[lev].get(),
                                                   lsm_data[lev], lsm_data_name, lsm_flux[lev], lsm_flux_name,
                                                   sst_lev[lev], tsk_lev[lev], lmask_lev[lev]);
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
    if (verbose) {
        amrex::Print() <<" REMAKING WITH NEW BA AT LEVEL " << lev << " " << ba << std::endl;
    }

    AMREX_ALWAYS_ASSERT(solverChoice.terrain_type != TerrainType::MovingFittedMesh);

    BoxArray            ba_old(vars_new[lev][Vars::cons].boxArray());
    DistributionMapping dm_old(vars_new[lev][Vars::cons].DistributionMap());

    if (verbose) {
        amrex::Print() <<"               OLD BA AT LEVEL " << lev << " " << ba_old << std::endl;
    }

    //
    // Re-define subdomain at this level within the domain such that
    // 1) all boxes in a given subdomain are "connected"
    // 2) no boxes in a subdomain touch any boxes in any other subdomain
    //
    if (solverChoice.anelastic[lev] == 1) {
        make_subdomains(ba.simplified_list(), subdomains[lev]);
    }

    int     ncomp_cons  = vars_new[lev][Vars::cons].nComp();
    IntVect ngrow_state = vars_new[lev][Vars::cons].nGrowVect();

    int ngrow_vels = ComputeGhostCells(solverChoice);

    Vector<MultiFab> temp_lev_new(Vars::NumTypes);
    Vector<MultiFab> temp_lev_old(Vars::NumTypes);
    MultiFab temp_base_state;

    std::unique_ptr<MultiFab> temp_zphys_nd;

    //********************************************************************************************
    // This allocates all kinds of things, including but not limited to: solution arrays,
    //      terrain arrays and metrics, and base state.
    // *******************************************************************************************
    init_stuff(lev, ba, dm, temp_lev_new, temp_lev_old, temp_base_state, temp_zphys_nd);

    // ********************************************************************************************
    // Build the data structures for terrain-related quantities
    // ********************************************************************************************
    if ( solverChoice.terrain_type == TerrainType::EB ||
         solverChoice.terrain_type == TerrainType::ImmersedForcing ||
         solverChoice.buildings_type == BuildingsType::ImmersedForcing)
    {
        const amrex::EB2::IndexSpace& ebis = amrex::EB2::IndexSpace::top();
        const EB2::Level& eb_level = ebis.getLevel(geom[lev]);
        if (solverChoice.terrain_type == TerrainType::EB) {
            eb[lev]->make_all_factories(lev, geom[lev], ba, dm, eb_level);
        } else if (solverChoice.terrain_type == TerrainType::ImmersedForcing ||
                   solverChoice.buildings_type == BuildingsType::ImmersedForcing) {
            eb[lev]->make_cc_factory(lev, geom[lev], ba, dm, eb_level);
        }
    }
    remake_zphys(lev, time, temp_zphys_nd);
    update_terrain_arrays(lev);

    // ********************************************************************************************
    // Make sure that detJ and z_phys_cc are the average of the data on a finer level if there is one
    // Note that this shouldn't be necessary because the fine grid is created by interpolation
    // from the coarse ... but just in case ...
    // ********************************************************************************************
    if ( (SolverChoice::mesh_type != MeshType::ConstantDz) && (solverChoice.coupling_type == CouplingType::TwoWay) ) {
        for (int crse_lev = lev-1; crse_lev >= 0; crse_lev--) {
            average_down(  *detJ_cc[crse_lev+1],   *detJ_cc[crse_lev], 0, 1, refRatio(crse_lev));
            average_down(*z_phys_cc[crse_lev+1], *z_phys_cc[crse_lev], 0, 1, refRatio(crse_lev));
        }
    }

    // ********************************************************************************************
    // Build the data structures for canopy model (depends upon z_phys)
    // ********************************************************************************************
    if (solverChoice.do_forest_drag) {
        m_forest_drag[lev]->define_drag_field(ba, dm, geom[lev], z_phys_cc[lev].get(), z_phys_nd[lev].get());
    }

    // *****************************************************************************************************
    // Create the physbcs objects (after initializing the terrain but before calling FillCoarsePatch
    // *****************************************************************************************************
    make_physbcs(lev);

    // ********************************************************************************************
    // Update the base state at this level by interpolation from coarser level AND copy
    //    from previous (pre-regrid) base_state array
    // ********************************************************************************************
    if (lev > 0) {
        Interpolater* mapper = &cell_cons_interp;

        Vector<MultiFab*> fmf = {&base_state[lev  ], &base_state[lev  ]};
        Vector<MultiFab*> cmf = {&base_state[lev-1], &base_state[lev-1]};
        Vector<Real> ftime    = {time, time};
        Vector<Real> ctime    = {time, time};

        // Call FillPatch which ASSUMES that all ghost cells at lev-1 have already been filled
        FillPatchTwoLevels(temp_base_state, temp_base_state.nGrowVect(), IntVect(0,0,0),
                           time, cmf, ctime, fmf, ftime,
                           0, 0, temp_base_state.nComp(), geom[lev-1], geom[lev],
                           refRatio(lev-1), mapper, domain_bcs_type,
                           BaseBCVars::rho0_bc_comp);

        // Impose bc's outside the domain
        (*physbcs_base[lev])(temp_base_state,0,temp_base_state.nComp(),base_state[lev].nGrowVect());

        // *************************************************************************************************
        // This will fill the temporary MultiFabs with data from vars_new
        // NOTE: the momenta here are only used as scratch space, the momenta themselves are not fillpatched
        // NOTE: we must create the new base state before calling FillPatch because we will
        //       interpolate perturbational quantities
        // *************************************************************************************************
        FillPatchFineLevel(lev, time, {&temp_lev_new[Vars::cons],&temp_lev_new[Vars::xvel],
                           &temp_lev_new[Vars::yvel],&temp_lev_new[Vars::zvel]},
                          {&temp_lev_new[Vars::cons],&rU_new[lev],&rV_new[lev],&rW_new[lev]},
                           base_state[lev], temp_base_state, false);
    } else {
        temp_base_state.ParallelCopy(base_state[lev],0,0,base_state[lev].nComp(),
                                     base_state[lev].nGrowVect(),base_state[lev].nGrowVect());
        temp_lev_new[Vars::cons].ParallelCopy(vars_new[lev][Vars::cons],0,0,ncomp_cons,ngrow_state,ngrow_state);
        temp_lev_new[Vars::xvel].ParallelCopy(vars_new[lev][Vars::xvel],0,0,         1,ngrow_vels,ngrow_vels);
        temp_lev_new[Vars::yvel].ParallelCopy(vars_new[lev][Vars::yvel],0,0,         1,ngrow_vels,ngrow_vels);

        temp_lev_new[Vars::zvel].setVal(0.);
        temp_lev_new[Vars::zvel].ParallelCopy(vars_new[lev][Vars::zvel],0,0,         1,
                                              IntVect(ngrow_vels,ngrow_vels,0),IntVect(ngrow_vels,ngrow_vels,0));
    }

    // Now swap the pointers since we needed both old and new in the FillPatch
    std::swap(temp_base_state, base_state[lev]);

    // ********************************************************************************************
    // Copy from new into old just in case
    // ********************************************************************************************
    MultiFab::Copy(temp_lev_old[Vars::cons],temp_lev_new[Vars::cons],0,0,ncomp_cons,ngrow_state);
    MultiFab::Copy(temp_lev_old[Vars::xvel],temp_lev_new[Vars::xvel],0,0,         1,ngrow_vels);
    MultiFab::Copy(temp_lev_old[Vars::yvel],temp_lev_new[Vars::yvel],0,0,         1,ngrow_vels);
    MultiFab::Copy(temp_lev_old[Vars::zvel],temp_lev_new[Vars::zvel],0,0,         1,IntVect(ngrow_vels,ngrow_vels,0));

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

    //********************************************************************************************
    // Radiation
    // *******************************************************************************************
    if (solverChoice.rad_type != RadiationType::None)
    {
        rad[lev]->Init(geom[lev], ba, &vars_new[lev][Vars::cons]);
    }

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

    // ********************************************************************************************
    // Update the SurfaceLayer arrays at this level
    // ********************************************************************************************
    if (phys_bc_type[Orientation(Direction::z,Orientation::low)] == ERF_BC::surface_layer) {
        int nlevs = finest_level+1;
        Vector<MultiFab*> mfv_old = {&vars_old[lev][Vars::cons], &vars_old[lev][Vars::xvel],
                                     &vars_old[lev][Vars::yvel], &vars_old[lev][Vars::zvel]};
        m_SurfaceLayer->make_SurfaceLayer_at_level(lev,nlevs,
                                                   mfv_old, Theta_prim[lev], Qv_prim[lev],
                                                   Qr_prim[lev], z_phys_nd[lev],
                                                   Hwave[lev].get(),Lwave[lev].get(),eddyDiffs_lev[lev].get(),
                                                   lsm_data[lev], lsm_data_name, lsm_flux[lev], lsm_flux_name,
                                                   sst_lev[lev], tsk_lev[lev], lmask_lev[lev]);
    }

    // These calls are done in AmrCore::regrid if this is a regrid at lev > 0
    // For a level 0 regrid we must explicitly do them here
    if (lev == 0) {
        // Define grids[lev] to be ba
        SetBoxArray(lev, ba);

        // Define dmap[lev] to be dm
        SetDistributionMap(lev, dm);
    }

    // Clear the 2D arrays
    if (sst_lev[lev][0]) {
        for (int n = 0; n < sst_lev[lev].size(); n++) {
            sst_lev[lev][n].reset();
        }
    }
    if (tsk_lev[lev][0]) {
        for (int n = 0; n < tsk_lev[lev].size(); n++) {
            tsk_lev[lev][n].reset();
        }
    }
    if (lat_m[lev]) {
        lat_m[lev].reset();
    }
    if (lon_m[lev]) {
        lon_m[lev].reset();
    }
    if (sinPhi_m[lev]) {
        sinPhi_m[lev].reset();
    }
    if (cosPhi_m[lev]) {
        cosPhi_m[lev].reset();
    }

    //
    // Interpolate the 2D arrays at the lower boundary. We assume that since we created
    //     them by interpolation it is ok just to recreate them by interpolation.
    // Note that ba2d is constructed already in init_stuff, but we have not yet defined dmap[lev]
    //     so we must explicitly pass dm.
    Interp2DArrays(lev,ba2d[lev],dm);

#ifdef ERF_USE_PARTICLES
    particleData.Redistribute();
#endif
}

//
// Delete level data (overrides the pure virtual function in AmrCore)
// NOTE: this is only called for levels that no longer exist
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

    if (lev > 0) {
        zmom_crse_rhs[lev].clear();
    }

    if ( (solverChoice.anelastic[lev] == 1) || (solverChoice.project_initial_velocity[lev] == 1) ) {
        pp_inc[lev].clear();
    }
    if (solverChoice.anelastic[lev] == 0) {
        lagged_delta_rt[lev].clear();
    }
    avg_xmom[lev].clear();
    avg_ymom[lev].clear();
    avg_zmom[lev].clear();

    // Clears the integrator memory
    mri_integrator_mem[lev].reset();

    // Clears the physical boundary condition routines
    physbcs_cons[lev].reset();
    physbcs_u[lev].reset();
    physbcs_v[lev].reset();
    physbcs_w[lev].reset();
    physbcs_base[lev].reset();

    // Clears the flux register array
    advflux_reg[lev]->reset();

    // Clears the 2D arrays
    if (sst_lev[lev][0]) {
        for (int n = 0; n < sst_lev[lev].size(); n++) {
            sst_lev[lev][n].reset();
        }
    }
    if (tsk_lev[lev][0]) {
        for (int n = 0; n < tsk_lev[lev].size(); n++) {
            tsk_lev[lev][n].reset();
        }
    }
    if (lat_m[lev]) {
        lat_m[lev].reset();
    }
    if (lon_m[lev]) {
        lon_m[lev].reset();
    }
    if (sinPhi_m[lev]) {
        sinPhi_m[lev].reset();
    }
    if (cosPhi_m[lev]) {
        cosPhi_m[lev].reset();
    }
}

void
ERF::init_thin_body (int lev, const BoxArray& ba, const DistributionMapping& dm)
{
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
        amrex::Print() << "Setting up thin immersed body for "
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
        amrex::Print() << "Setting up thin immersed body for "
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
}
