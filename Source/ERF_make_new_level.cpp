/**
 * \file ERF_make_new_level.cpp
 */

/**
 * Routines that make and fill new levels after initialization, restart or regridding
 * The routines here call common routines in ERF_make_new_arrays.cpp
*/

#include "ERF_prob_common.H"
#include <ERF.H>
#include <AMReX_buildInfo.H>
#include <ERF_Utils.H>
#include <memory>

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

    // Define grids[lev] to be ba
    SetBoxArray(lev, ba);

    // Define dmap[lev] to be dm
    SetDistributionMap(lev, dm);

    // amrex::Print() <<" BA FROM SCRATCH AT LEVEL " << lev << " " << ba << std::endl;

    if (lev == 0) init_bcs();

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
    //     so that we can go back in forth between velocity and momentum on all faces
    // int ngrow_state = ComputeGhostCells(solverChoice.advChoice, solverChoice.use_NumDiff) + 1;
    // int ngrow_vels  = ComputeGhostCells(solverChoice.advChoice, solverChoice.use_NumDiff);

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

    // ********************************************************************************************
    // Thin immersed body
    // *******************************************************************************************
    init_immersed_body(lev, ba, dm);

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
    if (restart_chkfile.empty()) {
        if ((init_type == "real") || (init_type == "metgrid")) {
            init_only(lev, start_time);
            init_zphys(lev, time);
            update_terrain_arrays(lev);
            make_physbcs(lev);
        } else {
            init_zphys(lev, time);
            update_terrain_arrays(lev);
            // Note that for init_type != real or metgrid,
            // make_physbcs is called inside init_only
            init_only(lev, start_time);
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
    // If we are making a new level then the FillPatcher for this level hasn't been allocated yet
    // ********************************************************************************************
    if (lev > 0 && cf_width >= 0) {
        Construct_ERFFillPatchers(lev);
           Define_ERFFillPatchers(lev);
    }

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
    //      terrain arrays, metric terms and base state.
    // *******************************************************************************************
    init_stuff(lev, ba, dm, vars_new[lev], vars_old[lev], base_state[lev], z_phys_nd[lev]);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    // ********************************************************************************************
    // Build the data structures for terrain-related quantities
    // ********************************************************************************************
    init_zphys(lev, time);
    update_terrain_arrays(lev);

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
    // Update the base state at this level by interpolation from coarser level
    // ********************************************************************************************
    //
    // NOTE: this interpolater assumes that ALL ghost cells of the coarse MultiFab
    //       have been pre-filled - this includes ghost cells both inside and outside
    //       the domain
    //
    InterpFromCoarseLevel(base_state[lev], base_state[lev].nGrowVect(),
                          IntVect(0,0,0), // do not fill ghost cells outside the domain
                          base_state[lev-1], 0, 0, 3,
                          geom[lev-1], geom[lev],
                          refRatio(lev-1), &cell_cons_interp,
                          domain_bcs_type, BCVars::cons_bc);

    initHSE(lev);

    // ********************************************************************************************
    // Build the data structures for calculating diffusive/turbulent terms
    // ********************************************************************************************
    update_diffusive_arrays(lev, ba, dm);

    // *****************************************************************************************************
    // Initialize the boundary conditions (after initializing the terrain but before calling FillCoarsePatch
    // *****************************************************************************************************
    make_physbcs(lev);

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

    AMREX_ALWAYS_ASSERT(lev > 0);
    AMREX_ALWAYS_ASSERT(solverChoice.terrain_type != TerrainType::Moving);

    BoxArray            ba_old(vars_new[lev][Vars::cons].boxArray());
    DistributionMapping dm_old(vars_new[lev][Vars::cons].DistributionMap());

    // amrex::Print() <<"               OLD BA AT LEVEL " << lev << " " << ba_old << std::endl;

    int     ncomp_cons  = vars_new[lev][Vars::cons].nComp();
    IntVect ngrow_state = vars_new[lev][Vars::cons].nGrowVect();

    // int ngrow_state = ComputeGhostCells(solverChoice.advChoice, solverChoice.use_NumDiff) + 1;
    int ngrow_vels  = ComputeGhostCells(solverChoice.advChoice, solverChoice.use_NumDiff);

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
    remake_zphys(lev, temp_zphys_nd);
    update_terrain_arrays(lev);

    //
    // Make sure that detJ and z_phys_cc are the average of the data on a finer level if there is one
    //
    if (solverChoice.use_terrain != 0) {
        for (int crse_lev = lev-1; crse_lev >= 0; crse_lev--) {
            average_down(  *detJ_cc[crse_lev+1],   *detJ_cc[crse_lev], 0, 1, refRatio(crse_lev));
            average_down(*z_phys_cc[crse_lev+1], *z_phys_cc[crse_lev], 0, 1, refRatio(crse_lev));
        }
    }

    // *****************************************************************************************************
    // Create the physbcs objects (after initializing the terrain but before calling FillCoarsePatch
    // *****************************************************************************************************
    make_physbcs(lev);

    // ********************************************************************************************
    // This will fill the temporary MultiFabs with data from vars_new
    // ********************************************************************************************
    FillPatch(lev, time, {&temp_lev_new[Vars::cons],&temp_lev_new[Vars::xvel],
                          &temp_lev_new[Vars::yvel],&temp_lev_new[Vars::zvel]},
                         {&temp_lev_new[Vars::cons],&rU_new[lev],&rV_new[lev],&rW_new[lev]},
                          false);

    // ********************************************************************************************
    // Update the base state at this level by interpolation from coarser level AND copy
    //    from previous (pre-regrid) base_state array
    // ********************************************************************************************
    if (lev > 0) {
        // Interp all three components: rho, p, pi
        int  icomp = 0; int bccomp = 0; int  ncomp = 3;

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
                           BCVars::base_bc);

        // Impose bc's outside the domain
        (*physbcs_base[lev])(temp_base_state,icomp,ncomp,base_state[lev].nGrowVect(),time,bccomp);

        std::swap(temp_base_state, base_state[lev]);
    }

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

    if (solverChoice.anelastic[lev] == 1) {
        pp_inc[lev].clear();
    }

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
}

void
ERF::init_immersed_body (int lev, const BoxArray& ba, const DistributionMapping& dm)
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
}
