/**
 * \file ERF_make_new_arrays.cpp
 */

/**
 * Worker routines for filling data at new levels after initialization, restart or regridding
*/

#include "prob_common.H"
#include <EOS.H>
#include <ERF.H>

#include <AMReX_buildInfo.H>

#include <Utils.H>
#include <TerrainMetrics.H>
#include <Utils/ParFunctions.H>
#include <memory>

#ifdef ERF_USE_MULTIBLOCK
#include <MultiBlockContainer.H>
#endif

using namespace amrex;

void
ERF::init_stuff (int lev, const BoxArray& ba, const DistributionMapping& dm,
                 Vector<MultiFab>& lev_new, Vector<MultiFab>& lev_old,
                 MultiFab& tmp_base_state,
                 std::unique_ptr<MultiFab>& tmp_zphys_nd)
{
    if (lev == 0) {
        min_k_at_level[lev] = 0;
        max_k_at_level[lev] = geom[lev].Domain().bigEnd(2);
    } else {
        // Start with unreasonable values so we compute the right min/max
        min_k_at_level[lev] = geom[lev].Domain().bigEnd(2);
        max_k_at_level[lev] = geom[lev].Domain().smallEnd(2);
        for (int n = 0; n < ba.size(); n++)
        {
            min_k_at_level[lev] = std::min(min_k_at_level[lev], ba[n].smallEnd(2));
            max_k_at_level[lev] = std::max(max_k_at_level[lev], ba[n].bigEnd(2));
        }
    }

    // ********************************************************************************************
    // Base state holds r_0, pres_0, pi_0 (in that order)
    // ********************************************************************************************
    tmp_base_state.define(ba,dm,3,1);
    tmp_base_state.setVal(0.);

    if (solverChoice.use_terrain && solverChoice.terrain_type != TerrainType::Static) {
        base_state_new[lev].define(ba,dm,3,1);
        base_state_new[lev].setVal(0.);
    }

    // ********************************************************************************************
    // Allocate terrain arrays
    // ********************************************************************************************
    if (solverChoice.use_terrain) {
        z_phys_cc[lev] = std::make_unique<MultiFab>(ba,dm,1,1);

        if (solverChoice.terrain_type != TerrainType::Static)
        {
            detJ_cc_new[lev] = std::make_unique<MultiFab>(ba,dm,1,1);
            detJ_cc_src[lev] = std::make_unique<MultiFab>(ba,dm,1,1);

            ax_src[lev] = std::make_unique<MultiFab>(convert(ba,IntVect(1,0,0)),dm,1,1);
            ay_src[lev] = std::make_unique<MultiFab>(convert(ba,IntVect(0,1,0)),dm,1,1);
            az_src[lev] = std::make_unique<MultiFab>(convert(ba,IntVect(0,0,1)),dm,1,1);

            ax_new[lev] = std::make_unique<MultiFab>(convert(ba,IntVect(1,0,0)),dm,1,1);
            ay_new[lev] = std::make_unique<MultiFab>(convert(ba,IntVect(0,1,0)),dm,1,1);
            az_new[lev] = std::make_unique<MultiFab>(convert(ba,IntVect(0,0,1)),dm,1,1);

            z_t_rk[lev] = std::make_unique<MultiFab>( convert(ba, IntVect(0,0,1)), dm, 1, 1 );
        }

        BoxArray ba_nd(ba);
        ba_nd.surroundingNodes();

        // We need this to be one greater than the ghost cells to handle levels > 0
        int ngrow = ComputeGhostCells(solverChoice.advChoice, solverChoice.use_NumDiff) + 2;
        tmp_zphys_nd = std::make_unique<MultiFab>(ba_nd,dm,1,IntVect(ngrow,ngrow,ngrow));

        if (solverChoice.terrain_type != TerrainType::Static) {
            z_phys_nd_new[lev] = std::make_unique<MultiFab>(ba_nd,dm,1,IntVect(ngrow,ngrow,ngrow));
            z_phys_nd_src[lev] = std::make_unique<MultiFab>(ba_nd,dm,1,IntVect(ngrow,ngrow,ngrow));
        }

    } else {
            z_phys_nd[lev] = nullptr;
            z_phys_cc[lev] = nullptr;

        z_phys_nd_new[lev] = nullptr;
          detJ_cc_new[lev] = nullptr;

        z_phys_nd_src[lev] = nullptr;
          detJ_cc_src[lev] = nullptr;

               z_t_rk[lev] = nullptr;
    }

   // We use these area arrays regardless of terrain, EB or none of the above
   detJ_cc[lev] = std::make_unique<MultiFab>(ba,dm,1,1);
        ax[lev] = std::make_unique<MultiFab>(convert(ba,IntVect(1,0,0)),dm,1,1);
        ay[lev] = std::make_unique<MultiFab>(convert(ba,IntVect(0,1,0)),dm,1,1);
        az[lev] = std::make_unique<MultiFab>(convert(ba,IntVect(0,0,1)),dm,1,1);

   detJ_cc[lev]->setVal(1.0);
        ax[lev]->setVal(1.0);
        ay[lev]->setVal(1.0);
        az[lev]->setVal(1.0);

    // ********************************************************************************************
    // These are the persistent containers for the old and new data
    // ********************************************************************************************
    int ncomp;
    if (lev > 0) {
        ncomp = vars_new[lev-1][Vars::cons].nComp();
    } else {
        int n_qstate   = micro->Get_Qstate_Size();
        ncomp = NVAR_max - (NMOIST_max - n_qstate);
    }

    // ********************************************************************************************
    // The number of ghost cells for density must be 1 greater than that for velocity
    //     so that we can go back in forth between velocity and momentum on all faces
    // ********************************************************************************************
    int ngrow_state = ComputeGhostCells(solverChoice.advChoice, solverChoice.use_NumDiff) + 1;
    int ngrow_vels  = ComputeGhostCells(solverChoice.advChoice, solverChoice.use_NumDiff);

    // ********************************************************************************************
    // New solution data containers
    // ********************************************************************************************
    lev_new[Vars::cons].define(ba, dm, ncomp, ngrow_state);
    lev_old[Vars::cons].define(ba, dm, ncomp, ngrow_state);

    lev_new[Vars::xvel].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);
    lev_old[Vars::xvel].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);

    lev_new[Vars::yvel].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);
    lev_old[Vars::yvel].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);

    // Note that we need the ghost cells in the z-direction if we are doing any
    // kind of domain decomposition in the vertical (at level 0 or above)
    lev_new[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, ngrow_vels);
    lev_old[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, ngrow_vels);

#ifdef ERF_USE_POISSON_SOLVE
    pp_inc[lev].define(ba, dm, 1, 1);
    pp_inc[lev].setVal(0.0);
#endif

    // ********************************************************************************************
    // These are just used for scratch in the time integrator but we might as well define them here
    // ********************************************************************************************
    rU_old[lev].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);
    rU_new[lev].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);

    rV_old[lev].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);
    rV_new[lev].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);

    rW_old[lev].define(convert(ba, IntVect(0,0,1)), dm, 1, ngrow_vels);
    rW_new[lev].define(convert(ba, IntVect(0,0,1)), dm, 1, ngrow_vels);

    // We do this here just so they won't be undefined in the initial FillPatch
    rU_new[lev].setVal(1.2e21);
    rV_new[lev].setVal(3.4e22);
    rW_new[lev].setVal(5.6e23);

    // ********************************************************************************************
    // These are just time averaged fields for diagnostics
    // ********************************************************************************************

    // NOTE: We are not completing a fillpach call on the time averaged data;
    //       which would copy on intersection and interpolate from coarse.
    //       Therefore, we are restarting the averaging when the ba changes,
    //       this may give poor statistics for dynamic mesh refinement.
    vel_t_avg[lev] = nullptr;
    if (solverChoice.time_avg_vel) {
        vel_t_avg[lev] = std::make_unique<MultiFab>(ba, dm, 4, 0); // Each vel comp and the mag
        vel_t_avg[lev]->setVal(0.0);
        t_avg_cnt[lev] = 0.0;
    }

    // ********************************************************************************************
    // Initialize flux registers whenever we create/re-create a level
    // ********************************************************************************************
    if (solverChoice.coupling_type == CouplingType::TwoWay) {
        if (lev == 0) {
            advflux_reg[0] = nullptr;
        } else {
            int ncomp_reflux = vars_new[0][Vars::cons].nComp();
            advflux_reg[lev] = new YAFluxRegister(ba       , grids[lev-1],
                                                  dm       ,  dmap[lev-1],
                                                  geom[lev],  geom[lev-1],
                                                  ref_ratio[lev-1], lev, ncomp_reflux);
        }
    }

    // ********************************************************************************************
    // Define Theta_prim storage if using MOST BC
    // ********************************************************************************************
    if (phys_bc_type[Orientation(Direction::z,Orientation::low)] == ERF_BC::MOST) {
        Theta_prim[lev] = std::make_unique<MultiFab>(ba,dm,1,IntVect(ngrow_state,ngrow_state,0));
        if (solverChoice.moisture_type != MoistureType::None) {
            Qv_prim[lev]    = std::make_unique<MultiFab>(ba,dm,1,IntVect(ngrow_state,ngrow_state,0));
            Qr_prim[lev]    = std::make_unique<MultiFab>(ba,dm,1,IntVect(ngrow_state,ngrow_state,0));
        } else {
            Qv_prim[lev]    = nullptr;
            Qr_prim[lev]    = nullptr;
        }
    } else {
        Theta_prim[lev] = nullptr;
        Qv_prim[lev]    = nullptr;
        Qr_prim[lev]    = nullptr;
    }

    // ********************************************************************************************
    // Map factors
    // ********************************************************************************************
    BoxList bl2d_mf = ba.boxList();
    for (auto& b : bl2d_mf) {
        b.setRange(2,0);
    }
    BoxArray ba2d_mf(std::move(bl2d_mf));

    mapfac_m[lev] = std::make_unique<MultiFab>(ba2d_mf,dm,1,3);
    mapfac_u[lev] = std::make_unique<MultiFab>(convert(ba2d_mf,IntVect(1,0,0)),dm,1,3);
    mapfac_v[lev] = std::make_unique<MultiFab>(convert(ba2d_mf,IntVect(0,1,0)),dm,1,3);
    if (solverChoice.test_mapfactor) {
        mapfac_m[lev]->setVal(0.5);
        mapfac_u[lev]->setVal(0.5);
        mapfac_v[lev]->setVal(0.5);
    }
    else {
        mapfac_m[lev]->setVal(1.);
        mapfac_u[lev]->setVal(1.);
        mapfac_v[lev]->setVal(1.);
    }

#if defined(ERF_USE_WINDFARM)
    //*********************************************************
    // Variables for Ftich model for windfarm parametrization
    //*********************************************************
    if (solverChoice.windfarm_type == WindFarmType::Fitch){
        vars_windfarm[lev].define(ba, dm, 5, ngrow_state); // V, dVabsdt, dudt, dvdt, dTKEdt
        Nturb[lev].define(ba, dm, 1, ngrow_state); // Number of turbines in a cell
    }
    if (solverChoice.windfarm_type == WindFarmType::EWP){
        vars_windfarm[lev].define(ba, dm, 3, ngrow_state); // dudt, dvdt, dTKEdt
        Nturb[lev].define(ba, dm, 1, ngrow_state); // Number of turbines in a cell
    }
    if (solverChoice.windfarm_type == WindFarmType::SimpleAD) {
        vars_windfarm[lev].define(ba, dm, 2, ngrow_state);// dudt, dvdt
        Nturb[lev].define(ba, dm, 1, ngrow_state); // Number of turbines in a cell
    }
#endif


#ifdef ERF_USE_WW3_COUPLING
    // create a new BoxArray and DistributionMapping for a MultiFab with 1 box
    BoxArray ba_onegrid(geom[lev].Domain());
    BoxList bl2d_onegrid = ba_onegrid.boxList();
    for (auto& b : bl2d_onegrid) {
        b.setRange(2,0);
    }
    BoxArray ba2d_onegrid(std::move(bl2d_onegrid));
    Vector<int> pmap;
    pmap.resize(1);
    pmap[0]=0;
    DistributionMapping dm_onegrid(ba2d_onegrid);
    dm_onegrid.define(pmap);

    Hwave_onegrid[lev] = std::make_unique<MultiFab>(ba2d_onegrid,dm_onegrid,1,IntVect(1,1,0));
    Lwave_onegrid[lev] = std::make_unique<MultiFab>(ba2d_onegrid,dm_onegrid,1,IntVect(1,1,0));

    BoxList bl2d_wave = ba.boxList();
    for (auto& b : bl2d_wave) {
        b.setRange(2,0);
    }
    BoxArray ba2d_wave(std::move(bl2d_wave));

    Hwave[lev] = std::make_unique<MultiFab>(ba2d_wave,dm,1,IntVect(3,3,0));
    Lwave[lev] = std::make_unique<MultiFab>(ba2d_wave,dm,1,IntVect(3,3,0));

    std::cout<<ba_onegrid<<std::endl;
    std::cout<<ba2d_onegrid<<std::endl;
    std::cout<<dm_onegrid<<std::endl;
#endif


#if defined(ERF_USE_RRTMGP)
    //*********************************************************
    // Radiation heating source terms
    //*********************************************************
    qheating_rates[lev] = std::make_unique<MultiFab>(ba, dm, 2, ngrow_state);
    qheating_rates[lev]->setVal(0.);

    //*********************************************************
    // Radiation fluxes for coupling to LSM
    //*********************************************************

    // NOTE: Finer levels do not need to coincide with the bottom domain boundary
    //       at k=0. We make slabs here with the kmin for a given box. Therefore,
    //       care must be taken before applying these fluxes to an LSM model. For

    // Radiative fluxes for LSM
    if (solverChoice.lsm_type != LandSurfaceType::None)
    {
        BoxList m_bl = ba.boxList();
        for (auto& b : m_bl) {
            int kmin = b.smallEnd(2);
            b.setRange(2,kmin);
        }
        BoxArray m_ba(std::move(m_bl));

        sw_lw_fluxes[lev] = std::make_unique<MultiFab>(m_ba, dm, 5, ngrow_state); // SW direct (2), SW diffuse (2), LW
        solar_zenith[lev] = std::make_unique<MultiFab>(m_ba, dm, 2, ngrow_state);

        sw_lw_fluxes[lev]->setVal(0.);
        solar_zenith[lev]->setVal(0.);
    }
#endif

    //*********************************************************
    // Turbulent perturbation region initialization
    //*********************************************************
    // TODO: Test perturbation on multiple levels
    if (solverChoice.pert_type == PerturbationType::perturbSource ||
        solverChoice.pert_type == PerturbationType::perturbDirect)
    {
        if (lev == 0) {
            turbPert.init_tpi(lev, geom[lev].Domain().bigEnd(), geom[lev].CellSizeArray(), ba, dm, ngrow_state);
        }
    }

    //
    // Define the land mask here and set it to all land by default
    // NOTE: the logic below will BREAK if we have any grids not touching the bottom boundary
    //
    {
    lmask_lev[lev].resize(1);
    auto ngv = lev_new[Vars::cons].nGrowVect(); ngv[2] = 0;
    BoxList bl2d_mask = ba.boxList();
    for (auto& b : bl2d_mask) {
        b.setRange(2,0);
    }
    BoxArray ba2d_mask(std::move(bl2d_mask));
    lmask_lev[lev][0] = std::make_unique<iMultiFab>(ba2d_mask,dm,1,ngv);
    lmask_lev[lev][0]->setVal(1);
    lmask_lev[lev][0]->FillBoundary(geom[lev].periodicity());
    }

    // Read in tables needed for windfarm simulations
    // fill in Nturb multifab - number of turbines in each mesh cell
    // write out the vtk files for wind turbine location and/or
    // actuator disks
    #ifdef ERF_USE_WINDFARM
        init_windfarm(lev);
    #endif
}

void
ERF::update_diffusive_arrays (int lev, const BoxArray& ba, const DistributionMapping& dm)
{
    // ********************************************************************************************
    // Diffusive terms
    // ********************************************************************************************
    bool l_use_terrain = solverChoice.use_terrain;
    bool l_use_diff    = ( (solverChoice.diffChoice.molec_diff_type != MolecDiffType::None) ||
                           (solverChoice.turbChoice[lev].les_type        !=       LESType::None) ||
                           (solverChoice.turbChoice[lev].pbl_type        !=       PBLType::None) );
    bool l_use_kturb   = ( (solverChoice.turbChoice[lev].les_type        != LESType::None)   ||
                           (solverChoice.turbChoice[lev].pbl_type        != PBLType::None) );
    bool l_use_ddorf   = (  solverChoice.turbChoice[lev].les_type        == LESType::Deardorff);
    bool l_use_moist   = (  solverChoice.moisture_type != MoistureType::None  );

    BoxArray ba12 = convert(ba, IntVect(1,1,0));
    BoxArray ba13 = convert(ba, IntVect(1,0,1));
    BoxArray ba23 = convert(ba, IntVect(0,1,1));

    if (l_use_diff) {
        //
        // NOTE: We require ghost cells in the vertical when allowing grids that don't
        //       cover the entire vertical extent of the domain at this level
        //
        Tau11_lev[lev] = std::make_unique<MultiFab>( ba  , dm, 1, IntVect(1,1,1) );
        Tau22_lev[lev] = std::make_unique<MultiFab>( ba  , dm, 1, IntVect(1,1,1) );
        Tau33_lev[lev] = std::make_unique<MultiFab>( ba  , dm, 1, IntVect(1,1,1) );
        Tau12_lev[lev] = std::make_unique<MultiFab>( ba12, dm, 1, IntVect(1,1,1) );
        Tau13_lev[lev] = std::make_unique<MultiFab>( ba13, dm, 1, IntVect(1,1,1) );
        Tau23_lev[lev] = std::make_unique<MultiFab>( ba23, dm, 1, IntVect(1,1,1) );
        if (l_use_terrain) {
            Tau21_lev[lev] = std::make_unique<MultiFab>( ba12, dm, 1, IntVect(1,1,1) );
            Tau31_lev[lev] = std::make_unique<MultiFab>( ba13, dm, 1, IntVect(1,1,1) );
            Tau32_lev[lev] = std::make_unique<MultiFab>( ba23, dm, 1, IntVect(1,1,1) );
        } else {
            Tau21_lev[lev] = nullptr;
            Tau31_lev[lev] = nullptr;
            Tau32_lev[lev] = nullptr;
        }
        SFS_hfx1_lev[lev] = std::make_unique<MultiFab>( convert(ba,IntVect(1,0,0)), dm, 1, IntVect(1,1,1) );
        SFS_hfx2_lev[lev] = std::make_unique<MultiFab>( convert(ba,IntVect(0,1,0)), dm, 1, IntVect(1,1,1) );
        SFS_hfx3_lev[lev] = std::make_unique<MultiFab>( convert(ba,IntVect(0,0,1)), dm, 1, IntVect(1,1,1) );
        SFS_diss_lev[lev] = std::make_unique<MultiFab>( ba  , dm, 1, IntVect(1,1,1) );
        SFS_hfx1_lev[lev]->setVal(0.);
        SFS_hfx2_lev[lev]->setVal(0.);
        SFS_hfx3_lev[lev]->setVal(0.);
        SFS_diss_lev[lev]->setVal(0.);
        if (l_use_moist) {
            SFS_q1fx3_lev[lev] = std::make_unique<MultiFab>( convert(ba,IntVect(0,0,1)), dm, 1, IntVect(1,1,1) );
            SFS_q2fx3_lev[lev] = std::make_unique<MultiFab>( convert(ba,IntVect(0,0,1)), dm, 1, IntVect(1,1,1) );
            SFS_q1fx3_lev[lev]->setVal(0.0);
            SFS_q2fx3_lev[lev]->setVal(0.0);
            if (solverChoice.use_rotate_most) {
                SFS_q1fx1_lev[lev] = std::make_unique<MultiFab>( convert(ba,IntVect(1,0,0)), dm, 1, IntVect(1,1,1) );
                SFS_q1fx2_lev[lev] = std::make_unique<MultiFab>( convert(ba,IntVect(0,1,0)), dm, 1, IntVect(1,1,1) );
                SFS_q1fx1_lev[lev]->setVal(0.0);
                SFS_q1fx2_lev[lev]->setVal(0.0);
            } else {
                SFS_q1fx1_lev[lev] = nullptr;
                SFS_q1fx2_lev[lev] = nullptr;
            }
        } else {
            SFS_q1fx1_lev[lev] = nullptr;
            SFS_q1fx2_lev[lev] = nullptr;
            SFS_q1fx3_lev[lev] = nullptr;
            SFS_q2fx3_lev[lev] = nullptr;
        }
    } else {
        Tau11_lev[lev] = nullptr; Tau22_lev[lev] = nullptr; Tau33_lev[lev] = nullptr;
        Tau12_lev[lev] = nullptr; Tau21_lev[lev] = nullptr;
        Tau13_lev[lev] = nullptr; Tau31_lev[lev] = nullptr;
        Tau23_lev[lev] = nullptr; Tau32_lev[lev] = nullptr;
        SFS_hfx1_lev[lev] = nullptr; SFS_hfx2_lev[lev] = nullptr; SFS_hfx3_lev[lev] = nullptr;
        SFS_diss_lev[lev] = nullptr;
    }

    if (l_use_kturb) {
        eddyDiffs_lev[lev] = std::make_unique<MultiFab>(ba, dm, EddyDiff::NumDiffs, 2);
        eddyDiffs_lev[lev]->setVal(0.0);
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
ERF::init_zphys (int lev, Real time)
{
    if (solverChoice.use_terrain) {
        //
        // First interpolate from coarser level if there is one
        //
        if (lev > 0) {
            Vector<MultiFab*> fmf = {z_phys_nd[lev].get(), z_phys_nd[lev].get()};
            Vector<Real> ftime    = {t_old[lev], t_new[lev]};
            Vector<MultiFab*> cmf = {z_phys_nd[lev-1].get(), z_phys_nd[lev-1].get()};
            Vector<Real> ctime    = {t_old[lev-1], t_new[lev-1]};

            //
            // First we fill z_phys_nd at lev>0 through interpolation
            //
            Interpolater* mapper = &node_bilinear_interp;
            PhysBCFunctNoOp null_bc;
            InterpFromCoarseLevel(*z_phys_nd[lev], time, *z_phys_nd[lev-1],
                                  0, 0, 1,
                                  geom[lev-1], geom[lev],
                                  null_bc, 0, null_bc, 0, refRatio(lev-1),
                                  mapper, domain_bcs_type, 0);
        }

        //
        // NOTE: we are NOT currently doing this as described -- we will simply use
        //       the interpolation from coarse done above
        // Then, if not using real/metgrid data,
        // 1) redefine the terrain at k=0 for every fine box which includes k=0
        // 2) recreate z_phys_nd at every fine node using
        // the data at the bottom of each fine grid
        // which has been either been interpolated from the coarse grid (k>0)
        // or set in init_custom_terrain (k=0)
        //
        if (lev == 0) {
            if (init_type != "real" && init_type != "metgrid")
            {
                prob->init_custom_terrain(geom[lev],*z_phys_nd[lev],time);

                Vector<Real> zmax(1); // only reduce at lev==0
                reduce_to_max_per_height(zmax, z_phys_nd[lev]);
                amrex::Print() << "Max terrain elevation = " << zmax[0] << std::endl;
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(zlevels_stag[zlevels_stag.size()-1] > zmax[0],
                    "Terrain is taller than domain top!");

                init_terrain_grid(lev,geom[lev],*z_phys_nd[lev],zlevels_stag,phys_bc_type);
            } // init_type
        } // lev == 0
    }
}

void
ERF::remake_zphys (int lev, Real time, std::unique_ptr<MultiFab>& temp_zphys_nd)
{
    if (solverChoice.use_terrain && lev > 0) {

        Vector<MultiFab*> fmf = {z_phys_nd[lev].get(), z_phys_nd[lev].get()};
        Vector<MultiFab*> cmf = {z_phys_nd[lev-1].get(), z_phys_nd[lev-1].get()};
        Vector<Real> ftime    = {time, time};
        Vector<Real> ctime    = {time, time};

        PhysBCFunctNoOp null_bc;
        Interpolater* mapper = &node_bilinear_interp;

        FillPatchTwoLevels(*temp_zphys_nd, time,
                           cmf, ctime, fmf, ftime,
                           0, 0, 1, geom[lev-1], geom[lev],
                           null_bc, 0, null_bc, 0, refRatio(lev-1),
                           mapper, domain_bcs_type, 0);

        std::swap(temp_zphys_nd, z_phys_nd[lev]);

    } // use_terrain && lev > 0
}

void
ERF::update_terrain_arrays (int lev)
{
    if (solverChoice.use_terrain) {
        make_J(geom[lev],*z_phys_nd[lev],*detJ_cc[lev]);
        make_areas(geom[lev],*z_phys_nd[lev],*ax[lev],*ay[lev],*az[lev]);
        make_zcc(geom[lev],*z_phys_nd[lev],*z_phys_cc[lev]);
    }
}

void
ERF::initialize_integrator (int lev, MultiFab& cons_mf, MultiFab& vel_mf)
{
    const BoxArray& ba(cons_mf.boxArray());
    const DistributionMapping& dm(cons_mf.DistributionMap());

    int ncomp_cons = cons_mf.nComp();

    // Initialize the integrator memory
    Vector<MultiFab> int_state; // integration state data structure example
    int_state.push_back(MultiFab(cons_mf, make_alias, 0, ncomp_cons));         // cons
    int_state.push_back(MultiFab(convert(ba,IntVect(1,0,0)), dm, 1, vel_mf.nGrow())); // xmom
    int_state.push_back(MultiFab(convert(ba,IntVect(0,1,0)), dm, 1, vel_mf.nGrow())); // ymom
    int_state.push_back(MultiFab(convert(ba,IntVect(0,0,1)), dm, 1, vel_mf.nGrow())); // zmom

    mri_integrator_mem[lev] = std::make_unique<MRISplitIntegrator<Vector<MultiFab> > >(int_state);
    mri_integrator_mem[lev]->setNoSubstepping(solverChoice.no_substepping);
    mri_integrator_mem[lev]->setIncompressible(solverChoice.incompressible[lev]);
    mri_integrator_mem[lev]->setNcompCons(ncomp_cons);
    mri_integrator_mem[lev]->setForceFirstStageSingleSubstep(solverChoice.force_stage1_single_substep);
}

void
ERF::initialize_bcs (int lev)
{
    // Dirichlet BC data only lives on level 0
    Real* u_bc_tmp(nullptr);
    Real* v_bc_tmp(nullptr);
    Real* w_bc_tmp(nullptr);
    if (lev==0) {
        u_bc_tmp = xvel_bc_data.data();
        v_bc_tmp = yvel_bc_data.data();
        w_bc_tmp = zvel_bc_data.data();
    }

    physbcs_cons[lev] = std::make_unique<ERFPhysBCFunct_cons> (lev, geom[lev], domain_bcs_type, domain_bcs_type_d,
                                                               m_bc_extdir_vals, m_bc_neumann_vals,
                                                               z_phys_nd[lev], use_real_bcs);
    physbcs_u[lev]    = std::make_unique<ERFPhysBCFunct_u> (lev, geom[lev], domain_bcs_type, domain_bcs_type_d,
                                                            m_bc_extdir_vals, m_bc_neumann_vals,
                                                            z_phys_nd[lev], use_real_bcs, u_bc_tmp);
    physbcs_v[lev]    = std::make_unique<ERFPhysBCFunct_v> (lev, geom[lev], domain_bcs_type, domain_bcs_type_d,
                                                            m_bc_extdir_vals, m_bc_neumann_vals,
                                                            z_phys_nd[lev], use_real_bcs, v_bc_tmp);
    physbcs_w[lev]    = std::make_unique<ERFPhysBCFunct_w> (lev, geom[lev], domain_bcs_type, domain_bcs_type_d,
                                                            m_bc_extdir_vals, m_bc_neumann_vals,
                                                            solverChoice.terrain_type, z_phys_nd[lev],
                                                            use_real_bcs, w_bc_tmp);
    physbcs_w_no_terrain[lev]    = std::make_unique<ERFPhysBCFunct_w_no_terrain>
                                                           (lev, geom[lev], domain_bcs_type, domain_bcs_type_d,
                                                            m_bc_extdir_vals, m_bc_neumann_vals, use_real_bcs,
                                                            w_bc_tmp);
}
