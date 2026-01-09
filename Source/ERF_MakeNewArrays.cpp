/**
 * \file ERF_MakeNewArrays.cpp
 */

/**
 * Worker routines for filling data at new levels after initialization, restart or regridding
*/

#include <memory>

#include "AMReX_buildInfo.H"

#include "ERF_ProbCommon.H"
#include "ERF_EOS.H"
#include "ERF.H"
#include "ERF_Utils.H"
#include "ERF_TerrainMetrics.H"
#include "ERF_ParFunctions.H"


using namespace amrex;

void
ERF::init_stuff (int lev, const BoxArray& ba, const DistributionMapping& dm,
                 Vector<MultiFab>& lev_new, Vector<MultiFab>& lev_old,
                 MultiFab& tmp_base_state,
                 std::unique_ptr<MultiFab>& tmp_zphys_nd)
{
    // ********************************************************************************************
    // Base state holds r_0, pres_0, pi_0, th_0 (in that order)
    //
    // Here is where we set the number of ghost cells for the base state!
    // ********************************************************************************************
    int ngb = (solverChoice.terrain_type == TerrainType::EB) ? 4 : 3;
    tmp_base_state.define(ba,dm,BaseState::num_comps,ngb);
    tmp_base_state.setVal(0.);

    if (solverChoice.terrain_type == TerrainType::MovingFittedMesh) {
        base_state_new[lev].define(ba,dm,BaseState::num_comps,base_state[lev].nGrowVect());
        base_state_new[lev].setVal(0.);
    }

    // ********************************************************************************************
    // Allocate terrain arrays
    // ********************************************************************************************

    BoxArray ba_nd(ba);
    ba_nd.surroundingNodes();

    // NOTE: this is where we actually allocate z_phys_nd -- but here it's called "tmp_zphys_nd"
    // We need this to be one greater than the ghost cells to handle levels > 0

    int ngrow = ComputeGhostCells(solverChoice) + 2;
    tmp_zphys_nd = std::make_unique<MultiFab>(ba_nd,dm,1,IntVect(ngrow,ngrow,ngrow));

    z_phys_cc[lev] = std::make_unique<MultiFab>(ba,dm,1,1);
    init_default_zphys(lev, geom[lev], *tmp_zphys_nd, *z_phys_cc[lev]);

    if (solverChoice.terrain_type == TerrainType::MovingFittedMesh)
    {
        detJ_cc_new[lev] = std::make_unique<MultiFab>(ba,dm,1,1);
        detJ_cc_src[lev] = std::make_unique<MultiFab>(ba,dm,1,1);

        ax_src[lev] = std::make_unique<MultiFab>(convert(ba,IntVect(1,0,0)),dm,1,1);
        ay_src[lev] = std::make_unique<MultiFab>(convert(ba,IntVect(0,1,0)),dm,1,1);
        az_src[lev] = std::make_unique<MultiFab>(convert(ba,IntVect(0,0,1)),dm,1,1);

        z_t_rk[lev] = std::make_unique<MultiFab>( convert(ba, IntVect(0,0,1)), dm, 1, 1 );

        z_phys_nd_new[lev] = std::make_unique<MultiFab>(ba_nd,dm,1,IntVect(ngrow,ngrow,ngrow));
        z_phys_nd_src[lev] = std::make_unique<MultiFab>(ba_nd,dm,1,IntVect(ngrow,ngrow,ngrow));
        z_phys_cc_src[lev] = std::make_unique<MultiFab>(ba,dm,1,1);
    }
    else
    {
        z_phys_nd_new[lev] = nullptr;
          detJ_cc_new[lev] = nullptr;

        z_phys_nd_src[lev] = nullptr;
        z_phys_cc_src[lev] = nullptr;
          detJ_cc_src[lev] = nullptr;

               z_t_rk[lev] = nullptr;
    }

    if (solverChoice.terrain_type == TerrainType::ImmersedForcing ||
        solverChoice.buildings_type == BuildingsType::ImmersedForcing)
    {
        terrain_blanking[lev] = std::make_unique<MultiFab>(ba,dm,1,ngrow);
        terrain_blanking[lev]->setVal(1.0);
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
    // Create wall distance array for RANS modeling
    // ********************************************************************************************
    if (solverChoice.turbChoice[lev].rans_type != RANSType::None) {
        walldist[lev] = std::make_unique<MultiFab>(ba,dm,1,1);
        walldist[lev]->setVal(1e23);
    } else {
        walldist[lev] = nullptr;
    }

    // ********************************************************************************************
    // These are the persistent containers for the old and new data
    // ********************************************************************************************
    int ncomp;
    if (lev > 0) {
        ncomp = vars_new[lev-1][Vars::cons].nComp();
    } else {
        int n_qstate   = micro->Get_Qstate_Size();
        ncomp = NDRY + NSCALARS + n_qstate;
    }

    // ********************************************************************************************
    // The number of ghost cells for density must be 1 greater than that for velocity
    //     so that we can go back in forth between velocity and momentum on all faces
    // ********************************************************************************************
    int ngrow_state = ComputeGhostCells(solverChoice) + 1;
    int ngrow_vels  = ComputeGhostCells(solverChoice);

    // ********************************************************************************************
    // New solution data containers
    // ********************************************************************************************
    if (solverChoice.terrain_type != TerrainType::EB) {
        lev_new[Vars::cons].define(ba, dm, ncomp, ngrow_state);
        lev_old[Vars::cons].define(ba, dm, ncomp, ngrow_state);
    } else {
        // EB: Define the MultiFabs with the EBFactory
        lev_new[Vars::cons].define(ba, dm, ncomp, ngrow_state, MFInfo(), EBFactory(lev));
        lev_old[Vars::cons].define(ba, dm, ncomp, ngrow_state, MFInfo(), EBFactory(lev));
    }
    lev_new[Vars::xvel].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);
    lev_old[Vars::xvel].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);

    lev_new[Vars::yvel].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);
    lev_old[Vars::yvel].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);

    gradp[lev][GpVars::gpx].define(convert(ba, IntVect(1,0,0)), dm, 1, 1); gradp[lev][GpVars::gpx].setVal(0.);
    gradp[lev][GpVars::gpy].define(convert(ba, IntVect(0,1,0)), dm, 1, 1); gradp[lev][GpVars::gpy].setVal(0.);
    gradp[lev][GpVars::gpz].define(convert(ba, IntVect(0,0,1)), dm, 1, 1); gradp[lev][GpVars::gpz].setVal(0.);

    // Note that we need the ghost cells in the z-direction if we are doing any
    // kind of domain decomposition in the vertical (at level 0 or above)
    lev_new[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, ngrow_vels);
    lev_old[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, ngrow_vels);

    if ( (solverChoice.anelastic[lev] == 1) || (solverChoice.project_initial_velocity[lev] == 1) ) {
        pp_inc[lev].define(ba, dm, 1, 1);
        pp_inc[lev].setVal(0.0);
    }

    // We use this in the fast substepping only
    if (solverChoice.anelastic[lev] == 0) {
        lagged_delta_rt[lev].define(ba, dm, 1, 1);
        lagged_delta_rt[lev].setVal(0.0);
    }

    // We use these for advecting the slow variables, whether anelastic or compressible
    avg_xmom[lev].define(convert(ba, IntVect(1,0,0)), dm, 1, 1);
    avg_ymom[lev].define(convert(ba, IntVect(0,1,0)), dm, 1, 1);
    avg_zmom[lev].define(convert(ba, IntVect(0,0,1)), dm, 1, 1);
    avg_xmom[lev].setVal(0.0); avg_ymom[lev].setVal(0.0); avg_zmom[lev].setVal(0.0);

    // ********************************************************************************************
    // These are just used for scratch in the time integrator but we might as well define them here
    // ********************************************************************************************
    if (solverChoice.terrain_type != TerrainType::EB) {
        rU_old[lev].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);
        rU_new[lev].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);

        rV_old[lev].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);
        rV_new[lev].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);

        rW_old[lev].define(convert(ba, IntVect(0,0,1)), dm, 1, ngrow_vels);
        rW_new[lev].define(convert(ba, IntVect(0,0,1)), dm, 1, ngrow_vels);
    } else {
        // EB: Define the MultiFabs with the EBFactory
        rU_old[lev].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels, MFInfo(), EBFactory(lev));
        rU_new[lev].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels, MFInfo(), EBFactory(lev));

        rV_old[lev].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels, MFInfo(), EBFactory(lev));
        rV_new[lev].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels, MFInfo(), EBFactory(lev));

        rW_old[lev].define(convert(ba, IntVect(0,0,1)), dm, 1, ngrow_vels, MFInfo(), EBFactory(lev));
        rW_new[lev].define(convert(ba, IntVect(0,0,1)), dm, 1, ngrow_vels, MFInfo(), EBFactory(lev));
    }

    if (lev > 0) {
        //xmom_crse_rhs[lev].define(convert(ba, IntVect(1,0,0)), dm, 1, IntVect{0});
        //ymom_crse_rhs[lev].define(convert(ba, IntVect(0,1,0)), dm, 1, IntVect{0});
        zmom_crse_rhs[lev].define(convert(ba, IntVect(0,0,1)), dm, 1, IntVect{0});
    }

    // We do this here just so they won't be undefined in the initial FillPatch
    rU_old[lev].setVal(1.2e21);
    rV_old[lev].setVal(3.4e22);
    rW_old[lev].setVal(5.6e23);
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
    // Define Theta_prim storage if using surface_layer BC
    // ********************************************************************************************
    if (phys_bc_type[Orientation(Direction::z,Orientation::low)] == ERF_BC::surface_layer) {
        Theta_prim[lev] = std::make_unique<MultiFab>(ba,dm,1,IntVect(ngrow_state,ngrow_state,1));
        if (solverChoice.moisture_type != MoistureType::None) {
            Qv_prim[lev]    = std::make_unique<MultiFab>(ba,dm,1,IntVect(ngrow_state,ngrow_state,1));
            Qr_prim[lev]    = std::make_unique<MultiFab>(ba,dm,1,IntVect(ngrow_state,ngrow_state,1));
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

    mapfac[lev].resize(MapFacType::num);
    mapfac[lev][MapFacType::m_x] = std::make_unique<MultiFab>(                        ba2d_mf,dm,1,IntVect(3,3,0));
    mapfac[lev][MapFacType::u_x] = std::make_unique<MultiFab>(convert(ba2d_mf,IntVect(1,0,0)),dm,1,IntVect(3,3,0));
    mapfac[lev][MapFacType::v_x] = std::make_unique<MultiFab>(convert(ba2d_mf,IntVect(0,1,0)),dm,1,IntVect(3,3,0));

#if 0
    // For now we comment this out to avoid CI failures but we will need to re-enable
    //     this if using non-conformal mappings
    if (MapFacType::m_y != MapFacType::m_x) {
        mapfac[lev][MapFacType::m_y] = std::make_unique<MultiFab>(ba2d_mf,dm,1,IntVect(3,3,0));
    }
    if (MapFacType::u_y != MapFacType::u_x) {
        mapfac[lev][MapFacType::u_y] = std::make_unique<MultiFab>(convert(ba2d_mf,IntVect(1,0,0)),dm,1,IntVect(3,3,0));
    }
    if (MapFacType::v_y != MapFacType::v_x) {
        mapfac[lev][MapFacType::v_y] = std::make_unique<MultiFab>(convert(ba2d_mf,IntVect(0,1,0)),dm,1,IntVect(3,3,0));
    }
#endif

    if (solverChoice.test_mapfactor) {
        for (int i = 0; i < 3; i++) {
            mapfac[lev][i]->setVal(0.5);
        }
        for (int i = 3; i < mapfac[lev].size(); i++) {
            mapfac[lev][i]->setVal(0.25);
        }
    } else {
        for (int i = 0; i < mapfac[lev].size(); i++) {
            mapfac[lev][i]->setVal(1.0);
        }
    }

    // ********************************************************************************************
    // Build 1D BA and 2D BA
    // ********************************************************************************************
    BoxList bl1d = ba.boxList();
    for (auto& b : bl1d) {
        b.setRange(0,0);
        b.setRange(1,0);
    }
    ba1d[lev]  = BoxArray(std::move(bl1d));

    // Build 2D BA
    BoxList bl2d = ba.boxList();
    for (auto& b : bl2d) {
        b.setRange(2,0);
    }
    ba2d[lev]  = BoxArray(std::move(bl2d));

    IntVect ng  = vars_new[lev][Vars::cons].nGrowVect();

    if (lev == 0) {
        mf_C1H = std::make_unique<MultiFab>(ba1d[lev],dm,1,IntVect(ng[0],ng[1],ng[2]));
        mf_C2H = std::make_unique<MultiFab>(ba1d[lev],dm,1,IntVect(ng[0],ng[1],ng[2]));
        mf_MUB = std::make_unique<MultiFab>(ba2d[lev],dm,1,IntVect(ng[0],ng[1],ng[2]));
    }

    mf_PSFC[lev] = std::make_unique<MultiFab>(ba2d[lev],dm,1,ng);

    //*********************************************************
    // Variables for Fitch model for windfarm parametrization
    //*********************************************************
#if defined(ERF_USE_WINDFARM)
    if (solverChoice.windfarm_type == WindFarmType::Fitch){
        vars_windfarm[lev].define(ba, dm, 5, ngrow_state); // V, dVabsdt, dudt, dvdt, dTKEdt
    }
    if (solverChoice.windfarm_type == WindFarmType::EWP){
        vars_windfarm[lev].define(ba, dm, 3, ngrow_state); // dudt, dvdt, dTKEdt
    }
    if (solverChoice.windfarm_type == WindFarmType::SimpleAD) {
        vars_windfarm[lev].define(ba, dm, 2, ngrow_state);// dudt, dvdt
    }
    if (solverChoice.windfarm_type == WindFarmType::GeneralAD) {
        vars_windfarm[lev].define(ba, dm, 3, ngrow_state);// dudt, dvdt, dwdt
    }
        Nturb[lev].define(ba, dm, 1, ngrow_state); // Number of turbines in a cell
        SMark[lev].define(ba, dm, 2, 1); // Free stream velocity/source term
                                                   // sampling marker in a cell - 2 components
#endif

    if(solverChoice.init_type == InitType::HindCast and
        solverChoice.hindcast_lateral_forcing) {

        int ncomp_extra = 2;
        int nvars = vars_new[lev].size();

        // Resize all containers
        forecast_state_1[lev].resize(nvars + 1);
        forecast_state_2[lev].resize(nvars + 1);
        forecast_state_interp[lev].resize(nvars + 1);

        // Define the "normal" components
        for (int comp = 0; comp < nvars; ++comp) {
            const MultiFab& src = vars_new[lev][comp];
            ncomp = src.nComp();
            ngrow = src.nGrow();

            forecast_state_1[lev][comp].define(ba, dm, ncomp, ng);
            forecast_state_2[lev][comp].define(ba, dm, ncomp, ng);
            forecast_state_interp[lev][comp].define(ba, dm, ncomp, ng);
        }

        // Define the "extra" component (last slot)
        {
            const MultiFab& src0 = vars_new[lev][0];
            ngrow = src0.nGrow();
            int idx = nvars;

            forecast_state_1[lev][idx].define(ba, dm, ncomp_extra, ngrow);
            forecast_state_2[lev][idx].define(ba, dm, ncomp_extra, ngrow);
            forecast_state_interp[lev][idx].define(ba, dm, ncomp_extra, ngrow);
        }
        bool regrid_forces_file_read = true;
        WeatherDataInterpolation(lev, t_new[0],z_phys_nd, regrid_forces_file_read);
    }


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


    //*********************************************************
    // Radiation heating source terms
    //*********************************************************
    if (solverChoice.rad_type != RadiationType::None)
    {
        qheating_rates[lev] = std::make_unique<MultiFab>(ba, dm, 2, 0);
        rad_fluxes[lev]     = std::make_unique<MultiFab>(ba, dm, 4, 0);
        qheating_rates[lev]->setVal(0.);
        rad_fluxes[lev]->setVal(0.);
    }

    //*********************************************************
    // Radiation fluxes for coupling to LSM
    //*********************************************************

    // NOTE: Finer levels do not need to coincide with the bottom domain boundary
    //       at k=0. We make slabs here with the kmin for a given box. Therefore,
    //       care must be taken before applying these fluxes to an LSM model. For

    // Radiative fluxes for LSM
    if (solverChoice.lsm_type != LandSurfaceType::None &&
        solverChoice.rad_type != RadiationType::None)
    {
        BoxList m_bl = ba.boxList();
        for (auto& b : m_bl) {
            int kmin = b.smallEnd(2);
            b.setRange(2,kmin);
        }
        BoxArray m_ba(std::move(m_bl));

        sw_lw_fluxes[lev] = std::make_unique<MultiFab>(m_ba, dm, 6, 0); // DIR/DIF VIS/NIR (4), NET SW (1), LW (1)
        solar_zenith[lev] = std::make_unique<MultiFab>(m_ba, dm, 1, 0);

        sw_lw_fluxes[lev]->setVal(0.);
        solar_zenith[lev]->setVal(0.);
    }

    //*********************************************************
    // Turbulent perturbation region initialization
    //*********************************************************
    if (solverChoice.pert_type == PerturbationType::Source ||
        solverChoice.pert_type == PerturbationType::Direct ||
        solverChoice.pert_type == PerturbationType::CPM)
    {
        amrex::Box bnd_bx = ba.minimalBox();
        turbPert.init_tpi_type(solverChoice.pert_type);
        turbPert.init_tpi(lev, bnd_bx.smallEnd(), bnd_bx.bigEnd(), geom[lev].CellSizeArray(),
                          ba, dm, ngrow_state, pp_prefix, refRatio(), max_level);
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

    land_type_lev[lev].resize(1);
    land_type_lev[lev][0] = std::make_unique<iMultiFab>(ba2d_mask,dm,1,ngv);
    land_type_lev[lev][0]->setVal(0);
    land_type_lev[lev][0]->FillBoundary(geom[lev].periodicity());

    soil_type_lev[lev].resize(1);
    soil_type_lev[lev][0] = std::make_unique<iMultiFab>(ba2d_mask,dm,1,ngv);
    soil_type_lev[lev][0]->setVal(0);
    soil_type_lev[lev][0]->FillBoundary(geom[lev].periodicity());

    urb_frac_lev[lev].resize(1);
    urb_frac_lev[lev][0] = std::make_unique<MultiFab>(ba2d_mask,dm,1,ngv);
    urb_frac_lev[lev][0]->setVal(1.0);
    urb_frac_lev[lev][0]->FillBoundary(geom[lev].periodicity());
    }

    // Read in tables needed for windfarm simulations
    // fill in Nturb multifab - number of turbines in each mesh cell
    // write out the vtk files for wind turbine location and/or
    // actuator disks
    #ifdef ERF_USE_WINDFARM
        //init_windfarm(lev);
    #endif

    if (lev > 0) {
        fine_mask[lev] = std::make_unique<MultiFab>(grids[lev-1], dmap[lev-1], 1, 0);
        build_fine_mask(lev, *fine_mask[lev].get());
    }
}

void
ERF::update_diffusive_arrays (int lev, const BoxArray& ba, const DistributionMapping& dm)
{
    // ********************************************************************************************
    // Diffusive terms
    // ********************************************************************************************
    bool l_use_terrain = (SolverChoice::terrain_type != TerrainType::None);
    bool l_use_kturb   = solverChoice.turbChoice[lev].use_kturb;
    bool l_use_diff    = ( (solverChoice.diffChoice.molec_diff_type != MolecDiffType::None) ||
                           l_use_kturb );
    bool l_need_SmnSmn = solverChoice.turbChoice[lev].use_keqn;
    bool l_use_moist   = (  solverChoice.moisture_type != MoistureType::None  );
    bool l_rotate      = (  solverChoice.use_rotate_surface_flux  );

    bool l_implicit_diff = (solverChoice.vert_implicit_fac[0] > 0 ||
                            solverChoice.vert_implicit_fac[1] > 0 ||
                            solverChoice.vert_implicit_fac[2] > 0);

    BoxArray ba12 = convert(ba, IntVect(1,1,0));
    BoxArray ba13 = convert(ba, IntVect(1,0,1));
    BoxArray ba23 = convert(ba, IntVect(0,1,1));

    Tau[lev].resize(9);
    Tau_corr[lev].resize(3);

    if (l_use_diff) {
        //
        // NOTE: We require ghost cells in the vertical when allowing grids that don't
        //       cover the entire vertical extent of the domain at this level
        //
        for (int i = 0; i < 3; i++) {
            Tau[lev][i] = std::make_unique<MultiFab>( ba  , dm, 1, IntVect(1,1,1) );
        }
        Tau[lev][TauType::tau12] = std::make_unique<MultiFab>( ba12, dm, 1, IntVect(1,1,1) );
        Tau[lev][TauType::tau13] = std::make_unique<MultiFab>( ba13, dm, 1, IntVect(1,1,1) );
        Tau[lev][TauType::tau23] = std::make_unique<MultiFab>( ba23, dm, 1, IntVect(1,1,1) );
        Tau[lev][TauType::tau12]->setVal(0.);
        Tau[lev][TauType::tau13]->setVal(0.);
        Tau[lev][TauType::tau23]->setVal(0.);
        if (l_use_terrain) {
            Tau[lev][TauType::tau21] = std::make_unique<MultiFab>( ba12, dm, 1, IntVect(1,1,1) );
            Tau[lev][TauType::tau31] = std::make_unique<MultiFab>( ba13, dm, 1, IntVect(1,1,1) );
            Tau[lev][TauType::tau32] = std::make_unique<MultiFab>( ba23, dm, 1, IntVect(1,1,1) );
            Tau[lev][TauType::tau21]->setVal(0.);
            Tau[lev][TauType::tau31]->setVal(0.);
            Tau[lev][TauType::tau32]->setVal(0.);
        } else if (l_implicit_diff) {
            Tau[lev][TauType::tau31] = std::make_unique<MultiFab>( ba13, dm, 1, IntVect(1,1,1) );
            Tau[lev][TauType::tau32] = std::make_unique<MultiFab>( ba23, dm, 1, IntVect(1,1,1) );
            Tau[lev][TauType::tau31]->setVal(0.);
            Tau[lev][TauType::tau32]->setVal(0.);
        } else {
            Tau[lev][TauType::tau21] = nullptr;
            Tau[lev][TauType::tau31] = nullptr;
            Tau[lev][TauType::tau32] = nullptr;
        }

        if (l_implicit_diff && solverChoice.implicit_momentum_diffusion)
        {
            Tau_corr[lev][0] = std::make_unique<MultiFab>( ba13, dm, 1, IntVect(1,1,1) ); // Tau31
            Tau_corr[lev][1] = std::make_unique<MultiFab>( ba23, dm, 1, IntVect(1,1,1) ); // Tau32
            Tau_corr[lev][0]->setVal(0.);
            Tau_corr[lev][1]->setVal(0.);
#ifdef ERF_IMPLICIT_W
            Tau_corr[lev][2] = std::make_unique<MultiFab>( ba  , dm, 1, IntVect(1,1,1) ); // Tau33
            Tau_corr[lev][2]->setVal(0.);
#else
            Tau_corr[lev][2] = nullptr;
#endif
        } else {
            Tau_corr[lev][0] = nullptr;
            Tau_corr[lev][1] = nullptr;
            Tau_corr[lev][2] = nullptr;
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
            if (l_rotate) {
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
        for (int i = 0; i < 9; i++) {
            Tau[lev][i] = nullptr;
        }
        SFS_hfx1_lev[lev] = nullptr; SFS_hfx2_lev[lev] = nullptr; SFS_hfx3_lev[lev] = nullptr;
        SFS_diss_lev[lev] = nullptr;
    }

    if (l_use_kturb) {
        eddyDiffs_lev[lev] = std::make_unique<MultiFab>(ba, dm, EddyDiff::NumDiffs, 2);
        eddyDiffs_lev[lev]->setVal(0.0);
        if(l_need_SmnSmn) {
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
    if (solverChoice.init_type != InitType::WRFInput && solverChoice.init_type != InitType::Metgrid)
    {
        if (lev > 0) {
            //
            // First interpolate from coarser level if there is one
            // NOTE: this interpolater assumes that ALL ghost cells of the coarse MultiFab
            //       have been pre-filled - this includes ghost cells both inside and outside
            //       the domain
            //
            InterpFromCoarseLevel(*z_phys_nd[lev], z_phys_nd[lev]->nGrowVect(),
                                  IntVect(0,0,0), // do NOT fill ghost cells outside the domain
                                  *z_phys_nd[lev-1], 0, 0, 1,
                                  geom[lev-1], geom[lev],
                                  refRatio(lev-1), &node_bilinear_interp,
                                  domain_bcs_type, BCVars::cons_bc);
        }

        int ngrow = ComputeGhostCells(solverChoice) + 2;
        Box bx(surroundingNodes(Geom(lev).Domain())); bx.grow(ngrow);
        FArrayBox terrain_fab(makeSlab(bx,2,0),1);

        //
        // If we are using fitted mesh then we use the surface as defined above
        // If we are not using fitted mesh but are using z_levels, we still need z_phys (for now)
        //    but we need to use a flat terrain for the mesh itself (the EB data has already been made
        //    from the correct terrain)
        //
        if (solverChoice.terrain_type != TerrainType::StaticFittedMesh &&
            solverChoice.terrain_type != TerrainType::MovingFittedMesh) {
                terrain_fab.template setVal<RunOn::Device>(0.0);
        } else {
            //
            // Fill the values of the terrain height at k=0 only
            //
            prob->init_terrain_surface(geom[lev],terrain_fab,time);
        }

        for (MFIter mfi(*z_phys_nd[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box isect = terrain_fab.box() & (*z_phys_nd[lev])[mfi].box();
            if (!isect.isEmpty()) {
                (*z_phys_nd[lev])[mfi].template copy<RunOn::Device>(terrain_fab,isect,0,isect,0,1);
            }
        }

        make_terrain_fitted_coords(lev,geom[lev],*z_phys_nd[lev],zlevels_stag[lev],phys_bc_type);

        z_phys_nd[lev]->FillBoundary(geom[lev].periodicity());

        if (lev == 0) {
            Real zmax = z_phys_nd[0]->max(0,0,false);
            Real rel_diff = (zmax - zlevels_stag[0][zlevels_stag[0].size()-1]) / zmax;
            if (rel_diff < 1.e-8) {
                amrex::Print() << "max of zphys_nd " << zmax << std::endl;
                amrex::Print() << "max of zlevels  " << zlevels_stag[0][zlevels_stag[0].size()-1] << std::endl;
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(rel_diff < 1.e-8, "Terrain is taller than domain top!");
            }
        } // lev == 0

    } else {
        // NOTE: If a WRFInput file is NOT provided for a finer level,
        //       we simply interpolate from the coarse. This is necessary
        //       since we average_down the terrain (ERF_MakeNewLevel.cpp L351).
        //       If a WRFInput file IS present, it overwrites the terrain data.
        if (lev > 0) {
            //
            // First interpolate from coarser level if there is one
            // NOTE: this interpolater assumes that ALL ghost cells of the coarse MultiFab
            //       have been pre-filled - this includes ghost cells both inside and outside
            //       the domain
            //
            InterpFromCoarseLevel(*z_phys_nd[lev], z_phys_nd[lev]->nGrowVect(),
                                  z_phys_nd[lev]->nGrowVect(), // DO fill ghost cells outside the domain
                                  *z_phys_nd[lev-1], 0, 0, 1,
                                  geom[lev-1], geom[lev],
                                  refRatio(lev-1), &node_bilinear_interp,
                                  domain_bcs_type, BCVars::cons_bc);
        }
    } // init_type

    if (solverChoice.terrain_type == TerrainType::ImmersedForcing ||
        solverChoice.buildings_type == BuildingsType::ImmersedForcing) {
        terrain_blanking[lev]->setVal(1.0);
        MultiFab::Subtract(*terrain_blanking[lev], EBFactory(lev).getVolFrac(), 0, 0, 1, ComputeGhostCells(solverChoice) + 2);
        terrain_blanking[lev]->FillBoundary(geom[lev].periodicity());
        init_immersed_forcing(lev); // needed for real cases
    }

    // Compute the min dz and pass to the micro model
    Real dzmin = get_dzmin_terrain(*z_phys_nd[lev]);
    micro->Set_dzmin(lev, dzmin);
}

void
ERF::remake_zphys (int lev, Real /*time*/, std::unique_ptr<MultiFab>& temp_zphys_nd)
{
    if (solverChoice.init_type != InitType::WRFInput && solverChoice.init_type != InitType::Metgrid)
    {
        if (lev > 0)
        {
            //
            // First interpolate from coarser level
            // NOTE: this interpolater assumes that ALL ghost cells of the coarse MultiFab
            //       have been pre-filled - this includes ghost cells both inside and outside
            //       the domain
            //
            InterpFromCoarseLevel(*temp_zphys_nd, z_phys_nd[lev]->nGrowVect(),
                                  IntVect(0,0,0), // do NOT fill ghost cells outside the domain
                                  *z_phys_nd[lev-1], 0, 0, 1,
                                  geom[lev-1], geom[lev],
                                  refRatio(lev-1), &node_bilinear_interp,
                                  domain_bcs_type, BCVars::cons_bc);

            // This recomputes the fine values using the bottom terrain at the fine resolution,
            //    and also fills values of z_phys_nd outside the domain
            make_terrain_fitted_coords(lev,geom[lev],*temp_zphys_nd,zlevels_stag[lev],phys_bc_type);

            std::swap(temp_zphys_nd, z_phys_nd[lev]);
        } // lev > 0
    } else {
        if (lev > 0)
        {
            //
            // First interpolate from coarser level
            // NOTE: this interpolater assumes that ALL ghost cells of the coarse MultiFab
            //       have been pre-filled - this includes ghost cells both inside and outside
            //       the domain
            //
            InterpFromCoarseLevel(*temp_zphys_nd, z_phys_nd[lev]->nGrowVect(),
                                  z_phys_nd[lev]->nGrowVect(), // DO fill ghost cells outside the domain
                                  *z_phys_nd[lev-1], 0, 0, 1,
                                  geom[lev-1], geom[lev],
                                  refRatio(lev-1), &node_bilinear_interp,
                                  domain_bcs_type, BCVars::cons_bc);

            std::swap(temp_zphys_nd, z_phys_nd[lev]);
        } // lev > 0
    }

    if (solverChoice.terrain_type == TerrainType::ImmersedForcing ||
        solverChoice.buildings_type == BuildingsType::ImmersedForcing) {
        //
        // This assumes we have already remade the EBGeometry
        //
        terrain_blanking[lev]->setVal(1.0);
        MultiFab::Subtract(*terrain_blanking[lev], EBFactory(lev).getVolFrac(), 0, 0, 1, z_phys_nd[lev]->nGrowVect());
    }

    // Compute the min dz and pass to the micro model
    Real dzmin = get_dzmin_terrain(*z_phys_nd[lev]);
    micro->Set_dzmin(lev, dzmin);
}

void
ERF::update_terrain_arrays (int lev)
{
    if (SolverChoice::mesh_type == MeshType::StretchedDz ||
        SolverChoice::mesh_type == MeshType::VariableDz) {
        make_J(geom[lev],*z_phys_nd[lev],*detJ_cc[lev]);
        make_areas(geom[lev],*z_phys_nd[lev],*ax[lev],*ay[lev],*az[lev]);
        make_zcc(geom[lev],*z_phys_nd[lev],*z_phys_cc[lev]);
    } else { // MeshType::ConstantDz
        if (SolverChoice::terrain_type == TerrainType::EB) {
            const auto& ebfact = *eb[lev]->get_const_factory();
            const MultiFab& volfrac = ebfact.getVolFrac();
            detJ_cc[lev] = std::make_unique<MultiFab>(volfrac, amrex::make_alias, 0, volfrac.nComp());
        }
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
    mri_integrator_mem[lev]->setNoSubstepping((solverChoice.substepping_type[lev] == SubsteppingType::None));
    mri_integrator_mem[lev]->setAnelastic(solverChoice.anelastic[lev]);
    mri_integrator_mem[lev]->setNcompCons(ncomp_cons);
    mri_integrator_mem[lev]->setForceFirstStageSingleSubstep(solverChoice.force_stage1_single_substep);
}

void
ERF::make_physbcs (int lev)
{
    if (SolverChoice::mesh_type == MeshType::VariableDz) {
        AMREX_ALWAYS_ASSERT(z_phys_nd[lev] != nullptr);
    }

    physbcs_cons[lev] = std::make_unique<ERFPhysBCFunct_cons> (lev, geom[lev], domain_bcs_type, domain_bcs_type_d,
                                                               m_bc_extdir_vals, m_bc_neumann_vals,
                                                               z_phys_nd[lev], solverChoice.use_real_bcs, th_bc_data[lev].data());
    physbcs_u[lev]    = std::make_unique<ERFPhysBCFunct_u> (lev, geom[lev], domain_bcs_type, domain_bcs_type_d,
                                                            m_bc_extdir_vals, m_bc_neumann_vals,
                                                            z_phys_nd[lev], solverChoice.use_real_bcs, xvel_bc_data[lev].data());
    physbcs_v[lev]    = std::make_unique<ERFPhysBCFunct_v> (lev, geom[lev], domain_bcs_type, domain_bcs_type_d,
                                                            m_bc_extdir_vals, m_bc_neumann_vals,
                                                            z_phys_nd[lev], solverChoice.use_real_bcs, yvel_bc_data[lev].data());
    physbcs_w[lev]    = std::make_unique<ERFPhysBCFunct_w> (lev, geom[lev], domain_bcs_type, domain_bcs_type_d,
                                                            m_bc_extdir_vals, m_bc_neumann_vals,
                                                            solverChoice.terrain_type, mapfac[lev], z_phys_nd[lev],
                                                            solverChoice.use_real_bcs, zvel_bc_data[lev].data());
    physbcs_base[lev] = std::make_unique<ERFPhysBCFunct_base> (lev, geom[lev], domain_bcs_type, domain_bcs_type_d, z_phys_nd[lev],
                                                               (solverChoice.terrain_type == TerrainType::MovingFittedMesh));
}
