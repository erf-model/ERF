#include <ERF.H>
#include <ERF_Utils.H>

#ifdef ERF_USE_WINDFARM
#include <ERF_WindFarm.H>
#endif

using namespace amrex;

/**
 * Function that advances the solution at one level for a single time step --
 * this does some preliminaries then calls erf_advance
 *
 * @param[in] lev level of refinement (coarsest level is 0)
 * @param[in] time start time for time advance
 * @param[in] dt_lev time step for this time advance
 */

void
ERF::Advance (int lev, Real time, Real dt_lev, int iteration, int /*ncycle*/)
{
    BL_PROFILE("ERF::Advance()");

    // We must swap the pointers so the previous step's "new" is now this step's "old"
    std::swap(vars_old[lev], vars_new[lev]);

    MultiFab& S_old = vars_old[lev][Vars::cons];
    MultiFab& S_new = vars_new[lev][Vars::cons];

    MultiFab& U_old = vars_old[lev][Vars::xvel];
    MultiFab& V_old = vars_old[lev][Vars::yvel];
    MultiFab& W_old = vars_old[lev][Vars::zvel];

    MultiFab& U_new = vars_new[lev][Vars::xvel];
    MultiFab& V_new = vars_new[lev][Vars::yvel];
    MultiFab& W_new = vars_new[lev][Vars::zvel];

    // We need to set these because otherwise in the first call to erf_advance we may
    //    read uninitialized data on ghost values in setting the bc's on the velocities
    U_new.setVal(1.e34,U_new.nGrowVect());
    V_new.setVal(1.e34,V_new.nGrowVect());
    W_new.setVal(1.e34,W_new.nGrowVect());

    //
    // NOTE: the momenta here are not fillpatched (they are only used as scratch space)
    // If lev == 0 we have already FillPatched this in ERF::TimeStep
    //
    if (lev > 0) {
        FillPatchFineLevel(lev, time, {&S_old, &U_old, &V_old, &W_old},
                           {&S_old, &rU_old[lev], &rV_old[lev], &rW_old[lev]},
                           base_state[lev], base_state[lev]);
    }

    //
    // So we must convert the fillpatched to momenta, including the ghost values
    //
    const MultiFab* c_vfrac = nullptr;
    if (solverChoice.terrain_type == TerrainType::EB) {
        c_vfrac = &((get_eb(lev).get_const_factory())->getVolFrac());
    }

    VelocityToMomentum(U_old, rU_old[lev].nGrowVect(),
                       V_old, rV_old[lev].nGrowVect(),
                       W_old, rW_old[lev].nGrowVect(),
                       S_old, rU_old[lev], rV_old[lev], rW_old[lev],
                       Geom(lev).Domain(),
                       domain_bcs_type, c_vfrac);

    // Update the inflow perturbation update time and amplitude
    if (solverChoice.pert_type == PerturbationType::Source ||
        solverChoice.pert_type == PerturbationType::Direct ||
        solverChoice.pert_type == PerturbationType::CPM)
    {
        turbPert.calc_tpi_update(lev, dt_lev, U_old, V_old, S_old);
    }

    // If PerturbationType::Direct or CPM is selected, directly add the computed perturbation
    // on the conserved field
    if (solverChoice.pert_type == PerturbationType::Direct ||
        solverChoice.pert_type == PerturbationType::CPM)
    {
        auto m_ixtype = S_old.boxArray().ixType(); // Conserved term
        for (MFIter mfi(S_old,TileNoZ()); mfi.isValid(); ++mfi) {
            Box bx  = mfi.tilebox();
            const Array4<Real> &cell_data  = S_old.array(mfi);
            const Array4<const Real> &pert_cell = turbPert.pb_cell[lev].array(mfi);
            turbPert.apply_tpi(lev, bx, RhoTheta_comp, m_ixtype, cell_data, pert_cell);
        }
    }

    // configure SurfaceLayer params if needed
    if (phys_bc_type[Orientation(Direction::z,Orientation::low)] == ERF_BC::surface_layer) {
        if (m_SurfaceLayer) {
            IntVect ng = Theta_prim[lev]->nGrowVect();
            MultiFab::Copy(  *Theta_prim[lev], S_old, RhoTheta_comp, 0, 1, ng);
            MultiFab::Divide(*Theta_prim[lev], S_old, Rho_comp     , 0, 1, ng);
            if (solverChoice.moisture_type != MoistureType::None) {
                ng = Qv_prim[lev]->nGrowVect();

                MultiFab::Copy(  *Qv_prim[lev], S_old, RhoQ1_comp, 0, 1, ng);
                MultiFab::Divide(*Qv_prim[lev], S_old, Rho_comp  , 0, 1, ng);

                if (solverChoice.moisture_indices.qr > -1) {
                    MultiFab::Copy(  *Qr_prim[lev], S_old, solverChoice.moisture_indices.qr, 0, 1, ng);
                    MultiFab::Divide(*Qr_prim[lev], S_old, Rho_comp  , 0, 1, ng);
                } else {
                    Qr_prim[lev]->setVal(0.0);
                }
            }
            // NOTE: std::swap above causes the field ptrs to be out of date.
            //       Reassign the field ptrs for MAC avg computation.
            m_SurfaceLayer->update_mac_ptrs(lev, vars_old, Theta_prim, Qv_prim, Qr_prim);
            m_SurfaceLayer->update_pblh(lev, vars_old, z_phys_cc[lev].get(),
                                        solverChoice.moisture_indices);
            m_SurfaceLayer->update_fluxes(lev, time, S_old, z_phys_nd[lev]);
        }
    }

#if defined(ERF_USE_WINDFARM)
    // **************************************************************************************
    // Update the windfarm sources
    // **************************************************************************************
    if (solverChoice.windfarm_type != WindFarmType::None) {
        advance_windfarm(Geom(lev), dt_lev, S_old,
                         U_old, V_old, W_old, vars_windfarm[lev],
                         Nturb[lev], SMark[lev], time);
    }

#endif

    // **************************************************************************************
    // Update the radiation sources with the "old" state
    // **************************************************************************************
    advance_radiation(lev, S_old, dt_lev);

#ifdef ERF_USE_SHOC
    // **************************************************************************************
    // Update the "old" state using SHOC
    // **************************************************************************************
    if (solverChoice.use_shoc) {
        // Get SFC fluxes from SurfaceLayer
        if (m_SurfaceLayer) {
            Vector<const MultiFab*> mfs = {&S_old, &U_old, &V_old, &W_old};
            m_SurfaceLayer->impose_SurfaceLayer_bcs(lev, mfs, Tau[lev],
                                                    SFS_hfx1_lev[lev].get() , SFS_hfx2_lev[lev].get() , SFS_hfx3_lev[lev].get(),
                                                    SFS_q1fx1_lev[lev].get(), SFS_q1fx2_lev[lev].get(), SFS_q1fx3_lev[lev].get(),
                                                    z_phys_nd[lev].get());
        }

        // Get Shoc tendencies and update the state
        Real* w_sub = (solverChoice.custom_w_subsidence) ? d_w_subsid[lev].data() : nullptr;
        compute_shoc_tendencies(lev, &S_old, &U_old, &V_old, &W_old, w_sub,
                                Tau[lev][TauType::tau13].get(), Tau[lev][TauType::tau23].get(),
                                SFS_hfx3_lev[lev].get()       , SFS_q1fx3_lev[lev].get()      ,
                                eddyDiffs_lev[lev].get()      , z_phys_nd[lev].get()          ,
                                dt_lev);
    }
#endif

    const BoxArray&            ba = S_old.boxArray();
    const DistributionMapping& dm = S_old.DistributionMap();

    int nvars = S_old.nComp();

    // Source array for conserved cell-centered quantities -- this will be filled
    //     in the call to make_sources in ERF_TI_slow_rhs_pre.H
    MultiFab cc_source(ba,dm,nvars,1); cc_source.setVal(0.0);

    // Source arrays for momenta -- these will be filled
    //     in the call to make_mom_sources in ERF_TI_slow_rhs_pre.H
    BoxArray ba_x(ba); ba_x.surroundingNodes(0);
    MultiFab xmom_source(ba_x,dm,1,1); xmom_source.setVal(0.0);

    BoxArray ba_y(ba); ba_y.surroundingNodes(1);
    MultiFab ymom_source(ba_y,dm,1,1); ymom_source.setVal(0.0);

    BoxArray ba_z(ba); ba_z.surroundingNodes(2);
    MultiFab zmom_source(ba_z,dm,1,1); zmom_source.setVal(0.0);
    MultiFab    buoyancy(ba_z,dm,1,1); buoyancy.setVal(0.0);

    amrex::Vector<MultiFab> state_old;
    amrex::Vector<MultiFab> state_new;

    // **************************************************************************************
    // Here we define state_old and state_new which are to be advanced
    // **************************************************************************************
    // Initial solution
    // Note that "old" and "new" here are relative to each RK stage.
    state_old.push_back(MultiFab(S_old      , amrex::make_alias, 0, nvars)); // cons
    state_old.push_back(MultiFab(rU_old[lev], amrex::make_alias, 0,     1)); // xmom
    state_old.push_back(MultiFab(rV_old[lev], amrex::make_alias, 0,     1)); // ymom
    state_old.push_back(MultiFab(rW_old[lev], amrex::make_alias, 0,     1)); // zmom

    // Final solution
    // state_new at the end of the last RK stage holds the t^{n+1} data
    state_new.push_back(MultiFab(S_new      , amrex::make_alias, 0, nvars)); // cons
    state_new.push_back(MultiFab(rU_new[lev], amrex::make_alias, 0,     1)); // xmom
    state_new.push_back(MultiFab(rV_new[lev], amrex::make_alias, 0,     1)); // ymom
    state_new.push_back(MultiFab(rW_new[lev], amrex::make_alias, 0,     1)); // zmom

    // **************************************************************************************
    // Tests on the reasonableness of the solution before the dycore
    // **************************************************************************************
    // Test for NaNs after dycore
    if (check_for_nans > 1) {
        if (verbose > 1) {
            amrex::Print() << "Testing old state and vels for NaNs before dycore" << std::endl;
        }
        check_state_for_nans(S_old);
        check_vels_for_nans(rU_old[lev],rV_old[lev],rW_old[lev]);
    }

    // We only test on low temp if we have a moisture model because we are protecting against
    //    the test on low temp inside the moisture models
    if (solverChoice.moisture_type != MoistureType::None) {
        if (verbose > 1) {
            amrex::Print() << "Testing on low temperature before dycore" << std::endl;
        }
        check_for_low_temp(S_old);
    } else {
        if (verbose > 1) {
            amrex::Print() << "Testing on negative temperature before dycore" << std::endl;
        }
        check_for_negative_theta(S_old);
    }

    // **************************************************************************************
    // Update the dycore
    // **************************************************************************************
    advance_dycore(lev, state_old, state_new,
                   U_old, V_old, W_old,
                   U_new, V_new, W_new,
                   cc_source, xmom_source, ymom_source, zmom_source, buoyancy,
                   Geom(lev), dt_lev, time);

    // **************************************************************************************
    // Tests on the reasonableness of the solution after the dycore
    // **************************************************************************************
    // Test for NaNs after dycore
    if (check_for_nans > 0) {
        if (verbose > 1) {
            amrex::Print() << "Testing new state and vels for NaNs after dycore" << std::endl;
        }
        check_state_for_nans(S_new);
        check_vels_for_nans(rU_new[lev],rV_new[lev],rW_new[lev]);
    }

    // We only test on low temp if we have a moisture model because we are protecting against
    //    the test on low temp inside the moisture models
    if (solverChoice.moisture_type != MoistureType::None) {
        if (verbose > 1) {
            amrex::Print() << "Testing on low temperature after dycore" << std::endl;
        }
        check_for_low_temp(S_new);
    } else {
        // Otherwise we will test on negative (rhotheta) coming out of the dycore
        if (verbose > 1) {
            amrex::Print() << "Testing on negative temperature after dycore" << std::endl;
        }
        check_for_negative_theta(S_new);
    }

    // **************************************************************************************
    // Update the microphysics (moisture)
    // **************************************************************************************
    if (!solverChoice.moisture_tight_coupling)
    {
        advance_microphysics(lev, S_new, dt_lev, iteration, time);

        // Test for NaNs after microphysics
        if (check_for_nans > 0) {
            amrex::Print() << "Testing new state for NaNs after advance_microphysics" << std::endl;
            check_state_for_nans(S_new);
        }
    }

    // **************************************************************************************
    // Update the land surface model
    // **************************************************************************************
    advance_lsm(lev, S_new, U_new, V_new, dt_lev);

#ifdef ERF_USE_PARTICLES
    // **************************************************************************************
    // Update the particle positions
    // **************************************************************************************
   evolveTracers( lev, dt_lev, vars_new, z_phys_nd );
#endif

    // ***********************************************************************************************
    // Impose domain boundary conditions here so that in FillPatching the fine data we won't
    // need to re-fill these
    // ***********************************************************************************************
    if (lev < finest_level) {
         IntVect ngvect_vels = vars_new[lev][Vars::xvel].nGrowVect();
         (*physbcs_cons[lev])(vars_new[lev][Vars::cons], vars_new[lev][Vars::xvel], vars_new[lev][Vars::yvel],
                              0,vars_new[lev][Vars::cons].nComp(),
                              vars_new[lev][Vars::cons].nGrowVect(),time,BCVars::cons_bc,true);
            (*physbcs_u[lev])(vars_new[lev][Vars::xvel], vars_new[lev][Vars::xvel], vars_new[lev][Vars::yvel],
                              ngvect_vels,time,BCVars::xvel_bc,true);
            (*physbcs_v[lev])(vars_new[lev][Vars::yvel], vars_new[lev][Vars::xvel], vars_new[lev][Vars::yvel],
                              ngvect_vels,time,BCVars::yvel_bc,true);
            (*physbcs_w[lev])(vars_new[lev][Vars::zvel], vars_new[lev][Vars::xvel], vars_new[lev][Vars::yvel],
                              ngvect_vels,time,BCVars::zvel_bc,true);
    }

    // **************************************************************************************
    // Register old and new coarse data if we are at a level less than the finest level
    // **************************************************************************************
    if (lev < finest_level) {
        if (cf_width > 0) {
            // We must fill the ghost cells of these so that the parallel copy works correctly
            state_old[IntVars::cons].FillBoundary(geom[lev].periodicity());
            state_new[IntVars::cons].FillBoundary(geom[lev].periodicity());
            FPr_c[lev].RegisterCoarseData({&state_old[IntVars::cons], &state_new[IntVars::cons]},
                                          {time, time+dt_lev});
        }

        if (cf_width >= 0) {
            // We must fill the ghost cells of these so that the parallel copy works correctly
            state_old[IntVars::xmom].FillBoundary(geom[lev].periodicity());
            state_new[IntVars::xmom].FillBoundary(geom[lev].periodicity());
            FPr_u[lev].RegisterCoarseData({&state_old[IntVars::xmom], &state_new[IntVars::xmom]},
                                          {time, time+dt_lev});

            state_old[IntVars::ymom].FillBoundary(geom[lev].periodicity());
            state_new[IntVars::ymom].FillBoundary(geom[lev].periodicity());
            FPr_v[lev].RegisterCoarseData({&state_old[IntVars::ymom], &state_new[IntVars::ymom]},
                                          {time, time+dt_lev});

            state_old[IntVars::zmom].FillBoundary(geom[lev].periodicity());
            state_new[IntVars::zmom].FillBoundary(geom[lev].periodicity());
            FPr_w[lev].RegisterCoarseData({&state_old[IntVars::zmom], &state_new[IntVars::zmom]},
                                          {time, time+dt_lev});
        }

            //
            // Now create a MultiFab that holds (S_new - S_old) / dt from the coarse level interpolated
            //     on to the coarse/fine boundary at the fine resolution
            //
            Interpolater* mapper_f = &face_cons_linear_interp;

            // PhysBCFunctNoOp null_bc;
            // MultiFab tempx(vars_new[lev+1][Vars::xvel].boxArray(),vars_new[lev+1][Vars::xvel].DistributionMap(),1,0);
            // tempx.setVal(0.0);
            // xmom_crse_rhs[lev+1].setVal(0.0);
            // FPr_u[lev].FillSet(tempx               , time       , null_bc, domain_bcs_type);
            // FPr_u[lev].FillSet(xmom_crse_rhs[lev+1], time+dt_lev, null_bc, domain_bcs_type);
            // MultiFab::Subtract(xmom_crse_rhs[lev+1],tempx,0,0,1,IntVect{0});
            // xmom_crse_rhs[lev+1].mult(1.0/dt_lev,0,1,0);

            // MultiFab tempy(vars_new[lev+1][Vars::yvel].boxArray(),vars_new[lev+1][Vars::yvel].DistributionMap(),1,0);
            // tempy.setVal(0.0);
            // ymom_crse_rhs[lev+1].setVal(0.0);
            // FPr_v[lev].FillSet(tempy               , time       , null_bc, domain_bcs_type);
            // FPr_v[lev].FillSet(ymom_crse_rhs[lev+1], time+dt_lev, null_bc, domain_bcs_type);
            // MultiFab::Subtract(ymom_crse_rhs[lev+1],tempy,0,0,1,IntVect{0});
            // ymom_crse_rhs[lev+1].mult(1.0/dt_lev,0,1,0);

            MultiFab temp_state(zmom_crse_rhs[lev+1].boxArray(),zmom_crse_rhs[lev+1].DistributionMap(),1,0);
            InterpFromCoarseLevel(temp_state,            IntVect{0}, IntVect{0}, state_old[IntVars::zmom], 0, 0, 1,
                                  geom[lev], geom[lev+1], refRatio(lev), mapper_f, domain_bcs_type, BCVars::zvel_bc);
            InterpFromCoarseLevel(zmom_crse_rhs[lev+1],  IntVect{0}, IntVect{0}, state_new[IntVars::zmom], 0, 0, 1,
                                  geom[lev], geom[lev+1], refRatio(lev), mapper_f, domain_bcs_type, BCVars::zvel_bc);
            MultiFab::Subtract(zmom_crse_rhs[lev+1],temp_state,0,0,1,IntVect{0});
            zmom_crse_rhs[lev+1].mult(1.0/dt_lev,0,1,0);
    }

    // ***********************************************************************************************
    // Update the time averaged velocities if they are requested
    // ***********************************************************************************************
    if (solverChoice.time_avg_vel) {
        Time_Avg_Vel_atCC(dt[lev], t_avg_cnt[lev], vel_t_avg[lev].get(), U_new, V_new, W_new);
    }
}
