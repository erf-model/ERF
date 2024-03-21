#include <ERF.H>
#include <TileNoZ.H>
#include <Utils.H>
#ifdef ERF_USE_WINDFARM
#include <Fitch.H>
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
ERF::Advance (int lev, Real time, Real dt_lev, int /*iteration*/, int /*ncycle*/)
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


    // configure ABLMost params if used MostWall boundary condition
    if (phys_bc_type[Orientation(Direction::z,Orientation::low)] == ERF_BC::MOST) {
        if (m_most) {
            IntVect ng = Theta_prim[lev]->nGrowVect();
            MultiFab::Copy(  *Theta_prim[lev], S_old, RhoTheta_comp, 0, 1, ng);
            MultiFab::Divide(*Theta_prim[lev], S_old, Rho_comp     , 0, 1, ng);
            if (solverChoice.moisture_type != MoistureType::None) {
                ng = Qv_prim[lev]->nGrowVect();
                MultiFab::Copy(  *Qv_prim[lev], S_old, RhoQ1_comp, 0, 1, ng);
                MultiFab::Divide(*Qv_prim[lev], S_old, Rho_comp  , 0, 1, ng);
            }
            // NOTE: std::swap above causes the field ptrs to be out of date.
            //       Reassign the field ptrs for MAC avg computation.
            m_most->update_mac_ptrs(lev, vars_old, Theta_prim, Qv_prim);
            m_most->update_fluxes(lev, time);
        }
    }

    // We need to set these because otherwise in the first call to erf_advance we may
    //    read uninitialized data on ghost values in setting the bc's on the velocities
    U_new.setVal(1.e34,U_new.nGrowVect());
    V_new.setVal(1.e34,V_new.nGrowVect());
    W_new.setVal(1.e34,W_new.nGrowVect());

    FillPatch(lev, time, {&S_old, &U_old, &V_old, &W_old},
                         {&S_old, &rU_old[lev], &rV_old[lev], &rW_old[lev]});

    if (solverChoice.moisture_type != MoistureType::None) {
        // TODO: This is only qv
        if (qmoist[lev].size() > 0) FillPatchMoistVars(lev, *(qmoist[lev][0]));
    }

#if defined(ERF_USE_WINDFARM)
    // Update with the Fitch source terms
    if (solverChoice.windfarm_type == WindFarmType::Fitch) {
        fitch_advance(lev, Geom(lev), dt_lev, S_old,
                      U_old, V_old, W_old, vars_fitch[lev], Nturb[lev]);
    }
#endif

    const BoxArray&            ba = S_old.boxArray();
    const DistributionMapping& dm = S_old.DistributionMap();

    int nvars = S_old.nComp();

    // Place-holder for source array -- for now just set to 0
    MultiFab source(ba,dm,nvars,1);
    source.setVal(0.0);

    // We don't need to call FillPatch on cons_mf because we have fillpatch'ed S_old above
    MultiFab cons_mf(ba,dm,nvars,S_old.nGrowVect());
    MultiFab::Copy(cons_mf,S_old,0,0,S_old.nComp(),S_old.nGrowVect());

    // Define Multifab for buoyancy term -- only added to vertical velocity
    MultiFab buoyancy(W_old.boxArray(),W_old.DistributionMap(),1,1);

    amrex::Vector<MultiFab> state_old;
    amrex::Vector<MultiFab> state_new;

    // **************************************************************************************
    // Here we define state_old and state_new which are to be advanced
    // **************************************************************************************
    // Initial solution
    // Note that "old" and "new" here are relative to each RK stage.
    state_old.push_back(MultiFab(cons_mf    , amrex::make_alias, 0, nvars)); // cons
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
    // Update the dycore
    // **************************************************************************************
    advance_dycore(lev, state_old, state_new,
                   U_old, V_old, W_old,
                   U_new, V_new, W_new,
                   source, buoyancy,
                   Geom(lev), dt_lev, time);

    // **************************************************************************************
    // Update the microphysics (moisture)
    // **************************************************************************************
    advance_microphysics(lev, S_new, dt_lev);

    // **************************************************************************************
    // Update the land surface model
    // **************************************************************************************
    advance_lsm(lev, S_new, dt_lev);

#if defined(ERF_USE_RRTMGP)
    // **************************************************************************************
    // Update the radiation
    // **************************************************************************************
    advance_radiation(lev, S_new, dt_lev);
#endif

#ifdef ERF_USE_PARTICLES
    // **************************************************************************************
    // Update the particle positions
    // **************************************************************************************
   evolveTracers( lev, dt_lev, vars_new, z_phys_nd );
#endif

    // **************************************************************************************
    // Register old and new coarse data if we are at a level less than the finest level
    // **************************************************************************************
    if (lev < finest_level)
    {
        if (cf_width > 0) {
            // We must fill the ghost cells of these so that the parallel copy works correctly
            state_old[IntVars::cons].FillBoundary(geom[lev].periodicity());
            state_new[IntVars::cons].FillBoundary(geom[lev].periodicity());
            FPr_c[lev].RegisterCoarseData({&state_old[IntVars::cons], &state_new[IntVars::cons]},
                                          {time, time + dt_lev});
        }

        if (cf_width >= 0) {
            // We must fill the ghost cells of these so that the parallel copy works correctly
            state_old[IntVars::xmom].FillBoundary(geom[lev].periodicity());
            state_new[IntVars::xmom].FillBoundary(geom[lev].periodicity());
            FPr_u[lev].RegisterCoarseData({&state_old[IntVars::xmom], &state_new[IntVars::xmom]},
                                          {time, time + dt_lev});

            state_old[IntVars::ymom].FillBoundary(geom[lev].periodicity());
            state_new[IntVars::ymom].FillBoundary(geom[lev].periodicity());
            FPr_v[lev].RegisterCoarseData({&state_old[IntVars::ymom], &state_new[IntVars::ymom]},
                                          {time, time + dt_lev});

            state_old[IntVars::zmom].FillBoundary(geom[lev].periodicity());
            state_new[IntVars::zmom].FillBoundary(geom[lev].periodicity());
            FPr_w[lev].RegisterCoarseData({&state_old[IntVars::zmom], &state_new[IntVars::zmom]},
                                          {time, time + dt_lev});
        }
    }
}
