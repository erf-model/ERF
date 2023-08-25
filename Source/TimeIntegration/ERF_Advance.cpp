#include <ERF.H>
#include <TileNoZ.H>
#include <Utils.H>

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
        amrex::IntVect ng = S_old.nGrowVect(); ng[2]=0;
        MultiFab::Copy(  *Theta_prim[lev], S_old, Cons::RhoTheta, 0, 1, ng);
        MultiFab::Divide(*Theta_prim[lev], S_old, Cons::Rho     , 0, 1, ng);
        // NOTE: std::swap above causes the field ptrs to be out of date.
        //       Reassign the field ptrs for MAC avg computation.
        m_most->update_mac_ptrs(lev, vars_old, Theta_prim);
        m_most->update_fluxes(lev);
      }
    }

    // We need to set these because otherwise in the first call to erf_advance we may
    //    read uninitialized data on ghost values in setting the bc's on the velocities
    U_new.setVal(1.e34,U_new.nGrowVect());
    V_new.setVal(1.e34,V_new.nGrowVect());
    W_new.setVal(1.e34,W_new.nGrowVect());

    FillPatch(lev, time, {&vars_old[lev][Vars::cons], &vars_old[lev][Vars::xvel],
                          &vars_old[lev][Vars::yvel], &vars_old[lev][Vars::zvel]});
#if defined(ERF_USE_MOISTURE)
    FillPatchMoistVars(lev, qmoist[lev]);
#endif

    MultiFab* S_crse;
    MultiFab rU_crse, rV_crse, rW_crse;

    if (lev > 0)
    {
        S_crse = &vars_old[lev-1][Vars::cons];

        MultiFab& U_crse = vars_old[lev-1][Vars::xvel];
        MultiFab& V_crse = vars_old[lev-1][Vars::yvel];
        MultiFab& W_crse = vars_old[lev-1][Vars::zvel];

        rU_crse.define(U_crse.boxArray(), U_crse.DistributionMap(), 1, U_crse.nGrow());
        rV_crse.define(V_crse.boxArray(), V_crse.DistributionMap(), 1, V_crse.nGrow());
        rW_crse.define(W_crse.boxArray(), W_crse.DistributionMap(), 1, W_crse.nGrow());

        MultiFab density(*S_crse, make_alias, Rho_comp, 1);
        VelocityToMomentum(U_crse, U_crse.nGrowVect(),
                           V_crse, V_crse.nGrowVect(),
                           W_crse, W_crse.nGrowVect(),
                           density ,rU_crse, rV_crse, rW_crse,
                           solverChoice.use_NumDiff);
    }

    // Do an error check
    if (solverChoice.pbl_type == PBLType::MYNN25 &&
        phys_bc_type[Orientation(Direction::z,Orientation::low)] != ERF_BC::MOST) {
        amrex::Error("Must use MOST BC for MYNN2.5 PBL model");
    }

    const auto& local_ref_ratio = (lev > 0) ? ref_ratio[lev-1] : IntVect(1,1,1);

    InterpFaceRegister ifr;
    if (lev > 0)
    {
        ifr.define(S_old.boxArray(), S_old.DistributionMap(), Geom(lev), local_ref_ratio);
    }

    const BoxArray&            ba = S_old.boxArray();
    const DistributionMapping& dm = S_old.DistributionMap();

    int nvars = S_old.nComp();

    // Place-holder for source array -- for now just set to 0
    MultiFab source(ba,dm,nvars,1);
    source.setVal(0.0);
#if defined(ERF_USE_WARM_NO_PRECIP)
    Real tau_cond = solverChoice.tau_cond;
    Real      c_p = solverChoice.c_p;
    condensation_source(source, S_new, tau_cond, c_p);
#endif

    // We don't need to call FillPatch on cons_mf because we have fillpatch'ed S_old above
    MultiFab cons_mf(ba,dm,nvars,S_old.nGrowVect());
    MultiFab::Copy(cons_mf,S_old,0,0,S_old.nComp(),S_old.nGrowVect());

    // Define Multifab for buoyancy term -- only added to vertical velocity
    MultiFab buoyancy(W_old.boxArray(),W_old.DistributionMap(),1,1);

    // Update the dycore
    advance_dycore(lev,
                  cons_mf, S_new,
                  U_old, V_old, W_old,
                  U_new, V_new, W_new,
                  rU_old[lev], rV_old[lev], rW_old[lev],
                  rU_new[lev], rV_new[lev], rW_new[lev],
                  rU_crse, rV_crse, rW_crse,
                  source, buoyancy,
                  Geom(lev), dt_lev, time, &ifr);

#if defined(ERF_USE_MOISTURE)
    // Update the microphysics
    advance_microphysics(lev, S_new, dt_lev);
#endif

#ifdef ERF_USE_PARTICLES
    // Update tracer particles on level 0
    if (lev == 0 && use_tracer_particles) {
        MultiFab* Umac = &vars_new[lev][Vars::xvel];
        tracer_particles->AdvectWithUmac(Umac, lev, dt_lev, *z_phys_nd[0]);
    }
#endif
}
