#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>

#include <TimeIntegration.H>
#include <utils.H>

#ifdef AMREX_USE_SUNDIALS
#include <arkode/arkode_erkstep.h>     /* prototypes for ERKStep fcts., consts */
#include <arkode/arkode_arkstep.h>     /* prototypes for ARKStep fcts., consts */
#include <arkode/arkode_mristep.h>     /* prototypes for MRIStep fcts., consts */
#include <nvector/nvector_manyvector.h>/* manyvector N_Vector types, fcts. etc */
#include <AMReX_NVector_MultiFab.H>    /* MultiFab N_Vector types, fcts., macros */
#include <AMReX_Sundials.H>    /* MultiFab N_Vector types, fcts., macros */
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h> /* access to SPGMR SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype, etc */
#endif

using namespace amrex;

#ifdef AMREX_USE_SUNDIALS
/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int f0(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int ProcessStage(realtype t, N_Vector y_data, void *user_data);
#endif

// TODO: Check if the order of applying BC on cell-centered state or face-centered mom makes any difference

void erf_advance(int level,
                 MultiFab& cons_old,  MultiFab& cons_new,
                 MultiFab& xvel_old, MultiFab& yvel_old, MultiFab& zvel_old,
                 MultiFab& xvel_new, MultiFab& yvel_new, MultiFab& zvel_new,
                 MultiFab& xmom_crse, MultiFab& ymom_crse, MultiFab& zmom_crse,
                 MultiFab& source,
                 std::array< MultiFab, AMREX_SPACEDIM>& advflux,
                 std::array< MultiFab, AMREX_SPACEDIM>& diffflux,
                 const amrex::Geometry crse_geom,
                 const amrex::Geometry fine_geom,
                 const amrex::IntVect ref_ratio,
                 const amrex::Real dt, const amrex::Real time,
                 amrex::InterpFaceRegister* ifr,
                 const SolverChoice& solverChoice,
                 const amrex::Real* dptr_dens_hse,
                 const amrex::Real* dptr_pres_hse,
                 const amrex::Real* dptr_rayleigh_tau,
                 const amrex::Real* dptr_rayleigh_ubar,
                 const amrex::Real* dptr_rayleigh_vbar,
                 const amrex::Real* dptr_rayleigh_thetabar)
{
    BL_PROFILE_VAR("erf_advance()",erf_advance);

    int nvars = cons_old.nComp();

    // Determine the number of ghost cells depending on the spatial order
    // **************************************************************************************
    //TODO: Check if this is the only place to specify the number of ghost cells
    //TODO: Also explore how 'ngrow' should be related to the spatial_order
    //int ngc = ComputeGhostCells(solverChoice.spatial_order);

    const BoxArray& ba            = cons_old.boxArray();
    const DistributionMapping& dm = cons_old.DistributionMap();

    // Initial momentum -- on faces
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > state_old;
    state_old.emplace_back(std::make_unique<amrex::MultiFab>(ba, dm, nvars, cons_old.nGrow())); // cons
    state_old.emplace_back(std::make_unique<amrex::MultiFab>(convert(ba,IntVect(1,0,0)), dm, 1, 1)); // xmom
    state_old.emplace_back(std::make_unique<amrex::MultiFab>(convert(ba,IntVect(0,1,0)), dm, 1, 1)); // ymom
    state_old.emplace_back(std::make_unique<amrex::MultiFab>(convert(ba,IntVect(0,0,1)), dm, 1, 1)); // zmom

    // Final momentum -- on faces
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > state_new;
    state_new.emplace_back(std::make_unique<amrex::MultiFab>(ba, dm, nvars, cons_old.nGrow())); // cons
    state_new.emplace_back(std::make_unique<amrex::MultiFab>(convert(ba,IntVect(1,0,0)), dm, 1, 1)); // xmom
    state_new.emplace_back(std::make_unique<amrex::MultiFab>(convert(ba,IntVect(0,1,0)), dm, 1, 1)); // ymom
    state_new.emplace_back(std::make_unique<amrex::MultiFab>(convert(ba,IntVect(0,0,1)), dm, 1, 1)); // zmom

    // **************************************************************************************
    // Prepare the old-time data for calling the integrator
    // **************************************************************************************

    // Apply BC on old state data at cells
    // **************************************************************************************
    cons_old.FillBoundary(fine_geom.periodicity());

    // We need to apply the boundary conditions here because we are converting from velocity to momentum
    //    which requires having set boundary conditions on density
    amrex::Vector<MultiFab*> vars_orig{&cons_old};
    ERF::applyBCs(fine_geom, vars_orig);

    MultiFab::Copy(*state_old[IntVar::cons], cons_old, 0, 0, cons_old.nComp(), cons_old.nGrow());

    // Convert old velocity available on faces to old momentum on faces to be used in time integration
    // **************************************************************************************
    VelocityToMomentum(xvel_old, yvel_old, zvel_old,
                       *state_old[IntVar::cons],
                       *state_old[IntVar::xmom],
                       *state_old[IntVar::ymom],
                       *state_old[IntVar::zmom],
                       solverChoice.spatial_order);

    // Apply BC on old momentum data on faces before integration
    // **************************************************************************************
    state_old[IntVar::xmom]->FillBoundary(fine_geom.periodicity());
    state_old[IntVar::ymom]->FillBoundary(fine_geom.periodicity());
    state_old[IntVar::zmom]->FillBoundary(fine_geom.periodicity());

    auto interpolate_coarse_fine_faces = [&](Vector<std::unique_ptr<MultiFab> >& S_data) {
        if (level > 0)
        {
            amrex::Array<const MultiFab*,3> cmf_const{&xmom_crse, &ymom_crse, &zmom_crse};
            amrex::Array<MultiFab*,3> fmf{S_data[IntVar::xmom].get(),
                                          S_data[IntVar::ymom].get(),
                                          S_data[IntVar::zmom].get()};

            // Interpolate from coarse faces to fine faces *only* on the coarse-fine boundary
            ifr->interp(fmf,cmf_const,0,1);

            amrex::Array<MultiFab*,3> cmf{&xmom_crse, &ymom_crse, &zmom_crse};

            int nGrow = 1;
            BoxArray fine_grids(cons_old.boxArray());

            // Interpolate from coarse faces on fine faces outside the fine region
            create_umac_grown(level, nGrow, fine_grids, crse_geom, fine_geom, cmf, fmf, ref_ratio);
        }
    };

    auto apply_bcs = [&](Vector<std::unique_ptr<MultiFab> >& S_data) {
        amrex::Vector<MultiFab*> state_p{S_data[IntVar::cons].get(),
                                         S_data[IntVar::xmom].get(),
                                         S_data[IntVar::ymom].get(),
                                         S_data[IntVar::zmom].get()};

        ERF::applyBCs(fine_geom, state_p);
    };

    // are these in the right order?
    interpolate_coarse_fine_faces(state_old);
    apply_bcs(state_old);

    // **************************************************************************************
    // Setup the integrator
    // **************************************************************************************
    TimeIntegrator<amrex::Vector<std::unique_ptr<amrex::MultiFab> > > integrator(state_old);

#ifdef AMREX_USE_SUNDIALS

    ////STEP ONE
    // Create SUNDIALS specific objects
    SUNNonlinearSolver NLS = NULL;    /* empty nonlinear solver object */
    void *inner_mem = NULL;      /* empty ARKode memory structure */
    void *mristep_mem = NULL;      /* empty ARKode memory structure */
    // Create an N_Vector wrapper for the solution MultiFab
    auto get_length = [&](int index) -> sunindextype {
        auto* p_mf = state_old[index].get();
        return p_mf->nComp() * (p_mf->boxArray()).numPts();
    };

    ////STEP TWO
    sunindextype length = get_length(IntVar::cons);
    sunindextype length_mx = get_length(IntVar::xmom);
    sunindextype length_my = get_length(IntVar::ymom);
    sunindextype length_mz = get_length(IntVar::zmom);
    // Testing data structures, this length may be different for "real" initial condition
    int NVar             = 4;
    //Arbitrary tolerances
    Real reltol          = 1e-4;
    Real abstol          = 1e-4;
    Real t               = time;
    Real tout            = time+dt;
    Real hfixed          = dt;
    Real m               = 10;
    Real hfixed_mri      = dt / m;
    N_Vector nv_cons     = N_VMake_MultiFab(length, state_old[IntVar::cons].get());
    N_Vector nv_xmom     = N_VMake_MultiFab(length_mx, state_old[IntVar::xmom].get());
    N_Vector nv_ymom     = N_VMake_MultiFab(length_my, state_old[IntVar::ymom].get());
    N_Vector nv_zmom     = N_VMake_MultiFab(length_mz, state_old[IntVar::zmom].get());
    N_Vector nv_many_arr[NVar];              /* vector array composed of cons, xmom, ymom, zmom component vectors */

    ////STEP THREE
    /* Create manyvector for solution */
    nv_many_arr[0] = nv_cons;
    nv_many_arr[1] = nv_xmom;
    nv_many_arr[2] = nv_ymom;
    nv_many_arr[3] = nv_zmom;
    N_Vector nv_S = N_VNew_ManyVector(NVar, nv_many_arr);
#endif

    //Create function lambdas
    bool l_lo_z_is_no_slip = ERF::lo_z_is_no_slip;
    bool l_hi_z_is_no_slip = ERF::hi_z_is_no_slip;

    auto rhs_fun = [&](Vector<std::unique_ptr<MultiFab> >& S_rhs, const Vector<std::unique_ptr<MultiFab> >& S_data, const Real time) {
        erf_rhs(level, S_rhs, S_data,
                source,
                advflux, diffflux,
                fine_geom, dt,
                ifr,
                solverChoice,
                l_lo_z_is_no_slip,
                l_hi_z_is_no_slip,
                dptr_dens_hse, dptr_pres_hse,
                dptr_rayleigh_tau, dptr_rayleigh_ubar,
                dptr_rayleigh_vbar, dptr_rayleigh_thetabar);
    };

    auto post_update_fun = [&](Vector<std::unique_ptr<MultiFab> >& S_data, const Real time) {
        // Apply BC on updated state and momentum data
        for (auto& mfp : S_data) {
            mfp->FillBoundary(fine_geom.periodicity());
        }

        // are these in the right order?
        apply_bcs(S_data);

        // TODO: we should interpolate coarse data in time first, so that this interplation
        // in space is at the correct time indicated by the `time` function argument.
        interpolate_coarse_fine_faces(S_data);
    };

    // define rhs and 'post update' utility function that is called after calculating
    // any state data (e.g. at RK stages or at the end of a timestep)
    integrator.set_rhs(rhs_fun);
    integrator.set_post_update(post_update_fun);

#ifdef AMREX_USE_SUNDIALS
    bool use_erk3 = false;
    ////STEP FOUR
    /* Call ARKStepCreate to initialize the inner ARK timestepper module and
    specify the right-hand side function in y'=f(t,y), the inital time
    T0, and the initial dependent variable vector y. */
    if(use_erk3)
      inner_mem = ARKStepCreate(f0, NULL, time, nv_S);
    else
      inner_mem = ARKStepCreate(NULL, f0, time, nv_S);

    ////STEP FIVE
    ARKStepSetFixedStep(inner_mem, hfixed_mri);            // Specify fixed time step size

    ARKodeButcherTable B = ARKodeButcherTable_Alloc(3, SUNFALSE);
    if(use_erk3)
    {
    // Use SSP-RK3
    B->A[1][0] = 1.0;
    B->A[2][0] = 0.25;
    B->A[2][1] = 0.25;
    B->b[0] = 1./6.;
    B->b[1] = 1./6.;
    B->b[2] = 2./3.;
    B->c[1] = 1.0;
    B->c[2] = 0.5;
    B->q=3;
    B->p=0;
    ARKStepSetTables(inner_mem, B->q, B->p, NULL, B);       // Specify Butcher table
    }
    else
    {
    B->A[1][0] = 1.0;
    B->A[2][0] = 0.5;
    B->A[2][2] = 0.5;
    B->b[0] = 0.5;
    B->b[2] = 0.5;
    B->c[1] = 1.0;
    B->c[2] = 1.0;
    B->q=2;
    ARKStepSetTables(inner_mem, B->q, B->p, B, NULL);       // Specify Butcher table
    }
    //Set table
    //    ERKStepSetTable(arkode_mem, B);
    ////STEP SIX
    mristep_mem = MRIStepCreate(f, time, nv_S, MRISTEP_ARKSTEP, inner_mem);
    ////STEP SEVEN
    MRIStepSetFixedStep(mristep_mem, hfixed);
    ////STEP 8.1
    /* Specify tolerances */
    MRIStepSStolerances(mristep_mem, reltol, abstol);
    ////STEP 8.2
    NLS = SUNNonlinSol_FixedPoint(nv_S, 50);
    ////STEP 8.3
    //    ARKStepSetNonlinearSolver(inner_mem, NLS);
    MRIStepSetNonlinearSolver(mristep_mem, NLS);
    ////STEP NINE
    MRIStepSetUserData(mristep_mem, (void *) &integrator);  /* Pass udata to user functions */
    MRIStepSetPostprocessStageFn(mristep_mem, ProcessStage);

    MRIStepSetTable(mristep_mem, B->q, B);

    // Free the Butcher table
    ARKodeButcherTable_Free(B);

    bool advance_erk=false;
    amrex::ParmParse pp("integration.sundials");

    pp.query("erk", advance_erk);

    if(advance_erk)
    {
    ////STEP ELEVEN
    // Use MRIStep to evolve state_old data (wrapped in nv_S) from t to tout=t+dt
    MRIStepEvolve(mristep_mem, tout, nv_S, &t, ARK_NORMAL);

    // Copy the result stored in nv_S to state_new
    for(int i=0; i<N_VGetNumSubvectors_ManyVector(nv_S); i++)
    {
    MultiFab::Copy(*state_new[i], *NV_MFAB(N_VGetSubvector_ManyVector(nv_S, i)), 0, 0, state_new[i]->nComp(), state_new[i]->nGrow());
    }
    }
    else
#endif
    {
    // **************************************************************************************
    // Integrate for a single timestep
    // **************************************************************************************
    integrator.advance(state_old, state_new, time, dt);
    }

#ifdef AMREX_USE_SUNDIALS
  ////STEP THIRTEEN
  N_VDestroy(nv_cons);
  N_VDestroy(nv_xmom);
  N_VDestroy(nv_ymom);
  N_VDestroy(nv_zmom);
  N_VDestroy(nv_S);
  ////STEP FOURTEEN
  MRIStepFree(&mristep_mem);
  ARKStepFree(&inner_mem);
  SUNNonlinSolFree(NLS);    // Free nonlinear solvers
#endif
    // **************************************************************************************
    // Convert updated momentum to updated velocity on faces after we have taken a timestep
    // **************************************************************************************
    MomentumToVelocity(xvel_new, yvel_new, zvel_new,
                       *state_new[IntVar::cons],
                       *state_new[IntVar::xmom],
                       *state_new[IntVar::ymom],
                       *state_new[IntVar::zmom],
                       0, solverChoice.spatial_order);

    // **************************************************************************************
    // Get the final cell centered variables after the step
    // (do this at the very end because its a swap not a copy)
    // **************************************************************************************
    std::swap(cons_new, *state_new[IntVar::cons]);

    // One final application of non-periodic BCs

    xvel_new.FillBoundary(fine_geom.periodicity());
    yvel_new.FillBoundary(fine_geom.periodicity());
    zvel_new.FillBoundary(fine_geom.periodicity());
    amrex::Vector<MultiFab*> vars{&cons_new, &xvel_new, &yvel_new, &zvel_new};

    ERF::applyBCs(fine_geom, vars);

}

#ifdef AMREX_USE_SUNDIALS
// f0 routine to compute a zero-valued ODE RHS function f(t,y).
static int f0(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  // Initialize ydot to zero and return
  N_VConst(0.0, ydot);
  return 0;
}

/* f routine to compute the ODE RHS function f(t,y). */
static int f(realtype t, N_Vector y_data, N_Vector y_rhs, void *user_data)
{
  TimeIntegrator<amrex::Vector<std::unique_ptr<amrex::MultiFab> > > *integrator = (TimeIntegrator<amrex::Vector<std::unique_ptr<amrex::MultiFab> > > *) user_data;
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > S_data;
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > S_rhs;

  const int num_vecs = N_VGetNumSubvectors_ManyVector(y_data);
  S_data.resize(num_vecs);
  S_rhs.resize(num_vecs);

  for(int i=0; i<num_vecs; i++)
  {
      S_data[i].reset(NV_MFAB(N_VGetSubvector_ManyVector(y_data, i)));
      S_rhs[i].reset(NV_MFAB(N_VGetSubvector_ManyVector(y_rhs, i)));
  }

  integrator->call_post_update(S_data, t);
  integrator->rhs(S_rhs, S_data, t);

  for(int i=0; i<num_vecs; i++)
  {
      S_data[i].release();
      S_rhs[i].release();
  }

  return 0;
}

static int ProcessStage(realtype t, N_Vector y_data, void *user_data)
{
  TimeIntegrator<amrex::Vector<std::unique_ptr<amrex::MultiFab> > > *integrator = (TimeIntegrator<amrex::Vector<std::unique_ptr<amrex::MultiFab> > > *) user_data;
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > S_data;

  const int num_vecs = N_VGetNumSubvectors_ManyVector(y_data);
  S_data.resize(num_vecs);

  for(int i=0; i<num_vecs; i++)
  {
      S_data[i].reset(NV_MFAB(N_VGetSubvector_ManyVector(y_data, i)));
  }

  integrator->call_post_update(S_data, t);

  for(int i=0; i<num_vecs; i++)
  {
      S_data[i].release();
  }

  return 0;
}
#endif
