#include <AMReX_BC_TYPES.H>
#include <SpatialStencils.H>
#include <prob_common.H>
#include <AMReX_TimeIntegrator.H>
#include <ERF_MRI.H>
#include <ERF_SRI.H>
#include <TerrainMetrics.H>
#include <TimeIntegration.H>
#include <ERF.H>

using namespace amrex;

void ERF::erf_advance(int level,
                      MultiFab& cons_old,  MultiFab& cons_new,
                      MultiFab& xvel_old, MultiFab& yvel_old, MultiFab& zvel_old,
                      MultiFab& xvel_new, MultiFab& yvel_new, MultiFab& zvel_new,
                      MultiFab& xmom_crse, MultiFab& ymom_crse, MultiFab& zmom_crse,
                      MultiFab& source,
                      std::array< MultiFab, AMREX_SPACEDIM>& flux,
                      const amrex::Geometry fine_geom,
                      const amrex::Real dt_advance, const amrex::Real old_time,
                      amrex::InterpFaceRegister* ifr,
                      MultiFab* r0, MultiFab* p0,
                      const amrex::Real* dptr_rayleigh_tau,
                      const amrex::Real* dptr_rayleigh_ubar,
                      const amrex::Real* dptr_rayleigh_vbar,
                      const amrex::Real* dptr_rayleigh_thetabar)
{
    BL_PROFILE_VAR("erf_advance()",erf_advance);
    if (verbose) amrex::Print() << "Starting advance at level " << level << std::endl;

    int nvars = cons_old.nComp();

    bool use_fluxes = (finest_level > 0);

    const BoxArray& ba            = cons_old.boxArray();
    const DistributionMapping& dm = cons_old.DistributionMap();

#include "TI_utils.H"

    MultiFab S_prim(ba, dm, NUM_PRIM, cons_old.nGrowVect());

    // **************************************************************************************
    // These are temporary arrays that we use to store the accumulation of the fluxes
    // **************************************************************************************
    std::array< MultiFab, AMREX_SPACEDIM >  advflux;
    std::array< MultiFab, AMREX_SPACEDIM > diffflux;

     advflux[0].define(convert(ba,IntVect(1,0,0)), dm, nvars, 0);
     advflux[1].define(convert(ba,IntVect(0,1,0)), dm, nvars, 0);
     advflux[2].define(convert(ba,IntVect(0,0,1)), dm, nvars, 0);

     advflux[0].setVal(0.);
     advflux[1].setVal(0.);
     advflux[2].setVal(0.);

    diffflux[0].define(convert(ba,IntVect(1,0,0)), dm, nvars, 0);
    diffflux[1].define(convert(ba,IntVect(0,1,0)), dm, nvars, 0);
    diffflux[2].define(convert(ba,IntVect(0,0,1)), dm, nvars, 0);

    diffflux[0].setVal(0.);
    diffflux[1].setVal(0.);
    diffflux[2].setVal(0.);

    // **************************************************************************************
    // Here we define state_old and state_new which are to be advanced
    // **************************************************************************************
    // Initial solution
    amrex::Vector<amrex::MultiFab> state_old;
    state_old.push_back(MultiFab(cons_old, amrex::make_alias, 0, nvars)); // cons
    state_old.push_back(MultiFab(convert(ba,IntVect(1,0,0)), dm, 1, xvel_old.nGrow())); // xmom
    state_old.push_back(MultiFab(convert(ba,IntVect(0,1,0)), dm, 1, yvel_old.nGrow())); // ymom
    state_old.push_back(MultiFab(convert(ba,IntVect(0,0,1)), dm, 1, zvel_old.nGrow())); // zmom
    if (use_fluxes) {
        state_old.push_back(MultiFab(convert(ba,IntVect(1,0,0)), dm, nvars, 1)); // x-fluxes
        state_old[IntVar::xflux].setVal(0.0_rt);

        state_old.push_back(MultiFab(convert(ba,IntVect(0,1,0)), dm, nvars, 1)); // y-fluxes
        state_old[IntVar::yflux].setVal(0.0_rt);

        state_old.push_back(MultiFab(convert(ba,IntVect(0,0,1)), dm, nvars, 1)); // z-fluxes
        state_old[IntVar::zflux].setVal(0.0_rt);
    }

    // Final solution
    amrex::Vector<amrex::MultiFab> state_new;
    state_new.push_back(MultiFab(cons_new, amrex::make_alias, 0, nvars)); // cons
    state_new.push_back(MultiFab(convert(ba,IntVect(1,0,0)), dm, 1, xvel_old.nGrow())); // xmom
    state_new.push_back(MultiFab(convert(ba,IntVect(0,1,0)), dm, 1, yvel_old.nGrow())); // ymom
    state_new.push_back(MultiFab(convert(ba,IntVect(0,0,1)), dm, 1, zvel_old.nGrow())); // zmom
    if (use_fluxes) {
        state_new.push_back(MultiFab(flux[0], amrex::make_alias, 0, nvars)); // x-fluxes
        state_new.push_back(MultiFab(flux[1], amrex::make_alias, 0, nvars)); // y-fluxes
        state_new.push_back(MultiFab(flux[2], amrex::make_alias, 0, nvars)); // z-fluxes
    }

    // ***********************************************************************************************
    // Convert old velocity available on faces to old momentum on faces to be used in time integration
    // ***********************************************************************************************

    VelocityToMomentum(xvel_old, xvel_old.nGrowVect(),
                       yvel_old, yvel_old.nGrowVect(),
                       zvel_old, zvel_old.nGrowVect(),
                       state_old[IntVar::cons],
                       state_old[IntVar::xmom],
                       state_old[IntVar::ymom],
                       state_old[IntVar::zmom]);

    bool fast_only = false;
    apply_bcs(state_old, old_time, state_old[IntVar::cons].nGrow(), state_old[IntVar::xmom].nGrow(), fast_only);
    cons_to_prim(state_old[IntVar::cons], S_prim, state_old[IntVar::cons].nGrow());

#include "TI_slow_rhs_fun.H"
#include "TI_fast_rhs_fun.H"

    // ***************************************************************************************
    // Setup the integrator and integrate for a single timestep
    // **************************************************************************************
    if (use_native_mri) {

        MRISplitIntegrator<Vector<MultiFab> >& mri_integrator = *mri_integrator_mem[level];

        // Define rhs and 'post update' utility function that is called after calculating
        // any state data (e.g. at RK stages or at the end of a timestep)
        mri_integrator.set_slow_rhs(slow_rhs_fun);
        mri_integrator.set_post_update(post_update_fun);

        mri_integrator.set_fast_rhs(fast_rhs_fun);
        mri_integrator.set_slow_fast_timestep_ratio(fixed_mri_dt_ratio > 0 ? fixed_mri_dt_ratio : dt_mri_ratio[level]);
        mri_integrator.set_post_substep(post_substep_fun);

        mri_integrator.advance(state_old, state_new, old_time, dt_advance);

    } else {
        SRIIntegrator<Vector<MultiFab> >& sri_integrator = *sri_integrator_mem[level];

        sri_integrator.set_rhs(slow_rhs_fun);
        sri_integrator.set_post_update(post_update_fun);

        sri_integrator.advance(state_old, state_new, old_time, dt_advance);
    }

    // **************************************************************************************
    // Convert updated momentum to updated velocity on faces after we have taken a timestep
    // **************************************************************************************
    MomentumToVelocity(xvel_new, IntVect::TheZeroVector(),
                       yvel_new, IntVect::TheZeroVector(),
                       zvel_new, IntVect::TheZeroVector(),
                       state_new[IntVar::cons],
                       state_new[IntVar::xmom],
                       state_new[IntVar::ymom],
                       state_new[IntVar::zmom]);

    if (verbose) Print() << "Done with advance at level " << level << std::endl;
}
