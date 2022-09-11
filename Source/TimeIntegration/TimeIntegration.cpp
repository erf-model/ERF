#include <AMReX_BC_TYPES.H>
#include <AMReX_TimeIntegrator.H>
#include <ERF_MRI.H>
#include <EddyViscosity.H>
#include <TerrainMetrics.H>
#include <TimeIntegration.H>
#include <ERF.H>
#include <EOS.H>

using namespace amrex;

void ERF::erf_advance(int level,
                      MultiFab& cons_old,  MultiFab& cons_new,
                      MultiFab& xvel_old, MultiFab& yvel_old, MultiFab& zvel_old,
                      MultiFab& xvel_new, MultiFab& yvel_new, MultiFab& zvel_new,
                      MultiFab& xmom_old, MultiFab& ymom_old, MultiFab& zmom_old,
                      MultiFab& xmom_new, MultiFab& ymom_new, MultiFab& zmom_new,
                      MultiFab& xmom_crse, MultiFab& ymom_crse, MultiFab& zmom_crse,
                      MultiFab& source,
                      const amrex::Geometry fine_geom,
                      const amrex::Real dt_advance, const amrex::Real old_time,
                      amrex::InterpFaceRegister* ifr,
                      MultiFab* r0, MultiFab* p0, MultiFab* pi0,
                      const amrex::Real* dptr_rayleigh_tau,
                      const amrex::Real* dptr_rayleigh_ubar,
                      const amrex::Real* dptr_rayleigh_vbar,
                      const amrex::Real* dptr_rayleigh_thetabar)
{
    BL_PROFILE_VAR("erf_advance()",erf_advance);
    if (verbose) amrex::Print() << "Starting advance at level " << level << std::endl;

    int nvars = cons_old.nComp();

    const BoxArray& ba            = cons_old.boxArray();
    const BoxArray& ba_z          = zvel_old.boxArray();
    const DistributionMapping& dm = cons_old.DistributionMap();

    MultiFab    S_prim  (ba  , dm, NUM_PRIM,          cons_old.nGrowVect());
    MultiFab  pi_stage  (ba  , dm,        1,          cons_old.nGrowVect());
    MultiFab fast_coeffs(ba_z, dm,        5,          0);
    MultiFab eddyDiffs  (ba  , dm, EddyDiff::NumDiffs,1);

    MultiFab Omega (zmom_old.boxArray(),dm,1,1);

#include "TI_utils.H"

    amrex::Vector<amrex::MultiFab> state_old;
    amrex::Vector<amrex::MultiFab> state_new;

    // **************************************************************************************
    // These are temporary arrays that we use to store the accumulation of the fluxes
    // **************************************************************************************
    std::array< MultiFab, AMREX_SPACEDIM > diffflux;

    {
    BL_PROFILE("erf_advance_part_1");
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
    state_old.push_back(MultiFab(cons_old, amrex::make_alias, 0, nvars)); // cons
    state_old.push_back(MultiFab(xmom_old, amrex::make_alias, 0,     1)); // xmom
    state_old.push_back(MultiFab(ymom_old, amrex::make_alias, 0,     1)); // ymom
    state_old.push_back(MultiFab(zmom_old, amrex::make_alias, 0,     1)); // zmom

    // Final solution
    state_new.push_back(MultiFab(cons_new, amrex::make_alias, 0, nvars)); // cons
    state_new.push_back(MultiFab(xmom_new, amrex::make_alias, 0,     1)); // xmom
    state_new.push_back(MultiFab(ymom_new, amrex::make_alias, 0,     1)); // ymom
    state_new.push_back(MultiFab(zmom_new, amrex::make_alias, 0,     1)); // zmom
    } // end profile

    // *************************************************************************
    // Calculate cell-centered eddy viscosity & diffusivities
    //
    // Notes -- we fill all the data in ghost cells before calling this so
    //    that we can fill the eddy viscosity in the ghost regions and
    //    not have to call a boundary filler on this data itself
    //
    // LES - updates both horizontal and vertical eddy viscosity components
    // PBL - only updates vertical eddy viscosity components so horizontal
    //       components come from the LES model or are left as zero.
    // *************************************************************************
    if ( (solverChoice.les_type        !=       LESType::None) ||
         (solverChoice.pbl_type        !=       PBLType::None) )
    {
        ComputeTurbulentViscosity(xvel_old, yvel_old, zvel_old, state_old[IntVar::cons],
                                  eddyDiffs, fine_geom, solverChoice, m_most, domain_bcs_type_d);
    }
    // *************************************************************************

    // ***********************************************************************************************
    // Convert old velocity available on faces to old momentum on faces to be used in time integration
    // ***********************************************************************************************

    {
    BL_PROFILE("pre_set_up_mri");
    VelocityToMomentum(xvel_old, xvel_old.nGrowVect(),
                       yvel_old, yvel_old.nGrowVect(),
                       zvel_old, zvel_old.nGrowVect(),
                       state_old[IntVar::cons],
                       state_old[IntVar::xmom],
                       state_old[IntVar::ymom],
                       state_old[IntVar::zmom]);

    Real time_mt   = old_time - 0.5*dt_advance; // Moving terrain
    bool fast_only = false;
    apply_bcs(state_old, old_time, time_mt, dt_advance,
              state_old[IntVar::cons].nGrow(), state_old[IntVar::xmom].nGrow(), fast_only);
    cons_to_prim(state_old[IntVar::cons], state_old[IntVar::cons].nGrow());
    }

#include "TI_no_substep_fun.H"
#include "TI_slow_rhs_fun.H"
#include "TI_fast_rhs_fun.H"

    // ***************************************************************************************
    // Setup the integrator and integrate for a single timestep
    // **************************************************************************************
    MRISplitIntegrator<Vector<MultiFab> >& mri_integrator = *mri_integrator_mem[level];

    {
    BL_PROFILE("set_up_mri_integrator");
    // Define rhs and 'post update' utility function that is called after calculating
    // any state data (e.g. at RK stages or at the end of a timestep)
    mri_integrator.set_slow_rhs_pre(slow_rhs_fun_pre);
    mri_integrator.set_slow_rhs_post(slow_rhs_fun_post);
    mri_integrator.set_pre_update (pre_update_fun);
    mri_integrator.set_post_update(post_update_fun);

    mri_integrator.set_fast_rhs(fast_rhs_fun);
    mri_integrator.set_slow_fast_timestep_ratio(fixed_mri_dt_ratio > 0 ? fixed_mri_dt_ratio : dt_mri_ratio[level]);
    mri_integrator.set_no_substep(no_substep_fun);
    }

    mri_integrator.advance(state_old, state_new, old_time, dt_advance);

    if (verbose) Print() << "Done with advance at level " << level << std::endl;
}
