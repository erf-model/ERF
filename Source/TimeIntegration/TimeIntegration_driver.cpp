#include <AMReX.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>
#include <SpatialStencils.H>
#include <AMReX_TimeIntegrator.H>
#include <ERF_MRI.H>
#include <ERF_SRI.H>
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

    const BoxArray& ba            = cons_old.boxArray();
    const DistributionMapping& dm = cons_old.DistributionMap();

    // **************************************************************************************
    // Temporary array that we use to store primitive advected quantities for the RHS
    // **************************************************************************************
    auto cons_to_prim = [&](const MultiFab& cons_state, MultiFab& prim_state, int ng) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(cons_state,TilingIfNotGPU()); mfi.isValid(); ++mfi) 
      {
          // const Box& gbx = mfi.growntilebox(cons_state.nGrowVect());
          const Box& gbx = mfi.growntilebox(ng);
          const Array4<const Real>& cons_arr = cons_state.array(mfi);
          const Array4<Real>& prim_arr = prim_state.array(mfi);

          amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int n = 0; n < NUM_PRIM; ++n) {
              prim_arr(i,j,k,PrimTheta_comp + n) = cons_arr(i,j,k,RhoTheta_comp + n) / cons_arr(i,j,k,Rho_comp);
            }
          });
      } // mfi
    };

    MultiFab S_prim(ba, dm, NUM_PRIM, cons_old.nGrowVect());

    // **************************************************************************************
    // These are temporary arrays that we use to store the accumulation of the fluxes
    // **************************************************************************************
    std::array< MultiFab, AMREX_SPACEDIM >  advflux;
    std::array< MultiFab, AMREX_SPACEDIM > diffflux;

     advflux[0].define(convert(ba,IntVect(1,0,0)), dm, nvars, 0);
    diffflux[0].define(convert(ba,IntVect(1,0,0)), dm, nvars, 0);

     advflux[1].define(convert(ba,IntVect(0,1,0)), dm, nvars, 0);
    diffflux[1].define(convert(ba,IntVect(0,1,0)), dm, nvars, 0);

     advflux[2].define(convert(ba,IntVect(0,0,1)), dm, nvars, 0);
    diffflux[2].define(convert(ba,IntVect(0,0,1)), dm, nvars, 0);

     advflux[0].setVal(0.);
     advflux[1].setVal(0.);
     advflux[2].setVal(0.);

    diffflux[0].setVal(0.);
    diffflux[1].setVal(0.);
    diffflux[2].setVal(0.);

    // **************************************************************************************
    // Here we define state_old and state_new which are to be advanced
    // **************************************************************************************
    // Initial solution
    amrex::Vector<amrex::MultiFab> state_old;
    state_old.push_back(MultiFab(ba, dm, nvars, cons_old.nGrow())); // cons
    state_old.push_back(MultiFab(convert(ba,IntVect(1,0,0)), dm, 1, xvel_old.nGrow())); // xmom
    state_old.push_back(MultiFab(convert(ba,IntVect(0,1,0)), dm, 1, yvel_old.nGrow())); // ymom
    state_old.push_back(MultiFab(convert(ba,IntVect(0,0,1)), dm, 1, zvel_old.nGrow())); // zmom
    state_old.push_back(MultiFab(convert(ba,IntVect(1,0,0)), dm, nvars, 1)); // x-fluxes
    state_old.push_back(MultiFab(convert(ba,IntVect(0,1,0)), dm, nvars, 1)); // y-fluxes
    state_old.push_back(MultiFab(convert(ba,IntVect(0,0,1)), dm, nvars, 1)); // z-fluxes

    // Final solution
    amrex::Vector<amrex::MultiFab> state_new;
    state_new.push_back(MultiFab(ba, dm, nvars, cons_old.nGrow())); // cons
    state_new.push_back(MultiFab(convert(ba,IntVect(1,0,0)), dm, 1, xvel_old.nGrow())); // xmom
    state_new.push_back(MultiFab(convert(ba,IntVect(0,1,0)), dm, 1, yvel_old.nGrow())); // ymom
    state_new.push_back(MultiFab(convert(ba,IntVect(0,0,1)), dm, 1, zvel_old.nGrow())); // zmom
    state_new.push_back(MultiFab(convert(ba,IntVect(1,0,0)), dm, nvars, 1)); // x-fluxes
    state_new.push_back(MultiFab(convert(ba,IntVect(0,1,0)), dm, nvars, 1)); // y-fluxes
    state_new.push_back(MultiFab(convert(ba,IntVect(0,0,1)), dm, nvars, 1)); // z-fluxes

    // ***********************************************************************************************
    // Prepare the old-time data for calling the integrator
    // Note that we have filled the ghost cells of cons_old and we are copying the ghost cell values
    //      so we don't need to enforce BCs again here before calling VelocityToMomentum
    // ***********************************************************************************************
    MultiFab::Copy(state_old[IntVar::cons], cons_old, 0, 0, cons_old.nComp(), cons_old.nGrow());

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

    // ***********************************************************************************************
    // Initialize the fluxes to zero
    // ***********************************************************************************************
    state_old[IntVar::xflux].setVal(0.0_rt);
    state_old[IntVar::yflux].setVal(0.0_rt);
    state_old[IntVar::zflux].setVal(0.0_rt);

    // ***************************************************************************************
    // This routine is called before the first step of the time integration, *and* in the case
    //  of a multi-stage method like RK3, this is called from "post_update_fun" which is called
    //  before every subsequent stage.  Since we advance the variables in conservative form,
    //  we must convert momentum to velocity before imposing the bcs.
    // ***************************************************************************************
    auto apply_bcs = [&](Vector<MultiFab>& S_data, const Real time_for_fp, int ng_cons, int ng_vel)
    {
        amrex::Array<const MultiFab*,3> cmf_const{&xmom_crse, &ymom_crse, &zmom_crse};
        amrex::Array<MultiFab*,3> fmf{&S_data[IntVar::xmom],
                                      &S_data[IntVar::ymom],
                                      &S_data[IntVar::zmom]};

        // ***************************************************************************************
        // Interpolate momentum from coarse faces to fine faces *only* on the coarse-fine boundary
        // ***************************************************************************************
        if (level > 0) {
            ifr->interp(fmf,cmf_const,0,1);
        }

        // ***************************************************************************************
        // Call FillPatch routines for the density only because we need it to convert between
        //      momentum and velocity
        // This fills ghost cells/faces from
        //     1) coarser level if lev > 0
        //     2) physical boundaries
        //     3) other grids at the same level
        // ***************************************************************************************
        bool rho_only = true;
        FillIntermediatePatch(level, time_for_fp, {S_data[IntVar::cons], xvel_new, yvel_new, zvel_new},
                              ng_cons, ng_vel, rho_only);

        // Here we don't use include any of the ghost region because we have only updated
        //      momentum on valid faces
        MomentumToVelocity(xvel_new, IntVect::TheZeroVector(),
                           yvel_new, IntVect::TheZeroVector(),
                           zvel_new, IntVect::TheZeroVector(),
                           S_data[IntVar::cons],
                           S_data[IntVar::xmom],
                           S_data[IntVar::ymom],
                           S_data[IntVar::zmom]);

        // ***************************************************************************************
        // Call FillPatch routines for all data
        // This fills ghost cells/faces from
        //     1) coarser level if lev > 0
        //     2) physical boundaries
        //     3) other grids at the same level
        // ***************************************************************************************
        FillIntermediatePatch(level, time_for_fp, {S_data[IntVar::cons], xvel_new, yvel_new, zvel_new}, 
                              ng_cons, ng_vel);

        // Now we can convert back to momentum on valid+ghost since we have
        //     filled the ghost regions for both velocity and density
        // VelocityToMomentum(xvel_new, xvel_new.nGrowVect(),
        //                    yvel_new, yvel_new.nGrowVect(),
        //                    zvel_new, zvel_new.nGrowVect(),
        VelocityToMomentum(xvel_new, IntVect(ng_vel,ng_vel,ng_vel),
                           yvel_new, IntVect(ng_vel,ng_vel,ng_vel),
                           zvel_new, IntVect(ng_vel,ng_vel,0),
                           S_data[IntVar::cons],
                           S_data[IntVar::xmom],
                           S_data[IntVar::ymom],
                           S_data[IntVar::zmom]);
    };

    apply_bcs(state_old, old_time, state_old[IntVar::cons].nGrow(), state_old[IntVar::xmom].nGrow());
    cons_to_prim(state_old[IntVar::cons], S_prim, state_old[IntVar::cons].nGrow());

    //Create function lambdas
    auto slow_rhs_fun = [&](      Vector<MultiFab>& S_rhs,
                            const Vector<MultiFab>& S_data,
                                  Vector<MultiFab>& S_scratch,
                            const Real time,
                            const int rhs_vars=RHSVar::all) {
        if (verbose) Print() << "Calling slow rhs at level " << level << ", time = " << time << std::endl;
        erf_slow_rhs(level, S_rhs, S_data, S_prim, S_scratch,
                     xvel_new, yvel_new, zvel_new,
                     source, advflux, diffflux,
                     fine_geom, ifr, solverChoice, m_most, domain_bcs_type_d,
                     z_phys_nd[level], detJ_cc[level], r0, p0,
                     dptr_rayleigh_tau, dptr_rayleigh_ubar,
                     dptr_rayleigh_vbar, dptr_rayleigh_thetabar,
                     rhs_vars);
    };

    auto fast_rhs_fun = [&](      Vector<MultiFab>& S_rhs,
                                  Vector<MultiFab>& S_slow_rhs,
                                  Vector<MultiFab>& S_stage_data,
                            const Vector<MultiFab>& S_data,
                                  Vector<MultiFab>& S_scratch,
                            const Real fast_dt, const Real inv_fac)
    {
        if (verbose) amrex::Print() << "Calling fast rhs at level " << level << " with dt = " << fast_dt << std::endl;
        erf_fast_rhs(level, S_rhs, S_slow_rhs, S_stage_data, S_prim,
                     S_data, S_scratch, advflux, fine_geom, ifr, solverChoice,
                     z_phys_nd[level], detJ_cc[level], r0, p0,
                     fast_dt, inv_fac);
    };

    auto post_update_fun = [&](Vector<MultiFab>& S_data, const Real time_for_fp, int ng_cons, int ng_vel)
    {
        apply_bcs(S_data, time_for_fp, ng_cons, ng_vel);
        cons_to_prim(S_data[IntVar::cons], S_prim, ng_cons);
    };

    auto post_substep_fun = [&](Vector<MultiFab>& S_data, const Real time_for_fp, int ng_cons, int ng_vel)
    {
        // TODO: we only need to apply bcs for the "fast" variables -- this would be an optimization
        //       but shouldn't affect correctness
        apply_bcs(S_data, time_for_fp, ng_cons, ng_vel);
    };

    // ***************************************************************************************
    // Setup the integrator
    // **************************************************************************************
    if (use_native_mri) {
      MRISplitIntegrator<Vector<MultiFab> > mri_integrator(state_old);

      // Define rhs and 'post update' utility function that is called after calculating
      // any state data (e.g. at RK stages or at the end of a timestep)
      mri_integrator.set_slow_rhs(slow_rhs_fun);
      mri_integrator.set_post_update(post_update_fun);

      mri_integrator.set_fast_rhs(fast_rhs_fun);
      mri_integrator.set_slow_fast_timestep_ratio(fixed_mri_dt_ratio > 0 ? fixed_mri_dt_ratio : dt_mri_ratio[level]);
      mri_integrator.set_post_substep(post_substep_fun);

      // **************************************************************************************
      // Integrate for a single timestep
      // **************************************************************************************
      mri_integrator.advance(state_old, state_new, old_time, dt_advance);
    } else {
      // TimeIntegrator<Vector<MultiFab> > lev_integrator(state_old);
      SRIIntegrator<Vector<MultiFab> > sri_integrator(state_old);

      sri_integrator.set_rhs(slow_rhs_fun);
      sri_integrator.set_post_update(post_update_fun);

      // **************************************************************************************
      // Integrate for a single timestep
      // **************************************************************************************
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

    // **************************************************************************************
    // Get the final cell centered variables after the step
    // (do this at the very end because its a swap not a copy)
    // **************************************************************************************
    std::swap(cons_new, state_new[IntVar::cons]);

    std::swap(flux[0], state_new[IntVar::xflux]);
    std::swap(flux[1], state_new[IntVar::yflux]);
    std::swap(flux[2], state_new[IntVar::zflux]);

    // One final BC fill
    amrex::Real new_time = old_time + dt_advance;
    FillIntermediatePatch(level, new_time, {cons_new, xvel_new, yvel_new, zvel_new}, 
                          cons_new.nGrow(), xvel_new.nGrow());
    if (verbose) Print() << "Done with advance at level " << level << std::endl;
}
