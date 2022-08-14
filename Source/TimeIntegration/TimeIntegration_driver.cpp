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

    // **************************************************************************************
    // Temporary array that we use to store primitive advected quantities for the RHS
    // **************************************************************************************
    auto cons_to_prim = [&](const MultiFab& cons_state, MultiFab& prim_state, int ng) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(cons_state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
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
        state_old.push_back(MultiFab(convert(ba,IntVect(0,1,0)), dm, nvars, 1)); // y-fluxes
        state_old.push_back(MultiFab(convert(ba,IntVect(0,0,1)), dm, nvars, 1)); // z-fluxes
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

    // ***********************************************************************************************
    // Initialize the fluxes to zero
    // ***********************************************************************************************
    if (use_fluxes) {
        state_old[IntVar::xflux].setVal(0.0_rt);
        state_old[IntVar::yflux].setVal(0.0_rt);
        state_old[IntVar::zflux].setVal(0.0_rt);
    }

    // ***************************************************************************************
    // This routine is called before the first step of the time integration, *and* in the case
    //  of a multi-stage method like RK3, this is called from "post_update_fun" which is called
    //  before every subsequent stage.  Since we advance the variables in conservative form,
    //  we must convert momentum to velocity before imposing the bcs.
    // ***************************************************************************************
    auto apply_bcs = [&](Vector<MultiFab>& S_data, const Real time_for_fp, int ng_cons, int ng_vel, bool fast_only)
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
        int scomp_cons = 0;
        int ncomp_cons = 1;
        bool cons_only = true;
        FillIntermediatePatch(level, time_for_fp, {S_data[IntVar::cons], xvel_new, yvel_new, zvel_new},
                              ng_cons, 0, cons_only, scomp_cons, ncomp_cons);

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
        scomp_cons = 0;
        if (fast_only) {
            ncomp_cons = 2; // rho and (rho theta)
        } else {
            ncomp_cons = S_data[IntVar::cons].nComp();
        }
        cons_only = false;
        FillIntermediatePatch(level, time_for_fp, {S_data[IntVar::cons], xvel_new, yvel_new, zvel_new},
                              ng_cons, ng_vel, cons_only, scomp_cons, ncomp_cons);

        // Now we can convert back to momentum on valid+ghost since we have
        //     filled the ghost regions for both velocity and density
        VelocityToMomentum(xvel_new, IntVect(ng_vel,ng_vel,ng_vel),
                           yvel_new, IntVect(ng_vel,ng_vel,ng_vel),
                           zvel_new, IntVect(ng_vel,ng_vel,0),
                           S_data[IntVar::cons],
                           S_data[IntVar::xmom],
                           S_data[IntVar::ymom],
                           S_data[IntVar::zmom]);
    };

    bool fast_only = false;
    apply_bcs(state_old, old_time, state_old[IntVar::cons].nGrow(), state_old[IntVar::xmom].nGrow(), fast_only);
    cons_to_prim(state_old[IntVar::cons], S_prim, state_old[IntVar::cons].nGrow());

    //Create function lambdas
    auto slow_rhs_fun = [&](      Vector<MultiFab>& S_rhs,
                            const Vector<MultiFab>& S_data,
                                  Vector<MultiFab>& S_scratch,
                            const Real step_time, const Real stage_time,
                            const int rhs_vars=RHSVar::all)
    {
        if (verbose) Print() << "Calling slow rhs at level " << level << " for advancing from "
                             << step_time << " to " << stage_time << std::endl;

        // Moving terrain
        if ( solverChoice.use_terrain &&  (solverChoice.terrain_type == 1) )
        {
            if (verbose) Print() << "Making new geometry at stage_time " << std::endl;

            const auto dz = fine_geom.CellSizeArray()[2];
            const Real old_time = step_time;

            init_custom_terrain(fine_geom,*z_phys_nd_new[level],stage_time); // This defines h(i,j,k=0,t)
            init_terrain_grid  (fine_geom,*z_phys_nd_new[level]);        // This defines z_phys for all k given h
            make_metrics       (fine_geom,*z_phys_nd_new[level], *z_phys_cc_new[level], *detJ_cc_new[level]);  // This defines the other quantities

            MultiFab* z_t = new MultiFab(S_data[IntVar::zmom].boxArray(), S_data[IntVar::zmom].DistributionMap(), 1, 1);

            for (MFIter mfi(*z_t,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box gbx = mfi.growntilebox(IntVect(1,1,0));

                // This is the top face
                int khi = gbx.bigEnd()[2];

                // Now make this a 2D box (after we extract khi)
                gbx.setRange(2,0);

                const Array4<const Real>& detJ_new  =   detJ_cc_new[level]->const_array(mfi);
                const Array4<const Real>& detJ_old  =       detJ_cc[level]->const_array(mfi);

                const Array4<      Real>& z_t_arr   =  z_t->array(mfi);

                Real delta_t  = (stage_time - step_time);

                const auto& dX = fine_geom.CellSizeArray();

                // Loop over horizontal plane
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    z_t_arr(i,j,0) = dhdt(i,j,dX,old_time);
                    for (int k = 1; k <= khi; k++)
                    {
                        Real Jinv_new = 1./detJ_new(i,j,k-1);
                        Real Jinv_old = 1./detJ_old(i,j,k-1);
                        z_t_arr(i,j,k) = z_t_arr(i,j,k-1) + ( (Jinv_new - Jinv_old) / delta_t ) * dz;
                    }
                });
            } // mfi

            erf_slow_rhs(level, S_rhs, S_data, S_prim, S_scratch,
                         xvel_new, yvel_new, zvel_new, z_t,
                         source, advflux, diffflux,
                         fine_geom, ifr, solverChoice, m_most, domain_bcs_type_d,
                         z_phys_nd_new[level], detJ_cc_new[level], r0, p0,
                         dptr_rayleigh_tau, dptr_rayleigh_ubar,
                         dptr_rayleigh_vbar, dptr_rayleigh_thetabar,
                         rhs_vars);

             delete z_t;

        } else {
            MultiFab* z_t = 0;
            erf_slow_rhs(level, S_rhs, S_data, S_prim, S_scratch,
                         xvel_new, yvel_new, zvel_new, z_t,
                         source, advflux, diffflux,
                         fine_geom, ifr, solverChoice, m_most, domain_bcs_type_d,
                         z_phys_nd[level], detJ_cc[level], r0, p0,
                         dptr_rayleigh_tau, dptr_rayleigh_ubar,
                         dptr_rayleigh_vbar, dptr_rayleigh_thetabar,
                         rhs_vars);
        }
    };

    auto fast_rhs_fun = [&](      Vector<MultiFab>& S_fast_rhs,
                                  Vector<MultiFab>& S_slow_rhs,
                                  Vector<MultiFab>& S_stage_data,
                            const Vector<MultiFab>& S_data,
                                  Vector<MultiFab>& S_scratch,
                            const Real fast_dt, const Real inv_fac)
    {
        if (verbose) amrex::Print() << "Calling fast rhs at level " << level << " with dt = " << fast_dt << std::endl;
        erf_fast_rhs(level, S_fast_rhs, S_slow_rhs, S_stage_data, S_prim,
                     S_data, S_scratch, advflux, fine_geom, ifr, solverChoice,
                     z_phys_nd[level], detJ_cc[level], r0, p0,
                     fast_dt, inv_fac);
    };

    auto post_update_fun = [&](Vector<MultiFab>& S_data, const Real time_for_fp, int ng_cons, int ng_vel)
    {
        bool fast_only = false;
        apply_bcs(S_data, time_for_fp, ng_cons, ng_vel, fast_only);
        cons_to_prim(S_data[IntVar::cons], S_prim, ng_cons);
    };

    auto post_substep_fun = [&](Vector<MultiFab>& S_data, const Real time_for_fp, int ng_cons, int ng_vel)
    {
        bool fast_only = true;
        apply_bcs(S_data, time_for_fp, ng_cons, ng_vel, fast_only);
    };

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
