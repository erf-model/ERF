#include <AMReX_BC_TYPES.H>
#include <AMReX_TimeIntegrator.H>
#include <ERF_MRI.H>
#include <EddyViscosity.H>
#include <EOS.H>
#include <ERF.H>
#include <TerrainMetrics.H>
//#include <TI_headers.H>
//#include <PlaneAverage.H>
#include <Diffusion.H>
#include <TileNoZ.H>
#include <Utils.H>

using namespace amrex;

/**
 * Function that advances the solution at one level for a single time step --
 * this sets up the multirate time integrator and calls the integrator's advance function
 *
 * @param[in] level level of refinement (coarsest level is 0)
 * @param[in] state_old old-time conserved variables
 * @param[in] state_new new-time conserved variables
 * @param[in] xvel_old old-time x-component of velocity
 * @param[in] yvel_old old-time y-component of velocity
 * @param[in] zvel_old old-time z-component of velocity
 * @param[in] xvel_new new-time x-component of velocity
 * @param[in] yvel_new new-time y-component of velocity
 * @param[in] zvel_new new-time z-component of velocity
 * @param[in] cc_src source term for conserved variables
 * @param[in] xmom_src source term for x-momenta
 * @param[in] ymom_src source term for y-momenta
 * @param[in] zmom_src source term for z-momenta
 * @param[in] fine_geom container for geometry information at current level
 * @param[in] dt_advance time step for this time advance
 * @param[in] old_time old time for this time advance
 */

void ERF::advance_dycore(int level,
                         Vector<MultiFab>& state_old,
                         Vector<MultiFab>& state_new,
                         MultiFab& xvel_old, MultiFab& yvel_old, MultiFab& zvel_old,
                         MultiFab& xvel_new, MultiFab& yvel_new, MultiFab& zvel_new,
                         MultiFab&   cc_src, MultiFab& xmom_src,
                         MultiFab& ymom_src, MultiFab& zmom_src,
                         const Geometry fine_geom,
                         const Real dt_advance, const Real old_time)
{
    BL_PROFILE_VAR("erf_advance_dycore()",erf_advance_dycore);

    const Box& domain = fine_geom.Domain();

    DiffChoice dc    = solverChoice.diffChoice;
    TurbChoice tc    = solverChoice.turbChoice[level];
    SpongeChoice sc  = solverChoice.spongeChoice;

    MultiFab r_hse (base_state[level], make_alias, 0, 1); // r_0 is first  component
    MultiFab p_hse (base_state[level], make_alias, 1, 1); // p_0 is second component
    MultiFab pi_hse(base_state[level], make_alias, 2, 1); // pi_0 is second component

    // These pointers are used in the MRI utility functions
    MultiFab* r0  = &r_hse;
    MultiFab* p0  = &p_hse;
    MultiFab* pi0 = &pi_hse;

    Real* dptr_rhotheta_src = solverChoice.custom_rhotheta_forcing ? d_rhotheta_src[level].data() : nullptr;
    Real* dptr_rhoqt_src    = solverChoice.custom_moisture_forcing ? d_rhoqt_src[level].data()    : nullptr;
    Real* dptr_wbar_sub     = solverChoice.custom_w_subsidence     ? d_w_subsid[level].data()     : nullptr;

    // Turbulent Perturbation Pointer
    //Real* dptr_rhotheta_src = solverChoice.pert_type ? d_rhotheta_src[level].data() : nullptr;

    Vector<Real*> d_rayleigh_ptrs_at_lev;
    d_rayleigh_ptrs_at_lev.resize(Rayleigh::nvars);
    bool rayleigh_damp_any = (solverChoice.rayleigh_damp_U ||solverChoice.rayleigh_damp_V ||
                              solverChoice.rayleigh_damp_W ||solverChoice.rayleigh_damp_T);
    d_rayleigh_ptrs_at_lev[Rayleigh::tau]      =              rayleigh_damp_any ? d_rayleigh_ptrs[level][Rayleigh::tau ].data() : nullptr;
    d_rayleigh_ptrs_at_lev[Rayleigh::ubar]     = solverChoice.rayleigh_damp_U   ? d_rayleigh_ptrs[level][Rayleigh::ubar].data() : nullptr;
    d_rayleigh_ptrs_at_lev[Rayleigh::vbar]     = solverChoice.rayleigh_damp_V   ? d_rayleigh_ptrs[level][Rayleigh::vbar].data() : nullptr;
    d_rayleigh_ptrs_at_lev[Rayleigh::wbar]     = solverChoice.rayleigh_damp_W   ? d_rayleigh_ptrs[level][Rayleigh::wbar].data() : nullptr;
    d_rayleigh_ptrs_at_lev[Rayleigh::thetabar] = solverChoice.rayleigh_damp_T   ? d_rayleigh_ptrs[level][Rayleigh::thetabar].data() : nullptr;

    Vector<Real*> d_sponge_ptrs_at_lev;
    if(sc.sponge_type=="input_sponge")
    {
        d_sponge_ptrs_at_lev.resize(Sponge::nvars_sponge);
        d_sponge_ptrs_at_lev[Sponge::ubar_sponge]  =  d_sponge_ptrs[level][Sponge::ubar_sponge].data();
        d_sponge_ptrs_at_lev[Sponge::vbar_sponge]  =  d_sponge_ptrs[level][Sponge::vbar_sponge].data();
    }

    bool l_use_terrain = solverChoice.use_terrain;
    bool l_use_diff    = ( (dc.molec_diff_type != MolecDiffType::None) ||
                           (tc.les_type        !=       LESType::None) ||
                           (tc.pbl_type        !=       PBLType::None) );
    bool l_use_kturb   = ( (tc.les_type != LESType::None)   ||
                           (tc.pbl_type != PBLType::None) );

    const bool use_most = (m_most != nullptr);
    const bool exp_most = (solverChoice.use_explicit_most);
    amrex::ignore_unused(use_most);

    const BoxArray& ba            = state_old[IntVars::cons].boxArray();
    const BoxArray& ba_z          = zvel_old.boxArray();
    const DistributionMapping& dm = state_old[IntVars::cons].DistributionMap();

    int num_prim = state_old[IntVars::cons].nComp() - 1;

    MultiFab    S_prim  (ba  , dm, num_prim,          state_old[IntVars::cons].nGrowVect());
    MultiFab  pi_stage  (ba  , dm,        1,          state_old[IntVars::cons].nGrowVect());
    MultiFab fast_coeffs(ba_z, dm,        5,          0);
    MultiFab* eddyDiffs = eddyDiffs_lev[level].get();
    MultiFab* SmnSmn    = SmnSmn_lev[level].get();

    // **************************************************************************************
    // Compute strain for use in slow RHS, Smagorinsky model, and MOST
    // **************************************************************************************
    {
    BL_PROFILE("erf_advance_strain");
    if (l_use_diff) {

        const BCRec* bc_ptr_h = domain_bcs_type.data();
        const GpuArray<Real, AMREX_SPACEDIM> dxInv = fine_geom.InvCellSizeArray();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(state_new[IntVars::cons],TileNoZ()); mfi.isValid(); ++mfi)
        {
            Box bxcc  = mfi.growntilebox(IntVect(1,1,0));
            Box tbxxy = mfi.tilebox(IntVect(1,1,0),IntVect(1,1,0));
            Box tbxxz = mfi.tilebox(IntVect(1,0,1),IntVect(1,1,0));
            Box tbxyz = mfi.tilebox(IntVect(0,1,1),IntVect(1,1,0));

            const Array4<const Real> & u = xvel_old.array(mfi);
            const Array4<const Real> & v = yvel_old.array(mfi);
            const Array4<const Real> & w = zvel_old.array(mfi);

            Array4<Real> tau11 = Tau11_lev[level].get()->array(mfi);
            Array4<Real> tau22 = Tau22_lev[level].get()->array(mfi);
            Array4<Real> tau33 = Tau33_lev[level].get()->array(mfi);
            Array4<Real> tau12 = Tau12_lev[level].get()->array(mfi);
            Array4<Real> tau13 = Tau13_lev[level].get()->array(mfi);
            Array4<Real> tau23 = Tau23_lev[level].get()->array(mfi);

            Array4<Real> tau21  = l_use_terrain ? Tau21_lev[level].get()->array(mfi) : Array4<Real>{};
            Array4<Real> tau31  = l_use_terrain ? Tau31_lev[level].get()->array(mfi) : Array4<Real>{};
            Array4<Real> tau32  = l_use_terrain ? Tau32_lev[level].get()->array(mfi) : Array4<Real>{};
            const Array4<const Real>& z_nd = l_use_terrain ? z_phys_nd[level]->const_array(mfi) : Array4<const Real>{};

            const Array4<const Real> mf_m = mapfac_m[level]->array(mfi);
            const Array4<const Real> mf_u = mapfac_u[level]->array(mfi);
            const Array4<const Real> mf_v = mapfac_v[level]->array(mfi);

            if (l_use_terrain) {
                ComputeStrain_T(bxcc, tbxxy, tbxxz, tbxyz, domain,
                                u, v, w,
                                tau11, tau22, tau33,
                                tau12, tau13,
                                tau21, tau23,
                                tau31, tau32,
                                z_nd, detJ_cc[level]->const_array(mfi), bc_ptr_h, dxInv,
                                mf_m, mf_u, mf_v);
            } else {
                ComputeStrain_N(bxcc, tbxxy, tbxxz, tbxyz, domain,
                                u, v, w,
                                tau11, tau22, tau33,
                                tau12, tau13, tau23,
                                bc_ptr_h, dxInv,
                                mf_m, mf_u, mf_v);
            }
        } // mfi
    } // l_use_diff
    } // profile

    MultiFab Omega (state_old[IntVars::zmom].boxArray(),dm,1,1);

#include "TI_utils.H"

    // Additional SFS quantities, calculated once per timestep
    MultiFab* Hfx1 = SFS_hfx1_lev[level].get();
    MultiFab* Hfx2 = SFS_hfx2_lev[level].get();
    MultiFab* Hfx3 = SFS_hfx3_lev[level].get();
    MultiFab* Diss = SFS_diss_lev[level].get();

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
    if (l_use_kturb)
    {
        // NOTE: state_new transfers to state_old for PBL (due to ptr swap in advance)
        const BCRec* bc_ptr_d = domain_bcs_type_d.data();
        ComputeTurbulentViscosity(xvel_old, yvel_old,
                                  *Tau11_lev[level].get(), *Tau22_lev[level].get(), *Tau33_lev[level].get(),
                                  *Tau12_lev[level].get(), *Tau13_lev[level].get(), *Tau23_lev[level].get(),
                                  state_old[IntVars::cons],
                                  *eddyDiffs, *Hfx1, *Hfx2, *Hfx3, *Diss, // to be updated
                                  fine_geom, *mapfac_u[level], *mapfac_v[level],
                                  z_phys_nd[level], tc, solverChoice.gravity,
                                  m_most, exp_most, level, bc_ptr_d);
    }

    // ***********************************************************************************************
    // Update user-defined source terms -- these are defined once per time step (not per RK stage)
    // ***********************************************************************************************
    if (solverChoice.custom_rhotheta_forcing) {
        prob->update_rhotheta_sources(old_time,
                                      h_rhotheta_src[level], d_rhotheta_src[level],
                                      fine_geom, z_phys_cc[level]);
    }

    if (solverChoice.custom_moisture_forcing) {
        prob->update_rhoqt_sources(old_time,
                                   h_rhoqt_src[level], d_rhoqt_src[level],
                                   fine_geom, z_phys_cc[level]);
    }

    if (solverChoice.custom_geostrophic_profile) {
        prob->update_geostrophic_profile(old_time,
                                   h_u_geos[level], d_u_geos[level],
                                   h_v_geos[level], d_v_geos[level],
                                   fine_geom, z_phys_cc[level]);
    }

    // ***********************************************************************************************
    // Convert old velocity available on faces to old momentum on faces to be used in time integration
    // ***********************************************************************************************
    MultiFab density(state_old[IntVars::cons], make_alias, Rho_comp, 1);

    //
    // This is an optimization since we won't need more than one ghost
    // cell of momentum in the integrator if not using NumDiff
    //
    IntVect ngu = (solverChoice.use_NumDiff) ? IntVect(1,1,1) : xvel_old.nGrowVect();
    IntVect ngv = (solverChoice.use_NumDiff) ? IntVect(1,1,1) : yvel_old.nGrowVect();
    IntVect ngw = (solverChoice.use_NumDiff) ? IntVect(1,1,0) : zvel_old.nGrowVect();
    VelocityToMomentum(xvel_old, ngu, yvel_old, ngv, zvel_old, ngw, density,
                       state_old[IntVars::xmom],
                       state_old[IntVars::ymom],
                       state_old[IntVars::zmom],
                       domain, domain_bcs_type);

    MultiFab::Copy(xvel_new,xvel_old,0,0,1,xvel_old.nGrowVect());
    MultiFab::Copy(yvel_new,yvel_old,0,0,1,yvel_old.nGrowVect());
    MultiFab::Copy(zvel_new,zvel_old,0,0,1,zvel_old.nGrowVect());

    bool fast_only = false;
    bool vel_and_mom_synced = true;

    apply_bcs(state_old, old_time,
              state_old[IntVars::cons].nGrow(), state_old[IntVars::xmom].nGrow(),
              fast_only, vel_and_mom_synced);
    cons_to_prim(state_old[IntVars::cons], state_old[IntVars::cons].nGrow());

#include "TI_no_substep_fun.H"
#include "TI_slow_rhs_fun.H"
#include "TI_fast_rhs_fun.H"

    // ***************************************************************************************
    // Setup the integrator and integrate for a single timestep
    // **************************************************************************************
    MRISplitIntegrator<Vector<MultiFab> >& mri_integrator = *mri_integrator_mem[level];

    // Define rhs and 'post update' utility function that is called after calculating
    // any state data (e.g. at RK stages or at the end of a timestep)
    mri_integrator.set_slow_rhs_pre(slow_rhs_fun_pre);
    mri_integrator.set_slow_rhs_post(slow_rhs_fun_post);
    mri_integrator.set_pre_update (pre_update_fun);
    mri_integrator.set_post_update(post_update_fun);

#ifdef ERF_USE_POISSON_SOLVE
    if (solverChoice.incompressible[level]) {
        mri_integrator.set_slow_rhs_inc(slow_rhs_fun_inc);
    }
#endif

    mri_integrator.set_fast_rhs(fast_rhs_fun);
    mri_integrator.set_slow_fast_timestep_ratio(fixed_mri_dt_ratio > 0 ? fixed_mri_dt_ratio : dt_mri_ratio[level]);
    mri_integrator.set_no_substep(no_substep_fun);

    mri_integrator.advance(state_old, state_new, old_time, dt_advance);

    if (verbose) Print() << "Done with advance_dycore at level " << level << std::endl;
}
