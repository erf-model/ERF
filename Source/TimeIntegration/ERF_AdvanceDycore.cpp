#include <AMReX_BC_TYPES.H>
#include <AMReX_TimeIntegrator.H>

#include <ERF.H>
#include <ERF_MRI.H>
#include <ERF_EddyViscosity.H>
#include <ERF_EOS.H>
#include <ERF_TerrainMetrics.H>
#include <ERF_Diffusion.H>
#include <ERF_TileNoZ.H>
#include <ERF_Utils.H>
#include <ERF_EBRedistribute.H>

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

void ERF::advance_dycore (int level,
                          Vector<MultiFab>& state_old,
                          Vector<MultiFab>& state_new,
                          MultiFab& xvel_old, MultiFab& yvel_old, MultiFab& zvel_old,
                          MultiFab& xvel_new, MultiFab& yvel_new, MultiFab& zvel_new,
                          MultiFab&   cc_src, MultiFab& xmom_src,
                          MultiFab& ymom_src, MultiFab& zmom_src,
                          MultiFab& buoyancy,
                          const Geometry fine_geom,
                          const Real dt_advance, const Real old_time)
{
    BL_PROFILE_VAR("erf_advance_dycore()",erf_advance_dycore);

    const Box& domain = fine_geom.Domain();

    DiffChoice dc    = solverChoice.diffChoice;
    TurbChoice tc    = solverChoice.turbChoice[level];
    SpongeChoice sc  = solverChoice.spongeChoice;

    MultiFab r_hse (base_state[level], make_alias, BaseState::r0_comp , 1);
    MultiFab p_hse (base_state[level], make_alias, BaseState::p0_comp , 1);
    MultiFab pi_hse(base_state[level], make_alias, BaseState::pi0_comp, 1);

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
    d_rayleigh_ptrs_at_lev[Rayleigh::ubar]     = solverChoice.dampingChoice.rayleigh_damp_U   ? d_rayleigh_ptrs[level][Rayleigh::ubar].data() : nullptr;
    d_rayleigh_ptrs_at_lev[Rayleigh::vbar]     = solverChoice.dampingChoice.rayleigh_damp_V   ? d_rayleigh_ptrs[level][Rayleigh::vbar].data() : nullptr;
    d_rayleigh_ptrs_at_lev[Rayleigh::wbar]     = solverChoice.dampingChoice.rayleigh_damp_W   ? d_rayleigh_ptrs[level][Rayleigh::wbar].data() : nullptr;
    d_rayleigh_ptrs_at_lev[Rayleigh::thetabar] = solverChoice.dampingChoice.rayleigh_damp_T   ? d_rayleigh_ptrs[level][Rayleigh::thetabar].data() : nullptr;

    bool use_rayleigh =
       (solverChoice.dampingChoice.rayleigh_damp_U ||solverChoice.dampingChoice.rayleigh_damp_V ||
        solverChoice.dampingChoice.rayleigh_damp_W ||solverChoice.dampingChoice.rayleigh_damp_T);
    Real* d_sinesq_at_lev      = (use_rayleigh)  ? d_sinesq_ptrs[level].data() : nullptr;
    Real* d_sinesq_stag_at_lev = (use_rayleigh)  ? d_sinesq_stag_ptrs[level].data() : nullptr;

    Vector<Real*> d_sponge_ptrs_at_lev;
    if(sc.sponge_type=="input_sponge")
    {
        d_sponge_ptrs_at_lev.resize(Sponge::nvars_sponge);
        d_sponge_ptrs_at_lev[Sponge::ubar_sponge]  =  d_sponge_ptrs[level][Sponge::ubar_sponge].data();
        d_sponge_ptrs_at_lev[Sponge::vbar_sponge]  =  d_sponge_ptrs[level][Sponge::vbar_sponge].data();
    }

    bool l_use_terrain_fitted_coords = (solverChoice.mesh_type != MeshType::ConstantDz);
    bool l_use_kturb   = tc.use_kturb;
    bool l_use_diff    = ( (dc.molec_diff_type != MolecDiffType::None) ||
                           l_use_kturb );

    const bool use_SurfLayer = (m_SurfaceLayer != nullptr);
    const MultiFab* z_0     = (use_SurfLayer) ? m_SurfaceLayer->get_z0(level) : nullptr;

    const BoxArray& ba            = state_old[IntVars::cons].boxArray();
    const BoxArray& ba_z          = zvel_old.boxArray();
    const DistributionMapping& dm = state_old[IntVars::cons].DistributionMap();

    int num_prim = state_old[IntVars::cons].nComp() - 1;

    MultiFab    S_prim  (ba  , dm, num_prim,          state_old[IntVars::cons].nGrowVect());
    MultiFab  pi_stage  (ba  , dm,        1,          1);
    MultiFab fast_coeffs(ba_z, dm,        5,          0);

    MultiFab* eddyDiffs = eddyDiffs_lev[level].get();
    MultiFab* SmnSmn    = SmnSmn_lev[level].get();

    // **************************************************************************************
    // Compute strain for use in slow RHS and Smagorinsky model
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

            if (bxcc.smallEnd(2) != domain.smallEnd(2)) {
                 bxcc.growLo(2,1);
                tbxxy.growLo(2,1);
                tbxxz.growLo(2,1);
                tbxyz.growLo(2,1);
            }

            if (bxcc.bigEnd(2) != domain.bigEnd(2)) {
                 bxcc.growHi(2,1);
                tbxxy.growHi(2,1);
                tbxxz.growHi(2,1);
                tbxyz.growHi(2,1);
            }

            const Array4<const Real> & u = xvel_old.array(mfi);
            const Array4<const Real> & v = yvel_old.array(mfi);
            const Array4<const Real> & w = zvel_old.array(mfi);

            Array4<Real> tau11 = Tau[level][TauType::tau11].get()->array(mfi);
            Array4<Real> tau22 = Tau[level][TauType::tau22].get()->array(mfi);
            Array4<Real> tau33 = Tau[level][TauType::tau33].get()->array(mfi);
            Array4<Real> tau12 = Tau[level][TauType::tau12].get()->array(mfi);
            Array4<Real> tau13 = Tau[level][TauType::tau13].get()->array(mfi);
            Array4<Real> tau23 = Tau[level][TauType::tau23].get()->array(mfi);

            Array4<Real> tau21 = l_use_terrain_fitted_coords ? Tau[level][TauType::tau21].get()->array(mfi) : Array4<Real>{};
            Array4<Real> tau31 = l_use_terrain_fitted_coords ? Tau[level][TauType::tau31].get()->array(mfi) : Array4<Real>{};
            Array4<Real> tau32 = l_use_terrain_fitted_coords ? Tau[level][TauType::tau32].get()->array(mfi) : Array4<Real>{};
            const Array4<const Real>& z_nd = z_phys_nd[level]->const_array(mfi);

            const Array4<const Real> mf_mx = mapfac[level][MapFacType::m_x]->const_array(mfi);
            const Array4<const Real> mf_ux = mapfac[level][MapFacType::u_x]->const_array(mfi);
            const Array4<const Real> mf_vx = mapfac[level][MapFacType::v_x]->const_array(mfi);
            const Array4<const Real> mf_my = mapfac[level][MapFacType::m_y]->const_array(mfi);
            const Array4<const Real> mf_uy = mapfac[level][MapFacType::u_y]->const_array(mfi);
            const Array4<const Real> mf_vy = mapfac[level][MapFacType::v_y]->const_array(mfi);

            // We update Tau_corr[level] in erf_make_tau_terms, not here
            Array4<Real> no_tau_corr_update_here{};

            if (solverChoice.mesh_type == MeshType::StretchedDz) {
                ComputeStrain_S(bxcc, tbxxy, tbxxz, tbxyz, domain,
                                u, v, w,
                                tau11, tau22, tau33,
                                tau12, tau21,
                                tau13, tau31,
                                tau23, tau32,
                                stretched_dz_d[level], dxInv,
                                mf_mx, mf_ux, mf_vx, mf_my, mf_uy, mf_vy, bc_ptr_h,
                                no_tau_corr_update_here, no_tau_corr_update_here);
            } else if (l_use_terrain_fitted_coords) {
                ComputeStrain_T(bxcc, tbxxy, tbxxz, tbxyz, domain,
                                u, v, w,
                                tau11, tau22, tau33,
                                tau12, tau21,
                                tau13, tau31,
                                tau23, tau32,
                                z_nd, detJ_cc[level]->const_array(mfi), dxInv,
                                mf_mx, mf_ux, mf_vx, mf_my, mf_uy, mf_vy, bc_ptr_h,
                                no_tau_corr_update_here, no_tau_corr_update_here);
            } else {
                ComputeStrain_N(bxcc, tbxxy, tbxxz, tbxyz, domain,
                                u, v, w,
                                tau11, tau22, tau33,
                                tau12, tau13, tau23,
                                dxInv,
                                mf_mx, mf_ux, mf_vx, mf_my, mf_uy, mf_vy, bc_ptr_h,
                                no_tau_corr_update_here, no_tau_corr_update_here);
            }
        } // mfi
    } // l_use_diff
    } // profile

#include "ERF_TI_utils.H"

    // Additional SFS quantities, calculated once per timestep
    MultiFab* Hfx1  = SFS_hfx1_lev[level].get();
    MultiFab* Hfx2  = SFS_hfx2_lev[level].get();
    MultiFab* Hfx3  = SFS_hfx3_lev[level].get();
    MultiFab* Q1fx1 = SFS_q1fx1_lev[level].get();
    MultiFab* Q1fx2 = SFS_q1fx2_lev[level].get();
    MultiFab* Q1fx3 = SFS_q1fx3_lev[level].get();
    MultiFab* Q2fx3 = SFS_q2fx3_lev[level].get();
    MultiFab* Diss  = SFS_diss_lev[level].get();

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
        bool l_use_moisture = ( solverChoice.moisture_type != MoistureType::None );
        const BCRec* bc_ptr_h = domain_bcs_type.data();
        ComputeTurbulentViscosity(dt_advance, xvel_old, yvel_old,Tau[level],
                                  state_old[IntVars::cons],
                                  *walldist[level].get(),
                                  *eddyDiffs, *Hfx1, *Hfx2, *Hfx3, *Diss, // to be updated
                                  fine_geom, mapfac[level],
                                  z_phys_nd[level], solverChoice,
                                  m_SurfaceLayer, z_0, l_use_terrain_fitted_coords,
                                  l_use_moisture, level,
                                  bc_ptr_h);
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

    if (solverChoice.custom_w_subsidence) {
        prob->update_w_subsidence(old_time,
                                  h_w_subsid[level], d_w_subsid[level],
                                  fine_geom, z_phys_nd[level]);
    }

    // ***********************************************************************************************
    // Convert old velocity available on faces to old momentum on faces to be used in time integration
    // ***********************************************************************************************
    MultiFab density(state_old[IntVars::cons], make_alias, Rho_comp, 1);

    //
    // This is an optimization since we won't need more than one ghost
    // cell of momentum in the integrator if not using numerical diffusion
    //
    IntVect ngu = (!solverChoice.use_num_diff) ? IntVect(1,1,1) : xvel_old.nGrowVect();
    IntVect ngv = (!solverChoice.use_num_diff) ? IntVect(1,1,1) : yvel_old.nGrowVect();
    IntVect ngw = (!solverChoice.use_num_diff) ? IntVect(1,1,0) : zvel_old.nGrowVect();

    const MultiFab* c_vfrac = nullptr;
    if (solverChoice.terrain_type == TerrainType::EB) {
        c_vfrac = &((get_eb(level).get_const_factory())->getVolFrac());
    }

    VelocityToMomentum(xvel_old, ngu, yvel_old, ngv, zvel_old, ngw, density,
                       state_old[IntVars::xmom],
                       state_old[IntVars::ymom],
                       state_old[IntVars::zmom],
                       domain, domain_bcs_type, c_vfrac);

    MultiFab::Copy(xvel_new,xvel_old,0,0,1,xvel_old.nGrowVect());
    MultiFab::Copy(yvel_new,yvel_old,0,0,1,yvel_old.nGrowVect());
    MultiFab::Copy(zvel_new,zvel_old,0,0,1,zvel_old.nGrowVect());

    bool fast_only = false;
    bool vel_and_mom_synced = true;

    apply_bcs(state_old, old_time,
              state_old[IntVars::cons].nGrow(), state_old[IntVars::xmom].nGrow(),
              fast_only, vel_and_mom_synced);
    cons_to_prim(state_old[IntVars::cons], state_old[IntVars::cons].nGrow());

    // ***********************************************************************************************
    // Define a new MultiFab that holds q_total and fill it by summing the moisture components --
    //      to be used in buoyancy calculation and as part of the inertial weighting in the
    // ***********************************************************************************************

    const bool l_eb_terrain = (solverChoice.terrain_type == TerrainType::EB);
    MultiFab qt(grids[level], dmap[level], 1, (l_eb_terrain) ? 2 : 1);
    qt.setVal(0.0);

#include "ERF_TI_no_substep_fun.H"
#include "ERF_TI_substep_fun.H"
#include "ERF_TI_slow_rhs_pre.H"
#include "ERF_TI_slow_rhs_post.H"

    // ***************************************************************************************
    // Setup the integrator and integrate for a single timestep
    // **************************************************************************************
    MRISplitIntegrator<Vector<MultiFab> >& mri_integrator = *mri_integrator_mem[level];

    // Define rhs and 'post update' utility function that is called after calculating
    // any state data (e.g. at RK stages or at the end of a timestep)
    mri_integrator.set_slow_rhs_pre(slow_rhs_fun_pre);
    mri_integrator.set_slow_rhs_post(slow_rhs_fun_post);

    mri_integrator.set_acoustic_substepping(acoustic_substepping_fun);
    mri_integrator.set_slow_fast_timestep_ratio(fixed_mri_dt_ratio > 0 ? fixed_mri_dt_ratio : dt_mri_ratio[level]);
    mri_integrator.set_no_substep(no_substep_fun);

    mri_integrator.advance(state_old, state_new, old_time, dt_advance);

    if (verbose) Print() << "Done with advance_dycore at level " << level << std::endl;
}
