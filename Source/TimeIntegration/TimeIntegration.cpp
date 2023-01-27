#include <AMReX_BC_TYPES.H>
#include <AMReX_TimeIntegrator.H>
#include <ERF_MRI.H>
#include <EddyViscosity.H>
#include <TerrainMetrics.H>
#include <TimeIntegration.H>
#include <PlaneAverage.H>
#include <Diffusion.H>
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
                      MultiFab& source, MultiFab& buoyancy,
#if defined(ERF_USE_MOISTURE)
                      const MultiFab& qvapor, const MultiFab& qcloud, const MultiFab& qice,
#endif
                      const amrex::Geometry fine_geom,
                      const amrex::Real dt_advance, const amrex::Real old_time,
                      amrex::InterpFaceRegister* ifr)
{
    BL_PROFILE_VAR("erf_advance()",erf_advance);
    if (verbose) amrex::Print() << "Starting advance at level " << level << std::endl;

    int nvars = cons_old.nComp();

    MultiFab r_hse (base_state[level], make_alias, 0, 1); // r_0 is first  component
    MultiFab p_hse (base_state[level], make_alias, 1, 1); // p_0 is second component
    MultiFab pi_hse(base_state[level], make_alias, 2, 1); // pi_0 is second component

    MultiFab* r0  = &r_hse;
    MultiFab* p0  = &p_hse;
    MultiFab* pi0 = &pi_hse;

    Real* dptr_rayleigh_tau      = solverChoice.use_rayleigh_damping ? d_rayleigh_tau[level].data() : nullptr;
    Real* dptr_rayleigh_ubar     = solverChoice.use_rayleigh_damping ? d_rayleigh_ubar[level].data() : nullptr;
    Real* dptr_rayleigh_vbar     = solverChoice.use_rayleigh_damping ? d_rayleigh_vbar[level].data() : nullptr;
    amrex::Real* dptr_rayleigh_thetabar = solverChoice.use_rayleigh_damping ? d_rayleigh_thetabar[level].data() : nullptr;

    bool l_use_terrain = solverChoice.use_terrain;
    bool l_use_diff    = ( (solverChoice.molec_diff_type != MolecDiffType::None) ||
                           (solverChoice.les_type        !=       LESType::None) ||
                           (solverChoice.pbl_type        !=       PBLType::None) );
    bool l_use_kturb   = ( (solverChoice.les_type != LESType::None)   ||
                           (solverChoice.pbl_type != PBLType::None) );

    const BoxArray& ba            = cons_old.boxArray();
    const BoxArray& ba_z          = zvel_old.boxArray();
    const DistributionMapping& dm = cons_old.DistributionMap();

    MultiFab    S_prim  (ba  , dm, NUM_PRIM,          cons_old.nGrowVect());
    MultiFab  pi_stage  (ba  , dm,        1,          cons_old.nGrowVect());
    MultiFab fast_coeffs(ba_z, dm,        5,          0);
    MultiFab* eddyDiffs;
    MultiFab* SmnSmn;
    if (l_use_kturb) {
      eddyDiffs = new MultiFab(ba, dm, EddyDiff::NumDiffs, 1);
      if(solverChoice.les_type == LESType::Deardorff) SmnSmn = new MultiFab(ba, dm, 1, 0);
    } else {
      eddyDiffs = nullptr;
      SmnSmn    = nullptr;
    }

    // **************************************************************************************
    // Compute strain for use in slow RHS, Smagorinsky model, and MOST
    // **************************************************************************************
    BoxArray ba12 = convert(ba, IntVect(1,1,0));
    BoxArray ba13 = convert(ba, IntVect(1,0,1));
    BoxArray ba23 = convert(ba, IntVect(0,1,1));
    MultiFab* Tau11 = nullptr;
    MultiFab* Tau22 = nullptr;
    MultiFab* Tau33 = nullptr;
    MultiFab* Tau12 = nullptr;
    MultiFab* Tau13 = nullptr;
    MultiFab* Tau23 = nullptr;
    MultiFab* Tau21 = nullptr;
    MultiFab* Tau31 = nullptr;
    MultiFab* Tau32 = nullptr;
    {
    BL_PROFILE("erf_advance_strain");
    if (l_use_diff) {
        Tau11 = new MultiFab(ba  , dm, 1, IntVect(1,1,0));
        Tau22 = new MultiFab(ba  , dm, 1, IntVect(1,1,0));
        Tau33 = new MultiFab(ba  , dm, 1, IntVect(1,1,0));
        Tau12 = new MultiFab(ba12, dm, 1, IntVect(1,1,0));
        Tau13 = new MultiFab(ba13, dm, 1, IntVect(1,1,0));
        Tau23 = new MultiFab(ba23, dm, 1, IntVect(1,1,0));
        if (l_use_terrain) {
            Tau21 = new MultiFab(ba12, dm, 1, IntVect(1,1,0));
            Tau31 = new MultiFab(ba13, dm, 1, IntVect(1,1,0));
            Tau32 = new MultiFab(ba23, dm, 1, IntVect(1,1,0));
        }

        const amrex::BCRec* bc_ptr_h = domain_bcs_type.data();
        const GpuArray<Real, AMREX_SPACEDIM> dxInv = fine_geom.InvCellSizeArray();

        for ( MFIter mfi(cons_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box bx    = mfi.tilebox();
            Box bxcc  = mfi.growntilebox(IntVect(1,1,0));
            Box tbxxy = bx; tbxxy.convert(IntVect(1,1,0));
            Box tbxxz = bx; tbxxz.convert(IntVect(1,0,1));
            Box tbxyz = bx; tbxyz.convert(IntVect(0,1,1));

            // Fill strain ghost cells for building K_turb
            tbxxy.growLo(0,1);tbxxy.growLo(1,1);
            tbxxz.growLo(0,1);tbxxz.growLo(1,1);
            tbxyz.growLo(0,1);tbxyz.growLo(1,1);
            tbxxy.growHi(0,1);tbxxy.growHi(1,1);
            tbxxz.growHi(0,1);tbxxz.growHi(1,1);
            tbxyz.growHi(0,1);tbxyz.growHi(1,1);

            const Array4<const Real> & u = xvel_old.array(mfi);
            const Array4<const Real> & v = yvel_old.array(mfi);
            const Array4<const Real> & w = zvel_old.array(mfi);

            Array4<Real> tau11 = Tau11->array(mfi);
            Array4<Real> tau22 = Tau22->array(mfi);
            Array4<Real> tau33 = Tau33->array(mfi);
            Array4<Real> tau12 = Tau12->array(mfi);
            Array4<Real> tau13 = Tau13->array(mfi);
            Array4<Real> tau23 = Tau23->array(mfi);

            Array4<Real> tau21  = l_use_terrain ? Tau21->array(mfi) : Array4<Real>{};
            Array4<Real> tau31  = l_use_terrain ? Tau31->array(mfi) : Array4<Real>{};
            Array4<Real> tau32  = l_use_terrain ? Tau32->array(mfi) : Array4<Real>{};
            const Array4<const Real>& z_nd = l_use_terrain ? z_phys_nd[level]->const_array(mfi) : Array4<const Real>{};

            const Array4<const Real> mf_m = mapfac_m[level]->array(mfi);
            const Array4<const Real> mf_u = mapfac_u[level]->array(mfi);
            const Array4<const Real> mf_v = mapfac_v[level]->array(mfi);

            if (l_use_terrain) {
                ComputeStrain_T(bxcc, tbxxy, tbxxz, tbxyz,
                                u, v, w,
                                tau11, tau22, tau33,
                                tau12, tau13,
                                tau21, tau23,
                                tau31, tau32,
                                z_nd, bc_ptr_h, dxInv,
                                mf_m, mf_u, mf_v);
            } else {
                ComputeStrain_N(bxcc, tbxxy, tbxxz, tbxyz,
                                u, v, w,
                                tau11, tau22, tau33,
                                tau12, tau13, tau23,
                                bc_ptr_h, dxInv,
                                mf_m, mf_u, mf_v);
            }
        } // mfi
    } // l_use_diff
    } // profile


    MultiFab Omega (zmom_old.boxArray(),dm,1,1);

#include "TI_utils.H"

    amrex::Vector<amrex::MultiFab> state_old;
    amrex::Vector<amrex::MultiFab> state_new;

    // **************************************************************************************
    // Here we define state_old and state_new which are to be advanced
    // **************************************************************************************
    {
    BL_PROFILE("erf_advance_part_1");
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
    } // profile

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
        ComputeTurbulentViscosity(xvel_old, yvel_old,
                                  *Tau11, *Tau22, *Tau33,
                                  *Tau12, *Tau13, *Tau23,
                                  state_old[IntVar::cons],
                                  *eddyDiffs, fine_geom, *mapfac_u[level], *mapfac_v[level],
                                  solverChoice, m_most);
    }

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

    MultiFab::Copy(xvel_new,xvel_old,0,0,1,xvel_old.nGrowVect());
    MultiFab::Copy(yvel_new,yvel_old,0,0,1,yvel_old.nGrowVect());
    MultiFab::Copy(zvel_new,zvel_old,0,0,1,zvel_old.nGrowVect());

    bool fast_only          = false;
    bool vel_and_mom_synced = true;
    apply_bcs(state_old, old_time,
              state_old[IntVar::cons].nGrow(), state_old[IntVar::xmom].nGrow(), fast_only,
              vel_and_mom_synced);
    cons_to_prim(state_old[IntVar::cons], state_old[IntVar::cons].nGrow());
    } // profile

#ifdef ERF_USE_MOISTURE
    PlaneAverage qv_ave(&qvapor, geom[level], solverChoice.ave_plane);
    PlaneAverage qc_ave(&qcloud, geom[level], solverChoice.ave_plane);
    PlaneAverage qi_ave(&qice, geom[level], solverChoice.ave_plane);

    // Compute plane averages
    qv_ave.compute_averages(ZDir(), qv_ave.field());
    qc_ave.compute_averages(ZDir(), qc_ave.field());
    qi_ave.compute_averages(ZDir(), qi_ave.field());

    // get plane averaged data
    int ncell = qv_ave.ncell_line();

    Gpu::HostVector  <Real> qv_h(ncell), qi_h(ncell), qc_h(ncell);
    Gpu::DeviceVector<Real> qv_d(ncell), qc_d(ncell), qi_d(ncell);

    // Fill the vectors with the line averages computed above
    qv_ave.line_average(0, qv_h);
    qi_ave.line_average(0, qi_h);
    qc_ave.line_average(0, qc_h);


    // Copy data to device
    Gpu::copyAsync(Gpu::hostToDevice, qv_h.begin(), qv_h.end(), qv_d.begin());
    Gpu::copyAsync(Gpu::hostToDevice, qi_h.begin(), qi_h.end(), qi_d.begin());
    Gpu::copyAsync(Gpu::hostToDevice, qc_h.begin(), qc_h.end(), qc_d.begin());
    Gpu::streamSynchronize();
#endif

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
    } // profile

    mri_integrator.advance(state_old, state_new, old_time, dt_advance);

    if (l_use_kturb) delete eddyDiffs;

    if(solverChoice.les_type == LESType::Deardorff) delete SmnSmn;

    if (l_use_diff) {
      delete Tau11;
      delete Tau22;
      delete Tau33;
      delete Tau12;
      delete Tau13;
      delete Tau23;
      if (l_use_terrain) {
        delete Tau21;
        delete Tau31;
        delete Tau32;
      }
    }

    if (verbose) Print() << "Done with advance at level " << level << std::endl;
}
