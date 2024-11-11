#include <ERF.H>
#include <ERF_PhysBCFunct.H>
#include <ERF_IndexDefines.H>
#include <ERF_TimeInterpolatedData.H>
#include <ERF_FillPatcher.H>
#include <ERF_Utils.H>

using namespace amrex;

/*
 * Fill valid and ghost data
 * This version fills mfs in valid regions with the values in "mfs" when it is passed in;
 * it is used only to compute ghost values for intermediate stages of a time integrator.
 *
 * @param[in]  lev            level of refinement at which to fill the data
 * @param[in]  time           time at which the data should be filled
 * @param[out] mfs_vel        Vector of MultiFabs to be filled containing, in order: cons, xvel, yvel, and zvel
 * @param[out] mfs_mom        Vector of MultiFabs to be filled containing, in order: cons, xmom, ymom, and zmom
 * @param[in]  ng_cons        number of ghost cells to be filled for conserved (cell-centered) variables
 * @param[in]  ng_vel         number of ghost cells to be filled for velocity components
 * @param[in]  cons_only      if 1 then only fill conserved variables
 * @param[in]  icomp_cons     starting component for conserved variables
 * @param[in]  ncomp_cons     number of components for conserved variables
 * @param[in]  eddyDiffs      diffusion coefficients for LES turbulence models
 * @param[in]  allow_most_bcs if true then use MOST bcs at the low boundary
 */
void
ERF::FillIntermediatePatch (int lev, Real time,
                            const Vector<MultiFab*>& mfs_vel,     // This includes cc quantities and VELOCITIES
                            const Vector<MultiFab*>& mfs_mom,     // This includes cc quantities and MOMENTA
                            int ng_cons, int ng_vel, bool cons_only,
                            int icomp_cons, int ncomp_cons,
                            bool allow_most_bcs)
{
    BL_PROFILE_VAR("FillIntermediatePatch()",FillIntermediatePatch);
    Interpolater* mapper;

    PhysBCFunctNoOp null_bc;

    //
    // ***************************************************************************
    // The first thing we do is interpolate the momenta on the "valid" faces of
    // the fine grids (where the interface is coarse/fine not fine/fine) -- this
    // will not be over-written by interpolation below because the FillPatch
    // operators see these as valid faces.  But we must have these interpolated
    // values in the fine data before we call FillPatchTwoLevels.
    //
    // Also -- note that we might be filling values by interpolation at physical boundaries
    //         here but that's ok because we will overwrite those values when we impose
    //         the physical bc's below
    // ***************************************************************************
    if (lev>0) {
        if (cf_set_width > 0) {
            // We note that mfs_vel[Vars::cons] and mfs_mom[Vars::cons] are in fact the same pointer
            FPr_c[lev-1].FillSet(*mfs_vel[Vars::cons], time, null_bc, domain_bcs_type);
        }
        if ( !cons_only && (cf_set_width >= 0) ) {
            FPr_u[lev-1].FillSet(*mfs_mom[IntVars::xmom], time, null_bc, domain_bcs_type);
            FPr_v[lev-1].FillSet(*mfs_mom[IntVars::ymom], time, null_bc, domain_bcs_type);
            FPr_w[lev-1].FillSet(*mfs_mom[IntVars::zmom], time, null_bc, domain_bcs_type);
        }
    }

    // amrex::Print() << "LEVEL " << lev << " CONS ONLY   " << cons_only <<
    //                   " ICOMP NCOMP " << icomp_cons << " " << ncomp_cons << " NGHOST " << ng_cons << std::endl;

    if (!cons_only) {
        AMREX_ALWAYS_ASSERT(mfs_mom.size() == IntVars::NumTypes);
        AMREX_ALWAYS_ASSERT(mfs_vel.size() == Vars::NumTypes);
    }

    // Enforce no penetration for thin immersed body
    if (!cons_only) {
        // Enforce no penetration for thin immersed body
        if (xflux_imask[lev]) {
            ApplyMask(*mfs_mom[IntVars::xmom], *xflux_imask[lev]);
        }
        if (yflux_imask[lev]) {
            ApplyMask(*mfs_mom[IntVars::ymom], *yflux_imask[lev]);
        }
        if (zflux_imask[lev]) {
            ApplyMask(*mfs_mom[IntVars::zmom], *zflux_imask[lev]);
        }
    }

    //
    // We now start working on conserved quantities + VELOCITY
    //
    if (lev == 0)
    {
        // We don't do anything here because we will call the physbcs routines below,
        // which calls FillBoundary and fills other domain boundary conditions
        // Physical boundaries will be filled below

        if (!cons_only)
        {
            // ***************************************************************************
            // We always come in to this call with updated momenta but we need to create updated velocity
            //    in order to impose the rest of the bc's
            // ***************************************************************************
            // This only fills VALID region of velocity
            MomentumToVelocity(*mfs_vel[Vars::xvel], *mfs_vel[Vars::yvel], *mfs_vel[Vars::zvel],
                               *mfs_vel[Vars::cons],
                               *mfs_mom[IntVars::xmom], *mfs_mom[IntVars::ymom], *mfs_mom[IntVars::zmom],
                                Geom(lev).Domain(), domain_bcs_type);
        }
    }
    else
    {
        //
        // We must fill a temporary then copy it back so we don't double add/subtract
        //
        MultiFab mf(mfs_vel[Vars::cons]->boxArray(),mfs_vel[Vars::cons]->DistributionMap(),
                    mfs_vel[Vars::cons]->nComp()   ,mfs_vel[Vars::cons]->nGrowVect());
        //
        // Set all components to 1.789e19, then copy just the density from *mfs_vel[Vars::cons]
        //
        mf.setVal(1.789e19);
        MultiFab::Copy(mf,*mfs_vel[Vars::cons],Rho_comp,Rho_comp,1,mf.nGrowVect());

        Vector<MultiFab*> fmf = {mfs_vel[Vars::cons],mfs_vel[Vars::cons]};
        Vector<MultiFab*> cmf = {&vars_old[lev-1][Vars::cons], &vars_new[lev-1][Vars::cons]};
        Vector<Real> ctime    = {t_old[lev-1], t_new[lev-1]};
        Vector<Real> ftime    = {time,time};

        // Impose physical bc's on fine data (note time and 0 are not used)
        bool do_fb = true; bool do_terrain_adjustment = false;
        (*physbcs_cons[lev])(*mfs_vel[Vars::cons],icomp_cons,ncomp_cons,IntVect{ng_cons},time,BCVars::cons_bc,
                             do_fb, do_terrain_adjustment);

        if ( (icomp_cons+ncomp_cons > 1) && (interpolation_type == StateInterpType::Perturbational) )
        {
            // Divide (rho theta) by rho to get theta
            MultiFab::Divide(*mfs_vel[Vars::cons],*mfs_vel[Vars::cons],Rho_comp,RhoTheta_comp,1,IntVect{0});

            // Subtract theta_0 from theta
            MultiFab::Subtract(*mfs_vel[Vars::cons],base_state[lev],BaseState::th0_comp,RhoTheta_comp,1,IntVect{0});

            if (!amrex::almostEqual(time,ctime[1])) {
                MultiFab::Divide(vars_old[lev-1][Vars::cons], vars_old[lev-1][Vars::cons],
                                 Rho_comp,RhoTheta_comp,1,vars_old[lev-1][Vars::cons].nGrowVect());
                MultiFab::Subtract(vars_old[lev-1][Vars::cons], base_state[lev-1],
                                   BaseState::th0_comp,RhoTheta_comp,1,vars_old[lev-1][Vars::cons].nGrowVect());
            }
            if (!amrex::almostEqual(time,ctime[0])) {
                MultiFab::Divide(vars_new[lev-1][Vars::cons], vars_new[lev-1][Vars::cons],
                                 Rho_comp,RhoTheta_comp,1,vars_new[lev-1][Vars::cons].nGrowVect());
                MultiFab::Subtract(vars_new[lev-1][Vars::cons], base_state[lev-1],
                                   BaseState::th0_comp,RhoTheta_comp,1,vars_new[lev-1][Vars::cons].nGrowVect());
            }
        }

        // Subtract rho_0 from rho before we interpolate -- note we only subtract
        //    on valid region of mf since the ghost cells will be filled below
        if (icomp_cons == 0 && (interpolation_type == StateInterpType::Perturbational))
        {
            MultiFab::Subtract(*mfs_vel[Vars::cons],base_state[lev],BaseState::r0_comp,Rho_comp,1,IntVect{0});

            if (!amrex::almostEqual(time,ctime[1])) {
                MultiFab::Subtract(vars_old[lev-1][Vars::cons], base_state[lev-1],
                                   BaseState::r0_comp,Rho_comp,1,vars_old[lev-1][Vars::cons].nGrowVect());
            }
            if (!amrex::almostEqual(time,ctime[0])) {
                MultiFab::Subtract(vars_new[lev-1][Vars::cons], base_state[lev-1],
                                   BaseState::r0_comp,Rho_comp,1,vars_new[lev-1][Vars::cons].nGrowVect());
            }
        }

        // Call FillPatchTwoLevels which ASSUMES that all ghost cells have already been filled
        mapper = &cell_cons_interp;
        FillPatchTwoLevels(mf, IntVect{ng_cons}, IntVect(0,0,0),
                           time, cmf, ctime, fmf, ftime,
                           icomp_cons, icomp_cons, ncomp_cons, geom[lev-1], geom[lev],
                           refRatio(lev-1), mapper, domain_bcs_type,
                           icomp_cons);

        if (icomp_cons == 0 && (interpolation_type == StateInterpType::Perturbational))
        {
            // Restore the coarse values to what they were
            if (!amrex::almostEqual(time,ctime[1])) {
                MultiFab::Add(vars_old[lev-1][Vars::cons], base_state[lev-1],
                              BaseState::r0_comp,Rho_comp,1,vars_old[lev-1][Vars::cons].nGrowVect());
            }
            if (!amrex::almostEqual(time,ctime[0])) {
                MultiFab::Add(vars_new[lev-1][Vars::cons], base_state[lev-1],
                              BaseState::r0_comp,Rho_comp,1,vars_new[lev-1][Vars::cons].nGrowVect());
            }

            // Set values in the cells outside the domain boundary so that we can do the Add
            //     without worrying about uninitialized values outside the domain -- these
            //     will be filled in the physbcs call
            mf.setDomainBndry(1.234e20,Rho_comp,1,geom[lev]);

            // Add rho_0 back to rho after we interpolate -- on all the valid + ghost region
            MultiFab::Add(mf, base_state[lev],BaseState::r0_comp,Rho_comp,1,IntVect{ng_cons});
        }

        if ( (icomp_cons+ncomp_cons > 1) && (interpolation_type == StateInterpType::Perturbational) )
        {
            // Add theta_0 to theta
            if (!amrex::almostEqual(time,ctime[1])) {
                MultiFab::Add(vars_old[lev-1][Vars::cons], base_state[lev-1],
                              BaseState::th0_comp,RhoTheta_comp,1,vars_old[lev-1][Vars::cons].nGrowVect());
                MultiFab::Multiply(vars_old[lev-1][Vars::cons], vars_old[lev-1][Vars::cons],
                                   Rho_comp,RhoTheta_comp,1,vars_old[lev-1][Vars::cons].nGrowVect());
            }
            if (!amrex::almostEqual(time,ctime[0])) {
                MultiFab::Add(vars_new[lev-1][Vars::cons], base_state[lev-1],
                              BaseState::th0_comp,RhoTheta_comp,1,vars_new[lev-1][Vars::cons].nGrowVect());
                MultiFab::Multiply(vars_new[lev-1][Vars::cons], vars_new[lev-1][Vars::cons],
                                   Rho_comp,RhoTheta_comp,1,vars_new[lev-1][Vars::cons].nGrowVect());
            }

            // Multiply theta by rho to get (rho theta)
            MultiFab::Multiply(*mfs_vel[Vars::cons],*mfs_vel[Vars::cons],Rho_comp,RhoTheta_comp,1,IntVect{0});

            // Add theta_0 to theta
            MultiFab::Add(*mfs_vel[Vars::cons],base_state[lev],BaseState::th0_comp,RhoTheta_comp,1,IntVect{0});

            // Add theta_0 back to theta
            MultiFab::Add(mf,base_state[lev],BaseState::th0_comp,RhoTheta_comp,1,IntVect{ng_cons});

            // Multiply (theta) by rho to get (rho theta)
            MultiFab::Multiply(mf,mf,Rho_comp,RhoTheta_comp,1,IntVect{ng_cons});
        }

        // Make sure to only copy back the components we worked on
        MultiFab::Copy(*mfs_vel[Vars::cons],mf,icomp_cons,icomp_cons,ncomp_cons,IntVect{ng_cons});

        // *****************************************************************************************

        if (!cons_only)
        {
            // ***************************************************************************
            // We always come in to this call with updated momenta but we need to create updated velocity
            //    in order to impose the rest of the bc's
            // ***************************************************************************
            // This only fills VALID region of velocity
            MomentumToVelocity(*mfs_vel[Vars::xvel], *mfs_vel[Vars::yvel], *mfs_vel[Vars::zvel],
                               *mfs_vel[Vars::cons],
                               *mfs_mom[IntVars::xmom], *mfs_mom[IntVars::ymom], *mfs_mom[IntVars::zmom],
                                Geom(lev).Domain(), domain_bcs_type);

            mapper = &face_cons_linear_interp;

            //
            // NOTE: All interpolation here happens on velocities not momenta;
            //       note we only do the interpolation and FillBoundary here,
            //       physical bc's are imposed later
            //
            // NOTE: This will only fill velocity from coarse grid *outside* the fine grids
            //       unlike the FillSet calls above which filled momenta on the coarse/fine bdy
            //

            MultiFab& mfu = *mfs_vel[Vars::xvel];

            fmf = {&mfu,&mfu};
            cmf = {&vars_old[lev-1][Vars::xvel], &vars_new[lev-1][Vars::xvel]};

            // Call FillPatchTwoLevels which ASSUMES that all ghost cells have already been filled
            FillPatchTwoLevels(mfu, IntVect{ng_vel}, IntVect(0,0,0),
                               time, cmf, ctime, fmf, ftime,
                               0, 0, 1, geom[lev-1], geom[lev],
                               refRatio(lev-1), mapper, domain_bcs_type,
                               BCVars::xvel_bc);

            // *****************************************************************************************

            MultiFab& mfv = *mfs_vel[Vars::yvel];

            fmf = {&mfv,&mfv};
            cmf = {&vars_old[lev-1][Vars::yvel], &vars_new[lev-1][Vars::yvel]};

            // Call FillPatchTwoLevels which ASSUMES that all ghost cells have already been filled
            FillPatchTwoLevels(mfv, IntVect{ng_vel}, IntVect(0,0,0),
                               time, cmf, ctime, fmf, ftime,
                               0, 0, 1, geom[lev-1], geom[lev],
                               refRatio(lev-1), mapper, domain_bcs_type,
                               BCVars::yvel_bc);

            // *****************************************************************************************

            MultiFab& mfw = *mfs_vel[Vars::zvel];

            fmf = {&mfw,&mfw};
            cmf = {&vars_old[lev-1][Vars::zvel], &vars_new[lev-1][Vars::zvel]};

            // Call FillPatchTwoLevels which ASSUMES that all ghost cells have already been filled
            FillPatchTwoLevels(mfw, IntVect{ng_vel}, IntVect(0,0,0),
                               time, cmf, ctime, fmf, ftime,
                               0, 0, 1, geom[lev-1], geom[lev],
                               refRatio(lev-1), mapper, domain_bcs_type,
                               BCVars::zvel_bc);
        } // !cons_only
    } // lev > 0

    // ***************************************************************************
    // Physical bc's at domain boundary
    // ***************************************************************************
    IntVect ngvect_cons = IntVect(ng_cons,ng_cons,ng_cons);
    IntVect ngvect_vels = IntVect(ng_vel ,ng_vel ,ng_vel);

    bool do_fb = true;

#ifdef ERF_USE_NETCDF
    // We call this here because it is an ERF routine
    if (use_real_bcs && (lev==0)) {
        fill_from_realbdy(mfs_vel,time,cons_only,icomp_cons,ncomp_cons,ngvect_cons,ngvect_vels);
        do_fb = false;
    }
#endif

    if (m_r2d) fill_from_bndryregs(mfs_vel,time);

    // We call this even if init_type == InitType::Real because this routine will fill the vertical bcs
    (*physbcs_cons[lev])(*mfs_vel[Vars::cons],icomp_cons,ncomp_cons,ngvect_cons,time,BCVars::cons_bc, do_fb);
    if (!cons_only) {
        (*physbcs_u[lev])(*mfs_vel[Vars::xvel],0,1,ngvect_vels,time,BCVars::xvel_bc, do_fb);
        (*physbcs_v[lev])(*mfs_vel[Vars::yvel],0,1,ngvect_vels,time,BCVars::yvel_bc, do_fb);
        (*physbcs_w[lev])(*mfs_vel[Vars::zvel],*mfs_vel[Vars::xvel],*mfs_vel[Vars::yvel],
                          ngvect_vels,time,BCVars::zvel_bc, do_fb);
    }
    // ***************************************************************************

    // MOST boundary conditions
    if (!(cons_only && ncomp_cons == 1) && m_most && allow_most_bcs) {
        m_most->impose_most_bcs(lev,mfs_vel,
                                Tau11_lev[lev].get(),
                                Tau22_lev[lev].get(),
                                Tau33_lev[lev].get(),
                                Tau12_lev[lev].get(), Tau21_lev[lev].get(),
                                Tau13_lev[lev].get(), Tau31_lev[lev].get(),
                                Tau23_lev[lev].get(), Tau32_lev[lev].get(),
                                SFS_hfx1_lev[lev].get(),
                                SFS_hfx2_lev[lev].get(),
                                SFS_hfx3_lev[lev].get(),
                                SFS_q1fx1_lev[lev].get(),
                                SFS_q1fx2_lev[lev].get(),
                                SFS_q1fx3_lev[lev].get(),
                                z_phys_nd[lev].get());
    }

    // We always come in to this call with momenta so we need to leave with momenta!
    // We need to make sure we convert back on all ghost cells/faces because this is
    // how velocity from fine-fine copies (as well as physical and interpolated bcs) will be filled
    if (!cons_only)
    {
        IntVect ngu = (!solverChoice.use_NumDiff) ? IntVect(1,1,1) : mfs_vel[Vars::xvel]->nGrowVect();
        IntVect ngv = (!solverChoice.use_NumDiff) ? IntVect(1,1,1) : mfs_vel[Vars::yvel]->nGrowVect();
        IntVect ngw = (!solverChoice.use_NumDiff) ? IntVect(1,1,0) : mfs_vel[Vars::zvel]->nGrowVect();

        VelocityToMomentum(*mfs_vel[Vars::xvel], ngu,
                           *mfs_vel[Vars::yvel], ngv,
                           *mfs_vel[Vars::zvel], ngw,
                           *mfs_vel[Vars::cons],
                           *mfs_mom[IntVars::xmom], *mfs_mom[IntVars::ymom], *mfs_mom[IntVars::zmom],
                           Geom(lev).Domain(),
                           domain_bcs_type);
    }

    mfs_mom[IntVars::cons]->FillBoundary(geom[lev].periodicity());
    mfs_mom[IntVars::xmom]->FillBoundary(geom[lev].periodicity());
    mfs_mom[IntVars::ymom]->FillBoundary(geom[lev].periodicity());
    mfs_mom[IntVars::zmom]->FillBoundary(geom[lev].periodicity());
}
