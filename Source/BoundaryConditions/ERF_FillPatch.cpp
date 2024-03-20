#include <ERF.H>
#include <ERF_PhysBCFunct.H>
#include <IndexDefines.H>
#include <TimeInterpolatedData.H>
#include <ERF_FillPatcher.H>
#include <Utils.H>

using namespace amrex;

PhysBCFunctNoOp null_bc;

/*
 * Fill valid and ghost data with the "state data" at the given time
 * NOTE: THIS OPERATES ON VELOCITY (MOMENTA ARE JUST TEMPORARIES)
 *
 * @param[in] lev  level of refinement at which to fill the data
 * @param[in] time time at which the data should be filled
 * @param[out] mfs_vel Vector of MultiFabs to be filled containing, in order: cons, xvel, yvel, and zvel
 * @param[out] mfs_mom Vector of MultiFabs to be filled containing, in order: cons, xmom, ymom, and zmom
 */
void
ERF::FillPatch (int lev, Real time,
                const Vector<MultiFab*>& mfs_vel,     // This includes cc quantities and VELOCITIES
                const Vector<MultiFab*>& mfs_mom,     // This includes cc quantities and MOMENTA
                bool fillset, bool cons_only)
{
    BL_PROFILE_VAR("ERF::FillPatch()",ERF_FillPatch);
    Interpolater* mapper = nullptr;

    //
    // ***************************************************************************
    // The first thing we do is interpolate the momenta on the "valid" faces of
    // the fine grids (where the interface is coarse/fine not fine/fine) -- this
    // will not be over-written below because the FillPatch operators see these as
    // valid faces.
    //
    // Note that we interpolate momentum not velocity, but all the other boundary
    // conditions are imposed on velocity, so we convert to momentum here then
    // convert back.
    // ***************************************************************************
    if (lev>0 && fillset) {
        if (cf_set_width > 0) {
            FPr_c[lev-1].FillSet(*mfs_vel[Vars::cons], time, null_bc, domain_bcs_type);
        }
        if (cf_set_width >= 0 && !cons_only) {
            VelocityToMomentum(*mfs_vel[Vars::xvel], mfs_vel[Vars::xvel]->nGrowVect(),
                               *mfs_vel[Vars::yvel], mfs_vel[Vars::yvel]->nGrowVect(),
                               *mfs_vel[Vars::zvel], mfs_vel[Vars::zvel]->nGrowVect(),
                               *mfs_vel[Vars::cons],
                               *mfs_mom[IntVars::xmom],
                               *mfs_mom[IntVars::ymom],
                               *mfs_mom[IntVars::zmom],
                               Geom(lev).Domain(),
                               domain_bcs_type,
                               solverChoice.use_NumDiff);
            FPr_u[lev-1].FillSet(*mfs_mom[IntVars::xmom], time, null_bc, domain_bcs_type);
            FPr_v[lev-1].FillSet(*mfs_mom[IntVars::ymom], time, null_bc, domain_bcs_type);
            FPr_w[lev-1].FillSet(*mfs_mom[IntVars::zmom], time, null_bc, domain_bcs_type);

            MomentumToVelocity(*mfs_vel[Vars::xvel], *mfs_vel[Vars::yvel], *mfs_vel[Vars::zvel],
                               *mfs_vel[Vars::cons],
                               *mfs_mom[IntVars::xmom],
                               *mfs_mom[IntVars::ymom],
                               *mfs_mom[IntVars::zmom],
                               Geom(lev).Domain(),
                               domain_bcs_type);
        }
    }

    IntVect ngvect_cons = mfs_vel[Vars::cons]->nGrowVect();
    IntVect ngvect_vels = mfs_vel[Vars::xvel]->nGrowVect();

    if (lev == 0)
    {
        const int icomp = 0;

        Vector<Real> ftime    = {t_old[lev], t_new[lev]};

        Vector<MultiFab*> fmf = {&vars_old[lev][Vars::cons], &vars_new[lev][Vars::cons]};
        const int  ncomp = mfs_vel[Vars::cons]->nComp();
        FillPatchSingleLevel(*mfs_vel[Vars::cons], time, fmf, ftime, icomp, icomp, ncomp,
                             geom[lev], *physbcs_cons[lev], BCVars::cons_bc+icomp);

        if (!cons_only) {
            fmf = {&vars_old[lev][Vars::xvel], &vars_new[lev][Vars::xvel]};
            const int ncomp_u = 1;
            FillPatchSingleLevel(*mfs_vel[Vars::xvel], time, fmf, ftime, icomp, icomp, ncomp_u,
                                 geom[lev], *physbcs_u[lev], BCVars::xvel_bc);

            fmf = {&vars_old[lev][Vars::yvel], &vars_new[lev][Vars::yvel]};
            const int ncomp_v = 1;
            FillPatchSingleLevel(*mfs_vel[Vars::yvel], time, fmf, ftime, icomp, icomp, ncomp_v,
                                 geom[lev], *physbcs_v[lev], BCVars::yvel_bc);

            fmf = {&vars_old[lev][Vars::zvel], &vars_new[lev][Vars::zvel]};
            const int ncomp_w = 1;
            FillPatchSingleLevel(*mfs_vel[Vars::zvel], time, fmf, ftime, icomp, icomp, ncomp_w,
                                 geom[lev], *physbcs_w_no_terrain[lev], BCVars::zvel_bc);
            (*physbcs_w[lev])(*mfs_vel[Vars::zvel],*mfs_vel[Vars::xvel],*mfs_vel[Vars::yvel],
                               ngvect_vels,time,BCVars::xvel_bc,BCVars::yvel_bc,BCVars::zvel_bc);
        } // !cons_only

    } else {

        Vector<Real> ftime    = {t_old[lev], t_new[lev]};
        Vector<Real> ctime    = {t_old[lev-1], t_new[lev-1]};

        Vector<MultiFab*> fmf = {&vars_old[lev  ][Vars::cons], &vars_new[lev  ][Vars::cons]};
        Vector<MultiFab*> cmf = {&vars_old[lev-1][Vars::cons], &vars_new[lev-1][Vars::cons]};
        MultiFab& mf_c = *mfs_vel[Vars::cons];
        mapper = &cell_cons_interp;
        FillPatchTwoLevels(mf_c, time, cmf, ctime, fmf, ftime,
                           0, 0, mf_c.nComp(), geom[lev-1], geom[lev],
                           null_bc, BCVars::cons_bc, null_bc, BCVars::cons_bc, refRatio(lev-1),
                           mapper, domain_bcs_type, BCVars::cons_bc);

        if (!cons_only) {
            mapper = &face_cons_linear_interp;

            MultiFab& mf_u = *mfs_vel[Vars::xvel];
            fmf = {&vars_old[lev  ][Vars::xvel], &vars_new[lev  ][Vars::xvel]};
            cmf = {&vars_old[lev-1][Vars::xvel], &vars_new[lev-1][Vars::xvel]};
            FillPatchTwoLevels(mf_u, time, cmf, ctime, fmf, ftime,
                               0, 0, 1, geom[lev-1], geom[lev],
                               null_bc, BCVars::xvel_bc, null_bc, BCVars::xvel_bc, refRatio(lev-1),
                               mapper, domain_bcs_type, BCVars::xvel_bc);

            MultiFab& mf_v = *mfs_vel[Vars::yvel];
            fmf = {&vars_old[lev  ][Vars::yvel], &vars_new[lev  ][Vars::yvel]};
            cmf = {&vars_old[lev-1][Vars::yvel], &vars_new[lev-1][Vars::yvel]};
            FillPatchTwoLevels(mf_v, time, cmf, ctime, fmf, ftime,
                               0, 0, 1, geom[lev-1], geom[lev],
                               null_bc, BCVars::yvel_bc, null_bc, BCVars::yvel_bc, refRatio(lev-1),
                               mapper, domain_bcs_type, BCVars::yvel_bc);

            MultiFab& mf_w = *mfs_vel[Vars::zvel];
            fmf = {&vars_old[lev  ][Vars::zvel], &vars_new[lev  ][Vars::zvel]};
            cmf = {&vars_old[lev-1][Vars::zvel], &vars_new[lev-1][Vars::zvel]};
            FillPatchTwoLevels(mf_w, time, cmf, ctime, fmf, ftime,
                               0, 0, 1, geom[lev-1], geom[lev],
                               null_bc, BCVars::zvel_bc, null_bc, BCVars::zvel_bc, refRatio(lev-1),
                               mapper, domain_bcs_type, BCVars::zvel_bc);
        } // !cons_only
    } // lev > 0

    // ***************************************************************************
    // Physical bc's at domain boundary
    // ***************************************************************************
    int icomp_cons = 0;
    int ncomp_cons = mfs_vel[Vars::cons]->nComp();

#ifdef ERF_USE_NETCDF
    // We call this here because it is an ERF routine
    if (use_real_bcs && (lev==0)) {
        fill_from_realbdy(mfs_vel,time,false,0,ncomp_cons);
    }
#endif

    if (m_r2d) fill_from_bndryregs(mfs_vel,time);

    // We call these even if init_type == real because these will fill the vertical bcs
    (*physbcs_cons[lev])(*mfs_vel[Vars::cons],icomp_cons,ncomp_cons,ngvect_cons,time,BCVars::cons_bc);
    if (!cons_only) {
        (*physbcs_u[lev])(*mfs_vel[Vars::xvel],0,1,ngvect_vels,time,BCVars::xvel_bc);
        (*physbcs_v[lev])(*mfs_vel[Vars::yvel],0,1,ngvect_vels,time,BCVars::yvel_bc);
        (*physbcs_w[lev])(*mfs_vel[Vars::zvel],*mfs_vel[Vars::xvel],*mfs_vel[Vars::yvel],
                           ngvect_vels,time,BCVars::xvel_bc,BCVars::yvel_bc,BCVars::zvel_bc);
    }
}

/*
 * Fill ghost cells of qmoist
 *
 * @param[in] lev  level of refinement at which to fill the data
 * @param[in] time time at which the data should be filled
 * @param[out] mf MultiFab to be filled (qmoist[lev])
 */
void
ERF::FillPatchMoistVars (int lev, MultiFab& mf)
{
    BL_PROFILE_VAR("ERF::FillPatchMoistVars()",ERF_FillPatchMoistVars);
    // ***************************************************************************
    // Physical bc's at domain boundary
    // ***************************************************************************
    int icomp_cons = 0;
    int ncomp_cons = 1; // We only fill qv, the first component

    // Note that we are filling qv, stored in qmoist[lev], with the input data (if there is any), stored
    // in RhoQ1_comp.

    if (!use_real_bcs) {
        Real time = Real(0.0);
        IntVect ngvect_cons = mf.nGrowVect();
        int bccomp_cons = BCVars::RhoQ1_bc_comp;

        (*physbcs_cons[lev])(mf,icomp_cons,ncomp_cons,ngvect_cons,time,bccomp_cons);
    }

    mf.FillBoundary(geom[lev].periodicity());
}

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
    int bccomp;
    Interpolater* mapper;

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

    AMREX_ALWAYS_ASSERT(mfs_mom.size() == IntVars::NumTypes);
    AMREX_ALWAYS_ASSERT(mfs_vel.size() == Vars::NumTypes);

    // Impose poor person's immersed boundary
    if (xflux_mask[lev]) {
        amrex::Multiply(*mfs_mom[IntVars::xmom], *xflux_mask[lev], 0,0,1,0);
    }
    if (yflux_mask[lev]) {
        amrex::Multiply(*mfs_mom[IntVars::ymom], *yflux_mask[lev], 0,0,1,0);
    }
    if (zflux_mask[lev]) {
        amrex::Multiply(*mfs_mom[IntVars::zmom], *zflux_mask[lev], 0,0,1,0);
    }

    // We always come in to this call with updated momenta but we need to create updated velocity
    //    in order to impose the rest of the bc's
    if (!cons_only) {
        // This only fills VALID region of velocity
        MomentumToVelocity(*mfs_vel[Vars::xvel], *mfs_vel[Vars::yvel], *mfs_vel[Vars::zvel],
                           *mfs_vel[Vars::cons],
                           *mfs_mom[IntVars::xmom], *mfs_mom[IntVars::ymom], *mfs_mom[IntVars::zmom],
                            Geom(lev).Domain(), domain_bcs_type);
    }

    // We now start working on VELOCITY
    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx)
    {
        if (cons_only && var_idx != Vars::cons) continue;

        MultiFab& mf = *mfs_vel[var_idx];

        IntVect ngvect;
        int icomp, ncomp;
        if (var_idx == Vars::cons)
        {
            bccomp = icomp_cons;
            mapper = &cell_cons_interp;
            ngvect = IntVect(ng_cons,ng_cons,ng_cons);
            icomp  = icomp_cons;
            ncomp  = ncomp_cons;
        }
        else if (var_idx == IntVars::xmom)
        {
            bccomp = BCVars::xvel_bc;
            mapper = &face_cons_linear_interp;
            ngvect = IntVect(ng_vel,ng_vel,ng_vel);
            icomp  = 0;
            ncomp  = 1;
        }
        else if (var_idx == IntVars::ymom)
        {
            bccomp = BCVars::yvel_bc;
            mapper = &face_cons_linear_interp;
            ngvect = IntVect(ng_vel,ng_vel,ng_vel);
            icomp  = 0;
            ncomp  = 1;
        }
        else if (var_idx == IntVars::zmom)
        {
            bccomp = BCVars::zvel_bc;
            mapper = &face_cons_linear_interp;
            ngvect = IntVect(ng_vel,ng_vel,0);
            icomp  = 0;
            ncomp  = 1;
        }

        if (lev == 0)
        {
            // This fills fine-fine ghost values of VELOCITY
            mf.FillBoundary(icomp,ncomp,ngvect,geom[lev].periodicity());
        }
        else
        {
            //
            // NOTE: All interpolation here happens on velocities not momenta;
            //       note we only do the interpolation and FillBoundary here,
            //       physical bc's are imposed later
            //
            // NOTE: This will only fill velocity from coarse grid *outside* the fine grids
            //       unlike the FillSet calls above which filled momenta on the coarse/fine bdy
            //
            Vector<MultiFab*> fmf = {&mf};
            Vector<MultiFab*> cmf = {&vars_old[lev-1][var_idx], &vars_new[lev-1][var_idx]};
            Vector<Real> ctime    = {t_old[lev-1], t_new[lev-1]};

            FillPatchTwoLevels(mf, time, cmf, ctime, fmf, {time},
                               icomp, icomp, ncomp, geom[lev-1], geom[lev],
                               null_bc, 0, null_bc, 0, refRatio(lev-1),
                               mapper, domain_bcs_type, bccomp);
        } // lev > 0
    } // var_idx

    // ***************************************************************************
    // Physical bc's at domain boundary
    // ***************************************************************************
    IntVect ngvect_cons = IntVect(ng_cons,ng_cons,ng_cons);
    IntVect ngvect_vels = IntVect(ng_vel ,ng_vel ,ng_vel);

#ifdef ERF_USE_NETCDF
    // We call this here because it is an ERF routine
    if (use_real_bcs && (lev==0)) {
        fill_from_realbdy(mfs_vel,time,false,0,ncomp_cons);
    }
#endif

    if (m_r2d) fill_from_bndryregs(mfs_vel,time);

    // We call this even if init_type == real because this routine will fill the vertical bcs
    (*physbcs_cons[lev])(*mfs_vel[Vars::cons],icomp_cons,ncomp_cons,ngvect_cons,time,BCVars::cons_bc);
    if (!cons_only) {
        (*physbcs_u[lev])(*mfs_vel[Vars::xvel],0,1,ngvect_vels,time,BCVars::xvel_bc);
        (*physbcs_v[lev])(*mfs_vel[Vars::yvel],0,1,ngvect_vels,time,BCVars::yvel_bc);
        (*physbcs_w[lev])(*mfs_vel[Vars::zvel],*mfs_vel[Vars::xvel],*mfs_vel[Vars::yvel],
                           ngvect_vels,time,BCVars::xvel_bc,BCVars::yvel_bc,BCVars::zvel_bc);
    }
    // ***************************************************************************

    // MOST boundary conditions
    if (!(cons_only && ncomp_cons == 1) && m_most && allow_most_bcs) {
        m_most->impose_most_bcs(lev,mfs_vel,
#ifdef ERF_EXPLICIT_MOST_STRESS
                                Tau13_lev[lev].get(), Tau31_lev[lev].get(),
                                Tau23_lev[lev].get(), Tau32_lev[lev].get(),
                                SFS_hfx3_lev[lev].get(),
#else
                                eddyDiffs_lev[lev].get(),
#endif
                                z_phys_nd[lev].get());
    }

    // We always come in to this call with momenta so we need to leave with momenta!
    // We need to make sure we convert back on all ghost cells/faces because this is
    // how velocity from fine-fine copies (as well as physical and interpolated bcs) will be filled
    if (!cons_only) {
        VelocityToMomentum(*mfs_vel[Vars::xvel], mfs_vel[Vars::xvel]->nGrowVect(),
                           *mfs_vel[Vars::yvel], mfs_vel[Vars::yvel]->nGrowVect(),
                           *mfs_vel[Vars::zvel], mfs_vel[Vars::zvel]->nGrowVect(),
                           *mfs_vel[Vars::cons],
                           *mfs_mom[IntVars::xmom], *mfs_mom[IntVars::ymom], *mfs_mom[IntVars::zmom],
                           Geom(lev).Domain(),
                           domain_bcs_type,
                           solverChoice.use_NumDiff);
    }
}

/*
 * Fill valid and ghost data.
 * This version fills an entire MultiFab by interpolating from the coarser level -- this is used
 * only when a new level of refinement is being created during a run (i.e not at initialization)
 * This will never be used with static refinement.
 *
 * @param[in]  lev            level of refinement at which to fill the data
 * @param[in]  time           time at which the data should be filled
 * @param[out] mfs            Vector of MultiFabs to be filled containing, in order: cons, xvel, yvel, and zvel data
 */
void
ERF::FillCoarsePatch (int lev, Real time)
{
    BL_PROFILE_VAR("FillCoarsePatch()",FillCoarsePatch);
    AMREX_ASSERT(lev > 0);

    int icomp = 0;

    int bccomp = 0;
    Interpolater* mapper = &cell_cons_interp;
    InterpFromCoarseLevel(vars_new[lev][Vars::cons], time, vars_new[lev-1][Vars::cons],
                          icomp, icomp, vars_new[lev][Vars::cons].nComp(),
                          geom[lev-1], geom[lev],
                          null_bc, 0, null_bc, 0, refRatio(lev-1),
                          mapper, domain_bcs_type, bccomp);


    mapper = &face_cons_linear_interp;

    for (int which_lev = lev-1; which_lev <= lev; which_lev++)
    {
        // First fill density so we can convert velocity to momenta
        // Note we can only do this because we first interpolated density above
        bool cons_only = true;
        FillPatch(which_lev, time, {&vars_new[which_lev][Vars::cons], &vars_new[which_lev][Vars::xvel],
                                    &vars_new[which_lev][Vars::yvel], &vars_new[which_lev][Vars::zvel]},
                                   {&vars_new[which_lev][Vars::cons],
                                    &rU_new[which_lev], &rV_new[which_lev], &rW_new[which_lev]},
                                    false, cons_only);

        VelocityToMomentum(vars_new[which_lev][Vars::xvel], IntVect(0,0,0),
                           vars_new[which_lev][Vars::yvel], IntVect(0,0,0),
                           vars_new[which_lev][Vars::zvel], IntVect(0,0,0),
                           vars_new[which_lev][Vars::cons],
                             rU_new[which_lev],
                             rV_new[which_lev],
                             rW_new[which_lev],
                           Geom(lev).Domain(),
                           domain_bcs_type,
                           true);
    }

    bccomp = BCVars::xvel_bc;
    InterpFromCoarseLevel(rU_new[lev], time, rU_new[lev-1],
                          0, 0, 1, geom[lev-1], geom[lev],
                          null_bc, 0, null_bc, 0, refRatio(lev-1),
                          mapper, domain_bcs_type, bccomp);

    bccomp = BCVars::yvel_bc;
    InterpFromCoarseLevel(rV_new[lev], time, rV_new[lev-1],
                          0, 0, 1, geom[lev-1], geom[lev],
                          null_bc, 0, null_bc, 0, refRatio(lev-1),
                          mapper, domain_bcs_type, bccomp);

    bccomp = BCVars::zvel_bc;
    InterpFromCoarseLevel(rW_new[lev], time, rW_new[lev-1],
                          0, 0, 1, geom[lev-1], geom[lev],
                          null_bc, 0, null_bc, 0, refRatio(lev-1),
                          mapper, domain_bcs_type, bccomp);

   for (int which_lev = lev-1; which_lev <= lev; which_lev++)
   {
        MomentumToVelocity(vars_new[which_lev][Vars::xvel],
                           vars_new[which_lev][Vars::yvel],
                           vars_new[which_lev][Vars::zvel],
                           vars_new[which_lev][Vars::cons],
                             rU_new[which_lev],
                             rV_new[which_lev],
                             rW_new[which_lev],
                           Geom(lev).Domain(),
                           domain_bcs_type);
    }

    vars_new[lev][Vars::cons].FillBoundary(geom[lev].periodicity());
    vars_new[lev][Vars::xvel].FillBoundary(geom[lev].periodicity());
    vars_new[lev][Vars::yvel].FillBoundary(geom[lev].periodicity());
    vars_new[lev][Vars::zvel].FillBoundary(geom[lev].periodicity());

    // ***************************************************************************
    // Physical bc's at domain boundary
    // ***************************************************************************
    IntVect ngvect_cons = vars_new[lev][Vars::cons].nGrowVect();
    IntVect ngvect_vels = vars_new[lev][Vars::xvel].nGrowVect();

    int ncomp_cons = vars_new[lev][Vars::cons].nComp();

    (*physbcs_cons[lev])(vars_new[lev][Vars::cons],0,ncomp_cons,ngvect_cons,time,BCVars::cons_bc);
    (   *physbcs_u[lev])(vars_new[lev][Vars::xvel],0,1         ,ngvect_vels,time,BCVars::xvel_bc);
    (   *physbcs_v[lev])(vars_new[lev][Vars::yvel],0,1         ,ngvect_vels,time,BCVars::yvel_bc);
    (   *physbcs_w[lev])(vars_new[lev][Vars::zvel],vars_new[lev][Vars::xvel],vars_new[lev][Vars::yvel],
                         ngvect_vels,time,BCVars::xvel_bc,BCVars::yvel_bc,BCVars::zvel_bc);

    // ***************************************************************************
    // Since lev > 0 here we don't worry about m_r2d or wrfbdy data
    // ***************************************************************************
}
