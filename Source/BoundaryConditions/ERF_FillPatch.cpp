#include <ERF.H>
#include <ERF_PhysBCFunct.H>
#include <IndexDefines.H>
#include <TimeInterpolatedData.H>

using namespace amrex;

PhysBCFunctNoOp null_bc;

/*
 * Fill valid and ghost data with the "state data" at the given time
 *
 * @param[in] lev  level of refinement at which to fill the data
 * @param[in] time time at which the data should be filled
 * @param[out] mfs Vector of MultiFabs to be filled containing, in order: cons, xvel, yvel, and zvel data
 */

void
ERF::FillPatch (int lev, Real time, const Vector<MultiFab*>& mfs)
{
    BL_PROFILE_VAR("ERF::FillPatch()",ERF_FillPatch);
    int bccomp;
    amrex::Interpolater* mapper = nullptr;

    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx) {
        MultiFab& mf = *mfs[var_idx];
        const int icomp = 0;
        const int ncomp = mf.nComp();

        if (var_idx == Vars::cons)
        {
            bccomp = BCVars::cons_bc + icomp;
            mapper = &cell_cons_interp;
        }
        else if (var_idx == Vars::xvel)
        {
            bccomp = BCVars::xvel_bc;
            mapper = &face_linear_interp;
        }
        else if (var_idx == Vars::yvel)
        {
            bccomp = BCVars::yvel_bc;
            mapper = &face_linear_interp;
        }
        else if (var_idx == Vars::zvel)
        {
            bccomp = BCVars::zvel_bc;
            mapper = &face_linear_interp;
        } else {
          amrex::Abort("Dont recognize this variable type in ERF_Fillpatch");
        }

        if (lev == 0)
        {
            Vector<MultiFab*> fmf = {&vars_old[lev][var_idx], &vars_new[lev][var_idx]};
            Vector<Real> ftime    = {t_old[lev], t_new[lev]};
            amrex::FillPatchSingleLevel(mf, time, fmf, ftime, icomp, icomp, ncomp,
                                        geom[lev], null_bc, bccomp);
        }
        else
        {
            Vector<MultiFab*> fmf = {&vars_old[lev][var_idx], &vars_new[lev][var_idx]};
            Vector<Real> ftime    = {t_old[lev], t_new[lev]};
            Vector<MultiFab*> cmf = {&vars_old[lev-1][var_idx], &vars_new[lev-1][var_idx]};
            Vector<Real> ctime    = {t_old[lev-1], t_new[lev-1]};

            amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                      0, icomp, ncomp, geom[lev-1], geom[lev],
                                      null_bc, bccomp, null_bc, bccomp, refRatio(lev-1),
                                      mapper, domain_bcs_type, bccomp);
        } // lev > 0
    } // var_idx

    // ***************************************************************************
    // Physical bc's at domain boundary
    // ***************************************************************************
    bool cons_only = false;
    int icomp_cons = 0;
    int ncomp_cons = mfs[Vars::cons]->nComp();

    IntVect ngvect_cons = mfs[Vars::cons]->nGrowVect();
    IntVect ngvect_vels = mfs[Vars::xvel]->nGrowVect();

#ifdef ERF_USE_NETCDF
    // We call this here because it is an ERF routine
    if (init_type == "real") fill_from_wrfbdy(mfs,time);
    if (init_type == "metgrid") fill_from_metgrid(mfs,time);
#endif

    if (m_r2d) fill_from_bndryregs(mfs,time);

    // We call this even if init_type == real because this routine will fill the vertical bcs
    (*physbcs[lev])(mfs,icomp_cons,ncomp_cons,ngvect_cons,ngvect_vels,init_type,cons_only,BCVars::cons_bc);
}

/*
 * Fill ghost cells of qmoist
 *
 * @param[in] lev  level of refinement at which to fill the data
 * @param[in] time time at which the data should be filled
 * @param[out] mf MultiFab to be filled (qmoist[lev])
 */

#ifdef ERF_USE_MOISTURE
void
ERF::FillPatchMoistVars (int lev, MultiFab& mf)
{
    BL_PROFILE_VAR("ERF::FillPatchMoistVars()",ERF_FillPatchMoistVars);
    // ***************************************************************************
    // Physical bc's at domain boundary
    // ***************************************************************************
    bool cons_only = true;
    int icomp_cons = 0;
    int ncomp_cons = 1; // We only fill qv, the first component

    // Note that we are filling qv, stored in qmoist[lev], with the input data (if there is any), stored
    // in RhoQt_comp.
    int bccomp_cons = BCVars::RhoQt_bc_comp;

    IntVect ngvect_cons = mf.nGrowVect();
    IntVect ngvect_vels = {0,0,0};

    if ((init_type != "real") and (init_type != "metgrid")) {
        (*physbcs[lev])({&mf},icomp_cons,ncomp_cons,ngvect_cons,ngvect_vels,init_type,cons_only,bccomp_cons);
    }

    mf.FillBoundary(geom[lev].periodicity());
}
#endif

/*
 * Fill valid and ghost data
 * This version fills mfs in valid regions with the values in "mfs" when it is passed in;
 * it is used only to compute ghost values for intermediate stages of a time integrator.
 *
 * @param[in]  lev            level of refinement at which to fill the data
 * @param[in]  time           time at which the data should be filled
 * @param[out] mfs            Vector of MultiFabs to be filled containing, in order: cons, xvel, yvel, and zvel data
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
                            const Vector<MultiFab*>& mfs,
                            int ng_cons, int ng_vel, bool cons_only,
                            int icomp_cons, int ncomp_cons,
                            MultiFab* eddyDiffs,
                            bool allow_most_bcs)
{
    BL_PROFILE_VAR("FillIntermediatePatch()",FillIntermediatePatch);
    int bccomp;
    amrex::Interpolater* mapper;

    // We should always pass cons, xvel, yvel, and zvel (in that order) in the mfs vector
    AMREX_ALWAYS_ASSERT(mfs.size() == Vars::NumTypes);

    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx)
    {
        if (cons_only && var_idx != Vars::cons) continue;

        MultiFab& mf = *mfs[var_idx];

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
        else if (var_idx == Vars::xvel)
        {
            bccomp = BCVars::xvel_bc;
            mapper = &face_linear_interp;
            ngvect = IntVect(ng_vel,ng_vel,ng_vel);
            icomp  = 0;
            ncomp  = 1;
        }
        else if (var_idx == Vars::yvel)
        {
            bccomp = BCVars::yvel_bc;
            mapper = &face_linear_interp;
            ngvect = IntVect(ng_vel,ng_vel,ng_vel);
            icomp  = 0;
            ncomp  = 1;
        }
        else if (var_idx == Vars::zvel)
        {
            bccomp = BCVars::zvel_bc;
            mapper = &face_linear_interp;
            ngvect = IntVect(ng_vel,ng_vel,0);
            icomp  = 0;
            ncomp  = 1;
        }

        if (lev == 0)
        {
            mf.FillBoundary(icomp,ncomp,ngvect,geom[lev].periodicity());
        }
        else
        {
            Vector<MultiFab*> fmf = {&mf};
            Vector<MultiFab*> cmf = {&vars_old[lev-1][var_idx], &vars_new[lev-1][var_idx]};
            Vector<Real> ctime    = {t_old[lev-1], t_new[lev-1]};

            amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, {time},
                                      icomp, icomp, ncomp, geom[lev-1], geom[lev],
                                      null_bc, 0, null_bc, 0, refRatio(lev-1),
                                      mapper, domain_bcs_type, bccomp);
        }
    }

    // ***************************************************************************
    // Physical bc's at domain boundary
    // ***************************************************************************
    IntVect ngvect_cons = IntVect(ng_cons,ng_cons,ng_cons);
    IntVect ngvect_vels = IntVect(ng_vel ,ng_vel ,ng_vel);

#ifdef ERF_USE_NETCDF
    // We call this here because it is an ERF routine
    if (init_type == "real") fill_from_wrfbdy(mfs,time);
    if (init_type == "metgrid") fill_from_metgrid(mfs,time);
#endif

    if (m_r2d) fill_from_bndryregs(mfs,time);

    // We call this even if init_type == real because this routine will fill the vertical bcs
    (*physbcs[lev])(mfs,icomp_cons,ncomp_cons,ngvect_cons,ngvect_vels,init_type,cons_only,BCVars::cons_bc);
    // ***************************************************************************

    //
    // It is important that we apply the MOST bcs after we have imposed all the others
    //    so that we have enough information in the ghost cells to calculate the viscosity
    //
    if (!(cons_only && ncomp_cons == 1) && m_most && allow_most_bcs)
        m_most->impose_most_bcs(lev,mfs,eddyDiffs);
}

//

/*
 * Fill valid and ghost data.
 * This version sill an entire MultiFab by interpolating from the coarser level -- this is used
 * only when a new level of refinement is being created during a run (i.e not at initialization)
 * This will never be used with static refinement.
 *
 * @param[in]  lev            level of refinement at which to fill the data
 * @param[in]  time           time at which the data should be filled
 * @param[out] mfs            Vector of MultiFabs to be filled containing, in order: cons, xvel, yvel, and zvel data
 */

void
ERF::FillCoarsePatch (int lev, Real time, const Vector<MultiFab*>& mfs)
{
    BL_PROFILE_VAR("FillCoarsePatch()",FillCoarsePatch);
    AMREX_ASSERT(lev > 0);

    int icomp = 0;
    int bccomp;
    amrex::Interpolater* mapper;

    for (int var_idx = 0; var_idx < mfs.size(); ++var_idx) {

        if (var_idx == Vars::cons)
        {
            bccomp = 0;
            bccomp = 0;
            mapper = &cell_cons_interp;
        }
        else if (var_idx == Vars::xvel || var_idx == Vars::xmom)
        {
            bccomp = NVAR;
            mapper = &face_linear_interp;
        }
        else if (var_idx == Vars::yvel || var_idx == Vars::ymom)
        {
            bccomp = NVAR+1;
            mapper = &face_linear_interp;
        }
        else if (var_idx == Vars::zvel || var_idx == Vars::zmom)
        {
            bccomp = NVAR+2;
            mapper = &face_linear_interp;
        }

        amrex::InterpFromCoarseLevel(*mfs[var_idx], time, vars_new[lev-1][var_idx],
                                     icomp, icomp, mfs[var_idx]->nComp(),
                                     geom[lev-1], geom[lev],
                                     null_bc, 0, null_bc, 0, refRatio(lev-1),
                                     mapper, domain_bcs_type, bccomp);
    }

    // ***************************************************************************
    // Physical bc's at domain boundary
    // ***************************************************************************
    IntVect ngvect_cons = mfs[Vars::cons]->nGrowVect();
    IntVect ngvect_vels = mfs[Vars::xvel]->nGrowVect();
    bool cons_only = false;

    (*physbcs[lev])(mfs,0,mfs[Vars::cons]->nComp(),ngvect_cons,ngvect_vels,init_type,cons_only,BCVars::cons_bc);

    // ***************************************************************************
    // Since lev > 0 here we don't worry about m_r2d or wrfbdy data
    // ***************************************************************************
}
