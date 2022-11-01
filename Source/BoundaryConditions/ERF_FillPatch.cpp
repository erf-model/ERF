#include <ERF.H>
#include <ERF_PhysBCFunct.H>
#include <IndexDefines.H>
#include <TimeInterpolatedData.H>

using namespace amrex;

PhysBCFunctNoOp null_bc;

//
// Fill valid and ghost data in the MultiFab "mf"
// This version fills the MultiFab mf in valid regions with the "state data" at the given time;
// values in mf when it is passed in are *not* used.
//
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
            bccomp = 0;
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

    (*physbcs[lev])(mfs,icomp_cons,ncomp_cons,ngvect_cons,ngvect_vels,time,cons_only);

    if (m_r2d) fill_from_bndryregs(mfs,time);
#ifdef ERF_USE_NETCDF
    if (init_type == "real") fill_from_wrfbdy(mfs,time);
#endif
}

//
// Fill valid and ghost data in the MultiFabs in "mfs"
// mfs is a Vector<std::reference_wrapper<MultiFab> > containing, in order: cons, xvel, yvel, and zvel data
// This version fills the MultiFabs mfs in valid regions with the values in "mfs" when it is passed in;
// it is used only to compute ghost values for intermediate stages of a time integrator.
//
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

    (*physbcs[lev])(mfs,icomp_cons,ncomp_cons,ngvect_cons,ngvect_vels,time,cons_only);
    // ***************************************************************************

    if (m_r2d) fill_from_bndryregs(mfs,time);
#ifdef ERF_USE_NETCDF
    if (init_type == "real") fill_from_wrfbdy(mfs,time);
#endif

    //
    // It is important that we apply the MOST bcs after we have imposed all the others
    //    so that we have enough information in the ghost cells to calculate the viscosity
    //
    if (!(cons_only && ncomp_cons == 1) && m_most && allow_most_bcs)
    {
        /*
        const int icomp = 0;
        for (MFIter mfi(*mfs[0]); mfi.isValid(); ++mfi)
        {
            const auto cons_arr = (*mfs[Vars::cons])[mfi].array();
            const auto velx_arr = (*mfs[Vars::xvel])[mfi].array();
            const auto vely_arr = (*mfs[Vars::yvel])[mfi].array();
            const auto  eta_arr = (*eddyDiffs)[mfi].array();

            for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx)
            {
                const Box& bx = (*mfs[var_idx])[mfi].box();
                auto dest_arr = (*mfs[var_idx])[mfi].array();
                int zlo = 0;
                m_most->impose_most_bcs(lev,bx,dest_arr,cons_arr,velx_arr,vely_arr,eta_arr,var_idx,icomp,zlo);
            } // var_idx
        } // mf
        */

        m_most->impose_most_bcs(lev,mfs,eddyDiffs);
    } // most
}

// Fill an entire multifab by interpolating from the coarser level -- this is used
//     only when a new level of refinement is being created during a run (i.e not at initialization)
//     This will never be used with static refinement.
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

    (*physbcs[lev])(mfs,0,mfs[Vars::cons]->nComp(),ngvect_cons,ngvect_vels,time,cons_only);

    // ***************************************************************************
    // Since lev > 0 here we don't worry about m_r2d or wrfbdy data
    // ***************************************************************************
}
