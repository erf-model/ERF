#include <ERF.H>

// utility to copy in data from old/new data into a struct that holds data for FillPatching
TimeInterpolatedData
ERF::GetDataAtTime (int lev, Real time)
{
    TimeInterpolatedData data;

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
        for (int i = 0; i < Vars::NumTypes; ++i) {
            data.add_var(&vars_new[lev][i], data.non_owning);
        }
        data.set_time(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
        for (int i = 0; i < Vars::NumTypes; ++i) {
            data.add_var(&vars_old[lev][i], data.non_owning);
        }
        data.set_time(t_old[lev]);
    }
    else if (time > t_old[lev] && time < t_new[lev])
    {
        // do first order interpolation in time between [t_old[lev], t_new[lev]]
        // time interpolation includes the ghost cells
        for (int i = 0; i < Vars::NumTypes; ++i) {
            MultiFab* mf_temp = new MultiFab(vars_new[lev][i].boxArray(), dmap[lev],
                                             vars_new[lev][i].nComp(), vars_new[lev][i].nGrowVect());
            mf_temp->setVal(0.0_rt);

            const Real dt_fraction = (time - t_old[lev]) / (t_new[lev] - t_old[lev]);
            MultiFab::Saxpy(*mf_temp, 1.0_rt - dt_fraction, vars_old[lev][i], 0, 0, mf_temp->nComp(), mf_temp->nGrowVect());
            MultiFab::Saxpy(*mf_temp,          dt_fraction, vars_new[lev][i], 0, 0, mf_temp->nComp(), mf_temp->nGrowVect());

            data.add_var(mf_temp, data.owning);
        }
        data.set_time(time);
    }
    else
    {
        amrex::Error("Requested data at a time outside the interval [t_old, t_new]");
    }

    return data;
}
//
// Fill valid and ghost data in the MultiFab "mf"
// This version fills the MultiFab mf in valid regions with the "state data" at the given time;
// values in mf when it is passed in are *not* used.
//
void
ERF::FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp, int var_idx)
{
    int bccomp;
    amrex::Interpolater* mapper = nullptr;

    if (var_idx == Vars::cons)
    {
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

    if (lev == 0)
    {
        TimeInterpolatedData sdata = GetDataAtTime(lev, time);
        Vector<MultiFab*> smf = {&sdata.get_var(var_idx)}; // todo: if var_idx is momentum, then all instances of get_var(var_idx) will have to get velocity & cons and then convert to momentum
        Vector<Real> stime = {sdata.get_time()};

        ERFPhysBCFunct physbc(geom[lev],domain_bcs_type,domain_bcs_type_d,var_idx,sdata,m_bc_extdir_vals,solverChoice,
#ifdef ERF_USE_TERRAIN
                              z_phys_nd[lev],detJ_cc[lev],
#endif
                              get_most(),m_r2d);
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                    geom[lev], physbc, bccomp);

    }
    else
    {
        TimeInterpolatedData cdata = GetDataAtTime(lev-1, time);
        TimeInterpolatedData fdata = GetDataAtTime(lev  , time);
        Vector<MultiFab*> cmf = {&cdata.get_var(var_idx)};
        Vector<MultiFab*> fmf = {&fdata.get_var(var_idx)};
        Vector<Real> ctime = {cdata.get_time()};
        Vector<Real> ftime = {fdata.get_time()};

        ERFPhysBCFunct cphysbc(geom[lev-1],domain_bcs_type,domain_bcs_type_d,var_idx,cdata,m_bc_extdir_vals,solverChoice,
#ifdef ERF_USE_TERRAIN
                               z_phys_nd[lev-1],detJ_cc[lev-1],
#endif
                               get_most(),m_r2d);
        ERFPhysBCFunct fphysbc(geom[lev],domain_bcs_type,domain_bcs_type_d,var_idx,fdata,m_bc_extdir_vals,solverChoice,
#ifdef ERF_USE_TERRAIN
                               z_phys_nd[lev],detJ_cc[lev],
#endif
                               get_most(),m_r2d);

        amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                  0, icomp, ncomp, geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                  mapper, domain_bcs_type, bccomp);
    }
}

//
// Fill valid and ghost data in the MultiFabs in "mfs"
// mfs is a Vector<std::reference_wrapper<MultiFab> > containing, in order: cons, xvel, yvel, and zvel data
// This version fills the MultiFabs mfs in valid regions with the values in "mfs" when it is passed in;
// it is used only to compute ghost values for interemediate stages of a time integrator.
//
void
ERF::FillIntermediatePatch (int lev, Real time, Vector<std::reference_wrapper<MultiFab> > mfs, int which)
{
    int bccomp;
    amrex::Interpolater* mapper;
    TimeInterpolatedData level_data;

    // We should always pass cons, xvel, yvel, and zvel (in that order) in the mfs vector
    AMREX_ALWAYS_ASSERT(mfs.size() == Vars::NumTypes);

    for (int imf = 0; imf < Vars::NumTypes; ++imf) {
        level_data.add_var(&mfs[imf].get(), level_data.non_owning);
    }
    level_data.set_time(time);

    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx) {
        if (which >= 0 && which != var_idx) continue;
        MultiFab& mf = mfs[var_idx].get();
        const int icomp = 0;
        const int ncomp = mf.nComp();

        if (var_idx == Vars::cons)
        {
            bccomp = 0;
            mapper = &cell_cons_interp;
        }
        else if (var_idx == Vars::xvel)
        {
            bccomp = NVAR;
            mapper = &face_linear_interp;
        }
        else if (var_idx == Vars::yvel)
        {
            bccomp = NVAR+1;
            mapper = &face_linear_interp;
        }
        else if (var_idx == Vars::zvel)
        {
            bccomp = NVAR+2;
            mapper = &face_linear_interp;
        }
        else
        {

        }

        if (lev == 0)
        {
            // on lev, use the mf data and time passed to FillIntermediatePatch().
            Vector<MultiFab*> smf { &mf };
            Vector<Real> stime { time };

            ERFPhysBCFunct physbc(geom[lev],domain_bcs_type,domain_bcs_type_d,var_idx,level_data,m_bc_extdir_vals,solverChoice,
#ifdef ERF_USE_TERRAIN
                                  z_phys_nd[lev],detJ_cc[lev],
#endif
                                  get_most(),m_r2d);
            amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                        geom[lev], physbc, bccomp);
        }
        else
        {
            MultiFab mf_temp(mf.boxArray(), mf.DistributionMap(), mf.nComp(), mf.nGrowVect());

            TimeInterpolatedData cdata = GetDataAtTime(lev-1, time);
            Vector<MultiFab*> cmf = {&cdata.get_var(var_idx)};
            Vector<MultiFab*> fmf = {&mf};
            Vector<Real> ctime = {cdata.get_time()};
            Vector<Real> ftime = {level_data.get_time()};

            ERFPhysBCFunct cphysbc(geom[lev-1],domain_bcs_type,domain_bcs_type_d,var_idx,cdata,m_bc_extdir_vals,solverChoice,
#ifdef ERF_USE_TERRAIN
                                   z_phys_nd[lev-1],detJ_cc[lev-1],
#endif
                                   get_most(),m_r2d);
            ERFPhysBCFunct fphysbc(geom[lev],domain_bcs_type,domain_bcs_type_d,var_idx,level_data,m_bc_extdir_vals,solverChoice,
#ifdef ERF_USE_TERRAIN
                                   z_phys_nd[lev],detJ_cc[lev],
#endif
                                   get_most(),m_r2d);

            amrex::FillPatchTwoLevels(mf_temp, time, cmf, ctime, fmf, ftime,
                                    0, icomp, ncomp, geom[lev-1], geom[lev],
                                    cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                    mapper, domain_bcs_type, bccomp);

            // Replace mf with mf_temp
            std::swap(mf_temp, mf);
        }
    }
}

// Fill an entire multifab by interpolating from the coarser level -- this is used
//     only when a new level of refinement is being created during a run (i.e not at initialization)
//     This will never be used with static refinement.
void
ERF::FillCoarsePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp, int var_idx)
{
    AMREX_ASSERT(lev > 0);

    int bccomp;
    amrex::Interpolater* mapper;

    if (var_idx == Vars::cons)
    {
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

    TimeInterpolatedData cdata = GetDataAtTime(lev-1, time);
    TimeInterpolatedData fdata = GetDataAtTime(lev  , time);
    Vector<MultiFab*> cmf = {&cdata.get_var(var_idx)};
    Vector<MultiFab*> fmf = {&fdata.get_var(var_idx)};
    Vector<Real> ctime = {cdata.get_time()};
    Vector<Real> ftime = {fdata.get_time()};

    ERFPhysBCFunct cphysbc(geom[lev-1],domain_bcs_type,domain_bcs_type_d,var_idx,cdata,m_bc_extdir_vals,solverChoice,
#ifdef ERF_USE_TERRAIN
                           z_phys_nd[lev-1],detJ_cc[lev-1],
#endif
                           get_most(),m_r2d);
    ERFPhysBCFunct fphysbc(geom[lev],domain_bcs_type,domain_bcs_type_d,var_idx,fdata,m_bc_extdir_vals,solverChoice,
#ifdef ERF_USE_TERRAIN
                           z_phys_nd[lev],detJ_cc[lev],
#endif
                           get_most(),m_r2d);

    amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
                                    cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                    mapper, domain_bcs_type, bccomp);
}

void
ERF::FillCoarsePatchAllVars (int lev, Real time, Vector<MultiFab>& vmf)
{
    for (int var_idx = 0; var_idx < vmf.size(); ++var_idx) {
        FillCoarsePatch(lev, time, vmf[var_idx], 0, vmf[var_idx].nComp(), var_idx);
    }
}
