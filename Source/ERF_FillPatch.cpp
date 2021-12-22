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
            MultiFab* mf_temp = new MultiFab(grids[lev], dmap[lev], vars_new[lev][i].nComp(), vars_new[lev][i].nGrowVect());
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

Vector<BCRec>
ERF::GetBCVarVector (int var_idx)
{
    Vector<BCRec> bcs_var;
    if (var_idx == Vars::cons) {
        bcs_var = bcs[BCVars::cons];
    } else if (var_idx == Vars::xvel) {
        bcs_var = {bcs[BCVars::vels][0]};
    } else if (var_idx == Vars::yvel) {
        bcs_var = {bcs[BCVars::vels][1]};
    } else if (var_idx == Vars::zvel) {
        bcs_var = {bcs[BCVars::vels][2]};
    } else {
        Error("Invalid var_idx passed to GetBCVarVector. It should be in the Vars enum.");
    }
    return bcs_var;
}

// compute a new multifab by copying in data from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
ERF::FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp, int var_idx)
{
    if (lev == 0)
    {
        TimeInterpolatedData sdata = GetDataAtTime(lev, time);
        Vector<MultiFab*> smf = {&sdata.get_var(var_idx)};
        Vector<Real> stime = {sdata.get_time()};

        Vector<BCRec> bcs_var = GetBCVarVector(var_idx);
        ERFPhysBCFunct physbc(geom[lev],bcs_var,var_idx,sdata,m_bc_cons,m_bc_vels);
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                    geom[lev], physbc, 0);

    }
    else
    {
        TimeInterpolatedData cdata = GetDataAtTime(lev-1, time);
        TimeInterpolatedData fdata = GetDataAtTime(lev  , time);
        Vector<MultiFab*> cmf = {&cdata.get_var(var_idx)};
        Vector<MultiFab*> fmf = {&fdata.get_var(var_idx)};
        Vector<Real> ctime = {cdata.get_time()};
        Vector<Real> ftime = {fdata.get_time()};

        Vector<BCRec> bcs_var = GetBCVarVector(var_idx);
        ERFPhysBCFunct cphysbc(geom[lev-1],bcs_var,var_idx,cdata,m_bc_cons,m_bc_vels);
        ERFPhysBCFunct fphysbc(geom[lev  ],bcs_var,var_idx,fdata,m_bc_cons,m_bc_vels);

        amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                  0, icomp, ncomp, geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                  mapper, bcs_var, 0);
    }
}

// Compute a new multifab by copying in data from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
// unlike FillPatch, FillIntermediatePatch will use the supplied multifab instead of fine level data.
// This is to support filling boundary cells at an intermediate time between old/new times
// on the fine level when valid data at a specific time is already available (such as
// at each RK stage when integrating between initial and final times at a given level).
void
ERF::FillIntermediatePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp, int var_idx)
{
    if (lev == 0)
    {
        // on lev, use the mf data and time passed to FillIntermediatePatch().
        TimeInterpolatedData sdata = GetDataAtTime(lev, time);
        Vector<MultiFab*> smf { &mf };
        Vector<Real> stime { time };

        Vector<BCRec> bcs_var = GetBCVarVector(var_idx);
        ERFPhysBCFunct physbc(geom[lev],bcs_var,var_idx,sdata,m_bc_cons,m_bc_vels);
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                    geom[lev], physbc, 0);
    }
    else
    {
        MultiFab mf_temp(mf.boxArray(), mf.DistributionMap(), mf.nComp(), mf.nGrowVect());

        TimeInterpolatedData cdata = GetDataAtTime(lev-1, time);
        TimeInterpolatedData fdata = GetDataAtTime(lev  , time);
        Vector<MultiFab*> cmf = {&cdata.get_var(var_idx)};
        Vector<MultiFab*> fmf = {&mf};
        Vector<Real> ctime = {cdata.get_time()};
        Vector<Real> ftime = {fdata.get_time()};

        Vector<BCRec> bcs_var = GetBCVarVector(var_idx);
        ERFPhysBCFunct cphysbc(geom[lev-1],bcs_var,var_idx,cdata,m_bc_cons,m_bc_vels);
        ERFPhysBCFunct fphysbc(geom[lev  ],bcs_var,var_idx,fdata,m_bc_cons,m_bc_vels);

        amrex::FillPatchTwoLevels(mf_temp, time, cmf, ctime, fmf, ftime,
                                  0, icomp, ncomp, geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                  mapper, bcs_var, 0);

        // Replace mf with mf_temp
        std::swap(mf_temp, mf);
    }
}

// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
void
ERF::FillCoarsePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp, int var_idx)
{
    AMREX_ASSERT(lev > 0);

    TimeInterpolatedData cdata = GetDataAtTime(lev-1, time);
    TimeInterpolatedData fdata = GetDataAtTime(lev  , time);
    Vector<MultiFab*> cmf = {&cdata.get_var(var_idx)};
    Vector<MultiFab*> fmf = {&fdata.get_var(var_idx)};
    Vector<Real> ctime = {cdata.get_time()};
    Vector<Real> ftime = {fdata.get_time()};

    Vector<BCRec> bcs_var = GetBCVarVector(var_idx);
    ERFPhysBCFunct cphysbc(geom[lev-1],bcs_var,var_idx,cdata,m_bc_cons,m_bc_vels);
    ERFPhysBCFunct fphysbc(geom[lev  ],bcs_var,var_idx,fdata,m_bc_cons,m_bc_vels);

    amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
                                    cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                    mapper, bcs_var, 0);
}

void
ERF::FillCoarsePatchAllVars (int lev, Real time, Vector<MultiFab>& vmf)
{
    for (int var_idx = 0; var_idx < vmf.size(); ++var_idx) {
        FillCoarsePatch(lev, time, vmf[var_idx], 0, vmf[var_idx].nComp(), var_idx);
    }
}
