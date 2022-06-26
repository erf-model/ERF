#include <ERF.H>
#include <ERF_PhysBCFunct.H>
#include <IndexDefines.H>
#include <TimeInterpolatedData.H>

using namespace amrex;

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
            MultiFab* mf_temp = new MultiFab(vars_new[lev][i].boxArray(),
                                             vars_new[lev][i].DistributionMap(),
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

    // We need to make sure to fill these before we compute the viscosity
    for (int i = 0; i < Vars::NumTypes; ++i) {
        data.get_var(Vars::xvel).FillBoundary(geom[lev].periodicity());
        data.get_var(Vars::yvel).FillBoundary(geom[lev].periodicity());
        data.get_var(Vars::zvel).FillBoundary(geom[lev].periodicity());
        data.get_var(Vars::cons).FillBoundary(geom[lev].periodicity());
    }

    return data;
}
//
// Fill valid and ghost data in the MultiFab "mf"
// This version fills the MultiFab mf in valid regions with the "state data" at the given time;
// values in mf when it is passed in are *not* used.
//
void
ERF::FillPatch (int lev, Real time, Vector<MultiFab>& mfs)
{
    int bccomp;
    amrex::Interpolater* mapper = nullptr;

    TimeInterpolatedData fdata = GetDataAtTime(lev, time);
    Vector<Real> ftime         = {fdata.get_time()};

    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx) {
        MultiFab& mf = mfs[var_idx];
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
            Vector<MultiFab*> smf = {&fdata.get_var(var_idx)};
            ERFPhysBCFunct physbc(lev,geom[lev],domain_bcs_type,domain_bcs_type_d,var_idx,fdata,
                                  m_bc_extdir_vals,
#ifdef ERF_USE_TERRAIN
                                  z_phys_nd[lev], detJ_cc[lev],
#endif
#ifdef ERF_USE_NETCDF
                                  bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi, bdy_time_interval,
#endif
                                  m_r2d);
            amrex::FillPatchSingleLevel(mf, time, smf, ftime, 0, icomp, ncomp,
                                        geom[lev], physbc, bccomp);
        }
        else
        {
            TimeInterpolatedData cdata = GetDataAtTime(lev-1, time);
            Vector<Real> ctime = {cdata.get_time()};
            Vector<MultiFab*> cmf = {&cdata.get_var(var_idx)};
            Vector<MultiFab*> fmf = {&fdata.get_var(var_idx)};

            ERFPhysBCFunct cphysbc(lev-1,geom[lev-1],domain_bcs_type,domain_bcs_type_d,var_idx,cdata,
                                   m_bc_extdir_vals,
#ifdef ERF_USE_TERRAIN
                                   z_phys_nd[lev-1],detJ_cc[lev-1],
#endif
#ifdef ERF_USE_NETCDF
                                   bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi, bdy_time_interval,
#endif
                                   m_r2d);
            ERFPhysBCFunct fphysbc(lev,geom[lev],domain_bcs_type,domain_bcs_type_d,var_idx,fdata,
                                   m_bc_extdir_vals,
#ifdef ERF_USE_TERRAIN
                                   z_phys_nd[lev],detJ_cc[lev],
#endif
#ifdef ERF_USE_NETCDF
                                   bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi, bdy_time_interval,
#endif
                                   m_r2d);

            amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                      0, icomp, ncomp, geom[lev-1], geom[lev],
                                      cphysbc, bccomp, fphysbc, bccomp, refRatio(lev-1),
                                      mapper, domain_bcs_type, bccomp);
        } // lev > 0
    } // var_idx

    //
    // It is important that we apply the MOST bcs after we have imposed all the others
    //    so that we have enough information in the ghost cells to calculate the viscosity
    // Note that we don't test on the BCRec's for these variables because those have been set
    //      to match those of no-slip-wall in order to fill values before computing viscosity
    //
    if (m_most)
    {
        MultiFab eddyDiffs(fdata.get_var(Vars::cons).boxArray(),
                           fdata.get_var(Vars::cons).DistributionMap(),
                           EddyDiff::NumDiffs,3);
        bool vert_only = true;
        ComputeTurbulentViscosity(fdata.get_var(Vars::xvel), fdata.get_var(Vars::yvel),
                                  fdata.get_var(Vars::zvel), fdata.get_var(Vars::cons),
                                  eddyDiffs, geom[lev], solverChoice, m_most, domain_bcs_type_d, vert_only);
        eddyDiffs.FillBoundary(geom[lev].periodicity());

        for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx)
        {
            MultiFab& mf = mfs[var_idx];
            const int icomp = 0;

            for (MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                const Box& bx       = mfs[var_idx][mfi].box();
                      auto dest_arr = mfs[var_idx][mfi].array();

                const auto velx_arr = fdata.get_var(Vars::xvel)[mfi].array();
                const auto vely_arr = fdata.get_var(Vars::yvel)[mfi].array();
                const auto cons_arr = fdata.get_var(Vars::cons)[mfi].array();
                const auto  eta_arr = eddyDiffs[mfi].array();

                int zlo = 0;
                m_most->impose_most_bcs(lev,bx,dest_arr,cons_arr,velx_arr,vely_arr,eta_arr,var_idx,icomp,zlo);
            } // mf
        } // var_idx
    } // most
}

//
// Fill valid and ghost data in the MultiFabs in "mfs"
// mfs is a Vector<std::reference_wrapper<MultiFab> > containing, in order: cons, xvel, yvel, and zvel data
// This version fills the MultiFabs mfs in valid regions with the values in "mfs" when it is passed in;
// it is used only to compute ghost values for intermediate stages of a time integrator.
//
void
ERF::FillIntermediatePatch (int lev, Real time, Vector<std::reference_wrapper<MultiFab> > mfs,
                            bool rho_only)
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

    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx)
    {
        if (rho_only && var_idx != Vars::cons) continue;

        MultiFab& mf = mfs[var_idx].get();
        const int icomp = 0;
        int ncomp = mf.nComp();

        if (rho_only && var_idx == Vars::cons) ncomp = 1;

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

        if (lev == 0)
        {
            // on lev, use the mf data and time passed to FillIntermediatePatch().
            Vector<MultiFab*> smf { &mf };
            Vector<Real> stime { time };

            ERFPhysBCFunct physbc(lev,geom[lev],domain_bcs_type,domain_bcs_type_d,var_idx,level_data,
                                  m_bc_extdir_vals,
#ifdef ERF_USE_TERRAIN
                                  z_phys_nd[lev],detJ_cc[lev],
#endif
#ifdef ERF_USE_NETCDF
                                   bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi, bdy_time_interval,
#endif
                                   m_r2d);

            amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                        geom[lev], physbc, bccomp);
        }
        else
        {
            MultiFab mf_temp(mf.boxArray(), mf.DistributionMap(), ncomp, mf.nGrowVect());

            TimeInterpolatedData cdata = GetDataAtTime(lev-1, time);
            Vector<MultiFab*> cmf = {&cdata.get_var(var_idx)};
            Vector<MultiFab*> fmf = {&mf};
            Vector<Real> ctime = {cdata.get_time()};
            Vector<Real> ftime = {level_data.get_time()};

            ERFPhysBCFunct cphysbc(lev-1,geom[lev-1],domain_bcs_type,domain_bcs_type_d,var_idx,cdata,
                                  m_bc_extdir_vals,
#ifdef ERF_USE_TERRAIN
                                   z_phys_nd[lev-1],detJ_cc[lev-1],
#endif
#ifdef ERF_USE_NETCDF
                                   bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi, bdy_time_interval,
#endif
                                   m_r2d);
            ERFPhysBCFunct fphysbc(lev,geom[lev],domain_bcs_type,domain_bcs_type_d,var_idx,level_data,
                                  m_bc_extdir_vals,
#ifdef ERF_USE_TERRAIN
                                   z_phys_nd[lev],detJ_cc[lev],
#endif
#ifdef ERF_USE_NETCDF
                                   bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi, bdy_time_interval,
#endif
                                   m_r2d);

            amrex::FillPatchTwoLevels(mf_temp, time, cmf, ctime, fmf, ftime,
                                    0, icomp, ncomp, geom[lev-1], geom[lev],
                                    cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                    mapper, domain_bcs_type, bccomp);

            // Replace mf with mf_temp
            if (ncomp == mf.nComp())
                std::swap(mf_temp, mf);
            else
                MultiFab::Copy(mf,mf_temp,0,0,1,mf.nGrowVect());
        }
    }

    //
    // It is important that we apply the MOST bcs after we have imposed all the others
    //    so that we have enough information in the ghost cells to calculate the viscosity
    //
    if (!rho_only && m_most)
    {
        MultiFab eddyDiffs(mfs[Vars::cons].get().boxArray(),
                           mfs[Vars::cons].get().DistributionMap(),
                           EddyDiff::NumDiffs,3);
        bool vert_only = true;
        ComputeTurbulentViscosity(mfs[Vars::xvel].get(), mfs[Vars::yvel].get(),
                                  mfs[Vars::zvel].get(), mfs[Vars::cons].get(),
                                  eddyDiffs, geom[lev], solverChoice, m_most, domain_bcs_type_d, vert_only);
        eddyDiffs.FillBoundary(geom[lev].periodicity());

        for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx)
        {
            MultiFab& mf = mfs[var_idx].get();
            const int icomp = 0;

            for (MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                const Box& bx       = mf[mfi].box();
                      auto dest_arr = mf[mfi].array();

                const auto velx_arr = level_data.get_var(Vars::xvel)[mfi].array();
                const auto vely_arr = level_data.get_var(Vars::yvel)[mfi].array();
                const auto cons_arr = level_data.get_var(Vars::cons)[mfi].array();
                const auto  eta_arr = eddyDiffs[mfi].array();

                int zlo = 0;
                m_most->impose_most_bcs(lev,bx,dest_arr,cons_arr,velx_arr,vely_arr,eta_arr,var_idx,icomp,zlo);
            } // mf
        } // var_idx
    } // most
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

    ERFPhysBCFunct cphysbc(lev-1,geom[lev-1],domain_bcs_type,domain_bcs_type_d,var_idx,cdata,
                           m_bc_extdir_vals,
#ifdef ERF_USE_TERRAIN
                           z_phys_nd[lev-1],detJ_cc[lev-1],
#endif
#ifdef ERF_USE_NETCDF
                           bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi, bdy_time_interval,
#endif
                           m_r2d);
    ERFPhysBCFunct fphysbc(lev,geom[lev],domain_bcs_type,domain_bcs_type_d,var_idx,fdata,
                           m_bc_extdir_vals,
#ifdef ERF_USE_TERRAIN
                           z_phys_nd[lev],detJ_cc[lev],
#endif
#ifdef ERF_USE_NETCDF
                           bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi, bdy_time_interval,
#endif
                           m_r2d);

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
