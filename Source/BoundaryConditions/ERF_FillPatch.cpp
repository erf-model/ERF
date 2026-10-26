#include <ERF.H>
#include "ERF_Constants.H"
#include <ERF_PhysBCFunct.H>
#include <ERF_IndexDefines.H>
#include <ERF_TimeInterpolatedData.H>
#include <ERF_FillPatcher.H>
#include <ERF_Utils.H>

using namespace amrex;

/**
 * Fill valid and ghost data with the "state data" at the given time
 * NOTE: THIS OPERATES ON VELOCITY (MOMENTA ARE JUST TEMPORARIES)
 *
 * @param[in]     lev            level of refinement at which to fill the data
 * @param[in]     time_d         time at which the data should be filled
 * @param[in,out] mfs_vel        Vector of MultiFabs to be filled containing, in order: cons, xvel, yvel, and zvel
 * @param[in,out] mfs_mom        Vector of MultiFabs to be filled containing, in order: cons, xmom, ymom, and zmom
 * @param[in]     old_base_state base-state data at the old time
 * @param[in]     new_base_state base-state data at the new time
 * @param[in]     fillset        whether to fill the coarse-fine set region before standard fillpatch
 * @param[in]     cons_only      whether to fill only conserved variables
 */
void
ERF::FillPatchFineLevel (int lev, double time_d,
                         const Vector<MultiFab*>& mfs_vel,     // This includes cc quantities and VELOCITIES
                         const Vector<MultiFab*>& mfs_mom,     // This includes cc quantities and MOMENTA
                         const MultiFab& old_base_state,
                         const MultiFab& new_base_state,
                         bool fillset, bool cons_only)
{
    BL_PROFILE_VAR("ERF::FillPatchFineLevel()",ERF_FillPatchFineLevel);

    AMREX_ALWAYS_ASSERT(lev > 0);

    Real time = static_cast<Real>(time_d);

    Interpolater* mapper = nullptr;

    PhysBCFunctNoOp null_bc;

    IntVect ngvect_cons = mfs_vel[Vars::cons]->nGrowVect();
    IntVect ngvect_vels = mfs_vel[Vars::xvel]->nGrowVect();

    // Which time level(s) to take on each level is decided in double precision
    const TimeLevelSel fsel = select_time_level(time_d, t_old[lev  ], t_new[lev  ]);
    const TimeLevelSel csel = select_time_level(time_d, t_old[lev-1], t_new[lev-1]);

    // AMReX computes its time-interpolation weights in Real from the times it is given,
    // so give it times measured from the coarse new time (formed in double) rather than
    // the elapsed times themselves, which can be large compared to dt.  The physical
    // bc's below still see the elapsed time.
    const double tref = t_new[lev-1];
    const Real time_rel = static_cast<Real>(time_d - tref);
    Vector<Real> ftime  = {static_cast<Real>(t_old[lev  ] - tref), static_cast<Real>(t_new[lev  ] - tref)};
    Vector<Real> ctime  = {static_cast<Real>(t_old[lev-1] - tref), static_cast<Real>(t_new[lev-1] - tref)};

    Vector<MultiFab*> fmf = select_time_data(fsel, vars_old[lev  ][Vars::cons], vars_new[lev  ][Vars::cons]);
    Vector<MultiFab*> cmf = select_time_data(csel, vars_old[lev-1][Vars::cons], vars_new[lev-1][Vars::cons]);

    //
    // ***************************************************************************
    // The momentum conversion below reads the density on both sides of every face
    // of the valid box -- one ghost cell in each face-normal direction -- and at
    // this point in the step nothing has filled those ghost cells yet: they still
    // hold whatever last wrote them.  That made the interface momenta depend on
    // the run's output history, because writing a 3D plotfile fills them (it calls
    // this routine with fillset false), so a run that wrote one got different
    // momenta on the step that followed than a run that did not, and no two-level
    // restart could be bit-exact (issue 4224).
    //
    // So fill exactly what the conversion reads first, and no more: Rho_comp, one
    // ghost cell.  The whole conserved state is filled below anyway, so anything
    // wider would be thrown away.  Perturbational interpolation is not the one
    // used here, and the density the conversion sees has to be the density the
    // state ends up with, so that mode fills the conserved state the way it will
    // be filled below -- at the cost of doing it twice -- rather than taking the
    // cheap path and being inconsistent with itself.
    // ***************************************************************************
    if (fillset && !cons_only) {
        //
        // Perturbational interpolation cannot reach this: it is only ever used from
        // here, and ERF::ReadParameters refuses it together with amr.max_level > 0,
        // which this routine requires (issue 4257).  Were it allowed, the density the
        // conversion needs could only be had by running the transform again, which
        // would perturb the live state a second time per fill -- so the mode has to
        // earn a deck and a test before it comes back, not be given a cheap path here.
        //
        AMREX_ALWAYS_ASSERT(interpolation_type != StateInterpType::Perturbational);

        {
            const IntVect ng_rho(1,1,1);
            MultiFab rho_mf(mfs_vel[Vars::cons]->boxArray(),
                            mfs_vel[Vars::cons]->DistributionMap(), 1, ng_rho);
            FillPatchTwoLevels(rho_mf, ng_rho, IntVect(0,0,0),
                               time_rel, cmf, ctime, fmf, ftime,
                               Rho_comp, 0, 1, geom[lev-1], geom[lev],
                               refRatio(lev-1), &cell_cons_interp, domain_bcs_type,
                               BCVars::cons_bc);
            MultiFab::Copy(*mfs_vel[Vars::cons], rho_mf, 0, Rho_comp, 1, ng_rho);
            //
            // A density ingested from boundary planes is not imposed by physbcs_cons, so
            // apply it here, in the same order as below, or the conversion would read a
            // stale ghost density at the ingested faces (issue 4265).  This writes the
            // ingested ghost cells of every conserved component; all of them are
            // overwritten when the conserved state is filled below.
            //
            if (m_r2d && !solverChoice.use_real_bcs) {
                fill_from_bndryregs(lev, mfs_vel, time_d, true);
            }
            //
            // FillPatchTwoLevels is told to fill nothing outside the domain, and
            // MomentumToVelocity does not clamp its z faces to the domain the way
            // VelocityToMomentum does, so the cells under and over the domain have
            // to come from the bcs.  The velocities passed here are only read at
            // the lateral domain faces for ext_dir_upwind, and those are valid
            // faces whose data is already correct at this point.
            //
            (*physbcs_cons[lev])(*mfs_vel[Vars::cons],
                                 *mfs_vel[Vars::xvel], *mfs_vel[Vars::yvel],
                                 Rho_comp, 1, ng_rho, time, BCVars::cons_bc, true);
        }
    }

    //
    // ***************************************************************************
    // We now interpolate the momenta on the "valid" faces of the fine grids (where
    // the interface is coarse/fine not fine/fine) -- this will not be over-written
    // below because the FillPatch operators see these as valid faces.
    //
    // Note that we interpolate momentum not velocity, but all the other boundary
    // conditions are imposed on velocity, so we convert to momentum here then
    // convert back.
    // ***************************************************************************
    if (fillset) {
        if (cf_set_width > 0) {
            FPr_c[lev-1].FillSet(*mfs_vel[Vars::cons], time_d, null_bc, domain_bcs_type);
        }
        if (cf_set_width >= 0 && !cons_only) {

            const MultiFab* c_vfrac = nullptr;
            if (solverChoice.terrain_type == TerrainType::EB) {
                c_vfrac = &((get_eb(lev).get_const_factory())->getVolFrac());
            }

            VelocityToMomentum(*mfs_vel[Vars::xvel], IntVect{0},
                               *mfs_vel[Vars::yvel], IntVect{0},
                               *mfs_vel[Vars::zvel], IntVect{0},
                               *mfs_vel[Vars::cons],
                               *mfs_mom[IntVars::xmom],
                               *mfs_mom[IntVars::ymom],
                               *mfs_mom[IntVars::zmom],
                               Geom(lev).Domain(),
                               domain_bcs_type, c_vfrac);

            FillSetMomentaCF(lev, time_d, *mfs_mom[IntVars::xmom], *mfs_mom[IntVars::ymom],
                             *mfs_mom[IntVars::zmom]);

            MomentumToVelocity(*mfs_vel[Vars::xvel], *mfs_vel[Vars::yvel], *mfs_vel[Vars::zvel],
                               *mfs_vel[Vars::cons],
                               *mfs_mom[IntVars::xmom],
                               *mfs_mom[IntVars::ymom],
                               *mfs_mom[IntVars::zmom],
                               Geom(lev).Domain(),
                               domain_bcs_type, c_vfrac);
        }
    }

    // We must fill a temporary then copy it back so we don't double add/subtract
    MultiFab mf_c(mfs_vel[Vars::cons]->boxArray(),mfs_vel[Vars::cons]->DistributionMap(),
                  mfs_vel[Vars::cons]->nComp()   ,mfs_vel[Vars::cons]->nGrowVect());

    mapper = &cell_cons_interp;

    if (interpolation_type == StateInterpType::Perturbational)
    {
        // Divide (rho theta) by rho to get theta (before we subtract rho0 from rho!)
        if (csel != TimeLevelSel::New) {
            MultiFab::Divide(vars_old[lev-1][Vars::cons],vars_old[lev-1][Vars::cons],
                             Rho_comp,RhoTheta_comp,1,ngvect_cons);
            MultiFab::Subtract(vars_old[lev-1][Vars::cons],base_state[lev-1],
                               BaseState::r0_comp,Rho_comp,1,ngvect_cons);
            MultiFab::Subtract(vars_old[lev-1][Vars::cons],base_state[lev-1],
                               BaseState::th0_comp,RhoTheta_comp,1,ngvect_cons);
        }
        if (csel != TimeLevelSel::Old) {
            MultiFab::Divide(vars_new[lev-1][Vars::cons],vars_new[lev-1][Vars::cons],
                             Rho_comp,RhoTheta_comp,1,ngvect_cons);
            MultiFab::Subtract(vars_new[lev-1][Vars::cons],base_state[lev-1],
                               BaseState::r0_comp,Rho_comp,1,ngvect_cons);
            MultiFab::Subtract(vars_new[lev-1][Vars::cons],base_state[lev-1],
                               BaseState::th0_comp,RhoTheta_comp,1,ngvect_cons);
        }

        if (fsel != TimeLevelSel::New) {
            MultiFab::Divide(vars_old[lev  ][Vars::cons],vars_old[lev  ][Vars::cons],
                             Rho_comp,RhoTheta_comp,1,IntVect{0});
            MultiFab::Subtract(vars_old[lev  ][Vars::cons],old_base_state,
                               BaseState::r0_comp,Rho_comp,1,IntVect{0});
            MultiFab::Subtract(vars_old[lev  ][Vars::cons],old_base_state,
                               BaseState::th0_comp,RhoTheta_comp,1,IntVect{0});
        }
        if (fsel != TimeLevelSel::Old) {
            MultiFab::Divide(vars_new[lev  ][Vars::cons],vars_new[lev  ][Vars::cons],
                             Rho_comp,RhoTheta_comp,1,IntVect{0});
            MultiFab::Subtract(vars_new[lev  ][Vars::cons],old_base_state,
                               BaseState::r0_comp,Rho_comp,1,IntVect{0});
            MultiFab::Subtract(vars_new[lev  ][Vars::cons],old_base_state,
                               BaseState::th0_comp,RhoTheta_comp,1,IntVect{0});
        }
    }

    // Call FillPatchTwoLevels which ASSUMES that all ghost cells at lev-1 have already been filled
    FillPatchTwoLevels(mf_c, ngvect_cons, IntVect(0,0,0),
                       time_rel, cmf, ctime, fmf, ftime,
                       0, 0, mf_c.nComp(), geom[lev-1], geom[lev],
                       refRatio(lev-1), mapper, domain_bcs_type,
                       BCVars::cons_bc);

    // FillPatchTwoLevels does not fill the cells outside the domain boundary, and all of
    // mf_c is copied into the state below; set them to a known value rather than leaving
    // uninitialized memory -- they will be filled in the physbcs call
    mf_c.setDomainBndry(bogus_large_value,geom[lev]);

    if (interpolation_type == StateInterpType::Perturbational)
    {
        // Restore the coarse values to what they were
        if (csel != TimeLevelSel::New) {
            MultiFab::Add(vars_old[lev-1][Vars::cons], base_state[lev-1],
                          BaseState::r0_comp,Rho_comp,1,ngvect_cons);
            MultiFab::Add(vars_old[lev-1][Vars::cons], base_state[lev-1],
                          BaseState::th0_comp,RhoTheta_comp,1,ngvect_cons);
            MultiFab::Multiply(vars_old[lev-1][Vars::cons], vars_old[lev-1][Vars::cons],
                               Rho_comp,RhoTheta_comp,1,ngvect_cons);
        }
        if (csel != TimeLevelSel::Old) {
            MultiFab::Add(vars_new[lev-1][Vars::cons], base_state[lev-1],
                          BaseState::r0_comp,Rho_comp,1,ngvect_cons);
            MultiFab::Add(vars_new[lev-1][Vars::cons], base_state[lev-1],
                          BaseState::th0_comp,RhoTheta_comp,1,ngvect_cons);
            MultiFab::Multiply(vars_new[lev-1][Vars::cons], vars_new[lev-1][Vars::cons],
                               Rho_comp,RhoTheta_comp,1,ngvect_cons);
        }

        if (fsel != TimeLevelSel::New) {
            MultiFab::Add(vars_old[lev][Vars::cons],base_state[lev  ],BaseState::r0_comp,Rho_comp,1,IntVect{0});
            MultiFab::Add(vars_old[lev][Vars::cons],base_state[lev  ],BaseState::th0_comp,RhoTheta_comp,1,IntVect{0});
            MultiFab::Multiply(vars_old[lev][Vars::cons], vars_old[lev][Vars::cons],
                               Rho_comp,RhoTheta_comp,1,IntVect{0});
        }
        if (fsel != TimeLevelSel::Old) {
            MultiFab::Add(vars_new[lev][Vars::cons], base_state[lev],BaseState::r0_comp,Rho_comp,1,IntVect{0});
            MultiFab::Add(vars_new[lev][Vars::cons], base_state[lev],BaseState::th0_comp,RhoTheta_comp,1,IntVect{0});
            MultiFab::Multiply(vars_new[lev][Vars::cons], vars_new[lev][Vars::cons],
                               Rho_comp,RhoTheta_comp,1,IntVect{0});
        }

        // Add rho_0 back to rho and theta_0 back to theta
        MultiFab::Add(mf_c, new_base_state,BaseState::r0_comp,Rho_comp,1,ngvect_cons);
        MultiFab::Add(mf_c, new_base_state,BaseState::th0_comp,RhoTheta_comp,1,ngvect_cons);

        // Multiply (theta) by rho to get (rho theta)
        MultiFab::Multiply(mf_c,mf_c,Rho_comp,RhoTheta_comp,1,ngvect_cons);
    }

    MultiFab::Copy(*mfs_vel[Vars::cons],mf_c,0,0,mf_c.nComp(),mf_c.nGrowVect());

    // ***************************************************************************************

    if (!cons_only)
    {
            mapper = &face_cons_linear_interp;

            MultiFab& mf_u = *mfs_vel[Vars::xvel];
            MultiFab& mf_v = *mfs_vel[Vars::yvel];
            MultiFab& mf_w = *mfs_vel[Vars::zvel];

            Vector<MultiFab*> fmf_u; Vector<MultiFab*> fmf_v; Vector<MultiFab*> fmf_w;
            Vector<MultiFab*> cmf_u; Vector<MultiFab*> cmf_v; Vector<MultiFab*> cmf_w;

            // **********************************************************************

            fmf_u = select_time_data(fsel, vars_old[lev  ][Vars::xvel], vars_new[lev  ][Vars::xvel]);
            fmf_v = select_time_data(fsel, vars_old[lev  ][Vars::yvel], vars_new[lev  ][Vars::yvel]);
            fmf_w = select_time_data(fsel, vars_old[lev  ][Vars::zvel], vars_new[lev  ][Vars::zvel]);
            cmf_u = select_time_data(csel, vars_old[lev-1][Vars::xvel], vars_new[lev-1][Vars::xvel]);
            cmf_v = select_time_data(csel, vars_old[lev-1][Vars::yvel], vars_new[lev-1][Vars::yvel]);
            cmf_w = select_time_data(csel, vars_old[lev-1][Vars::zvel], vars_new[lev-1][Vars::zvel]);

            // Call FillPatchTwoLevels which ASSUMES that all ghost cells at lev-1 have already been filled
            FillPatchTwoLevels(mf_u, ngvect_vels, IntVect(0,0,0),
                               time_rel, cmf_u, ctime, fmf_u, ftime,
                               0, 0, 1, geom[lev-1], geom[lev],
                               refRatio(lev-1), mapper, domain_bcs_type,
                               BCVars::xvel_bc);

            FillPatchTwoLevels(mf_v, ngvect_vels, IntVect(0,0,0),
                               time_rel, cmf_v, ctime, fmf_v, ftime,
                               0, 0, 1, geom[lev-1], geom[lev],
                               refRatio(lev-1), mapper, domain_bcs_type,
                               BCVars::yvel_bc);

            // We put these here because these may be used in constructing omega outside the
            // domain when fillpatching w
            bool do_fb = true;
            (*physbcs_u[lev])(*mfs_vel[Vars::xvel],*mfs_vel[Vars::xvel],*mfs_vel[Vars::yvel],
                              ngvect_vels,time,BCVars::xvel_bc, do_fb);
            (*physbcs_v[lev])(*mfs_vel[Vars::yvel],*mfs_vel[Vars::xvel],*mfs_vel[Vars::yvel],
                              ngvect_vels,time,BCVars::yvel_bc, do_fb);

            // **********************************************************************

            // Call FillPatchTwoLevels which ASSUMES that all ghost cells at lev-1 have already been filled
            FillPatchTwoLevels(mf_w, ngvect_vels, IntVect(0,0,0),
                               time_rel, cmf_w, ctime, fmf_w, ftime,
                               0, 0, 1, geom[lev-1], geom[lev],
                               refRatio(lev-1), mapper, domain_bcs_type,
                               BCVars::zvel_bc);
    } // !cons_only

    // ***************************************************************************
    // Physical bc's at domain boundary
    // ***************************************************************************
    int icomp_cons = 0;
    int ncomp_cons = mfs_vel[Vars::cons]->nComp();

    bool do_fb = true;

    if (m_r2d && !solverChoice.use_real_bcs) { fill_from_bndryregs(lev,mfs_vel,time_d); }

    // We call this even if use_real_bcs is true because these will fill the vertical bcs
    // Note that we call FillBoundary inside the physbcs call
    (*physbcs_cons[lev])(*mfs_vel[Vars::cons],*mfs_vel[Vars::xvel],*mfs_vel[Vars::yvel],
                         icomp_cons,ncomp_cons,ngvect_cons,time,BCVars::cons_bc, do_fb);
    if (!cons_only) {
        // Note that we need to fill u and v in the case of terrain because we will use
        //      these in the call of WFromOmega in lateral ghost cells of the fine grid
        // (*physbcs_u[lev])(*mfs_vel[Vars::xvel],*mfs_vel[Vars::xvel],*mfs_vel[Vars::yvel],
        //                   ngvect_vels,time,BCVars::xvel_bc, do_fb);
        // (*physbcs_v[lev])(*mfs_vel[Vars::yvel],*mfs_vel[Vars::xvel],*mfs_vel[Vars::yvel],
        //                   ngvect_vels,time,BCVars::yvel_bc, do_fb);
        (*physbcs_w[lev])(*mfs_vel[Vars::zvel],*mfs_vel[Vars::xvel],*mfs_vel[Vars::yvel],
                          ngvect_vels,time,BCVars::zvel_bc, do_fb);
    }
}

/**
 * Fill valid and ghost data on the coarse level with state data at the given time.
 *
 * @param[in]     lev       coarse level to fill
 * @param[in]     time_d    time at which the data should be filled
 * @param[in,out] mfs_vel   Vector of MultiFabs to be filled containing cons, xvel, yvel, and zvel
 * @param[in]     cons_only whether to fill only conserved variables
 */
void
ERF::FillPatchCrseLevel (int lev, double time_d,
                         const Vector<MultiFab*>& mfs_vel,     // This includes cc quantities and VELOCITIES
                         bool cons_only)
{
    BL_PROFILE_VAR("ERF::FillPatchCrseLevel()",ERF_FillPatchCrseLevel);

    AMREX_ALWAYS_ASSERT(lev == 0);

    Real time = static_cast<Real>(time_d);

    IntVect ngvect_cons = mfs_vel[Vars::cons]->nGrowVect();
    IntVect ngvect_vels = mfs_vel[Vars::xvel]->nGrowVect();

    const TimeLevelSel fsel = select_time_level(time_d, t_old[lev], t_new[lev]);

    // As in FillPatchFineLevel, hand AMReX times measured from the new time, formed in double
    const double tref = t_new[lev];
    const Real time_rel = static_cast<Real>(time_d - tref);
    Vector<Real> ftime  = {static_cast<Real>(t_old[lev] - tref), Real(0)};

    //
    // Below we call FillPatchSingleLevel which does NOT fill ghost cells outside the domain
    //

    Vector<MultiFab*> fmf;
    Vector<MultiFab*> fmf_u;
    Vector<MultiFab*> fmf_v;
    Vector<MultiFab*> fmf_w;

    fmf = select_time_data(fsel, vars_old[lev][Vars::cons], vars_new[lev][Vars::cons]);

    const int  ncomp = mfs_vel[Vars::cons]->nComp();

    FillPatchSingleLevel(*mfs_vel[Vars::cons], ngvect_cons, time_rel, fmf, IntVect(0,0,0), ftime,
                         0, 0, ncomp, geom[lev]);

    if (!cons_only) {
        fmf_u = select_time_data(fsel, vars_old[lev][Vars::xvel], vars_new[lev][Vars::xvel]);
        fmf_v = select_time_data(fsel, vars_old[lev][Vars::yvel], vars_new[lev][Vars::yvel]);
        fmf_w = select_time_data(fsel, vars_old[lev][Vars::zvel], vars_new[lev][Vars::zvel]);
        FillPatchSingleLevel(*mfs_vel[Vars::xvel], ngvect_vels, time_rel, fmf_u,
                             IntVect(0,0,0), ftime,  0, 0, 1, geom[lev]);

        FillPatchSingleLevel(*mfs_vel[Vars::yvel], ngvect_vels, time_rel, fmf_v,
                             IntVect(0,0,0), ftime,  0, 0, 1, geom[lev]);

        FillPatchSingleLevel(*mfs_vel[Vars::zvel], ngvect_vels, time_rel, fmf_w,
                             IntVect(0,0,0), ftime,  0, 0, 1, geom[lev]);
    } // !cons_only

    // ***************************************************************************
    // Physical bc's at domain boundary
    // ***************************************************************************
    int icomp_cons = 0;
    int ncomp_cons = mfs_vel[Vars::cons]->nComp();

    bool do_fb = true;

#ifdef ERF_USE_NETCDF
    if(solverChoice.use_real_bcs && (lev==0)) {
        fill_from_realbdy(mfs_vel,time_d,cons_only,icomp_cons,ncomp_cons,ngvect_cons,ngvect_vels);
        do_fb = false;
    }
#endif

    if (m_r2d && !solverChoice.use_real_bcs) { fill_from_bndryregs(lev,mfs_vel,time_d); }

    // We call this even if use_real_bcs is true because these will fill the vertical bcs
    // Note that we call FillBoundary inside the physbcs call
    (*physbcs_cons[lev])(*mfs_vel[Vars::cons],*mfs_vel[Vars::xvel],*mfs_vel[Vars::yvel],
                          icomp_cons,ncomp_cons,ngvect_cons,time,BCVars::cons_bc, do_fb);
    if (!cons_only) {
        (*physbcs_u[lev])(*mfs_vel[Vars::xvel],*mfs_vel[Vars::xvel],*mfs_vel[Vars::yvel],
                          ngvect_vels,time,BCVars::xvel_bc, do_fb);
        (*physbcs_v[lev])(*mfs_vel[Vars::yvel],*mfs_vel[Vars::xvel],*mfs_vel[Vars::yvel],
                          ngvect_vels,time,BCVars::yvel_bc, do_fb);
        (*physbcs_w[lev])(*mfs_vel[Vars::zvel],*mfs_vel[Vars::xvel],*mfs_vel[Vars::yvel],
                          ngvect_vels,time,BCVars::zvel_bc, do_fb);
    }
}
