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
 * @param[in]     fillset        whether to fill the coarse-fine set region before standard fillpatch
 * @param[in]     cons_only      whether to fill only conserved variables
 * @param[in]     allow_perturbational false to fill the conserved state from coarse in full
 *                               whatever erf.perturbational_interpolation, as for a destination
 *                               that is not on the layout of base_state[lev]
 */
void
ERF::FillPatchFineLevel (int lev, double time_d,
                         const Vector<MultiFab*>& mfs_vel,     // This includes cc quantities and VELOCITIES
                         const Vector<MultiFab*>& mfs_mom,     // This includes cc quantities and MOMENTA
                         bool fillset, bool cons_only, bool allow_perturbational)
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
    // wider would be thrown away.
    // ***************************************************************************
    if (fillset && !cons_only) {
        {
            const IntVect ng_rho(1,1,1);
            MultiFab rho_mf(mfs_vel[Vars::cons]->boxArray(),
                            mfs_vel[Vars::cons]->DistributionMap(), 1, ng_rho);
            FillPatchTwoLevelsCons(lev, rho_mf, ng_rho, time_rel, cmf, ctime, fmf, ftime,
                                   Rho_comp, 0, 1, BCVars::cons_bc, allow_perturbational);
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

    // Fill a temporary, then copy it back: the destination is often vars_new[lev] itself,
    // which is also one of the fine-level sources of the fill
    MultiFab mf_c(mfs_vel[Vars::cons]->boxArray(),mfs_vel[Vars::cons]->DistributionMap(),
                  mfs_vel[Vars::cons]->nComp()   ,mfs_vel[Vars::cons]->nGrowVect());

    // This ASSUMES that all ghost cells at lev-1 have already been filled
    FillPatchTwoLevelsCons(lev, mf_c, ngvect_cons, time_rel, cmf, ctime, fmf, ftime,
                           0, 0, mf_c.nComp(), BCVars::cons_bc, allow_perturbational);

    // FillPatchTwoLevels does not fill the cells outside the domain boundary, and all of
    // mf_c is copied into the state below; set them to a known value rather than leaving
    // uninitialized memory -- they will be filled in the physbcs call
    mf_c.setDomainBndry(bogus_large_value,geom[lev]);

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

namespace {
/**
 * Add sign * rho0 to component ir, sign * rho0 theta0 to component irt and sign * rho0 qv0 to
 * component iqv of mf (each may be -1, for absent), in the cells of mf's boxes grown by ng
 * that lie in region.  The base state bs may be on any layout of mf's level; every cell
 * touched must lie in its valid or ghost region.
 */
void
add_base_cons (MultiFab& mf, const IntVect& ng, const Box& region, const MultiFab& bs,
               const Periodicity& period, Real sign, int ir, int irt, int iqv)
{
    // The base state on mf's layout, starting from a value that cannot pass unnoticed
    //    in a cell bs does not cover
    MultiFab b(mf.boxArray(), mf.DistributionMap(), 3, ng);
    b.setVal(bogus_large_value);
    b.ParallelCopy(bs, BaseState::r0_comp , 0, 1, bs.nGrowVect(), ng, period);
    b.ParallelCopy(bs, BaseState::th0_comp, 1, 1, bs.nGrowVect(), ng, period);
    if (iqv >= 0) {
        b.ParallelCopy(bs, BaseState::qv0_comp, 2, 1, bs.nGrowVect(), ng, period);
    }

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box bx = amrex::grow(mfi.validbox(), ng) & region;
        if (bx.isEmpty()) { continue; }
        const auto& m_arr = mf.array(mfi);
        const auto& b_arr = b.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real r0 = b_arr(i,j,k,0);
            if (ir  >= 0) { m_arr(i,j,k,ir ) += sign * r0; }
            if (irt >= 0) { m_arr(i,j,k,irt) += sign * r0 * b_arr(i,j,k,1); }
            if (iqv >= 0) { m_arr(i,j,k,iqv) += sign * r0 * b_arr(i,j,k,2); }
        });
    }
}
} // namespace

/**
 * FillPatchTwoLevels for the cell-centered conserved state on level lev.
 *
 * With erf.perturbational_interpolation, the cells filled from the coarse level take
 * rho0 + I(rho - rho0), rho0 theta0 + I(rho theta - rho0 theta0) and, with moisture,
 * rho0 qv0 + I(rho qv - rho0 qv0), with I the
 * coarse-to-fine interpolation and rho0, theta0, qv0 the base state of the level the
 * data lives on.  An atmosphere that is in its base state on the coarse level is
 * then in the fine base state in those cells, whatever the two levels' heights.
 * The fine-level data, including the fine/fine ghost cells, is copied as is.
 *
 * @param[in]     lev    fine level
 * @param[in,out] mf     destination, on level lev
 * @param[in]     nghost number of ghost cells of mf to fill
 * @param[in]     time   time of the fill, measured as ctime and ftime are
 * @param[in]     cmf    coarse source data at ctime
 * @param[in]     fmf    fine source data at ftime
 * @param[in]     scomp  first source component
 * @param[in]     dcomp  first destination component
 * @param[in]     ncomp  number of components
 * @param[in]     bccomp first component of domain_bcs_type to use
 * @param[in]     perturbational false for a full-state fill whatever erf.perturbational_interpolation
 */
void
ERF::FillPatchTwoLevelsCons (int lev, MultiFab& mf, const IntVect& nghost, Real time,
                             const Vector<MultiFab*>& cmf, const Vector<Real>& ctime,
                             const Vector<MultiFab*>& fmf, const Vector<Real>& ftime,
                             int scomp, int dcomp, int ncomp, int bccomp,
                             bool perturbational)
{
    Interpolater* mapper = &cell_cons_interp;

    // Where Rho_comp and RhoTheta_comp sit in the patches the hooks see, which hold
    //    the source components scomp to scomp+ncomp-1 starting at 0
    const int ir  = Rho_comp      - scomp;
    const int irt = RhoTheta_comp - scomp;
    const int iqv = RhoQ1_comp    - scomp;
    const bool l_moist = (solverChoice.moisture_type != MoistureType::None);
    const bool has_r  = (ir  >= 0 && ir  < ncomp);
    const bool has_rt = (irt >= 0 && irt < ncomp);
    const bool has_qv = l_moist && (iqv >= 0 && iqv < ncomp);

    if (!perturbational_interp || !perturbational || !(has_r || has_rt || has_qv)) {
        FillPatchTwoLevels(mf, nghost, IntVect(0,0,0), time, cmf, ctime, fmf, ftime,
                           scomp, dcomp, ncomp, geom[lev-1], geom[lev],
                           refRatio(lev-1), mapper, domain_bcs_type, bccomp);
        return;
    }

    const MultiFab& bs_crse = base_state[lev-1];
    const MultiFab& bs_fine = base_state[lev  ];
    const Geometry& cgeom   = geom[lev-1];
    const Geometry& fgeom   = geom[lev  ];

    const int ir_p  = (has_r ) ? ir  : -1;
    const int irt_p = (has_rt) ? irt : -1;
    const int iqv_p = (has_qv) ? iqv : -1;

    // The patches have no ghost cells, and the boxes of each lie in the base state's
    //    valid or ghost region
    const Box everywhere(IntVect(std::numeric_limits<int>::lowest()/2),
                         IntVect(std::numeric_limits<int>::max()/2));
    auto pre_interp  = [&] (MultiFab& crse_patch, int, int) {
        add_base_cons(crse_patch, IntVect(0), everywhere, bs_crse, cgeom.periodicity(), -one, ir_p, irt_p, iqv_p);
    };
    auto post_interp = [&] (MultiFab& fine_patch, int, int) {
        add_base_cons(fine_patch, IntVect(0), everywhere, bs_fine, fgeom.periodicity(),  one, ir_p, irt_p, iqv_p);
    };

    PhysBCFunctUseCoarseGhost erfbc(*cmf[0], nghost, IntVect(0,0,0), refRatio(lev-1), mapper);
    FillPatchTwoLevels(mf, nghost, time, cmf, ctime, fmf, ftime, scomp, dcomp, ncomp,
                       cgeom, fgeom, erfbc, 0, erfbc, 0, refRatio(lev-1), mapper,
                       domain_bcs_type, bccomp, pre_interp, post_interp);
}

/**
 * InterpFromCoarseLevel of the whole cell-centered conserved state of level lev from the
 * coarser level, into the valid cells and the ghost cells inside the domain: as
 * FillPatchTwoLevelsCons, perturbational in Rho_comp, RhoTheta_comp and RhoQ1_comp with
 * erf.perturbational_interpolation.  For a level made from the coarser one, so that mf must
 * be on the layout of base_state[lev], which must already be built on it.
 *
 * @param[in]     lev    fine level
 * @param[in,out] mf     destination, on level lev
 * @param[in]     nghost number of ghost cells of mf to fill
 * @param[in]     cmf    coarse source data
 */
void
ERF::InterpFromCoarseLevelCons (int lev, MultiFab& mf, const IntVect& nghost, const MultiFab& cmf)
{
    Interpolater* mapper = &cell_cons_interp;
    const int ncomp = mf.nComp();
    const int iqv = (solverChoice.moisture_type != MoistureType::None) ? RhoQ1_comp : -1;

    if (!perturbational_interp) {
        InterpFromCoarseLevel(mf, nghost, IntVect(0,0,0), cmf, 0, 0, ncomp,
                              geom[lev-1], geom[lev], refRatio(lev-1), mapper,
                              domain_bcs_type, BCVars::cons_bc);
        return;
    }

    AMREX_ALWAYS_ASSERT(mf.boxArray() == base_state[lev].boxArray() &&
                        mf.DistributionMap() == base_state[lev].DistributionMap());

    const Geometry& cgeom = geom[lev-1];
    const Geometry& fgeom = geom[lev  ];
    const MultiFab& bs_crse = base_state[lev-1];

    const Box everywhere(IntVect(std::numeric_limits<int>::lowest()/2),
                         IntVect(std::numeric_limits<int>::max()/2));
    // With PhysBCFunctUseCoarseGhost, InterpFromCoarseLevel's coarse patch is coarsen(ba) with
    //    ghost cells, which the interpolation reads, so take the base state out of those too:
    //    left as the full state, the slopes across the edge of each patch mix it with the
    //    perturbation, and the fine cells near every box edge -- including the edges between
    //    a level's boxes -- then depend on how the level is chopped
    auto pre_interp = [&] (MultiFab& crse_patch, int, int) {
        add_base_cons(crse_patch, crse_patch.nGrowVect(), everywhere, bs_crse, cgeom.periodicity(),
                      -one, Rho_comp, RhoTheta_comp, iqv);
    };
    // InterpFromCoarseLevel's post hook sees one fab at a time, without its index, so the
    //    fine base state is added below instead
    auto post_interp = [] (FArrayBox&, const Box&, int, int) {};

    PhysBCFunctUseCoarseGhost erfbc(cmf, nghost, IntVect(0,0,0), refRatio(lev-1), mapper);
    InterpFromCoarseLevel(mf, nghost, Real(0.0), cmf, 0, 0, ncomp, cgeom, fgeom,
                          erfbc, 0, erfbc, 0, refRatio(lev-1), mapper,
                          domain_bcs_type, BCVars::cons_bc, pre_interp, post_interp);

    // The cells the interpolation filled: inside the domain, which extends across a
    //    periodic boundary
    Box fdomain_g(fgeom.Domain());
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (fgeom.isPeriodic(d)) { fdomain_g.grow(d, nghost[d]); }
    }
    add_base_cons(mf, nghost, fdomain_g, base_state[lev], fgeom.periodicity(),
                  one, Rho_comp, RhoTheta_comp, iqv);
}
