#include <ERF.H>
#include "ERF_Constants.H"
#include <ERF_PhysBCFunct.H>
#include <ERF_IndexDefines.H>
#include <ERF_TimeInterpolatedData.H>
#include <ERF_FillPatcher.H>
#include <ERF_Utils.H>

#include <algorithm>

using namespace amrex;

/**
 * Fill valid and ghost data
 * This version fills mfs in valid regions with the values in "mfs" when it is passed in;
 * it is used only to compute ghost values for intermediate stages of a time integrator.
 *
 * @param[in]     lev         level of refinement at which to fill the data
 * @param[in]     time_d      time at which the data should be filled
 * @param[in,out] mfs_vel     Vector of MultiFabs to be filled containing, in order: cons, xvel, yvel, and zvel
 * @param[in,out] mfs_mom     Vector of MultiFabs to be filled containing, in order: cons, xmom, ymom, and zmom
 * @param[in]     ng_cons     number of ghost cells to be filled for conserved variables
 * @param[in]     ng_vel      number of ghost cells to be filled for velocity components
 * @param[in]     cons_only   if true then only fill conserved variables
 * @param[in]     icomp_cons  starting component for conserved variables
 * @param[in]     ncomp_cons  number of components for conserved variables
 */
void
ERF::FillIntermediatePatch (int lev, double time_d,
                            const Vector<MultiFab*>& mfs_vel,     // This includes cc quantities and VELOCITIES
                            const Vector<MultiFab*>& mfs_mom,     // This includes cc quantities and MOMENTA
                            int ng_cons, int ng_vel, bool cons_only,
                            int icomp_cons, int ncomp_cons)
{
    BL_PROFILE_VAR("FillIntermediatePatch()",FillIntermediatePatch);
    Interpolater* mapper;

    Real time = static_cast<Real>(time_d);

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
            FPr_c[lev-1].FillSet(*mfs_vel[Vars::cons], time_d, null_bc, domain_bcs_type);
        }
        if ( !cons_only && (cf_set_width >= 0) ) {
            // The coarse/fine faces from the coarse level, those in z of a level that stops
            // inside the domain with the coarse Omega (FillSetMomentaCF)
            FillSetMomentaCF(lev, time_d, *mfs_mom[IntVars::xmom], *mfs_mom[IntVars::ymom],
                             *mfs_mom[IntVars::zmom]);
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
            const MultiFab* c_vfrac = nullptr;
            if (solverChoice.terrain_type == TerrainType::EB) {
                c_vfrac = &((get_eb(lev).get_const_factory())->getVolFrac());
            }

            // This only fills VALID region of velocity
            MomentumToVelocity(*mfs_vel[Vars::xvel], *mfs_vel[Vars::yvel], *mfs_vel[Vars::zvel],
                               *mfs_vel[Vars::cons],
                               *mfs_mom[IntVars::xmom], *mfs_mom[IntVars::ymom], *mfs_mom[IntVars::zmom],
                                Geom(lev).Domain(), domain_bcs_type, c_vfrac);
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
        // Set all components to bogus_large_value, then copy just the density from *mfs_vel[Vars::cons]
        //
        mf.setVal(bogus_large_value);
        MultiFab::Copy(mf,*mfs_vel[Vars::cons],Rho_comp,Rho_comp,1,mf.nGrowVect());

        // Which coarse time level(s) to take is decided in double precision
        const TimeLevelSel csel = select_time_level(time_d, t_old[lev-1], t_new[lev-1]);

        // AMReX computes its time-interpolation weights in Real from the times it is given,
        // so give it times measured from the coarse new time (formed in double) rather than
        // the elapsed times themselves, which can be large compared to dt.  The physical
        // bc's below still see the elapsed time.
        const double tref = t_new[lev-1];
        const Real time_rel = static_cast<Real>(time_d - tref);
        Vector<Real> ctime  = {static_cast<Real>(t_old[lev-1] - tref), Real(0)};
        Vector<Real> ftime  = {time_rel, time_rel};

        Vector<MultiFab*> fmf = {mfs_vel[Vars::cons],mfs_vel[Vars::cons]};
        Vector<MultiFab*> cmf = select_time_data(csel, vars_old[lev-1][Vars::cons], vars_new[lev-1][Vars::cons]);

        if (interpolation_type == StateInterpType::Perturbational)
        {
            if (icomp_cons+ncomp_cons > 1)
            {
                // Divide (rho theta) by rho to get theta
                MultiFab::Divide(*mfs_vel[Vars::cons],*mfs_vel[Vars::cons],Rho_comp,RhoTheta_comp,1,IntVect{0});

                // Subtract theta_0 from theta
                MultiFab::Subtract(*mfs_vel[Vars::cons],base_state[lev],BaseState::th0_comp,RhoTheta_comp,1,IntVect{0});

                if (csel != TimeLevelSel::New) {
                    MultiFab::Divide(vars_old[lev-1][Vars::cons], vars_old[lev-1][Vars::cons],
                                     Rho_comp,RhoTheta_comp,1,vars_old[lev-1][Vars::cons].nGrowVect());
                    MultiFab::Subtract(vars_old[lev-1][Vars::cons], base_state[lev-1],
                                       BaseState::th0_comp,RhoTheta_comp,1,vars_old[lev-1][Vars::cons].nGrowVect());
                }
                if (csel != TimeLevelSel::Old) {
                    MultiFab::Divide(vars_new[lev-1][Vars::cons], vars_new[lev-1][Vars::cons],
                                     Rho_comp,RhoTheta_comp,1,vars_new[lev-1][Vars::cons].nGrowVect());
                    MultiFab::Subtract(vars_new[lev-1][Vars::cons], base_state[lev-1],
                                       BaseState::th0_comp,RhoTheta_comp,1,vars_new[lev-1][Vars::cons].nGrowVect());
                }
            }

            // Subtract rho_0 from rho before we interpolate -- note we only subtract
            //    on valid region of mf since the ghost cells will be filled below
            if (icomp_cons == 0)
            {
                MultiFab::Subtract(*mfs_vel[Vars::cons],base_state[lev],BaseState::r0_comp,Rho_comp,1,IntVect{0});

                if (csel != TimeLevelSel::New) {
                    MultiFab::Subtract(vars_old[lev-1][Vars::cons], base_state[lev-1],
                                       BaseState::r0_comp,Rho_comp,1,vars_old[lev-1][Vars::cons].nGrowVect());
                }
                if (csel != TimeLevelSel::Old) {
                    MultiFab::Subtract(vars_new[lev-1][Vars::cons], base_state[lev-1],
                                       BaseState::r0_comp,Rho_comp,1,vars_new[lev-1][Vars::cons].nGrowVect());
                }
            }
        } // interpolation_type == StateInterpType::Perturbational

        // Call FillPatchTwoLevels which ASSUMES that all ghost cells at lev-1 have already been filled
        mapper = &cell_cons_interp;
        FillPatchTwoLevels(mf, IntVect{ng_cons}, IntVect(0,0,0),
                           time_rel, cmf, ctime, fmf, ftime,
                           icomp_cons, icomp_cons, ncomp_cons, geom[lev-1], geom[lev],
                           refRatio(lev-1), mapper, domain_bcs_type,
                           icomp_cons);

        if (interpolation_type == StateInterpType::Perturbational)
        {
            if (icomp_cons == 0)
            {
                // Restore the coarse values to what they were
                if (csel != TimeLevelSel::New) {
                    MultiFab::Add(vars_old[lev-1][Vars::cons], base_state[lev-1],
                                  BaseState::r0_comp,Rho_comp,1,vars_old[lev-1][Vars::cons].nGrowVect());
                }
                if (csel != TimeLevelSel::Old) {
                    MultiFab::Add(vars_new[lev-1][Vars::cons], base_state[lev-1],
                                  BaseState::r0_comp,Rho_comp,1,vars_new[lev-1][Vars::cons].nGrowVect());
                }

                // Set values in the cells outside the domain boundary so that we can do the Add
                //     without worrying about uninitialized values outside the domain -- these
                //     will be filled in the physbcs call
                mf.setDomainBndry(bogus_large_value,Rho_comp,1,geom[lev]);

                // Add rho_0 back to rho after we interpolate -- on all the valid + ghost region
                MultiFab::Add(mf, base_state[lev],BaseState::r0_comp,Rho_comp,1,IntVect{ng_cons});
            }

            if (icomp_cons+ncomp_cons > 1)
            {
                // Add theta_0 to theta
                if (csel != TimeLevelSel::New) {
                    MultiFab::Add(vars_old[lev-1][Vars::cons], base_state[lev-1],
                                  BaseState::th0_comp,RhoTheta_comp,1,vars_old[lev-1][Vars::cons].nGrowVect());
                    MultiFab::Multiply(vars_old[lev-1][Vars::cons], vars_old[lev-1][Vars::cons],
                                       Rho_comp,RhoTheta_comp,1,vars_old[lev-1][Vars::cons].nGrowVect());
                }
                if (csel != TimeLevelSel::Old) {
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
        } // interpolation_type == StateInterpType::Perturbational

        // Impose physical bc's on fine data (note time and 0 are not used)
        // Note that we do this after the FillPatch because imposing physical bc's on fine ghost
        //      cells that need to be filled from coarse requires that we have done the interpolation first
        // Quantities ingested from boundary planes are not imposed by physbcs_cons; the density
        // among them is also read by MomentumToVelocity below at the ingested faces
        if (m_r2d && !solverChoice.use_real_bcs) {
            fill_from_bndryregs(lev, {&mf, mfs_vel[Vars::xvel], mfs_vel[Vars::yvel], mfs_vel[Vars::zvel]},
                                time_d, true);
        }
        bool do_fb = true; bool do_terrain_adjustment = false;
        (*physbcs_cons[lev])(mf,*mfs_vel[Vars::xvel],*mfs_vel[Vars::yvel],
                             icomp_cons,ncomp_cons,IntVect{ng_cons},time,BCVars::cons_bc,
                             do_fb, do_terrain_adjustment);

        // Make sure to only copy back the components we worked on
        MultiFab::Copy(*mfs_vel[Vars::cons],mf,icomp_cons,icomp_cons,ncomp_cons,IntVect{ng_cons});

        // *****************************************************************************************

        if (!cons_only)
        {
            // ***************************************************************************
            // We always come in to this call with updated momenta but we need to create updated velocity
            //    in order to impose the rest of the bc's
            // ***************************************************************************
            const MultiFab* c_vfrac = nullptr;
            if (solverChoice.terrain_type == TerrainType::EB) {
                c_vfrac = &((get_eb(lev).get_const_factory())->getVolFrac());
            }

            // This only fills VALID region of velocity
            MomentumToVelocity(*mfs_vel[Vars::xvel], *mfs_vel[Vars::yvel], *mfs_vel[Vars::zvel],
                               *mfs_vel[Vars::cons],
                               *mfs_mom[IntVars::xmom], *mfs_mom[IntVars::ymom], *mfs_mom[IntVars::zmom],
                                Geom(lev).Domain(), domain_bcs_type, c_vfrac);

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
            cmf = select_time_data(csel, vars_old[lev-1][Vars::xvel], vars_new[lev-1][Vars::xvel]);

            // Call FillPatchTwoLevels which ASSUMES that all ghost cells at lev-1 have already been filled
            FillPatchTwoLevels(mfu, IntVect{ng_vel}, IntVect(0,0,0),
                               time_rel, cmf, ctime, fmf, ftime,
                               0, 0, 1, geom[lev-1], geom[lev],
                               refRatio(lev-1), mapper, domain_bcs_type,
                               BCVars::xvel_bc);

            // *****************************************************************************************

            MultiFab& mfv = *mfs_vel[Vars::yvel];

            fmf = {&mfv,&mfv};
            cmf = select_time_data(csel, vars_old[lev-1][Vars::yvel], vars_new[lev-1][Vars::yvel]);

            // Call FillPatchTwoLevels which ASSUMES that all ghost cells at lev-1 have already been filled
            FillPatchTwoLevels(mfv, IntVect{ng_vel}, IntVect(0,0,0),
                               time_rel, cmf, ctime, fmf, ftime,
                               0, 0, 1, geom[lev-1], geom[lev],
                               refRatio(lev-1), mapper, domain_bcs_type,
                               BCVars::yvel_bc);

            // *****************************************************************************************

            MultiFab& mfw = *mfs_vel[Vars::zvel];

            fmf = {&mfw,&mfw};
            cmf = select_time_data(csel, vars_old[lev-1][Vars::zvel], vars_new[lev-1][Vars::zvel]);

            // Call FillPatchTwoLevels which ASSUMES that all ghost cells at lev-1 have already been filled
            FillPatchTwoLevels(mfw, IntVect{ng_vel}, IntVect(0,0,0),
                               time_rel, cmf, ctime, fmf, ftime,
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
    if (solverChoice.use_real_bcs && (lev==0)) {
        fill_from_realbdy(mfs_vel,time_d,cons_only,icomp_cons,ncomp_cons,ngvect_cons,ngvect_vels);
        do_fb = false;
    }
#endif

    if (m_r2d && !solverChoice.use_real_bcs) { fill_from_bndryregs(lev,mfs_vel,time_d); }

    // We call this even if use_real_bcs is true because these will fill the vertical bcs
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
    // ***************************************************************************

    // We always come in to this call with momenta so we need to leave with momenta!
    // We need to make sure we convert back on all ghost cells/faces because this is
    // how velocity from fine-fine copies (as well as physical and interpolated bcs) will be filled
    if (!cons_only)
    {
        // One ghost cell of momentum is enough for the integrator unless we are using numerical
        // diffusion or EB -- the EB momentum fluxes and Omega are built on grown boxes and read
        // momenta two cells outside the valid region
        const bool all_mom_ghosts = solverChoice.use_num_diff ||
                                    (solverChoice.terrain_type == TerrainType::EB);
        IntVect ngu = (!all_mom_ghosts) ? IntVect(1,1,1) : mfs_vel[Vars::xvel]->nGrowVect();
        IntVect ngv = (!all_mom_ghosts) ? IntVect(1,1,1) : mfs_vel[Vars::yvel]->nGrowVect();
        IntVect ngw = (!all_mom_ghosts) ? IntVect(1,1,0) : mfs_vel[Vars::zvel]->nGrowVect();

        const MultiFab* c_vfrac = nullptr;
        if (solverChoice.terrain_type == TerrainType::EB) {
            c_vfrac = &((get_eb(lev).get_const_factory())->getVolFrac());
        }

        VelocityToMomentum(*mfs_vel[Vars::xvel], ngu,
                           *mfs_vel[Vars::yvel], ngv,
                           *mfs_vel[Vars::zvel], ngw,
                           *mfs_vel[Vars::cons],
                           *mfs_mom[IntVars::xmom], *mfs_mom[IntVars::ymom], *mfs_mom[IntVars::zmom],
                           Geom(lev).Domain(),
                           domain_bcs_type, c_vfrac);
    }

    // NOTE: There are not FillBoundary calls here for the following reasons:
    // Removal of the FillBoundary (FB) calls has bee completed for the following reasons:
    //
    // 1. physbc_cons is called before VelocityToMomentum and a FB is completed in that functor.
    //    Therefore, the conserved CC vars have their inter-rank ghost cells filled and then their
    //    domain ghost cells filled from the BC operations. We should not call FB on this MF again.
    //
    // 2. physbc_u/v/w is also called before VelocityToMomentum and a FB is completed those functors.
    //    Furthermore, VelocityToMomentum operates on a growntilebox so we exit that routine with momentum
    //    filled everywhere---i.e., physbc_u/v/w fills velocity ghost cells (inter-rank and domain)
    //    and then V2M does the conversion to momenta everywhere; so there is again no need to do a FB on momenta.
}

/**
 * Mask of the coarse/fine ghost cells of level lev (1 there, 0 in the valid region, in ghost
 * cells covered by another grid at this level and outside the domain), for the cell-centred
 * (which = 0) or face-centred (which = 1 + dir) data of FillCoarseFineGhostsForImplicit.  It
 * depends only on the grids, so it is kept and rebuilt only when they change.
 */
const iMultiFab&
ERF::get_cf_implicit_mask (int lev, int which, const BoxArray& ba,
                           const DistributionMapping& dm, const IntVect& ng)
{
    std::unique_ptr<iMultiFab>& mask = cf_implicit_mask[lev][which];
    if (!mask || mask->boxArray() != ba || mask->DistributionMap() != dm || mask->nGrowVect() != ng)
    {
        mask = std::make_unique<iMultiFab>(ba, dm, 1, ng);
        mask->BuildMask(geom[lev].Domain(), geom[lev].periodicity(), 0, 1, 0, 0);
    }
    return *mask;
}

/**
 * Fill the coarse/fine ghost cells of the state used by an implicit vertical solve with coarse
 * data at the time that solve advances to.
 *
 * The solves are done on a state whose valid region has been advanced to the new stage time
 * but whose ghost cells have not: the scratch state of ERF_TI_slow_rhs_pre.H (S_old + dt * S_rhs,
 * with the ghost values of S_old), and S_new in ERF_TI_slow_rhs_post.H (updated on the valid
 * region only, with the ghost values of the previous stage).  A column whose top or bottom is a
 * coarse/fine boundary closes its tridiagonal system with the value on the other side of that
 * boundary, which should be at the new stage time rather than lagged.
 *
 * Only ghost cells that are neither covered by another grid at this level nor outside the
 * domain are changed: those between grids at this level are filled by FillBoundary, and those
 * outside the domain carry the boundary data the physical-BC rows of the solve read.
 *
 * Only the ghost cells the solves read are filled; the others keep their old values, which
 * nothing reads before the state is next filled (the scratch state of slow_rhs_pre is used on
 * its valid region only, and apply_bcs fills all the ghost cells of S_new at the end of
 * slow_rhs_post).  Every column is closed in z, so a solved component is read at the cell above
 * and below the column, but not at its sides; the same holds for the momenta.  Density is read
 * at those cells too, and also by the face averages: at the neighbouring cell in x (y) all the
 * way up the column for the x- (y-) momentum, and at two cells above and below for the
 * z-momentum.
 *
 * The momenta are formed from the interpolated velocity and the refreshed density with the
 * same face average the solve divides by, so that it recovers the interpolated velocity.
 *
 * @param[in]    lev        level (must be > 0)
 * @param[in]    time_d     time the implicit solve advances to
 * @param[inout] cons       conserved state, with rho in Rho_comp and the solved components
 *                          at their usual indices
 * @param[in]    cons_comps components solved for, besides rho (which is always filled since
 *                          the solve divides by it)
 * @param[inout] mom        x-, y- and z-momenta (any entry may be nullptr)
 */
void
ERF::FillCoarseFineGhostsForImplicit (int lev, double time_d,
                                      MultiFab& cons,
                                      const Vector<int>& cons_comps,
                                      const Vector<MultiFab*>& mom)
{
    BL_PROFILE("FillCoarseFineGhostsForImplicit()");
    AMREX_ALWAYS_ASSERT(lev > 0);
    AMREX_ALWAYS_ASSERT(interpolation_type == StateInterpType::FullState);

    const TimeLevelSel csel = select_time_level(time_d, t_old[lev-1], t_new[lev-1]);

    // As in FillIntermediatePatch, give AMReX times measured from the coarse new time
    const double tref = t_new[lev-1];
    const Real time_rel = static_cast<Real>(time_d - tref);
    Vector<Real> ctime  = {static_cast<Real>(t_old[lev-1] - tref), Real(0)};
    Vector<Real> ftime  = {time_rel, time_rel};

    // Copy components [scomp, scomp+ncomp) of src into dst on the cells where mask is set
    auto copy_cf_ghosts = [&] (MultiFab& dst, const MultiFab& src, const iMultiFab& mask,
                               int scomp, int ncomp, const IntVect& ng)
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box gbx = mfi.growntilebox(ng);
            const Array4<const int > m = mask.const_array(mfi);
            const Array4<const Real> s = src.const_array(mfi);
            const Array4<      Real> d = dst.array(mfi);
            ParallelFor(gbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (m(i,j,k)) { d(i,j,k,scomp+n) = s(i,j,k,scomp+n); }
            });
        }
    };

    const bool fill_mom[AMREX_SPACEDIM] = {mom[0] != nullptr, mom[1] != nullptr, mom[2] != nullptr};

    // The ghost cells filled (see above): one in z for the solved components and the momenta;
    // for rho, also one in x (y) for the face average of the x- (y-) momentum, and two in z for
    // that of the z-momentum
    const IntVect ng_solved(0,0,1);
    const IntVect ng_rho(fill_mom[0] ? 1 : 0, fill_mom[1] ? 1 : 0, fill_mom[2] ? 2 : 1);
    const IntVect ngm(0,0,1);

    // The masks are kept with two ghost cells (one for faces) and used over the region filled
    const IntVect ngc(2);
    AMREX_ALWAYS_ASSERT(cons.nGrowVect().allGE(ngc));
    {
        Vector<int> comps = {Rho_comp};
        for (int c : cons_comps) {
            AMREX_ALWAYS_ASSERT(c > Rho_comp && c < cons.nComp());
            comps.push_back(c);
        }
        std::sort(comps.begin(), comps.end());
        comps.erase(std::unique(comps.begin(), comps.end()), comps.end());

        const iMultiFab& mask = get_cf_implicit_mask(lev, 0, cons.boxArray(), cons.DistributionMap(), ngc);

        MultiFab tmp(cons.boxArray(), cons.DistributionMap(), cons.nComp(), ngc);
        Vector<MultiFab*> fmf = {&tmp, &tmp};
        Vector<MultiFab*> cmf = select_time_data(csel, vars_old[lev-1][Vars::cons], vars_new[lev-1][Vars::cons]);

        // One fill per run of contiguous components, over the ghost cells any of them needs.
        // cell_cons_interp limits the slopes of each component separately, and the value it
        // gives a fine cell depends only on the coarse data around it, so this gives the same
        // result in each cell as one fill per component over all its ghost cells.
        const int ncomps = static_cast<int>(comps.size());
        for (int i0 = 0; i0 < ncomps; )
        {
            IntVect ngf = (comps[i0] == Rho_comp) ? ng_rho : ng_solved;
            int i1 = i0 + 1;
            while (i1 < ncomps && comps[i1] == comps[i1-1] + 1) {
                ngf.max((comps[i1] == Rho_comp) ? ng_rho : ng_solved);
                ++i1;
            }
            const int scomp = comps[i0];
            const int ncomp = i1 - i0;

            MultiFab::Copy(tmp, cons, scomp, scomp, ncomp, ngf);
            FillPatchTwoLevels(tmp, ngf, IntVect(0,0,0),
                               time_rel, cmf, ctime, fmf, ftime,
                               scomp, scomp, ncomp, geom[lev-1], geom[lev],
                               refRatio(lev-1), &cell_cons_interp, domain_bcs_type,
                               BCVars::cons_bc+scomp);
            copy_cf_ghosts(cons, tmp, mask, scomp, ncomp, ngf);

            i0 = i1;
        }
    }

    // Momenta: interpolate the velocity, then multiply by the face-averaged density
    const IntVect ngm_mask(1);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    {
        if (!fill_mom[dir]) { continue; }
        MultiFab& mom_dir = *mom[dir];
        AMREX_ALWAYS_ASSERT(mom_dir.nGrowVect().allGE(ngm_mask));

        MultiFab vel(mom_dir.boxArray(), mom_dir.DistributionMap(), 1, ngm);
        vel.setVal(zero);

        Vector<MultiFab*> fmf = {&vel, &vel};
        Vector<MultiFab*> cmf = select_time_data(csel, vars_old[lev-1][Vars::xvel+dir],
                                                       vars_new[lev-1][Vars::xvel+dir]);
        FillPatchTwoLevels(vel, ngm, IntVect(0,0,0),
                           time_rel, cmf, ctime, fmf, ftime,
                           0, 0, 1, geom[lev-1], geom[lev],
                           refRatio(lev-1), &face_cons_linear_interp, domain_bcs_type,
                           BCVars::xvel_bc+dir);

        const int ioff = (dir == 0) ? 1 : 0;
        const int joff = (dir == 1) ? 1 : 0;
        const int koff = (dir == 2) ? 1 : 0;

        MultiFab tmp(mom_dir.boxArray(), mom_dir.DistributionMap(), 1, ngm);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(tmp,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box gbx = mfi.growntilebox(ngm);
            const Array4<const Real> u   = vel.const_array(mfi);
            const Array4<const Real> rho = cons.const_array(mfi);
            const Array4<      Real> m   = tmp.array(mfi);
            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                m(i,j,k) = u(i,j,k) * myhalf * (rho(i,j,k,Rho_comp) + rho(i-ioff,j-joff,k-koff,Rho_comp));
            });
        }

        const iMultiFab& mask = get_cf_implicit_mask(lev, 1+dir, mom_dir.boxArray(),
                                                     mom_dir.DistributionMap(), ngm_mask);
        copy_cf_ghosts(mom_dir, tmp, mask, 0, 1, ngm);
    }
}
