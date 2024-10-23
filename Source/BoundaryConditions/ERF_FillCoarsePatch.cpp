#include <ERF.H>
#include <ERF_PhysBCFunct.H>
#include <ERF_IndexDefines.H>
#include <ERF_TimeInterpolatedData.H>
#include <ERF_FillPatcher.H>
#include <ERF_Utils.H>

using namespace amrex;

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

    //
    //****************************************************************************************************************
    // First fill velocities and density at the COARSE level so we can convert velocity to momenta at the COARSE level
    //****************************************************************************************************************
    //
    bool cons_only = false;
    FillPatch(lev-1, time, {&vars_new[lev-1][Vars::cons], &vars_new[lev-1][Vars::xvel],
                            &vars_new[lev-1][Vars::yvel], &vars_new[lev-1][Vars::zvel]},
                           {&vars_new[lev-1][Vars::cons],
                            &rU_new[lev-1], &rV_new[lev-1], &rW_new[lev-1]},
                            false, cons_only);

    //
    // ************************************************
    // Convert velocity to momentum at the COARSE level
    // ************************************************
    //
    VelocityToMomentum(vars_new[lev-1][Vars::xvel], IntVect{0},
                       vars_new[lev-1][Vars::yvel], IntVect{0},
                       vars_new[lev-1][Vars::zvel], IntVect{0},
                       vars_new[lev-1][Vars::cons],
                         rU_new[lev-1],
                         rV_new[lev-1],
                         rW_new[lev-1],
                       Geom(lev).Domain(),
                       domain_bcs_type);
    //
    // *****************************************************************
    // Interpolate all cell-centered variables from coarse to fine level
    // *****************************************************************
    //
    Interpolater* mapper_c = &cell_cons_interp;
    Interpolater* mapper_f = &face_cons_linear_interp;

    //
    //************************************************************************************************
    // Interpolate cell-centered data from coarse to fine level
    // with InterpFromCoarseLevel which ASSUMES that all ghost cells have already been filled
    // ************************************************************************************************
    IntVect ngvect_cons = vars_new[lev][Vars::cons].nGrowVect();
    int      ncomp_cons = vars_new[lev][Vars::cons].nComp();

    InterpFromCoarseLevel(vars_new[lev  ][Vars::cons], ngvect_cons, IntVect(0,0,0),
                          vars_new[lev-1][Vars::cons], 0, 0, ncomp_cons,
                          geom[lev-1], geom[lev],
                          refRatio(lev-1), mapper_c, domain_bcs_type, BCVars::cons_bc);

    //
    //************************************************************************************************
    // Interpolate x-momentum from coarse to fine level
    // with InterpFromCoarseLevel which ASSUMES that all ghost cells have already been filled
    // ************************************************************************************************
    //
    InterpFromCoarseLevel(rU_new[lev], IntVect{0}, IntVect{0}, rU_new[lev-1], 0, 0, 1,
                          geom[lev-1], geom[lev],
                          refRatio(lev-1), mapper_f, domain_bcs_type, BCVars::xvel_bc);

    //
    //************************************************************************************************
    // Interpolate y-momentum from coarse to fine level
    // with InterpFromCoarseLevel which ASSUMES that all ghost cells have already been filled
    // ************************************************************************************************
    //
    InterpFromCoarseLevel(rV_new[lev], IntVect{0}, IntVect{0}, rV_new[lev-1], 0, 0, 1,
                          geom[lev-1], geom[lev],
                          refRatio(lev-1), mapper_f, domain_bcs_type, BCVars::yvel_bc);

    //************************************************************************************************
    // Interpolate z-momentum from coarse to fine level
    // with InterpFromCoarseLevel which ASSUMES that all ghost cells have already been filled
    // ************************************************************************************************
    InterpFromCoarseLevel(rW_new[lev],  IntVect{0}, IntVect{0}, rW_new[lev-1], 0, 0, 1,
                          geom[lev-1], geom[lev],
                          refRatio(lev-1), mapper_f, domain_bcs_type, BCVars::zvel_bc);
    //
    // *********************************************************
    // After interpolation of momentum, convert back to velocity
    // *********************************************************
    //
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

    // ***************************************************************************
    // Physical bc's at domain boundary
    // ***************************************************************************
    IntVect ngvect_vels = vars_new[lev][Vars::xvel].nGrowVect();

    (*physbcs_cons[lev])(vars_new[lev][Vars::cons],0,ncomp_cons,ngvect_cons,time,BCVars::cons_bc,true);
    (   *physbcs_u[lev])(vars_new[lev][Vars::xvel],0,1         ,ngvect_vels,time,BCVars::xvel_bc,true);
    (   *physbcs_v[lev])(vars_new[lev][Vars::yvel],0,1         ,ngvect_vels,time,BCVars::yvel_bc,true);
    (   *physbcs_w[lev])(vars_new[lev][Vars::zvel],vars_new[lev][Vars::xvel],vars_new[lev][Vars::yvel],
                         ngvect_vels,time,BCVars::zvel_bc,true);

    // ***************************************************************************
    // Since lev > 0 here we don't worry about m_r2d or wrfbdy data
    // ***************************************************************************
}
