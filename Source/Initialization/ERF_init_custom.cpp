/**
 * \file ERF_init_custom.cpp
 */

#include <ERF.H>
#include <ERF_Constants.H>
#include <TileNoZ.H>
#include <prob_common.H>

using namespace amrex;

/**
 * Wrapper for custom problem-specific initialization routines that can be
 * defined by the user as they set up a new problem in ERF. This wrapper
 * handles all the overhead of defining both the background and perturbation
 * state as well as initializing the random seed.
 *
 * This wrapper calls a user function to customize initialization on a per-Fab
 * level inside an MFIter loop, so all the MultiFab operations are hidden from
 * the user.
 *
 * @param lev Integer specifying the current level
 */
void
ERF::init_custom (int lev)
{
    auto& lev_new = vars_new[lev];
#if defined(ERF_USE_MOISTURE)
    auto& qmoist_new  = qmoist[lev];
#endif
    MultiFab r_hse(base_state[lev], make_alias, 0, 1); // r_0 is first  component
    MultiFab p_hse(base_state[lev], make_alias, 1, 1); // p_0 is second component

    MultiFab cons_pert(lev_new[Vars::cons].boxArray(), lev_new[Vars::cons].DistributionMap(),
                       lev_new[Vars::cons].nComp()   , lev_new[Vars::cons].nGrow());
    MultiFab xvel_pert(lev_new[Vars::xvel].boxArray(), lev_new[Vars::xvel].DistributionMap(), 1, lev_new[Vars::xvel].nGrowVect());
    MultiFab yvel_pert(lev_new[Vars::yvel].boxArray(), lev_new[Vars::yvel].DistributionMap(), 1, lev_new[Vars::yvel].nGrowVect());
    MultiFab zvel_pert(lev_new[Vars::zvel].boxArray(), lev_new[Vars::zvel].DistributionMap(), 1, lev_new[Vars::zvel].nGrowVect());

    // Default all perturbations to zero
    cons_pert.setVal(0.);
    xvel_pert.setVal(0.);
    yvel_pert.setVal(0.);
    zvel_pert.setVal(0.);

#if defined(ERF_USE_MOISTURE)
    MultiFab qmoist_pert(qmoist[lev].boxArray(), qmoist[lev].DistributionMap(), 3, qmoist[lev].nGrow());
    qmoist_pert.setVal(0.);
    MultiFab qv_pert(qmoist_pert, amrex::make_alias, 0, 1);
    MultiFab qc_pert(qmoist_pert, amrex::make_alias, 1, 1);
    MultiFab qi_pert(qmoist_pert, amrex::make_alias, 2, 1);
#elif defined(ERF_USE_WARM_NO_PRECIP)
    MultiFab qmoist_pert(cons_pert.boxArray(), cons_pert.DistributionMap(), 2, cons_pert.nGrow());
    qmoist_pert.setVal(0.);
    MultiFab qv_pert(qmoist_pert, amrex::make_alias, 0, 1);
    MultiFab qc_pert(qmoist_pert, amrex::make_alias, 1, 1);
#endif

    int fix_random_seed = 0;
    ParmParse pp("erf"); pp.query("fix_random_seed", fix_random_seed);
    // Note that the value of 1024UL is not significant -- the point here is just to set the
    //     same seed for all MPI processes for the purpose of regression testing
    if (fix_random_seed) {
        amrex::Print() << "Fixing the random seed" << std::endl;
        amrex::InitRandom(1024UL);
    }

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(lev_new[Vars::cons], TileNoZ()); mfi.isValid(); ++mfi)
    {
        const Box &bx  = mfi.tilebox();
        const Box &xbx = mfi.tilebox(IntVect(1,0,0));
        const Box &ybx = mfi.tilebox(IntVect(0,1,0));
        const Box &zbx = mfi.tilebox(IntVect(0,0,1));

        const auto &cons_pert_arr = cons_pert.array(mfi);
        const auto &xvel_pert_arr = xvel_pert.array(mfi);
        const auto &yvel_pert_arr = yvel_pert.array(mfi);
        const auto &zvel_pert_arr = zvel_pert.array(mfi);

        Array4<Real const> z_nd_arr = (solverChoice.use_terrain) ? z_phys_nd[lev]->const_array(mfi) : Array4<Real const>{};
        Array4<Real const> z_cc_arr = (solverChoice.use_terrain) ? z_phys_cc[lev]->const_array(mfi) : Array4<Real const>{};

        Array4<Real const> mf_m     = mapfac_m[lev]->array(mfi);
        Array4<Real const> mf_u     = mapfac_m[lev]->array(mfi);
        Array4<Real const> mf_v     = mapfac_m[lev]->array(mfi);

        Array4<Real> r_hse_arr = r_hse.array(mfi);
        Array4<Real> p_hse_arr = p_hse.array(mfi);

#if defined(ERF_USE_MOISTURE)
        const auto &qv_pert_arr = qv_pert.array(mfi);
        const auto &qc_pert_arr = qc_pert.array(mfi);
        const auto &qi_pert_arr = qi_pert.array(mfi);
#elif defined(ERF_USE_WARM_NO_PRECIP)
        const auto &qv_pert_arr = qv_pert.array(mfi);
        const auto &qc_pert_arr = qc_pert.array(mfi);
#endif
        prob->init_custom_prob(bx, xbx, ybx, zbx, cons_pert_arr, xvel_pert_arr, yvel_pert_arr, zvel_pert_arr,
                         r_hse_arr, p_hse_arr, z_nd_arr, z_cc_arr,
#if defined(ERF_USE_MOISTURE)
                         qv_pert_arr, qc_pert_arr, qi_pert_arr,
#elif defined(ERF_USE_WARM_NO_PRECIP)
                         qv_pert_arr, qc_pert_arr,
#endif
                         geom[lev].data(), mf_m, mf_u, mf_v,
                         solverChoice);
    } //mfi

    // Add problem-specific perturbation to background flow
    MultiFab::Add(lev_new[Vars::cons], cons_pert, Rho_comp,      Rho_comp,      1, cons_pert.nGrow());
    MultiFab::Add(lev_new[Vars::cons], cons_pert, RhoTheta_comp, RhoTheta_comp, 1, cons_pert.nGrow());
    MultiFab::Add(lev_new[Vars::cons], cons_pert, RhoScalar_comp,RhoScalar_comp,1, cons_pert.nGrow());
    MultiFab::Add(lev_new[Vars::cons], cons_pert, RhoQKE_comp,   RhoQKE_comp,   1, cons_pert.nGrow());
#if defined(ERF_USE_MOISTURE)
    MultiFab::Add(lev_new[Vars::cons], cons_pert, RhoQt_comp,    RhoQt_comp,    1, cons_pert.nGrow());
    MultiFab::Add(lev_new[Vars::cons], cons_pert, RhoQp_comp,    RhoQp_comp,    1, cons_pert.nGrow());
    MultiFab::Add(         qmoist_new, qmoist_pert, 0,           0,             3, qmoist_pert.nGrow()); // qv, qc, qi
#elif defined(ERF_USE_WARM_NO_PRECIP)
    MultiFab::Add(lev_new[Vars::cons], cons_pert, RhoQv_comp,    RhoQv_comp,    1, cons_pert.nGrow());
    MultiFab::Add(lev_new[Vars::cons], cons_pert, RhoQc_comp,    RhoQc_comp,    1, cons_pert.nGrow());
#endif
    MultiFab::Add(lev_new[Vars::xvel], xvel_pert, 0,             0,             1, xvel_pert.nGrowVect());
    MultiFab::Add(lev_new[Vars::yvel], yvel_pert, 0,             0,             1, yvel_pert.nGrowVect());
    MultiFab::Add(lev_new[Vars::zvel], zvel_pert, 0,             0,             1, zvel_pert.nGrowVect());
}
