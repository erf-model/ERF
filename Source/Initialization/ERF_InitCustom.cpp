
/**
 * \file ERF_InitCustom.cpp
 */

#include <ERF.H>
#include <ERF_Constants.H>
#include <ERF_TileNoZ.H>
#include <ERF_ProbCommon.H>

using namespace amrex;

/**
 * Wrapper for custom problem-specific initialization routines that can be
 * defined by the user as they set up a new problem in ERF. This wrapper
 * handles all the overhead of defining the perturbation as well as initializing
 * the random seed if needed.
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

    MultiFab r_hse(base_state[lev], make_alias, BaseState::r0_comp, 1);
    MultiFab p_hse(base_state[lev], make_alias, BaseState::p0_comp, 1);

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

        Array4<Real const> cons_arr = lev_new[Vars::cons].const_array(mfi);
        Array4<Real const> z_nd_arr = (z_phys_nd[lev]) ? z_phys_nd[lev]->const_array(mfi) : Array4<Real const>{};
        Array4<Real const> z_cc_arr = (z_phys_cc[lev]) ? z_phys_cc[lev]->const_array(mfi) : Array4<Real const>{};

        // Here we arbitrarily choose the x-oriented map factor -- this should be generalized
        Array4<Real const> mf_m     = mapfac[lev][MapFacType::m_x]->const_array(mfi);
        Array4<Real const> mf_u     = mapfac[lev][MapFacType::u_x]->const_array(mfi);
        Array4<Real const> mf_v     = mapfac[lev][MapFacType::v_y]->const_array(mfi);

        Array4<Real> r_hse_arr = r_hse.array(mfi);
        Array4<Real> p_hse_arr = p_hse.array(mfi);

        prob->init_custom_pert(bx, xbx, ybx, zbx, cons_arr, cons_pert_arr,
                               xvel_pert_arr, yvel_pert_arr, zvel_pert_arr,
                               r_hse_arr, p_hse_arr, z_nd_arr, z_cc_arr,
                               geom[lev].data(), mf_m, mf_u, mf_v,
                               solverChoice, lev);
    } //mfi

    // Add problem-specific perturbation to background flow if not doing anelastic with fixed-in-time density
    if (!solverChoice.fixed_density[lev]) {
        MultiFab::Add(lev_new[Vars::cons], cons_pert, Rho_comp,      Rho_comp,             1, cons_pert.nGrow());
    }
    MultiFab::Add(lev_new[Vars::cons], cons_pert, RhoTheta_comp, RhoTheta_comp,        1, cons_pert.nGrow());
    MultiFab::Add(lev_new[Vars::cons], cons_pert, RhoScalar_comp,RhoScalar_comp,NSCALARS, cons_pert.nGrow());

    // RhoKE is relevant if using Deardorff with LES, k-equation for RANS, or MYNN with PBL
    if (solverChoice.turbChoice[lev].use_tke) {
        MultiFab::Add(lev_new[Vars::cons], cons_pert, RhoKE_comp,    RhoKE_comp,    1, cons_pert.nGrow());
    }

    if (solverChoice.moisture_type != MoistureType::None) {
        int qstate_size = micro->Get_Qstate_Size();
        for (int q_offset(0); q_offset<qstate_size; ++q_offset) {
            int q_idx = RhoQ1_comp+q_offset;
            MultiFab::Add(lev_new[Vars::cons], cons_pert, q_idx, q_idx, 1, cons_pert.nGrow());
        }
    }

    MultiFab::Add(lev_new[Vars::xvel], xvel_pert, 0,             0,             1, xvel_pert.nGrowVect());
    MultiFab::Add(lev_new[Vars::yvel], yvel_pert, 0,             0,             1, yvel_pert.nGrowVect());
    MultiFab::Add(lev_new[Vars::zvel], zvel_pert, 0,             0,             1, zvel_pert.nGrowVect());
}
