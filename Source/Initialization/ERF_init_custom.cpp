/**
 * \file ERF_init_custom.cpp
 */

#include <ERF.H>
#include <EOS.H>
#include <ERF_Constants.H>
#include <Utils.H>
#include <prob_common.H>

using namespace amrex;

void
ERF::init_custom(int lev)
{
    auto& lev_new = vars_new[lev];
#if defined(ERF_USE_MOISTURE)
    auto& qv_new  = qv[lev];
    auto& qc_new  = qc[lev];
    auto& qi_new  = qi[lev];
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
    MultiFab qv_pert(qv[lev].boxArray(), qv[lev].DistributionMap(), 1, qv[lev].nGrow());
    MultiFab qc_pert(qc[lev].boxArray(), qc[lev].DistributionMap(), 1, qc[lev].nGrow());
    MultiFab qi_pert(qi[lev].boxArray(), qi[lev].DistributionMap(), 1, qi[lev].nGrow());
    qv_pert.setVal(0.);
    qc_pert.setVal(0.);
    qi_pert.setVal(0.);
#elif defined(ERF_USE_WARM_NO_PRECIP)
    MultiFab qv_pert(cons_pert.boxArray(), cons_pert.DistributionMap(), 1, cons_pert.nGrow());
    MultiFab qc_pert(cons_pert.boxArray(), cons_pert.DistributionMap(), 1, cons_pert.nGrow());
    qv_pert.setVal(0.);
    qc_pert.setVal(0.);
#endif

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.tilebox();

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
        init_custom_prob(bx, cons_pert_arr, xvel_pert_arr, yvel_pert_arr, zvel_pert_arr,
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
    MultiFab::Add(             qv_new,   qv_pert, 0,             0,             1,   qv_pert.nGrow());
    MultiFab::Add(             qc_new,   qc_pert, 0,             0,             1,   qc_pert.nGrow());
    MultiFab::Add(             qi_new,   qi_pert, 0,             0,             1,   qi_pert.nGrow());
#elif defined(ERF_USE_WARM_NO_PRECIP)
    MultiFab::Add(lev_new[Vars::cons], cons_pert, RhoQv_comp,    RhoQv_comp,    1, cons_pert.nGrow());
    MultiFab::Add(lev_new[Vars::cons], cons_pert, RhoQc_comp,    RhoQc_comp,    1, cons_pert.nGrow());
#endif
    MultiFab::Add(lev_new[Vars::xvel], xvel_pert, 0,             0,             1, xvel_pert.nGrowVect());
    MultiFab::Add(lev_new[Vars::yvel], yvel_pert, 0,             0,             1, yvel_pert.nGrowVect());
    MultiFab::Add(lev_new[Vars::zvel], zvel_pert, 0,             0,             1, zvel_pert.nGrowVect());
}
