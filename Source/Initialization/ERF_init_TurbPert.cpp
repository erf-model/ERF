/**
 * \ERF_init_TurbPert.cpp
 */

//DUSTIN MA: May 28th, 2024

#include <ERF.H>
#include <AMReX_MultiFabUtil.H>
#include <TileNoZ.H>
#include <prob_common.H>

using namespace amrex;

void
ERF::turbPert_update (const int lev, const Real local_dt)
{
    // Grabbing data from velocity field
    auto& lev_new = vars_new[lev];

    // Accessing local data
    MultiFab cons_data(lev_new[Vars::cons].boxArray(), lev_new[Vars::cons].DistributionMap(), 1, lev_new[Vars::cons].nGrowVect());
    MultiFab xvel_data(lev_new[Vars::xvel].boxArray(), lev_new[Vars::xvel].DistributionMap(), 1, lev_new[Vars::xvel].nGrowVect());
    MultiFab yvel_data(lev_new[Vars::yvel].boxArray(), lev_new[Vars::yvel].DistributionMap(), 1, lev_new[Vars::yvel].nGrowVect());
    MultiFab::Copy (cons_data, lev_new[Vars::cons], 0, 0, 1, lev_new[Vars::cons].nGrowVect());
    MultiFab::Copy (xvel_data, lev_new[Vars::xvel], 0, 0, 1, lev_new[Vars::xvel].nGrowVect());
    MultiFab::Copy (yvel_data, lev_new[Vars::yvel], 0, 0, 1, lev_new[Vars::yvel].nGrowVect());

    // Computing perturbation update time
    turbPert.calc_tpi_update(lev, local_dt, xvel_data, yvel_data, cons_data);

    Print() << "Successfully initialized turbulent perturbation update time and amplitude\n";
}

// Calculate the perturbation region amplitude.
// This function heavily emmulates the ERF::init_custom ()
void
ERF::turbPert_amplitude (int lev)
{
    auto& lev_new = vars_new[lev];

    MultiFab r_hse(base_state[lev], make_alias, 0, 1); // r_0 is first  component
    MultiFab p_hse(base_state[lev], make_alias, 1, 1); // p_0 is second component

    MultiFab cons_pert(lev_new[Vars::cons].boxArray(), lev_new[Vars::cons].DistributionMap(),
                       lev_new[Vars::cons].nComp()   , lev_new[Vars::cons].nGrow());
    MultiFab xvel_pert(lev_new[Vars::xvel].boxArray(), lev_new[Vars::xvel].DistributionMap(), 1, lev_new[Vars::xvel].nGrowVect());
    MultiFab yvel_pert(lev_new[Vars::yvel].boxArray(), lev_new[Vars::yvel].DistributionMap(), 1, lev_new[Vars::yvel].nGrowVect());
    MultiFab zvel_pert(lev_new[Vars::zvel].boxArray(), lev_new[Vars::zvel].DistributionMap(), 1, lev_new[Vars::zvel].nGrowVect());

    // Only storing for conserved quantity for now
    auto m_ixtype = cons_pert.boxArray().ixType();

    // Default perturbations to zero
    cons_pert.setVal(0.);
    xvel_pert.setVal(0.);
    yvel_pert.setVal(0.);
    zvel_pert.setVal(0.);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(lev_new[Vars::cons], TileNoZ()); mfi.isValid(); ++mfi) {
        const Box &bx  = mfi.validbox();
        const auto &cons_pert_arr = cons_pert.array(mfi);
        const amrex::Array4<const amrex::Real> &pert_cell = turbPert.pb_cell.array(mfi);

        turbPert.apply_tpi(lev, bx, RhoTheta_comp, m_ixtype, cons_pert_arr, pert_cell);

        /*
        prob->init_custom_pert(bx, xbx, ybx, zbx, cons_arr, cons_pert_arr,
                               xvel_pert_arr, yvel_pert_arr, zvel_pert_arr,
                               r_hse_arr, p_hse_arr, z_nd_arr, z_cc_arr,
                               geom[lev].data(), mf_m, mf_u, mf_v,
                               solverChoice);
       */
    } // mfi
}
