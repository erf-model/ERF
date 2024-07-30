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
    // Accessing data
    auto& lev_new = vars_new[lev];

    // Creating local data
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
    // Accessing data
    auto& lev_new = vars_new[lev];

    // Creating local data
    MultiFab cons_pert(lev_new[Vars::cons].boxArray(), lev_new[Vars::cons].DistributionMap(),
                       lev_new[Vars::cons].nComp()   , lev_new[Vars::cons].nGrow());

    // Initializing perturbation to zero
    cons_pert.setVal(0.);

    // Defining BoxArray type
    auto m_ixtype = cons_pert.boxArray().ixType();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(lev_new[Vars::cons], TileNoZ()); mfi.isValid(); ++mfi) {
        const Box &bx  = mfi.validbox();
        const auto &cons_pert_arr = cons_pert.array(mfi); // Address of perturbation array
        const amrex::Array4<const amrex::Real> &pert_cell = turbPert.pb_cell.array(mfi); // per-cell perturbation stored in structure

        turbPert.apply_tpi(lev, bx, RhoTheta_comp, m_ixtype, cons_pert_arr, pert_cell);
    } // mfi
}
