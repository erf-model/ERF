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

    // Create aliases to state data to pass to calc_tpi_update
    int ncons = lev_new[Vars::cons].nComp();
    MultiFab cons_data(lev_new[Vars::cons], make_alias, 0, ncons);
    MultiFab xvel_data(lev_new[Vars::xvel], make_alias, 0, 1);
    MultiFab yvel_data(lev_new[Vars::yvel], make_alias, 0, 1);

    // This logic is done once then stored within TurbPertStruct.H
    turbPert.pt_type = -1;
    if (solverChoice.pert_type == PerturbationType::perturbSource) {
        turbPert.pt_type = 0;
    } else if (solverChoice.pert_type == PerturbationType::perturbDirect) {
        turbPert.pt_type = 1;
    }
    AMREX_ALWAYS_ASSERT(turbPert.pt_type >= 0);

    // Computing perturbation update time
    turbPert.calc_tpi_update(lev, local_dt, xvel_data, yvel_data, cons_data);

    Print() << "Successfully initialized turbulent perturbation update time and amplitude with type: "<< turbPert.pt_type <<"\n";
}

// Calculate the perturbation region amplitude.
// This function heavily emmulates the ERF::init_custom ()
void
ERF::turbPert_amplitude (int lev)
{
    // Accessing data
    auto& lev_new = vars_new[lev];

    // Creating local data
    int ncons = lev_new[Vars::cons].nComp();
    MultiFab cons_data(lev_new[Vars::cons], make_alias, 0, ncons);

    // Defining BoxArray type
    auto m_ixtype = cons_data.boxArray().ixType();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(lev_new[Vars::cons], TileNoZ()); mfi.isValid(); ++mfi) {
        const Box &bx  = mfi.validbox();
        const auto &cons_pert_arr = cons_data.array(mfi); // Address of perturbation array
        const amrex::Array4<const amrex::Real> &pert_cell = turbPert.pb_cell.array(mfi); // per-cell perturbation stored in structure

        turbPert.apply_tpi(lev, bx, RhoTheta_comp, m_ixtype, cons_pert_arr, pert_cell);
    } // mfi
}
