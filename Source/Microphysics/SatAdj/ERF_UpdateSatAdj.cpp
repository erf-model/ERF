#include "ERF_SatAdj.H"

using namespace amrex;

/**
 * Updates conserved and microphysics variables in the provided MultiFabs from
 * the internal MultiFabs that store Microphysics module data.
 *
 * @param[out] cons Conserved variables
 * @param[out] qmoist: qv, qc
 */
void SatAdj::Copy_Micro_to_State (MultiFab& cons)
{
    // Get the temperature, density, theta, qt and qp from input
    for (MFIter mfi(cons,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto& tbx = mfi.tilebox();

        auto states_arr = cons.array(mfi);

        auto rho_arr    = mic_fab_vars[MicVar_SatAdj::rho]->array(mfi);
        auto theta_arr  = mic_fab_vars[MicVar_SatAdj::theta]->array(mfi);
        auto qv_arr     = mic_fab_vars[MicVar_SatAdj::qv]->array(mfi);
        auto qc_arr     = mic_fab_vars[MicVar_SatAdj::qc]->array(mfi);

        // get potential total density, temperature, qt, qp
        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            states_arr(i,j,k,RhoTheta_comp) = rho_arr(i,j,k)*theta_arr(i,j,k);
            states_arr(i,j,k,RhoQ1_comp)    = rho_arr(i,j,k)*qv_arr(i,j,k);
            states_arr(i,j,k,RhoQ2_comp)    = rho_arr(i,j,k)*qc_arr(i,j,k);
        });
    }

    // Fill interior ghost cells and periodic boundaries
    cons.FillBoundary(m_geom.periodicity());
}

