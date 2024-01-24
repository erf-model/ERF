#include "MicrophysicsVariables.H"
#include "Kessler.H"
#include "IndexDefines.H"
#include "TileNoZ.H"

/**
 * Updates conserved and microphysics variables in the provided MultiFabs from
 * the internal MultiFabs that store Microphysics module data.
 *
 * @param[out] cons Conserved variables
 * @param[out] qmoist: qv, qc, qi, qr, qs, qg
 */
void Kessler::Copy_Micro_to_State (amrex::MultiFab& cons)
{
    // Get the temperature, density, theta, qt and qp from input
    for ( amrex::MFIter mfi(cons,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();

        auto states_arr = cons.array(mfi);

        auto rho_arr    = mic_fab_vars[MicVar_Kess::rho]->array(mfi);
        auto theta_arr  = mic_fab_vars[MicVar_Kess::theta]->array(mfi);
        auto qv_arr     = mic_fab_vars[MicVar_Kess::qv]->array(mfi);
        auto qc_arr     = mic_fab_vars[MicVar_Kess::qcl]->array(mfi);
        auto qp_arr     = mic_fab_vars[MicVar_Kess::qp]->array(mfi);

        // get potential total density, temperature, qt, qp
        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            states_arr(i,j,k,RhoTheta_comp) = rho_arr(i,j,k)*theta_arr(i,j,k);
            states_arr(i,j,k,RhoQ1_comp)    = rho_arr(i,j,k)*qv_arr(i,j,k);
            states_arr(i,j,k,RhoQ2_comp)    = rho_arr(i,j,k)*qc_arr(i,j,k);
            states_arr(i,j,k,RhoQ3_comp)    = rho_arr(i,j,k)*qp_arr(i,j,k);
        });
    }

    // Fill interior ghost cells and periodic boundaries
    cons.FillBoundary(m_geom.periodicity());
}

