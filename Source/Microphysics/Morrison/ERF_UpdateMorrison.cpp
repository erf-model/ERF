#include "ERF_Morrison.H"
#include "ERF_Constants.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;
#include <stdbool.h>

/**
 * Updates conserved and microphysics variables in the provided MultiFabs from
 * the internal MultiFabs that store Microphysics module data.
 *
 * @param[out] cons Conserved variables
 * @param[out] qmoist: qv, qc, qi, qr, qs, qg
 */
void
Morrison::Copy_Micro_to_State (MultiFab& cons)
{
    // Get the temperature, density, theta, qt and qp from input
    for ( MFIter mfi(cons,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();

        auto states_arr = cons.array(mfi);

        auto rho_arr   = mic_fab_vars[MicVar_Morr::rho]->array(mfi);
        auto theta_arr = mic_fab_vars[MicVar_Morr::theta]->array(mfi);

        auto qv_arr    = mic_fab_vars[MicVar_Morr::qv]->array(mfi);
        auto qc_arr    = mic_fab_vars[MicVar_Morr::qcl]->array(mfi);
        auto qi_arr    = mic_fab_vars[MicVar_Morr::qci]->array(mfi);

        auto qpr_arr    = mic_fab_vars[MicVar_Morr::qpr]->array(mfi);
        auto qps_arr    = mic_fab_vars[MicVar_Morr::qps]->array(mfi);
        auto qpg_arr    = mic_fab_vars[MicVar_Morr::qpg]->array(mfi);

        auto nc_arr     = mic_fab_vars[MicVar_Morr::nc]->array(mfi);
        auto ni_arr     = mic_fab_vars[MicVar_Morr::ni]->array(mfi);
        auto nr_arr     = mic_fab_vars[MicVar_Morr::nr]->array(mfi);
        auto ns_arr     = mic_fab_vars[MicVar_Morr::ns]->array(mfi);
        auto ng_arr     = mic_fab_vars[MicVar_Morr::ng]->array(mfi);

        // get potential total density, temperature, qt, qp
        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            states_arr(i,j,k,RhoTheta_comp) = rho_arr(i,j,k)*theta_arr(i,j,k);

            states_arr(i,j,k,RhoQ1_comp)    = rho_arr(i,j,k)*std::max(Real(0),qv_arr(i,j,k));
            states_arr(i,j,k,RhoQ2_comp)    = rho_arr(i,j,k)*std::max(Real(0),qc_arr(i,j,k));
            states_arr(i,j,k,RhoQ3_comp)    = rho_arr(i,j,k)*std::max(Real(0),qi_arr(i,j,k));

            states_arr(i,j,k,RhoQ4_comp)    = rho_arr(i,j,k)*std::max(Real(0),qpr_arr(i,j,k));
            states_arr(i,j,k,RhoQ5_comp)    = rho_arr(i,j,k)*std::max(Real(0),qps_arr(i,j,k));
            states_arr(i,j,k,RhoQ6_comp)    = rho_arr(i,j,k)*std::max(Real(0),qpg_arr(i,j,k));

            states_arr(i,j,k,RhoQ7_comp)    = rho_arr(i,j,k)*std::max(Real(0),nc_arr(i,j,k));
            states_arr(i,j,k,RhoQ8_comp)    = rho_arr(i,j,k)*std::max(Real(0),ni_arr(i,j,k));
            states_arr(i,j,k,RhoQ9_comp)    = rho_arr(i,j,k)*std::max(Real(0),nr_arr(i,j,k));
            states_arr(i,j,k,RhoQ10_comp)   = rho_arr(i,j,k)*std::max(Real(0),ns_arr(i,j,k));
            states_arr(i,j,k,RhoQ11_comp)   = rho_arr(i,j,k)*std::max(Real(0),ng_arr(i,j,k));
        });
    }

    // Fill interior ghost cells and periodic boundaries
    cons.FillBoundary(m_geom.periodicity());
}
