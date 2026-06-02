#include "ERF_SAM.H"
#include "ERF_IndexDefines.H"
#include "ERF_TileNoZ.H"

using namespace amrex;

/**
 * Updates conserved and microphysics variables in the provided MultiFabs from
 * the internal MultiFabs that store Microphysics module data.
 *
 * @param[out] cons Conserved variables
 * @param[out] qmoist: qv, qc, qi, qr, qs, qg
 */
void
SAM::Copy_Micro_to_State (MultiFab& cons)
{
    // Get the temperature, density, theta, qt and qp from input
    for ( MFIter mfi(cons,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();

        auto states_arr = cons.array(mfi);

        auto rho_arr    = mic_fab_vars[MicVar::rho]->array(mfi);
        auto theta_arr  = mic_fab_vars[MicVar::theta]->array(mfi);

        auto qv_arr     = mic_fab_vars[MicVar::qv]->array(mfi);
        auto qc_arr     = mic_fab_vars[MicVar::qcl]->array(mfi);
        auto qi_arr     = mic_fab_vars[MicVar::qci]->array(mfi);

        auto qpr_arr     = mic_fab_vars[MicVar::qpr]->array(mfi);
        auto qps_arr     = mic_fab_vars[MicVar::qps]->array(mfi);
        auto qpg_arr     = mic_fab_vars[MicVar::qpg]->array(mfi);

        // get potential total density, temperature, qt, qp
        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            SAMPrimitiveCell primitive{};
            primitive.rho = rho_arr(i,j,k);
            primitive.theta = theta_arr(i,j,k);
            primitive.qv = qv_arr(i,j,k);
            primitive.qcl = qc_arr(i,j,k);
            primitive.qci = qi_arr(i,j,k);
            primitive.qpr = qpr_arr(i,j,k);
            primitive.qps = qps_arr(i,j,k);
            primitive.qpg = qpg_arr(i,j,k);
            sam_primitive_to_cons(primitive, states_arr, i, j, k);
        });
    }

    // Fill interior ghost cells and periodic boundaries
    cons.FillBoundary(m_geom.periodicity());
}


