#include "ERF_SatAdj.H"

using namespace amrex;

/**
 * Compute Precipitation-related Microphysics quantities.
 */
void SatAdj::AdvanceSatAdj (const SolverChoice& /*solverChoice*/)
{
    if (!m_do_cond) { return; }

    auto tabs  = mic_fab_vars[MicVar_SatAdj::tabs];

    // Expose for GPU
    Real d_fac_cond = m_fac_cond;
    Real rdOcp      = m_rdOcp;

    // get the temperature, density, theta, qt and qc from input
    for ( MFIter mfi(*tabs,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        auto tbx = mfi.tilebox();

        auto qv_array    = mic_fab_vars[MicVar_SatAdj::qv]->array(mfi);
        auto qc_array    = mic_fab_vars[MicVar_SatAdj::qc]->array(mfi);
        auto tabs_array  = mic_fab_vars[MicVar_SatAdj::tabs]->array(mfi);
        auto theta_array = mic_fab_vars[MicVar_SatAdj::theta]->array(mfi);
        auto pres_array  = mic_fab_vars[MicVar_SatAdj::pres]->array(mfi);

        ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real T  = tabs_array(i,j,k);
            Real p  = pres_array(i,j,k);
            Real th = theta_array(i,j,k);
            Real qv = qv_array(i,j,k);
            Real qc = qc_array(i,j,k);

            AdjustSatAdjCell(d_fac_cond, rdOcp, T, p, th, qv, qc);

            tabs_array(i,j,k)  = T;
            theta_array(i,j,k) = th;
            qv_array(i,j,k)    = qv;
            qc_array(i,j,k)    = qc;
        });
    }
}
