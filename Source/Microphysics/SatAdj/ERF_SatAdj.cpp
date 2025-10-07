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

        const auto& tbx = mfi.tilebox();

        auto qv_array    = mic_fab_vars[MicVar_SatAdj::qv]->array(mfi);
        auto qc_array    = mic_fab_vars[MicVar_SatAdj::qc]->array(mfi);
        auto tabs_array  = mic_fab_vars[MicVar_SatAdj::tabs]->array(mfi);
        auto theta_array = mic_fab_vars[MicVar_SatAdj::theta]->array(mfi);
        auto pres_array  = mic_fab_vars[MicVar_SatAdj::pres]->array(mfi);

        ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            //qc_array(i,j,k) = std::max(0.0, qc_array(i,j,k));

            //------- Evaporation/condensation
            Real qsat;
            erf_qsatw(tabs_array(i,j,k), pres_array(i,j,k), qsat);

            // There is enough moisture to drive to equilibrium
            if ((qv_array(i,j,k)+qc_array(i,j,k)) > qsat) {
                Real qvprev = qv_array(i,j,k);
                Real qcprev = qc_array(i,j,k);

                // clip qc but maintain total water
                if (qc_array(i,j,k) < 0) {
                    qv_array(i,j,k) += qc_array(i,j,k);
                    qc_array(i,j,k)  = 0.0;
                }

                // Update temperature
                tabs_array(i,j,k) = NewtonIterSat(i, j, k   ,
                                                  d_fac_cond, tabs_array, pres_array,
                                                  qv_array  , qc_array  );

                Real qsatnew;
                erf_qsatw(tabs_array(i,j,k), pres_array(i,j,k), qsatnew);

                AMREX_ASSERT(std::abs(qv_array(i,j,k)-qsatnew) < 1e-12);

                amrex::ignore_unused(qvprev);
                amrex::ignore_unused(qcprev);
                AMREX_ASSERT(std::abs(qv_array(i,j,k)+qc_array(i,j,k)-qvprev-qcprev) < 1e-14);

                // Update theta (constant pressure)
                theta_array(i,j,k) = getThgivenTandP(tabs_array(i,j,k), 100.0*pres_array(i,j,k), rdOcp);

            //
            // We cannot blindly relax to qsat, but we can convert qc/qi -> qv.
            // The concept here is that if we put all the moisture into qv and modify
            // the temperature, we can then check if qv > qsat occurs (for final T/P/qv).
            // If the reduction in T/qsat and increase in qv does trigger the
            // aforementioned condition, we can do Newton iteration to drive qv = qsat.
            //
            } else {
                // Changes in each component
                Real delta_qc = qc_array(i,j,k);

                // Partition the change in non-precipitating q
                qv_array(i,j,k) += qc_array(i,j,k);
                qc_array(i,j,k)  = 0.0;

                // Update temperature (endothermic since we evap/sublime)
                tabs_array(i,j,k) -= d_fac_cond * delta_qc;

                // Update theta
                theta_array(i,j,k) = getThgivenTandP(tabs_array(i,j,k), 100.0*pres_array(i,j,k), rdOcp);

                // Verify assumption that qv > qsat does not occur
                erf_qsatw(tabs_array(i,j,k), pres_array(i,j,k), qsat);
                if (qv_array(i,j,k) > qsat) {
                    Real qvprev = qv_array(i,j,k);
                    Real qcprev = qc_array(i,j,k);
                    Real Tprev = tabs_array(i,j,k);

                    // Update temperature
                    tabs_array(i,j,k) = NewtonIterSat(i, j, k     ,
                                                      d_fac_cond  , tabs_array, pres_array,
                                                      qv_array    , qc_array  );

                    Real qsatnew;
                    erf_qsatw(tabs_array(i,j,k), pres_array(i,j,k), qsatnew);
                    amrex::ignore_unused(qvprev);
                    amrex::ignore_unused(qcprev);
                    amrex::ignore_unused(Tprev);
                    AMREX_ASSERT(qv_array(i,j,k) < qvprev);
                    AMREX_ASSERT(qc_array(i,j,k) > qcprev);
                    AMREX_ASSERT(tabs_array(i,j,k) > Tprev);
                    AMREX_ASSERT(std::abs(qv_array(i,j,k)-qsatnew) < 1e-14);
                    AMREX_ASSERT(std::abs(qv_array(i,j,k)+qc_array(i,j,k)-qvprev-qcprev) < 1e-14);

                    // Update theta
                    theta_array(i,j,k) = getThgivenTandP(tabs_array(i,j,k), 100.0*pres_array(i,j,k), rdOcp);

                }
            }
        });
    }
}
