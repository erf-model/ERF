#include "SAM.H"
#include "IndexDefines.H"
#include "TileNoZ.H"
#include "EOS.H"

using namespace amrex;

/**
 * Compute Cloud-related Microphysics quantities.
 */
void SAM::Cloud () {

    constexpr Real an = 1.0/(tbgmax-tbgmin);
    constexpr Real bn = tbgmin*an;
    constexpr Real ap = 1.0/(tprmax-tprmin);
    constexpr Real bp = tprmin*ap;
    constexpr Real ag = 1.0/(tgrmax-tgrmin);
    constexpr Real bg = tgrmin*ag;

    Real fac_cond = m_fac_cond;
    Real fac_sub  = m_fac_sub;
    Real fac_fus  = m_fac_fus;
    Real rdOcp    = m_rdOcp;

    Real tol = 1.0e-4;

    for ( MFIter mfi(*(mic_fab_vars[MicVar::tabs]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto  qt_array = mic_fab_vars[MicVar::qt]->array(mfi);
        auto  qn_array = mic_fab_vars[MicVar::qn]->array(mfi);
        auto  qv_array = mic_fab_vars[MicVar::qv]->array(mfi);
        auto qcl_array = mic_fab_vars[MicVar::qcl]->array(mfi);
        auto qci_array = mic_fab_vars[MicVar::qci]->array(mfi);

        auto  qp_array = mic_fab_vars[MicVar::qp]->array(mfi);
        auto qpr_array = mic_fab_vars[MicVar::qpr]->array(mfi);
        auto qps_array = mic_fab_vars[MicVar::qps]->array(mfi);
        auto qpg_array = mic_fab_vars[MicVar::qpg]->array(mfi);

        auto   rho_array = mic_fab_vars[MicVar::rho]->array(mfi);
        auto  tabs_array = mic_fab_vars[MicVar::tabs]->array(mfi);
        auto theta_array = mic_fab_vars[MicVar::theta]->array(mfi);
        auto  pres_array = mic_fab_vars[MicVar::pres]->array(mfi);

        const auto& box3d = mfi.tilebox();

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            int dbg = 0;

            // Initial guess for temperature assuming no cloud water/ice:
            Real tabs = tabs_array(i,j,k);
            Real pres = pres_array(i,j,k);

            // Saturation moisture fractions
            Real omn, omp, omg;
            Real qsattw, qsatti;
            Real dqsatw, dqsati;

            // Newton iteration vars
            int niter;
            Real fff, dfff, dtabs;
            Real lstarw, dlstarw;
            Real lstari, dlstari;
            Real lstarp, dlstarp;
            Real delta_qv, delta_qc, delta_qi;

            // HACK (remove precip from Newton iters)
            tabs += fac_cond*qp_array(i,j,k);

            /*
            // Remove precipitation term from Newton iteration if
            // the temperature is out of range.
            if(tabs > tbgmax) {
                tabs += fac_cond*qp_array(i,j,k);
            }
            // Ice cloud:
            else if(tabs <= tbgmin) {
                tabs += fac_sub*qp_array(i,j,k);
            }
            */

            //  Newton iteration for latent thermal sources
            niter = 0;
            dtabs = 1;
            do {
                // Hack only do non-precip latent heat
                lstarw  = fac_cond;
                dlstarw = 0.0;
                lstari  = 0.0; //fac_fus;
                dlstari = 0.0;

                erf_qsatw(tabs, pres, qsattw);
                erf_qsati(tabs, pres, qsatti);
                erf_dtqsatw(tabs, pres, dqsatw);
                erf_dtqsati(tabs, pres, dqsati);

                qsatti = qci_array(i,j,k);
                dqsati = 0.0;

                /*
                // Non-precipitating components
                if(tabs >= tbgmax) {       // Warm cloud (condensation)
                    omn     = 1.0;

                    lstarw  = fac_cond;
                    dlstarw = 0.0;
                    lstari  = fac_fus;
                    dlstari = 0.0;

                    erf_qsatw(tabs, pres, qsattw);
                    erf_dtqsatw(tabs, pres, dqsatw);
                    qsatti = 0.0;
                    dqsati = 0.0;

                    dbg = 0;
                }
                else if(tabs <= tbgmin) {  // Ice cloud (sublimation)
                    omn     = 0.0;

                    lstarw  = fac_cond;
                    dlstarw = 0.0;
                    lstari  = fac_fus; //fac_sub;
                    dlstari = 0.0;

                    qsattw = 0.0;
                    dqsatw = 0.0;
                    erf_qsati(tabs, pres, qsatti);
                    erf_dtqsati(tabs, pres, dqsati);

                    dbg = 1;
                }
                else {                      // Mixed-phase cloud (condensation + fusion)
                    omn     = an*tabs-bn;

                    lstarw  = fac_cond;
                    dlstarw = 0.0;
                    lstari  = fac_fus;
                    dlstari = 0.0;

                    erf_qsatw(tabs, pres, qsattw);
                    erf_qsati(tabs, pres, qsatti);
                    erf_dtqsatw(tabs, pres, dqsatw);
                    erf_dtqsati(tabs, pres, dqsati);

                    dbg = 2;
                }


                // Precipitating components
                if(tabs >= tprmax) {
                    omp     = 1.0;
                    lstarp  = fac_cond;
                    dlstarp = 0.0;
                }
                else if(tabs <= tprmin) {
                    omp     = 0.0;
                    lstarp  = fac_sub;
                    dlstarp = 0.0;
                }
                else {
                    omp     = ap*tabs-bp;
                    lstarp  = fac_cond+(1.0-omp)*fac_fus;
                    dlstarp = -ap*fac_fus;
                }
                */

                // Function for root finding:
                // 0 = -T_new + T_old + L_v/C_p * (qv - qsatw) + L_i/C_p * (qi - qsati) - L_p/C_p * (qp)
                fff   = -tabs + tabs_array(i,j,k)
                      +  lstarw*( qv_array(i,j,k)-qsattw)
                      +  lstari*(qci_array(i,j,k)-qsatti)
                      -  lstarp*qp_array(i,j,k);

                // Derivative of function (T_new iterated on)
                dfff  = -1.0
                      + ( dlstarw*( qv_array(i,j,k)-qsattw) - lstarw*dqsatw )
                      + ( dlstari*(qci_array(i,j,k)-qsatti) - lstari*dqsati )
                      - dlstarp*qp_array(i,j,k);

                // Update the temperature
                dtabs = -fff/dfff;
                tabs  = tabs+dtabs;
                niter = niter+1;
            } while(std::abs(dtabs) > tol && niter < 20);

            // Update qsat from last iteration (dq = dq/dt * dt)
            qsattw   = qsattw + dqsatw*dtabs;
            qsatti   = qsatti + dqsati*dtabs;
            delta_qv = qv_array(i,j,k)-qsattw;
            delta_qi = qci_array(i,j,k)-qsatti;
            delta_qc = delta_qv - delta_qi;

            // Partition the change in non-precipitating q
             qv_array(i,j,k)  = qsattw;
            qci_array(i,j,k)  = qsatti;
            qcl_array(i,j,k)  = std::max(0.0,qcl_array(i,j,k) + delta_qc);
             qn_array(i,j,k)  = qcl_array(i,j,k) + qci_array(i,j,k);
             qt_array(i,j,k)  =  qv_array(i,j,k) +  qn_array(i,j,k);

            // Partition the change in precipitating q
            qpr_array(i,j,k) = 0.0;
            qps_array(i,j,k) = 0.0;
            qpg_array(i,j,k) = 0.0;
             qp_array(i,j,k) = qpr_array(i,j,k) + qps_array(i,j,k) + qpg_array(i,j,k);

            // Update temperature
            tabs_array(i,j,k) = tabs;

            // Update pressure
            pres_array(i,j,k) = rho_array(i,j,k) * R_d * tabs_array(i,j,k)
                              * (1.0 + R_v/R_d * qv_array(i,j,k));

            // Update theta from temperature
            theta_array(i,j,k) = getThgivenPandT(tabs_array(i,j,k), pres_array(i,j,k), rdOcp);

            // Pressure unit conversion
            pres_array(i,j,k) /= 100.0;
        });
    }
}
