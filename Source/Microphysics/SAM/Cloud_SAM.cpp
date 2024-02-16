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
            int Tcase = 0;

            // Initial guess for temperature assuming no cloud water/ice:
            Real tabs = tabs_array(i,j,k);
            Real pres = pres_array(i,j,k);

            // Saturation moisture fractions
            Real omn, omp, omg;
            Real qsattw, qsatti, qsattc;
            Real dqsatw, dqsati, dqsatc;

            // Newton iteration vars
            int niter;
            Real fff, dfff, dtabs;
            Real lstarw, dlstarw;
            Real lstari, dlstari;
            Real lstarc, dlstarc;
            Real lstarp, dlstarp;
            Real delta_qn, delta_qv, delta_qc, delta_qi;

            //  Newton iteration for latent thermal sources
            niter = 0;
            dtabs = 1;
            do {
                // Latent heats and their derivatives wrt to T
                lstarw  = fac_cond;
                dlstarw = 0.0;
                lstari  = fac_fus;
                dlstari = 0.0;
                lstarc  = fac_fus;
                dlstarc = 0.0;

                // Cloud ice not permitted (condensation & melting)
                if(tabs >= tbgmax) {
                    erf_qsatw(  tabs, pres, qsattw);
                    erf_dtqsatw(tabs, pres, dqsatw);
                    qsattc = qcl_array(i,j,k);
                    dqsatc = 0.0;
                    qsatti = 0.0;
                    dqsati = 0.0;

                    Tcase = 0;
                }
                // Cloud water not permitted (sublimation & fusion)
                else if(tabs <= tbgmin) {
                    lstarw  = fac_sub;
                    dlstarw = 0.0;

                    erf_qsatw(  tabs, pres, qsattw);
                    erf_dtqsatw(tabs, pres, dqsatw);
                    qsattc = 0.0;
                    dqsatc = 0.0;
                    qsatti = qci_array(i,j,k);
                    dqsati = 0.0;

                    Tcase = 1;
                }
                // Mixed cloud phase (condensation & fusion)
                else {
                    erf_qsatw(  tabs, pres, qsattw);
                    erf_dtqsatw(tabs, pres, dqsatw);
                    qsattc = qcl_array(i,j,k);
                    dqsatc = 0.0;
                    erf_qsati(  tabs, pres, qsatti);
                    erf_dtqsati(tabs, pres, dqsati);

                    Tcase = 2;
                }

                // Function for root finding:
                // 0 = -T_new + T_old + L_v/C_p * (qv - qsatw) + L_f/C_p * (qi - qsati) - L_f/C_p * (qc - qsatc)
                fff   = -tabs + tabs_array(i,j,k)
                      +  lstarw*( qv_array(i,j,k)-qsattw)
                      +  lstarc*(qcl_array(i,j,k)-qsattc)
                      -  lstari*(qci_array(i,j,k)-qsatti);

                // Derivative of function (T_new iterated on)
                dfff  = -1.0
                      + ( dlstarw*( qv_array(i,j,k)-qsattw) - lstarw*dqsatw )
                      + ( dlstarc*(qcl_array(i,j,k)-qsattc) - lstarc*dqsatc )
                      - ( dlstari*(qci_array(i,j,k)-qsatti) - lstari*dqsati );

                // Update the temperature
                dtabs = -fff/dfff;
                tabs  = tabs+dtabs;
                niter = niter+1;

                /*
                if (i==0 && j==0) {
                    amrex::Print() << "Iter: " << Tcase << ' ' << niter << ' '
                                   << IntVect(i,j,k) << ' '
                                   << qv_array(i,j,k)  - qsattw << ' '
                                   << qcl_array(i,j,k) - qsattc << ' '
                                   << qci_array(i,j,k) - qsatti << ' '
                                   << qsattw << ' '
                                   << qsattc << ' '
                                   << qsatti << "\n";
                }
                */
            } while(std::abs(dtabs) > tol && niter < 20);

            // Update qsat from last iteration (dq = dq/dt * dt)
            qsattw   = qsattw + dqsatw*dtabs;
            qsattc   = qsattc + dqsatc*dtabs;
            qsatti   = qsatti + dqsati*dtabs;

            // Changes in each component
            delta_qv =  qv_array(i,j,k)-qsattw;
            delta_qc = qcl_array(i,j,k)-qsattc;
            delta_qi = qci_array(i,j,k)-qsatti;

            // The net change in q_i (d_qi or d_qc must be 0 for each case)
            delta_qn = (delta_qi + delta_qc) + delta_qv;

            // Partition the change in non-precipitating q
             qv_array(i,j,k) = qsattw;
            qcl_array(i,j,k) = (Tcase == 1) ? qsattc : (qcl_array(i,j,k)+delta_qn);
            qci_array(i,j,k) = (Tcase != 1) ? qsatti : (qci_array(i,j,k)+delta_qn);

             qv_array(i,j,k)  = std::max( qv_array(i,j,k),0.0);
            qcl_array(i,j,k) = std::max(qcl_array(i,j,k),0.0);
            qci_array(i,j,k) = std::max(qci_array(i,j,k),0.0);

             qn_array(i,j,k) = qcl_array(i,j,k) + qci_array(i,j,k);
             qt_array(i,j,k) =  qv_array(i,j,k) +  qn_array(i,j,k);

             /*
             if (i==0 && j==0) {
                 amrex::Print() << "Update: " << IntVect(i,j,k) << ' '
                                << getThgivenPandT(tabs, pres_array(i,j,k)*100, rdOcp) - theta_array(i,j,k) << ' '
                                << tabs - tabs_array(i,j,k) << "\n";
                 amrex::Print() << "dQ: " << Tcase << ' '
                                << delta_qv << ' '
                                << delta_qc << ' '
                                << delta_qi << ' '
                                << delta_qn << "\n";
                 amrex::Print() << "Q: " << Tcase << ' '
                                << qv_array(i,j,k) << ' '
                                << qcl_array(i,j,k) << ' '
                                << qci_array(i,j,k) << "\n";
                 amrex::Print() << "\n";
             }
             */

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
