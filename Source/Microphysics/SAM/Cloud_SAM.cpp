#include "Microphysics.H"
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

    auto pres1d_t = pres1d.table();

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

        const auto& box3d = mfi.tilebox() & m_gtoe[mfi.index()];

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            int dbg = 0;

            // Initial guess for temperature assuming no cloud water/ice:
            Real tabs1 = tabs_array(i,j,k);
            Real pres1 = pres_array(i,j,k);

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

            //  Newton iteration for latent thermal sources
            Real delta_qv, delta_qc, delta_qi;
            //if(qt_array(i,j,k) > qsatt) {
                niter = 0;
                dtabs = 1;
                do {
                    // Non-precipitating components
                    if(tabs1 >= tbgmax) {       // Warm cloud (condensation)
                        omn     = 1.0;

                        lstarw  = fac_cond;
                        dlstarw = 0.0;
                        lstari  = fac_fus;
                        dlstari = 0.0;

                        erf_qsatw(tabs1, pres1, qsattw);
                        erf_dtqsatw(tabs1, pres1, dqsatw);
                        qsatti = 0.0;
                        dqsati = 0.0;

                        dbg = 0;
                    }
                    else if(tabs1 <= tbgmin) {  // Ice cloud (sublimation)
                        omn     = 0.0;

                        lstarw  = fac_cond;
                        dlstarw = 0.0;
                        lstari  = fac_sub;
                        dlstari = 0.0;

                        qsattw = 0.0;
                        dqsatw = 0.0;
                        erf_qsati(tabs1, pres1, qsatti);
                        erf_dtqsati(tabs1, pres1, dqsati);

                        dbg = 1;
                    }
                    else {                      // Mixed-phase cloud (condensation + fusion)
                        omn     = an*tabs1-bn;

                        lstarw  = fac_cond;
                        dlstarw = 0.0;
                        lstari  = fac_fus;
                        dlstari = 0.0;

                        erf_qsatw(tabs1, pres1, qsattw);
                        erf_qsati(tabs1, pres1, qsatti);
                        erf_dtqsatw(tabs1, pres1, dqsatw);
                        erf_dtqsati(tabs1, pres1, dqsati);

                        dbg = 2;
                    }

                    // Precipitating components
                    if(tabs1 >= tprmax) {
                        omp     = 1.0;
                        lstarp  = fac_cond;
                        dlstarp = 0.0;
                    }
                    else if(tabs1 <= tprmin) {
                        omp     = 0.0;
                        lstarp  = fac_sub;
                        dlstarp = 0.0;
                    }
                    else {
                        omp     = ap*tabs1-bp;
                        lstarp  = fac_cond+(1.0-omp)*fac_fus;
                        dlstarp = -ap*fac_fus;
                    }

                    if (k>50 && k<80 && i==0 && j==0) {
                        amrex::Print() << "Iter: " << niter << ' '
                                       << dbg << ' '
                                       << qv_array(i,j,k) << ' '
                                       << qsattw << ' '
                                       << qci_array(i,j,k) << ' '
                                       << qsatti << ' '
                                       << dtabs << "\n";
                    }

                    // Graupel fraction
                    omg = std::max(0.0, std::min(1.0,(ag*tabs1-bg)));

                    // Function for root finding: T_new = T_old + L_cond/C_p * (qt - qsat) + L_fus/C_p * (qp)
                    fff   = -tabs1 + tabs_array(i,j,k)
                            + lstarw*( qv_array(i,j,k)-qsattw)
                            + lstari*(qci_array(i,j,k)-qsatti)
                            + lstarp*qp_array(i,j,k);

                    // Derivative of function (T_new iterated on)
                    dfff  = -1.0
                            + dlstarw*( qv_array(i,j,k)-qsattw) - lstarw*dqsatw
                            + dlstari*(qci_array(i,j,k)-qsatti) - lstari*dqsati
                            + dlstarp*qp_array(i,j,k);

                    // Update the temperature
                    dtabs = -fff/dfff;
                    tabs1 = tabs1+dtabs;
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
                 qci_array(i,j,k) = qsatti;
                qcl_array(i,j,k) += delta_qc;
                qn_array(i,j,k)   = qcl_array(i,j,k) + qci_array(i,j,k);
                 qt_array(i,j,k)  =  qv_array(i,j,k) +  qn_array(i,j,k);

                // Partition the change in precipitating q
                qpr_array(i,j,k) = std::max(0.0,qpr_array(i,j,k) + delta_qi*omp);
                qps_array(i,j,k) = std::max(0.0,qps_array(i,j,k) + delta_qi*(1.0-omp)*(1.0-omg));
                qpg_array(i,j,k) = std::max(0.0,qpg_array(i,j,k) + delta_qi*(1.0-omp)*omg);
                qp_array(i,j,k)  = qpr_array(i,j,k) + qps_array(i,j,k) + qpg_array(i,j,k);

                if (k>50 && k<80 && i==0 && j==0) { //dtabs > 0.5 || omn != 1.0) {
                    amrex::Print() << "WTF: " << IntVect(i,j,k) << ' '
                                   << dbg << ' '
                                   << delta_qv << ' '
                                   << qsattw << ' '
                                   << delta_qi << ' '
                                   << qsatti << ' '
                                   << delta_qc << ' '
                                   << tabs1 - tabs_array(i,j,k) << "\n";
                    amrex::Print() << "CHECK: "
                                   << qv_array(i,j,k) << ' '
                                   << qcl_array(i,j,k) << ' '
                                   << qci_array(i,j,k) << ' '
                                   << theta_array(i,j,k) << ' '
                                   << getThgivenPandT(tabs_array(i,j,k), pres_array(i,j,k)*100, rdOcp) << ' '
                                   << getThgivenPandT(tabs1, pres_array(i,j,k)*100, rdOcp) <<  "\n";
                    amrex::Print() << "\n";
                    //exit(0);
                }

                // Update temperature
                tabs_array(i,j,k) = tabs1;

                // Update pressure
                pres_array(i,j,k) =  rho_array(i,j,k) * R_d * tabs_array(i,j,k)
                                   * (1.0 + R_v/R_d * qv_array(i,j,k));

                // Update theta from temperature
                theta_array(i,j,k) = getThgivenPandT(tabs_array(i,j,k), pres_array(i,j,k)*100, rdOcp);
           //}
        });
    }
}
