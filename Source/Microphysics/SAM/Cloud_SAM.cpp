#include "SAM.H"
#include "IndexDefines.H"
#include "TileNoZ.H"
#include "EOS.H"

using namespace amrex;

/**
 * Split cloud components according to saturation pressures; source theta from latent heat.
 */
void SAM::Cloud () {

    constexpr Real an = 1.0/(tbgmax-tbgmin);
    constexpr Real bn = tbgmin*an;

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

        auto   rho_array = mic_fab_vars[MicVar::rho]->array(mfi);
        auto  tabs_array = mic_fab_vars[MicVar::tabs]->array(mfi);
        auto theta_array = mic_fab_vars[MicVar::theta]->array(mfi);
        auto  pres_array = mic_fab_vars[MicVar::pres]->array(mfi);

        const auto& box3d = mfi.tilebox(IntVect(0), IntVect(0,0,1));

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // Saturation moisture fractions
            Real omn, domn;
            Real qsat, dqsat;
            Real qsatw, dqsatw;
            Real qsati, dqsati;

            // Newton iteration vars
            int niter;
            Real fff, dfff, dtabs;
            Real lstar, dlstar;
            Real lstarw, lstari;
            Real delta_qv, delta_qc, delta_qi;

            // NOTE: Conversion before iterations is necessary to
            //       convert cloud water to ice or vice versa.
            //       This ensures the omn splitting is enforced
            //       before the Newton iteration, which assumes it is.

            // Cloud ice not permitted (melt to form water)
            if(tabs_array(i,j,k) >= tbgmax) {
                omn = 1.0;
                delta_qi = qci_array(i,j,k);
                qci_array(i,j,k)   = 0.0;
                qcl_array(i,j,k)  += delta_qi;
                tabs_array(i,j,k) -= fac_fus * delta_qi;
                pres_array(i,j,k)  = rho_array(i,j,k) * R_d * tabs_array(i,j,k)
                                     * (1.0 + R_v/R_d * qv_array(i,j,k));
                theta_array(i,j,k) = getThgivenPandT(tabs_array(i,j,k), pres_array(i,j,k), rdOcp);
                pres_array(i,j,k) *= 0.01;
            }
            // Cloud water not permitted (freeze to form ice)
            else if(tabs_array(i,j,k) <= tbgmin) {
                omn = 0.0;
                delta_qc = qcl_array(i,j,k);
                qcl_array(i,j,k)   = 0.0;
                qci_array(i,j,k)  += delta_qc;
                tabs_array(i,j,k) += fac_fus * delta_qc;
                pres_array(i,j,k)  = rho_array(i,j,k) * R_d * tabs_array(i,j,k)
                                     * (1.0 + R_v/R_d * qv_array(i,j,k));
                theta_array(i,j,k) = getThgivenPandT(tabs_array(i,j,k), pres_array(i,j,k), rdOcp);
                pres_array(i,j,k) *= 0.01;
            }
            // Mixed cloud phase (split according to omn)
            else {
                omn = an*tabs_array(i,j,k)-bn;
                delta_qc = qcl_array(i,j,k) - qn_array(i,j,k) * omn;
                delta_qi = qci_array(i,j,k) - qn_array(i,j,k) * (1.0 - omn);
                qcl_array(i,j,k)   = qn_array(i,j,k) * omn;
                qci_array(i,j,k)   = qn_array(i,j,k) * (1.0 - omn);
                tabs_array(i,j,k) += fac_fus * delta_qc;
                pres_array(i,j,k)  = rho_array(i,j,k) * R_d * tabs_array(i,j,k)
                                     * (1.0 + R_v/R_d * qv_array(i,j,k));
                theta_array(i,j,k) = getThgivenPandT(tabs_array(i,j,k), pres_array(i,j,k), rdOcp);
                pres_array(i,j,k) *= 0.01;
            }

            // Initial guess for temperature & pressure
            Real tabs = tabs_array(i,j,k);
            Real pres = pres_array(i,j,k);

            // Saturation moisture fractions
            erf_qsatw(tabs, pres, qsatw);
            erf_qsati(tabs, pres, qsati);
            qsat = omn * qsatw  + (1.0-omn) * qsati;

            // We have enough total moisture to relax to equilibrium
            if (qt_array(i,j,k) > qsat) {
                niter = 0;
                dtabs = 1;
                //==================================================
                // Newton iteration to qv=qsat (cloud phase only)
                //==================================================
                do {
                    // Latent heats and their derivatives wrt to T
                    lstarw  = fac_cond;
                    lstari  = fac_fus;
                    domn    = 0.0;

                    // Saturation moisture fractions
                    erf_qsatw(tabs, pres, qsatw);
                    erf_qsati(tabs, pres, qsati);
                    erf_dtqsatw(tabs, pres, dqsatw);
                    erf_dtqsati(tabs, pres, dqsati);

                    // Cloud ice not permitted (condensation & fusion)
                    if(tabs >= tbgmax) {
                        omn   = 1.0;
                    }
                    // Cloud water not permitted (sublimation & fusion)
                    else if(tabs <= tbgmin) {
                        omn    = 0.0;
                        lstarw = fac_sub;
                    }
                    // Mixed cloud phase (condensation & fusion)
                    else {
                        omn   = an*tabs-bn;
                        domn  = an;
                    }

                    // Linear combination of each component
                    qsat   =  omn * qsatw  + (1.0-omn) * qsati;
                    dqsat  =  omn * dqsatw + (1.0-omn) * dqsati
                           + domn *  qsatw -     domn  *  qsati;
                    lstar  =  omn * lstarw + (1.0-omn) * lstari;
                    dlstar = domn * lstarw -     domn  * lstari;

                    // Function for root finding:
                    // 0 = -T_new + T_old + L_eff/C_p * (qv - qsat)
                    fff   = -tabs + tabs_array(i,j,k) +  lstar*(qv_array(i,j,k) - qsat);

                    // Derivative of function (T_new iterated on)
                    dfff  = -1.0 + dlstar*(qv_array(i,j,k) - qsat) - lstar*dqsat;

                    // Update the temperature
                    dtabs = -fff/dfff;
                    tabs  = tabs+dtabs;

                    // Update the pressure
                    pres = rho_array(i,j,k) * R_d * tabs
                         * (1.0 + R_v/R_d * qsat) * 0.01;

                    // Update iteration
                    niter = niter+1;
                } while(std::abs(dtabs) > tol && niter < 20);

                // Update qsat from last iteration (dq = dq/dt * dt)
                qsat = qsat + dqsat*dtabs;

                // Changes in each component
                delta_qv = qv_array(i,j,k) - qsat;
                delta_qc = std::max(-qcl_array(i,j,k), delta_qv * omn);
                delta_qi = std::max(-qci_array(i,j,k), delta_qv * (1.0-omn));

                // Partition the change in non-precipitating q
                 qv_array(i,j,k)  = qsat;
                qcl_array(i,j,k) += delta_qc;
                qci_array(i,j,k) += delta_qi;
                 qn_array(i,j,k)  = qcl_array(i,j,k) + qci_array(i,j,k);
                 qt_array(i,j,k)  =  qv_array(i,j,k) +  qn_array(i,j,k);

                // Update temperature
                tabs_array(i,j,k) = tabs;

                // Update pressure
                pres_array(i,j,k) = rho_array(i,j,k) * R_d * tabs_array(i,j,k)
                                    * (1.0 + R_v/R_d * qv_array(i,j,k));

                // Update theta from temperature
                theta_array(i,j,k) = getThgivenPandT(tabs_array(i,j,k), pres_array(i,j,k), rdOcp);

                // Pressure unit conversion
                pres_array(i,j,k) *= 0.01;

            }
            // We cannot blindly relax to qsat, but we can convert qc/qi -> qv
            else {
                // Changes in each component
                delta_qv = qcl_array(i,j,k) + qci_array(i,j,k);
                delta_qc = qcl_array(i,j,k);
                delta_qi = qci_array(i,j,k);

                // Partition the change in non-precipitating q
                 qv_array(i,j,k) += delta_qv;
                qcl_array(i,j,k)  = 0.0;
                qci_array(i,j,k)  = 0.0;
                 qn_array(i,j,k)  = 0.0;
                 qt_array(i,j,k)  = qv_array(i,j,k);

                // Update temperature (endothermic since we evap/sublime)
                tabs_array(i,j,k) -= fac_cond * delta_qc + fac_sub * delta_qi;

                // Update pressure
                pres_array(i,j,k) = rho_array(i,j,k) * R_d * tabs_array(i,j,k)
                                    * (1.0 + R_v/R_d * qv_array(i,j,k));

                // Update theta from temperature
                theta_array(i,j,k) = getThgivenPandT(tabs_array(i,j,k), pres_array(i,j,k), rdOcp);

                // Pressure unit conversion
                pres_array(i,j,k) *= 0.01;
            }
        });
    } // mfi
}
