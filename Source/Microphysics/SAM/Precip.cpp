#include "SAM.H"

using namespace amrex;

/**
 * Compute Precipitation-related Microphysics quantities.
 */
void SAM::Precip () {

    Real powr1 = (3.0 + b_rain) / 4.0;
    Real powr2 = (5.0 + b_rain) / 8.0;
    Real pows1 = (3.0 + b_snow) / 4.0;
    Real pows2 = (5.0 + b_snow) / 8.0;
    Real powg1 = (3.0 + b_grau) / 4.0;
    Real powg2 = (5.0 + b_grau) / 8.0;

    auto accrrc_t  = accrrc.table();
    auto accrsc_t  = accrsc.table();
    auto accrsi_t  = accrsi.table();
    auto accrgc_t  = accrgc.table();
    auto accrgi_t  = accrgi.table();
    auto coefice_t = coefice.table();
    auto evapr1_t  = evapr1.table();
    auto evapr2_t  = evapr2.table();
    auto evaps1_t  = evaps1.table();
    auto evaps2_t  = evaps2.table();
    auto evapg1_t  = evapg1.table();
    auto evapg2_t  = evapg2.table();
    auto pres1d_t  = pres1d.table();

    auto qt   = mic_fab_vars[MicVar::qt];
    auto qp   = mic_fab_vars[MicVar::qp];
    auto qn   = mic_fab_vars[MicVar::qn];
    auto tabs = mic_fab_vars[MicVar::tabs];

    Real dtn = dt;

    // get the temperature, dentisy, theta, qt and qp from input
    for ( MFIter mfi(*tabs,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto tabs_array = mic_fab_vars[MicVar::tabs]->array(mfi);
        auto qn_array   = mic_fab_vars[MicVar::qn]->array(mfi);
        auto qt_array   = mic_fab_vars[MicVar::qt]->array(mfi);
        auto qp_array   = mic_fab_vars[MicVar::qp]->array(mfi);

        const auto& box3d = mfi.tilebox();

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            //------- Autoconversion/accretion
            Real omn, omp, omg, qcc, qii, autor, autos, accrr, qrr, accrcs, accris,
                qss, accrcg, accrig, tmp, qgg, dq, qsatt, qsat;

            if (qn_array(i,j,k)+qp_array(i,j,k) > 0.0) {
                omn = std::max(0.0,std::min(1.0,(tabs_array(i,j,k)-tbgmin)*a_bg));
                omp = std::max(0.0,std::min(1.0,(tabs_array(i,j,k)-tprmin)*a_pr));
                omg = std::max(0.0,std::min(1.0,(tabs_array(i,j,k)-tgrmin)*a_gr));

                if (qn_array(i,j,k) > 0.0) {
                    qcc = qn_array(i,j,k) * omn;
                    qii = qn_array(i,j,k) * (1.0-omn);

                    if (qcc > qcw0) {
                        autor = alphaelq;
                    } else {
                        autor = 0.0;
                    }

                    if (qii > qci0) {
                        autos = betaelq*coefice_t(k);
                    } else {
                        autos = 0.0;
                    }

                    accrr = 0.0;
                    if (omp > 0.001) {
                        qrr = qp_array(i,j,k) * omp;
                        accrr = accrrc_t(k) * std::pow(qrr, powr1);
                    }

                    accrcs = 0.0;
                    accris = 0.0;

                    if (omp < 0.999 && omg < 0.999) {
                        qss = qp_array(i,j,k) * (1.0-omp)*(1.0-omg);
                        tmp = pow(qss, pows1);
                        accrcs = accrsc_t(k) * tmp;
                        accris = accrsi_t(k) * tmp;
                    }
                    accrcg = 0.0;
                    accrig = 0.0;
                    if (omp < 0.999 && omg > 0.001) {
                        qgg = qp_array(i,j,k) * (1.0-omp)*omg;
                        tmp = pow(qgg, powg1);
                        accrcg = accrgc_t(k) * tmp;
                        accrig = accrgi_t(k) * tmp;
                    }
                    qcc = (qcc+dtn*autor*qcw0)/(1.0+dtn*(accrr+accrcs+accrcg+autor));
                    qii = (qii+dtn*autos*qci0)/(1.0+dtn*(accris+accrig+autos));
                    dq = dtn *(accrr*qcc + autor*(qcc-qcw0)+(accris+accrig)*qii + (accrcs+accrcg)*qcc + autos*(qii-qci0));
                    dq = std::min(dq,qn_array(i,j,k));
                    qt_array(i,j,k) = qt_array(i,j,k) - dq;
                    qp_array(i,j,k) = qp_array(i,j,k) + dq;
                    qn_array(i,j,k) = qn_array(i,j,k) - dq;

                } else if(qp_array(i,j,k) > qp_threshold && qn_array(i,j,k) == 0.0) {

                    qsatt = 0.0;
                    if(omn > 0.001) {
                        erf_qsatw(tabs_array(i,j,k),pres1d_t(k),qsat);
                        qsatt = qsatt + omn*qsat;
                    }
                    if(omn < 0.999) {
                        erf_qsati(tabs_array(i,j,k),pres1d_t(k),qsat);
                        qsatt = qsatt + (1.-omn)*qsat;
                    }
                    dq = 0.0;
                    if(omp > 0.001) {
                        qrr = qp_array(i,j,k) * omp;
                        dq = dq + evapr1_t(k)*sqrt(qrr) + evapr2_t(k)*pow(qrr,powr2);
                    }
                    if(omp < 0.999 && omg < 0.999) {
                        qss = qp_array(i,j,k) * (1.0-omp)*(1.0-omg);
                        dq = dq + evaps1_t(k)*sqrt(qss) + evaps2_t(k)*pow(qss,pows2);
                    }
                    if(omp < 0.999 && omg > 0.001) {
                        qgg = qp_array(i,j,k) * (1.0-omp)*omg;
                        dq = dq + evapg1_t(k)*sqrt(qgg) + evapg2_t(k)*pow(qgg,powg2);
                    }
                    dq = dq * dtn * (qt_array(i,j,k) / qsatt-1.0);
                    dq = std::max(-0.5*qp_array(i,j,k),dq);
                    qt_array(i,j,k) = qt_array(i,j,k) - dq;
                    qp_array(i,j,k) = qp_array(i,j,k) + dq;

                } else {
                    qt_array(i,j,k) = qt_array(i,j,k) + qp_array(i,j,k);
                    qp_array(i,j,k) = 0.0;
                }
            }
            dq = qp_array(i,j,k);
            qp_array(i,j,k) = max(0.0,qp_array(i,j,k));
            qt_array(i,j,k) = qt_array(i,j,k) + (dq-qp_array(i,j,k));
        });
    }
}


