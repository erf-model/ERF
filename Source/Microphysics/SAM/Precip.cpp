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
        auto pres_array = mic_fab_vars[MicVar::pres]->array(mfi);

        // Non-precipitating
        auto qv_array    = mic_fab_vars[MicVar::qv]->array(mfi);
        auto qcl_array   = mic_fab_vars[MicVar::qcl]->array(mfi);
        auto qci_array   = mic_fab_vars[MicVar::qci]->array(mfi);
        auto qn_array    = mic_fab_vars[MicVar::qn]->array(mfi);
        auto qt_array    = mic_fab_vars[MicVar::qt]->array(mfi);

        // Precipitating
        auto qpr_array   = mic_fab_vars[MicVar::qpr]->array(mfi);
        auto qps_array   = mic_fab_vars[MicVar::qps]->array(mfi);
        auto qpg_array   = mic_fab_vars[MicVar::qpg]->array(mfi);
        auto qp_array    = mic_fab_vars[MicVar::qp]->array(mfi);

        const auto& box3d = mfi.tilebox();

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            //------- Autoconversion/accretion
            Real tmp, qsatt;
            Real omn, omp, omg;

            Real qcc, qii, qpr, qps, qpg;
            Real dprc, dpsc, dpgc;
            Real dpri, dpsi, dpgi;

            Real dqc, dqi;
            Real dqpr, dqps, dqpg;

            Real autor, autos;
            Real accrr, accrcs, accris, accrcg, accrig;

            if (qn_array(i,j,k)+qp_array(i,j,k) > 0.0) {

                omn = std::max(0.0,std::min(1.0,(tabs_array(i,j,k)-tbgmin)*a_bg));
                omp = std::max(0.0,std::min(1.0,(tabs_array(i,j,k)-tprmin)*a_pr));
                omg = std::max(0.0,std::min(1.0,(tabs_array(i,j,k)-tgrmin)*a_gr));

                if (qn_array(i,j,k) > 0.0) {
                    qcc = qcl_array(i,j,k);
                    qii = qci_array(i,j,k);

                    qpr = qpr_array(i,j,k);
                    qps = qps_array(i,j,k);
                    qpg = qpg_array(i,j,k);

                    accrr  = 0.0;
                    accrcs = 0.0;
                    accris = 0.0;
                    accrcg = 0.0;
                    accrig = 0.0;

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

                    if (omp > 0.001) {
                        accrr = accrrc_t(k) * std::pow(qpr, powr1);
                    }

                    if (omp < 0.999 && omg < 0.999) {
                        tmp = pow(qps, pows1);
                        accrcs = accrsc_t(k) * tmp;
                        accris = accrsi_t(k) * tmp;
                    }

                    if (omp < 0.999 && omg > 0.001) {
                        tmp = pow(qpg, powg1);
                        accrcg = accrgc_t(k) * tmp;
                        accrig = accrgi_t(k) * tmp;
                    }

                    // Autoconversion & accretion
                    dprc = dtn * autor  * (qcc-qcw0);
                    dpsc = dtn * accrcs * qcc * std::pow(qps,(3.+b_snow)/4.);
                    dpgc = dtn * accrcg * qcc * std::pow(qpg,(3.+b_grau)/4.);

                    dpri = dtn * autos  * (qii-qci0);
                    dpsi = dtn * accris * qii * std::pow(qps,(3.+b_snow)/4.);
                    dpgi = dtn * accrig * qii * std::pow(qpg,(3.+b_grau)/4.);

                    dqpr = dprc + dpri;
                    dqps = dpsc + dpsi;
                    dqpg = dpgc + dpgi;

                    dqc  = dprc + dpsc + dpgc;
                    dqi  = dpri + dpsi + dpgi;

                    // Update the primitive state variables
                    qcl_array(i,j,k) = std::max(0.0,qcl_array(i,j,k) - dqc);
                    qci_array(i,j,k) = std::max(0.0,qci_array(i,j,k) - dqi);
                    qpr_array(i,j,k) = std::max(0.0,qpr_array(i,j,k) + dqpr);
                    qps_array(i,j,k) = std::max(0.0,qps_array(i,j,k) + dqps);
                    qpg_array(i,j,k) = std::max(0.0,qpg_array(i,j,k) + dqpg);

                    // Update the primitive derived vars
                    qn_array(i,j,k) = qcl_array(i,j,k) + qci_array(i,j,k);
                    qt_array(i,j,k) =  qv_array(i,j,k) +  qn_array(i,j,k);
                    qp_array(i,j,k) = qpr_array(i,j,k) + qps_array(i,j,k) + qpg_array(i,j,k);

                    // TODO: Reconcile the gcc operations from original formulation!

                } else if(qp_array(i,j,k) > qp_threshold && qn_array(i,j,k) == 0.0) {

                    erf_qsatw(tabs_array(i,j,k),pres_array(i,j,k),qsatt);

                    if(omp > 0.001) {
                        qpr  = qpr_array(i,j,k);
                        dqpr = evapr1_t(k)*sqrt(qpr) + evapr2_t(k)*pow(qpr,powr2);
                    }
                    if(omp < 0.999 && omg < 0.999) {
                        qps  = qps_array(i,j,k);
                        dqps = evaps1_t(k)*sqrt(qps) + evaps2_t(k)*pow(qps,pows2);
                    }
                    if(omp < 0.999 && omg > 0.001) {
                        qpg  = qpg_array(i,j,k);
                        dqpg = evapg1_t(k)*sqrt(qpg) + evapg2_t(k)*pow(qpg,powg2);
                    }

                    // Evaporation
                    dqpr *= dtn * (qv_array(i,j,k)/qsatt - 1.0);
                    dqps *= dtn * (qv_array(i,j,k)/qsatt - 1.0);
                    dqpg *= dtn * (qv_array(i,j,k)/qsatt - 1.0);

                    // Update the primitive state variables
                    qcl_array(i,j,k) = std::max(0.0,qcl_array(i,j,k) - dqpr);
                    qci_array(i,j,k) = std::max(0.0,qci_array(i,j,k) - dqps - dqpg);
                    qpr_array(i,j,k) = std::max(0.0,qpr_array(i,j,k) + dqpr);
                    qps_array(i,j,k) = std::max(0.0,qps_array(i,j,k) + dqps);
                    qpg_array(i,j,k) = std::max(0.0,qpg_array(i,j,k) + dqpg);

                    // Update the primitive derived vars
                    qn_array(i,j,k) = qcl_array(i,j,k) + qci_array(i,j,k);
                    qt_array(i,j,k) =  qv_array(i,j,k) +  qn_array(i,j,k);
                    qp_array(i,j,k) = qpr_array(i,j,k) + qps_array(i,j,k) + qpg_array(i,j,k);

                    // TODO: Latent heat source term for theta?

                } else {
                    // Update the primitive state variables
                    qcl_array(i,j,k) += qpr_array(i,j,k);
                    qci_array(i,j,k) += qps_array(i,j,k) + qpg_array(i,j,k);
                    qpr_array(i,j,k)  = 0.0;
                    qps_array(i,j,k)  = 0.0;
                    qpg_array(i,j,k)  = 0.0;

                    // Update the primitive derived vars
                    qn_array(i,j,k) = qcl_array(i,j,k) + qci_array(i,j,k);
                    qt_array(i,j,k) =  qv_array(i,j,k) +  qn_array(i,j,k);
                    qp_array(i,j,k) = qpr_array(i,j,k) + qps_array(i,j,k) + qpg_array(i,j,k);

                    // TODO: Latent heat source term for theta?

                }
            }
        });
    }
}


