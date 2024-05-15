#include "SAM.H"
#include "EOS.H"

using namespace amrex;

/**
 * Autoconversion (A30), Accretion (A28), Evaporation (A24)
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

    Real fac_cond = m_fac_cond;
    Real fac_sub  = m_fac_sub;
    Real fac_fus  = m_fac_fus;
    Real rdOcp    = m_rdOcp;

    Real eps = std::numeric_limits<Real>::epsilon();

    Real dtn = dt;

    // get the temperature, dentisy, theta, qt and qp from input
    for ( MFIter mfi(*(mic_fab_vars[MicVar::tabs]),TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto rho_array   = mic_fab_vars[MicVar::rho]->array(mfi);
        auto theta_array = mic_fab_vars[MicVar::theta]->array(mfi);
        auto tabs_array  = mic_fab_vars[MicVar::tabs]->array(mfi);
        auto pres_array  = mic_fab_vars[MicVar::pres]->array(mfi);

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
            Real omn, omp, omg;
            Real qsat, qsatw, qsati;

            Real qcc, qii, qpr, qps, qpg;
            Real dprc, dpsc, dpgc;
            Real dpsi, dpgi;

            Real dqc, dqca, dqi, dqia, dqp;
            Real dqpr, dqps, dqpg;

            Real autor, autos;
            Real accrcr, accrcs, accris, accrcg, accrig;

            // Work to be done for autoc/accr or evap
            if (qn_array(i,j,k)+qp_array(i,j,k) > 0.0) {
                omn = std::max(0.0,std::min(1.0,(tabs_array(i,j,k)-tbgmin)*a_bg));
                omp = std::max(0.0,std::min(1.0,(tabs_array(i,j,k)-tprmin)*a_pr));
                omg = std::max(0.0,std::min(1.0,(tabs_array(i,j,k)-tgrmin)*a_gr));

                qcc = qcl_array(i,j,k);
                qii = qci_array(i,j,k);

                qpr = qpr_array(i,j,k);
                qps = qps_array(i,j,k);
                qpg = qpg_array(i,j,k);

                //==================================================
                // Autoconversion (A30/A31) and accretion (A27)
                //==================================================
                if (qn_array(i,j,k) > 0.0) {
                    accrcr = 0.0;
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
                        accrcr = accrrc_t(k);
                    }

                    if (omp < 0.999 && omg < 0.999) {
                        accrcs = accrsc_t(k);
                        accris = accrsi_t(k);
                    }

                    if (omp < 0.999 && omg > 0.001) {
                        accrcg = accrgc_t(k);
                        accrig = accrgi_t(k);
                    }

                    // Autoconversion & accretion (sink for cloud comps)
                    dqca = dtn * autor  * (qcc-qcw0);
                    dprc = dtn * accrcr * qcc * std::pow(qpr, powr1);
                    dpsc = dtn * accrcs * qcc * std::pow(qps, pows1);
                    dpgc = dtn * accrcg * qcc * std::pow(qpg, powg1);

                    dqia = dtn * autos  * (qii-qci0);
                    dpsi = dtn * accris * qii * std::pow(qps, pows1);
                    dpgi = dtn * accrig * qii * std::pow(qpg, powg1);

                    // Rescale sinks to avoid negative cloud fractions
                    dqc  = dqca + dprc + dpsc + dpgc;
                    dqi  = dqia + dpsi + dpgi;
                    Real scalec = std::min(qcl_array(i,j,k),dqc) / (dqc + eps);
                    Real scalei = std::min(qci_array(i,j,k),dqi) / (dqi + eps);
                    dqca *= scalec; dprc *= scalec; dpsc *= scalec; dpgc *= scalec;
                    dqia *= scalei; dpsi *= scalei; dpgi *= scalei;
                    dqc   = dqca + dprc + dpsc + dpgc;
                    dqi   = dqia + dpsi + dpgi;

                    // NOTE: Autoconversion of cloud water and ice are sources
                    //       to qp, while accretion is a source to an individual
                    //       precipitating component (e.g., qpr/qps/qpg). So we
                    //       only split autoconversion with omega. The omega
                    //       splitting does imply a latent heat source.

                    // Partition formed precip componentss
                    dqpr = (dqca + dqia) * omp + dprc;
                    dqps = (dqca + dqia) * (1.0 - omp) * (1.0 - omg) + dpsc + dpsi;
                    dqpg = (dqca + dqia) * (1.0 - omp) * omg         + dpgc + dpgi;

                    // Update the primitive state variables
                    qcl_array(i,j,k) -= dqc;
                    qci_array(i,j,k) -= dqi;
                    qpr_array(i,j,k) += dqpr;
                    qps_array(i,j,k) += dqps;
                    qpg_array(i,j,k) += dqpg;

                    // Update the primitive derived vars
                    qn_array(i,j,k) = qcl_array(i,j,k) + qci_array(i,j,k);
                    qt_array(i,j,k) =  qv_array(i,j,k) +  qn_array(i,j,k);
                    qp_array(i,j,k) = qpr_array(i,j,k) + qps_array(i,j,k) + qpg_array(i,j,k);

                    // Latent heat source for theta
                    tabs_array(i,j,k) += fac_fus * ( dqca * (1.0 - omp) - dqia * omp );
                    pres_array(i,j,k)  = rho_array(i,j,k) * R_d * tabs_array(i,j,k)
                                       * (1.0 + R_v/R_d * qv_array(i,j,k));
                    theta_array(i,j,k) = getThgivenPandT(tabs_array(i,j,k), pres_array(i,j,k), rdOcp);
                    pres_array(i,j,k) /= 100.0;
                }

                //==================================================
                // Evaporation (A24)
                //==================================================
                erf_qsatw(tabs_array(i,j,k),pres_array(i,j,k),qsatw);
                erf_qsati(tabs_array(i,j,k),pres_array(i,j,k),qsati);
                qsat = qsatw * omn + qsati * (1.0-omn);
                if((qp_array(i,j,k) > 0.0) && (qv_array(i,j,k) < qsat)) {

                    dqpr = evapr1_t(k)*sqrt(qpr) + evapr2_t(k)*pow(qpr,powr2);
                    dqps = evaps1_t(k)*sqrt(qps) + evaps2_t(k)*pow(qps,pows2);
                    dqpg = evapg1_t(k)*sqrt(qpg) + evapg2_t(k)*pow(qpg,powg2);

                    // NOTE: This is always a sink for precipitating comps
                    //       since qv<qsat and thus (qv/qsat-1)<0. If we are
                    //       in a super-saturated state (qv>qsat) the Newton
                    //       iterations in Cloud() will have handled condensation.
                    dqpr *= -dtn * (qv_array(i,j,k)/qsat - 1.0);
                    dqps *= -dtn * (qv_array(i,j,k)/qsat - 1.0);
                    dqpg *= -dtn * (qv_array(i,j,k)/qsat - 1.0);

                    // Limit to avoid negative moisture fractions
                    dqpr = std::min(qpr_array(i,j,k),dqpr);
                    dqps = std::min(qps_array(i,j,k),dqps);
                    dqpg = std::min(qpg_array(i,j,k),dqpg);
                    dqp  = dqpr + dqps + dqpg;

                    // Update the primitive state variables
                     qv_array(i,j,k) += dqp;
                    qpr_array(i,j,k) -= dqpr;
                    qps_array(i,j,k) -= dqps;
                    qpg_array(i,j,k) -= dqpg;

                    // Update the primitive derived vars
                    qt_array(i,j,k) =  qv_array(i,j,k) +  qn_array(i,j,k);
                    qp_array(i,j,k) = qpr_array(i,j,k) + qps_array(i,j,k) + qpg_array(i,j,k);

                    // Latent heat source for theta
                    tabs_array(i,j,k) += fac_cond * dqpr + fac_sub * (dqps + dqpg);
                    pres_array(i,j,k)  = rho_array(i,j,k) * R_d * tabs_array(i,j,k)
                                       * (1.0 + R_v/R_d * qv_array(i,j,k));
                    theta_array(i,j,k) = getThgivenPandT(tabs_array(i,j,k), pres_array(i,j,k), rdOcp);
                    pres_array(i,j,k) /= 100.0;
                }
            }
        });
    }
}


