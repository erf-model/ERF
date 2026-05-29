#include "ERF_SAM.H"
#include "ERF_EOS.H"
#include "ERF_SAMUtils.H"

using namespace amrex;

/**
 * Autoconversion (A30), Accretion (A28), Evaporation (A24)
 */
void
SAM::Precip (const SolverChoice& sc)
{
    if (sam_is_no_precip(sc.moisture_type)) return;

    Real powr1 = (three + b_rain) / Real(4.0);
    Real powr2 = (Real(5.0) + b_rain) / Real(8.0);
    Real pows1 = (three + b_snow) / Real(4.0);
    Real pows2 = (Real(5.0) + b_snow) / Real(8.0);
    Real powg1 = (three + b_grau) / Real(4.0);
    Real powg2 = (Real(5.0) + b_grau) / Real(8.0);

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

    int SAM_moisture_type = 1;
    if (sc.moisture_type == MoistureType::SAM_NoIce) {
        SAM_moisture_type = 2;
    }

    // get the temperature, density, theta, qt and qp from input
    for ( MFIter mfi(*(mic_fab_vars[MicVar::tabs]),TilingIfNotGPU()); mfi.isValid(); ++mfi) {
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

        auto tbx = mfi.tilebox();

        ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            //------- Autoconversion/accretion
            Real omn, omp, omg;
            Real qsat, qsatw, qsati;

            Real qcc, qii, qpr, qps, qpg;
            Real dqc, dqca, dqi, dqia, dqp;
            Real dqpr, dqps, dqpg;

            // Work to be done for autoc/accr or evap
            if (qn_array(i,j,k)+qp_array(i,j,k) > zero) {
                if (SAM_moisture_type == 2) {
                    omn = Real(1);
                    omp = Real(1);
                    omg = Real(0);
                } else {
                    omn = sam_cloud_liquid_fraction(SAM_moisture_type,
                                                    tabs_array(i,j,k), a_bg, tbgmin * a_bg);
                    omp = sam_precip_rain_fraction(SAM_moisture_type,
                                                   tabs_array(i,j,k));
                    omg = sam_graupel_fraction(SAM_moisture_type,
                                               tabs_array(i,j,k));
                }

                qcc = qcl_array(i,j,k);
                qii = qci_array(i,j,k);

                qpr = qpr_array(i,j,k);
                qps = qps_array(i,j,k);
                qpg = qpg_array(i,j,k);

                //==================================================
                // Autoconversion (A30/A31) and accretion (A27)
                //==================================================
                if (qn_array(i,j,k) > zero) {
                    SAMPrecipSources source_terms =
                        sam_autoconversion_rates(dtn, qcc, qii, coefice_t(k));
                    const SAMPrecipSources accretion_terms =
                        sam_accretion_rates(dtn, qcc, qii, qpr, qps, qpg,
                                            powr1, pows1, powg1,
                                            omp, omg,
                                            accrrc_t(k), accrsc_t(k), accrsi_t(k),
                                            accrgc_t(k), accrgi_t(k));
                    source_terms.dprc = accretion_terms.dprc;
                    source_terms.dpsc = accretion_terms.dpsc;
                    source_terms.dpgc = accretion_terms.dpgc;
                    source_terms.dpsi = accretion_terms.dpsi;
                    source_terms.dpgi = accretion_terms.dpgi;
                    source_terms = sam_rescale_cloud_sinks(qcl_array(i,j,k), qci_array(i,j,k),
                                                           eps, source_terms);

                    // NOTE: Autoconversion of cloud water and ice are sources
                    //       to qp, while accretion is a source to an individual
                    //       precipitating component (e.g., qpr/qps/qpg). So we
                    //       only split autoconversion with omega. The omega
                    //       splitting does imply a latent heat source.

                    // Partition formed precip componentss
                    source_terms = sam_partition_autoconverted_precip(source_terms, omp, omg);
                    dqca = source_terms.dqca;
                    dqia = source_terms.dqia;
                    dqc = source_terms.dqc;
                    dqi = source_terms.dqi;
                    dqpr = source_terms.dqpr;
                    dqps = source_terms.dqps;
                    dqpg = source_terms.dqpg;

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

                    // Update temperature
                    tabs_array(i,j,k) += fac_fus * ( dqca * (one - omp) - dqia * omp );

                    // Update theta
                    theta_array(i,j,k) = sam_theta_from_stored_mbar_converted_to_pa(tabs_array(i,j,k),
                                                                                     pres_array(i,j,k), rdOcp);
                }

                //==================================================
                // Evaporation (A24)
                //==================================================
                erf_qsatw(tabs_array(i,j,k),pres_array(i,j,k),qsatw);
                erf_qsati(tabs_array(i,j,k),pres_array(i,j,k),qsati);
                qsat = sam_mixed_qsat(omn, qsatw, qsati);
                if((qp_array(i,j,k) > zero) && (qv_array(i,j,k) < qsat)) {

                    SAMPrecipSources evaporation_terms =
                        sam_precip_evaporation_rates(qpr, qps, qpg,
                                                     powr2, pows2, powg2,
                                                     evapr1_t(k), evapr2_t(k),
                                                     evaps1_t(k), evaps2_t(k),
                                                     evapg1_t(k), evapg2_t(k));

                    // NOTE: This is always a sink for precipitating comps
                    //       since qv<qsat and thus (1 - qv/qsat)>zero If we are
                    //       in a super-saturated state (qv>qsat) the Newton
                    //       iterations in Cloud() will have handled condensation.
                    evaporation_terms.dqpr *= dtn * (one - qv_array(i,j,k)/qsat);
                    evaporation_terms.dqps *= dtn * (one - qv_array(i,j,k)/qsat);
                    evaporation_terms.dqpg *= dtn * (one - qv_array(i,j,k)/qsat);

                    // Limit to avoid negative moisture fractions
                    evaporation_terms = sam_apply_precip_evaporation_limiter(qpr_array(i,j,k),
                                                                             qps_array(i,j,k),
                                                                             qpg_array(i,j,k),
                                                                             evaporation_terms);
                    dqpr = evaporation_terms.dqpr;
                    dqps = evaporation_terms.dqps;
                    dqpg = evaporation_terms.dqpg;
                    dqp  = evaporation_terms.dqp;

                    // Update the primitive state variables
                     qv_array(i,j,k) += dqp;
                    qpr_array(i,j,k) -= dqpr;
                    qps_array(i,j,k) -= dqps;
                    qpg_array(i,j,k) -= dqpg;

                    // Update the primitive derived vars
                    qt_array(i,j,k) =  qv_array(i,j,k) +  qn_array(i,j,k);
                    qp_array(i,j,k) = qpr_array(i,j,k) + qps_array(i,j,k) + qpg_array(i,j,k);

                    // Update temperature
                    tabs_array(i,j,k) -= fac_cond * dqpr + fac_sub * (dqps + dqpg);

                    // Update theta
                    theta_array(i,j,k) = sam_theta_from_stored_mbar_converted_to_pa(tabs_array(i,j,k),
                                                                                     pres_array(i,j,k), rdOcp);
                }
            }
        });
    }
}
