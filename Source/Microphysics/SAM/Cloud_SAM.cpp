#include "SAM.H"
#include "IndexDefines.H"
#include "TileNoZ.H"
#include "EOS.H"

using namespace amrex;

/**
 * Split cloud components according to saturation pressures; source theta from latent heat.
 */
void
SAM::Cloud (const SolverChoice& sc)
{

    constexpr Real an = 1.0/(tbgmax-tbgmin);
    constexpr Real bn = tbgmin*an;

    Real fac_cond = m_fac_cond;
    Real fac_sub  = m_fac_sub;
    Real fac_fus  = m_fac_fus;
    Real rdOcp    = m_rdOcp;

    int SAM_moisture_type = 1;
    if (sc.moisture_type == MoistureType::SAM_NoIce ||
        sc.moisture_type == MoistureType::SAM_NoPrecip_NoIce) {
        SAM_moisture_type = 2;
    }

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
            Real omn;
            Real qsat;
            Real qsatw;
            Real qsati;

            // Newton iteration vars
            Real delta_qv, delta_qc, delta_qi;

            // NOTE: Conversion before iterations is necessary to
            //       convert cloud water to ice or vice versa.
            //       This ensures the omn splitting is enforced
            //       before the Newton iteration, which assumes it is.

            omn = 1.0;
            if (SAM_moisture_type == 1){
                // Cloud ice not permitted (melt to form water)
                if (tabs_array(i,j,k) >= tbgmax) {
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
                else if (tabs_array(i,j,k) <= tbgmin) {
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
            }
            else if (SAM_moisture_type == 2)
            {
                // No ice. ie omn = 1.0
                delta_qc = qcl_array(i,j,k) - qn_array(i,j,k);
                delta_qi = 0.0;
                qcl_array(i,j,k)   = qn_array(i,j,k);
                qci_array(i,j,k)   = 0.0;
                tabs_array(i,j,k) += fac_cond * delta_qc;
                pres_array(i,j,k)  = rho_array(i,j,k) * R_d * tabs_array(i,j,k)
                                     * (1.0 + R_v/R_d * qv_array(i,j,k));
                theta_array(i,j,k) = getThgivenPandT(tabs_array(i,j,k), pres_array(i,j,k), rdOcp);
                pres_array(i,j,k) *= 0.01;
            }

            // Saturation moisture fractions
            erf_qsatw(tabs_array(i,j,k), pres_array(i,j,k), qsatw);
            erf_qsati(tabs_array(i,j,k), pres_array(i,j,k), qsati);
            qsat = omn * qsatw  + (1.0-omn) * qsati;

            // We have enough total moisture to relax to equilibrium
            if (qt_array(i,j,k) > qsat) {

                // Update temperature
                tabs_array(i,j,k) = NewtonIterSat(i, j, k   , SAM_moisture_type   ,
                                                  fac_cond  , fac_fus   , fac_sub ,
                                                  an        , bn        ,
                                                  tabs_array, pres_array,
                                                  qv_array  , qcl_array  , qci_array,
                                                  qn_array  , qt_array);

                // Update theta
                theta_array(i,j,k) = getThgivenPandT(tabs_array(i,j,k), 100.0*pres_array(i,j,k), rdOcp);

            //
            // We cannot blindly relax to qsat, but we can convert qc/qi -> qv.
            // The concept here is that if we put all the moisture into qv and modify
            // the temperature, we can then check if qv > qsat occurs (for final T/P/qv).
            // If the reduction in T/qsat and increase in qv does trigger the
            // aforementioned condition, we can do Newton iteration to drive qv = qsat.
            //
            } else {

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

                // NOTE: delta_qc & delta_qi are positive!
                // Update temperature (endothermic since we evap/sublime)
                tabs_array(i,j,k) -= fac_cond * delta_qc + fac_sub * delta_qi;

                // Update theta
                theta_array(i,j,k) = getThgivenPandT(tabs_array(i,j,k), 100.0*pres_array(i,j,k), rdOcp);

                // Verify assumption that qv > qsat does not occur
                erf_qsatw(tabs_array(i,j,k), pres_array(i,j,k), qsatw);
                erf_qsati(tabs_array(i,j,k), pres_array(i,j,k), qsati);
                qsat = omn * qsatw  + (1.0-omn) * qsati;
                if (qt_array(i,j,k) > qsat) {

                    // Update temperature
                    tabs_array(i,j,k) = NewtonIterSat(i, j, k   , SAM_moisture_type   ,
                                                      fac_cond  , fac_fus   , fac_sub ,
                                                      an        , bn        ,
                                                      tabs_array, pres_array,
                                                      qv_array  , qcl_array  , qci_array,
                                                      qn_array  , qt_array);

                    // Update theta
                    theta_array(i,j,k) = getThgivenPandT(tabs_array(i,j,k), 100.0*pres_array(i,j,k), rdOcp);

                }
            }
        });
    } // mfi
}
