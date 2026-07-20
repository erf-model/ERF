#include "ERF_SAM.H"
#include "ERF_IndexDefines.H"
#include "ERF_SAMUtils.H"
#include "ERF_TileNoZ.H"
#include "ERF_EOS.H"

using namespace amrex;

/**
 * Split cloud components according to saturation pressures; source theta from latent heat.
 */
void
SAM::Cloud (const SolverChoice& sc)
{
    if (!m_do_cond) { return; }

    constexpr Real an = one/(tbgmax-tbgmin);
    constexpr Real bn = tbgmin*an;

    Real fac_cond = m_fac_cond;
    Real fac_sub  = m_fac_sub;
    Real fac_fus  = m_fac_fus;
    Real rdOcp    = m_rdOcp;

    int SAM_moisture_type = 1;
    if (sam_is_no_ice(sc.moisture_type)) {
        SAM_moisture_type = 2;
    }

    for ( MFIter mfi(*(mic_fab_vars[MicVar::tabs]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto  qt_array = mic_fab_vars[MicVar::qt]->array(mfi);
        auto  qn_array = mic_fab_vars[MicVar::qn]->array(mfi);
        auto  qv_array = mic_fab_vars[MicVar::qv]->array(mfi);
        auto qcl_array = mic_fab_vars[MicVar::qcl]->array(mfi);
        auto qci_array = mic_fab_vars[MicVar::qci]->array(mfi);

        auto  tabs_array = mic_fab_vars[MicVar::tabs]->array(mfi);
        auto theta_array = mic_fab_vars[MicVar::theta]->array(mfi);
        auto  pres_array = mic_fab_vars[MicVar::pres]->array(mfi);

        auto tbx = mfi.tilebox();

        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
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

            const SAMCloudPhaseChange phase_change =
                sam_partition_cloud_phase(SAM_moisture_type,
                                          tabs_array(i,j,k),
                                          qn_array(i,j,k),
                                          qcl_array(i,j,k),
                                          qci_array(i,j,k),
                                          fac_cond, fac_fus,
                                          an, bn);
            omn = phase_change.omn;
            delta_qc = phase_change.delta_qc;
            delta_qi = phase_change.delta_qi;
            qcl_array(i,j,k) = phase_change.qcl;
            qci_array(i,j,k) = phase_change.qci;
            tabs_array(i,j,k) = phase_change.tabs;
            // Cloud adjustment updates tabs under the held-pressure SAM source
            // convention, then refreshes theta from the same stored pressure.
            theta_array(i,j,k) = sam_theta_from_stored_mbar_converted_to_pa(tabs_array(i,j,k),
                                                                             pres_array(i,j,k), rdOcp);

            // Saturation moisture fractions
            erf_qsatw(tabs_array(i,j,k), pres_array(i,j,k), qsatw);
            erf_qsati(tabs_array(i,j,k), pres_array(i,j,k), qsati);
            qsat = sam_mixed_qsat(omn, qsatw, qsati);

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
                theta_array(i,j,k) = sam_theta_from_stored_mbar_converted_to_pa(tabs_array(i,j,k),
                                                                                 pres_array(i,j,k), rdOcp);

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
                qcl_array(i,j,k)  = zero;
                qci_array(i,j,k)  = zero;
                 qn_array(i,j,k)  = zero;
                 qt_array(i,j,k)  = qv_array(i,j,k);

                // Update temperature (endothermic since we evap/sublime)
                tabs_array(i,j,k) -= fac_cond * delta_qc + fac_sub * delta_qi;

                // Update theta
                theta_array(i,j,k) = sam_theta_from_stored_mbar_converted_to_pa(tabs_array(i,j,k),
                                                                                 pres_array(i,j,k), rdOcp);

                // Verify assumption that qv > qsat does not occur
                erf_qsatw(tabs_array(i,j,k), pres_array(i,j,k), qsatw);
                erf_qsati(tabs_array(i,j,k), pres_array(i,j,k), qsati);
                qsat = sam_mixed_qsat(omn, qsatw, qsati);
                if (qt_array(i,j,k) > qsat) {

                    // Update temperature
                    tabs_array(i,j,k) = NewtonIterSat(i, j, k   , SAM_moisture_type   ,
                                                      fac_cond  , fac_fus   , fac_sub ,
                                                      an        , bn        ,
                                                      tabs_array, pres_array,
                                                      qv_array  , qcl_array  , qci_array,
                                                      qn_array  , qt_array);

                    // Update theta
                    theta_array(i,j,k) = sam_theta_from_stored_mbar_converted_to_pa(tabs_array(i,j,k),
                                                                                     pres_array(i,j,k), rdOcp);

                }
            }
        });
    } // mfi
}
