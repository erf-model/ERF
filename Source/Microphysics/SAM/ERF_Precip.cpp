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

    const SAMPrecipConfig precip_config{
        SAM_moisture_type,
        true,
        dtn,
        rdOcp,
        fac_cond,
        fac_fus,
        fac_sub,
        eps,
        powr1,
        pows1,
        powg1,
        powr2,
        pows2,
        powg2};
    // get the temperature, density, theta, qt and qp from input
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

        auto tbx = mfi.tilebox();

        ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            SAMCellState state{};
            state.rho = rho_array(i,j,k);
            state.theta = theta_array(i,j,k);
            state.tabs = tabs_array(i,j,k);
            state.pres_mbar = pres_array(i,j,k);
            state.qv = qv_array(i,j,k);
            state.qcl = qcl_array(i,j,k);
            state.qci = qci_array(i,j,k);
            state.qn = qn_array(i,j,k);
            state.qt = qt_array(i,j,k);
            state.qpr = qpr_array(i,j,k);
            state.qps = qps_array(i,j,k);
            state.qpg = qpg_array(i,j,k);
            state.qp = qp_array(i,j,k);

            const SAMCoefficientRow coeffs{
                accrrc_t(k),
                accrsi_t(k),
                accrsc_t(k),
                coefice_t(k),
                evaps1_t(k),
                evaps2_t(k),
                accrgi_t(k),
                accrgc_t(k),
                evapg1_t(k),
                evapg2_t(k),
                evapr1_t(k),
                evapr2_t(k)};

            // The scalar helper owns the local held-pressure source update;
            // this public kernel handles state staging and coefficient lookup.
            state = sam_precip_cell_update(state, coeffs, precip_config);

            theta_array(i,j,k) = state.theta;
            tabs_array(i,j,k) = state.tabs;
            qv_array(i,j,k) = state.qv;
            qcl_array(i,j,k) = state.qcl;
            qci_array(i,j,k) = state.qci;
            qn_array(i,j,k) = state.qn;
            qt_array(i,j,k) = state.qt;
            qpr_array(i,j,k) = state.qpr;
            qps_array(i,j,k) = state.qps;
            qpg_array(i,j,k) = state.qpg;
            qp_array(i,j,k) = state.qp;
        });
    }
}
