#ifndef ERF_SAM_UTILS_H_
#define ERF_SAM_UTILS_H_

#include <cmath>
#include <limits>

#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>

#include "ERF_Constants.H"
#include "ERF_DataStruct.H"
#include "ERF_EOS.H"
#include "ERF_IndexDefines.H"
#include "ERF_MicrophysicsUtils.H"

struct SAMPhaseFractions {
    amrex::Real omn{one};
    amrex::Real omp{one};
    amrex::Real omg{zero};
};

struct SAMCloudPhaseChange {
    amrex::Real omn{one};
    amrex::Real delta_qv{zero};
    amrex::Real delta_qc{zero};
    amrex::Real delta_qi{zero};
    amrex::Real qv{zero};
    amrex::Real qcl{zero};
    amrex::Real qci{zero};
    amrex::Real qn{zero};
    amrex::Real qt{zero};
    amrex::Real tabs{zero};
};

struct SAMPrecipSources {
    amrex::Real dqca{zero};
    amrex::Real dprc{zero};
    amrex::Real dpsc{zero};
    amrex::Real dpgc{zero};
    amrex::Real dqia{zero};
    amrex::Real dpsi{zero};
    amrex::Real dpgi{zero};
    amrex::Real dqc{zero};
    amrex::Real dqi{zero};
    amrex::Real dqpr{zero};
    amrex::Real dqps{zero};
    amrex::Real dqpg{zero};
    amrex::Real dqp{zero};
};

struct SAMPrimitiveCell {
    amrex::Real rho{zero};
    amrex::Real theta{zero};
    amrex::Real tabs{zero};
    amrex::Real pres_mbar{zero};
    amrex::Real qv{zero};
    amrex::Real qcl{zero};
    amrex::Real qci{zero};
    amrex::Real qn{zero};
    amrex::Real qt{zero};
    amrex::Real qpr{zero};
    amrex::Real qps{zero};
    amrex::Real qpg{zero};
    amrex::Real qp{zero};
};

struct SAMCoefficientRow {
    amrex::Real accrrc{zero};
    amrex::Real accrsi{zero};
    amrex::Real accrsc{zero};
    amrex::Real coefice{zero};
    amrex::Real evaps1{zero};
    amrex::Real evaps2{zero};
    amrex::Real accrgi{zero};
    amrex::Real accrgc{zero};
    amrex::Real evapg1{zero};
    amrex::Real evapg2{zero};
    amrex::Real evapr1{zero};
    amrex::Real evapr2{zero};
};

struct SAMPrecipFaceState {
    amrex::Real rho_avg{zero};
    amrex::Real tabs_avg{zero};
    amrex::Real qp_avg{zero};
    amrex::Real omp{one};
    amrex::Real omg{zero};
    amrex::Real qrr{zero};
    amrex::Real qss{zero};
    amrex::Real qgg{zero};
};

struct SAMFaceState {
    amrex::Real rho_avg{zero};
    amrex::Real tabs_avg{zero};
    amrex::Real qci_avg{zero};
    amrex::Real qp_avg{zero};
};

struct SAMSurfaceAccumulation {
    amrex::Real rain{zero};
    amrex::Real snow{zero};
    amrex::Real graupel{zero};
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_mbar_to_pa (const amrex::Real pres_mbar) noexcept
{
    // SAM stores pressure in mbar. EOS and theta helpers expect Pa.
    return amrex::Real(100.0) * pres_mbar;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_pa_to_mbar (const amrex::Real pres_pa) noexcept
{
    return amrex::Real(0.01) * pres_pa;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_theta_from_stored_mbar_converted_to_pa (const amrex::Real tabs,
                                                        const amrex::Real pres_mbar,
                                                        const amrex::Real rdOcp) noexcept
{
    // SAM stores pressure in mbar on this path. EOS helpers expect Pa.
    return getThgivenTandP(tabs, sam_mbar_to_pa(pres_mbar), rdOcp);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_theta_from_stored_pressure_as_pa (const amrex::Real tabs,
                                                  const amrex::Real stored_pres,
                                                  const amrex::Real rdOcp) noexcept
{
    // This SAM path passes the stored pressure value as though it were already Pa.
    return getThgivenTandP(tabs, stored_pres, rdOcp);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMPrimitiveCell sam_cons_to_primitive (const amrex::Real rho,
                                        const amrex::Real rho_theta,
                                        const amrex::Real rho_qv,
                                        const amrex::Real rho_qcl,
                                        const amrex::Real rho_qci,
                                        const amrex::Real rho_qpr,
                                        const amrex::Real rho_qps,
                                        const amrex::Real rho_qpg) noexcept
{
    // Copy-in clips moisture and stores pressure in mbar.
    SAMPrimitiveCell result{};
    result.rho = rho;
    result.theta = rho_theta / rho;
    result.qv = amrex::max(amrex::Real(0.0), rho_qv / rho);
    result.qcl = amrex::max(amrex::Real(0.0), rho_qcl / rho);
    result.qci = amrex::max(amrex::Real(0.0), rho_qci / rho);
    result.qn = result.qcl + result.qci;
    result.qt = result.qv + result.qn;
    result.qpr = amrex::max(amrex::Real(0.0), rho_qpr / rho);
    result.qps = amrex::max(amrex::Real(0.0), rho_qps / rho);
    result.qpg = amrex::max(amrex::Real(0.0), rho_qpg / rho);
    result.qp = result.qpr + result.qps + result.qpg;
    result.tabs = getTgivenRandRTh(rho, rho_theta, result.qv);
    result.pres_mbar = sam_pa_to_mbar(getPgivenRTh(rho_theta, result.qv));
    return result;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void sam_primitive_to_cons (const SAMPrimitiveCell& primitive,
                            const amrex::Array4<amrex::Real>& states_arr,
                            const int i,
                            const int j,
                            const int k) noexcept
{
    // Copy-out clips moisture before writing conserved state.
    states_arr(i,j,k,RhoTheta_comp) = primitive.rho * primitive.theta;
    states_arr(i,j,k,RhoQ1_comp) = primitive.rho * amrex::max(amrex::Real(0.0), primitive.qv);
    states_arr(i,j,k,RhoQ2_comp) = primitive.rho * amrex::max(amrex::Real(0.0), primitive.qcl);
    states_arr(i,j,k,RhoQ3_comp) = primitive.rho * amrex::max(amrex::Real(0.0), primitive.qci);
    states_arr(i,j,k,RhoQ4_comp) = primitive.rho * amrex::max(amrex::Real(0.0), primitive.qpr);
    states_arr(i,j,k,RhoQ5_comp) = primitive.rho * amrex::max(amrex::Real(0.0), primitive.qps);
    states_arr(i,j,k,RhoQ6_comp) = primitive.rho * amrex::max(amrex::Real(0.0), primitive.qpg);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool sam_is_no_ice (const MoistureType moisture_type) noexcept
{
    return moisture_type == MoistureType::SAM_NoIce ||
           moisture_type == MoistureType::SAM_NoPrecip_NoIce;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool sam_is_no_precip (const MoistureType moisture_type) noexcept
{
    return moisture_type == MoistureType::SAM_NoPrecip_NoIce;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_cloud_liquid_fraction (const int SAM_moisture_type,
                                       const amrex::Real tabs,
                                       const amrex::Real an,
                                       const amrex::Real bn) noexcept
{
    if (SAM_moisture_type == 2) {
        return one;
    }
    if (tabs >= tbgmax) {
        return one;
    }
    if (tabs <= tbgmin) {
        return zero;
    }
    return an * tabs - bn;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_precip_rain_fraction (const int SAM_moisture_type,
                                      const amrex::Real tabs) noexcept
{
    if (SAM_moisture_type == 2) {
        return one;
    }
    return amrex::max(amrex::Real(0.0), amrex::min(amrex::Real(1.0), (tabs - tprmin) * a_pr));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_graupel_fraction (const int SAM_moisture_type,
                                  const amrex::Real tabs) noexcept
{
    if (SAM_moisture_type == 2) {
        return zero;
    }
    return amrex::max(amrex::Real(0.0), amrex::min(amrex::Real(1.0), (tabs - tgrmin) * a_gr));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMCloudPhaseChange sam_partition_cloud_phase (const int SAM_moisture_type,
                                               const amrex::Real tabs,
                                               const amrex::Real qn,
                                               const amrex::Real qcl,
                                               const amrex::Real qci,
                                               const amrex::Real fac_cond,
                                               const amrex::Real fac_fus,
                                               const amrex::Real an,
                                               const amrex::Real bn) noexcept
{
    SAMCloudPhaseChange result{};
    result.qcl = qcl;
    result.qci = qci;
    result.tabs = tabs;
    result.omn = sam_cloud_liquid_fraction(SAM_moisture_type, tabs, an, bn);

    if (SAM_moisture_type == 1) {
        if (tabs >= tbgmax) {
            result.omn = one;
            result.delta_qi = qci;
            result.qci = zero;
            result.qcl += result.delta_qi;
            result.tabs -= fac_fus * result.delta_qi;
        } else if (tabs <= tbgmin) {
            result.omn = zero;
            result.delta_qc = qcl;
            result.qcl = zero;
            result.qci += result.delta_qc;
            result.tabs += fac_fus * result.delta_qc;
        } else {
            result.delta_qc = qcl - qn * result.omn;
            result.delta_qi = qci - qn * (one - result.omn);
            result.qcl = qn * result.omn;
            result.qci = qn * (one - result.omn);
            result.tabs += fac_fus * result.delta_qc;
        }
    } else {
        result.omn = one;
        result.delta_qc = qcl - qn;
        result.qcl = qn;
        result.qci = zero;
        result.tabs += fac_cond * result.delta_qc;
    }

    result.qn = result.qcl + result.qci;
    return result;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_mixed_qsat (const amrex::Real omn,
                            const amrex::Real qsatw,
                            const amrex::Real qsati) noexcept
{
    return omn * qsatw + (one - omn) * qsati;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_mixed_dqsat_dT (const amrex::Real omn,
                                const amrex::Real domn,
                                const amrex::Real qsatw,
                                const amrex::Real qsati,
                                const amrex::Real dqsatw,
                                const amrex::Real dqsati) noexcept
{
    return omn * dqsatw + (one - omn) * dqsati + domn * qsatw - domn * qsati;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_newton_residual (const amrex::Real tabs_new,
                                 const amrex::Real tabs_old,
                                 const amrex::Real lstar,
                                 const amrex::Real qv,
                                 const amrex::Real qsat) noexcept
{
    return -tabs_new + tabs_old + lstar * (qv - qsat);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_newton_residual_derivative (const amrex::Real lstar,
                                            const amrex::Real dlstar,
                                            const amrex::Real qv,
                                            const amrex::Real qsat,
                                            const amrex::Real dqsat) noexcept
{
    return -one + dlstar * (qv - qsat) - lstar * dqsat;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMCloudPhaseChange sam_apply_condensate_limiter (const amrex::Real qv,
                                                  const amrex::Real qsat,
                                                  const amrex::Real qc,
                                                  const amrex::Real qi,
                                                  const amrex::Real omn) noexcept
{
    SAMCloudPhaseChange result{};
    result.delta_qv = qv - qsat;
    result.delta_qc = amrex::max(-qc, result.delta_qv * omn);
    result.delta_qi = amrex::max(-qi, result.delta_qv * (one - omn));
    result.qv = qsat;
    result.qcl = qc + result.delta_qc;
    result.qci = qi + result.delta_qi;
    result.qn = result.qcl + result.qci;
    result.qt = result.qv + result.qn;
    return result;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMPrecipSources sam_autoconversion_rates (const amrex::Real dtn,
                                           const amrex::Real qcc,
                                           const amrex::Real qii,
                                           const amrex::Real coefice_k) noexcept
{
    SAMPrecipSources result{};
    // MUST MATCH: SAM rain autoconversion threshold.
    const amrex::Real auto_r = (qcc > qcw0) ? alphaelq : zero;
    const amrex::Real autos = (qii > qci0) ? betaelq * coefice_k : zero;
    result.dqca = dtn * auto_r * (qcc - qcw0);
    result.dqia = dtn * autos * (qii - qci0);
    return result;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMPrecipSources sam_accretion_rates (const amrex::Real dtn,
                                      const amrex::Real qcc,
                                      const amrex::Real qii,
                                      const amrex::Real qpr,
                                      const amrex::Real qps,
                                      const amrex::Real qpg,
                                      const amrex::Real powr1,
                                      const amrex::Real pows1,
                                      const amrex::Real powg1,
                                      const amrex::Real omp,
                                      const amrex::Real omg,
                                      const amrex::Real accrrc_k,
                                      const amrex::Real accrsc_k,
                                      const amrex::Real accrsi_k,
                                      const amrex::Real accrgc_k,
                                      const amrex::Real accrgi_k) noexcept
{
    SAMPrecipSources result{};
    amrex::Real accrcr = zero;
    amrex::Real accrcs = zero;
    amrex::Real accris = zero;
    amrex::Real accrcg = zero;
    amrex::Real accrig = zero;

    // MUST MATCH: SAM precipitation phase-fraction activation thresholds.
    if (omp > amrex::Real(0.001)) {
        accrcr = accrrc_k;
    }
    if (omp < amrex::Real(0.999) && omg < amrex::Real(0.999)) {
        accrcs = accrsc_k;
        accris = accrsi_k;
    }
    if (omp < amrex::Real(0.999) && omg > amrex::Real(0.001)) {
        accrcg = accrgc_k;
        accrig = accrgi_k;
    }

    result.dprc = dtn * accrcr * qcc * std::pow(qpr, powr1);
    result.dpsc = dtn * accrcs * qcc * std::pow(qps, pows1);
    result.dpgc = dtn * accrcg * qcc * std::pow(qpg, powg1);
    result.dpsi = dtn * accris * qii * std::pow(qps, pows1);
    result.dpgi = dtn * accrig * qii * std::pow(qpg, powg1);
    return result;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMPrecipSources sam_rescale_cloud_sinks (const amrex::Real qcl,
                                          const amrex::Real qci,
                                          const amrex::Real eps,
                                          SAMPrecipSources result) noexcept
{
    result.dqc = result.dqca + result.dprc + result.dpsc + result.dpgc;
    result.dqi = result.dqia + result.dpsi + result.dpgi;

    const amrex::Real scalec = amrex::min(qcl, result.dqc) / (result.dqc + eps);
    const amrex::Real scalei = amrex::min(qci, result.dqi) / (result.dqi + eps);

    result.dqca *= scalec;
    result.dprc *= scalec;
    result.dpsc *= scalec;
    result.dpgc *= scalec;
    result.dqia *= scalei;
    result.dpsi *= scalei;
    result.dpgi *= scalei;
    result.dqc = result.dqca + result.dprc + result.dpsc + result.dpgc;
    result.dqi = result.dqia + result.dpsi + result.dpgi;

    return result;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMPrecipSources sam_partition_autoconverted_precip (SAMPrecipSources result,
                                                     const amrex::Real omp,
                                                     const amrex::Real omg) noexcept
{
    result.dqpr = (result.dqca + result.dqia) * omp + result.dprc;
    result.dqps = (result.dqca + result.dqia) * (one - omp) * (one - omg) + result.dpsc + result.dpsi;
    result.dqpg = (result.dqca + result.dqia) * (one - omp) * omg + result.dpgc + result.dpgi;
    return result;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMPrecipSources sam_precip_evaporation_rates (const amrex::Real qpr,
                                               const amrex::Real qps,
                                               const amrex::Real qpg,
                                               const amrex::Real powr2,
                                               const amrex::Real pows2,
                                               const amrex::Real powg2,
                                               const amrex::Real evapr1_k,
                                               const amrex::Real evapr2_k,
                                               const amrex::Real evaps1_k,
                                               const amrex::Real evaps2_k,
                                               const amrex::Real evapg1_k,
                                               const amrex::Real evapg2_k) noexcept
{
    SAMPrecipSources result{};
    result.dqpr = evapr1_k * std::sqrt(qpr) + evapr2_k * std::pow(qpr, powr2);
    result.dqps = evaps1_k * std::sqrt(qps) + evaps2_k * std::pow(qps, pows2);
    result.dqpg = evapg1_k * std::sqrt(qpg) + evapg2_k * std::pow(qpg, powg2);
    return result;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMPrecipSources sam_apply_precip_evaporation_limiter (const amrex::Real qpr,
                                                       const amrex::Real qps,
                                                       const amrex::Real qpg,
                                                       SAMPrecipSources result) noexcept
{
    result.dqpr = amrex::min(qpr, result.dqpr);
    result.dqps = amrex::min(qps, result.dqps);
    result.dqpg = amrex::min(qpg, result.dqpg);
    result.dqp = result.dqpr + result.dqps + result.dqpg;
    return result;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMCoefficientRow sam_compute_coefficient_row (const amrex::Real rho,
                                               const amrex::Real tabs,
                                               const amrex::Real gamr1,
                                               const amrex::Real gamr2,
                                               const amrex::Real gams1,
                                               const amrex::Real gams2,
                                               const amrex::Real gamg1,
                                               const amrex::Real gamg2) noexcept
{
    // Preserve the SAM baseline formula ordering in these closure rows.
    SAMCoefficientRow row{};
    amrex::Real prefactor;
    amrex::Real pratio = std::sqrt(amrex::Real(1.29) / rho);
    amrex::Real estw = amrex::Real(100.0) * erf_esatw(tabs);
    amrex::Real esti = amrex::Real(100.0) * erf_esati(tabs);

    amrex::Real coef1 = fourth * PI * nzeros * a_snow * gams1 * pratio
        / std::pow((PI * rhos * nzeros / rho), ((three + b_snow) / amrex::Real(4.0)));
    amrex::Real coef2 = std::exp(amrex::Real(0.025) * (tabs - amrex::Real(273.15)));
    row.accrsi = coef1 * coef2 * esicoef;
    row.accrsc = coef1 * esccoef;
    row.coefice = coef2;

    coef1 = (lsub / (tabs * R_v) - one) * lsub / (therco * tabs);
    coef2 = R_v * R_d / (diffelq * esti);
    prefactor = two * PI * nzeros / (rho * (coef1 + coef2));
    prefactor *= (two / PI);
    row.evaps1 = prefactor * amrex::Real(0.65) * std::sqrt(rho / (PI * rhos * nzeros));
    row.evaps2 = prefactor * amrex::Real(0.44) * std::sqrt(a_snow * rho / muelq) * gams2
        * std::sqrt(pratio) * std::pow(rho / (PI * rhos * nzeros), ((amrex::Real(5.0) + b_snow) / amrex::Real(8.0)));

    coef1 = fourth * PI * nzerog * a_grau * gamg1 * pratio
        / std::pow((PI * rhog * nzerog / rho), ((three + b_grau) / amrex::Real(4.0)));
    coef2 = std::exp(amrex::Real(0.025) * (tabs - amrex::Real(273.15)));
    row.accrgi = coef1 * coef2 * egicoef;
    row.accrgc = coef1 * egccoef;

    coef1 = (lsub / (tabs * R_v) - one) * lsub / (therco * tabs);
    coef2 = R_v * R_d / (diffelq * esti);
    prefactor = two * PI * nzerog / (rho * (coef1 + coef2));
    row.evapg1 = prefactor * amrex::Real(0.78) * std::sqrt(rho / (PI * rhog * nzerog));
    row.evapg2 = prefactor * amrex::Real(0.31) * std::sqrt(a_grau * rho / muelq) * gamg2
        * std::sqrt(pratio) * std::pow(rho / (PI * rhog * nzerog), ((amrex::Real(5.0) + b_grau) / amrex::Real(8.0)));

    row.accrrc = fourth * PI * nzeror * a_rain * gamr1 * pratio
        / std::pow((PI * rhor * nzeror / rho), ((three + b_rain) / amrex::Real(4.0))) * erccoef;

    coef1 = (lcond / (tabs * R_v) - one) * lcond / (therco * tabs);
    coef2 = R_v * R_d / (diffelq * estw);
    prefactor = two * PI * nzeror / (rho * (coef1 + coef2));
    row.evapr1 = prefactor * amrex::Real(0.78) * std::sqrt(rho / (PI * rhor * nzeror));
    row.evapr2 = prefactor * amrex::Real(0.31) * std::sqrt(a_rain * rho / muelq) * gamr2
        * std::sqrt(pratio) * std::pow(rho / (PI * rhor * nzeror), ((amrex::Real(5.0) + b_rain) / amrex::Real(8.0)));
    return row;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_cloud_ice_terminal_velocity (const amrex::Real qci_avg) noexcept
{
    // MUST MATCH: SAM cloud-ice terminal-velocity coefficient and exponent.
    return amrex::min(amrex::Real(0.4),
                      amrex::Real(8.66) * std::pow((amrex::max(amrex::Real(0.0), qci_avg) + amrex::Real(1.e-10)), amrex::Real(0.24)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMFaceState sam_face_average_state (const int k,
                                     const int k_lo,
                                     const int k_hi,
                                     const amrex::Real rho_km1,
                                     const amrex::Real rho_k,
                                     const amrex::Real tabs_km1,
                                     const amrex::Real tabs_k,
                                     const amrex::Real qci_km1,
                                     const amrex::Real qci_k,
                                     const amrex::Real qp_km1,
                                     const amrex::Real qp_k) noexcept
{
    SAMFaceState result{};
    if (k == k_lo) {
        result.rho_avg = rho_k;
        result.tabs_avg = tabs_k;
        result.qci_avg = qci_k;
        result.qp_avg = qp_k;
    } else if (k == k_hi + 1) {
        result.rho_avg = rho_km1;
        result.tabs_avg = tabs_km1;
        result.qci_avg = qci_km1;
        result.qp_avg = qp_km1;
    } else {
        result.rho_avg = myhalf * (rho_km1 + rho_k);
        result.tabs_avg = myhalf * (tabs_km1 + tabs_k);
        result.qci_avg = myhalf * (qci_km1 + qci_k);
        result.qp_avg = myhalf * (qp_km1 + qp_k);
    }
    return result;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int sam_substep_count_from_reduced_flux (const amrex::Real reduced_flux,
                                         const amrex::Real dt,
                                         const amrex::Real dz) noexcept
{
    // Substeps follow the legacy reduced-flux rule. Property tests will audit it.
    // MUST MATCH: SAM::CFL_MAX.
    return static_cast<int>(std::ceil(reduced_flux * (dt / dz) / myhalf));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMPrecipFaceState sam_precip_face_state (
    const int SAM_moisture_type,
    const amrex::Array4<const amrex::Real>& rho_array,
    const amrex::Array4<const amrex::Real>& tabs_array,
    const amrex::Array4<const amrex::Real>& qp_array,
    const int i,
    const int j,
    const int k,
    const int k_lo,
    const int k_hi) noexcept
{
    // MUST MATCH: SAM PrecipFall bottom/interior/top face rules.
    SAMPrecipFaceState result{};
    if (k == k_lo) {
        result.rho_avg = rho_array(i,j,k);
        result.tabs_avg = tabs_array(i,j,k);
        result.qp_avg = qp_array(i,j,k);
    } else if (k == k_hi + 1) {
        result.rho_avg = rho_array(i,j,k-1);
        result.tabs_avg = tabs_array(i,j,k-1);
        result.qp_avg = qp_array(i,j,k-1);
    } else {
        result.rho_avg = myhalf * (rho_array(i,j,k-1) + rho_array(i,j,k));
        result.tabs_avg = myhalf * (tabs_array(i,j,k-1) + tabs_array(i,j,k));
        result.qp_avg = myhalf * (qp_array(i,j,k-1) + qp_array(i,j,k));
    }
    result.omp = sam_precip_rain_fraction(SAM_moisture_type, result.tabs_avg);
    result.omg = sam_graupel_fraction(SAM_moisture_type, result.tabs_avg);
    result.qrr = result.omp * result.qp_avg;
    result.qss = (one - result.omp) * (one - result.omg) * result.qp_avg;
    result.qgg = (one - result.omp) * result.omg * result.qp_avg;
    return result;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_precip_flux_from_face_state (const SAMPrecipFaceState& face_state,
                                             const amrex::Real vrain,
                                             const amrex::Real vsnow,
                                             const amrex::Real vgrau) noexcept
{
    if (face_state.qp_avg <= qp_threshold) {
        return zero;
    }
    return face_state.omp * vrain * std::pow(face_state.rho_avg * face_state.qrr, one + crain)
        + (one - face_state.omp)
            * ((one - face_state.omg) * vsnow * std::pow(face_state.rho_avg * face_state.qss, one + csnow)
               + face_state.omg * vgrau * std::pow(face_state.rho_avg * face_state.qgg, one + cgrau));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_precip_flux_density_corrected (const amrex::Real precip_flux,
                                               const amrex::Real rho_0,
                                               const amrex::Real rho_avg) noexcept
{
    return precip_flux * std::sqrt(rho_0 / rho_avg);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_sedimentation_tendency (const amrex::Real fz_hi,
                                        const amrex::Real fz_lo,
                                        const amrex::Real rho,
                                        const amrex::Real dJinv,
                                        const amrex::Real coef) noexcept
{
    return dJinv * (one / rho) * (fz_hi - fz_lo) * coef;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMSurfaceAccumulation sam_surface_accumulation (const amrex::Real rho_avg,
                                                 const amrex::Real qp_avg,
                                                 const amrex::Real omp,
                                                 const amrex::Real omg,
                                                 const amrex::Real vrain,
                                                 const amrex::Real vsnow,
                                                 const amrex::Real vgrau,
                                                 const amrex::Real dtn) noexcept
{
    SAMSurfaceAccumulation result{};
    // Surface accumulation preserves the SAM baseline formula. Property tests will audit it.
    result.rain = rho_avg * (omp * qp_avg) * vrain * dtn / rhor * amrex::Real(1000.0);
    result.snow = rho_avg * (one - omp) * (one - omg) * qp_avg * vsnow * dtn / rhos * amrex::Real(1000.0);
    result.graupel = rho_avg * (one - omp) * omg * qp_avg * vgrau * dtn / rhog * amrex::Real(1000.0);
    return result;
}

#endif