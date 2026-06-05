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

using SAMCellState = SAMPrimitiveCell;

struct SAMPrecipConfig {
    int sam_moisture_type{1};
    bool enable_precip{true};
    amrex::Real dtn{zero};
    amrex::Real rdOcp{zero};
    amrex::Real fac_cond{zero};
    amrex::Real fac_fus{zero};
    amrex::Real fac_sub{zero};
    amrex::Real eps{zero};
    amrex::Real powr1{zero};
    amrex::Real pows1{zero};
    amrex::Real powg1{zero};
    amrex::Real powr2{zero};
    amrex::Real pows2{zero};
    amrex::Real powg2{zero};
};

struct SAMPrecipCellDiagnostics {
    SAMPrecipSources autoconversion{};
    SAMPrecipSources accretion{};
    SAMPrecipSources limited_sources{};
    SAMPrecipSources partitioned_sources{};
    SAMPrecipSources evaporation{};
    amrex::Real qsat{zero};
    amrex::Real omn{one};
    amrex::Real omp{one};
    amrex::Real omg{zero};
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

struct SAMPrecipComponentFaceState {
    amrex::Real rho_avg{zero};
    amrex::Real tabs_avg{zero};
    amrex::Real qpr_avg{zero};
    amrex::Real qps_avg{zero};
    amrex::Real qpg_avg{zero};
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

struct SAMPrecipFluxComponents {
    amrex::Real rain{zero};
    amrex::Real snow{zero};
    amrex::Real graupel{zero};
};

// SAM local source updates are treated as constant-pressure microphysics
// updates. Copy-in reconstructs tabs and pres_mbar from the conserved state;
// source helpers hold pres_mbar fixed while latent heating/cooling changes
// tabs; theta is then refreshed from T and the held pressure. Pressure is
// stored in mbar on this path, EOS/theta helpers expect Pa, and the saturation
// helpers below consume mbar where SAM currently calls them that way.
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
    // Do not replace this with getThgivenRandT(rho,T,qv): that would enforce a
    // fixed-density EOS projection rather than the SAM source-step
    // thermodynamic contract.
    return getThgivenTandP(tabs, sam_mbar_to_pa(pres_mbar), rdOcp);
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
    // Copy-in reconstructs the SAM primitive cell state from the conserved ERF
    // state, clips negative moisture defensively, and diagnoses tabs and
    // pres_mbar from the incoming rho*theta/rho*q state. That diagnosed mbar
    // pressure is the held pressure for the following source updates.
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
    // Copy-out writes rho*theta and clipped moisture species back to the
    // conserved state. Pressure is not stored here; it is diagnosed again on
    // the next copy-in.
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
    // Empirical cloud phase-fraction closure for omega_n. No-ice mode forces
    // liquid fraction to one; the warm branch is all liquid, the cold branch
    // is all ice, and the mixed branch varies linearly with temperature.
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
    // Empirical precipitating rain fraction omega_p. No-ice mode forces rain
    // fraction to one, and the result is clipped to [0, 1].
    if (SAM_moisture_type == 2) {
        return one;
    }
    return amrex::max(amrex::Real(0.0), amrex::min(amrex::Real(1.0), (tabs - tprmin) * a_pr));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_graupel_fraction (const int SAM_moisture_type,
                                  const amrex::Real tabs) noexcept
{
    // Empirical graupel fraction omega_g. No-ice mode forces graupel fraction
    // to zero, and the result is clipped to [0, 1].
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
    // Repartition cloud condensate between liquid and ice using omega_n while
    // preserving total cloud condensate qn. The temperature increment accounts
    // for fusion heating/cooling under the held-pressure source-step
    // convention; no-ice mode collapses condensate into liquid.
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
    // Mixed liquid/ice saturation relation used by the SAM cloud adjustment.
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
    // Temperature derivative of mixed saturation, including the mixed-branch
    // temperature derivative of the liquid/ice phase fraction.
    return omn * dqsatw + (one - omn) * dqsati + domn * qsatw - domn * qsati;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_newton_residual (const amrex::Real tabs_new,
                                 const amrex::Real tabs_old,
                                 const amrex::Real lstar,
                                 const amrex::Real qv,
                                 const amrex::Real qsat) noexcept
{
    // Residual for the local constant-pressure saturation adjustment solve.
    return -tabs_new + tabs_old + lstar * (qv - qsat);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_newton_residual_derivative (const amrex::Real lstar,
                                            const amrex::Real dlstar,
                                            const amrex::Real qv,
                                            const amrex::Real qsat,
                                            const amrex::Real dqsat) noexcept
{
    // Derivative for the local constant-pressure saturation adjustment solve.
    return -one + dlstar * (qv - qsat) - lstar * dqsat;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMCloudPhaseChange sam_apply_condensate_limiter (const amrex::Real qv,
                                                  const amrex::Real qsat,
                                                  const amrex::Real qc,
                                                  const amrex::Real qi,
                                                  const amrex::Real omn) noexcept
{
    // Prevents the adjustment from evaporating or depositing more condensate
    // than is locally available.
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
    // Thresholded cloud-liquid and cloud-ice autoconversion increments over
    // the microphysics time step dtn. The strict qcw0/qci0 thresholds are
    // intentional: equality leaves the source inactive.
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

    // Accretion increments over dtn. Liquid-origin sinks feed rain, snow, and
    // graupel; ice-origin sinks feed snow and graupel. These strict empirical
    // selector thresholds are implementation rules, not fundamental laws:
    // rain requires omp > 0.001, snow requires omp < 0.999 and omg < 0.999,
    // and graupel requires omp < 0.999 and omg > 0.001. Equality is inactive.
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
    // Proportional donor limiter. If requested liquid or ice sinks exceed
    // available qcl/qci, all competing sinks from the same donor phase share
    // one scale factor. This preserves source ratios while preventing negative
    // qcl or qci. The eps denominator intentionally leaves only an epsilon-
    // scale undershoot.
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
    // Partitions limited autoconversion sources into rain/snow/graupel using
    // the precipitating phase fractions and adds direct accretion
    // contributions. dqpr + dqps + dqpg equals the limited cloud sink total.
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
    // Positive-domain evaporation/sublimation rate form for precipitating
    // species: C1*sqrt(q_m) + C2*q_m^p. The derivative is singular at q_m = 0,
    // so derivative checks only make sense away from zero.
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
    // Limits rain evaporation and snow/graupel sublimation by available
    // qpr/qps/qpg so the source update cannot make precipitating species
    // negative.
    result.dqpr = amrex::min(qpr, result.dqpr);
    result.dqps = amrex::min(qps, result.dqps);
    result.dqpg = amrex::min(qpg, result.dqpg);
    result.dqp = result.dqpr + result.dqps + result.dqpg;
    return result;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMCellState sam_precip_cell_update (SAMCellState state,
                                     const SAMCoefficientRow& coeffs,
                                     const SAMPrecipConfig& config,
                                     SAMPrecipCellDiagnostics* diagnostics = nullptr) noexcept
{
    // Local SAM precipitation source update. This is a cell-local,
    // held-pressure microphysics update; sedimentation and surface
    // accumulation are handled in PrecipFall. Before sedimentation, the source
    // terms exchange water among qv/qcl/qci/qpr/qps/qpg, conserve total water,
    // and conserve the constant-pressure latent proxy
    // tabs + fac_cond*qv - fac_fus*(qci + qps + qpg)
    // up to limiter and roundoff effects.
    if (diagnostics != nullptr) {
        *diagnostics = {};
    }

    state.qn = state.qcl + state.qci;
    state.qt = state.qv + state.qn;
    state.qp = state.qpr + state.qps + state.qpg;

    if (!config.enable_precip || state.qn + state.qp <= zero) {
        return state;
    }

    amrex::Real omn;
    amrex::Real omp;
    amrex::Real omg;
    if (config.sam_moisture_type == 2) {
        omn = one;
        omp = one;
        omg = zero;
    } else {
        omn = sam_cloud_liquid_fraction(config.sam_moisture_type,
                                        state.tabs, a_bg, tbgmin * a_bg);
        omp = sam_precip_rain_fraction(config.sam_moisture_type,
                                       state.tabs);
        omg = sam_graupel_fraction(config.sam_moisture_type,
                                   state.tabs);
    }

    if (diagnostics != nullptr) {
        diagnostics->omn = omn;
        diagnostics->omp = omp;
        diagnostics->omg = omg;
    }

    if (state.qn > zero) {
        SAMPrecipSources source_terms =
            sam_autoconversion_rates(config.dtn, state.qcl, state.qci, coeffs.coefice);
        const SAMPrecipSources accretion_terms =
            sam_accretion_rates(config.dtn, state.qcl, state.qci, state.qpr, state.qps, state.qpg,
                                config.powr1, config.pows1, config.powg1,
                                omp, omg,
                                coeffs.accrrc, coeffs.accrsc, coeffs.accrsi,
                                coeffs.accrgc, coeffs.accrgi);
        if (diagnostics != nullptr) {
            diagnostics->autoconversion = source_terms;
            diagnostics->accretion = accretion_terms;
        }

        source_terms.dprc = accretion_terms.dprc;
        source_terms.dpsc = accretion_terms.dpsc;
        source_terms.dpgc = accretion_terms.dpgc;
        source_terms.dpsi = accretion_terms.dpsi;
        source_terms.dpgi = accretion_terms.dpgi;
        source_terms = sam_rescale_cloud_sinks(state.qcl, state.qci,
                                               config.eps, source_terms);
        if (diagnostics != nullptr) {
            diagnostics->limited_sources = source_terms;
        }

        source_terms = sam_partition_autoconverted_precip(source_terms, omp, omg);
        if (diagnostics != nullptr) {
            diagnostics->partitioned_sources = source_terms;
        }

        const amrex::Real dqca = source_terms.dqca;
        const amrex::Real dqia = source_terms.dqia;

        state.qcl -= source_terms.dqc;
        state.qci -= source_terms.dqi;
        state.qpr += source_terms.dqpr;
        state.qps += source_terms.dqps;
        state.qpg += source_terms.dqpg;

        state.qn = state.qcl + state.qci;
        state.qt = state.qv + state.qn;
        state.qp = state.qpr + state.qps + state.qpg;

        // Cloud-liquid accretion onto snow and graupel increases frozen mass
        // and must contribute the same fusion-heating term to the local latent
        // proxy as the autoconverted frozen branch.
        state.tabs += config.fac_fus * (dqca * (one - omp) - dqia * omp
                          + source_terms.dpsc + source_terms.dpgc);
        state.theta = sam_theta_from_stored_mbar_converted_to_pa(state.tabs,
                                                                 state.pres_mbar,
                                                                 config.rdOcp);
    }

    amrex::Real qsatw;
    amrex::Real qsati;
    erf_qsatw(state.tabs, state.pres_mbar, qsatw);
    erf_qsati(state.tabs, state.pres_mbar, qsati);
    const amrex::Real qsat = sam_mixed_qsat(omn, qsatw, qsati);
    if (diagnostics != nullptr) {
        diagnostics->qsat = qsat;
    }

    if (state.qp > zero && state.qv < qsat) {
        SAMPrecipSources evaporation_terms =
            sam_precip_evaporation_rates(state.qpr, state.qps, state.qpg,
                                         config.powr2, config.pows2, config.powg2,
                                         coeffs.evapr1, coeffs.evapr2,
                                         coeffs.evaps1, coeffs.evaps2,
                                         coeffs.evapg1, coeffs.evapg2);

        const amrex::Real supersat_factor = config.dtn * (one - state.qv / qsat);
        evaporation_terms.dqpr *= supersat_factor;
        evaporation_terms.dqps *= supersat_factor;
        evaporation_terms.dqpg *= supersat_factor;
        evaporation_terms = sam_apply_precip_evaporation_limiter(state.qpr,
                                                                 state.qps,
                                                                 state.qpg,
                                                                 evaporation_terms);
        if (diagnostics != nullptr) {
            diagnostics->evaporation = evaporation_terms;
        }

        state.qv += evaporation_terms.dqp;
        state.qpr -= evaporation_terms.dqpr;
        state.qps -= evaporation_terms.dqps;
        state.qpg -= evaporation_terms.dqpg;

        state.qt = state.qv + state.qn;
        state.qp = state.qpr + state.qps + state.qpg;

        // Rain evaporation cools with Lv/cp; snow and graupel sublimation cool
        // with Ls/cp, preserving the same local latent proxy.
        state.tabs -= config.fac_cond * evaporation_terms.dqpr
                    + config.fac_sub * (evaporation_terms.dqps + evaporation_terms.dqpg);
        state.theta = sam_theta_from_stored_mbar_converted_to_pa(state.tabs,
                                                                 state.pres_mbar,
                                                                 config.rdOcp);
    }

    return state;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
// Closure coefficients for the SAM source parameterizations at one vertical
// level. The formula ordering is kept close to the SAM baseline to avoid
// unintended numerical drift.
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
int sam_precip_face_donor_k (const int face_k,
                             const int k_lo,
                             const int k_hi) noexcept
{
    if (face_k <= k_lo) {
        return k_lo;
    }
    if (face_k >= k_hi + 1) {
        return k_hi;
    }
    return face_k;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int sam_substep_count_from_reduced_flux (const amrex::Real reduced_flux,
                                         const amrex::Real dt,
                                         const amrex::Real dz) noexcept
{
    // Substeps follow the legacy reduced-flux rule based on the maximum
    // density-corrected precip flux, not directly on a terminal fall speed.
    // Property tests characterize this behavior before any intended change.
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
// Component PrecipFall face state. Uses the same bottom/interior/top face
// selection rules, but reads qpr/qps/qpg directly so each precipitating
// species can be sedimented and limited conservatively.
SAMPrecipComponentFaceState sam_precip_component_face_state (
    const amrex::Array4<const amrex::Real>& rho_array,
    const amrex::Array4<const amrex::Real>& tabs_array,
    const amrex::Array4<const amrex::Real>& qpr_array,
    const amrex::Array4<const amrex::Real>& qps_array,
    const amrex::Array4<const amrex::Real>& qpg_array,
    const int i,
    const int j,
    const int k,
    const int k_lo,
    const int k_hi) noexcept
{
    // MUST MATCH: SAM PrecipFall bottom/interior/top face rules.
    SAMPrecipComponentFaceState result{};
    if (k == k_lo) {
        result.rho_avg = rho_array(i,j,k);
        result.tabs_avg = tabs_array(i,j,k);
        result.qpr_avg = qpr_array(i,j,k);
        result.qps_avg = qps_array(i,j,k);
        result.qpg_avg = qpg_array(i,j,k);
    } else if (k == k_hi + 1) {
        result.rho_avg = rho_array(i,j,k-1);
        result.tabs_avg = tabs_array(i,j,k-1);
        result.qpr_avg = qpr_array(i,j,k-1);
        result.qps_avg = qps_array(i,j,k-1);
        result.qpg_avg = qpg_array(i,j,k-1);
    } else {
        result.rho_avg = myhalf * (rho_array(i,j,k-1) + rho_array(i,j,k));
        result.tabs_avg = myhalf * (tabs_array(i,j,k-1) + tabs_array(i,j,k));
        result.qpr_avg = myhalf * (qpr_array(i,j,k-1) + qpr_array(i,j,k));
        result.qps_avg = myhalf * (qps_array(i,j,k-1) + qps_array(i,j,k));
        result.qpg_avg = myhalf * (qpg_array(i,j,k-1) + qpg_array(i,j,k));
    }
    return result;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
// Legacy total-qp precip flux helper. It partitions a total qp into rain/
// snow/graupel by temperature phase fractions and applies the qp_threshold
// branch. The conservative production PrecipFall path uses the component
// helper below.
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
// Legacy total-qp companion returning the partitioned rain/snow/graupel face
// fluxes. The conservative production PrecipFall path uses the component
// helper below instead of repartitioning total qp.
SAMPrecipFluxComponents sam_precip_flux_components_from_face_state (const SAMPrecipFaceState& face_state,
                                                                    const amrex::Real vrain,
                                                                    const amrex::Real vsnow,
                                                                    const amrex::Real vgrau) noexcept
{
    SAMPrecipFluxComponents result{};
    if (face_state.qp_avg <= qp_threshold) {
        return result;
    }

    result.rain = face_state.omp * vrain * std::pow(face_state.rho_avg * face_state.qrr, one + crain);
    result.snow = (one - face_state.omp) * (one - face_state.omg)
        * vsnow * std::pow(face_state.rho_avg * face_state.qss, one + csnow);
    result.graupel = (one - face_state.omp) * face_state.omg
        * vgrau * std::pow(face_state.rho_avg * face_state.qgg, one + cgrau);
    return result;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
// Component PrecipFall fluxes. Unlike the legacy total-qp helper, these
// fluxes are zero only when the corresponding component mass is nonpositive;
// tiny positive qpr/qps/qpg produce tiny positive fluxes and are donor-limited
// later.
SAMPrecipFluxComponents sam_precip_component_fluxes_from_face_state (
    const SAMPrecipComponentFaceState& face_state,
    const amrex::Real vrain,
    const amrex::Real vsnow,
    const amrex::Real vgrau) noexcept
{
    SAMPrecipFluxComponents result{};
    if (face_state.qpr_avg > zero) {
        result.rain = vrain * std::pow(face_state.rho_avg * face_state.qpr_avg, one + crain);
    }
    if (face_state.qps_avg > zero) {
        result.snow = vsnow * std::pow(face_state.rho_avg * face_state.qps_avg, one + csnow);
    }
    if (face_state.qpg_avg > zero) {
        result.graupel = vgrau * std::pow(face_state.rho_avg * face_state.qpg_avg, one + cgrau);
    }
    return result;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_precip_flux_density_corrected (const amrex::Real precip_flux,
                                               const amrex::Real rho_0,
                                               const amrex::Real rho_avg) noexcept
{
    return precip_flux * std::sqrt(rho_0 / rho_avg);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMPrecipFluxComponents sam_precip_flux_components_density_corrected (
    const SAMPrecipFluxComponents& precip_fluxes,
    const amrex::Real rho_0,
    const amrex::Real rho_avg) noexcept
{
    const amrex::Real density_correction = std::sqrt(rho_0 / rho_avg);
    return {
        precip_fluxes.rain * density_correction,
        precip_fluxes.snow * density_correction,
        precip_fluxes.graupel * density_correction};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
// Donor-cap limiter for conservative component sedimentation. The limited flux
// cannot export more component mass over one substep than the donor cell
// contains, including detJ for terrain-following cell mass.
amrex::Real sam_limit_precip_component_flux (
    const amrex::Real raw_flux,
    const amrex::Real rho_donor,
    const amrex::Real q_donor,
    const amrex::Real detJ_donor,
    const amrex::Real coef) noexcept
{
    const amrex::Real available_flux =
        rho_donor * amrex::max(amrex::Real(0.0), q_donor) * detJ_donor / coef;
    return amrex::min(raw_flux, available_flux);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
// Converts face-flux divergence into a mixing-ratio tendency. With the same
// limited face fluxes used for accumulation, column component budgets
// telescope.
amrex::Real sam_sedimentation_tendency (const amrex::Real fz_hi,
                                        const amrex::Real fz_lo,
                                        const amrex::Real rho,
                                        const amrex::Real dJinv,
                                        const amrex::Real coef) noexcept
{
    return dJinv * (one / rho) * (fz_hi - fz_lo) * coef;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
// Converts the already-limited bottom component fluxes to accumulated surface
// depths. These are the same fluxes used in the in-column sedimentation
// update.
SAMSurfaceAccumulation sam_surface_accumulation_from_component_fluxes (
    const SAMPrecipFluxComponents& component_fluxes,
    const amrex::Real dtn) noexcept
{
    SAMSurfaceAccumulation result{};
    result.rain = component_fluxes.rain * dtn / rhor * amrex::Real(1000.0);
    result.snow = component_fluxes.snow * dtn / rhos * amrex::Real(1000.0);
    result.graupel = component_fluxes.graupel * dtn / rhog * amrex::Real(1000.0);
    return result;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMSurfaceAccumulation sam_surface_accumulation (const SAMPrecipFaceState& face_state,
                                                 const amrex::Real rho_0,
                                                 const amrex::Real vrain,
                                                 const amrex::Real vsnow,
                                                 const amrex::Real vgrau,
                                                 const amrex::Real dtn) noexcept
{
    const SAMPrecipFluxComponents face_fluxes =
        sam_precip_flux_components_from_face_state(face_state, vrain, vsnow, vgrau);

    return sam_surface_accumulation_from_component_fluxes(
        sam_precip_flux_components_density_corrected(face_fluxes, rho_0, face_state.rho_avg),
        dtn);
}

#endif