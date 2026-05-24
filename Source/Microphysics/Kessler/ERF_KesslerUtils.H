#ifndef ERF_KESSLER_UTILS_H_
#define ERF_KESSLER_UTILS_H_

#include <algorithm>
#include <cmath>
#include <limits>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <ERF_Constants.H>
#include <ERF_MicrophysicsUtils.H>

struct KesslerSourceTerms {
    amrex::Real dq_vapor_to_cloud;
    amrex::Real dq_cloud_to_vapor;
    amrex::Real dq_cloud_to_rain;
    amrex::Real dq_rain_to_vapor;
};

struct KesslerSaturationAdjustment {
    amrex::Real dq_vapor_to_cloud;
    amrex::Real dq_cloud_to_vapor;
};

struct KesslerFaceState {
    amrex::Real rho;
    amrex::Real qp;
};

static constexpr amrex::Real kessler_sedimentation_cfl_max = myhalf;
static constexpr amrex::Real kessler_water_depth_m_per_kg_m2 = amrex::Real(1.0) / rhor;
static constexpr amrex::Real kessler_millimeters_per_meter = amrex::Real(1000.0);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int kessler_face_donor_k (const int k,
                          const int k_hi) noexcept
{
    return (k == k_hi + 1) ? k - 1 : k;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real kessler_terminal_velocity (const amrex::Real rho,
                                       const amrex::Real qp) noexcept
{
    // MUST MATCH: current AdvanceKessler terminal-velocity coefficients.
    return amrex::Real(36.34)
        * std::pow(rho * amrex::Real(0.001) * qp, amrex::Real(0.1346))
        * std::pow(rho / amrex::Real(1.16), -myhalf);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real kessler_precip_flux (const amrex::Real rho,
                                 const amrex::Real terminal_velocity,
                                 const amrex::Real qp) noexcept
{
    return rho * terminal_velocity * qp;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int kessler_num_sedimentation_substeps (const amrex::Real current_reduced_value,
                                        const amrex::Real dt,
                                        const amrex::Real dzmin) noexcept
{
    // MUST MATCH: current AdvanceKessler CFL_MAX precipitation substep formula.
    return static_cast<int>(std::ceil((current_reduced_value + std::numeric_limits<amrex::Real>::epsilon())
                                      * (dt / dzmin) / kessler_sedimentation_cfl_max));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real kessler_rain_accumulation_increment (const amrex::Real precip_mass_per_area) noexcept
{
    return precip_mass_per_area
        * kessler_water_depth_m_per_kg_m2
        * kessler_millimeters_per_meter;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool kessler_is_small_sedimentation_value (const amrex::Real value) noexcept
{
    // MUST MATCH: current AdvanceKessler sedimentation zero threshold literal.
    return std::fabs(value) < 1e-14;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
KesslerSaturationAdjustment kessler_saturation_adjustment (const amrex::Real qv,
                                                           const amrex::Real qc,
                                                           const amrex::Real qsat,
                                                           const amrex::Real dtqsat,
                                                           const bool do_cond,
                                                           const amrex::Real latent_over_cp) noexcept
{
    KesslerSaturationAdjustment source_terms{amrex::Real(0), amrex::Real(0)};

    const amrex::Real fac = latent_over_cp * dtqsat;

    // MUST MATCH: current AdvanceKessler condensation branch behavior.
    if ((qv > qsat) && do_cond) {
        source_terms.dq_vapor_to_cloud = std::min(qv, (qv - qsat) / (amrex::Real(1) + fac));
    }

    // MUST MATCH: current AdvanceKessler cloud-evaporation branch behavior.
    if ((qv < qsat) && (qc > amrex::Real(0)) && do_cond) {
        source_terms.dq_cloud_to_vapor = std::min(qc, (qsat - qv) / (amrex::Real(1) + fac));
    }

    return source_terms;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
KesslerFaceState kessler_face_state (const int k,
                                     const int k_hi,
                                     const amrex::Real rho_km1,
                                     const amrex::Real rho_k,
                                     const amrex::Real qp_km1,
                                     const amrex::Real qp_k) noexcept
{
    KesslerFaceState face_state{amrex::Real(0), amrex::Real(0)};

    // Downward sedimentation uses the upper-cell donor state at interior faces.
    // The top boundary face uses the top cell below the boundary.
    if (k == k_hi + 1) {
        face_state.rho = rho_km1;
        face_state.qp = qp_km1;
    } else {
        face_state.rho = rho_k;
        face_state.qp = qp_k;
    }

    // MUST MATCH: current AdvanceKessler nonnegative face qp clipping.
    face_state.qp = std::max(amrex::Real(0), face_state.qp);
    return face_state;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real kessler_sedimentation_tendency (const amrex::Real fz_hi,
                                            const amrex::Real fz_lo,
                                            const amrex::Real rho,
                                            const amrex::Real dJinv,
                                            const amrex::Real coef) noexcept
{
    // MUST MATCH: current AdvanceKessler sedimentation tendency arithmetic order.
    return dJinv * (amrex::Real(1) / rho) * (fz_hi - fz_lo) * coef;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
KesslerSourceTerms kessler_warm_rain_sources (const amrex::Real qv,
                                              const amrex::Real qc,
                                              const amrex::Real qp,
                                              const amrex::Real rho,
                                              const amrex::Real pressure_current_units,
                                              const amrex::Real qsat,
                                              const amrex::Real dtqsat,
                                              const amrex::Real dt,
                                              const bool do_cond,
                                              const amrex::Real latent_over_cp) noexcept
{
    // Units must match current AdvanceKessler: q* are mass fractions, rho is in
    // the production density units currently carried by Kessler, pressure is in
    // the current production mbar units, and dt is in seconds.
    const KesslerSaturationAdjustment saturation_adjustment =
        kessler_saturation_adjustment(qv, qc, qsat, dtqsat, do_cond, latent_over_cp);
    KesslerSourceTerms source_terms{
        saturation_adjustment.dq_vapor_to_cloud,
        saturation_adjustment.dq_cloud_to_vapor,
        amrex::Real(0),
        amrex::Real(0)};

    const amrex::Real capacity = amrex::max(
        amrex::Real(0), (qsat - qv) / (amrex::Real(1) + latent_over_cp * dtqsat));
    const amrex::Real remaining_capacity = amrex::max(
        amrex::Real(0), capacity - saturation_adjustment.dq_cloud_to_vapor);
    // Saturation adjustment returns at most one nonzero phase-change direction:
    // condensation reduces vapor in supersaturation, while cloud evaporation
    // increases vapor in subsaturation. Rain evaporation uses the vapor state
    // after that cloud adjustment.
    const amrex::Real qv_eff = qv
        - saturation_adjustment.dq_vapor_to_cloud
        + saturation_adjustment.dq_cloud_to_vapor;

    if ((qp > amrex::Real(0)) && (qsat > amrex::Real(0))
        && (pressure_current_units > amrex::Real(0))
        && (qv_eff < qsat) && (remaining_capacity > amrex::Real(0))) {
        // MUST MATCH: current AdvanceKessler rain-evaporation coefficients and exponents.
        const amrex::Real coeff = amrex::Real(1.6)
            + amrex::Real(124.9) * std::pow(amrex::Real(0.001) * rho * qp, amrex::Real(0.2046));
        const amrex::Real raw_rain_to_vapor = amrex::Real(1) / (amrex::Real(0.001) * rho)
            * (amrex::Real(1) - qv_eff / qsat)
            * coeff
            * std::pow(amrex::Real(0.001) * rho * qp, amrex::Real(0.525))
            / (amrex::Real(5.4e5) + amrex::Real(2.55e6) / (pressure_current_units * qsat))
            * dt;
        source_terms.dq_rain_to_vapor = amrex::min(qp, amrex::min(raw_rain_to_vapor, remaining_capacity));
    }

    if (qc > amrex::Real(0)) {
        const amrex::Real qcc = qc;
        const amrex::Real available_cloud = amrex::max(
            amrex::Real(0), qc + saturation_adjustment.dq_vapor_to_cloud - saturation_adjustment.dq_cloud_to_vapor);
        amrex::Real auto_r = amrex::Real(0);
        if (qcc > qcw0) {
            // MUST MATCH: current AdvanceKessler autoconversion trigger.
            auto_r = alphaelq;
        }

        // MUST MATCH: current AdvanceKessler accretion coefficient and exponent.
        amrex::Real accrr = amrex::Real(0);
        accrr = amrex::Real(2.2) * std::pow(qp, amrex::Real(0.875));
        source_terms.dq_cloud_to_rain = dt * (accrr * qcc + auto_r * (qcc - qcw0));
        source_terms.dq_cloud_to_rain = amrex::min(source_terms.dq_cloud_to_rain, available_cloud);
    }

    return source_terms;
}

#endif