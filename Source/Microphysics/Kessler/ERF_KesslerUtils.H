#ifndef ERF_KESSLER_UTILS_H_
#define ERF_KESSLER_UTILS_H_

#include <algorithm>
#include <cmath>
#include <limits>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <ERF_NumericalConstants.H>
#include <ERF_MicrophysicsConstants.H>
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
// More rain substeps than any healthy state asks for (AdvanceKessler aborts above this)
static constexpr int kessler_max_sedimentation_substeps = 1000;
static constexpr amrex::Real kessler_water_depth_m_per_kg_m2 = amrex::Real(1.0) / rhor;
static constexpr amrex::Real kessler_millimeters_per_meter = amrex::Real(1000.0);
// Klemp & Wilhelmson (1978) rain fall speed, as in WRF module_mp_kessler.F and DCMIP2016 kessler.f90:
// V [m/s] = 36.34 * (0.001*rho*qp)^0.1364 * (rho/1.16)^(-1/2), rho in kg/m^3 (0.001*rho is in g/cm^3).
static constexpr amrex::Real kessler_vt_coeff = amrex::Real(36.34);
static constexpr amrex::Real kessler_vt_rho_to_cgs = amrex::Real(0.001);
static constexpr amrex::Real kessler_vt_exponent = amrex::Real(0.1364);
static constexpr amrex::Real kessler_vt_rho0 = amrex::Real(1.16);
static constexpr amrex::Real kessler_vt_density_exponent = -myhalf;

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
    // This function, with the kessler_vt_* constants, defines the Kessler rain fall speed.
    return kessler_vt_coeff
        * std::pow(rho * kessler_vt_rho_to_cgs * qp, kessler_vt_exponent)
        * std::pow(rho / kessler_vt_rho0, kessler_vt_density_exponent);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real kessler_precip_flux (const amrex::Real rho,
                                 const amrex::Real terminal_velocity,
                                 const amrex::Real qp) noexcept
{
    return rho * terminal_velocity * qp;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real kessler_limit_precip_flux (const amrex::Real flux,
                                       const amrex::Real donor_rho,
                                       const amrex::Real donor_qp,
                                       const amrex::Real donor_detJ,
                                       const amrex::Real coef) noexcept
{
    // The donor cap: a face cannot export more rain in one substep than its donor cell
    // holds, F * dt / dz <= rho_donor * qp_donor * detJ_donor, with coef = dt / dz.
    return amrex::min(flux, donor_rho * donor_qp * donor_detJ / coef);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int kessler_num_sedimentation_substeps (const amrex::Real max_courant) noexcept
{
    // max_courant is the largest fall Courant number over the faces for the full step,
    // V*dt over the height of the cell the face drains: enough substeps that rain falls
    // at most CFL_MAX of that cell per substep, and always at least one.  A NaN, infinite
    // or huge max_courant returns kessler_max_sedimentation_substeps + 1 rather than casting
    // it to an int; AdvanceKessler aborts before it gets here.
    const amrex::Real ratio = max_courant / kessler_sedimentation_cfl_max;
    if (!(ratio <= amrex::Real(kessler_max_sedimentation_substeps))) {
        return kessler_max_sedimentation_substeps + 1;
    }
    return amrex::max(1, static_cast<int>(std::ceil(ratio)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real kessler_fall_courant (const amrex::Real terminal_velocity,
                                  const amrex::Real dt,
                                  const amrex::Real dz,
                                  const amrex::Real detJ_donor) noexcept
{
    // The fall Courant number at a face: the distance fallen in dt over the physical height
    // dz*detJ of the donor cell, the one the face drains, as the donor-cap limiter uses.
    return terminal_velocity * dt / (dz * detJ_donor);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real kessler_sizing_detJ (const int k,
                                 const int k_lo,
                                 const int k_hi,
                                 const amrex::Real detJ_donor,
                                 const amrex::Real detJ_below) noexcept
{
    // The detJ that sizes the substeps at face k: the smaller of the donor's (the cell the
    // face drains) and that of the cell below it (the one its rain falls into), for faces
    // inside the domain.  The substep count is fixed at the start of the step, while the
    // cells below the rain front are still dry; looking one cell down covers the first cell
    // the front enters.  The surface face drains out of the domain and the top face carries
    // no rain, so both use the donor alone.
    return (k > k_lo && k <= k_hi) ? amrex::min(detJ_donor, detJ_below) : detJ_donor;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real kessler_rain_accumulation_increment (const amrex::Real precip_mass_per_area) noexcept
{
    // Convert precipitation mass per unit area [kg m^-2] to liquid-water depth
    // [mm]: divide by water density [kg m^-3], then multiply by mm per meter.
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
                                     const amrex::Real /*qp_km1*/,
                                     const amrex::Real qp_k) noexcept
{
    KesslerFaceState face_state{amrex::Real(0), amrex::Real(0)};

    // Sedimentation is one-way downward transport. Use a face-centered density
    // for the z-face flux, but upwind the transported rain mixing ratio from the
    // donor cell. A centered qp can create flux from downstream rain when the
    // donor cell is dry. The top boundary face carries no rain: nothing lies
    // above the domain, so its donor is empty. It keeps the top cell's rho so
    // the fall speed stays defined.
    if (k == k_hi + 1) {
        face_state.rho = rho_km1;
        face_state.qp = amrex::Real(0);
    } else if (k == 0) {
        face_state.rho = rho_k;
        face_state.qp = qp_k;
    } else {
        face_state.rho = amrex::Real(0.5) * (rho_km1 + rho_k);
        face_state.qp = qp_k;
    }

    // MUST MATCH: current AdvanceKessler nonnegative donor-qp clip before flux formation.
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
    // Kessler's mbar / hPa path for qsat helpers, and dt is in seconds.
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
    // Cloud evaporation gets first claim on the linearized subsaturation
    // capacity. Rain evaporation can use only the remaining capacity.
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
        // Cloud-to-rain conversion is capped by cloud water available after
        // saturation adjustment, not by the original qc.
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
