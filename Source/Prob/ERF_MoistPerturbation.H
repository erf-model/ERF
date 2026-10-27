#ifndef ERF_MOIST_PERTURBATION_H_
#define ERF_MOIST_PERTURBATION_H_

/**
 * Initial temperature and vapour perturbations of a moist column, applied together.
 *
 * Shared by the BOMEX, GATE and RICO setups (and DYCOMS, which uses the RICO one).
 * rho*theta is held, so the temperature perturbation is stored as a density change;
 * the vapour mixing ratio after the perturbation is qv + dqv, whatever the density.
 */

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include "ERF_EOS.H"

namespace erf_moist_perturbation {

/**
 * Perturbation of rho and rho*qv for a temperature change dT and a vapour change dqv.
 *
 * With perturb_T, the new density gives T = T_old + dT at the pressure of the new state,
 * P(rho*theta, qv + dqv); with dqv = 0 that is the old pressure. Without perturb_T the
 * density is left as it is. rho*qv then holds qv + dqv on the new density.
 *
 * @param[in]  rho       dry density of the unperturbed state
 * @param[in]  rhotheta  rho*theta of the unperturbed state (held)
 * @param[in]  rhoqv     rho*qv of the unperturbed state (0 without moisture)
 * @param[in]  dT        temperature change [K]
 * @param[in]  dqv       vapour mixing ratio change
 * @param[in]  perturb_T whether the temperature is perturbed
 * @param[in]  moist     whether a moisture model is active
 * @param[in]  rdOcp     R_d / c_p
 * @param[out] drho      perturbation of rho
 * @param[out] drhoqv    perturbation of rho*qv (0 without moisture)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
perturb_temperature_and_vapour (const amrex::Real rho, const amrex::Real rhotheta,
                                const amrex::Real rhoqv, const amrex::Real dT,
                                const amrex::Real dqv, const bool perturb_T,
                                const bool moist, const amrex::Real rdOcp,
                                amrex::Real& drho, amrex::Real& drhoqv) noexcept
{
    const amrex::Real qv     = moist ? rhoqv / rho : amrex::Real(0.0);
    const amrex::Real qv_new = moist ? qv + dqv : amrex::Real(0.0);

    amrex::Real rho_new = rho;
    if (perturb_T) {
        const amrex::Real T_old     = getTgivenRandRTh(rho, rhotheta, qv);
        const amrex::Real P         = getPgivenRTh(rhotheta, qv_new);
        const amrex::Real theta_new = getThgivenTandP(T_old + dT, P, rdOcp);
        rho_new = getRhogivenThetaPress(theta_new, P, rdOcp, qv_new);
    }

    drho   = rho_new - rho;
    drhoqv = moist ? rho_new * qv_new - rho * qv : amrex::Real(0.0);
}

} // namespace erf_moist_perturbation

#endif
