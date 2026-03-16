#ifndef ERF_MOISTUTILS_H_
#define ERF_MOISTUTILS_H_

#include "ERF_EOS.H"
#include "ERF_MicrophysicsUtils.H"

// =============================================================================
// MOIST TURBULENCE FUNCTIONS
// =============================================================================

/**
 * Extract moisture variables and partition into liquid/ice phases
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void GetMoistureVars (int i, int j, int k,
                      const amrex::Array4<amrex::Real const>& cell_data,
                      amrex::Real& qv,
                      amrex::Real& qc_liquid,
                      amrex::Real& qc_ice,
                      const MoistureComponentIndices& moisture_indices)
{
    qv = 0.0; qc_liquid = 0.0; qc_ice = 0.0;

    // Water vapor
    if (moisture_indices.qv >= 0) {
        qv = cell_data(i,j,k,moisture_indices.qv)/cell_data(i,j,k,Rho_comp);
    }

    // Cloud liquid water
    if (moisture_indices.qc >= 0) {
        qc_liquid = cell_data(i,j,k,moisture_indices.qc)/cell_data(i,j,k,Rho_comp);
    }

    // Cloud ice (only if separate ice species exists)
    if (moisture_indices.qi >= 0) {
        qc_ice = cell_data(i,j,k,moisture_indices.qi)/cell_data(i,j,k,Rho_comp);
    }
    // No temperature-based partitioning - respect the microphysics scheme's decisions

    // Add precipitating species
    if (moisture_indices.qr >= 0) {  // Rain (always liquid)
        qc_liquid += cell_data(i,j,k,moisture_indices.qr)/cell_data(i,j,k,Rho_comp);
    }
    if (moisture_indices.qs >= 0) {  // Snow (always ice)
        qc_ice += cell_data(i,j,k,moisture_indices.qs)/cell_data(i,j,k,Rho_comp);
    }
    if (moisture_indices.qg >= 0) {  // Graupel (always ice)
        qc_ice += cell_data(i,j,k,moisture_indices.qg)/cell_data(i,j,k,Rho_comp);
    }
}

/**
 * Compute virtual potential temperature with moisture loading effects
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real ComputeVirtualPotentialTemperature (amrex::Real theta,
                                                amrex::Real qv,
                                                amrex::Real qc_liquid,
                                                amrex::Real qc_ice)
{
    amrex::Real qc_total = qc_liquid + qc_ice;
    return theta * (1.0 + 0.61*qv - qc_total);
}

/**
 * Wrapper around ComputeVirtualPotentialTemperature
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real GetThetav (int i, int j, int k,
                       const amrex::Array4<amrex::Real const>& cell_data,
                       const MoistureComponentIndices& moisture_indices)
{
    amrex::Real theta = cell_data(i, j, k, RhoTheta_comp) / cell_data(i, j, k, Rho_comp);

    amrex::Real qv, qcl, qci;
    GetMoistureVars(i, j, k, cell_data, qv, qcl, qci, moisture_indices);

    return ComputeVirtualPotentialTemperature(theta, qv, qcl, qci);
}

/**
 * Compute linearized liquid-water potential temperature
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real ComputeLiquidWaterPotentialTemperature (amrex::Real theta,
                                                    amrex::Real T,
                                                    amrex::Real qc_liquid)
{
    return theta - (theta / T) * (L_v / Cp_d) * qc_liquid;
}

/**
 * Wrapper around ComputeLiquidWaterPotentialTemperature
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real GetThetal (int i, int j, int k,
                       const amrex::Array4<amrex::Real const>& cell_data,
                       const MoistureComponentIndices& moisture_indices)
{
    amrex::Real qv, qcl, qci;
    GetMoistureVars(i, j, k, cell_data, qv, qcl, qci, moisture_indices);

    amrex::Real theta = cell_data(i, j, k, RhoTheta_comp) / cell_data(i, j, k, Rho_comp);

    amrex::Real T = getTgivenRandRTh(cell_data(i, j, k, Rho_comp),
                                     cell_data(i, j, k, RhoTheta_comp),
                                     qv);

    return ComputeLiquidWaterPotentialTemperature(theta, T, qcl);
}


/**
 * Compute moist stratification accounting for conditional instability
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real ComputeMoistStratification (int i, int j, int k,
                                        const amrex::Array4<amrex::Real const>& cell_data,
                                        amrex::Real dzInv,
                                        amrex::Real abs_g,
                                        amrex::Real inv_theta0,
                                        const MoistureComponentIndices& moisture_indices)
{
    // Get moisture variables partitioned into phases
    amrex::Real qv, qc_liquid, qc_ice;
    GetMoistureVars(i, j, k, cell_data, qv, qc_liquid, qc_ice, moisture_indices);

    // Compute virtual potential temperature gradient
    amrex::Real theta_v_upper = GetThetav(i, j, k+1, cell_data, moisture_indices);
    amrex::Real theta_v_lower = GetThetav(i, j, k-1, cell_data, moisture_indices);

    amrex::Real dthetav_dz = 0.5 * (theta_v_upper - theta_v_lower) * dzInv;
    amrex::Real stratification = abs_g * dthetav_dz * inv_theta0;

    // Apply conditional instability correction
    amrex::Real total_condensate = qc_liquid + qc_ice;
    if (total_condensate > 1e-8) {
        amrex::Real T_current = getTgivenRandRTh(cell_data(i,j,k,Rho_comp),
                                                 cell_data(i,j,k,RhoTheta_comp), qv);

        // Phase-weighted effective latent heat
        amrex::Real liquid_fraction = qc_liquid / (total_condensate + 1e-12);
        amrex::Real ice_fraction = 1.0 - liquid_fraction;
        amrex::Real L_eff = liquid_fraction * L_v + ice_fraction * (L_v + lat_ice);

        // Phase-weighted saturation mixing ratio
        amrex::Real qsat_liquid = 0.0, qsat_ice = 0.0;
        amrex::Real pres_current = getPgivenRTh(cell_data(i,j,k,RhoTheta_comp), qv) * 0.01;

        erf_qsatw(T_current, pres_current, qsat_liquid);
        erf_qsati(T_current, pres_current, qsat_ice);

        amrex::Real qsat_eff = liquid_fraction * qsat_liquid + ice_fraction * qsat_ice;

        // Moist adiabatic lapse rate correction
        amrex::Real gamma_moist_factor = (1.0 + L_eff*qsat_eff/(R_d*T_current)) /
                                         (1.0 + L_eff*L_eff*qsat_eff/(Cp_d*R_v*T_current*T_current));

        stratification *= gamma_moist_factor;
    }

    return stratification;
}

/**
 * Compute stratification for Smagorinsky scheme (moist or dry)
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real ComputeStratificationForSmagorinsky (int i, int j, int k,
                                                 const amrex::Array4<amrex::Real const>& cell_data,
                                                 amrex::Real dzInv,
                                                 amrex::Real abs_g,
                                                 amrex::Real inv_theta0,
                                                 bool use_moisture,
                                                 int rho_qv_comp,
                                                 const MoistureComponentIndices& moisture_indices)
{
    if (use_moisture && (rho_qv_comp >= 0)) {
        // Moist stratification with virtual temperature and conditional instability
        return ComputeMoistStratification(i, j, k, cell_data, dzInv, abs_g, inv_theta0, moisture_indices);
    } else {
        // Dry stratification (original approach)
        amrex::Real theta_hi = cell_data(i,j,k+1,RhoTheta_comp)/cell_data(i,j,k+1,Rho_comp);
        amrex::Real theta_lo = cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp);
        amrex::Real dtheta_dz = 0.5 * (theta_hi - theta_lo) * dzInv;
        return abs_g * dtheta_dz * inv_theta0;
    }
}

#endif
