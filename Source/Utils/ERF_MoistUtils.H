#ifndef ERF_MOISTUTILS_H_
#define ERF_MOISTUTILS_H_

#include "ERF_EOS.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_MicrophysicsConstants.H"
#include "ERF_Constants.H"

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
    qv = zero; qc_liquid = zero; qc_ice = zero;

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
 * Extract water vapor and the cloud condensate alone: cloud liquid qc and cloud ice qi, without
 * the precipitating species (rain, snow, graupel) that GetMoistureVars adds in
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void GetCloudMoistureVars (int i, int j, int k,
                           const amrex::Array4<amrex::Real const>& cell_data,
                           amrex::Real& qv,
                           amrex::Real& qc_liquid,
                           amrex::Real& qc_ice,
                           const MoistureComponentIndices& moisture_indices)
{
    const amrex::Real rho = cell_data(i,j,k,Rho_comp);
    qv        = (moisture_indices.qv >= 0) ? cell_data(i,j,k,moisture_indices.qv) / rho : zero;
    qc_liquid = (moisture_indices.qc >= 0) ? cell_data(i,j,k,moisture_indices.qc) / rho : zero;
    qc_ice    = (moisture_indices.qi >= 0) ? cell_data(i,j,k,moisture_indices.qi) / rho : zero;
}

/**
 * Compute virtual potential temperature with moisture loading effects
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real ComputeVirtualPotentialTemperature (const amrex::Real& theta,
                                                const amrex::Real& qv,
                                                const amrex::Real& qc_liquid,
                                                const amrex::Real& qc_ice)
{
    amrex::Real qc_total = qc_liquid + qc_ice;
    return theta * (one + epsv * qv - qc_total);
}

/**
 * Wrapper around ComputeVirtualPotentialTemperature
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real GetThetav (const int& i, const int& j, const int& k,
                       const amrex::Array4<amrex::Real const>& cell_data,
                       const MoistureComponentIndices& moisture_indices)
{
    amrex::Real theta = cell_data(i, j, k, RhoTheta_comp) / cell_data(i, j, k, Rho_comp);

    amrex::Real qv, qcl, qci;
    GetMoistureVars(i, j, k, cell_data, qv, qcl, qci, moisture_indices);

    return ComputeVirtualPotentialTemperature(theta, qv, qcl, qci);
}

/**
 * Compute linearized liquid-ice water potential temperature,
 * theta_l = theta - (theta / T) (L_v q_c + L_s q_i) / c_p
 * from the cloud liquid q_c and cloud ice q_i (WRF's thlix in YSU)
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real ComputeLiquidWaterPotentialTemperature (const amrex::Real& theta,
                                                    const amrex::Real& T,
                                                    const amrex::Real& qc_liquid,
                                                    const amrex::Real& qc_ice)
{
    return theta - (theta / T) * (L_v * qc_liquid + lsub * qc_ice) / Cp_d;
}

/**
 * Wrapper around ComputeLiquidWaterPotentialTemperature.  theta_l is the variable conserved under
 * reversible condensation, so it uses the cloud condensate only: precipitation falls out of the
 * parcel and does not enter the latent-heat term (erf-model/ERF#4892)
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real GetThetal (const int& i, const int& j, const int& k,
                       const amrex::Array4<amrex::Real const>& cell_data,
                       const MoistureComponentIndices& moisture_indices)
{
    amrex::Real qv, qcl, qci;
    GetCloudMoistureVars(i, j, k, cell_data, qv, qcl, qci, moisture_indices);

    amrex::Real theta = cell_data(i, j, k, RhoTheta_comp) / cell_data(i, j, k, Rho_comp);

    amrex::Real T = getTgivenRandRTh(cell_data(i, j, k, Rho_comp),
                                     cell_data(i, j, k, RhoTheta_comp),
                                     qv);

    return ComputeLiquidWaterPotentialTemperature(theta, T, qcl, qci);
}

#endif
