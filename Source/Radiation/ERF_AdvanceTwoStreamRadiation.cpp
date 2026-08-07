#include <ERF.H>
#include <ERF_RadStruct.H>
#include <ERF_RadiationDiagnostics.H>
#include <ERF_TwoStreamSW.H>
#include <ERF_TwoStreamLW.H>
#include <AMReX_Print.H>
#include <cmath>

using namespace amrex;

/**
 * @file ERF_AdvanceTwoStreamRadiation.cpp
 * @brief Phase 1 two-stream radiation driver for diagnostic computation.
 *
 * Computes SW/LW fluxes and heating rates using the Phase 1 two-stream model
 * and writes diagnostics to CSV file.
 *
 * This is a simplified implementation for Phase 1 that:
 * - Does NOT inject heating into RhoTheta (that's Phase 8)
 * - Only computes diagnostic flux values
 * - Writes to CSV file via RadiationDiagnostics
 *
 * Note: The actual heating injection and proper time-stepping coupling
 * will be implemented in later phases.
 */

void ERF::compute_twostream_radiation_diagnostics(
    int lev,
    int nstep,
    amrex::Real time_step)
{
    const auto& rad_choice = solverChoice.radChoice;

    // Only proceed if TwoStream radiation is enabled
    if (rad_choice.rad_type != RadType::TwoStream) {
        return;
    }

    // Create RadiationDiagnostics instance for this level
    // Note: Phase 1 only uses level 0
    RadiationDiagnostics rad_diag(rad_choice.verbosity, rad_choice.diag_file, lev);

    // ========================================
    // Phase 1 stub: Compute diagnostic values
    // ========================================
    
    // For Phase 1, we use simple analytical values based on input parameters
    // Later phases will integrate with actual atmospheric state and thermodynamics
    
    amrex::Real SW_surface = 0.0;
    amrex::Real SW_TOA = 0.0;
    amrex::Real F_up_surface = 0.0;
    amrex::Real F_down_toa = 0.0;
    amrex::Real heating_rate_max = 0.0;

    // Shortwave calculation (Phase 1)
    if (rad_choice.sw_enabled) {
        // Convert solar zenith angle to radians
        amrex::Real zenith_rad = rad_choice.solar_zenith_deg * M_PI / 180.0;
        amrex::Real cos_zenith = std::cos(zenith_rad);

        // TOA downwelling SW flux (solar constant × geometric factor)
        SW_TOA = rad_choice.S0 * std::max(0.0, cos_zenith);

        // For Phase 1 test: assume single layer with tau_per_layer
        // Surface flux from Beer-Lambert: F(surface) = S0 * cos(zenith) * exp(-tau/cos(zenith))
        // Using simplified tau_cumulative = tau_per_layer (single layer)
        if (cos_zenith > 0.0) {
            //amrex::Real tau_cum = rad_choice.tau_per_layer;
            int n_layers = geom[lev].Domain().length(2);
            amrex::Real tau_cum = rad_choice.tau_per_layer * static_cast<amrex::Real>(n_layers);
            SW_surface = rad_choice.S0 * cos_zenith * 
                        std::exp(-tau_cum / cos_zenith);
            // Rough estimate of max heating rate (simplified)
            heating_rate_max = std::abs(SW_TOA - SW_surface) / 1000.0; // dummy calc
        }
    }

    // Longwave calculation (Phase 1)
    if (rad_choice.lw_enabled) {
        // For Phase 1 isothermal test: if isothermal_test is true,
        // F_up and F_down at surface/TOA are computed from Stefan-Boltzmann
        if (rad_choice.isothermal_test) {
            // Stefan-Boltzmann constant: sigma = 5.670374419e-8 W/(m^2·K^4)
            amrex::Real sigma = 5.670374419e-8;
            amrex::Real T = rad_choice.T_iso_K;
            amrex::Real I_thermal = sigma * T * T * T * T;
            
            // In isothermal case: F_up(surf) = F_down(surf) = I_thermal everywhere
            // Net flux is zero (no heating)
            F_up_surface = I_thermal;
            F_down_toa = I_thermal;
            // heating_rate_max remains 0 for isothermal case
        }
    }

    // Append diagnostics to CSV file
    rad_diag.append(nstep, time_step, SW_surface, SW_TOA,
                    F_up_surface, F_down_toa, heating_rate_max);
}
