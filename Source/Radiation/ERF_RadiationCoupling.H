#ifndef ERF_RADIATION_COUPLING_H_
#define ERF_RADIATION_COUPLING_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <ERF_RadStruct.H>

/**
 * @file ERF_RadiationCoupling.H
 * @brief Interface for coupling radiation forcing to the atmospheric state.
 *
 * Phase 1 provides only a stub implementation (no-op). Full coupling
 * to RhoTheta fields is deferred to Phase 8.
 *
 * This file will be extended in later phases to:
 * - Inject SW/LW heating rates into prognostic equations
 * - Interface with moisture microphysics
 * - Handle surface flux boundary conditions
 */

/**
 * @brief Inject radiation heating into the prognostic state (Phase 1 stub).
 *
 * In Phase 1, this function does nothing (no-op placeholder).
 * Later phases will use this to:
 * - Update RhoTheta with radiative heating rates
 * - Apply surface LW/SW flux boundary conditions
 * - Couple to microphysics (cloud optical properties)
 *
 * @param[in,out] rho_theta  RhoTheta field (density × potential temperature).
 *                           Will remain unchanged in Phase 1.
 * @param[in]     rad_choice Radiation parameters and control flags.
 * @param[in]     geom       Grid geometry.
 * @param[in]     lev        AMR level index.
 *
 * @note This is a no-op in Phase 1. Integration into the main timestep
 *       loop is deferred to Phase 8.
 */
inline void
inject_radiation_heating(amrex::MultiFab& /*rho_theta*/,
                          const RadChoice& /*rad_choice*/,
                          const amrex::Geometry& /*geom*/,
                          int /*lev*/)
{
  // Phase 1 stub: no heating injection yet
  // Will be extended in Phase 8 to actually couple heating to RhoTheta
}

#endif // ERF_RADIATION_COUPLING_H_
