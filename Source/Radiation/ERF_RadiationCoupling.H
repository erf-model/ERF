#ifndef ERF_RADIATION_COUPLING_H_
#define ERF_RADIATION_COUPLING_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <ERF_RadStruct.H>

/**
 * @file ERF_RadiationCoupling.H
 * @brief Interface for coupling radiation forcing to the atmospheric state.
 *
 * Provides a placeholder interface for coupling radiation forcing to the prognostic state.
 * Full coupling to RhoTheta fields and microphysics is planned for future work.
 *
 * This file will be extended to:
 * - Inject SW/LW heating rates into prognostic equations
 * - Interface with moisture microphysics
 * - Handle surface flux boundary conditions
 */

/**
 * @brief Inject radiation heating into the prognostic state.
 *
 * Currently a no-op placeholder. When fully implemented, this function will:
 * - Update RhoTheta with radiative heating rates
 * - Apply surface LW/SW flux boundary conditions
 * - Couple to microphysics (cloud optical properties)
 *
 * @param[in,out] rho_theta  RhoTheta field (density × potential temperature).
 * @param[in]     rad_choice Radiation parameters and control flags.
 * @param[in]     geom       Grid geometry.
 * @param[in]     lev        AMR level index.
 *
 * @note Currently a no-op. Integration into the main timestep loop is
 *       planned for future work.
 */
inline void
inject_radiation_heating(amrex::MultiFab& /*rho_theta*/,
                          const RadChoice& /*rad_choice*/,
                          const amrex::Geometry& /*geom*/,
                          int /*lev*/)
{
  // No-op stub: heating injection not yet implemented
  // Will be extended in future work to couple heating to RhoTheta
}

#endif // ERF_RADIATION_COUPLING_H_
