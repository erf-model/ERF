#ifndef ERF_FIRE_SMOKE_EMISSION_H
#define ERF_FIRE_SMOKE_EMISSION_H

#ifdef ERF_ENABLE_FIRE

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

/**
 * @brief Inject smoke tracer from fire heat flux into cc_source at k=0.
 *
 * Converts fire heat flux [W/m2] to smoke mass flux [kg/m2/s] via:
 *   smoke_flux = emission_factor * Q_fire / heat_of_combustion
 * then injects into cc_source[RhoSmoke_comp] at k=0 divided by dz.
 *
 * The fire heat flux must already be coarsened to the ATM grid BA
 * by the caller (same ParallelCopy pattern as Phase 3 lofting).
 *
 * @param[in,out] cc_source           ATM cell-centred source MultiFab
 * @param[in]     fire_heat_atm       Fire heat flux coarsened to ATM grid [W/m2]
 * @param[in]     z_phys_cc           Cell-centre terrain height [m]
 * @param[in]     geom_atm            ATM geometry
 * @param[in]     emission_factor     kg smoke per kg fuel burned [-]
 * @param[in]     heat_of_combustion  J per kg fuel burned [J/kg]
 * @param[in]     smoke_comp          Component index in cc_source for smoke
 * @param[in]     fire_debug          Print debug output if true
 * @param[in]     step                Current step for debug output
 */
void inject_smoke_from_fire(
    amrex::MultiFab&       cc_source,
    const amrex::MultiFab& fire_heat_atm,
    const amrex::MultiFab& z_phys_cc,
    const amrex::Geometry& geom_atm,
    amrex::Real            emission_factor,
    amrex::Real            heat_of_combustion,
    int                    smoke_comp,
    bool                   fire_debug,
    int                    step);

#endif // ERF_ENABLE_FIRE
#endif // ERF_FIRE_SMOKE_EMISSION_H
