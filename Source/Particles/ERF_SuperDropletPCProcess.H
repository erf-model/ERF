/*! \file ERF_SuperDropletPCProcess.H
 *  \brief Policy-based template infrastructure for SuperDroplet physical processes
 *
 *  This file provides a framework to reduce code duplication across physical
 *  process implementations by:
 *  1. Defining compile-time traits for process requirements
 *  2. Bundling commonly-used pointers into a single struct
 *  3. Providing template-based pointer setup that only fetches needed data
 */

#ifndef ERF_SUPERDROPLET_PC_PROCESS_H_
#define ERF_SUPERDROPLET_PC_PROCESS_H_

#ifdef ERF_USE_PARTICLES

#include "ERF_SuperDropletPCDefinitions.H"
#include <AMReX_GpuContainers.H>
#include <AMReX_Particles.H>

namespace SDProcess {

using namespace amrex;
using namespace SDPCDefn;

//===========================================================================
// PART 1: Traits System - Compile-time feature requirements
//===========================================================================

/*! \brief Base traits with default values
 *
 * Process policies inherit from this and override what they need.
 * Using static constexpr bool allows if constexpr optimization.
 */
struct DefaultTraits {
    // Core pointer requirements
    static constexpr bool needs_mass         = true;   // Almost always needed
    static constexpr bool needs_radius       = true;
    static constexpr bool needs_active       = true;
    static constexpr bool needs_velocity     = false;
    static constexpr bool needs_multiplicity = false;
    static constexpr bool needs_term_vel     = false;

    // Ice-related pointers
    static constexpr bool needs_ice_axes     = false;  // a, c pointers
    static constexpr bool needs_ice_rime     = false;  // mrime, nmono
    static constexpr bool needs_ice_Tfz      = false;  // freezing temperature

    // Material property requirements
    static constexpr bool needs_density      = true;   // sp_rho, ae_rho
    static constexpr bool needs_solubility   = true;   // sp_sol, ae_sol
    static constexpr bool needs_ionization   = false;
    static constexpr bool needs_mol_weight   = false;
    static constexpr bool needs_INP_flags    = false;
};

//===========================================================================
// Intermediate Traits - Common combinations for ice processes
//===========================================================================

/*! \brief Traits for processes needing ice freezing temperature only
 *  Used by: MassChangeLVTraits
 */
struct IceTfzTraits : DefaultTraits {
    static constexpr bool needs_ice_Tfz = true;
};

/*! \brief Traits for processes needing ice axes (a, c) and Tfz
 *  Used by: MassChangeSLTraits, AdvectionTraits (for IceBohm terminal velocity)
 */
struct IceAxesTraits : DefaultTraits {
    static constexpr bool needs_ice_Tfz  = true;
    static constexpr bool needs_ice_axes = true;
};

/*! \brief Traits for processes needing full ice properties (axes, rime, Tfz)
 *  Used by: MassChangeSVTraits, CoalescenceTraits
 */
struct IceFullTraits : DefaultTraits {
    static constexpr bool needs_ice_Tfz  = true;
    static constexpr bool needs_ice_axes = true;
    static constexpr bool needs_ice_rime = true;
};

/*! \brief Traits for processes needing velocity and terminal velocity
 *  Used by: AdvectionTraits, CoalescenceTraits
 */
struct VelocityTraits : DefaultTraits {
    static constexpr bool needs_velocity = true;
    static constexpr bool needs_term_vel = true;
};

//===========================================================================
// PART 2: Unified Pointer Bundle
//===========================================================================

/*! \brief Bundle of all particle attribute pointers
 *
 * This struct holds all commonly-used pointers to particle data.
 * Pointers that aren't needed (based on traits) remain nullptr.
 * The struct is designed to be captured by GPU lambdas.
 */
struct ParticlePointers {
    // Number of particles in this tile
    int num_particles{0};

    // Core SOA pointers
    ParticleReal* mass_ptr{nullptr};
    ParticleReal* radius_ptr{nullptr};
    int* active_ptr{nullptr};

    // Velocity components
    GpuArray<ParticleReal*, AMREX_SPACEDIM> v_ptr{nullptr, nullptr, nullptr};
    ParticleReal* vterm_ptr{nullptr};
    ParticleReal* mult_ptr{nullptr};

    // Ice properties
    ParticleReal* Tfz_ptr{nullptr};
    ParticleReal* a_ptr{nullptr};     // equatorial radius
    ParticleReal* c_ptr{nullptr};     // polar radius
    ParticleReal* mrime_ptr{nullptr};
    ParticleReal* nmono_ptr{nullptr};

    // Species and aerosol mass arrays
    SDSpeciesMassArr sp_mass_ptrs{};
    SDAerosolMassArr ae_mass_ptrs{};

    // Material properties (read-only)
    const ParticleReal* sp_rho_arr{nullptr};
    const int* sp_sol_arr{nullptr};
    const ParticleReal* sp_ion_arr{nullptr};
    const ParticleReal* sp_mw_arr{nullptr};
    const int* sp_INP_arr{nullptr};

    const ParticleReal* ae_rho_arr{nullptr};
    const int* ae_sol_arr{nullptr};
    const ParticleReal* ae_ion_arr{nullptr};
    const ParticleReal* ae_mw_arr{nullptr};
    const int* ae_INP_arr{nullptr};
};

//===========================================================================
// PART 3: Process Context - Non-pointer data needed across processes
//===========================================================================

/*! \brief Context data that doesn't change per-tile */
struct ProcessContext {
    // Geometry
    GpuArray<Real, AMREX_SPACEDIM> plo;
    GpuArray<Real, AMREX_SPACEDIM> phi;
    GpuArray<Real, AMREX_SPACEDIM> dxi;
    GpuArray<Real, AMREX_SPACEDIM> dx;
    GpuArray<int, AMREX_SPACEDIM> is_periodic;
    Box domain;
    Real cell_volume{0.0};

    // Species/aerosol counts and indices
    int num_species{0};
    int num_aerosols{0};
    int idx_water{-1};
    int idx_ice{-1};

    // Material densities (commonly used)
    ParticleReal rho_water{0.0};
    ParticleReal rho_ice{0.0};
};

} // namespace SDProcess

#endif // ERF_USE_PARTICLES
#endif // ERF_SUPERDROPLET_PC_PROCESS_H_
