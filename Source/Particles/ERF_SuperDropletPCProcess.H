/*! \file ERF_SuperDropletPCProcess.H
 *  \brief Common data structures for SuperDroplet physical processes
 *
 *  This file provides shared infrastructure for physical process implementations:
 *  1. ParticlePointers - Bundle of particle attribute pointers for GPU kernels
 *  2. ProcessContext - Non-pointer data (geometry, species info) shared across processes
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

/*! \brief Bundle of all particle attribute pointers
 *
 * This struct holds pointers to particle data for use in GPU kernels.
 * All pointers are populated by setupParticlePointers(); processes use
 * only the pointers they need.
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

    // Species and aerosol mass arrays
    SDSpeciesMassArr sp_mass_ptrs{};
    SDAerosolMassArr ae_mass_ptrs{};

    // Material properties (read-only)
    const ParticleReal* sp_rho_arr{nullptr};
    const int* sp_sol_arr{nullptr};
    const ParticleReal* sp_ion_arr{nullptr};
    const ParticleReal* sp_mw_arr{nullptr};

    const ParticleReal* ae_rho_arr{nullptr};
    const int* ae_sol_arr{nullptr};
    const ParticleReal* ae_ion_arr{nullptr};
    const ParticleReal* ae_mw_arr{nullptr};
};

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

    // Material densities (commonly used)
    ParticleReal rho_water{0.0};
};

} // namespace SDProcess

#endif // ERF_USE_PARTICLES
#endif // ERF_SUPERDROPLET_PC_PROCESS_H_
