#ifndef SUPERDROPLET_PC_DEFN_H_
#define SUPERDROPLET_PC_DEFN_H_

#ifdef ERF_USE_PARTICLES

#include <AMReX_Enum.H>
#include <AMReX_Array.H>
#include <AMReX_REAL.H>
#include <AMReX_Particles.H>
#include "ERF_IndexDefines.H"
#include "ERF_Constants.H"
#include "ERFPC.H"
#include "ERF_SDInitialization.H"

namespace SDPCDefn {
    using SDSpeciesMassArr = amrex::GpuArray<amrex::ParticleReal*,SupDropInit::num_species_max>;
    using SDAerosolMassArr = amrex::GpuArray<amrex::ParticleReal*,SupDropInit::num_aerosols_max>;

    /*! \brief Real-type struct-of-array attributes of a super-droplet */
    using SuperDropletsRealIdx = ERFParticlesRealIdx;

    /*! \brief Int-type struct-of-array attributes of a super-droplet */
    using SuperDropletsIntIdx = ERFParticlesIntIdx;

    /*! \brief Int-type struct-of-array attributes added during runtime */
    struct SuperDropletsIntIdxSoA_RT
    {
        enum {
            active = 0,     /*!< Active/inactive */
            ncomps
        };
    };

    /*! \brief Real-type struct-of-array attributes added during runtime */
    struct SuperDropletsRealIdxSoA_RT
    {
        enum {
            radius = 0,     /*!< Radius of physical particles */
            multiplicity,   /*!< Number of droplets that this super-droplet represents
                                 (Real type to handle large values) */
            term_vel,       /*!< Terminal velocity */
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
            cond_tendency,  /*!< Condensation/evaporation tendency */
#endif
            uid,            /*!< unique ID */
            ncomps
        };
    };
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int ridx_a(const int a_i, /*!< Species index */
                  const int a_num_a, /*!< Number of aerosols */
                  const int a_num_s  /*!< Number of species */)
{
    amrex::ignore_unused(a_num_a);
    amrex::ignore_unused(a_num_s);
    return    SDPCDefn::SuperDropletsRealIdxSoA_RT::ncomps
            + a_i;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int ridx_s(const int a_i, /*!< Species index */
                  const int a_num_a, /*!< Number of aerosols */
                  const int a_num_s  /*!< Number of species */)
{
    amrex::ignore_unused(a_num_s);
    return    SDPCDefn::SuperDropletsRealIdxSoA_RT::ncomps
            + a_num_a
            + a_i;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int idx_a(const int a_i, /*!< Species index */
                 const int a_num_a, /*!< Number of aerosols */
                 const int a_num_s  /*!< Number of species */)
{
    return SDPCDefn::SuperDropletsRealIdx::ncomps + ridx_a(a_i,a_num_a,a_num_s);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int idx_s(const int a_i, /*!< Species index */
                 const int a_num_a, /*!< Number of aerosols */
                 const int a_num_s  /*!< Number of species */)
{
    return SDPCDefn::SuperDropletsRealIdx::ncomps + ridx_s(a_i,a_num_a,a_num_s);
}

/*! \brief Compute total mass of a droplet from its components */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static amrex::ParticleReal SD_total_mass( const int a_idx, /*!< Particle index */
                           const int a_num_sp,   /*!< number of species */
                           const int a_num_ae,   /*!< number of aerosols */
                           const SDPCDefn::SDSpeciesMassArr& a_sp_mass, /*!< species masses */
                           const SDPCDefn::SDAerosolMassArr& a_ae_mass  /*!< aerosol masses */ )
{
    amrex::ParticleReal retval = zero;
    for (int j = 0; j < a_num_sp; j++) {
        retval += a_sp_mass[j][a_idx];
    }
    for (int j = 0; j < a_num_ae; j++) {
        retval += a_ae_mass[j][a_idx];
    }
    return retval;
}

/*! \brief Compute effective radius of a droplet from its components */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static amrex::ParticleReal SD_effective_radius( const int a_idx, /*!< Particle index */
                                 const int a_idx_w, /*!< Species index of water */
                                 amrex::ParticleReal a_rho_w,  /*!< density of water */
                                 const int a_num_sp,   /*!< number of species */
                                 const int a_num_ae,   /*!< number of aerosols */
                                 const int* const a_sp_sol_arr, /*!< solubility of species */
                                 const int* const a_ae_sol_arr, /*!< solubility of aerosols */
                                 const SDPCDefn::SDSpeciesMassArr& a_sp_mass, /*!< species masses */
                                 const SDPCDefn::SDAerosolMassArr& a_ae_mass, /*!< aerosol masses */
                                 const amrex::ParticleReal* const a_sp_rho, /*!< species densities */
                                 const amrex::ParticleReal* const a_ae_rho  /*!< aerosol densities */ )
{
    amrex::ParticleReal m_w = a_sp_mass[a_idx_w][a_idx];
    amrex::ParticleReal m_s = zero;
    amrex::ParticleReal m_p = zero;
    amrex::ParticleReal rho_p = zero;
    for (int j = 0; j < a_num_sp; j++) {
        if (j != a_idx_w) {
            if (a_sp_sol_arr[j]) {
                m_s += a_sp_mass[j][a_idx];
            } else {
                m_p += a_sp_mass[j][a_idx];
                rho_p += a_sp_rho[j]*a_sp_mass[j][a_idx];
            }
        }
    }
    for (int j = 0; j < a_num_ae; j++) {
        if (a_ae_sol_arr[j]) {
            m_s += a_ae_mass[j][a_idx];
        } else {
            m_p += a_ae_mass[j][a_idx];
            rho_p += a_ae_rho[j]*a_ae_mass[j][a_idx];
        }
    }
    if (m_p > zero) { rho_p /= m_p; }
    else           { rho_p = one; }
    auto m_t = m_w + m_s + (a_rho_w/rho_p)*m_p;
    auto r_eff = std::cbrt(m_t / (four_thirds_pi*a_rho_w));
    return r_eff;
}

/*! \brief Compute "dry" radius of a droplet from its insoluble components */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static amrex::ParticleReal SD_dry_radius( const int a_idx, /*!< Particle index */
                           const int a_num_sp,   /*!< number of species */
                           const int a_num_ae,   /*!< number of aerosols */
                           const int* const a_sp_sol_arr, /*!< solubility of species */
                           const int* const a_ae_sol_arr, /*!< solubility of aerosols */
                           const SDPCDefn::SDSpeciesMassArr& a_sp_mass, /*!< species masses */
                           const SDPCDefn::SDAerosolMassArr& a_ae_mass, /*!< aerosol masses */
                           const amrex::ParticleReal* const a_sp_rho, /*!< species densities */
                           const amrex::ParticleReal* const a_ae_rho  /*!< aerosol densities */ )
{
    amrex::ParticleReal m_p = zero;
    amrex::ParticleReal rho_p = zero;
    for (int j = 0; j < a_num_sp; j++) {
        if (!a_sp_sol_arr[j]) {
            m_p += a_sp_mass[j][a_idx];
            rho_p += a_sp_rho[j]*a_sp_mass[j][a_idx];
        }
    }
    for (int j = 0; j < a_num_ae; j++) {
        if (!a_ae_sol_arr[j]) {
            m_p += a_ae_mass[j][a_idx];
            rho_p += a_ae_rho[j]*a_ae_mass[j][a_idx];
        }
    }
    if (m_p > zero) { rho_p /= m_p; }
    else           { rho_p = one; }
    auto r_dry = std::cbrt(m_p / (four_thirds_pi*rho_p));
    return r_dry;
}

/*! \brief Types of terminal velocity models */
AMREX_ENUM(SDTerminalVelocityType,
    RogersYau,
    AtlasUlbrich,
    CloudRainShima
);

/*! \brief Types of coalescence kernels */
enum struct SDCoalescenceKernelType {
    golovin,
    sedimentation,
    Longs,
    Halls
};

/*! \brief Types of phase change ODE time integrator */
enum struct SDMassChangeTIMethod {
    RK3BS, RK4, BE, CN, DIRK2
};

#endif
#endif
