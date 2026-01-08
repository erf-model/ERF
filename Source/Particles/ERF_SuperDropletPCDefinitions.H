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

using SDSpeciesMassArr = amrex::GpuArray<amrex::ParticleReal*,SupDropInit::num_species_max>;
using SDAerosolMassArr = amrex::GpuArray<amrex::ParticleReal*,SupDropInit::num_aerosols_max>;

/*! \brief Real-type array-of-struct attributes of a super-droplet */
using SuperDropletsRealIdxAoS = ERFParticlesRealIdxAoS;

/*! \brief Int-type array-of-struct attributes of a super-droplet */
using SuperDropletsIntIdxAoS = ERFParticlesIntIdxAoS;

/*! \brief Real-type struct-of-array attributes of a super-droplet */
using SuperDropletsRealIdxSoA = ERFParticlesRealIdxSoA;

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

/*! \brief Real-type struct-of-array attributes added during runtime for ice */
struct SDIceRealIdxSoA_RT
{
    enum {
        T_fz = 0,   /*!< freezing temperature */
        a,          /*!< equatorial radius */
        c,          /*!< porous radius */
        m_rime,     /*!< rime mass */
        n_mono,     /*!< number of primary ice crystals */
        ncomps
    };
};

/*! \brief Int-type struct-of-array attributes of a super-droplet */
using SuperDropletsIntIdxSoA = ERFParticlesIntIdxSoA;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int ridx_a(const int a_i, /*!< Species index */
                 const int a_num_a, /*!< Number of aerosols */
                 const int a_num_s  /*!< Number of species */)
{
    amrex::ignore_unused(a_num_a);
    amrex::ignore_unused(a_num_s);
    return    SuperDropletsRealIdxSoA_RT::ncomps
            + a_i;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int ridx_s(const int a_i, /*!< Species index */
                 const int a_num_a, /*!< Number of aerosols */
                 const int a_num_s  /*!< Number of species */)
{
    amrex::ignore_unused(a_num_s);
    return    SuperDropletsRealIdxSoA_RT::ncomps
            + a_num_a
            + a_i;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int idx_a(const int a_i, /*!< Species index */
                 const int a_num_a, /*!< Number of aerosols */
                 const int a_num_s  /*!< Number of species */)
{
    return SuperDropletsRealIdxSoA::ncomps + ridx_a(a_i,a_num_a,a_num_s);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int idx_s(const int a_i, /*!< Species index */
                 const int a_num_a, /*!< Number of aerosols */
                 const int a_num_s  /*!< Number of species */)
{
    return SuperDropletsRealIdxSoA::ncomps + ridx_s(a_i,a_num_a,a_num_s);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int ridx_ice_Tfz(const int a_num_a, /*!< Number of aerosols */
                        const int a_num_s  /*!< Number of species */)
{
    return   SuperDropletsRealIdxSoA_RT::ncomps
           + a_num_a
           + a_num_s
           + SDIceRealIdxSoA_RT::T_fz;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int ridx_ice_a(const int a_num_a, /*!< Number of aerosols */
                      const int a_num_s  /*!< Number of species */)
{
    return   SuperDropletsRealIdxSoA_RT::ncomps
           + a_num_a
           + a_num_s
           + SDIceRealIdxSoA_RT::a;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int ridx_ice_c(const int a_num_a, /*!< Number of aerosols */
                      const int a_num_s  /*!< Number of species */)
{
    return   SuperDropletsRealIdxSoA_RT::ncomps
           + a_num_a
           + a_num_s
           + SDIceRealIdxSoA_RT::c;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int ridx_ice_mrime(const int a_num_a, /*!< Number of aerosols */
                          const int a_num_s  /*!< Number of species */)
{
    return   SuperDropletsRealIdxSoA_RT::ncomps
           + a_num_a
           + a_num_s
           + SDIceRealIdxSoA_RT::m_rime;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int ridx_ice_nmono(const int a_num_a, /*!< Number of aerosols */
                          const int a_num_s  /*!< Number of species */)
{
    return   SuperDropletsRealIdxSoA_RT::ncomps
           + a_num_a
           + a_num_s
           + SDIceRealIdxSoA_RT::n_mono;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int idx_ice_Tfz(const int a_num_a, /*!< Number of aerosols */
                       const int a_num_s  /*!< Number of species */)
{
    return SuperDropletsRealIdxSoA::ncomps + ridx_ice_Tfz(a_num_a,a_num_s);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int idx_ice_a(const int a_num_a, /*!< Number of aerosols */
                     const int a_num_s  /*!< Number of species */)
{
    return SuperDropletsRealIdxSoA::ncomps + ridx_ice_a(a_num_a,a_num_s);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int idx_ice_c(const int a_num_a, /*!< Number of aerosols */
                     const int a_num_s  /*!< Number of species */)
{
    return SuperDropletsRealIdxSoA::ncomps + ridx_ice_c(a_num_a,a_num_s);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int idx_ice_mrime(const int a_num_a, /*!< Number of aerosols */
                         const int a_num_s  /*!< Number of species */)
{
    return SuperDropletsRealIdxSoA::ncomps + ridx_ice_mrime(a_num_a,a_num_s);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int idx_ice_nmono(const int a_num_a, /*!< Number of aerosols */
                         const int a_num_s  /*!< Number of species */)
{
    return SuperDropletsRealIdxSoA::ncomps + ridx_ice_nmono(a_num_a,a_num_s);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static auto ice_rho(const amrex::ParticleReal a_a, /*!< polar radius */
                    const amrex::ParticleReal a_c, /*!< equatorial radius */
                    const amrex::ParticleReal a_mass /*!< mass */ )
{
    auto vol = (4.0/3.0)*PI*a_a*a_a*a_c;
    if (vol == 0.0) {
        AMREX_ALWAYS_ASSERT(a_mass == 0.0);
        return 0.0;
    } else {
        return a_mass/vol;
    }
}

/*! \brief Compute total mass of a droplet from its components */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static amrex::ParticleReal SD_total_mass( const int a_idx, /*!< Particle index */
                           const int a_num_sp,   /*!< number of species */
                           const int a_num_ae,   /*!< number of aerosols */
                           const SDSpeciesMassArr& a_sp_mass, /*!< species masses */
                           const SDAerosolMassArr& a_ae_mass  /*!< aerosol masses */ )
{
    amrex::ParticleReal retval = 0.0;
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
                                 const SDSpeciesMassArr& a_sp_mass, /*!< species masses */
                                 const SDAerosolMassArr& a_ae_mass, /*!< aerosol masses */
                                 const amrex::ParticleReal* const a_sp_rho, /*!< species densities */
                                 const amrex::ParticleReal* const a_ae_rho  /*!< aerosol densities */ )
{
    amrex::ParticleReal m_w = a_sp_mass[a_idx_w][a_idx];
    amrex::ParticleReal m_s = 0.0;
    amrex::ParticleReal m_p = 0.0;
    amrex::ParticleReal rho_p = 0.0;
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
    if (m_p > 0.0) { rho_p /= m_p; }
    else           { rho_p = 1.0; }
    auto m_t = m_w + m_s + (a_rho_w/rho_p)*m_p;
    auto r_eff = std::cbrt(m_t / (4.0/3.0*PI*a_rho_w));
    return r_eff;
}

/*! \brief Compute "dry" radius of a droplet from its insoluble components */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static amrex::ParticleReal SD_dry_radius( const int a_idx, /*!< Particle index */
                           const int a_num_sp,   /*!< number of species */
                           const int a_num_ae,   /*!< number of aerosols */
                           const int* const a_sp_sol_arr, /*!< solubility of species */
                           const int* const a_ae_sol_arr, /*!< solubility of aerosols */
                           const SDSpeciesMassArr& a_sp_mass, /*!< species masses */
                           const SDAerosolMassArr& a_ae_mass, /*!< aerosol masses */
                           const amrex::ParticleReal* const a_sp_rho, /*!< species densities */
                           const amrex::ParticleReal* const a_ae_rho  /*!< aerosol densities */ )
{
    amrex::ParticleReal m_p = 0.0;
    amrex::ParticleReal rho_p = 0.0;
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
    if (m_p > 0.0) { rho_p /= m_p; }
    else           { rho_p = 1.0; }
    auto r_dry = std::cbrt(m_p / (4.0/3.0*PI*rho_p));
    return r_dry;
}

/* \brief SD phase (ice or water) */
AMREX_ENUM (SDPhase,
    ice,
    water
);

/*! \brief Is the SD ice or water? */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static auto SD_phase( const int a_p, /*!< Particle index */
                      const int a_i_w, /*!< Species index of water */
                      const int a_i_i, /*!< Species index of ice */
                      const SDSpeciesMassArr& a_sp_mass  /*!< species masses */)
{
    if (a_i_i < 0) {
        return SDPhase::water; // not modeling ice
    } else {
        auto m_water = a_sp_mass[a_i_w][a_p];
        auto m_ice   = a_sp_mass[a_i_i][a_p];
        if (m_ice > 0.0) {
            AMREX_ALWAYS_ASSERT(m_water == 0.0); // must be pure ice
            return SDPhase::ice;
        } else {
            AMREX_ALWAYS_ASSERT(m_ice == 0.0); // must be pure water
            return SDPhase::water;
        }
    }
}

/*! \brief Types of terminal velocity models */
AMREX_ENUM(SDTerminalVelocityType,
    RogersYau,
    AtlasUlbrich,
    CloudRainShima,
    IceBohm
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

