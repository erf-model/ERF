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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int ridx_ice_Tfz(const int a_num_a, /*!< Number of aerosols */
                        const int a_num_s  /*!< Number of species */)
{
    return   SDPCDefn::SuperDropletsRealIdxSoA_RT::ncomps
           + a_num_a
           + a_num_s
           + SDPCDefn::SDIceRealIdxSoA_RT::T_fz;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int ridx_ice_a(const int a_num_a, /*!< Number of aerosols */
                      const int a_num_s  /*!< Number of species */)
{
    return   SDPCDefn::SuperDropletsRealIdxSoA_RT::ncomps
           + a_num_a
           + a_num_s
           + SDPCDefn::SDIceRealIdxSoA_RT::a;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int ridx_ice_c(const int a_num_a, /*!< Number of aerosols */
                      const int a_num_s  /*!< Number of species */)
{
    return   SDPCDefn::SuperDropletsRealIdxSoA_RT::ncomps
           + a_num_a
           + a_num_s
           + SDPCDefn::SDIceRealIdxSoA_RT::c;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int ridx_ice_mrime(const int a_num_a, /*!< Number of aerosols */
                          const int a_num_s  /*!< Number of species */)
{
    return   SDPCDefn::SuperDropletsRealIdxSoA_RT::ncomps
           + a_num_a
           + a_num_s
           + SDPCDefn::SDIceRealIdxSoA_RT::m_rime;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int ridx_ice_nmono(const int a_num_a, /*!< Number of aerosols */
                          const int a_num_s  /*!< Number of species */)
{
    return   SDPCDefn::SuperDropletsRealIdxSoA_RT::ncomps
           + a_num_a
           + a_num_s
           + SDPCDefn::SDIceRealIdxSoA_RT::n_mono;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int idx_ice_Tfz(const int a_num_a, /*!< Number of aerosols */
                       const int a_num_s  /*!< Number of species */)
{
    return SDPCDefn::SuperDropletsRealIdx::ncomps + ridx_ice_Tfz(a_num_a,a_num_s);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int idx_ice_a(const int a_num_a, /*!< Number of aerosols */
                     const int a_num_s  /*!< Number of species */)
{
    return SDPCDefn::SuperDropletsRealIdx::ncomps + ridx_ice_a(a_num_a,a_num_s);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int idx_ice_c(const int a_num_a, /*!< Number of aerosols */
                     const int a_num_s  /*!< Number of species */)
{
    return SDPCDefn::SuperDropletsRealIdx::ncomps + ridx_ice_c(a_num_a,a_num_s);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int idx_ice_mrime(const int a_num_a, /*!< Number of aerosols */
                         const int a_num_s  /*!< Number of species */)
{
    return SDPCDefn::SuperDropletsRealIdx::ncomps + ridx_ice_mrime(a_num_a,a_num_s);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int idx_ice_nmono(const int a_num_a, /*!< Number of aerosols */
                         const int a_num_s  /*!< Number of species */)
{
    return SDPCDefn::SuperDropletsRealIdx::ncomps + ridx_ice_nmono(a_num_a,a_num_s);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static auto ice_rho(const amrex::ParticleReal a_a, /*!< polar radius */
                    const amrex::ParticleReal a_c, /*!< equatorial radius */
                    const amrex::ParticleReal a_mass /*!< mass */ )
{
    auto vol = static_cast<amrex::ParticleReal>(four_thirds_pi)*a_a*a_a*a_c;
    if (vol == amrex::ParticleReal(0)) {
        AMREX_ALWAYS_ASSERT(a_mass == amrex::ParticleReal(0));
        return amrex::ParticleReal(0);
    } else {
        return a_mass/vol;
    }
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

/*! \brief Compute total surface area of INP (ice-nucleating) components
 *
 *  The following code is adapted from SCALE-SDM:
 *  https://github.com/Shima-Lab/SCALE-SDM_mixed-phase_Shima2019
 *  Copyright (c) 2012-2015, Team SCALE
 *  All rights reserved.
 *  BSD 2-Clause License
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static amrex::ParticleReal SD_INP_surface_area(
    const int a_idx,                              /*!< Particle index */
    const int a_num_sp,                           /*!< number of species */
    const int a_num_ae,                           /*!< number of aerosols */
    const int* const a_sp_is_INP,                 /*!< INP flag for species */
    const int* const a_ae_is_INP,                 /*!< INP flag for aerosols */
    const SDPCDefn::SDSpeciesMassArr& a_sp_mass,  /*!< species masses */
    const SDPCDefn::SDAerosolMassArr& a_ae_mass,  /*!< aerosol masses */
    const amrex::ParticleReal* const a_sp_rho,    /*!< species densities */
    const amrex::ParticleReal* const a_ae_rho     /*!< aerosol densities */ )
{
    amrex::ParticleReal total_area = 0.0;

    // Species contribution
    for (int j = 0; j < a_num_sp; j++) {
        if (a_sp_is_INP[j] && a_sp_mass[j][a_idx] > 0.0) {
            auto vol = a_sp_mass[j][a_idx] / a_sp_rho[j];
            auto r = std::cbrt(vol * 3.0 / (4.0 * PI));
            total_area += 4.0 * PI * r * r;
        }
    }

    // Aerosol contribution
    for (int j = 0; j < a_num_ae; j++) {
        if (a_ae_is_INP[j] && a_ae_mass[j][a_idx] > 0.0) {
            auto vol = a_ae_mass[j][a_idx] / a_ae_rho[j];
            auto r = std::cbrt(vol * 3.0 / (4.0 * PI));
            total_area += 4.0 * PI * r * r;
        }
    }

    return total_area;
}

/*! \brief INAS density parameterization for immersion freezing
 *
 *  Implements the Ice Nucleation Active Site (INAS) density parameterization
 *  of Niemand et al. (2012) for heterogeneous ice nucleation.
 *
 *  Reference:
 *  - Niemand, M., et al. (2012): A Particle-Surface-Area-Based Parameterization
 *    of Immersion Freezing on Desert Dust Particles. J. Atmos. Sci., 69, 3077-3092.
 *
 *  The following code is adapted from SCALE-SDM:
 *  https://github.com/Shima-Lab/SCALE-SDM_mixed-phase_Shima2019
 *  Copyright (c) 2012-2015, Team SCALE
 *  All rights reserved.
 *  BSD 2-Clause License
 */
struct INAS_Niemand2012 {
    amrex::Real a0 = 1.0;          /*!< INAS a0 parameter [m^-2] */
    amrex::Real a1 = 0.517;        /*!< INAS a1 parameter [degC^-1] */
    amrex::Real a2 = 8.934;        /*!< INAS a2 parameter [-] */

    amrex::Real Tfz_min = -36.0;   /*!< coldest heterogeneous freezing [degC] */
    amrex::Real Tfz_max = -12.0;   /*!< warmest heterogeneous freezing [degC] */
    amrex::Real Tfz_homogeneous = -38.0; /*!< homogeneous freezing limit [degC] */

    /*! \brief Compute INAS density at temperature T
     *  n_s(T) = a0 * exp(-a1*T + a2)  [active sites per m^2]
     *  where T is in degrees Celsius */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real n_s(amrex::Real T_degC) const {
        return a0 * std::exp(-a1 * T_degC + a2);
    }

    /*! \brief Compute cumulative freezing probability at temperature T
     *  P(T) = 1 - exp(-A * n_s(T))
     *  where A is the INP surface area */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real P_freeze(amrex::Real T_degC, amrex::Real A) const {
        return 1.0 - std::exp(-A * n_s(T_degC));
    }

    /*! \brief Sample freezing temperature using inverse CDF method
     *
     *  Given uniform random number u in [0,1] and INP surface area A,
     *  returns Tfz in Kelvin.
     *
     *  Based on the singular hypothesis: each particle has a characteristic
     *  freezing temperature determined by its coldest active site.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real sample_Tfz(amrex::Real u, amrex::Real A) const {
        // No INP surface area -> homogeneous freezing
        if (A <= 0.0) {
            return Tfz_homogeneous + 273.15;
        }

        // Check if particle is IN-inactive (no active sites even at Tfz_min)
        // This corresponds to particles that won't nucleate ice heterogeneously
        auto P_max = P_freeze(Tfz_min, A);
        if (u > P_max) {
            // Particle has no active sites -> homogeneous freezing
            return Tfz_homogeneous + amrex::Real(273.15);
        }

        // Inverse CDF: solve P(T) = u for T
        // P(T) = 1 - exp(-A * a0 * exp(-a1*T + a2))
        // Let y = -ln(1-u) = A * a0 * exp(-a1*T + a2)
        // -a1*T + a2 = ln(y / (A*a0))
        // T = (a2 - ln(y/(A*a0))) / a1
        auto y = -std::log(amrex::Real(1) - u);
        auto T_degC = (a2 - std::log(y / (A * a0))) / a1;

        // Clamp to valid range
        T_degC = std::max(T_degC, Tfz_min);
        T_degC = std::min(T_degC, Tfz_max);

        return T_degC + amrex::Real(273.15);  // Convert to Kelvin
    }
};

/* \brief SD phase (pure ice, pure water, or mixed ice-water) */
AMREX_ENUM (SDPhase,
    ice,
    water,
    mixed
);

/*! \brief Liquid mass fraction of an SD: m_water/(m_water+m_ice).
 *  Returns 1 for pure water and 0 for pure ice. */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static auto SD_liquid_fraction( const int a_p, /*!< Particle index */
                                const int a_i_w, /*!< Species index of water */
                                const int a_i_i, /*!< Species index of ice */
                                const SDPCDefn::SDSpeciesMassArr& a_sp_mass  /*!< species masses */)
{
    if (a_i_i < 0) { return amrex::ParticleReal(1); } // not modeling ice
    auto m_water = a_sp_mass[a_i_w][a_p];
    auto m_ice   = a_sp_mass[a_i_i][a_p];
    auto m_cond  = m_water + m_ice;
    return (m_cond > amrex::ParticleReal(0)) ? (m_water/m_cond) : amrex::ParticleReal(1);
}

/*! \brief Is the SD pure ice, pure water, or a mixed ice-water particle? */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static auto SD_phase( const int a_p, /*!< Particle index */
                      const int a_i_w, /*!< Species index of water */
                      const int a_i_i, /*!< Species index of ice */
                      const SDPCDefn::SDSpeciesMassArr& a_sp_mass  /*!< species masses */)
{
    if (a_i_i < 0) {
        return SDPhase::water; // not modeling ice
    } else {
        auto m_water = a_sp_mass[a_i_w][a_p];
        auto m_ice   = a_sp_mass[a_i_i][a_p];
        if (m_ice > 0.0) {
            return (m_water > 0.0) ? SDPhase::mixed : SDPhase::ice;
        } else {
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

/*! \brief Types of kernel relative velocity calculation: terminal_velocity
 *  (only terminal fall speed; default), absolute_velocity (flow velocity minus
 *  terminal), radial_velocity (radial component of the relative velocity along
 *  the pair separation vector). */
AMREX_ENUM(SDKernelRelativeVelocityType,
    terminal_velocity,
    absolute_velocity,
    radial_velocity
);

/*! \brief Types of phase change ODE time integrator */
enum struct SDMassChangeTIMethod {
    RK3BS, RK4, BE, CN, DIRK2
};

#endif
#endif

