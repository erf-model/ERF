#include "ERF_Constants.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_SuperDropletPC.H"
#include "ERF_SuperDropletPCCoalescence.H"
#include "ERF_SuperDropletPCMassChange.H"
#include <AMReX_TracerParticle_mod_K.H>
#include "ERF_InterpolationUtils.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;
using namespace SDMassChangeUtils_SV;
using namespace SDPCDefn;

namespace {
    /*! \brief Update common attributes for predator particle (remainder > 0 case)
     *  Updates Tfz, species masses, and aerosol masses for particle i */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void update_common_positive_rmndr(
        const int a_i, const int a_j,
        const ParticleReal a_gamma,
        ParticleReal* const a_Tfz,
        const int a_n_sp, const SDSpeciesMassArr& a_sp_m,
        const int a_n_ae, const SDAerosolMassArr& a_ae_m)
    {
        a_Tfz[a_i] = std::max(a_Tfz[a_j], a_Tfz[a_i]);
        for (int n = 0; n < a_n_sp; n++) { a_sp_m[n][a_i] += a_gamma * a_sp_m[n][a_j]; }
        for (int n = 0; n < a_n_ae; n++) { a_ae_m[n][a_i] += a_gamma * a_ae_m[n][a_j]; }
    }

    /*! \brief Update common attributes for both particles (remainder == 0 case)
     *  Splits multiplicity and updates Tfz, species masses, aerosol masses for both */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void update_common_zero_rmndr(
        const int a_i, const int a_j,
        const ParticleReal a_gamma,
        ParticleReal* const a_mult,
        ParticleReal* const a_Tfz,
        const int a_n_sp, const SDSpeciesMassArr& a_sp_m,
        const int a_n_ae, const SDAerosolMassArr& a_ae_m)
    {
        // Split multiplicity
        ParticleReal dm = std::floor(a_mult[a_j]/2);
        a_mult[a_i] = dm;
        a_mult[a_j] -= dm;
        // Update Tfz
        a_Tfz[a_i] = a_Tfz[a_j] = std::max(a_Tfz[a_j], a_Tfz[a_i]);
        // Update species masses
        for (int n = 0; n < a_n_sp; n++) {
            a_sp_m[n][a_j] += a_gamma * a_sp_m[n][a_i];
            a_sp_m[n][a_i] = a_sp_m[n][a_j];
        }
        // Update aerosol masses
        for (int n = 0; n < a_n_ae; n++) {
            a_ae_m[n][a_j] += a_gamma * a_ae_m[n][a_i];
            a_ae_m[n][a_i] = a_ae_m[n][a_j];
        }
    }
}

/*! \brief Compute dynamic viscosity */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static auto viscCoeff ( const ParticleReal a_T /*!< temperature */ )
{
    auto T_degC = a_T - ParticleReal(tmelt); // [K] => [degC]
    ParticleReal visc_coeff = ParticleReal(zero);
    if( T_degC >= ParticleReal(zero) ) {
        visc_coeff = ( ParticleReal(1.7180) + ParticleReal(4.9E-3)*T_degC ) * ParticleReal(1.E-5);
    } else {
        visc_coeff = ( ParticleReal(1.7180) + ParticleReal(4.9E-3)*T_degC - ParticleReal(1.2E-5)*T_degC*T_degC ) * ParticleReal(1.E-5);
    }
    return visc_coeff;
}

/*! \brief Compute coalescence rate between two superdroplets */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static ParticleReal coalescence_rate ( const RandomEngine& a_rnd_eng, /*!< random engine */
                                       const Real a_p /*!< probability */ )
{
    ParticleReal p_int = std::floor(a_p);
    ParticleReal gamma = p_int;
    if (Random(a_rnd_eng) < (a_p-p_int)) { gamma += 1; }
    return gamma;
}

/*! \brief Binary coalescence between two superdroplets
 *
 *  The following Tfz update logic is adapted from SCALE-SDM:
 *  https://github.com/Shima-Lab/SCALE-SDM_mixed-phase_Shima2019
 *  Copyright (c) 2012-2015, Team SCALE
 *  All rights reserved.
 *  BSD 2-Clause License
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static void coal_update_attribs(const int a_i, /*!< index of particle */
                                const int a_j, /*!< index of coalescence partner */
                                const int* const a_prey, /*!< prey/predator */
                                const ParticleReal* const a_gamma, /*!< coalescence rate */
                                const ParticleReal* const a_rmndr, /*!< coalescence remainder*/
                                ParticleReal* const a_mass, /*!< mass */
                                ParticleReal* const a_radius, /*!< radius */
                                ParticleReal* const a_Tfz, /*!< freezing temperature */
                                ParticleReal* const a_mult, /*!< multiplicity */
                                const int a_n_sp, /*!< number of species */
                                const SDSpeciesMassArr& a_sp_m, /*!< species masses*/
                                const int a_n_ae, /*!< number of aerosols */
                                const SDAerosolMassArr& a_ae_m /*!< aerosol masses*/)
{
    int i = a_i;
    int j = a_j;
    auto gamma = a_gamma[i];
    AMREX_ALWAYS_ASSERT(gamma == a_gamma[j]);
    AMREX_ALWAYS_ASSERT(a_rmndr[a_i] >= ParticleReal(zero));

    if ( a_rmndr[a_i] > 0 ) {

        if (a_prey[i]) {
            a_mult[i] -= gamma * a_mult[j];
        } else {
            auto r3 = gamma*a_radius[j]*a_radius[j]*a_radius[j]
                          + a_radius[i]*a_radius[i]*a_radius[i];
            a_radius[i] = std::cbrt(r3);
            a_mass[i] += gamma * a_mass[j];
            update_common_positive_rmndr(i, j, gamma, a_Tfz, a_n_sp, a_sp_m, a_n_ae, a_ae_m);
        }

    } else if ( a_rmndr[a_i] == 0 ) {

        if (a_prey[i]) {
            auto r3 = gamma*a_radius[i]*a_radius[i]*a_radius[i]
                          + a_radius[j]*a_radius[j]*a_radius[j];
            a_radius[i] = a_radius[j] = std::cbrt(r3);
            a_mass[j] += gamma * a_mass[i];
            a_mass[i] = a_mass[j];
            update_common_zero_rmndr(i, j, gamma, a_mult, a_Tfz, a_n_sp, a_sp_m, a_n_ae, a_ae_m);
        }

    }

}

/*! \brief Binary aggregation between two superdroplets */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static void aggr_update_attribs(const int a_i, /*!< index of particle */
                                const int a_j, /*!< index of coalescence partner */
                                const int a_sp_idx_i, /*!< species index of ice */
                                const Real a_rho_ice, /*!< true ice density */
                                const Real a_rho_min, /*!< minimum ice density */
                                const int* const a_prey, /*!< prey/predator */
                                const ParticleReal* const a_gamma, /*!< coalescence rate */
                                const ParticleReal* const a_rmndr, /*!< coalescence remainder*/
                                ParticleReal* const a_Tfz, /*!< freezing temperature */
                                ParticleReal* const a_a, /*!< equatorial radius */
                                ParticleReal* const a_c, /*!< polar radius */
                                ParticleReal* const a_mrime, /*!< rime mass */
                                ParticleReal* const a_nmono, /*!< number of monomers */
                                ParticleReal* const a_mult, /*!< multiplicity */
                                const int a_n_sp, /*!< number of species */
                                const SDSpeciesMassArr& a_sp_m, /*!< species masses*/
                                const int a_n_ae, /*!< number of aerosols */
                                const SDAerosolMassArr& a_ae_m /*!< aerosol masses*/)
{
    AMREX_ALWAYS_ASSERT(a_gamma[a_i] == a_gamma[a_j]);
    AMREX_ALWAYS_ASSERT(a_rmndr[a_i] >= ParticleReal(zero));

    auto gamma = a_gamma[a_i];
    ParticleReal a_new = ParticleReal(zero);
    ParticleReal c_new = ParticleReal(zero);
    ParticleReal m_new = ParticleReal(zero);

    {
        // i will be the prey particle index, j is the partner
        int i = -1, j = -1;
        if (a_prey[a_i]) {
            AMREX_ALWAYS_ASSERT(!a_prey[a_j]);
            i = a_i;
            j = a_j;
        } else {
            AMREX_ALWAYS_ASSERT(a_prey[a_j]);
            i = a_j;
            j = a_i;
        }

        m_new = gamma * a_sp_m[a_sp_idx_i][i] + a_sp_m[a_sp_idx_i][j];
        auto V_new_min = ParticleReal(four_thirds_pi) * (gamma*a_a[i]*a_a[i]*a_c[i] + a_a[j]*a_a[j]*a_c[j]);
        auto rhoi_bar = m_new / V_new_min;

        auto maxR_i = std::max(a_a[i], a_c[i]);
        auto maxR_j = std::max(a_a[j], a_c[j]);

        if (maxR_i > maxR_j) {
            if (a_a[i] > a_c[i]) {
                a_new = a_a[i];
                auto c_new_min = V_new_min / (ParticleReal(four_thirds_pi)*a_new*a_new);
                auto c_new_max = a_c[i]*gamma + std::min(a_a[j],a_c[j]);
                c_new = (a_rho_ice-a_rho_min) / ((a_rho_ice-rhoi_bar)/c_new_min + (rhoi_bar-a_rho_min)/c_new_max);
            } else {
                c_new = a_c[i];
                auto a_new_min = std::sqrt(V_new_min / (ParticleReal(four_thirds_pi)*c_new));
                auto a_new_max = std::sqrt( std::max(std::max(a_a[i],a_a[j]),a_c[j])
                                            * (a_a[i]*gamma + std::min(a_a[j],a_c[j])) );
                a_new = std::sqrt( (a_rho_ice - a_rho_min)
                                    / (   (a_rho_ice-rhoi_bar)/(a_new_min*a_new_min)
                                        + (rhoi_bar-a_rho_min)/(a_new_max*a_new_max)) );
            }
        } else {
            if (a_a[j] > a_c[j]) {
                a_new = a_a[j];
                auto c_new_min = V_new_min / (ParticleReal(four_thirds_pi)*a_new*a_new);
                auto c_new_max = a_c[j] + std::min(a_a[i],a_c[i])*gamma;
                c_new = (a_rho_ice-a_rho_min) / ((a_rho_ice-rhoi_bar)/c_new_min + (rhoi_bar-a_rho_min)/c_new_max);
            } else {
                c_new = a_c[j];
                auto a_new_min = std::sqrt(V_new_min / (ParticleReal(four_thirds_pi)*c_new));
                auto a_new_max = std::sqrt( std::max(std::max(a_a[j],a_a[i]),a_c[i])
                                            * (a_a[j] + gamma*std::min(a_a[i],a_c[i])) );
                a_new = std::sqrt( (a_rho_ice - a_rho_min)
                                    / (   (a_rho_ice-rhoi_bar)/(a_new_min*a_new_min)
                                        + (rhoi_bar-a_rho_min)/(a_new_max*a_new_max)) );
            }
        }
    }

    if ( a_rmndr[a_i] > 0 ) {
        if (a_prey[a_i]) {
            a_mult[a_i] -= gamma * a_mult[a_j];
        } else {
            a_a[a_i] = a_new;
            a_c[a_i] = c_new;
            a_mrime[a_i] += gamma * a_mrime[a_j];
            a_nmono[a_i] += gamma * a_nmono[a_j];
            update_common_positive_rmndr(a_i, a_j, gamma, a_Tfz, a_n_sp, a_sp_m, a_n_ae, a_ae_m);
            AMREX_ALWAYS_ASSERT(m_new == a_sp_m[a_sp_idx_i][a_i]);
        }
    } else {
        if (a_prey[a_i]) {
            a_a[a_i] = a_a[a_j] = a_new;
            a_c[a_i] = a_c[a_j] = c_new;
            a_mrime[a_j] += gamma * a_mrime[a_i]; a_mrime[a_i] = a_mrime[a_j];
            a_nmono[a_j] += gamma * a_nmono[a_i]; a_nmono[a_i] = a_nmono[a_j];
            update_common_zero_rmndr(a_i, a_j, gamma, a_mult, a_Tfz, a_n_sp, a_sp_m, a_n_ae, a_ae_m);
            AMREX_ALWAYS_ASSERT(m_new == a_sp_m[a_sp_idx_i][a_i]);
        }
    }
}

/*! \brief Impact velocity ratio from Rasmussen and Heymsfield (1985)
 *
 *  Computes the ratio of impact velocity to relative velocity for ice-droplet
 *  collisions as a function of Reynolds and Stokes numbers.
 *
 *  Reference: Rasmussen, R. M., and A. J. Heymsfield, 1985: A generalized
 *  form for impact velocities used to determine graupel accretional densities.
 *  J. Atmos. Sci., 42, 2275-2279.
 *  https://doi.org/10.1175/1520-0469(1985)042<2275:AGFFIV>2.0.CO;2
 *
 *  Polynomial fit: v_impact/v_rel = A0 + A1*w + A2*w^2 + A3*w^3 + A4*w^4
 *  where w = log10(St)
 *
 *  Coefficients from Table 1 for different Re ranges:
 *  - Re < 20:    (0.1701, 0.7246, 0.2257, -1.13, 0.5756), asymptote 0.57
 *  - 20 ≤ Re < 65:  (0.2927, 0.5085, -0.03453, -0.2184, 0.03595), asymptote 0.59
 *  - 65 ≤ Re < 200: (0.3272, 0.4907, -0.09452, -0.1906, 0.07105), asymptote 0.61
 *  - Re ≥ 200:   (0.356, 0.4738, -0.1233, -0.1618, 0.08087), asymptote 0.63
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static auto impactVelocity_RasmussenHeymsfield1985( const ParticleReal a_Re, /*!< Reynolds number */
                                                    const ParticleReal a_St  /*!< Stokes number */ )
{
    auto w = std::log10(a_St);
    auto w2 = w*w;
    auto w3 = w2*w;
    auto w4 = w3*w;
    ParticleReal retval = ParticleReal(zero);
    if(a_Re < (ParticleReal(10.0)+ParticleReal(30.0))*ParticleReal(myhalf)) {
        if (a_St < ParticleReal(0.4)) {
            retval = ParticleReal(zero);
        } else if (a_St<ParticleReal(10.0)) {
            retval = ParticleReal(0.1701) + ParticleReal(0.7246)*w + ParticleReal(0.2257)*w2 - ParticleReal(1.13)*w3 + ParticleReal(0.5756)*w4;
        } else {
            retval = ParticleReal(0.57);
        }
    } else if (a_Re < (ParticleReal(30.0)+ParticleReal(100.0))*ParticleReal(myhalf)) {
        if (a_St < ParticleReal(0.1)) {
            retval = ParticleReal(zero);
        } else if (a_St < ParticleReal(10.0)) {
            retval = ParticleReal(0.2927) + ParticleReal(0.5085)*w - ParticleReal(0.03453)*w2 - ParticleReal(0.2184)*w3 + ParticleReal(0.03595)*w4;
        } else {
            retval = ParticleReal(0.59);
        }
    } else if (a_Re < (ParticleReal(100.0)+ParticleReal(300.0))*ParticleReal(myhalf)) {
        if (a_St < ParticleReal(0.1)) {
          retval = ParticleReal(zero);
        } else if (a_St < ParticleReal(10.0)) {
            retval = ParticleReal(0.3272) + ParticleReal(0.4907)*w - ParticleReal(0.09452)*w2 - ParticleReal(0.1906)*w3 + ParticleReal(0.07105)*w4;
        } else {
            retval = ParticleReal(0.61);
        }
    } else {
        if (a_St < ParticleReal(0.1)) {
            retval = ParticleReal(zero);
        } else if (a_St < ParticleReal(10.0)) {
            retval = ParticleReal(0.356) + ParticleReal(0.4738)*w - ParticleReal(0.1233)*w2 - ParticleReal(0.1618)*w3 + ParticleReal(0.08087)*w4;
        } else {
            retval = ParticleReal(0.63);
        }
    }
    retval = std::max(retval,ParticleReal(zero));
    return retval;
}

/*! \brief Ice surface temperature */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static auto iceSurfaceTemperature( const ParticleReal a_T, /*!< temperature */
                                   const ParticleReal a_P, /*!< pressure */
                                   const ParticleReal a_qv, /*!< vapour fraction */
                                   const ParticleReal a_D, /*!< diffusivity coeff */
                                   const dMdt<ParticleReal>& a_dmdt /*!< mass change utilities */)
{
    Real qsat = Real(zero); erf_qsati(a_T, a_P, qsat);
    auto sup_sat  = a_qv/qsat - Real(one);
    auto drho  = sup_sat / (a_D * a_dmdt.Fk_plus_Fd(a_T, erf_esati(a_T), a_D));
    auto dT = (a_dmdt.L * a_D / a_dmdt.K) * drho;
    return a_T + dT;
}

/*! \brief Rime density parameterization from Heymsfield and Pflaum (1985)
 *
 *  Computes the density of rime accreted onto ice particles based on the
 *  dimensionless impact parameter Y = -r[μm] * v_impact / T_surf[°C].
 *
 *  Reference: Heymsfield, A. J., and J. C. Pflaum, 1985: A quantitative
 *  assessment of the accuracy of techniques for calculating graupel growth.
 *  J. Atmos. Sci., 42, 2264-2274.
 *  https://doi.org/10.1175/1520-0469(1985)042<2264:AQAOTA>2.0.CO;2
 *
 *  Two regimes based on surface temperature and Y parameter:
 *  - Low-density (T_surf ≤ -5°C or Y < 1.6): ρ = (0.30*Y)^0.44  [Eq. 8]
 *  - High-density (otherwise): ρ = exp(-0.03115 - 1.703*Y + 0.9116*Y² - 0.1224*Y³) [Eq. 9]
 *
 *  Output bounded to 0.1 ≤ ρ_rime ≤ 0.91 g/cm³ (100-910 kg/m³)
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static auto rimeDensity_HeymsfieldPflaum1985( const ParticleReal a_radius,      /*!< droplet radius [m] */
                                              const ParticleReal a_a,           /*!< ice particle equatorial radius [m] */
                                              const ParticleReal a_c,           /*!< ice particle polar radius [m] */
                                              const ParticleReal a_vz_w,        /*!< vertical velocity of water droplet [m/s] */
                                              const ParticleReal a_vz_i,        /*!< vertical velocity of ice particle [m/s] */
                                              const ParticleReal a_rho_w,       /*!< water density [kg/m³] */
                                              const ParticleReal a_D,           /*!< diffusivity coefficient [m²/s] */
                                              const dMdt<ParticleReal>& a_dmdt, /*!< mass change utilities */
                                              const ParticleReal a_T,           /*!< temperature [K] */
                                              const ParticleReal a_rhom,        /*!< moist air density [kg/m³] */
                                              const ParticleReal a_P,           /*!< pressure [Pa] */
                                              const ParticleReal a_qv           /*!< water vapor mixing ratio [kg/kg] */)
{
    auto r_um = a_radius * ParticleReal(1.0e6);
    auto mu = viscCoeff(a_T);

    auto maxD = ParticleReal(two) * std::max(a_a,a_c);
    auto Re = a_rhom * maxD * std::abs(a_vz_i) / mu;
    auto eqr_i = std::cbrt(a_a*a_a*a_c);
    auto St = ParticleReal(two)*std::abs(a_vz_i)*a_radius*a_radius*a_rho_w/(ParticleReal(9.0)*mu*eqr_i);

    auto v_impact_ratio = impactVelocity_RasmussenHeymsfield1985(Re,St);
    auto v_impact = std::abs(a_vz_i - a_vz_w) * v_impact_ratio;

    auto T_surf = std::min(ParticleReal(-0.01), iceSurfaceTemperature(a_T,a_P,a_qv,a_D,a_dmdt) - ParticleReal(tmelt));
    auto var_Y = -r_um * v_impact/T_surf;

    ParticleReal rho_rime = ParticleReal(zero);
    if ((T_surf <= ParticleReal(-5.0)) || (var_Y < ParticleReal(1.6))) {
        rho_rime = ParticleReal(0.30) * std::pow(var_Y, ParticleReal(0.44));
    } else {
        var_Y = std::min(var_Y,ParticleReal(3.5));
        rho_rime = std::exp(ParticleReal(-0.03115) - ParticleReal(1.7030)*var_Y + ParticleReal(0.9116)*var_Y*var_Y - ParticleReal(0.1224)*var_Y*var_Y*var_Y);
    }
    rho_rime = std::min(ParticleReal(0.91),std::max(rho_rime,ParticleReal(0.1))) * ParticleReal(1000.0);
    return rho_rime;
}

/*! \brief Binary riming between a water droplet and an ice particle */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static void rime_update_attribs(const int a_i, /*!< index of particle */
                                const int a_j, /*!< index of coalescence partner */
                                const int a_sp_idx_w, /*!< species index of water */
                                const int a_sp_idx_i, /*!< species index of ice */
                                const Real a_rho_water, /*!< water density */
                                const Real a_rho_ice, /*!< true ice density */
                                const int* const a_prey, /*!< prey/predator */
                                const ParticleReal* const a_gamma, /*!< coalescence rate */
                                const ParticleReal* const a_rmndr, /*!< coalescence remainder*/
                                const SDPhase a_phase_i, /*!< phase of particle i */
                                const SDPhase a_phase_j, /*!< phase of particle j (partner) */
                                const GpuArray<ParticleReal*,AMREX_SPACEDIM>& a_vel, /*!< velocity */
                                ParticleReal* const a_vterm, /*!< terminal velocity */
                                ParticleReal* const a_radius, /*!< radius */
                                ParticleReal* const a_Tfz, /*!< freezing temperature */
                                ParticleReal* const a_a, /*!< equatorial radius */
                                ParticleReal* const a_c, /*!< polar radius */
                                ParticleReal* const a_mrime, /*!< rime mass */
                                ParticleReal* const a_nmono, /*!< number of monomers */
                                ParticleReal* const a_mult, /*!< multiplicity */
                                const ParticleReal a_T, /*!< temperature */
                                const ParticleReal a_rhom, /*!< moist density */
                                const ParticleReal a_P, /*!< pressure */
                                const ParticleReal a_qv, /*!< qv */
                                const ParticleReal a_D, /*!< diffusivity coeff */
                                const dMdt<ParticleReal>& a_dmdt, /*!< mass change utilities */
                                const int a_n_sp, /*!< number of species */
                                const SDSpeciesMassArr& a_sp_m, /*!< species masses*/
                                const int a_n_ae, /*!< number of aerosols */
                                const SDAerosolMassArr& a_ae_m /*!< aerosol masses*/)
{
    AMREX_ALWAYS_ASSERT(a_gamma[a_i] == a_gamma[a_j]);
    AMREX_ALWAYS_ASSERT(a_rmndr[a_i] >= ParticleReal(0.0));
    auto gamma = a_gamma[a_i];

    // Identify ice and water particle indices
    int id_ice, id_water;
    if (a_phase_i == SDPhase::ice) {
        AMREX_ALWAYS_ASSERT(a_phase_j == SDPhase::water);
        id_ice = a_i;
        id_water = a_j;
    } else {
        AMREX_ALWAYS_ASSERT(a_phase_j == SDPhase::ice);
        id_ice = a_j;
        id_water = a_i;
    }

    // Determine gamma multipliers based on prey/predator status
    // gamma copies of prey merge with 1 predator
    ParticleReal gamma_ice, gamma_water;
    if (a_prey[id_ice]) {
        // ice is prey (higher multiplicity)
        gamma_ice = gamma;
        gamma_water = ParticleReal(one);
    } else {
        // water is prey (higher multiplicity)
        gamma_ice = ParticleReal(one);
        gamma_water = gamma;
    }

    // Compute new ice mass
    ParticleReal mi_new = gamma_ice * a_sp_m[a_sp_idx_i][id_ice]
                        + gamma_water * a_sp_m[a_sp_idx_w][id_water];

    // Compute new rime mass
    ParticleReal mrime_new = gamma_ice * a_mrime[id_ice]
                           + gamma_water * a_sp_m[a_sp_idx_w][id_water];

    // Compute new number of monomers
    ParticleReal nmono_new = gamma_ice * a_nmono[id_ice];

    // Compute new dimensions
    ParticleReal a_new = ParticleReal(zero);
    ParticleReal c_new = ParticleReal(zero);

    if (a_radius[id_water] > std::max(a_a[id_ice], a_c[id_ice])) {
        // Water droplet is larger than ice particle - form spherical ice
        a_new = c_new = std::cbrt(mi_new / (ParticleReal(four_thirds_pi)*a_rho_ice));
    } else {
        // Normal riming - ice particle collects water droplet
        // Pass terminal velocities for Re/St calculation (velocity relative to air)
        auto rho_rime = rimeDensity_HeymsfieldPflaum1985( a_radius[id_water],
                                                          a_a[id_ice], a_c[id_ice],
                                                          a_vterm[id_water], a_vterm[id_ice],
                                                          a_rho_water,
                                                          a_D, a_dmdt,
                                                          a_T, a_rhom, a_P, a_qv );
        auto V_new = ParticleReal(four_thirds_pi) * ( gamma_ice * a_a[id_ice]*a_a[id_ice]*a_c[id_ice]
                                    + gamma_water * a_radius[id_water]*a_radius[id_water]*a_radius[id_water]
                                                  * (a_rho_water/rho_rime) );
        auto phi = a_c[id_ice] / a_a[id_ice];
        // Following Shima et al. (2019) sdm_outcome_riming:
        // - When ice is prey: rime dimension doesn't include gamma
        // - When water is prey: rime dimension scales with gamma^(1/3)
        if ((phi < ParticleReal(0.8)) || ((phi < ParticleReal(1.25)) && (phi >= ParticleReal(one)))) {
            // Oblate or quasi-spherical prolate: constrain equatorial radius
            if (a_prey[id_ice]) {
                a_new = std::max(a_a[id_ice], a_radius[id_water]*std::cbrt(a_rho_water/rho_rime));
            } else {
                // Water is prey - include gamma inside cube root (Fortran line 1445)
                a_new = std::max(a_a[id_ice], a_radius[id_water]*std::cbrt((a_rho_water/rho_rime)*gamma));
            }
            c_new = V_new / (ParticleReal(four_thirds_pi)*a_new*a_new);
        } else {
            // Prolate or quasi-spherical oblate: constrain polar radius
            if (a_prey[id_ice]) {
                c_new = std::max(a_c[id_ice], a_radius[id_water]*std::cbrt(a_rho_water/rho_rime));
            } else {
                // Water is prey - include gamma inside cube root (Fortran line 1460)
                c_new = std::max(a_c[id_ice], a_radius[id_water]*std::cbrt((a_rho_water/rho_rime)*gamma));
            }
            a_new = std::sqrt(V_new / (ParticleReal(four_thirds_pi)*c_new));
        }
    }

    if ( a_rmndr[a_i] > 0 ) {
        if (a_prey[a_i]) {
            a_mult[a_i] -= gamma * a_mult[a_j];
        } else {
            a_a[a_i] = a_new;
            a_c[a_i] = c_new;
            a_mrime[a_i] = mrime_new;
            a_nmono[a_i] = nmono_new;
            a_Tfz[a_i] = std::max(a_Tfz[a_j], a_Tfz[a_i]);
            // Riming-specific species update: water->0, ice->mi_new
            for (int n = 0; n < a_n_sp; n++) {
                if (n == a_sp_idx_w) { a_sp_m[n][a_i] = ParticleReal(zero); }
                else if (n == a_sp_idx_i) { a_sp_m[n][a_i] = mi_new; }
                else { a_sp_m[n][a_i] += gamma*a_sp_m[n][a_j]; }
            }
            for (int n = 0; n < a_n_ae; n++) { a_ae_m[n][a_i] += gamma * a_ae_m[n][a_j]; }
            AMREX_ALWAYS_ASSERT(mi_new == a_sp_m[a_sp_idx_i][a_i]);
        }
    } else {
        if (a_prey[a_i]) {
            ParticleReal dm = std::floor(a_mult[a_j]/2);
            a_mult[a_i] = dm;
            a_mult[a_j] -= dm;
            a_a[a_i] = a_a[a_j] = a_new;
            a_c[a_i] = a_c[a_j] = c_new;
            a_mrime[a_j] = a_mrime[a_i] = mrime_new;
            a_nmono[a_j] = a_nmono[a_i] = nmono_new;
            a_Tfz[a_i] = a_Tfz[a_j] = std::max(a_Tfz[a_j], a_Tfz[a_i]);
            // Riming-specific species update: water->0, ice->mi_new
            for (int n = 0; n < a_n_sp; n++) {
                if (n == a_sp_idx_w) {
                    a_sp_m[n][a_i] = a_sp_m[n][a_j] = 0.0;
                } else if (n == a_sp_idx_i) {
                    a_sp_m[n][a_i] = a_sp_m[n][a_j] = mi_new;
                } else {
                    a_sp_m[n][a_j] += gamma*a_sp_m[n][a_i];
                    a_sp_m[n][a_i] = a_sp_m[n][a_j];
                }
            }
            for (int n = 0; n < a_n_ae; n++) {
                a_ae_m[n][a_j] += gamma * a_ae_m[n][a_i];
                a_ae_m[n][a_i] = a_ae_m[n][a_j];
            }
            AMREX_ALWAYS_ASSERT(mi_new == a_sp_m[a_sp_idx_i][a_i]);
        }
    }
}

/*! Compute the coalescence of superdroplets in each time step */
void SuperDropletPC::Coalescence( int   a_lev,
                                  Real  a_dt,
                                  const MultiFab& a_pressure,
                                  const MultiFab& a_moist_density,
                                  const MultiFab& a_temperature,
                                  const MultiFab& a_qv,
                                  const Vector<std::unique_ptr<MultiFab>>& a_z_phys_nd )
{
    BL_PROFILE("SuperDropletPC::Coalescence()");
    AMREX_ASSERT( a_lev == m_lev );

    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const auto domain = geom.Domain();

    const auto is_periodic = geom.isPeriodic();
    auto is_periodic_z = is_periodic[2];

    const std::unique_ptr<MultiFab>& z_height = a_z_phys_nd[a_lev];

    const auto rho_water = m_species_mat[m_idx_w]->m_density;
    const auto rho_ice = (m_idx_i >= 0 ? m_species_mat[m_idx_i]->m_density : 0.0);
    const auto sp_idx_w = m_idx_w;
    const auto sp_idx_i = m_idx_i;
    const auto mat_prop(*(m_species_mat[(sp_idx_i>=0?sp_idx_i:sp_idx_w)]));

    const int num_ae = m_num_aerosols;
    const int num_sp  = m_num_species;
    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
    const ParticleReal inv_bin_size
        = ParticleReal(one) / (  static_cast<ParticleReal>(m_coalescence_bin_size[0])
                 * static_cast<ParticleReal>(m_coalescence_bin_size[1])
                 * static_cast<ParticleReal>(m_coalescence_bin_size[2]) );
    const ParticleReal inv_bin_volume = inv_cell_volume*inv_bin_size;

    Real num_collisions = 0;
    const auto& gvec = a_temperature.nGrowVect();

    auto kernel_choice = m_coalescence_kernel;
    auto include_brownian_coalescence = m_include_brownian_coalescence;

    // Build process context
    const auto ctx = buildProcessContext(a_lev);

    // Use serial (non-OMP) iteration because DenseBins is not thread-safe
    forEachParticleTileSerial(a_lev, ctx,
        [&](ParIterType& /*pti*/, int grid, ParticleType* pstruct_ptr,
            const SDProcess::ParticlePointers& ptrs,
            const SDProcess::ProcessContext& /*ctx*/)
    {
        const size_t np = static_cast<size_t>(ptrs.num_particles);
        Box box = a_temperature[grid].box(); box.grow(-gvec);
        int ntiles = numTilesInBox(box, true, m_coalescence_bin_size);
        auto binner = GetParticleBin{plo, dxi, domain, m_coalescence_bin_size, box};
        DenseBins<ParticleType> bins;
        bins.build( np, pstruct_ptr, ntiles, binner);
        AMREX_ALWAYS_ASSERT(np == static_cast<size_t>(bins.numItems()));
        AMREX_ALWAYS_ASSERT(bins.numBins() >= 0);
        auto inds = bins.permutationPtr();
        auto offsets = bins.offsetsPtr();

        BL_PROFILE_VAR("SuperDropletPC::Coalescence::MCShuffle", mcshuffle);

#ifdef AMREX_USE_GPU
        {
            // get the max bin size
            Gpu::Buffer<unsigned int> max_np_bin_d({0});
            auto max_np_bin_ptr = max_np_bin_d.data();
            ParallelFor( bins.numBins(), [=] AMREX_GPU_DEVICE (int i_bin) noexcept
            {
                auto bin_start = offsets[i_bin];
                auto bin_stop = offsets[i_bin+1];
                auto np_bin = bin_stop - bin_start;
                Gpu::Atomic::Max(max_np_bin_ptr, static_cast<unsigned int>(np_bin));
            });
            Gpu::synchronize();
            auto max_np_bin = *(max_np_bin_d.copyToHost());
            // create a stencil vector with integers [0, max_np_bin-1]
            Vector<unsigned int> stencil_vec(max_np_bin);
            for (unsigned int i = 0; i < max_np_bin; i++) { stencil_vec[i] = i; }
            // now shuffle it
            std::shuffle ( stencil_vec.begin(),stencil_vec.end(), m_rndeng );
            // Copy to device
            Gpu::DeviceVector<unsigned int> stencil_vec_d;
            stencil_vec_d.resize(max_np_bin);
            Gpu::copy(  Gpu::hostToDevice,
                        stencil_vec.begin(),
                        stencil_vec.end(),
                        stencil_vec_d.begin() );
            Gpu::synchronize();
            // use this stencil to shuffle each bin
            Gpu::DeviceVector<unsigned int> inds_tmp;
            inds_tmp.resize(np);
            auto stencil_vec_ptr = stencil_vec_d.data();
            auto inds_tmp_ptr = inds_tmp.data();
            ParallelFor( bins.numBins(), [=] AMREX_GPU_DEVICE (int i_bin) noexcept
            {
                auto bin_start = offsets[i_bin];
                auto bin_stop = offsets[i_bin+1];
                auto np_bin = static_cast<unsigned int>(bin_stop - bin_start);
                if (np_bin <= 1) { return; }
                AMREX_ALWAYS_ASSERT(np_bin <= max_np_bin);

                unsigned int j = 0;
                for (unsigned int i = 0; i < max_np_bin; i++) {
                    if (stencil_vec_ptr[i] < np_bin) {
                        inds_tmp_ptr[bin_start+j] = inds[bin_start+stencil_vec_ptr[i]];
                        j++;
                    }
                }
                AMREX_ALWAYS_ASSERT(j == np_bin);
                for (unsigned int i = 0; i < np_bin; i++) {
                    inds[bin_start + i] = inds_tmp_ptr[bin_start + i];
                }
            });
            Gpu::synchronize();
        }
#else
        for (int i_bin = 0; i_bin < bins.numBins(); i_bin++) {
            std::shuffle( inds + offsets[i_bin],
                          inds + offsets[i_bin+1],
                          m_rndeng );
        }
#endif

        BL_PROFILE_VAR_STOP(mcshuffle);

        const auto& pressure_arr = a_pressure[grid].const_array();
        const auto& temperature_arr = a_temperature[grid].const_array();
        const auto& moist_density_arr = a_moist_density[grid].const_array();
        const auto& qv_arr = a_qv[grid].const_array();
        auto zheight = (*z_height)[grid].array();

        CollisionKernel<ParticleReal,AMREX_SPACEDIM> ckernel{};

        Gpu::Buffer<Real> particle_collisions({0});
        auto particle_collisions_ptr = particle_collisions.data();

        // initialize coalescence quantities
        Gpu::DeviceVector<int> coll_partner_idx, flag_prey, num_particles_bin;
        Gpu::DeviceVector<Real> coll_rate, coll_rmndr;
        num_particles_bin.resize(np);
        coll_partner_idx.resize(np);
        flag_prey.resize(np);
        coll_rate.resize(np);
        coll_rmndr.resize(np);
        auto np_bin_ptr = num_particles_bin.data();
        auto partner_idx_ptr = coll_partner_idx.data();
        auto flag_prey_ptr = flag_prey.data();
        auto coll_rate_ptr = coll_rate.data();
        auto coll_rmndr_ptr = coll_rmndr.data();
        ParallelFor( np, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            np_bin_ptr[i] = 0;
            partner_idx_ptr[i] = -1;
            flag_prey_ptr[i] = -1;
            coll_rate_ptr[i] = ParticleReal(zero);
            coll_rmndr_ptr[i] = ParticleReal(zero);
        });
        Gpu::synchronize();

        BL_PROFILE_VAR("SuperDropletPC::Coalescence::MCPairing", mcpairing);

        // create pairs: note that the SD with larger multiplicity is the "prey"
        ParallelFor( bins.numBins(), [=] AMREX_GPU_DEVICE (int i_bin) noexcept
        {
            auto bin_start = offsets[i_bin];
            auto bin_stop = offsets[i_bin+1];
            auto np_bin = bin_stop - bin_start;
            if (np_bin <= 1) { return; }

            for (int p = 0; p < np_bin/2; p++) {
                auto pi = inds[bin_start+p];
                auto pj = inds[bin_stop-1-p];

                if (pi == pj) { continue; }
                if (ptrs.active_ptr[pi] == 0) { continue; }
                if (ptrs.active_ptr[pj] == 0) { continue; }
                if (ptrs.mult_ptr[pi] == 0) { continue; }
                if (ptrs.mult_ptr[pj] == 0) { continue; }

                np_bin_ptr[pi] = np_bin_ptr[pj] = np_bin;
                partner_idx_ptr[pi] = pj;
                partner_idx_ptr[pj] = pi;

                int i = -1, j = -1;
                if (ptrs.mult_ptr[pi] >= ptrs.mult_ptr[pj]) { i = pi; j = pj; }
                else                              { i = pj; j = pi; }
                flag_prey_ptr[i] = 1;
                flag_prey_ptr[j] = 0;
            }
        } );
        Gpu::synchronize();

        BL_PROFILE_VAR_STOP(mcpairing);

        BL_PROFILE_VAR("SuperDropletPC::Coalescence::CollisionCalc", coalescence_calc);

        // calculate collision efficiencies for each pair
        ParallelForRNG( np, [=] AMREX_GPU_DEVICE (int i, RandomEngine const& rnd_eng) noexcept
        {
            if (partner_idx_ptr[i] < 0) { return; }
            if (!flag_prey_ptr[i]) { return; }

            int pi = i; // prey - higher multiplicity
            int pj = partner_idx_ptr[i]; // predator - lower multiplicity
            AMREX_ALWAYS_ASSERT(ptrs.mult_ptr[pi] >= ptrs.mult_ptr[pj]);

            // get phases for the two particles
            auto phase_i = SD_phase(pi, sp_idx_w, sp_idx_i, ptrs.sp_mass_ptrs);
            auto phase_j = SD_phase(pj, sp_idx_w, sp_idx_i, ptrs.sp_mass_ptrs);

            ParticleReal k_val = ParticleReal(zero);
            if ((phase_i == SDPhase::water) && (phase_j == SDPhase::water)) {

                // coalescence between two water droplets

                if (kernel_choice == SDCoalescenceKernelType::golovin) {

                    k_val = ckernel.golovin(ptrs.radius_ptr[pi],ptrs.radius_ptr[pj]);

                } else {

                    ParticleReal v_i[AMREX_SPACEDIM], v_j[AMREX_SPACEDIM];
                    for (int d = 0; d < AMREX_SPACEDIM; d++) {
                        v_i[d] = ptrs.v_ptr[d][pi];
                        v_j[d] = ptrs.v_ptr[d][pj];
                    }
                    v_i[AMREX_SPACEDIM-1] -= ptrs.vterm_ptr[pi];
                    v_j[AMREX_SPACEDIM-1] -= ptrs.vterm_ptr[pj];

                    if (kernel_choice == SDCoalescenceKernelType::sedimentation) {
                        k_val = ckernel.sedimentation(ptrs.radius_ptr[pi],ptrs.radius_ptr[pj],v_i,v_j);
                    } else if (kernel_choice == SDCoalescenceKernelType::Longs) {
                        k_val = ckernel.Longs(ptrs.radius_ptr[pi],ptrs.radius_ptr[pj],v_i,v_j);
                    } else if (kernel_choice == SDCoalescenceKernelType::Halls) {
                        k_val = ckernel.Halls(ptrs.radius_ptr[pi],ptrs.radius_ptr[pj],v_i,v_j);
                    }

                    if (k_val < 0.0) {
                        amrex::Abort("Invalid value for k_val");
                    }
                }

            } else if ((phase_i == SDPhase::ice) && (phase_j == SDPhase::ice) && (sp_idx_i >= 0)) {

                // aggregation between two ice particles

                // ice particle 1
                auto vz_i = ptrs.v_ptr[AMREX_SPACEDIM-1][pi] - ptrs.vterm_ptr[pi];
                auto a_i = ptrs.a_ptr[pi];
                auto c_i = ptrs.c_ptr[pi];
                auto rhoi_i = ice_rho(a_i, c_i, ptrs.sp_mass_ptrs[sp_idx_i][pi]);
                auto maxR_i = std::max(a_i, c_i);
                auto k_i = std::exp(-ckernel.k_coeff * c_i/a_i);
                auto area_i = PI * a_i * maxR_i * std::pow(rhoi_i/rho_ice, k_i);

                // ice particle 2
                auto vz_j = ptrs.v_ptr[AMREX_SPACEDIM-1][pj] - ptrs.vterm_ptr[pj];
                auto a_j = ptrs.a_ptr[pj];
                auto c_j = ptrs.c_ptr[pj];
                auto rhoi_j = ice_rho(a_j, c_j, ptrs.sp_mass_ptrs[sp_idx_i][pj]);
                auto maxR_j = std::max(a_j, c_j);
                auto k_j = std::exp(-ckernel.k_coeff * c_j/a_j);
                auto area_j = PI * a_j * maxR_j * std::pow(rhoi_j/rho_ice, k_j);

                // velocity difference
                auto dvz = std::sqrt((vz_i-vz_j)*(vz_i-vz_j));

                // collision efficiency
                k_val = 0.1 * (std::sqrt(area_i)+std::sqrt(area_j))*(std::sqrt(area_i)+std::sqrt(area_j)) * dvz;

            } else if (sp_idx_i >= 0) {

                // riming between a water droplet and an ice particle
                AMREX_ALWAYS_ASSERT(    ((phase_i == SDPhase::ice) && (phase_j == SDPhase::water))
                                     || ((phase_i == SDPhase::water) && (phase_j == SDPhase::ice)) );

                // figure out which one is water, and which is ice
                int id_w = -1, id_i = -1;
                if (phase_i == SDPhase::water) { id_w = pi; id_i = pj; }
                else { id_w = pj; id_i = pi; }

                // water droplet
                auto vz_w = ptrs.v_ptr[AMREX_SPACEDIM-1][id_w] - ptrs.vterm_ptr[id_w];
                auto r_w = ptrs.radius_ptr[id_w];

                // ice particle
                auto vz_i = ptrs.v_ptr[AMREX_SPACEDIM-1][id_i] - ptrs.vterm_ptr[id_i];
                auto a_i = ptrs.a_ptr[id_i];
                auto c_i = ptrs.c_ptr[id_i];
                auto rhoi_i = ice_rho(a_i, c_i, ptrs.sp_mass_ptrs[sp_idx_i][id_i]);
                auto eqr_i = std::cbrt(a_i*a_i*c_i);
                auto maxR_i = std::max(a_i, c_i);

                // interpolate flow quantities at ice particle location
                ParticleType& p_ice = pstruct_ptr[id_i];
                constexpr int nf = static_cast<int>(InterpFieldsTR::NUM_FIELDS);
                ParticleReal fv[nf];
                const Array4<const Real> fa[nf] = {
                    temperature_arr, moist_density_arr
                };
                ERF::Interpolation::interpolateFields(
                    p_ice, plo, dxi, fa, fv, nf,
                    is_periodic_z ? 1 : 0, is_periodic_z ? nullptr : &zheight
                );
                const auto temperature   = fv[static_cast<int>(InterpFieldsTR::temperature)];
                const auto moist_density = fv[static_cast<int>(InterpFieldsTR::moist_density)];

                // dynamic viscosity
                auto mu = viscCoeff(temperature);

                // velocity difference
                auto dvz = std::sqrt((vz_i-vz_w)*(vz_i-vz_w));

                if (std::abs(vz_w) > std::abs(vz_i)) {

                    // collector is water droplet

                    auto size_ratio = eqr_i / r_w;
                    auto Re = moist_density * (2.0*r_w) * std::abs(vz_w) / mu;
                    auto Kn = ckernel.lambda_mfp / (2.0*eqr_i);
                    // Cunningham slip correction factor
                    auto C_sc = 1.0 + Kn * (ckernel.para_Asc + ckernel.para_Bsc * std::exp(-ckernel.para_Csc/Kn));
                    // Stokes number
                    auto St = size_ratio*size_ratio * rhoi_i * Re * C_sc / (9.0*moist_density);

                    k_val = ckernel.BeardGrover1974(size_ratio, Re, St);

                    // collision efficiency
                    k_val *= (PI * (r_w+a_i) * (r_w+maxR_i) * dvz);

                } else {

                    // collector is ice particle

                    auto size_ratio = r_w / eqr_i;
                    auto Re = moist_density * 2*maxR_i * std::abs(vz_i) / mu;
                    auto mixFr = std::abs((vz_i - vz_w) * vz_w) / (maxR_i * CONST_GRAV);

                    k_val = ckernel.BeardGrover1974(size_ratio, Re, mixFr);

                    auto phi_i = c_i / a_i;
                    if (phi_i < 1.0) {
                        // oblate
                        auto tmp = ckernel.ErfaniMitchell2017_Plate(Re, mixFr);
                        k_val = phi_i*k_val + (1.0-phi_i)*tmp;
                    } else {
                        // prolate
                        auto Re_col = moist_density * (2*a_i) * std::abs(vz_i) / mu;
                        auto tmp = ckernel.ErfaniMitchell2017_Column(Re_col, mixFr);
                        k_val = (1.0/phi_i)*k_val + (1.0-(1.0/phi_i))*tmp;
                    }

                    // circumscribed ellipse area of the particle projected to the flow direction
                    auto area_i_ce = PI * a_i * maxR_i;
                    // area of the particle projected to the flow direction
                    auto k_i = std::exp(-ckernel.k_coeff * phi_i);
                    auto area_i = area_i_ce * std::pow(rhoi_i/rho_ice, k_i);

                    // collision-riming kernel: E_rime*A_g*|vj-vk|
                    k_val *= ( (PI*(r_w+a_i)*(r_w+maxR_i) - (area_i_ce-area_i)) * dvz);
                }

            }

            if (include_brownian_coalescence) {

                ParticleReal pressure = ParticleReal(zero), temperature = ParticleReal(zero);
                {
                    ParticleType& par_1 = pstruct_ptr[pi];
                    auto iv = getParticleCell(par_1, plo, dxi, domain);
                    pressure = pressure_arr(iv[0],iv[1],iv[2],0);
                    temperature = temperature_arr(iv[0],iv[1],iv[2],0);
                }

                ParticleReal sd_mass_1 = ParticleReal(zero),
                             sd_mass_2 = ParticleReal(zero);
                for (int ia = 0; ia < num_ae; ia++) {
                    sd_mass_1 += ptrs.ae_mass_ptrs[ia][pi];
                    sd_mass_2 += ptrs.ae_mass_ptrs[ia][pj];
                }
                for (int ia = 0; ia < num_sp; ia++) {
                    sd_mass_1 += ptrs.sp_mass_ptrs[ia][pi];
                    sd_mass_2 += ptrs.sp_mass_ptrs[ia][pj];
                }

                auto r_eff_1 = SD_effective_radius( pi, sp_idx_w,
                                                    rho_water,
                                                    num_sp, num_ae,
                                                    ptrs.sp_sol_arr, ptrs.ae_sol_arr,
                                                    ptrs.sp_mass_ptrs, ptrs.ae_mass_ptrs,
                                                    ptrs.sp_rho_arr, ptrs.ae_rho_arr );
                auto r_eff_2 = SD_effective_radius( pj, sp_idx_w,
                                                    rho_water,
                                                    num_sp, num_ae,
                                                    ptrs.sp_sol_arr, ptrs.ae_sol_arr,
                                                    ptrs.sp_mass_ptrs, ptrs.ae_mass_ptrs,
                                                    ptrs.sp_rho_arr, ptrs.ae_rho_arr );

                auto k_brown = ckernel.Brownian_SeinfeldPandis( r_eff_1,
                                                                r_eff_2,
                                                                sd_mass_1,
                                                                sd_mass_2,
                                                                pressure,
                                                                temperature );
                if (k_brown < ParticleReal(zero)) {
                    amrex::Abort("Invalid value for k_brown");
                }

                k_val += k_brown;
            }


            auto prob_ij = k_val*inv_bin_volume;
            auto prob_sd_ij = std::max(ptrs.mult_ptr[pi],ptrs.mult_ptr[pj])*prob_ij;

            auto ns = static_cast<ParticleReal>(np_bin_ptr[i]);
            auto scaling_factor = myhalf*ns*(ns-1)/std::floor(myhalf*ns);
            auto scaled_prob = prob_sd_ij * scaling_factor;

            auto gamma = coalescence_rate ( rnd_eng, (scaled_prob*a_dt) );
            if (gamma > 0) {
                amrex::Gpu::Atomic::Add(particle_collisions_ptr, gamma);
                coll_rate_ptr[pi] = std::min(gamma,std::floor(ptrs.mult_ptr[pi]/ptrs.mult_ptr[pj]));
                coll_rate_ptr[pj] = coll_rate_ptr[pi];
                coll_rmndr_ptr[pi] = ptrs.mult_ptr[pi] - coll_rate_ptr[pi]*ptrs.mult_ptr[pj];
                coll_rmndr_ptr[pj] = coll_rmndr_ptr[pi];
            } else {
                partner_idx_ptr[pi] = -1;
                partner_idx_ptr[pj] = -1;
            }

        } );
        Gpu::synchronize();
        num_collisions = *(particle_collisions.copyToHost());

        dMdt<ParticleReal> dmdt{ mat_prop.m_lat_vap,
                                 therco, /* ERF_Constants.H */
                                 mat_prop.m_Rv,
                                 mat_prop.m_density };

        ParallelFor( np, [=] AMREX_GPU_DEVICE (int i)
        {
            const auto j = partner_idx_ptr[i];
            if (j < 0) { return; }

            // get phases for the two particles
            auto phase_i = SD_phase(i, sp_idx_w, sp_idx_i, ptrs.sp_mass_ptrs);
            auto phase_j = SD_phase(j, sp_idx_w, sp_idx_i, ptrs.sp_mass_ptrs);

            if ((phase_i == SDPhase::water) && (phase_j == SDPhase::water)) {

                // coalescence between two water droplets
                coal_update_attribs( i, j,
                                     flag_prey_ptr,
                                     coll_rate_ptr,
                                     coll_rmndr_ptr,
                                     ptrs.mass_ptr,
                                     ptrs.radius_ptr,
                                     ptrs.Tfz_ptr,
                                     ptrs.mult_ptr,
                                     num_sp,
                                     ptrs.sp_mass_ptrs,
                                     num_ae,
                                     ptrs.ae_mass_ptrs );

            } else if ((phase_i == SDPhase::ice) && (phase_j == SDPhase::ice) && (sp_idx_i >= 0)) {

                // aggregation between two ice particles
                aggr_update_attribs( i, j,
                                     sp_idx_i,
                                     rho_ice,
                                     10.0,
                                     flag_prey_ptr,
                                     coll_rate_ptr,
                                     coll_rmndr_ptr,
                                     ptrs.Tfz_ptr,
                                     ptrs.a_ptr,
                                     ptrs.c_ptr,
                                     ptrs.mrime_ptr,
                                     ptrs.nmono_ptr,
                                     ptrs.mult_ptr,
                                     num_sp,
                                     ptrs.sp_mass_ptrs,
                                     num_ae,
                                     ptrs.ae_mass_ptrs );

            } else if (sp_idx_i >= 0) {

                // riming between a water droplet and an ice particle
                AMREX_ALWAYS_ASSERT(    ((phase_i == SDPhase::ice) && (phase_j == SDPhase::water))
                                     || ((phase_i == SDPhase::water) && (phase_j == SDPhase::ice)) );

                ParticleType& p_ice = pstruct_ptr[i];
                constexpr int nf = static_cast<int>(InterpFieldsFull::NUM_FIELDS);
                ParticleReal fv[nf];
                const Array4<const Real> fa[nf] = {
                    temperature_arr, pressure_arr, moist_density_arr, qv_arr
                };
                ERF::Interpolation::interpolateFields(
                    p_ice, plo, dxi, fa, fv, nf,
                    is_periodic_z ? 1 : 0, is_periodic_z ? nullptr : &zheight
                );
                const auto temperature   = fv[static_cast<int>(InterpFieldsFull::temperature)];
                const auto pressure      = fv[static_cast<int>(InterpFieldsFull::pressure)];
                const auto moist_density = fv[static_cast<int>(InterpFieldsFull::moist_density)];
                const auto qv            = fv[static_cast<int>(InterpFieldsFull::qv)];
                auto coeff_moldiff = mat_prop.coeffMolecularDiffusion(temperature, pressure);

                rime_update_attribs( i, j,
                                     sp_idx_w, sp_idx_i,
                                     rho_water, rho_ice,
                                     flag_prey_ptr,
                                     coll_rate_ptr,
                                     coll_rmndr_ptr,
                                     phase_i, phase_j,
                                     ptrs.v_ptr,
                                     ptrs.vterm_ptr,
                                     ptrs.radius_ptr,
                                     ptrs.Tfz_ptr,
                                     ptrs.a_ptr,
                                     ptrs.c_ptr,
                                     ptrs.mrime_ptr,
                                     ptrs.nmono_ptr,
                                     ptrs.mult_ptr,
                                     temperature,
                                     moist_density,
                                     pressure,
                                     qv,
                                     coeff_moldiff,
                                     dmdt,
                                     num_sp,
                                     ptrs.sp_mass_ptrs,
                                     num_ae,
                                     ptrs.ae_mass_ptrs );
            }

        } );
        Gpu::synchronize();

        // update particle radius and total mass
        ParallelFor( np, [=] AMREX_GPU_DEVICE (int i)
        {
            SuperDropletPC::updateParticleAttributes(
                i, ptrs.radius_ptr, ptrs.mass_ptr, sp_idx_w, rho_water,
                num_sp, num_ae, ptrs.sp_sol_arr, ptrs.ae_sol_arr,
                ptrs.sp_mass_ptrs, ptrs.ae_mass_ptrs, ptrs.sp_rho_arr, ptrs.ae_rho_arr);
        } );
        Gpu::synchronize();

        BL_PROFILE_VAR_STOP(coalescence_calc);
    }); // end forEachParticleTileSerial

    ParallelDescriptor::ReduceRealSum(  &num_collisions,
                                        1,
                                        ParallelDescriptor::IOProcessorNumber() );

    Print() << "SuperDropletPC(" << m_name << "): "
            << "number of collisions = " << num_collisions << "\n";
}

#endif
