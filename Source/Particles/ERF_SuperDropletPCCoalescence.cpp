#ifndef _WIN32
#include <sys/time.h>
#endif
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

/*! \brief Compute dynamic viscosity */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static auto viscCoeff ( const ParticleReal a_T /*!< temperature */ )
{
    auto T_degC = a_T - 273.15; // [K] => [degC]
    ParticleReal visc_coeff = 0.0;
    if( T_degC >= 0.0 ) {
        visc_coeff = ( 1.7180 + 4.9E-3*T_degC ) * 1.E-5;
    } else {
        visc_coeff = ( 1.7180 + 4.9E-3*T_degC -1.2E-5*T_degC*T_degC ) * 1.E-5;
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
    AMREX_ALWAYS_ASSERT(a_rmndr[a_i] >= 0.0);

    if ( a_rmndr[a_i] > 0 ) {

        if (a_prey[i]) {
            a_mult[i] -= gamma*a_mult[j];
        } else {
            auto r3 = gamma*a_radius[j]*a_radius[j]*a_radius[j]
                            + a_radius[i]*a_radius[i]*a_radius[i];
            a_radius[i] = std::cbrt(r3);
            a_mass[i] += gamma*a_mass[j];
            // Update Tfz: combined droplet inherits the warmest (least negative) Tfz
            // from either parent since both INP surfaces are now present
            a_Tfz[i] = std::max(a_Tfz[j], a_Tfz[i]);
            for (int n = 0; n < a_n_sp; n++) {
                a_sp_m[n][i] += gamma*a_sp_m[n][j];
            }
            for (int n = 0; n < a_n_ae; n++) {
                a_ae_m[n][i] += gamma*a_ae_m[n][j];
            }
        }

    } else if ( a_rmndr[a_i] == 0 ) {

        if (a_prey[i]) {
            ParticleReal dm = std::floor(a_mult[j]/2);
            a_mult[i] = dm;
            a_mult[j] -= dm;
            ParticleReal r3 = gamma*a_radius[i]*a_radius[i]*a_radius[i]
                                    + a_radius[j]*a_radius[j]*a_radius[j];
            a_radius[i] = a_radius[j] = std::cbrt(r3);
            a_mass[j] += gamma*a_mass[i];
            a_mass[i] = a_mass[j];
            // Update Tfz: both particles get the warmest (least negative) Tfz
            a_Tfz[i] = a_Tfz[j] = std::max(a_Tfz[j], a_Tfz[i]);
            for (int n = 0; n < a_n_sp; n++) {
                a_sp_m[n][j] += gamma*a_sp_m[n][i];
                a_sp_m[n][i] = a_sp_m[n][j];
            }
            for (int n = 0; n < a_n_ae; n++) {
                a_ae_m[n][j] += gamma*a_ae_m[n][i];
                a_ae_m[n][i] = a_ae_m[n][j];
            }
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
    AMREX_ALWAYS_ASSERT(a_rmndr[a_i] >= 0.0);

    auto gamma = a_gamma[a_i];
    ParticleReal a_new = 0.0;
    ParticleReal c_new = 0.0;
    ParticleReal m_new = 0.0;

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
        auto V_new_min = (4.0*PI/3.0) * (gamma*a_a[i]*a_a[i]*a_c[i] + a_a[j]*a_a[j]*a_c[j]);
        auto rhoi_bar = m_new / V_new_min;

        auto maxR_i = std::max(a_a[i], a_c[i]);
        auto maxR_j = std::max(a_a[j], a_c[j]);

        if (maxR_i > maxR_j) {
            if (a_a[i] > a_c[i]) {
                a_new = a_a[i];
                auto c_new_min = V_new_min / ((4.0*PI/3.0)*a_new*a_new);
                auto c_new_max = a_c[i]*gamma + std::min(a_a[j],a_c[j]);
                c_new = (a_rho_ice-a_rho_min) / ((a_rho_ice-rhoi_bar)/c_new_min + (rhoi_bar-a_rho_min)/c_new_max);
            } else {
                c_new = a_c[i];
                auto a_new_min = std::sqrt(V_new_min / ((4.0*PI/3.0)*c_new));
                auto a_new_max = std::sqrt( std::max(std::max(a_a[i],a_a[j]),a_c[j])
                                            * (a_a[i]*gamma + std::min(a_a[j],a_c[j])) );
                a_new = std::sqrt( (a_rho_ice - a_rho_min)
                                    / (   (a_rho_ice-rhoi_bar)/(a_new_min*a_new_min)
                                        + (rhoi_bar-a_rho_min)/(a_new_max*a_new_max)) );
            }
        } else {
            if (a_a[j] > a_c[j]) {
                a_new = a_a[j];
                auto c_new_min = V_new_min / ((4.0*PI/3.0)*a_new*a_new);
                auto c_new_max = a_c[j] + std::min(a_a[i],a_c[i])*gamma;
                c_new = (a_rho_ice-a_rho_min) / ((a_rho_ice-rhoi_bar)/c_new_min + (rhoi_bar-a_rho_min)/c_new_max);
            } else {
                c_new = a_c[j];
                auto a_new_min = std::sqrt(V_new_min / ((4.0*PI/3.0)*c_new));
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
            a_mult[a_i] -= gamma*a_mult[a_j];
        } else {
            a_a[a_i] = a_new;
            a_c[a_i] = c_new;
            a_mrime[a_i] += gamma * a_mrime[a_j];
            a_nmono[a_i] += gamma * a_nmono[a_j];
            a_Tfz[a_i] = std::max(a_Tfz[a_j], a_Tfz[a_i]);
            for (int n = 0; n < a_n_sp; n++) { a_sp_m[n][a_i] += gamma*a_sp_m[n][a_j]; }
            for (int n = 0; n < a_n_ae; n++) { a_ae_m[n][a_i] += gamma*a_ae_m[n][a_j]; }
            AMREX_ALWAYS_ASSERT(m_new == a_sp_m[a_sp_idx_i][a_i]);
        }
    } else {
        if (a_prey[a_i]) {
            auto dm = std::floor(a_mult[a_j]/2);
            a_mult[a_i] = dm;
            a_mult[a_j] -= dm;
            a_a[a_i] = a_a[a_j] = a_new;
            a_c[a_i] = a_c[a_j] = c_new;
            a_mrime[a_j] += gamma * a_mrime[a_i]; a_mrime[a_i] = a_mrime[a_j];
            a_nmono[a_j] += gamma * a_nmono[a_i]; a_nmono[a_i] = a_nmono[a_j];
            a_Tfz[a_i] = a_Tfz[a_j] = std::max(a_Tfz[a_j], a_Tfz[a_i]);
            for (int n = 0; n < a_n_sp; n++) {
                a_sp_m[n][a_j] += gamma*a_sp_m[n][a_i];
                a_sp_m[n][a_i] = a_sp_m[n][a_j];
            }
            for (int n = 0; n < a_n_ae; n++) {
                a_ae_m[n][a_j] += gamma*a_ae_m[n][a_i];
                a_ae_m[n][a_i] = a_ae_m[n][a_j];
            }
            AMREX_ALWAYS_ASSERT(m_new == a_sp_m[a_sp_idx_i][a_i]);
        }
    }
}

/*! \brief Rasmussen and Heymsfield (1985) for impact velocity */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static auto impactVelocity_RasmussenHeymsfield1985( const ParticleReal a_Re, /*!< Reynolds number */
                                                    const ParticleReal a_St  /*!< Stokes number */ )
{
    auto w = std::log10(a_St);
    auto w2 = w*w;
    auto w3 = w2*w;
    auto w4 = w3*w;
    ParticleReal retval = 0.0;;
    if(a_Re < (10.0+30.0)*0.5) {
        if (a_St < 0.4) {
            retval = 0.0;
        } else if (a_St<10.0) {
            retval = 0.1701 + 0.7246*w + 0.2257*w2 - 1.13*w3 + 0.5756*w4;
        } else {
            retval = 0.57;
        }
    } else if (a_Re < (30.0+100.0)*0.5) {
        if (a_St < 0.1) {
            retval = 0.0;
        } else if (a_St < 10.0) {
            retval = 0.2927 + 0.5085*w - 0.03453*w2 - 0.2184*w3 + 0.03595*w4;
        } else {
            retval = 0.59;
        }
    } else if (a_Re < (100.0+300.0)*0.5) {
        if (a_St < 0.1) {
          retval = 0.0;
        } else if (a_St < 10.0) {
            retval = 0.3272 + 0.4907*w - 0.09452*w2 - 0.1906*w3 + 0.07105*w4;
        } else {
            retval = 0.61;
        }
    } else {
        if (a_St < 0.1) {
            retval = 0.0;
        } else if (a_St < 10.0) {
            retval = 0.356 + 0.4738*w - 0.1233*w2 - 0.1618*w3 + 0.08087*w4;
        } else {
            retval = 0.63;
        }
    }
    retval = std::max(retval,0.0);
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
    Real qsat = 0.0; erf_qsati(a_T, a_P, qsat);
    auto sup_sat  = a_qv/qsat - 1.0;
    auto drho  = sup_sat / (a_D * a_dmdt.Fk_plus_Fd(a_T, erf_esati(a_T), a_D));
    auto dT = (a_dmdt.L * a_D / a_dmdt.K) * drho;
    return a_T + dT;
}

/*! \brief Heymsfield and Pflaum (1985) formula for rime density */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static auto rimeDensity_HeymsfieldPflaum1985( const ParticleReal a_radius,      /*!< droplet radius */
                                              const ParticleReal a_a,           /*!< ice particle equatorial radius */
                                              const ParticleReal a_c,           /*!< ice particle polar radius */
                                              const ParticleReal a_vz_w,        /*!< vertical velocity of water droplet */
                                              const ParticleReal a_vz_i,        /*!< vertical velocity of ice particle */
                                              const ParticleReal a_rho_w,       /*!< water density */
                                              const ParticleReal a_D, /*!< diffusivity coeff */
                                              const dMdt<ParticleReal>& a_dmdt, /*!< mass change utilities */
                                              const ParticleReal a_T,           /*!< temperature */
                                              const ParticleReal a_rhom,        /*!< moist density */
                                              const ParticleReal a_P,           /*!< pressure */
                                              const ParticleReal a_qv           /*!< vapour fraction */)
{
    auto r_um = a_radius * 1.0e6;
    auto mu = viscCoeff(a_T);

    auto maxD = 2.0 * std::max(a_a,a_c);
    auto Re = a_rhom * maxD * std::abs(a_vz_i) / mu;
    auto eqr_i = std::cbrt(a_a*a_a*a_c);
    auto St = 2.0*std::abs(a_vz_i)*a_radius*a_radius*a_rho_w/(9.0*mu*eqr_i);

    auto v_impact_ratio = impactVelocity_RasmussenHeymsfield1985(Re,St);
    auto v_impact = std::abs(a_vz_i - a_vz_w) * v_impact_ratio;

    auto T_surf = std::min(-0.01, iceSurfaceTemperature(a_T,a_P,a_qv,a_D,a_dmdt) - 273.15);
    auto var_Y = -r_um * v_impact/T_surf;

    ParticleReal rho_rime = 0.0;
    if ((T_surf <= -5.0) || (var_Y < 1.6)) {
        rho_rime = std::exp(std::log(0.30*var_Y)*0.44);
    } else {
        var_Y = std::min(var_Y,3.5);
        rho_rime = std::exp(-0.03115 - 1.7030*var_Y + 0.9116*var_Y*var_Y - 0.1224*var_Y*var_Y*var_Y);
    }
    rho_rime = std::min(0.91,std::max(rho_rime,0.1)) * 1000.0;
    return rho_rime;
}

/*! \brief Binary aggregation between two superdroplets */
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
                                const Array<ParticleReal*,AMREX_SPACEDIM>& a_vel, /*!< velocity */
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
    AMREX_ALWAYS_ASSERT(a_rmndr[a_i] >= 0.0);
    auto gamma = a_gamma[a_i];

    ParticleReal mi_new = 0.0;
    ParticleReal a_new = 0.0;
    ParticleReal c_new = 0.0;
    ParticleReal mrime_new = 0.0;
    ParticleReal nmono_new = 0.0;

    if (a_phase_i == SDPhase::ice) {
        AMREX_ALWAYS_ASSERT(a_phase_j == SDPhase::water);
        mi_new = gamma*a_sp_m[a_sp_idx_i][a_i] + a_sp_m[a_sp_idx_w][a_j];
        if (a_radius[a_j] > std::max(a_a[a_i], a_c[a_i])) {
            a_new = c_new = std::cbrt(mi_new / ((4.0*PI/3.0)*a_rho_ice));
        } else {
            auto vz_w = a_vel[AMREX_SPACEDIM-1][a_j] - a_vterm[a_j];
            auto vz_i = a_vel[AMREX_SPACEDIM-1][a_i] - a_vterm[a_i];
            auto rho_rime = rimeDensity_HeymsfieldPflaum1985( a_radius[a_j],
                                                              a_a[a_i], a_c[a_i],
                                                              vz_w, vz_i,
                                                              a_rho_water,
                                                              a_D, a_dmdt,
                                                              a_T, a_rhom, a_P, a_qv );
            auto V_new = (4.0*PI/3.0) * (   a_a[a_i]*a_a[a_i]*a_c[a_i] * gamma
                                          + a_radius[a_j]*a_radius[a_j]*a_radius[a_j]*(a_rho_water/rho_rime) );
            auto phi = a_c[a_i] / a_a[a_i];
            if ((phi < 0.8) || ((phi < 1.25) && (phi >= 1.0))) {
                a_new = std::max(a_a[a_i], a_radius[a_j]*std::cbrt(a_rho_water/rho_rime));
                c_new = V_new / ((4.0*PI/3.0)*a_new*a_new);
            } else {
                c_new = std::max(a_c[a_i], a_radius[a_j]*std::cbrt(a_rho_water/rho_rime));
                a_new = std::sqrt(V_new / ((4.0*PI/3.0)*c_new));
            }
        }
        mrime_new = gamma*a_mrime[a_i] + a_sp_m[a_sp_idx_w][a_j];
        nmono_new = gamma*a_nmono[a_i];
    } else {
        AMREX_ALWAYS_ASSERT(a_phase_j == SDPhase::ice);
        mi_new = a_sp_m[a_sp_idx_i][a_i] + gamma*a_sp_m[a_sp_idx_w][a_j];
        if (a_radius[a_i] > std::max(a_a[a_j], a_c[a_j])) {
            a_new = c_new = std::cbrt(mi_new / ((4.0*PI/3.0)*a_rho_ice));
        } else {
            auto vz_w = a_vel[AMREX_SPACEDIM-1][a_i] - a_vterm[a_i];
            auto vz_i = a_vel[AMREX_SPACEDIM-1][a_j] - a_vterm[a_j];
            auto rho_rime = rimeDensity_HeymsfieldPflaum1985( a_radius[a_i],
                                                              a_a[a_j], a_c[a_j],
                                                              vz_w, vz_i,
                                                              a_rho_water,
                                                              a_D, a_dmdt,
                                                              a_T, a_rhom, a_P, a_qv );
            auto V_new = (4.0*PI/3.0) * (   a_a[a_j]*a_a[a_j]*a_c[a_j]
                                          + a_radius[a_i]*a_radius[a_i]*a_radius[a_i]*(a_rho_water/rho_rime)*gamma );
            auto phi = a_c[a_j] / a_a[a_j];
            if ((phi < 0.8) || ((phi < 1.25) && (phi >= 1.0))) {
                a_new = std::max(a_a[a_j], a_radius[a_i]*std::cbrt(a_rho_water/rho_rime));
                c_new = V_new / ((4.0*PI/3.0)*a_new*a_new);
            } else {
                c_new = std::max(a_c[a_j], a_radius[a_i]*std::cbrt(a_rho_water/rho_rime));
                a_new = std::sqrt(V_new / ((4.0*PI/3.0)*c_new));
            }
        }
        mrime_new = a_mrime[a_j] + gamma*a_sp_m[a_sp_idx_w][a_i];
        nmono_new = a_nmono[a_j];
    }

    if ( a_rmndr[a_i] > 0 ) {
        if (a_prey[a_i]) {
            a_mult[a_i] -= gamma*a_mult[a_j];
        } else {
            a_a[a_i] = a_new;
            a_c[a_i] = c_new;
            a_mrime[a_i] = mrime_new;
            a_nmono[a_i] = nmono_new;
            a_Tfz[a_i] = std::max(a_Tfz[a_j], a_Tfz[a_i]);
            for (int n = 0; n < a_n_sp; n++) {
                if (n == a_sp_idx_w) { a_sp_m[n][a_i] = 0.0; }
                else if (n == a_sp_idx_i) { a_sp_m[n][a_i] = mi_new; }
                else { a_sp_m[n][a_i] += gamma*a_sp_m[n][a_j]; }
            }
            for (int n = 0; n < a_n_ae; n++) { a_ae_m[n][a_i] += gamma*a_ae_m[n][a_j]; }
            AMREX_ALWAYS_ASSERT(mi_new == a_sp_m[a_sp_idx_i][a_i]);
        }
    } else {
        if (a_prey[a_i]) {
            auto dm = std::floor(a_mult[a_j]/2);
            a_mult[a_i] = dm;
            a_mult[a_j] -= dm;
            a_a[a_i] = a_a[a_j] = a_new;
            a_c[a_i] = a_c[a_j] = c_new;
            a_mrime[a_j] = a_mrime[a_i] = mrime_new;
            a_nmono[a_j] = a_nmono[a_i] = nmono_new;
            a_Tfz[a_i] = a_Tfz[a_j] = std::max(a_Tfz[a_j], a_Tfz[a_i]);
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
                a_ae_m[n][a_j] += gamma*a_ae_m[n][a_i];
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
#ifndef _WIN32
    struct timeval total_start, total_end;
    gettimeofday(&total_start, NULL);
#endif

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
    const auto mat_prop(*(m_species_mat[sp_idx_i]));

    const int num_ae = m_num_aerosols;
    const int num_sp  = m_num_species;
    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
    const ParticleReal inv_bin_size
        = 1.0 / (  static_cast<ParticleReal>(m_coalescence_bin_size[0])
                 * static_cast<ParticleReal>(m_coalescence_bin_size[1])
                 * static_cast<ParticleReal>(m_coalescence_bin_size[2]) );
    const ParticleReal inv_bin_volume = inv_cell_volume*inv_bin_size;

    Real num_collisions = 0;
    const auto& gvec = a_temperature.nGrowVect();

    auto kernel_choice = m_coalescence_kernel;
    auto include_brownian_coalescence = m_include_brownian_coalescence;

    Real mcshuffle_wtime_sec = 0.0;
    Real mcpairing_wtime_sec = 0.0;
    Real coalescence_wtime_sec = 0.0;

// Do NOT add OpenMP here; building DenseBins is not thread-safe.
    for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti) {

        auto& ptile = ParticlesAt(a_lev, pti);
        auto& aos = ptile.GetArrayOfStructs();
        auto& soa = ptile.GetStructOfArrays();
        const size_t np = aos.numParticles();
        auto pstruct_ptr = aos().dataPtr();

        /* SoA attributes */
        auto* mass_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::mass).data();
        Array<ParticleReal*,AMREX_SPACEDIM> v_ptr;
        v_ptr[0] = soa.GetRealData(SuperDropletsRealIdxSoA::vx).data();
        v_ptr[1] = soa.GetRealData(SuperDropletsRealIdxSoA::vy).data();
        v_ptr[2] = soa.GetRealData(SuperDropletsRealIdxSoA::vz).data();

        /* Runtime-added SoA attributes */
        int rtoff_i = SuperDropletsIntIdxSoA::ncomps;
        auto* active_ptr = soa.GetIntData(rtoff_i+SuperDropletsIntIdxSoA_RT::active).data();
        int rtoff_r = SuperDropletsRealIdxSoA::ncomps;
        auto* radius_ptr = soa.GetRealData(rtoff_r+SuperDropletsRealIdxSoA_RT::radius).data();
        auto* mult_ptr = soa.GetRealData(rtoff_r+SuperDropletsRealIdxSoA_RT::multiplicity).data();
        auto* vterm_ptr = soa.GetRealData(rtoff_r+SuperDropletsRealIdxSoA_RT::term_vel).data();
        auto* Tfz_ptr = soa.GetRealData(idx_ice_Tfz(num_ae,num_sp)).data();
        auto* a_ptr = soa.GetRealData(idx_ice_a(num_ae,num_sp)).data();
        auto* c_ptr = soa.GetRealData(idx_ice_c(num_ae,num_sp)).data();
        auto* mrime_ptr = soa.GetRealData(idx_ice_mrime(num_ae,num_sp)).data();
        auto* nmono_ptr = soa.GetRealData(idx_ice_nmono(num_ae,num_sp)).data();

        /* species and aerosol masses */
        SDSpeciesMassArr sp_mass_ptrs;
        SDAerosolMassArr ae_mass_ptrs;
        setupMassPointers(soa, sp_mass_ptrs, ae_mass_ptrs);

        // Get pointers to persistent device data
        const ParticleReal* sp_rho_arr = getSpeciesDensitiesDevice();
        const int* sp_sol_arr = getSpeciesSolubilitiesDevice();
        const ParticleReal* ae_rho_arr = getAerosolDensitiesDevice();
        const int* ae_sol_arr = getAerosolSolubilitiesDevice();

        int grid = pti.index();
        Box box = a_temperature[grid].box(); box.grow(-gvec);
        int ntiles = numTilesInBox(box, true, m_coalescence_bin_size);
        auto binner = GetParticleBin{plo, dxi, domain, m_coalescence_bin_size, box};
        DenseBins<ParticleType> bins;
        bins.build( np, pstruct_ptr, ntiles, binner);
        AMREX_ALWAYS_ASSERT(np == static_cast<size_t>(bins.numItems()));
        AMREX_ALWAYS_ASSERT(bins.numBins() >= 0);
        auto inds = bins.permutationPtr();
        auto offsets = bins.offsetsPtr();

#ifndef _WIN32
        struct timeval mcshuffle_start, mcshuffle_end;
        gettimeofday(&mcshuffle_start, NULL);
#endif

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

#ifndef _WIN32
        gettimeofday(&mcshuffle_end,NULL);
        long long mcshuffle_wtime;
        mcshuffle_wtime = (   (mcshuffle_end.tv_sec   * 1000000 + mcshuffle_end.tv_usec  )
                            - (mcshuffle_start.tv_sec * 1000000 + mcshuffle_start.tv_usec) );
        mcshuffle_wtime_sec += (double) mcshuffle_wtime / 1000000.0;
#endif

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
            coll_rate_ptr[i] = 0.0;
            coll_rmndr_ptr[i] = 0.0;
        });
        Gpu::synchronize();

#ifndef _WIN32
        struct timeval mcpairing_start, mcpairing_end;
        gettimeofday(&mcpairing_start, NULL);
#endif

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
                if (active_ptr[pi] == 0) { continue; }
                if (active_ptr[pj] == 0) { continue; }
                if (mult_ptr[pi] == 0) { continue; }
                if (mult_ptr[pj] == 0) { continue; }

                np_bin_ptr[pi] = np_bin_ptr[pj] = np_bin;
                partner_idx_ptr[pi] = pj;
                partner_idx_ptr[pj] = pi;

                int i = -1, j = -1;
                if (mult_ptr[pi] >= mult_ptr[pj]) { i = pi; j = pj; }
                else                              { i = pj; j = pi; }
                flag_prey_ptr[i] = 1;
                flag_prey_ptr[j] = 0;
            }
        } );
        Gpu::synchronize();

#ifndef _WIN32
        gettimeofday(&mcpairing_end,NULL);
        long long mcpairing_wtime;
        mcpairing_wtime = (   (mcpairing_end.tv_sec   * 1000000 + mcpairing_end.tv_usec  )
                            - (mcpairing_start.tv_sec * 1000000 + mcpairing_start.tv_usec) );
        mcpairing_wtime_sec += (double) mcpairing_wtime / 1000000.0;

        struct timeval coalescence_start, coalescence_end;
        gettimeofday(&coalescence_start, NULL);
#endif

        // We've already declared these pointers above, no need to redeclare

        // calculate collision efficiencies for each pair
        ParallelForRNG( np, [=] AMREX_GPU_DEVICE (int i, RandomEngine const& rnd_eng) noexcept
        {
            if (partner_idx_ptr[i] < 0) { return; }
            if (!flag_prey_ptr[i]) { return; }

            int pi = i; // prey - higher multiplicity
            int pj = partner_idx_ptr[i]; // predator - lower multiplicity
            AMREX_ALWAYS_ASSERT(mult_ptr[pi] >= mult_ptr[pj]);

            // get phases for the two particles
            auto phase_i = SD_phase(pi, sp_idx_w, sp_idx_i, sp_mass_ptrs);
            auto phase_j = SD_phase(pj, sp_idx_w, sp_idx_i, sp_mass_ptrs);

            ParticleReal k_val = 0.0;
            if ((phase_i == SDPhase::water) && (phase_j == SDPhase::water)) {

                // coalescence between two water droplets

                if (kernel_choice == SDCoalescenceKernelType::golovin) {

                    k_val = ckernel.golovin(radius_ptr[pi],radius_ptr[pj]);

                } else {

                    ParticleReal v_i[AMREX_SPACEDIM], v_j[AMREX_SPACEDIM];
                    for (int d = 0; d < AMREX_SPACEDIM; d++) {
                        v_i[d] = v_ptr[d][pi];
                        v_j[d] = v_ptr[d][pj];
                    }
                    v_i[AMREX_SPACEDIM-1] -= vterm_ptr[pi];
                    v_j[AMREX_SPACEDIM-1] -= vterm_ptr[pj];

                    if (kernel_choice == SDCoalescenceKernelType::sedimentation) {
                        k_val = ckernel.sedimentation(radius_ptr[pi],radius_ptr[pj],v_i,v_j);
                    } else if (kernel_choice == SDCoalescenceKernelType::Longs) {
                        k_val = ckernel.Longs(radius_ptr[pi],radius_ptr[pj],v_i,v_j);
                    } else if (kernel_choice == SDCoalescenceKernelType::Halls) {
                        k_val = ckernel.Halls(radius_ptr[pi],radius_ptr[pj],v_i,v_j);
                    }

                    if (k_val < 0.0) {
                        amrex::Abort("Invalid value for k_val");
                    }
                }

            } else if ((phase_i == SDPhase::ice) && (phase_j == SDPhase::ice)) {

                // aggregation between two ice particles

                // ice particle 1
                auto vz_i = v_ptr[AMREX_SPACEDIM-1][pi] - vterm_ptr[pi];
                auto a_i = a_ptr[pi];
                auto c_i = c_ptr[pi];
                auto rhoi_i = ice_rho(a_i, c_i, sp_mass_ptrs[sp_idx_i][pi]);
                auto maxR_i = std::max(a_i, c_i);
                auto k_i = std::exp(-ckernel.k_coeff * c_i/a_i);
                auto area_i = PI * a_i * maxR_i * std::exp(k_i*std::log(rhoi_i/rho_ice));

                // ice particle 2
                auto vz_j = v_ptr[AMREX_SPACEDIM-1][pj] - vterm_ptr[pj];
                auto a_j = a_ptr[pj];
                auto c_j = c_ptr[pj];
                auto rhoi_j = ice_rho(a_j, c_j, sp_mass_ptrs[sp_idx_i][pj]);
                auto maxR_j = std::max(a_j, c_j);
                auto k_j = std::exp(-ckernel.k_coeff * c_j/a_j);
                auto area_j = PI * a_j * maxR_j * std::exp(k_j*std::log(rhoi_j/rho_ice));

                // velocity difference
                auto dvz = std::sqrt((vz_i-vz_j)*(vz_i-vz_j));

                // collision efficiency
                k_val = 0.1 * (std::sqrt(area_i)+std::sqrt(area_j))*(std::sqrt(area_i)+std::sqrt(area_j)) * dvz;
                // if (std::min(rhoi_i, rhoi_j) < 10.0) { k_val = 0.0; }

            } else {

                // riming between a water droplet and an ice particle
                AMREX_ALWAYS_ASSERT(    ((phase_i == SDPhase::ice) && (phase_j == SDPhase::water))
                                     || ((phase_i == SDPhase::water) && (phase_j == SDPhase::ice)) );

                // figure out which one is water, and which is ice
                int id_w = -1, id_i = -1;
                if (phase_i == SDPhase::water) { id_w = pi; id_i = pj; }
                else { id_w = pj; id_i = pi; }

                // water droplet
                auto vz_w = v_ptr[AMREX_SPACEDIM-1][id_w] - vterm_ptr[id_w];
                auto r_w = radius_ptr[id_w];

                // ice particle
                auto vz_i = v_ptr[AMREX_SPACEDIM-1][id_i] - vterm_ptr[id_i];
                auto a_i = a_ptr[id_i];
                auto c_i = c_ptr[id_i];
                auto rhoi_i = ice_rho(a_i, c_i, sp_mass_ptrs[sp_idx_i][id_i]);
                auto eqr_i = std::exp((1.0/3.0)*std::log(a_i*a_i*c_i));
                auto maxR_i = std::max(a_i, c_i);

                // interpolate flow quantities at ice particle location
                ParticleType& p_ice = pstruct_ptr[id_i];
                ParticleReal temperature, moist_density;
                // Define field values array to store interpolated results
                ParticleReal field_values[2]; // temperature, moist_density

                // Define array of field arrays to interpolate from
                const Array4<const Real> field_arrays[2] = {
                    temperature_arr,
                    moist_density_arr
                };

                // Use the interpolation helper function to interpolate all fields at once
                ERF::Interpolation::interpolateFields(
                    p_ice, plo, dxi, field_arrays, field_values, 2,
                    is_periodic_z ? 1 : 0, is_periodic_z ? nullptr : &zheight
                );

                // Extract interpolated values
                temperature = field_values[0];
                moist_density = field_values[1];

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
                    auto area_i = area_i_ce * std::exp(k_i*std::log(rhoi_i/rho_ice));

                    // collision-riming kernel: E_rime*A_g*|vj-vk|
                    k_val *= ( (PI*(r_w+a_i)*(r_w+maxR_i) - (area_i_ce-area_i)) * dvz);

                }

            }

            if (include_brownian_coalescence) {

                ParticleReal pressure = 0.0, temperature = 0.0;
                {
                    ParticleType& par_1 = pstruct_ptr[pi];
                    auto iv = getParticleCell(par_1, plo, dxi, domain);
                    pressure = pressure_arr(iv[0],iv[1],iv[2],0);
                    temperature = temperature_arr(iv[0],iv[1],iv[2],0);
                }

                ParticleReal sd_mass_1 = 0.0,
                             sd_mass_2 = 0.0;
                for (int ia = 0; ia < num_ae; ia++) {
                    sd_mass_1 += ae_mass_ptrs[ia][pi];
                    sd_mass_2 += ae_mass_ptrs[ia][pj];
                }
                for (int ia = 0; ia < num_sp; ia++) {
                    sd_mass_1 += sp_mass_ptrs[ia][pi];
                    sd_mass_2 += sp_mass_ptrs[ia][pj];
                }

                auto r_eff_1 = SD_effective_radius( pi, sp_idx_w,
                                                    rho_water,
                                                    num_sp, num_ae,
                                                    sp_sol_arr, ae_sol_arr,
                                                    sp_mass_ptrs, ae_mass_ptrs,
                                                    sp_rho_arr, ae_rho_arr );
                auto r_eff_2 = SD_effective_radius( pj, sp_idx_w,
                                                    rho_water,
                                                    num_sp, num_ae,
                                                    sp_sol_arr, ae_sol_arr,
                                                    sp_mass_ptrs, ae_mass_ptrs,
                                                    sp_rho_arr, ae_rho_arr );

                auto k_brown = ckernel.Brownian_SeinfeldPandis( r_eff_1,
                                                                r_eff_2,
                                                                sd_mass_1,
                                                                sd_mass_2,
                                                                pressure,
                                                                temperature );
                if (k_brown < 0.0) {
                    amrex::Abort("Invalid value for k_brown");
                }

                k_val += k_brown;
            }


            auto prob_ij = k_val*inv_bin_volume;
            auto prob_sd_ij = std::max(mult_ptr[pi],mult_ptr[pj])*prob_ij;

            auto ns = static_cast<ParticleReal>(np_bin_ptr[i]);
            auto scaling_factor = 0.5*ns*(ns-1)/std::floor(0.5*ns);
            auto scaled_prob = prob_sd_ij * scaling_factor;

            auto gamma = coalescence_rate ( rnd_eng, (scaled_prob*a_dt) );
            if (gamma > 0) {
                amrex::Gpu::Atomic::Add(particle_collisions_ptr, gamma);
                coll_rate_ptr[pi] = std::min(gamma,std::floor(mult_ptr[pi]/mult_ptr[pj]));
                coll_rate_ptr[pj] = coll_rate_ptr[pi];
                coll_rmndr_ptr[pi] = mult_ptr[pi] - coll_rate_ptr[pi]*mult_ptr[pj];
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
            auto phase_i = SD_phase(i, sp_idx_w, sp_idx_i, sp_mass_ptrs);
            auto phase_j = SD_phase(j, sp_idx_w, sp_idx_i, sp_mass_ptrs);

            if ((phase_i == SDPhase::water) && (phase_j == SDPhase::water)) {

                // coalescence between two water droplets
                coal_update_attribs( i, j,
                                     flag_prey_ptr,
                                     coll_rate_ptr,
                                     coll_rmndr_ptr,
                                     mass_ptr,
                                     radius_ptr,
                                     Tfz_ptr,
                                     mult_ptr,
                                     num_sp,
                                     sp_mass_ptrs,
                                     num_ae,
                                     ae_mass_ptrs );

            } else if ((phase_i == SDPhase::ice) && (phase_j == SDPhase::ice)) {

                // aggregation between two ice particles
                aggr_update_attribs( i, j,
                                     sp_idx_i,
                                     rho_ice,
                                     10.0,
                                     flag_prey_ptr,
                                     coll_rate_ptr,
                                     coll_rmndr_ptr,
                                     Tfz_ptr,
                                     a_ptr,
                                     c_ptr,
                                     mrime_ptr,
                                     nmono_ptr,
                                     mult_ptr,
                                     num_sp,
                                     sp_mass_ptrs,
                                     num_ae,
                                     ae_mass_ptrs );

            } else {

                // riming between a water droplet and an ice particle
                AMREX_ALWAYS_ASSERT(    ((phase_i == SDPhase::ice) && (phase_j == SDPhase::water))
                                     || ((phase_i == SDPhase::water) && (phase_j == SDPhase::ice)) );

                ParticleType& p_ice = pstruct_ptr[i];
                ParticleReal temperature, pressure, moist_density, qv;
                // Define field values array to store interpolated results
                ParticleReal field_values[4]; // temperature, pressure, moist_density, qv

                // Define array of field arrays to interpolate from
                const Array4<const Real> field_arrays[4] = {
                    temperature_arr,
                    pressure_arr,
                    moist_density_arr,
                    qv_arr
                };

                // Use the interpolation helper function to interpolate all fields at once
                ERF::Interpolation::interpolateFields(
                    p_ice, plo, dxi, field_arrays, field_values, 4,
                    is_periodic_z ? 1 : 0, is_periodic_z ? nullptr : &zheight
                );

                // Extract interpolated values
                temperature = field_values[0];
                pressure = field_values[1];
                moist_density = field_values[2];
                qv = field_values[3];
                auto coeff_moldiff = mat_prop.coeffMolecularDiffusion(temperature, pressure);

                rime_update_attribs( i, j,
                                     sp_idx_w, sp_idx_i,
                                     rho_water, rho_ice,
                                     flag_prey_ptr,
                                     coll_rate_ptr,
                                     coll_rmndr_ptr,
                                     phase_i, phase_j,
                                     v_ptr,
                                     vterm_ptr,
                                     radius_ptr,
                                     Tfz_ptr,
                                     a_ptr,
                                     c_ptr,
                                     mrime_ptr,
                                     nmono_ptr,
                                     mult_ptr,
                                     temperature,
                                     moist_density,
                                     pressure,
                                     qv,
                                     coeff_moldiff,
                                     dmdt,
                                     num_sp,
                                     sp_mass_ptrs,
                                     num_ae,
                                     ae_mass_ptrs );
            }

        } );
        Gpu::synchronize();

        // update particle radius and total mass
        ParallelFor( np, [=] AMREX_GPU_DEVICE (int i)
        {
            SuperDropletPC::updateParticleAttributes(
                i, radius_ptr, mass_ptr, sp_idx_w, rho_water,
                num_sp, num_ae, sp_sol_arr, ae_sol_arr,
                sp_mass_ptrs, ae_mass_ptrs, sp_rho_arr, ae_rho_arr);
        } );
        Gpu::synchronize();

#ifndef _WIN32
        gettimeofday(&coalescence_end,NULL);
        long long coalescence_wtime;
        coalescence_wtime = (  (coalescence_end.tv_sec   * 1000000 + coalescence_end.tv_usec  )
                             - (coalescence_start.tv_sec * 1000000 + coalescence_start.tv_usec) );
        coalescence_wtime_sec += (double) coalescence_wtime / 1000000.0;
#endif
    }

    ParallelDescriptor::ReduceRealSum(  &num_collisions,
                                        1,
                                        ParallelDescriptor::IOProcessorNumber() );

#ifndef _WIN32
    gettimeofday(&total_end,NULL);
    long long total_wtime;
    total_wtime = (   (total_end.tv_sec   * 1000000 + total_end.tv_usec  )
                   -  (total_start.tv_sec * 1000000 + total_start.tv_usec) );
    Real total_wtime_sec = (double) total_wtime / 1000000.0;

    ParallelDescriptor::ReduceRealMax( &mcshuffle_wtime_sec,
                                       1,
                                       ParallelDescriptor::IOProcessorNumber() );
    ParallelDescriptor::ReduceRealMax( &mcpairing_wtime_sec,
                                       1,
                                       ParallelDescriptor::IOProcessorNumber() );
    ParallelDescriptor::ReduceRealMax( &coalescence_wtime_sec,
                                       1,
                                       ParallelDescriptor::IOProcessorNumber() );
    ParallelDescriptor::ReduceRealMax( &total_wtime_sec,
                                       1,
                                       ParallelDescriptor::IOProcessorNumber() );
#else
    Real total_wtime_sec = 0.0;
#endif

    Print() << "SuperDropletPC(" << m_name << "): "
            << "number of collisions = " << num_collisions << "\n"
            << "    "
            << "wall time (seconds) = " << total_wtime_sec << " (total), "
            << mcshuffle_wtime_sec << " (MC shuffle), "
            << mcpairing_wtime_sec << " (MC pairing), "
            << coalescence_wtime_sec << " (coalescence)"
            << "\n";
}

#endif

