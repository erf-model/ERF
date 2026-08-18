#ifndef SUPERDROPLET_PC_RIMING_H_
#define SUPERDROPLET_PC_RIMING_H_

#include <cmath>
#include <algorithm>

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>

#include "ERF_Constants.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_SuperDropletPCMassChange.H"

#ifdef ERF_USE_PARTICLES

/*! \brief Riming microphysics helpers (dynamic viscosity, Rasmussen-Heymsfield
 *  impact velocity, ice surface temperature, and Heymsfield-Pflaum rime density).
 *  Factored out of ERF_SuperDropletPCCoalescence.cpp into this header, and made
 *  host+device, so they can be exercised by unit tests. */
namespace SDRiming
{
    using amrex::Real;
    using amrex::ParticleReal;

    /*! \brief Compute dynamic viscosity */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto viscCoeff ( const ParticleReal a_T /*!< temperature */ )
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
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto impactVelocity_RasmussenHeymsfield1985( const ParticleReal a_Re, /*!< Reynolds number */
                                                 const ParticleReal a_St  /*!< Stokes number */ )
    {
        // the fits only apply for St >= 0.1; log10(0) raises FE_DIVBYZERO under fpe_trap_zero
        ParticleReal w(zero), w2(zero), w3(zero), w4(zero);
        if (a_St >= ParticleReal(0.1)) {
            w  = std::log10(a_St);
            w2 = w*w;
            w3 = w2*w;
            w4 = w3*w;
        }
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
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ParticleReal iceSurfaceTemperature( const ParticleReal a_T, /*!< temperature */
                                        const ParticleReal a_P, /*!< pressure */
                                        const ParticleReal a_qv, /*!< vapour fraction */
                                        const ParticleReal a_D, /*!< diffusivity coeff */
                                        const SDMassChangeUtils_SV::dMdt<ParticleReal>& a_dmdt /*!< mass change utilities */)
    {
        // erf_qsati takes pressure in mbar and erf_esati returns mbar; a_P arrives in Pa
        Real qsat_r = Real(zero); erf_qsati(a_T, a_P/Real(100.0), qsat_r);
        ParticleReal qsat = ParticleReal(qsat_r);
        auto sup_sat  = a_qv/qsat - ParticleReal(one);
        auto drho  = sup_sat / (a_D * a_dmdt.Fk_plus_Fd(a_T, ParticleReal(erf_esati(a_T)*Real(100.0)), a_D));
        auto dT = (a_dmdt.L * a_D / a_dmdt.K) * drho;
        return a_T + dT;
    }

    /*! \brief Rime density parameterization from Heymsfield and Pflaum (1985)
     *
     *  Computes the density of rime accreted onto ice particles based on the
     *  dimensionless impact parameter Y = -r[um] * v_impact / T_surf[degC].
     *
     *  Reference: Heymsfield, A. J., and J. C. Pflaum, 1985: A quantitative
     *  assessment of the accuracy of techniques for calculating graupel growth.
     *  J. Atmos. Sci., 42, 2264-2274.
     *
     *  Two regimes based on surface temperature and Y parameter:
     *  - Low-density (T_surf <= -5 degC or Y < 1.6): rho = 0.30*Y^0.44  [Eq. 8]
     *  - High-density (otherwise): rho = exp(-0.03115 - 1.703*Y + 0.9116*Y^2 - 0.1224*Y^3) [Eq. 9]
     *
     *  Output bounded to 0.1 <= rho_rime <= 0.91 g/cm^3 (100-910 kg/m^3)
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto rimeDensity_HeymsfieldPflaum1985( const ParticleReal a_radius,      /*!< droplet radius [m] */
                                           const ParticleReal a_a,           /*!< ice particle equatorial radius [m] */
                                           const ParticleReal a_c,           /*!< ice particle polar radius [m] */
                                           const ParticleReal a_vz_w,        /*!< vertical velocity of water droplet [m/s] */
                                           const ParticleReal a_vz_i,        /*!< vertical velocity of ice particle [m/s] */
                                           const ParticleReal a_rho_w,       /*!< water density [kg/m^3] */
                                           const ParticleReal a_D,           /*!< diffusivity coefficient [m^2/s] */
                                           const SDMassChangeUtils_SV::dMdt<ParticleReal>& a_dmdt, /*!< mass change utilities */
                                           const ParticleReal a_T,           /*!< temperature [K] */
                                           const ParticleReal a_rhom,        /*!< moist air density [kg/m^3] */
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

} // namespace SDRiming

#endif
#endif
