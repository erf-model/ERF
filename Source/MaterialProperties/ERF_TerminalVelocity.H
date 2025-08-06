#ifndef TERMINALVELOCITY_H
#define TERMINALVELOCITY_H

#include <cmath>
#include <AMReX_Gpu.H>
#include "ERF_Constants.H"

template <typename RT /*!< real-type */>
struct TerminalVelocity
{
    RT m_rho_w; /*!< density of water */
    RT m_rho_i; /*!< density of ice */

    /*! \brief Compute viscosity coefficient - Pruppacher & Klett, 1997*/
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto viscCoeff( const RT a_T /*!< temperature [K] */ ) const noexcept
    {
        RT T_degC = a_T - 273.15; // [K] -> [deg Celsius]
        RT viscosity = 0.0;
        if (T_degC >= 0) {
            viscosity = (1.7180 + 4.9e-3*T_degC) * 1.0e-4;
        } else {
            viscosity = (  1.7180 + 4.9e-3*T_degC
                         - 1.2e-5*T_degC*T_degC     ) * 1.0e-4;
        }
        return viscosity;
    }

    /*! \brief Terminal velocity - Rogers & Yau, 1989 */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT RogersYau( const RT a_r) const noexcept
    {
        static const RT k1 = 1.233e8; // m^{-1} s^{-1}
        return (k1 * a_r * a_r);
    };

    /*! \brief Terminal velocity - Atlas & Ulbrich, 1977 */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT AtlasUlbrich( const RT a_r) const noexcept
    {
        static const RT a = 3.778;
        static const RT b = 0.67;
        RT D = 2*a_r*1000; // meter -> mm
        return a*(std::exp(b*std::log(D)));
    };

    /*! \brief Terminal velocity - cloud/rain droplets
        Shima et al., 2009, "The super-droplet method for the numerical simulation of clouds and
        precipitation: A particle-based and probabilistic microphysics model coupled with a non-
        hydrostatic model." Quart. J. Roy. Meteorol. Soc., 135: 1307-1320
        https://github.com/Shima-Lab/SCALE-SDM_BOMEX_Sato2018/blob/master/contrib/SDM/sdm_motion.f90 */
    AMREX_GPU_HOST_DEVICE
    RT CloudRainShima ( const RT a_r,      /*!< radius */
                        const RT a_rho,    /*!< density */
                        const RT a_P,      /*!< pressure */
                        const RT a_T       /*!< temperature */) const noexcept
    {
        RT lb4l = 6.62E-6;    // base length[cm] for mean free path  of air molecules
        RT visb4l = 1.818E-4; // base dynamic viscosity[g/(cm*s)] for mean free path of air molecules
        RT pb4l = 1013.25;    // base pressure[hPa] for mean free path of air molecules
        RT tb4l = 293.15;     // base temperarure[20C] for mean free path of air molecules

        RT P_hPa = a_P / 100.0;           // [Pa] -> [hPa]
        RT diameter_cm = 2.0*a_r*100.0;   // [m] -> [cm]
        RT rho_mat_cgs = m_rho_w/1000.0;    // [kg/m^3] -> [g/cm^3]
        RT rho_air_cgs = a_rho/1000.0;    // [kg/m^3] -> [g/cm^3]
        RT grav_cgs = CONST_GRAV * 100;   // [m/s^2] -> [cm/s^2]

        RT gxdrow = grav_cgs * ( rho_mat_cgs - rho_air_cgs );
        RT viscosity = viscCoeff(a_T);

        RT v_terminal = 0.0;
        if (diameter_cm < 1.9e-3) {
            // small cloud droplets
            RT sd_l = lb4l * (viscosity/visb4l) * (pb4l/P_hPa) * std::sqrt(a_T/tb4l);
            RT c1  = gxdrow / (18.0*viscosity);
            RT csc = 1.0 + 2.510 * (sd_l/diameter_cm);
            v_terminal = c1 * csc * diameter_cm * diameter_cm;
        } else if (diameter_cm < 1.07e-1) {
            // large cloud droplets and small raindrops
            RT sd_l = lb4l * (viscosity/visb4l) * (pb4l/P_hPa) * std::sqrt(a_T/tb4l);
            RT c2 = rho_air_cgs * (4.0*gxdrow)/(3.0*viscosity*viscosity);
            RT nda = c2 * diameter_cm * diameter_cm * diameter_cm;
            RT csc = 1.0 + 2.510 * ( sd_l / diameter_cm );
            RT nre = csc * std::exp(vz_b_x(std::log(nda)));
            v_terminal = (viscosity*nre)/(rho_air_cgs*diameter_cm);
        } else {
            // large raindrops
            RT tau = 1.0 - a_T/647.0960;
            RT sigma = 0.2358 * std::exp( 1.2560*std::log(tau) )
                                * ( 1.0 - 0.6250*tau );
            if( a_T<267.5 ) {
                tau = std::tanh( (a_T-243.9)/35.35 );
                sigma = sigma - 2.854E-3 * tau + 1.666E-3;
            }
            sigma = sigma * 1.E+3; // N/m => g/s^2
            RT c3 = (4.0*gxdrow)/(3.0*sigma);
            RT bond = c3 * diameter_cm * diameter_cm;
            RT npp =    (rho_air_cgs*rho_air_cgs) * (sigma*sigma*sigma)
                        / (gxdrow*viscosity*viscosity*viscosity*viscosity);
            npp = std::exp( std::log(npp)/6.0 );
            RT nre = npp * std::exp(vz_c_x(std::log(bond*npp)));
            v_terminal = (viscosity*nre)/(rho_air_cgs*diameter_cm);
        }

        v_terminal /= 100.0; // [cm/s] -> [m/s]
        return v_terminal;
    }

    /*! \brief Compute vz_b(x) */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT vz_b_x ( const RT a_x /*!< x */) const noexcept
    {
        RT x1 = a_x;
        RT x2 = x1 * x1;
        RT x3 = x1 * x2;
        RT y = - 0.318657E+1
               + 0.9926960    * x1
               - 0.153193E-2  * x2
               - 0.987059E-3  * x3
               - 0.578878e-3  * x2*x2
               + 0.855176E-4  * x3*x2
               - 0.327815E-5  * x3*x3;
        return y;
    }

    /*! \brief Compute vc_b(x) */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT vz_c_x ( const RT a_x /*!< x */) const noexcept
    {
        RT x1 = a_x;
        RT x2 = x1 * x1;
        RT x3 = x1 * x2;
        RT y = - 0.500015E+1
               + 0.523778E+1  * x1
               - 0.204914E+1  * x2
               + 0.47529400   * x3
               - 0.542819E-1  * x2*x2
               + 0.238449E-2  * x3*x2;
        return y;
    }

    /*! \brief Terminal velocity - ice particles (Bohm's theory)
     *  Ref: Shima et al., 2020, "Predicting the morphology of ice particles
     *  in deep convection using the super-droplet method: development and
     *  evaluation of SCALE-SDM 0.2.5-2.2.0, -2.2.1, and -2.2.2.",
     *  Geosci. Model Dev., 13, 4107–4157, 2020
     *
     *  and
     *
     *  Code: https://doi.org/10.5281/zenodo.3483650
     *  File: [...]/contrib/SDM/sdm_motion.f90 */
    AMREX_GPU_HOST_DEVICE
    RT IceBohm ( const RT a_m,      /*!< mass */
                 const RT a_a,      /*!< ice particle equatorial radius */
                 const RT a_c,      /*!< ice particle polar radius */
                 const RT a_rhoi,   /*!< ice particle density */
                 const RT a_rho,    /*!< density */
                 const RT a_T       /*!< temperature */) const noexcept
    {
        RT mu = viscCoeff(a_T);

        auto phi = a_c / a_a; // aspect ratio
        auto d = 2.0 * a_a; // characteristic length
        auto q = std::exp( std::exp(-phi) * std::log(a_rhoi/m_rho_i) ); // area ratio

        RT X0 = 2.8e06;
        auto t3 = PI * mu*mu * std::max(phi,1.0) * std::max(std::sqrt(std::sqrt(q)), q);
        auto X = (8.0*a_m*CONST_GRAV*a_rho) / t3;
        auto Xp = X * (1.0 + (X/X0)*(X/X0)) / (1.0 + 1.6*(X/X0)*(X/X0));

        auto k = std::min( std::max(RT(0.82)+0.18*phi, RT(0.85)),
                           std::min( RT(0.37) + 0.63/std::sqrt(phi),
                                     RT(1.33)/(std::max(std::log(phi),0.0)+RT(1.19)) ));

        auto phi_poly = 3.76-8.41*phi+9.18*phi*phi-3.53*phi*phi*phi;
        auto TVGamma = std::max( RT(1.0), std::min(RT(1.98), phi_poly) );

        auto C_DPS = std::max(0.292*k*TVGamma, 0.492 - 0.2/std::sqrt(phi));
        auto C_DP  = std::max( RT(1.0), q*(1.46*q-0.46) ) * C_DPS;
        auto C_DO = 4.5 * k*k * std::max(phi, 1.0);

        auto t2 = RT(1.0) + C_DP/(6.0*k) * std::sqrt(Xp/C_DP);
        auto beta = std::sqrt(t2) - 1.0;
        auto TVgamma = (C_DO-C_DP) / (4*C_DP);
        auto t1 = RT(1.0) + (2*beta*std::exp(-beta*TVgamma)) / ((2.0+beta)*(1.0+beta));
        auto N_Re = (6.0*k/C_DP) * beta*beta * t1;

        auto vterm = (mu*N_Re)/(a_rho*d);
        return vterm;
    }

};

#endif
