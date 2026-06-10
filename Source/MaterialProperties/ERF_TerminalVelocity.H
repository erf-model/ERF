#ifndef TERMINALVELOCITY_H
#define TERMINALVELOCITY_H

#include <cmath>
#include <AMReX_Gpu.H>
#include "ERF_Constants.H"

template <typename RT /*!< real-type */>
struct TerminalVelocity
{
    RT m_rho_w; /*!< density of water */

    /*! \brief Compute viscosity coefficient - Pruppacher & Klett, 1997*/
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto viscCoeff( const RT a_T /*!< temperature [K] */ ) const noexcept
    {
        RT T_degC = a_T - RT(273.15); // [K] -> [deg Celsius]
        RT viscosity = RT(0.0);
        if (T_degC >= RT(0)) {
            viscosity = (RT(1.7180) + RT(4.9e-3)*T_degC) * RT(1.0e-4);
        } else {
            viscosity = (  RT(1.7180) + RT(4.9e-3)*T_degC
                         - RT(1.2e-5)*T_degC*T_degC     ) * RT(1.0e-4);
        }
        return viscosity;
    }

    /*! \brief Terminal velocity - Rogers & Yau, 1989 */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT RogersYau( const RT a_r) const noexcept
    {
        static const RT k1 = RT(1.233e8); // m^{-1} s^{-1}
        return (k1 * a_r * a_r);
    }

    /*! \brief Terminal velocity - Atlas & Ulbrich, 1977 */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT AtlasUlbrich( const RT a_r) const noexcept
    {
        static const RT a = RT(3.778);
        static const RT b = RT(0.67);
        RT D = RT(2)*a_r*RT(1000); // meter -> mm
        return a*(std::exp(b*std::log(D)));
    }

    /* The following code is adapted from SCALE-SDM:
     * Copyright (c) 2012-2014, Team SCALE, All rights reserved. */
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
        RT lb4l = RT(6.62E-6);    // base length[cm] for mean free path  of air molecules
        RT visb4l = RT(1.818E-4); // base dynamic viscosity[g/(cm*s)] for mean free path of air molecules
        RT pb4l = RT(1013.25);    // base pressure[hPa] for mean free path of air molecules
        RT tb4l = RT(293.15);     // base temperature[20C] for mean free path of air molecules

        RT P_hPa = a_P / RT(100.0);           // [Pa] -> [hPa]
        RT diameter_cm = RT(2.0)*a_r*RT(100.0);   // [m] -> [cm]
        RT rho_mat_cgs = m_rho_w/RT(1000.0);  // [kg/m^3] -> [g/cm^3]
        RT rho_air_cgs = a_rho/RT(1000.0);    // [kg/m^3] -> [g/cm^3]
        RT grav_cgs = CONST_GRAV * RT(100);   // [m/s^2] -> [cm/s^2]

        RT gxdrow = grav_cgs * ( rho_mat_cgs - rho_air_cgs );

        RT viscosity = viscCoeff(a_T);

        RT v_terminal = RT(0.0);
        if (diameter_cm < RT(1.9e-3)) {
            // small cloud droplets
            RT sd_l = lb4l * (viscosity/visb4l) * (pb4l/P_hPa) * std::sqrt(a_T/tb4l);
            RT c1  = gxdrow / (RT(18.0)*viscosity);
            RT csc = RT(1.0) + RT(2.510) * (sd_l/diameter_cm);
            v_terminal = c1 * csc * diameter_cm * diameter_cm;
        } else if (diameter_cm < RT(1.07e-1)) {
            // large cloud droplets and small raindrops
            RT sd_l = lb4l * (viscosity/visb4l) * (pb4l/P_hPa) * std::sqrt(a_T/tb4l);
            RT c2 = rho_air_cgs * (RT(4.0)*gxdrow)/(RT(three)*viscosity*viscosity);
            RT nda = c2 * diameter_cm * diameter_cm * diameter_cm;
            RT csc = RT(1.0) + RT(2.510) * ( sd_l / diameter_cm );
            RT nre = csc * std::exp(vz_b_x(std::log(nda)));
            v_terminal = (viscosity*nre)/(rho_air_cgs*diameter_cm);
        } else {
            // large raindrops
            RT tau = RT(1.0) - a_T/RT(647.0960);
            RT sigma = RT(0.2358) * std::exp( RT(1.2560)*std::log(tau) )
                                * ( RT(1.0) - RT(0.6250)*tau );
            if( a_T<RT(267.5) ) {
                tau = std::tanh( (a_T-RT(243.9))/RT(35.35) );
                sigma = sigma - RT(2.854E-3) * tau + RT(1.666E-3);
            }
            sigma = sigma * RT(1.E+3); // N/m => g/s^2
            RT c3 = (RT(4.0)*gxdrow)/(RT(three)*sigma);
            RT bond = c3 * diameter_cm * diameter_cm;
            RT npp =    (rho_air_cgs*rho_air_cgs) * (sigma*sigma*sigma)
                        / (gxdrow*viscosity*viscosity*viscosity*viscosity);
            npp = std::exp( std::log(npp)/RT(6.0) );
            RT nre = npp * std::exp(vz_c_x(std::log(bond*npp)));
            v_terminal = (viscosity*nre)/(rho_air_cgs*diameter_cm);
        }

        v_terminal /= RT(100.0); // [cm/s] -> [m/s]
        return v_terminal;
    }

    /* The following code is adapted from SCALE-SDM:
     * Copyright (c) 2012-2014, Team SCALE, All rights reserved. */
    /*! \brief Compute vz_b(x) */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT vz_b_x ( const RT a_x /*!< x */) const noexcept
    {
        RT x1 = a_x;
        RT x2 = x1 * x1;
        RT x3 = x1 * x2;
        RT y = - RT(0.318657E+1)
               + RT(0.9926960)    * x1
               - RT(0.153193E-2)  * x2
               - RT(0.987059E-3)  * x3
               - RT(0.578878e-3)  * x2*x2
               + RT(0.855176E-4)  * x3*x2
               - RT(0.327815E-5)  * x3*x3;
        return y;
    }

    /*! \brief Compute vc_b(x) */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT vz_c_x ( const RT a_x /*!< x */) const noexcept
    {
        RT x1 = a_x;
        RT x2 = x1 * x1;
        RT x3 = x1 * x2;
        RT y = - RT(0.500015E+1)
               + RT(0.523778E+1)  * x1
               - RT(0.204914E+1)  * x2
               + RT(0.47529400)   * x3
               - RT(0.542819E-1)  * x2*x2
               + RT(0.238449E-2)  * x3*x2;
        return y;
    }

};

#endif
