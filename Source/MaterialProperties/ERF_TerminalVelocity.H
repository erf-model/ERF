#ifndef TERMINALVELOCITY_H
#define TERMINALVELOCITY_H

#include <cmath>
#include <AMReX_Gpu.H>
#include "ERF_Constants.H"

template <typename RT /*!< real-type */>
struct TerminalVelocity
{
    RT m_rho;

    /*! \brief Terminal velocity - Rogers & Yau, 1989 */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT RogersYau( const RT a_r) const noexcept
    {
        static const RT k1 = 1.233e8; // m^{-1} s^{-1}
        return (k1 * a_r * a_r);
    }

    /*! \brief Terminal velocity - Atlas & Ulbrich, 1977 */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT AtlasUlbrich( const RT a_r) const noexcept
    {
        static const RT a = 3.778;
        static const RT b = 0.67;
        RT D = 2*a_r*1000; // meter -> mm
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
        RT lb4l = 6.62E-6;    // base length[cm] for mean free path  of air molecules
        RT visb4l = 1.818E-4; // base dynamic viscosity[g/(cm*s)] for mean free path of air molecules
        RT pb4l = 1013.25;    // base pressure[hPa] for mean free path of air molecules
        RT tb4l = 293.15;     // base temperature[20C] for mean free path of air molecules

        RT P_hPa = a_P / 100.0;           // [Pa] -> [hPa]
        RT T_degC = a_T - 273.15;         // [K] -> [deg Celsius]
        RT diameter_cm = 2.0*a_r*100.0;   // [m] -> [cm]
        RT rho_mat_cgs = m_rho/1000.0;    // [kg/m^3] -> [g/cm^3]
        RT rho_air_cgs = a_rho/1000.0;    // [kg/m^3] -> [g/cm^3]
        RT grav_cgs = CONST_GRAV * 100;   // [m/s^2] -> [cm/s^2]

        RT gxdrow = grav_cgs * ( rho_mat_cgs - rho_air_cgs );

        RT viscosity = 0.0;
        if (T_degC >= 0) {
            viscosity = (1.7180 + 4.9e-3*T_degC) * 1.0e-4;
        } else {
            viscosity = (  1.7180 + 4.9e-3*T_degC
                         - 1.2e-5*T_degC*T_degC     ) * 1.0e-4;
        }

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

    /* The following code is adapted from SCALE-SDM:
     * Copyright (c) 2012-2014, Team SCALE, All rights reserved. */
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

};

#endif
