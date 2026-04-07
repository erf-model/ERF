#ifndef SUPERDROPLET_PC_COALESCENCE_H_
#define SUPERDROPLET_PC_COALESCENCE_H_

#include <algorithm>
#include "ERF_Constants.H"

#ifdef ERF_USE_PARTICLES

/*! \brief Coalescence kernels

    Reference for Long's and Hall's kernels:
    + Shima et al., 2009, "The super-droplet method for the numerical simulation of clouds and
      precipitation: A particle-based and probabilistic microphysics model coupled with a non-
      hydrostatic model." Quart. J. Roy. Meteorol. Soc., 135: 1307-1320
    + https://github.com/Shima-Lab/SCALE-SDM_BOMEX_Sato2018/blob/master/contrib/SDM/sdm_motion.f90
*/
template<typename RT, int SPACEDIM>
struct CollisionKernel
{
    /*! \brief Golovin kernel */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT golovin ( const RT a_r1, /*!< radius of droplet 1 */
                 const RT a_r2  /*!< radius of droplet 2 */ ) const noexcept
    {
        RT b = amrex::Real(1.5e03);
        auto X_1 = (amrex::Real(4.0)*PI/three) * a_r1*a_r1*a_r1;
        auto X_2 = (amrex::Real(4.0)*PI/three) * a_r2*a_r2*a_r2;
        return b * (X_1 + X_2);
    }

    /* The following code is adapted from SCALE-SDM:
     * Copyright (c) 2012-2014, Team SCALE, All rights reserved. */
    /*! \brief Calculate collision kernel for gravitational sedimentation
     *
     * This function computes the collision kernel for collisions driven by
     * differential gravitational settling. It uses the terminal velocities
     * of particles to determine collision probability.
     *
     * \param[in] a_r1 Radius of particle 1 (m)
     * \param[in] a_r2 Radius of particle 2 (m)
     * \param[in] a_v1 Pointer to velocity vector of particle 1 (m/s)
     * \param[in] a_v2 Pointer to velocity vector of particle 2 (m/s)
     * \return Collision kernel value (m³/s)
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT sedimentation (  const RT a_r1,
                        const RT a_r2,
                        const RT* const a_v1,
                        const RT* const a_v2 ) const noexcept
    {
        auto p = std::min(a_r1,a_r2) / std::max(a_r1,a_r2);
        auto E = myhalf*p*p / ((one+p)*(one+p));
        RT dv = 0;
        for (int d = 0; d < SPACEDIM; d++) {
            dv += (a_v1[d]-a_v2[d])*(a_v1[d]-a_v2[d]);
        }
        dv = std::sqrt(dv);
        return PI*(a_r1+a_r2)*(a_r1+a_r2)*E*dv;
    }

    /* The following code is adapted from SCALE-SDM:
     * Copyright (c) 2012-2014, Team SCALE, All rights reserved. */
    /*! \brief Calculate collision kernel using Long's formulation
     *
     * Long's kernel is an empirical formulation for cloud droplet collision-coalescence
     * that accounts for different collision rates based on droplet size. It is particularly
     * effective for smaller cloud droplets where hydrodynamic effects are important.
     *
     * Reference: Long, A. B., 1974: Solutions to the droplet collection equation
     * for polynomial kernels. J. Atmos. Sci., 31, 1040-amrex::Real(1052.)
     *
     * \param[in] a_r1 Radius of particle 1 (m)
     * \param[in] a_r2 Radius of particle 2 (m)
     * \param[in] a_v1 Pointer to velocity vector of particle 1 (m/s)
     * \param[in] a_v2 Pointer to velocity vector of particle 2 (m/s)
     * \return Collision kernel value (m³/s)
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT Longs ( const RT a_r1,
               const RT a_r2,
               const RT* const a_v1,
               const RT* const a_v2 ) const noexcept
    {
        auto r1 = std::max( a_r1, a_r2 );  // large
        auto r2 = std::min( a_r1, a_r2 );  // small
        auto sumr = r1 + r2;

        RT c_rate = one;
        if( r1 <= amrex::Real(5.E-5) ) {
            c_rate = amrex::Real(4.5E+8) * ( r1*r1 ) * ( one - amrex::Real(3.E-6)/(std::max(amrex::Real(3.01E-6),r1)) );
        }
        c_rate *= (PI*sumr*sumr);

        RT dv = 0;
        for (int d = 0; d < SPACEDIM; d++) {
            dv += (a_v1[d]-a_v2[d])*(a_v1[d]-a_v2[d]);
        }
        dv = std::sqrt(dv);

        return c_rate*dv;
    }

    /* The following code is adapted from SCALE-SDM:
     * Copyright (c) 2012-2014, Team SCALE, All rights reserved. */
    /*! \brief r0col data for Hall's kernel */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT Halls_data_r0col (const int a_i /*!< index */) const noexcept
    {
        static const RT vec[15] = {  amrex::Real(6.0),
                                     amrex::Real(8.0),
                                    amrex::Real(10.0),
                                    amrex::Real(15.0),
                                    amrex::Real(20.0),
                                    amrex::Real(25.0),
                                    amrex::Real(30.0),
                                    amrex::Real(40.0),
                                    amrex::Real(50.0),
                                    amrex::Real(60.0),
                                    amrex::Real(70.0),
                                   amrex::Real(100.0),
                                   amrex::Real(150.0),
                                   amrex::Real(200.0),
                                   amrex::Real(300.0)};
        return vec[a_i];
    }

    /* The following code is adapted from SCALE-SDM:
     * Copyright (c) 2012-2014, Team SCALE, All rights reserved. */
    /*! \brief r0col data for Hall's kernel */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int Halls_data_r0col ( const RT a_r /*!< radius */) const noexcept
    {
        auto radius_microns = a_r*amrex::Real(1.0e6); // [m] -> [mu-m]
        int i;
        for (i = 0; i < 15; i++) {
            if (radius_microns <= Halls_data_r0col(i)) { return i; }
        }
        return i;
    }

    /* The following code is adapted from SCALE-SDM:
     * Copyright (c) 2012-2014, Team SCALE, All rights reserved. */
    /*! \brief ratcol data for Hall's kernel */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT Halls_data_ratcol ( const int a_i /*!< index */) const noexcept
    {
        static const RT vec[21] = {amrex::Real(0.00),
                                   amrex::Real(0.05),
                                   amrex::Real(0.10),
                                   amrex::Real(0.15),
                                   amrex::Real(0.20),
                                   fourth,
                                   amrex::Real(0.30),
                                   amrex::Real(0.35),
                                   amrex::Real(0.40),
                                   amrex::Real(0.45),
                                   amrex::Real(0.50),
                                   amrex::Real(0.55),
                                   amrex::Real(0.60),
                                   amrex::Real(0.65),
                                   amrex::Real(0.70),
                                   amrex::Real(0.75),
                                   amrex::Real(0.80),
                                   amrex::Real(0.85),
                                   amrex::Real(0.90),
                                   amrex::Real(0.95),
                                   amrex::Real(1.00)};
        return vec[a_i];
    }

    /* The following code is adapted from SCALE-SDM:
     * Copyright (c) 2012-2014, Team SCALE, All rights reserved. */
    /*! \brief ratcol data for Hall's kernel */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int Halls_data_ratcol ( const RT a_r /*!< ratio */) const noexcept
    {
        int i;
        for (i = 1; i < 20; i++) {
            if (a_r <= Halls_data_ratcol(i)) { return i; }
        }
        return i;
    }

    /* The following code is adapted from SCALE-SDM:
     * Copyright (c) 2012-2014, Team SCALE, All rights reserved. */
    /*! \brief ecoll data for Hall's kernel */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT Halls_data_ecoll ( const int a_i, /*!< row index */
                          const int a_j /*!< column index */) const noexcept
    {
        const RT ecoll[21][15] =
            { {amrex::Real(0.0010),amrex::Real(0.0010),amrex::Real(0.0010),amrex::Real(0.0010),amrex::Real(0.0010),amrex::Real(0.0010),amrex::Real(0.0010),amrex::Real(0.0010),amrex::Real(0.0010),amrex::Real(0.0010),amrex::Real(0.0010),amrex::Real(0.0010),amrex::Real(0.0010),amrex::Real(0.0010),amrex::Real(0.0010)},
              {amrex::Real(0.0030),amrex::Real(0.0030),amrex::Real(0.0030),amrex::Real(0.0040),amrex::Real(0.0050),amrex::Real(0.0050),amrex::Real(0.0050),amrex::Real(0.0100),amrex::Real(0.1000),amrex::Real(0.0500),amrex::Real(0.2000),amrex::Real(0.5000),amrex::Real(0.7700),amrex::Real(0.8700),amrex::Real(0.9700)},
              {amrex::Real(0.0070),amrex::Real(0.0070),amrex::Real(0.0070),amrex::Real(0.0080),amrex::Real(0.0090),amrex::Real(0.0100),amrex::Real(0.0100),amrex::Real(0.0700),amrex::Real(0.4000),amrex::Real(0.4300),amrex::Real(0.5800),amrex::Real(0.7900),amrex::Real(0.9300),amrex::Real(0.9600),amrex::Real(1.0000)},
              {amrex::Real(0.0090),amrex::Real(0.0090),amrex::Real(0.0090),amrex::Real(0.0120),amrex::Real(0.0150),amrex::Real(0.0100),amrex::Real(0.0200),amrex::Real(0.2800),amrex::Real(0.6000),amrex::Real(0.6400),amrex::Real(0.7500),amrex::Real(0.9100),amrex::Real(0.9700),amrex::Real(0.9800),amrex::Real(1.0000)},
              {amrex::Real(0.0140),amrex::Real(0.0140),amrex::Real(0.0140),amrex::Real(0.0150),amrex::Real(0.0160),amrex::Real(0.0300),amrex::Real(0.0600),amrex::Real(0.5000),amrex::Real(0.7000),amrex::Real(0.7700),amrex::Real(0.8400),amrex::Real(0.9500),amrex::Real(0.9700),amrex::Real(1.0000),amrex::Real(1.0000)},
              {amrex::Real(0.0170),amrex::Real(0.0170),amrex::Real(0.0170),amrex::Real(0.0200),amrex::Real(0.0220),amrex::Real(0.0600),amrex::Real(0.1000),amrex::Real(0.6200),amrex::Real(0.7800),amrex::Real(0.8400),amrex::Real(0.8800),amrex::Real(0.9500),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000)},
              {amrex::Real(0.0300),amrex::Real(0.0300),amrex::Real(0.0240),amrex::Real(0.0220),amrex::Real(0.0320),amrex::Real(0.0620),amrex::Real(0.2000),amrex::Real(0.6800),amrex::Real(0.8300),amrex::Real(0.8700),amrex::Real(0.9000),amrex::Real(0.9500),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000)},
              {amrex::Real(0.0250),amrex::Real(0.0250),amrex::Real(0.0250),amrex::Real(0.0360),amrex::Real(0.0430),amrex::Real(0.1300),amrex::Real(0.2700),amrex::Real(0.7400),amrex::Real(0.8600),amrex::Real(0.8900),amrex::Real(0.9200),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000)},
              {amrex::Real(0.0270),amrex::Real(0.0270),amrex::Real(0.0270),amrex::Real(0.0400),amrex::Real(0.0520),amrex::Real(0.2000),amrex::Real(0.4000),amrex::Real(0.7800),amrex::Real(0.8800),amrex::Real(0.9000),amrex::Real(0.9400),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000)},
              {amrex::Real(0.0300),amrex::Real(0.0300),amrex::Real(0.0300),amrex::Real(0.0470),amrex::Real(0.0640),amrex::Real(0.2500),amrex::Real(0.5000),amrex::Real(0.8000),amrex::Real(0.9000),amrex::Real(0.9100),amrex::Real(0.9500),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000)},
              {amrex::Real(0.0400),amrex::Real(0.0400),amrex::Real(0.0330),amrex::Real(0.0370),amrex::Real(0.0680),amrex::Real(0.2400),amrex::Real(0.5500),amrex::Real(0.8000),amrex::Real(0.9000),amrex::Real(0.9100),amrex::Real(0.9500),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000)},
              {amrex::Real(0.0350),amrex::Real(0.0350),amrex::Real(0.0350),amrex::Real(0.0550),amrex::Real(0.0790),amrex::Real(0.2900),amrex::Real(0.5800),amrex::Real(0.8000),amrex::Real(0.9000),amrex::Real(0.9100),amrex::Real(0.9500),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000)},
              {amrex::Real(0.0370),amrex::Real(0.0370),amrex::Real(0.0370),amrex::Real(0.0620),amrex::Real(0.0820),amrex::Real(0.2900),amrex::Real(0.5900),amrex::Real(0.7800),amrex::Real(0.9000),amrex::Real(0.9100),amrex::Real(0.9500),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000)},
              {amrex::Real(0.0370),amrex::Real(0.0370),amrex::Real(0.0370),amrex::Real(0.0600),amrex::Real(0.0800),amrex::Real(0.2900),amrex::Real(0.5800),amrex::Real(0.7700),amrex::Real(0.8900),amrex::Real(0.9100),amrex::Real(0.9500),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000)},
              {amrex::Real(0.0370),amrex::Real(0.0370),amrex::Real(0.0370),amrex::Real(0.0410),amrex::Real(0.0750),amrex::Real(0.2500),amrex::Real(0.5400),amrex::Real(0.7600),amrex::Real(0.8800),amrex::Real(0.9200),amrex::Real(0.9500),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000)},
              {amrex::Real(0.0370),amrex::Real(0.0370),amrex::Real(0.0370),amrex::Real(0.0520),amrex::Real(0.0670),amrex::Real(0.2500),amrex::Real(0.5100),amrex::Real(0.7700),amrex::Real(0.8800),amrex::Real(0.9300),amrex::Real(0.9700),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000)},
              {amrex::Real(0.0370),amrex::Real(0.0370),amrex::Real(0.0370),amrex::Real(0.0470),amrex::Real(0.0570),amrex::Real(0.2500),amrex::Real(0.4900),amrex::Real(0.7700),amrex::Real(0.8900),amrex::Real(0.9500),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000),amrex::Real(1.0000)},
              {amrex::Real(0.0360),amrex::Real(0.0360),amrex::Real(0.0360),amrex::Real(0.0420),amrex::Real(0.0480),amrex::Real(0.2300),amrex::Real(0.4700),amrex::Real(0.7800),amrex::Real(0.9200),amrex::Real(1.0000),amrex::Real(1.0200),amrex::Real(1.0200),amrex::Real(1.0200),amrex::Real(1.0200),amrex::Real(1.0200)},
              {amrex::Real(0.0400),amrex::Real(0.0400),amrex::Real(0.0350),amrex::Real(0.0330),amrex::Real(0.0400),amrex::Real(0.1120),amrex::Real(0.4500),amrex::Real(0.7900),amrex::Real(1.0100),amrex::Real(1.0300),amrex::Real(1.0400),amrex::Real(1.0400),amrex::Real(1.0400),amrex::Real(1.0400),amrex::Real(1.0400)},
              {amrex::Real(0.0330),amrex::Real(0.0330),amrex::Real(0.0330),amrex::Real(0.0330),amrex::Real(0.0330),amrex::Real(0.1190),amrex::Real(0.4700),amrex::Real(0.9500),amrex::Real(1.3000),amrex::Real(1.7000),amrex::Real(2.3000),amrex::Real(2.3000),amrex::Real(2.3000),amrex::Real(2.3000),amrex::Real(2.3000)},
              {amrex::Real(0.0270),amrex::Real(0.0270),amrex::Real(0.0270),amrex::Real(0.0270),amrex::Real(0.0270),amrex::Real(0.1250),amrex::Real(0.5200),amrex::Real(1.4000),amrex::Real(2.3000),amrex::Real(3.0000),amrex::Real(4.0000),amrex::Real(4.0000),amrex::Real(4.0000),amrex::Real(4.0000),amrex::Real(4.0000)} };
        return ecoll[a_i][a_j];
    }


    /* The following code is adapted from SCALE-SDM:
     * Copyright (c) 2012-2014, Team SCALE, All rights reserved. */
    /*! \brief Hall's kernel */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT Halls ( const RT a_r1, /*!< radius of particle 1 */
               const RT a_r2, /*!< radius of particle 2 */
               const RT* const a_v1, /*!< velocity of particle 1 */
               const RT* const a_v2  /*!< velocity of particle 2 */) const noexcept
    {
        auto r1 = std::max( a_r1, a_r2 );  // large
        auto r2 = std::min( a_r1, a_r2 );  // small
        auto ratio = r2 / r1;
        auto sumr = r1 + r2;

        int irr = Halls_data_r0col(r1);
        int iqq = Halls_data_ratcol(ratio);

        RT c_rate = one;

        if( irr >= 15 ) {

            RT q  =   ( ratio                  - Halls_data_ratcol(iqq-1) )
                    / ( Halls_data_ratcol(iqq) - Halls_data_ratcol(iqq-1) );

            RT ek =   (one-q) * Halls_data_ecoll(iqq-1,14)
                    + q       * Halls_data_ecoll(iqq  ,14);

            c_rate = std::min( ek, amrex::Real(1.0) );

        } else if ((irr >= 1) && (irr < 15)) {

            RT p =   ( r1*amrex::Real(1.0e6)              - Halls_data_r0col(irr-1) )
                   / ( Halls_data_r0col(irr) - Halls_data_r0col(irr-1) );

            RT q =   ( ratio                  - Halls_data_ratcol(iqq-1) )
                   / ( Halls_data_ratcol(iqq) - Halls_data_ratcol(iqq-1) );

            c_rate =   (one-p) * (one-q) * Halls_data_ecoll(iqq-1,irr-1)
                     + p       * (one-q) * Halls_data_ecoll(iqq-1,irr  )
                     + (one-p) * q       * Halls_data_ecoll(iqq  ,irr-1)
                     + p       * q       * Halls_data_ecoll(iqq  ,irr  );

        } else {

            RT q =  ( ratio                  - Halls_data_ratcol(iqq-1) )
                  / ( Halls_data_ratcol(iqq) - Halls_data_ratcol(iqq-1) );

            c_rate =   (one-q) * Halls_data_ecoll(iqq-1,0)
                     + q       * Halls_data_ecoll(iqq  ,0);

        }
        c_rate *= (PI*sumr*sumr);

        RT dv = 0;
        for (int d = 0; d < SPACEDIM; d++) {
            dv += (a_v1[d]-a_v2[d])*(a_v1[d]-a_v2[d]);
        }
        dv = std::sqrt(dv);

        return c_rate*dv;
    }

    /* The following code is adapted from SCALE-SDM:
     * Copyright (c) 2012-2014, Team SCALE, All rights reserved. */
    /*! \brief Brownian kernel */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT Brownian_SeinfeldPandis ( const RT a_r1, /*!< radius of particle 1 */
                                 const RT a_r2, /*!< radius of particle 2 */
                                 const RT a_mass_1, /*!< total mass of particle 1 */
                                 const RT a_mass_2, /*!< total mass of particle 2 */
                                 const RT a_pressure, /*!< pressure */
                                 const RT a_temperature /*!< temperature */ ) const noexcept
    {
        // diameter of droplets
        RT diameter_1 = 2*a_r1;
        RT diameter_2 = 2*a_r2;

        RT p_crs = a_pressure; // [Pa]
        RT t_crs = a_temperature; // [K]

        // dynamic viscosity [Pa*s] (Pruppacher & Klett,1997)
        RT T_degC = t_crs - amrex::Real(273.15); // [K] => [degC]
        RT vis_crs = zero;
        if( T_degC >= zero ) {
            vis_crs = ( amrex::Real(1.7180) + amrex::Real(4.9E-3)*T_degC ) * amrex::Real(1.E-5);
        } else {
            vis_crs = ( amrex::Real(1.7180) + amrex::Real(4.9E-3)*T_degC -amrex::Real(1.2E-5)*T_degC*T_degC ) * amrex::Real(1.E-5);
        }

        //air mean free path [m]
        RT lambda_crs = (two*vis_crs) / (p_crs*std::sqrt(amrex::Real(8.0)*mwdair*amrex::Real(1.E-3)/(PI*R_d*t_crs)));

        // temporary var
        RT dtmp = zero;

        // slip correction of droplets [-]
        dtmp   = amrex::Real(1.2570) + amrex::Real(0.40) * exp(-amrex::Real(0.550)*diameter_1/lambda_crs);
        RT slip_corr_1 = one + (two*lambda_crs*dtmp)/(diameter_1);
        dtmp   = amrex::Real(1.2570) + amrex::Real(0.40) * exp(-amrex::Real(0.550)*diameter_2/lambda_crs);
        RT slip_corr_2 = one + (two*lambda_crs*dtmp)/(diameter_2);

        // diffusion term [m*m/s]
        dtmp = (boltz*t_crs)/(three*PI*vis_crs);
        RT diff_term_1 = dtmp * (slip_corr_1/diameter_1);
        RT diff_term_2 = dtmp * (slip_corr_2/diameter_2);

        // velocity term [m/s]
        dtmp = (amrex::Real(8.0)*boltz*t_crs)/PI;
        RT vel_term_1 = std::sqrt(dtmp/a_mass_1);
        RT vel_term_2 = std::sqrt(dtmp/a_mass_2);

        // mean free path of droplets [m]
        dtmp = amrex::Real(8.0)/PI;
        RT lambda_1 = dtmp * (diff_term_1/vel_term_1);
        RT lambda_2 = dtmp * (diff_term_2/vel_term_2);

        //length term [m]
        dtmp = (diameter_1+lambda_1)*(diameter_1+lambda_1)*(diameter_1+lambda_1)
               - std::exp(amrex::Real(1.5)*std::log(diameter_1*diameter_1+lambda_1*lambda_1));
        RT length_term_1 = dtmp/(three*diameter_1*lambda_1) - diameter_1;
        dtmp = (diameter_2+lambda_2)*(diameter_2+lambda_2)*(diameter_2+lambda_2)
               - std::exp(amrex::Real(1.5)*std::log(diameter_2*diameter_2+lambda_2*lambda_2));
        RT length_term_2 = dtmp/(three*diameter_2*lambda_2) - diameter_2;

        // Brownian Coagulation Coefficient K12 [m3/s]
        RT sumdia = diameter_1 + diameter_2;
        RT sumd   = diff_term_1   + diff_term_2;
        RT sumc   = std::sqrt( vel_term_1*vel_term_1 + vel_term_2*vel_term_2 );
        RT sumg   = std::sqrt( two*length_term_1*length_term_1 + two*length_term_2*length_term_2 );

        dtmp = sumdia/(sumdia+two*sumg) + (amrex::Real(8.0)*sumd)/(sumdia*sumc);
        RT k12 = two*PI * sumdia*sumd/dtmp;

        return k12;
    }

};

#endif
#endif

