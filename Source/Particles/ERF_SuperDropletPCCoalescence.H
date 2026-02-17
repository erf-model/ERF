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
        RT b = 1.5e03;
        auto X_1 = (4.0*PI/3.0) * a_r1*a_r1*a_r1;
        auto X_2 = (4.0*PI/3.0) * a_r2*a_r2*a_r2;
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
        auto E = 0.5*p*p / ((1.0+p)*(1.0+p));
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
     * for polynomial kernels. J. Atmos. Sci., 31, 1040-1052.
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

        RT c_rate = 1.0;
        if( r1 <= 5.E-5 ) {
            c_rate = 4.5E+8 * ( r1*r1 ) * ( 1.0 - 3.E-6/(std::max(3.01E-6,r1)) );
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
        static const RT vec[15] = {  6.0,
                                     8.0,
                                    10.0,
                                    15.0,
                                    20.0,
                                    25.0,
                                    30.0,
                                    40.0,
                                    50.0,
                                    60.0,
                                    70.0,
                                   100.0,
                                   150.0,
                                   200.0,
                                   300.0};
        return vec[a_i];
    }

    /* The following code is adapted from SCALE-SDM:
     * Copyright (c) 2012-2014, Team SCALE, All rights reserved. */
    /*! \brief r0col data for Hall's kernel */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int Halls_data_r0col ( const RT a_r /*!< radius */) const noexcept
    {
        auto radius_microns = a_r*1.0e6; // [m] -> [mu-m]
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
        static const RT vec[21] = {0.00,
                                   0.05,
                                   0.10,
                                   0.15,
                                   0.20,
                                   0.25,
                                   0.30,
                                   0.35,
                                   0.40,
                                   0.45,
                                   0.50,
                                   0.55,
                                   0.60,
                                   0.65,
                                   0.70,
                                   0.75,
                                   0.80,
                                   0.85,
                                   0.90,
                                   0.95,
                                   1.00};
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
            { {0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010},
              {0.0030,0.0030,0.0030,0.0040,0.0050,0.0050,0.0050,0.0100,0.1000,0.0500,0.2000,0.5000,0.7700,0.8700,0.9700},
              {0.0070,0.0070,0.0070,0.0080,0.0090,0.0100,0.0100,0.0700,0.4000,0.4300,0.5800,0.7900,0.9300,0.9600,1.0000},
              {0.0090,0.0090,0.0090,0.0120,0.0150,0.0100,0.0200,0.2800,0.6000,0.6400,0.7500,0.9100,0.9700,0.9800,1.0000},
              {0.0140,0.0140,0.0140,0.0150,0.0160,0.0300,0.0600,0.5000,0.7000,0.7700,0.8400,0.9500,0.9700,1.0000,1.0000},
              {0.0170,0.0170,0.0170,0.0200,0.0220,0.0600,0.1000,0.6200,0.7800,0.8400,0.8800,0.9500,1.0000,1.0000,1.0000},
              {0.0300,0.0300,0.0240,0.0220,0.0320,0.0620,0.2000,0.6800,0.8300,0.8700,0.9000,0.9500,1.0000,1.0000,1.0000},
              {0.0250,0.0250,0.0250,0.0360,0.0430,0.1300,0.2700,0.7400,0.8600,0.8900,0.9200,1.0000,1.0000,1.0000,1.0000},
              {0.0270,0.0270,0.0270,0.0400,0.0520,0.2000,0.4000,0.7800,0.8800,0.9000,0.9400,1.0000,1.0000,1.0000,1.0000},
              {0.0300,0.0300,0.0300,0.0470,0.0640,0.2500,0.5000,0.8000,0.9000,0.9100,0.9500,1.0000,1.0000,1.0000,1.0000},
              {0.0400,0.0400,0.0330,0.0370,0.0680,0.2400,0.5500,0.8000,0.9000,0.9100,0.9500,1.0000,1.0000,1.0000,1.0000},
              {0.0350,0.0350,0.0350,0.0550,0.0790,0.2900,0.5800,0.8000,0.9000,0.9100,0.9500,1.0000,1.0000,1.0000,1.0000},
              {0.0370,0.0370,0.0370,0.0620,0.0820,0.2900,0.5900,0.7800,0.9000,0.9100,0.9500,1.0000,1.0000,1.0000,1.0000},
              {0.0370,0.0370,0.0370,0.0600,0.0800,0.2900,0.5800,0.7700,0.8900,0.9100,0.9500,1.0000,1.0000,1.0000,1.0000},
              {0.0370,0.0370,0.0370,0.0410,0.0750,0.2500,0.5400,0.7600,0.8800,0.9200,0.9500,1.0000,1.0000,1.0000,1.0000},
              {0.0370,0.0370,0.0370,0.0520,0.0670,0.2500,0.5100,0.7700,0.8800,0.9300,0.9700,1.0000,1.0000,1.0000,1.0000},
              {0.0370,0.0370,0.0370,0.0470,0.0570,0.2500,0.4900,0.7700,0.8900,0.9500,1.0000,1.0000,1.0000,1.0000,1.0000},
              {0.0360,0.0360,0.0360,0.0420,0.0480,0.2300,0.4700,0.7800,0.9200,1.0000,1.0200,1.0200,1.0200,1.0200,1.0200},
              {0.0400,0.0400,0.0350,0.0330,0.0400,0.1120,0.4500,0.7900,1.0100,1.0300,1.0400,1.0400,1.0400,1.0400,1.0400},
              {0.0330,0.0330,0.0330,0.0330,0.0330,0.1190,0.4700,0.9500,1.3000,1.7000,2.3000,2.3000,2.3000,2.3000,2.3000},
              {0.0270,0.0270,0.0270,0.0270,0.0270,0.1250,0.5200,1.4000,2.3000,3.0000,4.0000,4.0000,4.0000,4.0000,4.0000} };
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

        RT c_rate = 1.0;

        if( irr >= 15 ) {

            RT q  =   ( ratio                  - Halls_data_ratcol(iqq-1) )
                    / ( Halls_data_ratcol(iqq) - Halls_data_ratcol(iqq-1) );

            RT ek =   (1.0-q) * Halls_data_ecoll(iqq-1,14)
                    + q       * Halls_data_ecoll(iqq  ,14);

            c_rate = std::min( ek, 1.0 );

        } else if ((irr >= 1) && (irr < 15)) {

            RT p =   ( r1*1.0e6              - Halls_data_r0col(irr-1) )
                   / ( Halls_data_r0col(irr) - Halls_data_r0col(irr-1) );

            RT q =   ( ratio                  - Halls_data_ratcol(iqq-1) )
                   / ( Halls_data_ratcol(iqq) - Halls_data_ratcol(iqq-1) );

            c_rate =   (1.0-p) * (1.0-q) * Halls_data_ecoll(iqq-1,irr-1)
                     + p       * (1.0-q) * Halls_data_ecoll(iqq-1,irr  )
                     + (1.0-p) * q       * Halls_data_ecoll(iqq  ,irr-1)
                     + p       * q       * Halls_data_ecoll(iqq  ,irr  );

        } else {

            RT q =  ( ratio                  - Halls_data_ratcol(iqq-1) )
                  / ( Halls_data_ratcol(iqq) - Halls_data_ratcol(iqq-1) );

            c_rate =   (1.0-q) * Halls_data_ecoll(iqq-1,0)
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
        RT T_degC = t_crs - 273.15; // [K] => [degC]
        RT vis_crs = 0.0;
        if( T_degC >= 0.0 ) {
            vis_crs = ( 1.7180 + 4.9E-3*T_degC ) * 1.E-5;
        } else {
            vis_crs = ( 1.7180 + 4.9E-3*T_degC -1.2E-5*T_degC*T_degC ) * 1.E-5;
        }

        //air mean free path [m]
        RT lambda_crs = (2.0*vis_crs) / (p_crs*std::sqrt(8.0*mwdair*1.E-3/(PI*R_d*t_crs)));

        // temporary var
        RT dtmp = 0.0;

        // slip correction of droplets [-]
        dtmp   = 1.2570 + 0.40 * exp(-0.550*diameter_1/lambda_crs);
        RT slip_corr_1 = 1.0 + (2.0*lambda_crs*dtmp)/(diameter_1);
        dtmp   = 1.2570 + 0.40 * exp(-0.550*diameter_2/lambda_crs);
        RT slip_corr_2 = 1.0 + (2.0*lambda_crs*dtmp)/(diameter_2);

        // diffusion term [m*m/s]
        dtmp = (boltz*t_crs)/(3.0*PI*vis_crs);
        RT diff_term_1 = dtmp * (slip_corr_1/diameter_1);
        RT diff_term_2 = dtmp * (slip_corr_2/diameter_2);

        // velocity term [m/s]
        dtmp = (8.0*boltz*t_crs)/PI;
        RT vel_term_1 = std::sqrt(dtmp/a_mass_1);
        RT vel_term_2 = std::sqrt(dtmp/a_mass_2);

        // mean free path of droplets [m]
        dtmp = 8.0/PI;
        RT lambda_1 = dtmp * (diff_term_1/vel_term_1);
        RT lambda_2 = dtmp * (diff_term_2/vel_term_2);

        //length term [m]
        dtmp = (diameter_1+lambda_1)*(diameter_1+lambda_1)*(diameter_1+lambda_1)
               - std::exp(1.5*std::log(diameter_1*diameter_1+lambda_1*lambda_1));
        RT length_term_1 = dtmp/(3.0*diameter_1*lambda_1) - diameter_1;
        dtmp = (diameter_2+lambda_2)*(diameter_2+lambda_2)*(diameter_2+lambda_2)
               - std::exp(1.5*std::log(diameter_2*diameter_2+lambda_2*lambda_2));
        RT length_term_2 = dtmp/(3.0*diameter_2*lambda_2) - diameter_2;

        // Brownian Coagulation Coefficient K12 [m3/s]
        RT sumdia = diameter_1 + diameter_2;
        RT sumd   = diff_term_1   + diff_term_2;
        RT sumc   = std::sqrt( vel_term_1*vel_term_1 + vel_term_2*vel_term_2 );
        RT sumg   = std::sqrt( 2.0*length_term_1*length_term_1 + 2.0*length_term_2*length_term_2 );

        dtmp = sumdia/(sumdia+2.0*sumg) + (8.0*sumd)/(sumdia*sumc);
        RT k12 = 2.0*PI * sumdia*sumd/dtmp;

        return k12;
    }

};

#endif
#endif

