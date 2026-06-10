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
        RT b = RT(1.5e03);
        auto X_1 = RT(four_thirds_pi) * a_r1*a_r1*a_r1;
        auto X_2 = RT(four_thirds_pi) * a_r2*a_r2*a_r2;
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
        auto E = RT(0.5)*p*p / ((RT(1.0)+p)*(RT(1.0)+p));
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
     * for polynomial kernels. J. Atmos. Sci., 31, 1040-RT(1052.)
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

        RT c_rate = RT(1.0);
        if( r1 <= RT(5.E-5) ) {
            c_rate = RT(4.5E+8) * ( r1*r1 ) * ( RT(1.0) - RT(3.E-6)/(std::max(RT(3.01E-6),r1)) );
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
        static const RT vec[15] = {  RT(6.0),
                                     RT(8.0),
                                    RT(10.0),
                                    RT(15.0),
                                    RT(20.0),
                                    RT(25.0),
                                    RT(30.0),
                                    RT(40.0),
                                    RT(50.0),
                                    RT(60.0),
                                    RT(70.0),
                                   RT(100.0),
                                   RT(150.0),
                                   RT(200.0),
                                   RT(300.0)};
        return vec[a_i];
    }

    /* The following code is adapted from SCALE-SDM:
     * Copyright (c) 2012-2014, Team SCALE, All rights reserved. */
    /*! \brief r0col data for Hall's kernel */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int Halls_data_r0col ( const RT a_r /*!< radius */) const noexcept
    {
        auto radius_microns = a_r*RT(1.0e6); // [m] -> [mu-m]
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
        static const RT vec[21] = {RT(0.00),
                                   RT(0.05),
                                   RT(0.10),
                                   RT(0.15),
                                   RT(0.20),
                                   fourth,
                                   RT(0.30),
                                   RT(0.35),
                                   RT(0.40),
                                   RT(0.45),
                                   RT(0.50),
                                   RT(0.55),
                                   RT(0.60),
                                   RT(0.65),
                                   RT(0.70),
                                   RT(0.75),
                                   RT(0.80),
                                   RT(0.85),
                                   RT(0.90),
                                   RT(0.95),
                                   RT(1.00)};
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
            { {RT(0.0010),RT(0.0010),RT(0.0010),RT(0.0010),RT(0.0010),RT(0.0010),RT(0.0010),RT(0.0010),RT(0.0010),RT(0.0010),RT(0.0010),RT(0.0010),RT(0.0010),RT(0.0010),RT(0.0010)},
              {RT(0.0030),RT(0.0030),RT(0.0030),RT(0.0040),RT(0.0050),RT(0.0050),RT(0.0050),RT(0.0100),RT(0.1000),RT(0.0500),RT(0.2000),RT(0.5000),RT(0.7700),RT(0.8700),RT(0.9700)},
              {RT(0.0070),RT(0.0070),RT(0.0070),RT(0.0080),RT(0.0090),RT(0.0100),RT(0.0100),RT(0.0700),RT(0.4000),RT(0.4300),RT(0.5800),RT(0.7900),RT(0.9300),RT(0.9600),RT(1.0000)},
              {RT(0.0090),RT(0.0090),RT(0.0090),RT(0.0120),RT(0.0150),RT(0.0100),RT(0.0200),RT(0.2800),RT(0.6000),RT(0.6400),RT(0.7500),RT(0.9100),RT(0.9700),RT(0.9800),RT(1.0000)},
              {RT(0.0140),RT(0.0140),RT(0.0140),RT(0.0150),RT(0.0160),RT(0.0300),RT(0.0600),RT(0.5000),RT(0.7000),RT(0.7700),RT(0.8400),RT(0.9500),RT(0.9700),RT(1.0000),RT(1.0000)},
              {RT(0.0170),RT(0.0170),RT(0.0170),RT(0.0200),RT(0.0220),RT(0.0600),RT(0.1000),RT(0.6200),RT(0.7800),RT(0.8400),RT(0.8800),RT(0.9500),RT(1.0000),RT(1.0000),RT(1.0000)},
              {RT(0.0300),RT(0.0300),RT(0.0240),RT(0.0220),RT(0.0320),RT(0.0620),RT(0.2000),RT(0.6800),RT(0.8300),RT(0.8700),RT(0.9000),RT(0.9500),RT(1.0000),RT(1.0000),RT(1.0000)},
              {RT(0.0250),RT(0.0250),RT(0.0250),RT(0.0360),RT(0.0430),RT(0.1300),RT(0.2700),RT(0.7400),RT(0.8600),RT(0.8900),RT(0.9200),RT(1.0000),RT(1.0000),RT(1.0000),RT(1.0000)},
              {RT(0.0270),RT(0.0270),RT(0.0270),RT(0.0400),RT(0.0520),RT(0.2000),RT(0.4000),RT(0.7800),RT(0.8800),RT(0.9000),RT(0.9400),RT(1.0000),RT(1.0000),RT(1.0000),RT(1.0000)},
              {RT(0.0300),RT(0.0300),RT(0.0300),RT(0.0470),RT(0.0640),RT(0.2500),RT(0.5000),RT(0.8000),RT(0.9000),RT(0.9100),RT(0.9500),RT(1.0000),RT(1.0000),RT(1.0000),RT(1.0000)},
              {RT(0.0400),RT(0.0400),RT(0.0330),RT(0.0370),RT(0.0680),RT(0.2400),RT(0.5500),RT(0.8000),RT(0.9000),RT(0.9100),RT(0.9500),RT(1.0000),RT(1.0000),RT(1.0000),RT(1.0000)},
              {RT(0.0350),RT(0.0350),RT(0.0350),RT(0.0550),RT(0.0790),RT(0.2900),RT(0.5800),RT(0.8000),RT(0.9000),RT(0.9100),RT(0.9500),RT(1.0000),RT(1.0000),RT(1.0000),RT(1.0000)},
              {RT(0.0370),RT(0.0370),RT(0.0370),RT(0.0620),RT(0.0820),RT(0.2900),RT(0.5900),RT(0.7800),RT(0.9000),RT(0.9100),RT(0.9500),RT(1.0000),RT(1.0000),RT(1.0000),RT(1.0000)},
              {RT(0.0370),RT(0.0370),RT(0.0370),RT(0.0600),RT(0.0800),RT(0.2900),RT(0.5800),RT(0.7700),RT(0.8900),RT(0.9100),RT(0.9500),RT(1.0000),RT(1.0000),RT(1.0000),RT(1.0000)},
              {RT(0.0370),RT(0.0370),RT(0.0370),RT(0.0410),RT(0.0750),RT(0.2500),RT(0.5400),RT(0.7600),RT(0.8800),RT(0.9200),RT(0.9500),RT(1.0000),RT(1.0000),RT(1.0000),RT(1.0000)},
              {RT(0.0370),RT(0.0370),RT(0.0370),RT(0.0520),RT(0.0670),RT(0.2500),RT(0.5100),RT(0.7700),RT(0.8800),RT(0.9300),RT(0.9700),RT(1.0000),RT(1.0000),RT(1.0000),RT(1.0000)},
              {RT(0.0370),RT(0.0370),RT(0.0370),RT(0.0470),RT(0.0570),RT(0.2500),RT(0.4900),RT(0.7700),RT(0.8900),RT(0.9500),RT(1.0000),RT(1.0000),RT(1.0000),RT(1.0000),RT(1.0000)},
              {RT(0.0360),RT(0.0360),RT(0.0360),RT(0.0420),RT(0.0480),RT(0.2300),RT(0.4700),RT(0.7800),RT(0.9200),RT(1.0000),RT(1.0200),RT(1.0200),RT(1.0200),RT(1.0200),RT(1.0200)},
              {RT(0.0400),RT(0.0400),RT(0.0350),RT(0.0330),RT(0.0400),RT(0.1120),RT(0.4500),RT(0.7900),RT(1.0100),RT(1.0300),RT(1.0400),RT(1.0400),RT(1.0400),RT(1.0400),RT(1.0400)},
              {RT(0.0330),RT(0.0330),RT(0.0330),RT(0.0330),RT(0.0330),RT(0.1190),RT(0.4700),RT(0.9500),RT(1.3000),RT(1.7000),RT(2.3000),RT(2.3000),RT(2.3000),RT(2.3000),RT(2.3000)},
              {RT(0.0270),RT(0.0270),RT(0.0270),RT(0.0270),RT(0.0270),RT(0.1250),RT(0.5200),RT(1.4000),RT(2.3000),RT(3.0000),RT(4.0000),RT(4.0000),RT(4.0000),RT(4.0000),RT(4.0000)} };
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

            c_rate = std::min( ek, RT(1.0) );

        } else if ((irr >= 1) && (irr < 15)) {

            RT p =   ( r1*RT(1.0e6)              - Halls_data_r0col(irr-1) )
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
        RT T_degC = t_crs - RT(273.15); // [K] => [degC]
        RT vis_crs = zero;
        if( T_degC >= zero ) {
            vis_crs = ( RT(1.7180) + RT(4.9E-3)*T_degC ) * RT(1.E-5);
        } else {
            vis_crs = ( RT(1.7180) + RT(4.9E-3)*T_degC -RT(1.2E-5)*T_degC*T_degC ) * RT(1.E-5);
        }

        //air mean free path [m]
        RT lambda_crs = (two*vis_crs) / (p_crs*std::sqrt(RT(8.0)*mwdair*RT(1.E-3)/(PI*R_d*t_crs)));

        // temporary var
        RT dtmp = zero;

        // slip correction of droplets [-]
        dtmp   = RT(1.2570) + RT(0.40) * exp(-RT(0.550)*diameter_1/lambda_crs);
        RT slip_corr_1 = one + (two*lambda_crs*dtmp)/(diameter_1);
        dtmp   = RT(1.2570) + RT(0.40) * exp(-RT(0.550)*diameter_2/lambda_crs);
        RT slip_corr_2 = one + (two*lambda_crs*dtmp)/(diameter_2);

        // diffusion term [m*m/s]
        dtmp = (boltz*t_crs)/(three*PI*vis_crs);
        RT diff_term_1 = dtmp * (slip_corr_1/diameter_1);
        RT diff_term_2 = dtmp * (slip_corr_2/diameter_2);

        // velocity term [m/s]
        dtmp = (RT(8.0)*boltz*t_crs)/PI;
        RT vel_term_1 = std::sqrt(dtmp/a_mass_1);
        RT vel_term_2 = std::sqrt(dtmp/a_mass_2);

        // mean free path of droplets [m]
        dtmp = RT(8.0)/PI;
        RT lambda_1 = dtmp * (diff_term_1/vel_term_1);
        RT lambda_2 = dtmp * (diff_term_2/vel_term_2);

        //length term [m]
        dtmp = (diameter_1+lambda_1)*(diameter_1+lambda_1)*(diameter_1+lambda_1)
               - std::exp(RT(1.5)*std::log(diameter_1*diameter_1+lambda_1*lambda_1));
        RT length_term_1 = dtmp/(three*diameter_1*lambda_1) - diameter_1;
        dtmp = (diameter_2+lambda_2)*(diameter_2+lambda_2)*(diameter_2+lambda_2)
               - std::exp(RT(1.5)*std::log(diameter_2*diameter_2+lambda_2*lambda_2));
        RT length_term_2 = dtmp/(three*diameter_2*lambda_2) - diameter_2;

        // Brownian Coagulation Coefficient K12 [m3/s]
        RT sumdia = diameter_1 + diameter_2;
        RT sumd   = diff_term_1   + diff_term_2;
        RT sumc   = std::sqrt( vel_term_1*vel_term_1 + vel_term_2*vel_term_2 );
        RT sumg   = std::sqrt( two*length_term_1*length_term_1 + two*length_term_2*length_term_2 );

        dtmp = sumdia/(sumdia+two*sumg) + (RT(8.0)*sumd)/(sumdia*sumc);
        RT k12 = two*PI * sumdia*sumd/dtmp;

        return k12;
    }

};

#endif
#endif

