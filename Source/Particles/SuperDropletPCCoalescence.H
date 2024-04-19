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
                 const RT a_r2, /*!< radius of droplet 2 */
                 const RT* const /*a_v1*/,
                 const RT* const /*a_v2*/ ) const noexcept
    {
        RT b = 1.5e03;
        auto X_1 = (4.0*PI/3.0) * a_r1*a_r1*a_r1;
        auto X_2 = (4.0*PI/3.0) * a_r2*a_r2*a_r2;
        return b * (X_1 + X_2);
    }

    /*! \brief Sedimentation kernel */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT sedimentation (  const RT a_r1, /*!< radius of particle 1 */
                        const RT a_r2, /*!< radius of particle 2 */
                        const RT* const a_v1, /*!< velocity of particle 1 */
                        const RT* const a_v2  /*!< velocity of particle 2 */) const noexcept
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

    /*! \brief Long's kernel */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT Longs ( const RT a_r1, /*!< radius of particle 1 */
               const RT a_r2, /*!< radius of particle 2 */
               const RT* const a_v1, /*!< velocity of particle 1 */
               const RT* const a_v2  /*!< velocity of particle 2 */) const noexcept
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

            RT p =   ( a_r1*1.0e6            - Halls_data_r0col(irr-1) )
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

};

#endif
#endif

