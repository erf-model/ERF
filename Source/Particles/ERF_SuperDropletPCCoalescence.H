#ifndef SUPERDROPLET_PC_COALESCENCE_H_
#define SUPERDROPLET_PC_COALESCENCE_H_

#include <algorithm>
#include "ERF_Constants.H"
#include "ERF_SuperDropletPCDefinitions.H"

#ifdef ERF_USE_PARTICLES

/*! \brief Field indices for 2-field interpolation (ice-water collision efficiency) */
AMREX_ENUM(InterpFieldsTR,
    temperature,
    moist_density,
    NUM_FIELDS
);

/*! \brief Field indices for 4-field interpolation (riming) */
AMREX_ENUM(InterpFieldsFull,
    temperature,
    pressure,
    moist_density,
    qv,
    NUM_FIELDS
);

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
    const RT lambda_mfp = RT(6.62e-8); // mean free path of air
    const RT para_Asc = RT(2.51);
    const RT para_Bsc = RT(0.8);
    const RT para_Csc = RT(0.55);
    const RT k_coeff = RT(1.0);

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
     * \param[in] a_x1 Pointer to position vector of particle 1 (m)
     * \param[in] a_x2 Pointer to position vector of particle 2 (m)
     * \param[in] a_vel_type Kernel relative velocity type
     * \return Collision kernel value (m³/s)
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT sedimentation (  const RT a_r1,
                        const RT a_r2,
                        const RT* const a_v1,
                        const RT* const a_v2,
                        const RT* const a_x1,
                        const RT* const a_x2,
                        const SDKernelRelativeVelocityType a_vel_type ) const noexcept
    {
        auto p = std::min(a_r1,a_r2) / std::max(a_r1,a_r2);
        auto E = RT(0.5)*p*p / ((RT(1.0)+p)*(RT(1.0)+p));

        RT dv = 0;
        if (a_vel_type == SDKernelRelativeVelocityType::radial_velocity) {
            RT dr = 0, dv2 = 0;
            for (int d = 0; d < SPACEDIM; d++) {
                const RT ddv = a_v1[d]-a_v2[d];
                const RT ddx = a_x1[d]-a_x2[d];
                dv  += ddv*ddx;
                dr  += ddx*ddx;
                dv2 += ddv*ddv;
            }
            // radial component; fall back to the full relative speed for coincident pairs
            dv = (dr > RT(0)) ? std::abs(dv) / std::sqrt(dr) : std::sqrt(dv2);
        } else {
            for (int d = 0; d < SPACEDIM; d++) {
                dv += (a_v1[d]-a_v2[d])*(a_v1[d]-a_v2[d]);
            }
            dv = std::sqrt(dv);
        }

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
     * \param[in] a_x1 Pointer to position vector of particle 1 (m)
     * \param[in] a_x2 Pointer to position vector of particle 2 (m)
     * \param[in] a_vel_type Kernel relative velocity type
     * \return Collision kernel value (m³/s)
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT Longs ( const RT a_r1,
               const RT a_r2,
               const RT* const a_v1,
               const RT* const a_v2,
               const RT* const a_x1,
               const RT* const a_x2,
               const SDKernelRelativeVelocityType a_vel_type ) const noexcept
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
        if (a_vel_type == SDKernelRelativeVelocityType::radial_velocity) {
            RT dr = 0, dv2 = 0;
            for (int d = 0; d < SPACEDIM; d++) {
                const RT ddv = a_v1[d]-a_v2[d];
                const RT ddx = a_x1[d]-a_x2[d];
                dv  += ddv*ddx;
                dr  += ddx*ddx;
                dv2 += ddv*ddv;
            }
            // radial component; fall back to the full relative speed for coincident pairs
            dv = (dr > RT(0)) ? std::abs(dv) / std::sqrt(dr) : std::sqrt(dv2);
        } else {
            for (int d = 0; d < SPACEDIM; d++) {
                dv += (a_v1[d]-a_v2[d])*(a_v1[d]-a_v2[d]);
            }
            dv = std::sqrt(dv);
        }

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
               const RT* const a_v2, /*!< velocity of particle 2 */
               const RT* const a_x1, /*!< position of particle 1 */
               const RT* const a_x2, /*!< position of particle 2 */
               const SDKernelRelativeVelocityType a_vel_type /*!< velocity type */) const noexcept
    {
        auto r1 = std::max( a_r1, a_r2 );  // large
        auto r2 = std::min( a_r1, a_r2 );  // small
        auto ratio = r2 / r1;
        auto sumr = r1 + r2;

        int irr = Halls_data_r0col(r1);
        int iqq = Halls_data_ratcol(ratio);

        RT c_rate = RT(1.0);

        if( irr >= 15 ) {

            RT q  =   ( ratio                  - Halls_data_ratcol(iqq-1) )
                    / ( Halls_data_ratcol(iqq) - Halls_data_ratcol(iqq-1) );

            RT ek =   (RT(1.0)-q) * Halls_data_ecoll(iqq-1,14)
                    + q       * Halls_data_ecoll(iqq  ,14);

            c_rate = std::min( ek, RT(1.0) );

        } else if ((irr >= 1) && (irr < 15)) {

            RT p =   ( r1*RT(1.0e6)              - Halls_data_r0col(irr-1) )
                   / ( Halls_data_r0col(irr) - Halls_data_r0col(irr-1) );

            RT q =   ( ratio                  - Halls_data_ratcol(iqq-1) )
                   / ( Halls_data_ratcol(iqq) - Halls_data_ratcol(iqq-1) );

            c_rate =   (RT(1.0)-p) * (RT(1.0)-q) * Halls_data_ecoll(iqq-1,irr-1)
                     + p       * (RT(1.0)-q) * Halls_data_ecoll(iqq-1,irr  )
                     + (RT(1.0)-p) * q       * Halls_data_ecoll(iqq  ,irr-1)
                     + p       * q       * Halls_data_ecoll(iqq  ,irr  );

        } else {

            RT q =  ( ratio                  - Halls_data_ratcol(iqq-1) )
                  / ( Halls_data_ratcol(iqq) - Halls_data_ratcol(iqq-1) );

            c_rate =   (RT(1.0)-q) * Halls_data_ecoll(iqq-1,0)
                     + q       * Halls_data_ecoll(iqq  ,0);

        }
        c_rate *= (PI*sumr*sumr);

        RT dv = 0;
        if (a_vel_type == SDKernelRelativeVelocityType::radial_velocity) {
            RT dr = 0, dv2 = 0;
            for (int d = 0; d < SPACEDIM; d++) {
                const RT ddv = a_v1[d]-a_v2[d];
                const RT ddx = a_x1[d]-a_x2[d];
                dv  += ddv*ddx;
                dr  += ddx*ddx;
                dv2 += ddv*ddv;
            }
            // radial component; fall back to the full relative speed for coincident pairs
            dv = (dr > RT(0)) ? std::abs(dv) / std::sqrt(dr) : std::sqrt(dv2);
        } else {
            for (int d = 0; d < SPACEDIM; d++) {
                dv += (a_v1[d]-a_v2[d])*(a_v1[d]-a_v2[d]);
            }
            dv = std::sqrt(dv);
        }

        return c_rate*dv;
    }

    /* The following code is adapted from SCALE-SDM:
     * Copyright (c) 2012-2014, Team SCALE, All rights reserved. */
    /*! \brief Brownian coagulation kernel from Seinfeld and Pandis
     *
     *  Computes the Brownian coagulation coefficient K12 [m³/s] for two particles
     *  using the interpolation formula that accounts for both continuum and
     *  free-molecular regimes.
     *
     *  Reference: Seinfeld, J. H., and S. N. Pandis, 2016: Atmospheric Chemistry
     *  and Physics: From Air Pollution to Climate Change, 3rd ed. Wiley, 1152 pp.
     *  Chapter 13, Eq. 13.56 (Fuchs interpolation formula)
     *
     *  Dynamic viscosity from Pruppacher & Klett (1997), Eq. 3-11.
     *  Slip correction coefficients: 1.257, 0.40, 0.55 (Table 13.1)
     */
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
        RT vis_crs = RT(0.0);
        if( T_degC >= RT(0.0) ) {
            vis_crs = ( RT(1.7180) + RT(4.9E-3)*T_degC ) * RT(1.E-5);
        } else {
            vis_crs = ( RT(1.7180) + RT(4.9E-3)*T_degC -RT(1.2E-5)*T_degC*T_degC ) * RT(1.E-5);
        }

        //air mean free path [m]
        RT lambda_crs = (two*vis_crs) / (p_crs*std::sqrt(RT(8.0)*mwdair/(PI*boltz*avogadro*t_crs)));

        // temporary var
        RT dtmp = RT(0.0);

        // slip correction of droplets [-]
        dtmp   = RT(1.2570) + RT(0.40) * exp(-RT(0.550)*diameter_1/lambda_crs);
        RT slip_corr_1 = RT(1.0) + (RT(2.0)*lambda_crs*dtmp)/(diameter_1);
        dtmp   = RT(1.2570) + RT(0.40) * exp(-RT(0.550)*diameter_2/lambda_crs);
        RT slip_corr_2 = RT(1.0) + (RT(2.0)*lambda_crs*dtmp)/(diameter_2);

        // diffusion term [m*m/s]
        dtmp = (boltz*t_crs)/(RT(3.0)*PI*vis_crs);
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
        RT d1l1_sq = diameter_1*diameter_1+lambda_1*lambda_1;
        dtmp = (diameter_1+lambda_1)*(diameter_1+lambda_1)*(diameter_1+lambda_1)
               - d1l1_sq*std::sqrt(d1l1_sq);  // x^1.5 = x*sqrt(x)
        RT length_term_1 = dtmp/(RT(3.0)*diameter_1*lambda_1) - diameter_1;
        RT d2l2_sq = diameter_2*diameter_2+lambda_2*lambda_2;
        dtmp = (diameter_2+lambda_2)*(diameter_2+lambda_2)*(diameter_2+lambda_2)
               - d2l2_sq*std::sqrt(d2l2_sq);  // x^1.5 = x*sqrt(x)
        RT length_term_2 = dtmp/(RT(3.0)*diameter_2*lambda_2) - diameter_2;

        // Brownian Coagulation Coefficient K12 [m3/s]
        RT sumdia = diameter_1 + diameter_2;
        RT sumd   = diff_term_1   + diff_term_2;
        RT sumc   = std::sqrt( vel_term_1*vel_term_1 + vel_term_2*vel_term_2 );
        RT sumg   = std::sqrt( RT(2.0)*length_term_1*length_term_1 + RT(2.0)*length_term_2*length_term_2 );

        dtmp = sumdia/(sumdia+RT(2.0)*sumg) + (RT(8.0)*sumd)/(sumdia*sumc);
        RT k12 = RT(2.0)*PI * sumdia*sumd/dtmp;

        return k12;
    }

    /* The following code is adapted from SCALE-SDM:
     * Copyright (c) 2012-2014, Team SCALE, All rights reserved. */
    /*! \brief Collision efficiency for riming from Beard and Grover (1974)
     *
     *  Computes the collision efficiency between ice particles and cloud droplets
     *  based on Reynolds number, Froude number, and size ratio.
     *
     *  Reference: Beard, K. V., and S. N. Grover, 1974: Numerical collision
     *  efficiencies for small raindrops colliding with micron size particles.
     *  J. Atmos. Sci., 31, 543-550.
     *  https://doi.org/10.1175/1520-0469(1974)031<0543:NCEFSR>2.0.CO;2
     *
     *  Coefficients:
     *  - K0 fit: A0=-0.1007, A1=-0.358, A2=0.0261 (Re ≤ 400), K0=0.21 (Re > 400)
     *  - yc0 fit: B0=0.1465, B1=1.302, B2=-0.607, B3=0.293
     *  - E = ((yc0+p)/(1+p))² where p = size ratio
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT BeardGrover1974 ( const RT a_sr, /*!< Size ratio (r_droplet/r_ice) */
                         const RT a_Re, /*!< Reynolds number of ice particle */
                         const RT a_Fr  /*!< Froude number */ ) const noexcept
    {
        // riming efficiency of Beard and Grover (1974)
        const RT A0 = -RT(0.1007);
        const RT A1 = -RT(0.358);
        const RT A2 =  RT(0.0261);
        const RT B0 =  RT(0.1465);
        const RT B1 =  RT(1.302);
        const RT B2 = -RT(0.607);
        const RT B3 =  RT(0.293);

        RT K0 = RT(0.0);
        if ( a_Re > RT(400.0)) {
            K0 = RT(0.21);
        } else {
            auto F = std::log(a_Re);
            auto G = A0 + A1*F + A2*F*F;
            K0 = std::exp(G);
        }

        auto z = std::log(a_Fr/K0);
        auto H = std::max(B0 + B1*z + B2*z*z + B3*z*z*z,RT(0.0));
        auto yc0 = (RT(2.0)/PI)*std::atan(H);

        auto retval = (yc0+a_sr)*(yc0+a_sr) / ((RT(1.0)+a_sr)*(RT(1.0)+a_sr));
        return retval;
    }

    /* The following code is adapted from SCALE-SDM:
     * Copyright (c) 2012-2014, Team SCALE, All rights reserved. */
    /*! \brief Collision efficiency for plate-like ice crystals from Erfani and Mitchell (2017)
     *
     *  Computes the collision efficiency between plate-like ice crystals and
     *  cloud droplets based on Reynolds and Froude numbers.
     *
     *  Reference: Erfani, E., and D. L. Mitchell, 2017: Growth of ice particle
     *  mass and projected area during riming. Atmos. Chem. Phys., 17, 1241-1257.
     *  https://doi.org/10.5194/acp-17-1241-2017
     *
     *  Valid for Re ≥ 1, with polynomial fits for K_t (transition Froude number)
     *  and K_c (critical Froude number). Three regimes based on Fr vs K_t.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT ErfaniMitchell2017_Plate ( const RT a_Re, /*!< Reynolds number of ice crystal */
                                  const RT a_Fr  /*!< Froude number */ ) const noexcept
    {
        if (a_Re < RT(1.0)) {
            return RT(0.0);
        } else {

            auto Re_l = std::min(a_Re,RT(120.0));
            auto Fr_l = std::min(a_Fr,RT(35.0));
            auto Re_l_2 = Re_l*Re_l;
            auto Re_l_3 = Re_l_2*Re_l;
            auto Re_l_4 = Re_l_2*Re_l_2;
            auto Re_l_5 = Re_l_4*Re_l;
            auto K_t = -RT(5.07e-10)*Re_l_5
                       +RT(1.73e-7) *Re_l_4
                       -RT(2.17e-5) *Re_l_3
                       +RT(0.0013)  *Re_l_2
                       -RT(0.037)   *Re_l
                       +RT(0.8355);

            RT K_c = RT(0.0);
            if      (Re_l <=  RT(10.0)) { K_c = RT(1.250) * std::pow(Re_l, RT(-0.350)); }
            else if (Re_l <=  RT(40.0)) { K_c = RT(1.072) * std::pow(Re_l, RT(-0.301)); }
            else if (Re_l <= RT(120.0)) { K_c = RT(0.356) * std::pow(Re_l, RT(-0.003)); }

            auto retval = RT(0.0);
            if (Fr_l <= RT(0.35)) {
                retval = RT(0.787) * std::pow(Fr_l, RT(0.988)) * std::max(RT(0.263)*std::log(Re_l)-RT(0.264), RT(0.0));
            } else if(Fr_l <= K_t) {
                retval = (RT(0.7475)*std::log10(Fr_l)+RT(0.620)) * std::max(RT(0.263)*std::log(Re_l)-RT(0.264), RT(0.0));
            } else {
                auto E_1 = (RT(0.7475)*std::log10(K_t)+RT(0.620)) * std::max(RT(0.263)*std::log(Re_l)-RT(0.264), RT(0.0));
                auto E_2 = std::sqrt( std::max(RT(1.0) - RT(0.2)*(std::log10(Fr_l/K_c)-std::sqrt(RT(5.0)))
                                                        *(std::log10(Fr_l/K_c)-std::sqrt(RT(5.0))),
                                               RT(0.0)) );
                retval = std::max(E_1,E_2);
            }
            return retval;
        }
    }

    /* The following code is adapted from SCALE-SDM:
     * Copyright (c) 2012-2014, Team SCALE, All rights reserved. */
    /*! \brief Collision efficiency for columnar ice crystals from Erfani and Mitchell (2017)
     *
     *  Computes the collision efficiency between columnar ice crystals and
     *  cloud droplets based on Reynolds and Froude numbers.
     *
     *  Reference: Erfani, E., and D. L. Mitchell, 2017: Growth of ice particle
     *  mass and projected area during riming. Atmos. Chem. Phys., 17, 1241-1257.
     *  https://doi.org/10.5194/acp-17-1241-2017
     *
     *  Valid for Re ≥ 0.2. Uses different polynomial fits than the plate version
     *  for K_t (transition Froude) and K_c (critical Froude). Two regimes based on Fr vs K_t.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT ErfaniMitchell2017_Column ( const RT a_Re, /*!< Reynolds number of ice crystal */
                                   const RT a_Fr  /*!< Froude number */ ) const noexcept
    {
        if (a_Re < RT(0.2)) {
            return RT(0.0);
        } else {

            auto Re_l = std::min(a_Re,RT(20.0));
            auto Fr_l = std::min(a_Fr,RT(20.0));

            auto K_t = RT(0.0);
            if (Re_l <= RT(2.0)) {
                K_t =   RT(0.0251)*Re_l*Re_l
                      - RT(0.0144)*Re_l
                      + RT(0.811);
            } else {
                K_t = -RT(0.0003)*Re_l*Re_l*Re_l
                      +RT(0.0124)*Re_l*Re_l
                      -RT(0.1634)*Re_l
                      +RT(1.0075);
            }

            auto para_r = RT(0.0), K_c = RT(0.0);
            if (Re_l <= RT(1.7)) {
                para_r = RT(0.7422) * std::pow(Re_l, RT(0.2111));
                K_c = RT(0.7797) * std::pow(Re_l, RT(-0.009));
            } else {
                para_r = RT(0.8025) * std::pow(Re_l, RT(0.0604));
                K_c = RT(1.0916) * std::pow(Re_l, RT(-0.635));
            }

            auto retval = RT(0.0);
            if (Fr_l <= K_t) {
                if (Re_l <= RT(3.0)) {
                    retval= RT(0.787) * std::pow(Fr_l, RT(0.988)) * (-RT(0.0121)*Re_l*Re_l + RT(0.1297)*Re_l + RT(0.0598));
                } else {
                    retval = RT(0.787) * std::pow(Fr_l, RT(0.988)) * (-RT(0.0005)*Re_l*Re_l + RT(0.1028)*Re_l + RT(0.0359));
                }
            } else {
                auto tmp = std::log10(Fr_l/K_c)-std::sqrt(RT(3.5));
                retval = para_r * std::sqrt( std::max(RT(1.0)-tmp*tmp/RT(3.5), RT(0.0)) );
            }
            return retval;
        }
    }

};

#endif
#endif
