#ifndef SUPERDROPLET_PC_COALESCENCE_H_
#define SUPERDROPLET_PC_COALESCENCE_H_

#include <algorithm>
#include "ERF_Constants.H"

#ifdef ERF_USE_PARTICLES

/*! \brief Coalescence kernels */
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
};

#endif
#endif

