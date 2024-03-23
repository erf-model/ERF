#ifndef SUPERDROPLET_PC_COALESCENCE_H_
#define SUPERDROPLET_PC_COALESCENCE_H_

#include <algorithm>
#include "ERF_Constants.H"

#ifdef ERF_USE_PARTICLES

/*! \brief Compute collision cross-section for particles undergoing sedimentation*/
template<typename RT>
struct CollisionKernel_Golovin
{
    RT m_b;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT operator() ( const RT a_r1,
                    const RT a_r2,
                    const RT* const /*a_v1*/,
                    const RT* const /*a_v2*/ ) const noexcept
    {
        auto X_1 = (4.0*PI/3.0) * a_r1*a_r1*a_r1;
        auto X_2 = (4.0*PI/3.0) * a_r2*a_r2*a_r2;
        return m_b * (X_1 + X_2);
    }
};

/*! \brief Compute collision cross-section for particles undergoing sedimentation*/
template<typename RT,int SPACEDIM>
struct CollisionKernel_Sedimentation
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT operator() ( const RT a_r1,
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
};

#endif
#endif

