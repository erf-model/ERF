#ifndef SUPERDROPLET_PC_COALESCENCE_H_
#define SUPERDROPLET_PC_COALESCENCE_H_

#include <algorithm>
#include "ERF_Constants.H"

#ifdef ERF_USE_PARTICLES

/*! \brief Compute collision cross-section for particles undergoing sedimentation*/
template<typename RT>
struct CollisionCS_Sedimentation
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT operator() ( const RT a_r1,
                    const RT a_r2 ) const noexcept
    {
        auto p = std::min(a_r1,a_r2) / std::max(a_r1,a_r2);
        auto E = 0.5*p*p / ((1.0+p)*(1.0+p));
        return PI*(a_r1+a_r2)*(a_r1+a_r2)*E;
    }
};

#endif
#endif

