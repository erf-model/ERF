#ifndef ERF_AUXILIARY_DIFFUSION_FACE_H_
#define ERF_AUXILIARY_DIFFUSION_FACE_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace erf_auxiliary {

/**
 * Diffusive semantics of one nonperiodic domain face. A closed face carries no
 * diffusive transfer. A Dirichlet face exchanges with the state held in the
 * ghost cell beyond it, which is the value on the face, as ERF's ext_dir. A
 * DirichletWhereInward face does so where the dry-air carrier is inward or
 * zero and is closed where it is outward, as ERF's ext_dir_upwind.
 */
enum class DiffusionFaceKind : int { Closed, Dirichlet, DirichletWhereInward };

/** Whether one cell of a domain face exchanges; `carrier` is positive along +dir. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
DiffusionFaceIsOpen (const DiffusionFaceKind kind, const bool low_face,
                     const amrex::Real carrier) noexcept
{
    if (kind == DiffusionFaceKind::Dirichlet) { return true; }
    if (kind == DiffusionFaceKind::DirichletWhereInward) {
        return low_face ? carrier >= amrex::Real(0.0) : carrier <= amrex::Real(0.0);
    }
    return false;
}

} // namespace erf_auxiliary

#endif
