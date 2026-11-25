#ifndef ERF_AUXILIARY_DIFFUSION_FACE_H_
#define ERF_AUXILIARY_DIFFUSION_FACE_H_

#include <AMReX_Array4.H>
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

/**
 * Whether the domain face cell (i, j, k), on the low or high face of a
 * direction, exchanges. `carrier` is that direction's dry-air carrier on the
 * face layout, read only on a DirichletWhereInward face. The low- and
 * high-order diffusive rates both take their open faces from here, so a face
 * the limiter blends is open, or closed, in both.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
DiffusionDomainFaceIsOpen (const bool low_face, const DiffusionFaceKind kind_lo,
                           const DiffusionFaceKind kind_hi,
                           const amrex::Array4<const amrex::Real>& carrier,
                           const int i, const int j, const int k) noexcept
{
    const DiffusionFaceKind kind = low_face ? kind_lo : kind_hi;
    const amrex::Real sign = kind == DiffusionFaceKind::DirichletWhereInward
                                 ? carrier(i, j, k, 0) : amrex::Real(0.0);
    return DiffusionFaceIsOpen(kind, low_face, sign);
}

} // namespace erf_auxiliary

#endif
