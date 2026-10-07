#ifndef ERF_TERRAIN_IMPLICIT_METRIC_H_
#define ERF_TERRAIN_IMPLICIT_METRIC_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

/**
 * Terrain-metric part of the vertical diffusion on a terrain-fitted mesh, for
 * erf.implicit_terrain_metric.
 *
 * The flux through a zeta face projects the horizontal fluxes onto the face
 * normal.  Their h*d/dz parts add a vertical-diffusion term K_h * M * d/dz to the
 * zeta flux, with the factors below.  The explicit operator forms it from
 * averaged horizontal fluxes; with the option on, the same term in compact form
 * (face differences in zeta) is also put into the implicit vertical solve and
 * subtracted explicitly at the state the explicit terms see, so the stages
 * with an implicit solve treat the K_h * M * d2/dz2 stiffness implicitly.  The
 * explicit remainder (averaged minus compact) is small for smooth fields on a
 * uniform slope; at the grid scale it is the whole compact term with the
 * opposite sign, cancelled within the stage by the implicit solve, which is why
 * every stage needs the solve (SolverChoice::check_params).
 *
 * h_xi and h_eta are the terrain slopes dz/dxi and dz/deta at the face, mx and my
 * the map factors there.
 */

/**
 * Momentum: the zeta flux of u holds -0.5 * mu_h_tot * M_u * du/dz with
 * M_u = 2 (h_xi mx)^2 + (h_eta my)^2 (from h_xi S11 and h_eta S12); for v,
 * M_v = (h_xi mx)^2 + 2 (h_eta my)^2 (from h_xi S12 and h_eta S22).
 *
 * @param[in] stagdir 0 for u, 1 for v
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
TerrainMetricMomFactor (const int stagdir,
                        const amrex::Real h_xi, const amrex::Real h_eta,
                        const amrex::Real mx, const amrex::Real my) noexcept
{
    const amrex::Real ax = h_xi * mx;
    const amrex::Real ay = h_eta * my;
    return (stagdir == 0) ? amrex::Real(2.0)*ax*ax + ay*ay
                          : ax*ax + amrex::Real(2.0)*ay*ay;
}

/**
 * Scalars: G_z = F_z - h_xi mx Fx_bar - h_eta my Fy_bar holds
 * -K_h * M_s * dphi/dz with M_s = (h_xi mx)^2 + (h_eta my)^2.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
TerrainMetricScalarFactor (const amrex::Real h_xi, const amrex::Real h_eta,
                           const amrex::Real mx, const amrex::Real my) noexcept
{
    const amrex::Real ax = h_xi * mx;
    const amrex::Real ay = h_eta * my;
    return ax*ax + ay*ay;
}

/**
 * The split applies on the faces the implicit solve treats with a coefficient:
 * every zeta face strictly inside the domain (dom_lo < kface <= dom_hi, with
 * dom_hi the top cell).  The bottom and top faces carry the boundary condition
 * (wall, surface layer, Neumann) and stay explicit.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
TerrainMetricImplicitFace (const int kface, const int dom_lo, const int dom_hi) noexcept
{
    return (kface > dom_lo) && (kface <= dom_hi);
}

#endif
