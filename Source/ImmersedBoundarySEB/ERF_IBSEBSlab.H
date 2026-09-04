#ifndef ERF_IBSEB_SLAB_H
#define ERF_IBSEB_SLAB_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

/**
 * \file ERF_IBSEBSlab.H
 * \brief One-dimensional heat conduction through a wall or roof slab.
 *
 * Adapted from the ERF-SLUCM branch's ERF_UCMSlabConduction.H (implicit
 * Euler, Thomas algorithm, all-plus tridiagonal convention). The one change
 * is the top boundary: the SLUCM solver takes a flux at the top, this one
 * takes the skin temperature, which is what the face balance holds, and
 * returns the conduction into the slab that results. The bottom boundary
 * is the interior temperature of the building.
 *
 * Layers are uniform, ``dz = L / N``, with layer centres at
 * ``(l + 1/2) dz`` below the skin; the skin sits half a layer above the
 * first centre and the interior half a layer below the last, so both
 * boundary fluxes use the half-layer spacing ``2 k / dz``.
 */
namespace ibseb {

/// Largest number of layers the stack arrays allow.
constexpr int SLAB_MAX_LAYERS = 256;

/**
 * Advance one slab column by one step and return the conduction into it.
 *
 * Solves, for layers ``l = 0 .. N-1`` with implicit Euler,
 * ``rho_cp dz (T_l^{n+1} - T_l^n) / dt = flux_in - flux_out``, the fluxes
 * between neighbouring centres being ``k (T_{l-1} - T_l) / dz``, the flux
 * from the skin ``2 k (T_skin - T_0) / dz`` and the flux to the interior
 * ``2 k (T_{N-1} - T_int) / dz``. Unconditionally stable.
 *
 * @param[in,out] T       Layer temperatures [K], ``N`` values, top first.
 * @param[in]  T_skin     Skin temperature at the top [K] (Dirichlet).
 * @param[in]  T_int      Interior temperature at the bottom [K] (Dirichlet).
 * @param[in]  k          Thermal conductivity [W/m/K].
 * @param[in]  rho_cp     Volumetric heat capacity [J/m3/K].
 * @param[in]  dz         Layer thickness [m].
 * @param[in]  dt         Time step [s].
 * @param[in]  N          Number of layers, at most SLAB_MAX_LAYERS.
 * @return Conduction into the slab through the skin at the new time,
 *         ``2 k (T_skin - T_0^{n+1}) / dz`` [W/m2], positive into the slab.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real advance_slab_dirichlet (amrex::Real* T, amrex::Real T_skin, amrex::Real T_int,
                                    amrex::Real k, amrex::Real rho_cp, amrex::Real dz,
                                    amrex::Real dt, int N)
{
    if (N < 1 || N > SLAB_MAX_LAYERS) return 0.0;
    const amrex::Real Fo = k * dt / (rho_cp * dz * dz);   // interior Fourier number
    // Coefficients of a l T_{l-1} + b_l T_l + c_l T_{l+1} = d_l (all-plus form).
    amrex::Real alpha[SLAB_MAX_LAYERS], gamma[SLAB_MAX_LAYERS];
    // Row 0: skin above at half spacing (coefficient 2 Fo), layer 1 below.
    {
        const amrex::Real b0 = 1.0 + 2.0 * Fo + ((N > 1) ? Fo : 2.0 * Fo);
        const amrex::Real c0 = (N > 1) ? -Fo : 0.0;
        const amrex::Real d0 = T[0] + 2.0 * Fo * T_skin + ((N > 1) ? 0.0 : 2.0 * Fo * T_int);
        alpha[0] = c0 / b0;
        gamma[0] = d0 / b0;
    }
    // Interior rows.
    for (int l = 1; l < N - 1; ++l) {
        const amrex::Real a = -Fo, b = 1.0 + 2.0 * Fo, c = -Fo;
        const amrex::Real bm = b - a * alpha[l - 1];
        gamma[l] = (T[l] - a * gamma[l - 1]) / bm;
        alpha[l] = c / bm;
    }
    // Last row (N > 1): interior below at half spacing (coefficient 2 Fo).
    if (N > 1) {
        const int l = N - 1;
        const amrex::Real a = -Fo, b = 1.0 + Fo + 2.0 * Fo;
        const amrex::Real bm = b - a * alpha[l - 1];
        gamma[l] = (T[l] + 2.0 * Fo * T_int - a * gamma[l - 1]) / bm;
        alpha[l] = 0.0;
    }
    // Back substitution.
    T[N - 1] = gamma[N - 1];
    for (int l = N - 2; l >= 0; --l) { T[l] = gamma[l] - alpha[l] * T[l + 1]; }
    return 2.0 * k * (T_skin - T[0]) / dz;
}

} // namespace ibseb

#endif
