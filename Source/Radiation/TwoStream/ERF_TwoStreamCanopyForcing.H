#ifndef ERF_TWO_STREAM_CANOPY_FORCING_H_
#define ERF_TWO_STREAM_CANOPY_FORCING_H_

#include <AMReX_Array4.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <limits>

/**
 * \file ERF_TwoStreamCanopyForcing.H
 * \brief What the two-stream sweep hands the building faces of the
 * immersed-boundary surface energy balance (erf.ibseb.radiation = two_stream),
 * and how a face reads it.
 *
 * The sweep already writes the interface fluxes of every column into ERF's
 * ``rad_fluxes`` (RRTMGP's level layout: index k the lower interface of layer k,
 * component 0 SW up, 1 SW down, 2 LW up, 3 LW down, after the clear/cloudy
 * blend). What it does not keep is the direct beam, which only its scratch holds,
 * and the cosine of the zenith it used. On the levels where ERF::init_ibseb()
 * asks for them (TwoStreamRadiation::supply_canopy_forcing()) every sweep also
 * writes:
 *  - the direct beam on a horizontal surface at the interfaces m = 0 .. m_top of
 *    every column, counted from the surface, m_top the top of the tallest
 *    building of the level (a field on the level's grids flattened in z and then
 *    given those interfaces, index k = m), blended clear/cloudy like the fluxes;
 *  - the cosine of the zenith angle, floored at zero (a 2D field).
 * A face reads the column of its fluid cell at its own height
 * (two_stream_canopy_sample()): the buildings are height maps, so no building
 * stands above a face in its own column, and the downwelling fluxes the column
 * carries there are the light and the sky that reach it, whatever the column
 * does below, where it passes through the buildings.
 */

/**
 * @brief The fields a level's faces read for one step: the beam and the cosine
 * of the zenith the sweep of that step wrote, and that sweep's interface fluxes.
 * All null when the level has not swept this step (the initial report).
 */
struct TwoStreamCanopyView
{
    const amrex::MultiFab* beam = nullptr;        ///< Direct beam on a horizontal surface, interfaces 0 .. m_top [W/m^2]
    const amrex::MultiFab* cos_zenith = nullptr;  ///< Cosine of the zenith of the sweep, floored at zero (2D) [1]
    const amrex::MultiFab* fluxes = nullptr;      ///< The level's rad_fluxes (SW up, SW down, LW up, LW down) [W/m^2]
    [[nodiscard]] bool valid () const noexcept { return beam != nullptr && cos_zenith != nullptr && fluxes != nullptr; }
};

/**
 * @brief Add w times the direct beam on a horizontal surface at the interfaces
 * 0 .. m_top of column (i,j), from the scratch of the evaluation that has just
 * run there, into @p beam (w = 1 with @p set for a clear sky).
 *
 * vertical_two_stream_sweep() keeps the Beer-Lambert beam of every interface in
 * its scratch (TwoStreamScratch::F_DIR), but writes it only when the shortwave is
 * on and the sun is up; otherwise the scratch holds the other evaluation's beam or
 * nothing at all, and the beam is zero. Call it right after each evaluation, before
 * the next one overwrites the scratch.
 *
 * @param scratch    the sweep's scratch, as vertical_two_stream_sweep() takes it
 * @param f_dir_comp the scratch component of the beam (TwoStreamScratch::F_DIR)
 * @param kmin       lowest cell of the column (its surface interface)
 * @param m_top      highest interface of the beam field, at most the column's top
 * @param sw_on      whether the sweep ran its shortwave
 * @param cos_zenith cosine of the zenith angle the sweep used (unfloored)
 * @param w          weight of this evaluation (1 clear only; 1 - f and f with a cloud fraction f)
 * @param set        true for the first evaluation (overwrite), false to add
 * @param beam       the beam field (component 0, index k = m)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void two_stream_canopy_beam (int i, int j, int kmin, int m_top,
                             const amrex::Array4<const amrex::Real>& scratch, int f_dir_comp,
                             bool sw_on, amrex::Real cos_zenith, amrex::Real w, bool set,
                             const amrex::Array4<amrex::Real>& beam)
{
    const bool lit = sw_on && (cos_zenith > amrex::Real(0.0));
    for (int m = 0; m <= m_top; ++m) {
        const amrex::Real b = lit ? w * scratch(i, j, kmin + m, f_dir_comp) : amrex::Real(0.0);
        beam(i, j, m) = set ? b : beam(i, j, m) + b;
    }
}

/**
 * @brief Sample interface of a face and the weight of the one above: a roof (a
 * z face, which with height-map buildings always has its solid below) takes the
 * interface at its own height, the bottom of its fluid cell (m = k - k0); a wall
 * the mean of the two interfaces of its cell (m and m + 1), its mean height.
 *
 * @param k    fluid cell of the face (the level's index)
 * @param k0   the level's domain bottom
 * @param dir  face direction (0 x, 1 y, 2 z)
 * @param w_up weight of interface m + 1
 * @return the lower sample interface m
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int two_stream_canopy_sample (int k, int k0, int dir, amrex::Real& w_up)
{
    w_up = (dir == 2) ? amrex::Real(0.0) : amrex::Real(0.5);
    return k - k0;
}

/** @brief What one face takes from its column (all W/m^2). */
struct TwoStreamFaceSky
{
    amrex::Real dni = 0.0;      ///< Direct-normal irradiance at the face's height
    amrex::Real diffuse = 0.0;  ///< Diffuse shortwave down at the face's height
    amrex::Real lw_down = 0.0;  ///< Longwave down at the face's height
    amrex::Real sw_up_g = 0.0;  ///< Shortwave up at the ground
    amrex::Real lw_up_g = 0.0;  ///< Longwave up at the ground
};

/**
 * @brief The column radiation a face at fluid cell (i,j,k) of direction @p dir
 * takes: at its sample (two_stream_canopy_sample()), the direct-normal
 * irradiance (the beam over the column's cosine of the zenith), the diffuse
 * light (the shortwave down less the beam, floored at zero against round-off on
 * a transparent column) and the longwave down; at the ground (index @p k0 of the
 * fluxes), the shortwave and longwave up.
 *
 * @param beam    beam field (index k = m)
 * @param cosz    cosine of the zenith, floored at zero (2D)
 * @param fluxes  interface fluxes in RRTMGP's level layout (absolute k)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
TwoStreamFaceSky two_stream_face_sky (int i, int j, int k, int k0, int dir,
                                      const amrex::Array4<const amrex::Real>& beam,
                                      const amrex::Array4<const amrex::Real>& cosz,
                                      const amrex::Array4<const amrex::Real>& fluxes)
{
    amrex::Real w_up = 0.0;
    const int ms = two_stream_canopy_sample(k, k0, dir, w_up);
    const int mu = ms + ((w_up > amrex::Real(0.0)) ? 1 : 0);
    const amrex::Real w_lo = amrex::Real(1.0) - w_up;
    const amrex::Real beam_h = w_lo * beam(i, j, ms) + w_up * beam(i, j, mu);
    const amrex::Real sw_dn  = w_lo * fluxes(i, j, k0 + ms, 1) + w_up * fluxes(i, j, k0 + mu, 1);
    TwoStreamFaceSky s;
    const amrex::Real cz = cosz(i, j, 0);
    // Both operands valid in both arms: an optimiser may evaluate the division
    // before the selection, and 0 / 0 at night would raise the invalid flag.
    s.dni     = (cz > amrex::Real(0.0)) ? beam_h / (cz + std::numeric_limits<amrex::Real>::min()) : amrex::Real(0.0);
    s.diffuse = (sw_dn - beam_h > amrex::Real(0.0)) ? sw_dn - beam_h : amrex::Real(0.0);
    s.lw_down = w_lo * fluxes(i, j, k0 + ms, 3) + w_up * fluxes(i, j, k0 + mu, 3);
    s.sw_up_g = fluxes(i, j, k0, 0);
    s.lw_up_g = fluxes(i, j, k0, 2);
    return s;
}

#endif // ERF_TWO_STREAM_CANOPY_FORCING_H_
