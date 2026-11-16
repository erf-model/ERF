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
 *    every column, counted from the surface, m_top the highest interface a face
 *    of the level reads (a field on the level's grids flattened in z and then
 *    given those interfaces, index k = m), blended clear/cloudy like the fluxes;
 *  - the cosine of the zenith angle, floored at zero (a 2D field).
 * A face reads the column of its fluid cell at its own height
 * (two_stream_canopy_sample()): the downwelling fluxes there are the sun and
 * the sky that reach it (for height-map buildings no building stands above a
 * face in its own column; a face under the overhang of an STL building sees only
 * building in IBFaceSet's view fractions, which take each column as solid up to
 * its top, so its column's fluxes carry no weight), and the
 * upwelling fluxes there are the ground as seen from that height, through the air
 * in between.
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
    /// The direct-beam albedo the sweep used at each column's lower boundary (2D) [1], and, as a
    /// second component, the interface that boundary is at, counted from the ground (0, or a
    /// roof's with erf.ibseb.two_stream_roofs; two_stream_column_base()). Read only with
    /// erf.ibseb.ground_view_shadow, which checks for it; valid() does not ask for it.
    const amrex::MultiFab* ground_albedo_dir = nullptr;
    [[nodiscard]] bool valid () const noexcept { return beam != nullptr && cos_zenith != nullptr && fluxes != nullptr; }
};

/**
 * @brief The interface a column starts on, counted from the ground: 0, or the roof's
 * with erf.ibseb.two_stream_roofs. It is the second component of the field of
 * TwoStreamCanopyView::ground_albedo_dir, which the sweep writes with the albedo of the
 * column's lower boundary.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int two_stream_column_base (int i, int j, const amrex::Array4<const amrex::Real>& alb_dir) noexcept
{
    return static_cast<int>(alb_dir(i, j, 0, 1) + amrex::Real(0.5));
}

/**
 * @brief Write the direct beam on a horizontal surface at the interfaces
 * 0 .. m_top of column (i,j), from the scratch of the evaluation that has just
 * run there, into @p beam: the clear-sky evaluation's beam as it is, and the
 * cloudy evaluation's blended in with the cloud fraction, by the same expression
 * as the fluxes ((1 - f) clear + f cloudy), so the two round alike in single
 * precision and the shortwave down less the beam is the diffuse light.
 *
 * vertical_two_stream_sweep() keeps the Beer-Lambert beam of every interface in
 * its scratch (TwoStreamScratch::F_DIR), but writes it only when the shortwave is
 * on and the sun is up; otherwise the scratch holds the other evaluation's beam or
 * nothing at all, and the beam is zero. Call it right after each evaluation, before
 * the next one overwrites the scratch.
 *
 * @param scratch        the sweep's scratch, as vertical_two_stream_sweep() takes it
 * @param f_dir_comp     the scratch component of the beam (TwoStreamScratch::F_DIR)
 * @param kmin           lowest cell of the column (its surface interface)
 * @param m_top          highest interface of the beam field, at most the column's top
 * @param sw_on          whether the sweep ran its shortwave
 * @param cos_zenith     cosine of the zenith angle the sweep used (unfloored)
 * @param cloud_fraction the cloud fraction (used only when @p blend)
 * @param blend          false for the clear-sky evaluation (overwrite), true for the cloudy one
 * @param beam           the beam field (component 0, index k = m)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void two_stream_canopy_beam (int i, int j, int kmin, int m_top,
                             const amrex::Array4<const amrex::Real>& scratch, int f_dir_comp,
                             bool sw_on, amrex::Real cos_zenith, amrex::Real cloud_fraction, bool blend,
                             const amrex::Array4<amrex::Real>& beam)
{
    const bool lit = sw_on && (cos_zenith > amrex::Real(0.0));
    for (int m = 0; m <= m_top; ++m) {
        const amrex::Real b = lit ? scratch(i, j, kmin + m, f_dir_comp) : amrex::Real(0.0);
        beam(i, j, m) = blend ? (1.0 - cloud_fraction) * beam(i, j, m) + cloud_fraction * b : b;
    }
}

/**
 * @brief Sample interface of a face and the weight of the one above. A z face
 * takes the interface it lies on: a roof (solid below, side -1) the bottom of its
 * fluid cell (m = k - k0), a ceiling (solid above, side +1, as under an overhang
 * of a building from an STL file) the top of it (m = k - k0 + 1). A wall takes the
 * mean of the two interfaces of its cell (m and m + 1), its mean height.
 *
 * @param k    fluid cell of the face (the level's index)
 * @param k0   the level's domain bottom
 * @param dir  face direction (0 x, 1 y, 2 z)
 * @param side side of the fluid cell the solid is on (-1 or +1)
 * @param w_up weight of interface m + 1
 * @return the lower sample interface m
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int two_stream_canopy_sample (int k, int k0, int dir, int side, amrex::Real& w_up)
{
    w_up = (dir == 2) ? amrex::Real(0.0) : amrex::Real(0.5);
    return k - k0 + ((dir == 2 && side > 0) ? 1 : 0);
}

/** @brief What one face takes from its column (all W/m^2). */
struct TwoStreamFaceSky
{
    amrex::Real dni = 0.0;      ///< Direct-normal irradiance at the face's height
    amrex::Real diffuse = 0.0;  ///< Diffuse shortwave down at the face's height
    amrex::Real lw_down = 0.0;  ///< Longwave down at the face's height
    amrex::Real sw_up = 0.0;  ///< Shortwave up at the face's height (the ground's reflection, seen from there)
    amrex::Real lw_up = 0.0;  ///< Longwave up at the face's height (the ground's emission and the air's below)
};

/**
 * @brief The column radiation a face at fluid cell (i,j,k) of direction @p dir,
 * with its solid on side @p side, takes: at its sample (two_stream_canopy_sample()), the direct-normal
 * irradiance (the beam over the column's cosine of the zenith), the diffuse
 * light (the shortwave down less the beam, floored at zero against round-off on
 * a transparent column), the longwave down, and the shortwave and longwave up,
 * which the face takes through its view of the ground.
 *
 * @param beam    beam field (index k = m)
 * @param cosz    cosine of the zenith, floored at zero (2D)
 * @param fluxes  interface fluxes in RRTMGP's level layout (absolute k)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
TwoStreamFaceSky two_stream_face_sky (int i, int j, int k, int k0, int dir, int side,
                                      const amrex::Array4<const amrex::Real>& beam,
                                      const amrex::Array4<const amrex::Real>& cosz,
                                      const amrex::Array4<const amrex::Real>& fluxes)
{
    amrex::Real w_up = 0.0;
    const int ms = two_stream_canopy_sample(k, k0, dir, side, w_up);
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
    s.sw_up = w_lo * fluxes(i, j, k0 + ms, 0) + w_up * fluxes(i, j, k0 + mu, 0);
    s.lw_up = w_lo * fluxes(i, j, k0 + ms, 2) + w_up * fluxes(i, j, k0 + mu, 2);
    return s;
}

/**
 * @brief The share of column (i,j)'s shortwave up that is its ground's reflection of
 * the direct beam: (direct-beam albedo x beam at the surface) over the shortwave up at
 * the surface, which is that reflection plus the ground's reflection of the diffuse
 * light; in [0, 1], and 0 when nothing comes up. The building faces take this share of
 * the shortwave up at their own height as the reflected beam (erf.ibseb.ground_view_shadow):
 * exact where the air between the ground and the face only absorbs, as both parts are
 * then dimmed alike on the way up; with scattering the light the air sends up from the
 * beam and from the sky's diffuse light is split in the surface's proportion too.
 *
 * @param k0      the level's domain bottom (the surface interface of fluxes)
 * @param beam    beam field (index k = m; m = 0 the surface)
 * @param fluxes  interface fluxes in RRTMGP's level layout (absolute k; component 0 SW up)
 * @param alb_dir the ground's direct-beam albedo (2D)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real two_stream_ground_beam_share (int i, int j, int k0,
                                          const amrex::Array4<const amrex::Real>& beam,
                                          const amrex::Array4<const amrex::Real>& fluxes,
                                          const amrex::Array4<const amrex::Real>& alb_dir)
{
    // The column's lower boundary: the ground, or the roof it starts on.
    const int m_base = two_stream_column_base(i, j, alb_dir);
    const amrex::Real up0  = fluxes(i, j, k0 + m_base, 0);
    const amrex::Real refl = alb_dir(i, j, 0) * beam(i, j, m_base);
    // Both operands valid in both arms (see two_stream_face_sky()).
    const amrex::Real share = (up0 > amrex::Real(0.0) && refl > amrex::Real(0.0))
                            ? refl / (up0 + std::numeric_limits<amrex::Real>::min()) : amrex::Real(0.0);
    return (share < amrex::Real(1.0)) ? share : amrex::Real(1.0);
}

/// What two_stream_ground_wall_terms() adds to the ground balance's net shortwave and longwave [W/m^2].
struct GroundWallTerms
{
    amrex::Real d_sw = 0.0;
    amrex::Real d_lw = 0.0;
};

/**
 * @brief The change of the open ground's net radiation beside walls
 * (erf.ibseb.ground_sky_view). The walls take the share @p wall_view of the ground's
 * view of the sky. Through it the ground loses the sky's diffuse light and longwave and
 * gains the walls':
 *   d_sw = c [(1 - alb_dif) sw_walls - wall_view (sw_net - sw_beam)],
 *   d_lw = c emis (lw_walls - wall_view lw_down).
 * The share comes from a count of rays and can exceed 1 in a column closed in by walls;
 * a ground cannot see more than its whole view, so there c = 1 / wall_view scales the
 * loss and the gain alike (the ground then sees no sky and the walls' mean); c = 1
 * otherwise. The direct beam is not touched (the shadow takes it).
 *
 * @param sw_net     the sweep's absorbed shortwave at the ground, beam and diffuse [W/m^2]
 * @param sw_beam    its absorbed direct beam [W/m^2]
 * @param wall_view  the walls' share of the ground's view [1]
 * @param sw_walls   the walls' reflected shortwave arriving on the ground [W/m^2]
 * @param lw_walls   the walls' longwave arriving on the ground [W/m^2]
 * @param lw_down    the sky's longwave down at the ground [W/m^2]
 * @param absorb_dif the ground's diffuse absorptance, 1 - alb_dif [1]
 * @param emis       the ground's longwave emissivity [1]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GroundWallTerms two_stream_ground_wall_terms (amrex::Real sw_net, amrex::Real sw_beam, amrex::Real wall_view,
                                              amrex::Real sw_walls, amrex::Real lw_walls, amrex::Real lw_down,
                                              amrex::Real absorb_dif, amrex::Real emis) noexcept
{
    // The divisor is at least 1 whatever the share: an optimiser that forms both sides of
    // a choice would otherwise divide by a share of zero, which traps where traps are on.
    const amrex::Real whole = (wall_view > amrex::Real(1.0)) ? wall_view : amrex::Real(1.0);
    const amrex::Real c = amrex::Real(1.0) / whole;
    GroundWallTerms t;
    t.d_sw = c * (absorb_dif * sw_walls - wall_view * (sw_net - sw_beam));
    t.d_lw = c * emis * (lw_walls - wall_view * lw_down);
    return t;
}

#endif // ERF_TWO_STREAM_CANOPY_FORCING_H_
