#ifndef ERF_ACTUATOR_SAMPLING_H
#define ERF_ACTUATOR_SAMPLING_H

/**
 * Interpolation of the flow and of the terrain height at arbitrary points, and the check that
 * points lie on a level's grids.
 */

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include "ERF_ActuatorDeviceBuffers.H"

namespace erf_actuator {

/// Abort, naming who, unless every rank passes the same number of points (one collective).
void check_same_points_everywhere (int npts, const char* who);

/**
 * Sample the velocity at a set of points from the face-centred velocity fields.
 *
 * Each component is interpolated trilinearly from its own staggered grid: u from the x faces,
 * v from the y faces, w from the z faces, so a field linear in x, y and z is reproduced
 * exactly. The vertical weight uses the physical face heights (z_phys_nd on a stretched or
 * terrain-following mesh, prob_lo + k dz otherwise), so the same holds on those meshes as long
 * as the field is linear in the physical height.
 *
 * Every rank interpolates the points whose containing cell lies in one of its valid boxes,
 * reading the neighbours from the ghost cells, and the results are summed across ranks; the
 * velocity fields must therefore have at least one filled ghost cell, and z_phys_nd one ghost
 * node, on the velocities' boxes and distribution. The ghost cells across the domain's lateral
 * boundaries are read too (a point within half a cell of a non-periodic face takes up to half
 * its value from them), so the fields must be filled there as well (FillPatch, not only
 * FillBoundary). In a periodic direction a point is moved by whole periods into the domain; a
 * point on the domain's upper x, y or z face belongs to the last cell. A point outside the
 * domain otherwise, or one found by no rank or by more than one, aborts
 * with its coordinates: the caller owns the points and must keep them inside the domain. Every
 * rank must pass the same points: the reduction is collective and the call stops when the ranks
 * disagree on their number. A model whose ranks hold different points (those of the bodies
 * each owns) calls sample_velocity_local() instead, which has no such check.
 *
 * @param[in] U,V,W face-centred velocities on one level, ghost cells filled (m/s)
 * @param[in] z_phys_nd nodal physical heights (m), or nullptr on a uniform-dz mesh
 * @param[in] geom the level's geometry
 * @param[in] pos 3*n point coordinates (x0,y0,z0, x1,y1,z1, ...) (m, ERF domain coordinates)
 * @param[out] vel 3*n sampled velocities (m/s), the same on every rank on return
 */
void sample_velocity (const amrex::MultiFab& U,
                      const amrex::MultiFab& V,
                      const amrex::MultiFab& W,
                      const amrex::MultiFab* z_phys_nd,
                      const amrex::Geometry& geom,
                      const std::vector<amrex::Real>& pos,
                      std::vector<amrex::Real>& vel);

/**
 * The part of sample_velocity() that a rank does on its own: the velocity at each of the points that
 * one of this rank's boxes holds (cnt 1), zero with cnt 0 at the others. No reduction and no check:
 * the caller combines the ranks' answers. The points must already be moved into the domain across
 * its periodic seams (wrap_periodic).
 * @param[in] wpos 3*n point coordinates inside the domain (m)
 * @param[out] vel 3*n velocities (m/s), zero at the points no box of this rank holds
 * @param[out] cnt n counts: how many of this rank's boxes hold each point (0 or 1)
 * @param buffers device memory kept from call to call (the call uses its own when nullptr); the points
 *        go to the device and the answers come back with asynchronous copies and one synchronisation
 */
void sample_velocity_local (const amrex::MultiFab& U,
                            const amrex::MultiFab& V,
                            const amrex::MultiFab& W,
                            const amrex::MultiFab* z_phys_nd,
                            const amrex::Geometry& geom,
                            const std::vector<amrex::Real>& wpos,
                            std::vector<amrex::Real>& vel,
                            std::vector<int>& cnt,
                            ActuatorDeviceBuffers* buffers = nullptr);

/**
 * The points (3 coordinates each, m) with their coordinates in every periodic direction moved by
 * whole periods into [prob_lo, prob_hi): the image the boxes hold, so that a point that has crossed
 * a periodic seam is sampled on the far side.
 */
std::vector<amrex::Real> wrap_periodic (const std::vector<amrex::Real>& pos, const amrex::Geometry& geom);

/**
 * Sample one component of a cell-centred field at a set of points: bilinear between the four
 * cell centres around the point in the horizontal and linear in the physical height between
 * the two cell centres that bracket it in the column (cell-centre heights from z_phys_nd on a
 * stretched or terrain-following mesh); below the first and above the last cell centre the
 * nearest pair extrapolates linearly. Same ownership rules and reduction as sample_velocity();
 * the field needs one filled ghost cell.
 *
 * @param[in] mf the cell-centred field
 * @param[in] comp the component of mf (0-based)
 * @param[in] pos 3*n point coordinates (m, ERF domain coordinates)
 * @param[out] val n sampled values, the same on every rank on return
 */
void sample_cell_scalar (const amrex::MultiFab& mf, int comp,
                         const amrex::MultiFab* z_phys_nd,
                         const amrex::Geometry& geom,
                         const std::vector<amrex::Real>& pos,
                         std::vector<amrex::Real>& val);

/**
 * Return the height of the terrain surface (the k = 0 node plane of z_phys_nd, bilinear between
 * the four mesh nodes around each point's (x, y)) under a set of points; prob_lo z everywhere on
 * a uniform-dz mesh (z_phys_nd null), where no point is checked. Otherwise the same ownership
 * rules, periodic images and reduction as sample_velocity(): the points must lie inside the domain
 * horizontally, a point on the domain's upper x or y face belonging to the last cell; z is ignored.
 *
 * @param[in] pos 3*n positions (m)
 * @param[out] h n terrain heights (m), the same on every rank on return
 */
void terrain_heights (const amrex::MultiFab* z_phys_nd, const amrex::Geometry& geom,
                      const std::vector<amrex::Real>& pos, std::vector<amrex::Real>& h);

/**
 * Return whether every point, with the cells within `reach` metres around it, lies inside the
 * union of boxes of a level (the anchor level's grids): the sampler needs the neighbouring cells
 * and the spreader the kernel's reach, and on a refined level that is not the whole domain. In a
 * periodic direction the needed cells wrap across the seam to the domain's other side; in a
 * non-periodic one they are clamped to the domain. Returns false and names the first point
 * outside in `first_outside`. On a terrain-following or stretched mesh a point's k is not its
 * height over dz, so z is checked as `z` says: Reach (the cells within reach of the height over dz,
 * for a flat, uniform mesh), Footprint (the reach across only, at any height: the sampler then finds
 * a point above or below the grids by itself), or Column (for the spreading, whose kernel would
 * otherwise be cut where the grids stop: the reach across over every cell the heights within reach
 * may lie in, over the point's own columns from `cells` (column_cell_ranges, one range per point),
 * else whatever the column from the mesh's bounds `zb` (mesh_z_bounds), or over the whole height
 * without either). With `cells` a point whose range is empty (kb > kt) is not covered.
 * Not with a periodic z.
 */
enum class CoverZ { Reach, Footprint, Column };

/**
 * The lowest and highest ground and the smallest and largest spacing of the levels (m) of a mesh,
 * and, over one rank's boxes (mesh_z_bounds with over_all_ranks false) or the level's (mesh_z_levels), a table of each node level's
 * lowest and highest height and of the cells those boxes hold, which bounds the cells a height may
 * lie in whether or not the boxes reach the ground (a column split in z, a level aloft).
 */
struct ZBounds {
    amrex::Real ground_lo{0.0}, ground_hi{0.0}, dz_lo{0.0}, dz_hi{0.0};
    bool all_ground{true};   //!< every box of the level reaches the ground (else the bounds hold for those that do)
    int k0{0};                                  //!< the node level of node_lo[0] and the cell of has_cell[0]
    std::vector<amrex::Real> node_lo, node_hi;  //!< node level k0 + n: the lowest and highest height (empty: no table)
    std::vector<char> has_cell;                 //!< cell k0 + n is in one of the boxes

    /** Whether the level's grids reach the ground (a refined level aloft does not, and has no bounds). */
    bool has_ground () const { return ground_lo <= ground_hi && dz_lo > amrex::Real(0.0); }
    /** Whether the bounds come with the table of node heights: of one rank's boxes (mesh_z_bounds) or of
     *  every rank's (mesh_z_levels). */
    bool has_table () const { return !has_cell.empty(); }

    /**
     * The cells k that heights from zl to zh may lie in, whatever the column: [kb, kt], clamped to
     * [klo - 1, khi + 1] (so a range outside the domain stays outside it); [klo, khi] without bounds.
     * With the table: exactly the boxes' cells k whose lowest bottom node is not above zh and whose
     * highest top node is not below zl, or kb > kt when there is none. zl or zh may be infinite (not NaN):
     * heights from -infinity take every cell up to zh's, and a range at one infinity takes no cell of
     * the domain.
     */
    void k_range (amrex::Real zl, amrex::Real zh, int klo, int khi, int& kb, int& kt) const
    {
        kb = klo;
        kt = khi;
        if (has_table()) {
            kb = khi + 1;
            kt = klo - 1;
            // a relative 1e-5 more each way: the caller's zl and zh are rounded sums (a point's height and
            // its reach), and a face at a node level's own height h just past them still has weight
            // From the finite ends only: an infinite end needs none. An infinite end is replaced by 0 before
            // the difference is taken, not skipped by a test, so that inf - inf is never computed (a compiler
            // may evaluate both sides of a test, and the invalid operation would trap)
            const amrex::Real fl = std::isinf(zl) ? amrex::Real(0.0) : zl;
            const amrex::Real fh = std::isinf(zh) ? amrex::Real(0.0) : zh;
            const amrex::Real slack = amrex::Real(1.0e-5) * std::max(std::abs(fl), std::max(std::abs(fh), std::abs(fh - fl)));
            const amrex::Real lo = zl - slack, hi = zh + slack;
            for (std::size_t n = 0; n < has_cell.size(); ++n) {
                if (has_cell[n] == 0 || node_lo[n] > hi || node_hi[n+1] < lo) { continue; }
                const int k = k0 + static_cast<int>(n);
                kb = std::min(kb, k);
                kt = std::max(kt, k);
            }
            return;
        }
        if (!has_ground()) { return; }
        // counted from the ground's cell klo. A cell that only touches the range at a node is not always
        // within [kb, kt]: the spreader takes one face level more each way (reach_box_k())
        const double b = klo + std::floor((static_cast<double>(zl) - ground_hi) / dz_hi);
        const double t = klo + std::floor((static_cast<double>(zh) - ground_lo) / dz_lo);
        kb = static_cast<int>(std::min<double>(khi + 1, std::max<double>(klo - 1, b)));
        kt = static_cast<int>(std::min<double>(khi + 1, std::max<double>(klo - 1, t)));
    }
};

/**
 * The bounds of a terrain-following or stretched level's heights, from its z_phys_nd (nodal, one
 * component): a height Z lies in a cell k with floor((Z - ground_hi) / dz_hi) <= k <= floor((Z -
 * ground_lo) / dz_lo). Over every rank's boxes, the same on every rank (a collective), or with
 * over_all_ranks false over this rank's boxes only, which bounds the cells of those boxes alone.
 * With over_all_ranks false and a box of this rank that does not reach the ground (a column split
 * in z, a level aloft), where that formula does not hold, also the table of the boxes' node heights
 * (ZBounds::has_table). The table is this rank's: points_covered_by takes the global bounds only.
 */
ZBounds mesh_z_bounds (const amrex::MultiFab& z_phys_nd, const amrex::Geometry& geom, bool over_all_ranks = true);

/**
 * mesh_z_bounds over every rank's boxes, with the table of the level's node heights over every
 * rank's boxes too (ZBounds::has_table): each node level's lowest and highest height and the cells
 * the level holds, the same on every rank (a collective). ZBounds::k_range then gives the cells of
 * the level a height range may lie in whether or not the boxes reach the ground. For
 * region_touches_grids; points_covered_by takes the bounds without a table.
 */
ZBounds mesh_z_levels (const amrex::MultiFab& z_phys_nd, const amrex::Geometry& geom);

/**
 * mesh_z_levels of one level, kept between calls: the heights of a level change only when its boxes do
 * (a regrid) or when the terrain moves, so a caller that needs the table at every step builds it once.
 */
class LevelHeights {
public:
    /**
     * The table of z_phys_nd's level. Built at the first call, and again when z_phys_nd's boxes or geom's
     * domain differ from those of the last build, or with rebuild true (for heights that change in time on
     * the same boxes: a moving terrain). A collective when it builds; every rank builds at the same calls,
     * as the boxes, the domain and rebuild are the same on every rank.
     */
    const ZBounds& get (const amrex::MultiFab& z_phys_nd, const amrex::Geometry& geom, bool rebuild = false);
    /** How many times get() has built the table. */
    int builds () const { return m_builds; }
private:
    ZBounds m_zb;
    amrex::BoxArray m_ba;   //!< z_phys_nd's boxes at the last build
    amrex::Box m_domain;    //!< geom's domain at the last build
    int m_builds{0};
};

/**
 * The cells of a terrain-following or stretched level that the heights within reach of each of a set of
 * points may lie in, over the point's own columns (column_cell_ranges): point p needs the cells kb[p] to
 * kt[p]; kb > kt when its reach touches no node of the level.
 */
struct ColumnCells {
    std::vector<int> kb, kt;
    /** The ranges of n points from point `begin` (as a copy). */
    ColumnCells slice (std::size_t begin, std::size_t n) const
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(begin + n <= kb.size() && kb.size() == kt.size(), "ColumnCells::slice: past the ranges");
        ColumnCells c;
        c.kb.assign(kb.begin() + static_cast<std::ptrdiff_t>(begin), kb.begin() + static_cast<std::ptrdiff_t>(begin + n));
        c.kt.assign(kt.begin() + static_cast<std::ptrdiff_t>(begin), kt.begin() + static_cast<std::ptrdiff_t>(begin + n));
        return c;
    }
};

/**
 * The cells a point's heights within reach may lie in, over the point's own columns: for point p
 * (pos[3p .. 3p+2], reach[p] across and in height, in metres), the lowest and highest cell k such that,
 * in some node column the reach across touches (a periodic x or y wrapping across the seam, as
 * points_covered_by takes it), the cell's spacing [z(k), z(k+1)] overlaps [z - reach, z + reach] (to a
 * relative 1e-5, as ZBounds::k_range), from every rank's boxes of z_phys_nd. Where those heights go below
 * a box's lowest node in a column (the box's bottom face, unless it is the domain's floor), the cell under
 * that box counts too, and where they go above a box's highest node (unless it is the domain's top), the
 * cell above it: in a column split in z between boxes those cells are held by the next box, and where the
 * level stops they are not, so a kernel the level cuts fails points_covered_by, whatever the boxes. The
 * range is the hull [min k, max k] over the columns: points_covered_by never accepts a kernel some cell
 * of which the level lacks, and accepts every kernel whose cells the level holds when, over the columns
 * within reach, those cells form one run and no column has the heights below its ground under a box
 * that starts above the floor. Elsewhere it is conservative: it also asks for the cells of a gap between
 * two runs, or the cell under such a box (a cliff next to the point on a level aloft). The coverage
 * answer is the same whatever the boxes; the hull itself may differ by cells the level holds in the same
 * run. A reach partly below the floor or above the domain top needs nothing there, as with the mesh's
 * bounds; one wholly beyond them touches no node and is not covered. kb > kt for a point whose reach
 * touches no node of the level. A collective, the same on every rank: every rank passes the same points.
 *
 * @param[in] pos 3*n positions (m)
 * @param[in] reach n reaches (m)
 * @param[out] cells the ranges of the n points
 */
void column_cell_ranges (const amrex::MultiFab& z_phys_nd, const amrex::Geometry& geom,
                         const std::vector<amrex::Real>& pos, const std::vector<amrex::Real>& reach, ColumnCells& cells);

bool points_covered_by (const amrex::BoxArray& ba, const amrex::Geometry& geom,
                        const std::vector<amrex::Real>& pos, amrex::Real reach, std::string& first_outside,
                        CoverZ z = CoverZ::Reach, const ZBounds* zb = nullptr, const ColumnCells* cells = nullptr);

/**
 * A level's boxes as points_covered_by tests them, built once for every test of a check: the boxes
 * cell-centred, and flattened onto the domain's bottom layer (for CoverZ::Footprint).
 */
struct CoverageGrids
{
    CoverageGrids (const amrex::BoxArray& ba, const amrex::Geometry& geom);
    amrex::BoxArray cells;        //!< the boxes, cell-centred
    amrex::BoxArray footprints;   //!< each box flattened onto the bottom layer (they may overlap)
};

/** points_covered_by on boxes built once (CoverageGrids of the same level and geometry). */
bool points_covered_by (const CoverageGrids& grids, const amrex::Geometry& geom,
                        const std::vector<amrex::Real>& pos, amrex::Real reach, std::string& first_outside,
                        CoverZ z = CoverZ::Reach, const ZBounds* zb = nullptr, const ColumnCells* cells = nullptr);

/**
 * Whether a level's boxes may hold a cell of a region (lo and hi in metres, lo <= hi in every
 * direction, no NaN): the opposite question to points_covered_by, for a caller that must stay clear
 * of a level. x and y use the cells the region's ends fall in (points_covered_by's own); in a
 * periodic direction every image of the region counts, however many periods away, and a region a
 * period or more long takes every cell. z must not be periodic. A region that overlaps a cell of
 * the level is never missed; one that only meets a cell at a face may not count. In z:
 *  - on a uniform-dz mesh (zb null) the cells are the region's own, from prob_lo and the spacing;
 *  - with the level's table of node heights (mesh_z_levels), the cells k of the level some bottom
 *    node of which is not above the region and some top node not below it (ZBounds::k_range): right
 *    wherever the level's boxes are, on the ground or aloft, and no cell where the level has none
 *    at those heights;
 *  - with the level's bounds alone (mesh_z_bounds over all ranks) and every box on the ground,
 *    the cells between the ground and spacing extremes allow, which is wider on steep terrain;
 *  - with the bounds alone and a box that does not stand on the ground (ZBounds::all_ground false:
 *    a level or a box aloft, or a column split in z between boxes), the whole column: the answer is
 *    then yes wherever a box lies over or under the region's footprint.
 * The cell indices are bounded before they become integers, so a region that reaches far outside
 * the domain (or to infinity) is safe. The answer reads no field data and is the same on every
 * rank when zb is: pass the table or the bounds over all ranks, not one rank's (which this function
 * cannot tell apart).
 *
 * @param[in] ba the level's boxes (any index type)
 * @param[in] geom the level's geometry
 * @param[in] region the region (m)
 * @param[in] zb the level's table or bounds on a stretched or terrain-following mesh, or nullptr
 */
bool region_touches_grids (const amrex::BoxArray& ba, const amrex::Geometry& geom, const amrex::RealBox& region,
                           const ZBounds* zb = nullptr);

} // namespace erf_actuator

#endif
