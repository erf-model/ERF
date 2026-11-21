#ifndef ERF_PLANAR_BOUNDARY_H_
#define ERF_PLANAR_BOUNDARY_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Periodicity.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <memory>

/**
 * Ghost-cell fill for planar (2D) MultiFabs built on the collapse of a 3D BoxArray
 * along a surface-normal direction.
 *
 * The surface-layer fields and the MOST averages are defined on one 2D box per 3D box.
 * When the 3D BoxArray is split in the normal direction, several 3D boxes collapse onto
 * the same 2D box, so the planar BoxArray holds duplicate boxes with overlapping valid
 * regions, and only the copy that belongs to the 3D box touching the surface is ever
 * computed.  A FillBoundary on such a MultiFab may fill a ghost cell from an
 * uncomputed copy.
 *
 * define() records which planar boxes are those surface copies.  fill() copies them into
 * a buffer without duplicates and ParallelCopies the buffer onto every planar box, valid
 * region and ghost cells alike.  Unlike FillBoundary it therefore also overwrites the
 * valid region of the uncomputed copies, so that any copy can be read afterwards.  Without
 * duplicates fill() is a plain FillBoundary.
 *
 * The buffer is kept for each target layout, index type, and number of components, so only
 * the first fill of each kind allocates.  For a face-centered MultiFab (the MOST velocity averages)
 * neighboring surface boxes share a face in the buffer.  Both should hold the same value
 * there, but nothing enforces it, so fill() gives the shared faces a single value before
 * it copies: the lowest global box index wins, which does not depend on the decomposition
 * or the rank count, and every copy of the field then agrees.
 *
 * A level none of whose boxes reaches the surface (a fine level away from the surface) has
 * no computed copy to fill from, and fill() leaves the MultiFab as it is.
 */
class PlanarBoundary
{
public:
    /**
     * Record the surface copies of a planar BoxArray.
     *
     * @param[in] ba3d          3D BoxArray the planar BoxArray was collapsed from
     * @param[in] ba2d          planar BoxArray, one box per box of ba3d and in the same order
     * @param[in] dm            DistributionMapping shared by ba3d and ba2d
     * @param[in] surface_index index of the surface cell in the normal direction
     * @param[in] is_low        whether the surface is the low side of the domain
     * @param[in] normal_dir    normal direction of the surface
     */
    void define (const amrex::BoxArray& ba3d,
                 const amrex::BoxArray& ba2d,
                 const amrex::DistributionMapping& dm,
                 int surface_index,
                 bool is_low = true,
                 int normal_dir = 2);

    /**
     * Fill every copy of a planar MultiFab, valid region and ghost cells, from the
     * surface copies.
     *
     * @param[in,out] mf     planar MultiFab on ba2d, of any index type and number of components
     * @param[in]     period periodicity of the level
     */
    void fill (amrex::MultiFab& mf, const amrex::Periodicity& period);

    /**
     * Give every copy of a ghost cell outside a non-periodic lateral domain edge the value
     * of the surface copy that holds the column it faces.
     *
     * Where a box seam meets a non-periodic edge, the ghost cells beyond the edge next to
     * the seam lie in the grown boxes of both neighbors, and no FillBoundary reconciles
     * them: each box fills its copy from its own valid region.  A cell beyond the edge is
     * here assigned to the surface copy whose valid region holds the cell's projection into
     * the domain (its index clamped into the domain in each non-periodic direction), and that
     * copy's value is written into every other copy, periodic images included.  The result
     * then does not depend on the decomposition, and a checkpoint, which keeps one copy,
     * restarts from the same values.  Only ghost cells beyond a non-periodic edge are
     * written; interior, periodic and coarse/fine ghost cells are left as they are.
     *
     * On a fine level the column a cell beyond the edge faces need not be covered by any
     * grid: where a coarse/fine boundary meets the domain edge, or two grids along the
     * edge are separated by a gap narrower than twice the ghost width.  Such a cell has no
     * owner by projection, and it takes instead the value of the surface copy with the
     * lowest global box index whose ghost region holds it (OverrideSync's precedence in
     * fill()), so it too has one value that does not depend on the decomposition.  Only
     * surface copies take part: a cell held solely by the ghost region of an uncomputed
     * copy (a 3D box above the surface whose footprint differs from those below it) is
     * left as it is.
     *
     * Call it after fill(), on a cell-centered planar MultiFab; other index types are left
     * as they are.
     *
     * @param[in,out] mf   planar MultiFab on the planar BoxArray given to define
     * @param[in]     geom geometry of the level
     */
    void reconcile_outside_domain (amrex::MultiFab& mf, const amrex::Geometry& geom);

    /**
     * Copy the computed surface copies of a planar MultiFab into a MultiFab without
     * duplicates.
     *
     * define() gives each surface copy the rank of the planar box it is a copy of, so every
     * box of dst is filled from a FAB of mf that lives on the same rank and the gather needs
     * no communication.  Reading the surface copies directly, rather than any copy of the
     * field, means the caller need not have filled the duplicates with fill().
     *
     * @param[in]  mf    planar MultiFab on the planar BoxArray given to define
     * @param[out] dst   MultiFab on surface_boxes() converted to the index type of mf, and on
     *                   surface_dm()
     * @param[in]  scomp first component to read from mf
     * @param[in]  dcomp first component to write in dst
     * @param[in]  ncomp number of components
     */
    void gather_surface (const amrex::MultiFab& mf, amrex::MultiFab& dst,
                         int scomp, int dcomp, int ncomp) const;

    //! Cell-centered surface copies of the planar boxes, without duplicates
    [[nodiscard]] const amrex::BoxArray& surface_boxes () const { return m_ba_sfc; }

    //! Ranks owning the surface copies: each is on the rank of the planar box it is a copy of
    [[nodiscard]] const amrex::DistributionMapping& surface_dm () const { return m_dm_sfc; }

    //! Index of each surface copy in the planar BoxArray
    [[nodiscard]] const amrex::Vector<int>& surface_index () const { return m_src_index; }

    //! Whether a planar-box index is the computed copy that touches the surface
    [[nodiscard]] AMREX_GPU_HOST
    bool is_surface_copy (const int planar_index) const
    {
        return std::binary_search(m_src_index.begin(), m_src_index.end(), planar_index);
    }

private:
    amrex::MultiFab& buffer (const amrex::MultiFab& mf);

    struct Buffer
    {
        amrex::IndexType ixtype;
        int ncomp;
        amrex::BoxArray ba;
        amrex::DistributionMapping dm;
        amrex::Vector<const amrex::MultiFab*> sources;
        std::unique_ptr<amrex::MultiFab> mf;
    };

    // The strips beyond the non-periodic domain edges that each surface copy owns, for
    // one planar layout (BoxArray and DistributionMapping) and number of lateral ghost cells
    struct EdgeStrips
    {
        amrex::BoxArray planar_ba;
        amrex::DistributionMapping planar_dm;
        amrex::IntVect ngrow;
        amrex::BoxArray ba;              // the strips, disjoint, inside the domain in periodic directions
        amrex::DistributionMapping dm;   // each on the rank of its owner
        amrex::Vector<int> owner;        // planar index of the surface copy owning each strip
        amrex::Vector<amrex::IntVect> shift; // owner's cells are at the strip's shifted by this
        amrex::Vector<std::unique_ptr<amrex::MultiFab>> bufs; // copy buffers, one per ncomp
    };

    EdgeStrips& edge_strips (const amrex::MultiFab& mf, const amrex::Geometry& geom);

    int m_nplanar{0};                    // Number of boxes in the planar BoxArray
    amrex::BoxArray m_ba_sfc;            // Surface copies
    amrex::DistributionMapping m_dm_sfc; // Ranks owning them
    amrex::Vector<int> m_src_index;      // Their index in the planar BoxArray
    amrex::Vector<Buffer> m_buffers;     // Gather buffers
    amrex::Vector<EdgeStrips> m_strips;  // Edge strips, per planar layout and ghost width
};

/**
 * The ground plane of a field on a level's 3D grids, on one planar box per 3D box.
 *
 * Box ib of the result is the footprint of box ib of src, cut down to the plane k = klo in
 * the index type of src (so a nodal-in-z field keeps the nodes klo and klo+1), on the rank
 * of box ib.  It therefore lines up, box for box, with any MultiFab on the collapsed
 * BoxArray (ba2d) and the level's DistributionMapping, and a kernel that reads src at klo
 * can read it instead from any box, including one that starts above the ground.
 *
 * Each planar box first takes the lowest plane of its own 3D box, moved down to klo; the
 * ground plane is then copied in from the boxes that hold it.  Where no box of the level
 * reaches the ground (a refined level aloft), the planar box keeps that lowest plane: a
 * finite stand-in, not the ground, which callers must not report as the ground.
 *
 * @param[in] src field on the level's 3D grids (any index type, any number of components)
 * @param[in] klo index of the lowest cell of the domain
 */
std::unique_ptr<amrex::MultiFab>
copy_ground_plane (const amrex::MultiFab& src, int klo);

/**
 * Finish a column sum computed box by box on a planar MultiFab.
 *
 * A column sum over a level's 3D grids, accumulated into the planar box of each 3D box,
 * holds in each planar box only the part of the column inside its own 3D box.  When the
 * 3D BoxArray is split in z, several planar boxes cover one column; this adds their parts
 * and gives every one of them the whole sum.  A planar BoxArray without overlapping boxes
 * is left as it is.
 *
 * @param[in,out] mf    planar MultiFab, one box per 3D box, cell-centered, no ghost cells used
 * @param[in]     comp  first component holding a column sum
 * @param[in]     ncomp number of components
 */
void
sum_over_stacked_boxes (amrex::MultiFab& mf, int comp, int ncomp);

#endif
