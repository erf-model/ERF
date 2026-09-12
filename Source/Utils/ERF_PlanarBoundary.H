#ifndef ERF_PLANAR_BOUNDARY_H_
#define ERF_PLANAR_BOUNDARY_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_IndexType.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Periodicity.H>
#include <AMReX_Vector.H>

#include <memory>

/**
 * Ghost-cell fill for planar (2D) MultiFabs built on the z-collapse of a 3D BoxArray.
 *
 * The surface-layer fields and the MOST averages are defined on one 2D box per 3D box.
 * When the 3D BoxArray is split in z, several 3D boxes collapse onto the same 2D box, so
 * the planar BoxArray holds duplicate boxes with overlapping valid regions, and only the
 * copy that belongs to the 3D box touching the surface is ever computed.  A FillBoundary
 * on such a MultiFab may fill a ghost cell from an uncomputed copy.
 *
 * define() records which planar boxes are those surface copies.  fill() copies them into
 * a buffer without duplicates and ParallelCopies the buffer onto every planar box, valid
 * region and ghost cells alike.  Unlike FillBoundary it therefore also overwrites the
 * valid region of the uncomputed copies, so that any copy can be read afterwards.  Without
 * duplicates fill() is a plain FillBoundary.
 *
 * The buffer is kept for each index type and number of components, so only the first
 * fill of each kind allocates.  For a face-centered MultiFab (the MOST velocity averages)
 * neighbouring surface boxes share a face in the buffer; which of the two supplies that
 * face is up to ParallelCopy, as it was up to FillBoundary before, and both should hold
 * the same value.
 *
 * A level none of whose boxes reaches the surface (a fine level aloft) has no computed
 * copy to fill from, and fill() leaves the MultiFab as it is.
 */
class PlanarBoundary
{
public:
    /**
     * Record the surface copies of a planar BoxArray.
     *
     * @param[in] ba3d 3D BoxArray the planar BoxArray was collapsed from
     * @param[in] ba2d planar BoxArray, one box per box of ba3d and in the same order
     * @param[in] dm   DistributionMapping shared by ba3d and ba2d
     * @param[in] klo  k index of the lowest cell in the domain
     */
    void define (const amrex::BoxArray& ba3d,
                 const amrex::BoxArray& ba2d,
                 const amrex::DistributionMapping& dm,
                 int klo);

    /**
     * Fill every copy of a planar MultiFab, valid region and ghost cells, from the
     * surface copies.
     *
     * @param[in,out] mf     planar MultiFab on ba2d, of any index type and number of components
     * @param[in]     period periodicity of the level
     */
    void fill (amrex::MultiFab& mf, const amrex::Periodicity& period);

    //! Cell-centered surface copies of the planar boxes, without duplicates
    [[nodiscard]] const amrex::BoxArray& surface_boxes () const { return m_ba_sfc; }

    //! Index of each surface copy in the planar BoxArray
    [[nodiscard]] const amrex::Vector<int>& surface_index () const { return m_src_index; }

private:
    amrex::MultiFab& buffer (amrex::IndexType ixtype, int ncomp);

    struct Buffer
    {
        amrex::IndexType ixtype;
        int ncomp;
        std::unique_ptr<amrex::MultiFab> mf;
    };

    int m_nplanar{0};                    // Number of boxes in the planar BoxArray
    amrex::BoxArray m_ba_sfc;            // Surface copies
    amrex::DistributionMapping m_dm_sfc; // Ranks owning them
    amrex::Vector<int> m_src_index;      // Their index in the planar BoxArray
    amrex::Vector<Buffer> m_buffers;     // Gather buffers, one per index type and ncomp
};

#endif
