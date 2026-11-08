#ifndef ERF_READ_CHECKPOINT_MF_H_
#define ERF_READ_CHECKPOINT_MF_H_

#include <string>

#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Periodicity.H>
#include <AMReX_VisMF.H>

namespace erf_checkpoint {

/**
 * The periodicity of geom in x and y only. For planar fields, and for the terrain heights,
 * which are not periodic in z even when the domain is: in a domain periodic in z, the top
 * node of a column is the image of the bottom one but its height is not.
 */
inline amrex::Periodicity
lateral_periodicity (const amrex::Geometry& geom)
{
    const amrex::Box& dom = geom.Domain();
    return amrex::Periodicity(amrex::IntVect(AMREX_D_DECL(geom.isPeriodic(0) ? dom.length(0) : 0,
                                                          geom.isPeriodic(1) ? dom.length(1) : 0,
                                                          0)));
}

/**
 * Read one checkpoint MultiFab and deposit it into dst, which may live on different grids.
 *
 * The source is read into an undefined MultiFab, which VisMF::Read then defines on the
 * BoxArray recorded in the file, so no caller needs to know the checkpoint's grids -- which
 * is what lets a restart land a component's state on level-0 grids that moved (issue 4225).
 *
 * Only the VALID cells of the source are offered by default: a file's ghost cells are not
 * guaranteed to agree with its valid data, and ParallelCopy intersects the grown source with
 * the grown destination without any ordering guarantee between overlapping sources, so
 * offering them can land stale values in valid cells. Pass src_ghosts = true for a field
 * whose halo the checkpoint carries and whose ordinary restart has always restored it; the
 * valid cells are then laid over the top in a second pass, so they win wherever they reach.
 *
 * @param[inout] dst         destination, on the grids the run uses
 * @param[in]    mf_name     the checkpoint MultiFab, as for VisMF::Read
 * @param[in]    scomp       first component of the source to read
 * @param[in]    dcomp       first component of dst to fill
 * @param[in]    ncomp       number of components
 * @param[in]    dst_ghosts  ghost cells of dst to fill (capped at what dst holds)
 * @param[in]    src_ghosts  offer the file's ghost cells as well as its valid cells
 * @param[in]    period      periodicity of the domain dst lives in
 */
inline void
read_checkpoint_multifab (amrex::MultiFab& dst,
                          const std::string& mf_name,
                          int scomp, int dcomp, int ncomp,
                          const amrex::IntVect& dst_ghosts,
                          bool src_ghosts,
                          const amrex::Periodicity& period)
{
    //
    // Undefined on purpose: VisMF::Read defines it on the BoxArray the file records, which is
    // the decomposition the checkpoint was written with and need not be this run's.
    //
    amrex::MultiFab src;
    amrex::VisMF::Read(src, mf_name);

    //
    // A file written with fewer ghost cells than the destination holds cannot fill more than
    // it has, and the destination cannot take more than it holds.
    //
    const amrex::IntVect ng_dst = amrex::min(dst_ghosts, dst.nGrowVect());

    //
    // The overwhelmingly common case: this level's grids ARE the grids the checkpoint was
    // written on, so nothing has to move between boxes. Copy box to box, exactly as the
    // reader did before the grids could differ, so that an ordinary restart stays
    // byte-for-byte what it always was -- no periodic images, and no overlapping sources to
    // be ambiguous about.
    //
    // The distribution still has to be matched: VisMF::Read gives the MultiFab it defines a
    // DistributionMapping of its own, which is not this run's even when the BoxArray is.
    // Redistribute is the box-to-box move for exactly that case.
    //
    if (src.boxArray() == dst.boxArray())
    {
        const amrex::IntVect ng_copy = src_ghosts ? amrex::min(ng_dst, src.nGrowVect())
                                                  : amrex::IntVect(0);

        if (src.DistributionMap() == dst.DistributionMap()) {
            amrex::MultiFab::Copy(dst, src, scomp, dcomp, ncomp, ng_copy);
        } else {
            amrex::MultiFab aligned(src.boxArray(), dst.DistributionMap(), src.nComp(), src.nGrowVect());
            aligned.Redistribute(src, 0, 0, src.nComp(), src.nGrowVect());
            amrex::MultiFab::Copy(dst, aligned, scomp, dcomp, ncomp, ng_copy);
        }
        return;
    }

    if (src_ghosts) {
        //
        // Two passes. The first offers the file's halo so that destination ghost cells OUTSIDE
        // the domain -- which no valid cell can reach -- get the values the checkpoint recorded
        // for them. But a file's halo is not guaranteed to agree with its valid data, and
        // offering it alone would let a stale ghost cell land in a valid cell through a box
        // overlap or a periodic image, with no ordering guarantee between overlapping sources.
        // So the second pass copies the valid cells on top, and they win wherever they reach.
        //
        const amrex::IntVect ng_src = amrex::min(ng_dst, src.nGrowVect());
        dst.ParallelCopy(src, scomp, dcomp, ncomp, ng_src,            ng_dst, period);
        dst.ParallelCopy(src, scomp, dcomp, ncomp, amrex::IntVect(0), ng_dst, period);
    } else {
        dst.ParallelCopy(src, scomp, dcomp, ncomp, amrex::IntVect(0), ng_dst, period);
    }
}

/**
 * As read_checkpoint_multifab, for an integer field the checkpoint stores as a MultiFab,
 * valid cells only. dst is (re)defined on ba/dm, which are the grids the run uses.
 */
inline void
read_checkpoint_imultifab (amrex::iMultiFab& dst,
                           const amrex::BoxArray& ba,
                           const amrex::DistributionMapping& dm,
                           const std::string& mf_name,
                           const amrex::Periodicity& period)
{
    amrex::MultiFab tmp(ba, dm, 1, 0);
    read_checkpoint_multifab(tmp, mf_name, 0, 0, 1, amrex::IntVect(0), false, period);
    dst = amrex::cast<amrex::iMultiFab>(tmp);
}

} // namespace erf_checkpoint

#endif
