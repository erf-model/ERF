/**
 * \file ERF_ImmersedBuildingMasks.H
 *
 * How the buildings immersed forcing reads a momentum face and which forcing it applies there,
 * shared by the source kernels (ImmersedForcingBuildings_{X,Y,Z}mom) and the implicit projection
 * (ERF::make_if_projection_sigma, erf.if_implicit_projection, erf-model/ERF#4887), so the two
 * agree face by face on the blanking, on which faces carry the MOST wall law and on which the
 * linear drag toward zero.
 */
#ifndef ERF_IMMERSED_BUILDING_MASKS_H_
#define ERF_IMMERSED_BUILDING_MASKS_H_

#include <AMReX_Array4.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

#include "ERF_NumericalConstants.H"

namespace if_bld {

/** Cell blanking snapped to solid (1) or fluid (0) at half (erf.if_snap_partial_cells), else as is */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real snapped (amrex::Real b, bool snap) noexcept
{
    return snap ? ((b >= myhalf) ? amrex::Real(one) : amrex::Real(zero)) : b;
}

/**
 * The face blanking the explicit path zeroes on after the projection (ERF_SlowRhsPost.cpp): the
 * face-centred fraction when the build has it, else the mean of the two cells the face joins.
 * (io, jo, ko) is the unit vector of the face's direction.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real raw_face_blank (const amrex::Array4<const amrex::Real>& t, const amrex::Array4<const amrex::Real>& tf,
                            int i, int j, int k, int io, int jo, int ko) noexcept
{
    return (tf) ? tf(i,j,k) : myhalf * (t(i,j,k) + t(i-io,j-jo,k-ko));
}

/**
 * The face blanking as the building kernels read it: with the snap, from the snapped cells (solid
 * when either is), so both EB builds give the same staircase; otherwise raw_face_blank; then
 * thresholded at small (eb2.small_volfrac), below which the face counts as fluid.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real face_blank (const amrex::Array4<const amrex::Real>& t, const amrex::Array4<const amrex::Real>& tf,
                        int i, int j, int k, int io, int jo, int ko, bool snap, amrex::Real small) noexcept
{
    const amrex::Real b = snap ? amrex::max(snapped(t(i,j,k), snap), snapped(t(i-io,j-jo,k-ko), snap))
                               : raw_face_blank(t, tf, i, j, k, io, jo, ko);
    return (b < small) ? amrex::Real(zero) : b;
}

/**
 * A horizontal-momentum face (x or y) and its neighbours as the kernels read them. lo/hi are the
 * neighbouring faces across the other horizontal direction: south/north for an x-face, west/east
 * for a y-face. raw is the unthresholded blanking of the face (raw_face_blank); normal is true
 * when, with the snap, the face joins a solid and a fluid cell (wall-normal).
 */
struct HorizFace {
    amrex::Real t, below, above, lo, hi, raw;
    bool normal;
};

/**
 * @param t     Cell blanking (with at least one ghost cell in each direction)
 * @param tf    Face-centred blanking of this direction, or an empty Array4
 * @param klo   Lowest cell index of the domain (the face below the domain bottom counts as fluid)
 * @param snap  erf.if_snap_partial_cells
 * @param small eb2.small_volfrac
 */
template <int dir>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
HorizFace horiz_face (const amrex::Array4<const amrex::Real>& t, const amrex::Array4<const amrex::Real>& tf,
                      int i, int j, int k, int klo, bool snap, amrex::Real small) noexcept
{
    static_assert(dir == 0 || dir == 1, "horiz_face: x- or y-faces only");
    constexpr int io = (dir == 0) ? 1 : 0;
    constexpr int jo = (dir == 1) ? 1 : 0;
    HorizFace f;
    f.t      = face_blank(t, tf, i, j, k, io, jo, 0, snap, small);
    f.below  = (k == klo) ? amrex::Real(zero) : face_blank(t, tf, i, j, k-1, io, jo, 0, snap, small);
    f.above  = face_blank(t, tf, i, j, k+1, io, jo, 0, snap, small);
    f.lo     = (dir == 0) ? face_blank(t, tf, i, j-1, k, io, jo, 0, snap, small)
                          : face_blank(t, tf, i-1, j, k, io, jo, 0, snap, small);
    f.hi     = (dir == 0) ? face_blank(t, tf, i, j+1, k, io, jo, 0, snap, small)
                          : face_blank(t, tf, i+1, j, k, io, jo, 0, snap, small);
    f.raw    = raw_face_blank(t, tf, i, j, k, io, jo, 0);
    f.normal = snap && (snapped(t(i,j,k), snap) != snapped(t(i-io,j-jo,k), snap));
    return f;
}

/** A z-face and its neighbours as the kernels read them (see HorizFace) */
struct VertFace {
    amrex::Real t, above, south, north, west, east, raw;
    bool normal;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
VertFace vert_face (const amrex::Array4<const amrex::Real>& t, const amrex::Array4<const amrex::Real>& tf,
                    int i, int j, int k, bool snap, amrex::Real small) noexcept
{
    VertFace f;
    f.t      = face_blank(t, tf, i, j, k,   0, 0, 1, snap, small);
    f.above  = face_blank(t, tf, i, j, k+1, 0, 0, 1, snap, small);
    f.south  = face_blank(t, tf, i, j-1, k, 0, 0, 1, snap, small);
    f.north  = face_blank(t, tf, i, j+1, k, 0, 0, 1, snap, small);
    f.west   = face_blank(t, tf, i-1, j, k, 0, 0, 1, snap, small);
    f.east   = face_blank(t, tf, i+1, j, k, 0, 0, 1, snap, small);
    f.raw    = raw_face_blank(t, tf, i, j, k, 0, 0, 1);
    f.normal = snap && (snapped(t(i,j,k), snap) != snapped(t(i,j,k-1), snap));
    return f;
}

/**
 * Masks of a horizontal-momentum face (x or y). lo/hi are the neighbouring faces across the
 * other horizontal direction: south/north for an x-face, west/east for a y-face.
 */
struct HorizMasks {
    amrex::Real roof, lo, hi, wall, side, interior;
    //! faces carrying the MOST wall law (a target velocity from the log law)
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real most () const noexcept { return roof + lo + hi; }
    //! weight of the linear drag CdM |u| toward zero on this face, with t its blanking
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real linear (amrex::Real t) const noexcept { return (wall + side) * t + interior; }
};

/**
 * @param t           Blanking of the face
 * @param below       Blanking of the face below (zero at the domain bottom)
 * @param above       Blanking of the face above
 * @param lo          Blanking of the neighbouring face on the low side across the other direction
 * @param hi          Blanking of the neighbouring face on the high side across the other direction
 * @param normal_face With the snap, a face joining a solid and a fluid cell (wall-normal)
 * @param use_most    erf.if_use_most
 * @param snap        erf.if_snap_partial_cells
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
HorizMasks horiz_masks (amrex::Real t, amrex::Real below, amrex::Real above,
                        amrex::Real lo, amrex::Real hi,
                        bool normal_face, bool use_most, bool snap) noexcept
{
    HorizMasks m{};
    // With the snap every solid face has t = 1: a roof face lies in the top solid row
    // (t <= below, the face above fluid), and a face carrying a wall law is not also an
    // interior face. A wall-normal face gets the interior drag (no penetration), not a wall law.
    m.roof = (!normal_face && t > zero && (snap ? t <= below : t < below) && above == zero && use_most) ? one : zero;
    m.lo   = (!normal_face && t > zero && t <= hi && lo == zero && use_most) ? one : zero;
    m.hi   = (!normal_face && t > zero && t <= lo && hi == zero && use_most) ? one : zero;
    m.wall = (t > zero && t < one && !use_most) ? one : zero;                       // all walls without MOST
    m.side = (t > zero && t < one && use_most && m.most() == zero) ? one : zero;    // partial faces not under MOST
    m.interior = (t == one && !(snap && m.most() > zero)) ? one : zero;
    return m;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
HorizMasks horiz_masks (const HorizFace& f, bool use_most, bool snap) noexcept
{
    return horiz_masks(f.t, f.below, f.above, f.lo, f.hi, f.normal, use_most, snap);
}

/**
 * Masks of a z-face. The roof face takes the linear drag (weighted by its blanking), not a
 * wall law; the four sides carry the wall law of the vertical walls.
 */
struct VertMasks {
    amrex::Real south, north, west, east, wall, roof, interior;
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real most () const noexcept { return south + north + west + east; }
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real linear (amrex::Real t) const noexcept { return (wall + roof) * t + interior; }
};

/**
 * @param t           Blanking of the face
 * @param above       Blanking of the face above
 * @param south,north,west,east Blanking of the neighbouring z-faces
 * @param above_floor The face is not the domain bottom (k >= 1)
 * @param normal_face With the snap, a face joining a solid and a fluid cell (wall-normal)
 * @param use_most    erf.if_use_most
 * @param snap        erf.if_snap_partial_cells
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
VertMasks vert_masks (amrex::Real t, amrex::Real above,
                      amrex::Real south, amrex::Real north, amrex::Real west, amrex::Real east,
                      bool above_floor, bool normal_face, bool use_most, bool snap) noexcept
{
    VertMasks m{};
    m.south = (!normal_face && t > zero && t <= north && south == zero && use_most && above_floor) ? one : zero;
    m.north = (!normal_face && t > zero && t <= south && north == zero && use_most && above_floor) ? one : zero;
    m.west  = (!normal_face && t > zero && t <= east  && west  == zero && use_most && above_floor) ? one : zero;
    m.east  = (!normal_face && t > zero && t <= west  && east  == zero && use_most && above_floor) ? one : zero;
    m.wall  = (t > zero && t < one && !use_most) ? one : zero;
    m.roof  = (!normal_face && t > zero && above == zero && use_most) ? one : zero;
    // with the snap a face carrying a wall law or the roof drag is not also an interior face
    m.interior = (t == one && !(snap && (m.most() + m.roof) > zero)) ? one : zero;
    return m;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
VertMasks vert_masks (const VertFace& f, bool above_floor, bool use_most, bool snap) noexcept
{
    return vert_masks(f.t, f.above, f.south, f.north, f.west, f.east, above_floor, f.normal, use_most, snap);
}

/**
 * Weight of the drag the implicit projection applies on a face (erf.if_implicit_projection):
 * the linear-drag weight of the kernel, and the full rate on every face the explicit path zeroes
 * after the projection (raw blanking 1, which the kernels leave without a wall law). The latter
 * also caps the vertical roof face of an exact box without the snap, where the explicit kernel
 * adds the roof and the interior drag together (weight 2 t), at one.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real implicit_weight (amrex::Real linear, amrex::Real raw) noexcept
{
    return (raw == one) ? amrex::Real(one) : linear;
}

} // namespace if_bld

#endif
