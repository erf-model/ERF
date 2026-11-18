// Gaussian force spreading with discrete normalisation, and the source integral.

#include "ERF_ActuatorSpreading.H"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_OpenMP.H>
#include <AMReX_ParallelDescriptor.H>

#include "ERF_ActuatorDeviceBuffers.H"
#include "ERF_ActuatorGeometry.H"
#include "ERF_ActuatorSampling.H"

using namespace amrex;

namespace erf_actuator {

namespace {

// The kernel weight exp(-r^2/eps^2) of face (i,j,k) of grid dir for the point (px,py,pz) (one
// periodic image) of width eps, 0 beyond the cut-off radius sqrt(reach2) = 3 eps
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real kernel_weight (int dir, int i, int j, int k, Real px, Real py, Real pz,
                    const GpuArray<Real,AMREX_SPACEDIM>& plo,
                    const GpuArray<Real,AMREX_SPACEDIM>& dx,
                    Array4<Real const> const& znd, bool has_znd,
                    Real inv_eps2, Real reach2)
{
    const Real ddx = plo[0] + (i + ((dir == 0) ? Real(0.0) : Real(0.5))) * dx[0] - px;
    const Real ddy = plo[1] + (j + ((dir == 1) ? Real(0.0) : Real(0.5))) * dx[1] - py;
    const Real ddz = face_height(dir, i, j, k, znd, has_znd, plo[2], dx[2]) - pz;
    const Real r2 = ddx*ddx + ddy*ddy + ddz*ddz;
    return (r2 > reach2) ? Real(0.0) : std::exp(-r2 * inv_eps2);
}

// The faces of grid dir that a box owns: its cells' low faces, so a face shared with the next
// box (or, at a periodic end, with the first box's image) is counted once. On a non-periodic
// boundary the low domain face carries no source; the high one is never owned.
Box owned_faces (int dir, const Box& valid_cells)
{
    Box fbx = amrex::convert(valid_cells, IntVect::TheDimensionVector(dir));
    fbx.growHi(dir, -1);
    return fbx;
}

// whether face index f of grid dir lies on a non-periodic domain boundary
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool boundary_face (int dir, int i, int j, int k, const GpuArray<int,AMREX_SPACEDIM>& per,
                    const GpuArray<int,AMREX_SPACEDIM>& dlo, const GpuArray<int,AMREX_SPACEDIM>& dhi)
{
    if (per[dir] != 0) { return false; }
    const int f = (dir == 0) ? i : ((dir == 1) ? j : k);
    return (f == dlo[dir]) || (f == dhi[dir] + 1);
}

} // namespace

void
spread_forces (const std::vector<Real>& pos, const std::vector<Real>& force, Real epsilon,
               const MultiFab* z_phys_nd, const MultiFab* detJ_cc, const Geometry& geom,
               MultiFab& src_x, MultiFab& src_y, MultiFab& src_z)
{
    // also without points, where the per-point check of the widths has nothing to refuse
    if (!(epsilon > 0.0) || !std::isfinite(epsilon)) {
        Abort("spread_forces: the kernel width " + std::to_string(epsilon) + " m must be positive and finite");
    }
    spread_forces(pos, force, std::vector<Real>(pos.size() / 3, epsilon), z_phys_nd, detJ_cc, geom, src_x, src_y, src_z);
}

std::string
spread_points_problem (const std::vector<Real>& pos, const std::vector<Real>& force, const std::vector<Real>& epsilon)
{
    if (!(pos.size() == force.size() && pos.size() % 3 == 0)) { return "spread_forces: one x,y,z force per x,y,z point"; }
    if (3 * epsilon.size() != pos.size()) { return "spread_forces: one kernel width per point"; }
    for (std::size_t p = 0; p < epsilon.size(); ++p) {
        if (!(epsilon[p] > 0.0) || !std::isfinite(epsilon[p])) {
            return "spread_forces: point " + std::to_string(p) + " (0-based) has the kernel width " + std::to_string(epsilon[p]) +
                   " m; it must be positive and finite";
        }
    }
    for (std::size_t q = 0; q < pos.size(); ++q) {
        if (!std::isfinite(pos[q]) || !std::isfinite(force[q])) {
            const std::size_t p = q / 3;
            return "spread_forces: point " + std::to_string(p) + " (0-based) has a non-finite position or force: (" +
                   std::to_string(pos[3*p]) + ", " + std::to_string(pos[3*p+1]) + ", " + std::to_string(pos[3*p+2]) + ") m, (" +
                   std::to_string(force[3*p]) + ", " + std::to_string(force[3*p+1]) + ", " + std::to_string(force[3*p+2]) + ") N";
        }
    }
    return {};
}

void
check_spread_points (const std::vector<Real>& pos, const std::vector<Real>& force, const std::vector<Real>& epsilon)
{
    const std::string problem = spread_points_problem(pos, force, epsilon);
    if (!problem.empty()) { Abort(problem); }
}

void
check_spread_layout (const MultiFab* z_phys_nd, const MultiFab* detJ_cc, const MultiFab& src_x, const MultiFab& src_y,
                     const MultiFab& src_z)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(src_x.ixType() == IndexType(IntVect(1,0,0)) && src_y.ixType() == IndexType(IntVect(0,1,0)) &&
                                     src_z.ixType() == IndexType(IntVect(0,0,1)),
                                     "spread_forces: src_x, src_y and src_z must be face-centred in x, y and z");
    for (const MultiFab* mf : {&src_y, &src_z, z_phys_nd, detJ_cc}) {
        if (mf == nullptr) { continue; }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mf->DistributionMap() == src_x.DistributionMap() && mf->boxArray().CellEqual(src_x.boxArray()),
                                         "spread_forces: the sources, z_phys_nd and detJ_cc must share boxes and distribution");
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(detJ_cc == nullptr || detJ_cc->nGrow() >= 1,
                                     "spread_forces: detJ_cc needs one ghost cell (the face volumes read the cell below each low face)");
}

namespace {

// What the two passes over the faces share: the level's spacing, its periodicity, the bounds of a
// fitted mesh's heights over this rank's boxes, and the faces within reach of a point.
struct KernelFaces
{
    GpuArray<Real,AMREX_SPACEDIM> plo{}, dx{}, period{{0.0, 0.0, 0.0}};
    GpuArray<int,AMREX_SPACEDIM> per{{0, 0, 0}}, dlo{{0, 0, 0}}, dhi{{0, 0, 0}};
    Real dxdydz{0.0};
    bool has_znd{false}, has_detj{false};
    ZBounds zb;

    KernelFaces (const Geometry& geom, const MultiFab* z_phys_nd, const MultiFab* detJ_cc, const ZBounds* zb_local)
    {
        plo = geom.ProbLoArray();
        dx = geom.CellSizeArray();
        dxdydz = dx[0] * dx[1] * dx[2];
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            per[d] = geom.isPeriodic(d) ? 1 : 0;
            period[d] = geom.isPeriodic(d) ? geom.ProbLength(d) : Real(0.0);
            dlo[d] = geom.Domain().smallEnd(d);
            dhi[d] = geom.Domain().bigEnd(d);
        }
        has_znd = (z_phys_nd != nullptr);
        has_detj = (detJ_cc != nullptr);
        // over this rank's boxes only: the faces this rank visits are in them, so their bounds hold. With the
        // table of its node heights they hold for any box (one aloft, a column split in z); without it the
        // bounds from the boxes that reach the ground do not hold for a box aloft over higher ground
        if (has_znd) { zb = (zb_local != nullptr) ? *zb_local : mesh_z_bounds(*z_phys_nd, geom, false); }
        if (!zb.all_ground && !zb.has_table()) { zb = ZBounds{}; }
    }

    // The index boxes of the faces of the three face grids within reach of a point at (px,py,pz); with
    // z_phys_nd (terrain-following or stretched) a face's height is not its index times dz, so a box takes
    // every k a height within reach may lie in, from the mesh's bounds (with the table of the rank's node
    // heights where its boxes do not all reach the ground), and kernel_weight's cut-off at 3 eps, measured
    // in physical heights, picks the faces. The cells in z are found once for the three grids
    std::array<Box,3> reach_boxes (Real px, Real py, Real pz, Real eps) const
    {
        int kb = dlo[2], kt = dhi[2];
        if (has_znd) { zb.k_range(pz - Real(3.0) * eps, pz + Real(3.0) * eps, dlo[2], dhi[2], kb, kt); }
        return {{reach_box_k(0, px, py, pz, eps, kb, kt), reach_box_k(1, px, py, pz, eps, kb, kt),
                 reach_box_k(2, px, py, pz, eps, kb, kt)}};
    }
    // The box of grid dir, with the cells [kb, kt] a height within reach may lie in (used with z_phys_nd only)
    Box reach_box_k (int dir, Real px, Real py, Real pz, Real eps, int kb, int kt) const
    {
        const Real pc[3] = {px, py, pz};
        IntVect lo, hi;
        for (int d = 0; d < 3; ++d) {
            const int reach_cells = static_cast<int>(std::ceil(Real(3.0) * eps / dx[d])) + 1;
            const Real fi = (pc[d] - plo[d]) / dx[d] - ((d == dir) ? Real(0.0) : Real(0.5));
            lo[d] = static_cast<int>(std::floor(fi)) - reach_cells;
            hi[d] = static_cast<int>(std::floor(fi)) + reach_cells + 1;
        }
        if (has_znd) {
            if (kb > kt) {
                // none of this rank's cells holds a height within reach: no face of its boxes is
                lo[2] = dhi[2] + 2;
                hi[2] = dhi[2] + 1;
                return Box(lo, hi, IntVect::TheDimensionVector(dir));
            }
            // the faces bounding those cells, one more each way for the faces of a cell's neighbours
            lo[2] = std::max(kb - 1, dlo[2]);
            hi[2] = std::min(kt + 2, dhi[2] + 1);
        }
        return Box(lo, hi, IntVect::TheDimensionVector(dir));
    }

    // The periodic images of a point that can reach the domain, appended to images: shifts of
    // -n L .. +n L per periodic direction, with n the number of domain widths its kernel reaches
    // (one unless the domain is narrower than 3 epsilon_p, as in a unit test)
    void point_images (Real px, Real py, Real pz, Real eps, std::vector<std::array<Real,3>>& images) const
    {
        int n_img[2] = {0, 0};
        for (int d = 0; d < 2; ++d) {
            if (per[d] != 0) { n_img[d] = static_cast<int>(std::ceil(Real(3.0) * eps / period[d])); }
        }
        for (int sx = -n_img[0]; sx <= n_img[0]; ++sx) {
            for (int sy = -n_img[1]; sy <= n_img[1]; ++sy) {
                images.push_back({{px + sx * period[0], py + sy * period[1], pz}});
            }
        }
    }
};

// The side of the face tiles the second pass bins the kernels to (faces per direction)
constexpr int tile_size = 8;

// The periodic images of every point, made once per call (not per box and face grid) and sent to the
// device in one copy each, in point order and, for a point, in point_images' order, which is the order
// in which the passes add their terms: the image's position, its point, and the faces of each grid
// within its reach (reach_boxes). Also each point's 1 / eps^2 and squared cut-off radius.
struct KernelImages
{
    std::vector<std::array<Real,3>> pos;
    std::vector<int> point;
    std::vector<std::array<Box,3>> reach;

    KernelImages (const KernelFaces& kf, const std::vector<Real>& ppos, const std::vector<Real>& epsilon,
                  ActuatorDeviceBuffers& buf)
    {
        const int npts = static_cast<int>(epsilon.size());
        for (int p = 0; p < npts; ++p) {
            const std::size_t first = pos.size();
            kf.point_images(ppos[3*p], ppos[3*p+1], ppos[3*p+2], epsilon[p], pos);
            point.resize(pos.size(), p);
            for (std::size_t m = first; m < pos.size(); ++m) {
                reach.push_back(kf.reach_boxes(pos[m][0], pos[m][1], pos[m][2], epsilon[p]));
            }
        }
        buf.image_pos.host.clear();
        buf.image_reach.host.clear();
        for (std::size_t m = 0; m < pos.size(); ++m) {
            buf.image_pos.host.insert(buf.image_pos.host.end(), pos[m].begin(), pos[m].end());
            for (int dir = 0; dir < 3; ++dir) {
                const Box& rb = reach[m][dir];
                buf.image_reach.host.insert(buf.image_reach.host.end(), {rb.smallEnd(0), rb.smallEnd(1), rb.smallEnd(2),
                                                                         rb.bigEnd(0), rb.bigEnd(1), rb.bigEnd(2)});
            }
        }
        buf.image_point.host.assign(point.begin(), point.end());
        buf.point_eps.host.resize(2 * static_cast<std::size_t>(npts));
        for (int p = 0; p < npts; ++p) {
            const Real eps = epsilon[p];
            buf.point_eps.host[2*p]   = Real(1.0) / (eps * eps);
            buf.point_eps.host[2*p+1] = Real(9.0) * eps * eps;
        }
        buf.image_pos.upload();
        buf.image_point.upload();
        buf.image_reach.upload();
        buf.point_eps.upload();
    }
};

// Where a (box, face grid)'s bin starts in the bins buffer, and how many segments (points) or tiles it has
struct Bin
{
    int start{0};
    int count{0};
};

// The first pass's bin of one box and face grid, appended to bins: the points one of whose images
// reaches the faces fbx, each with the images that do. Laid out as the n points, then n + 1 offsets
// (into bins) of their images' lists, then the lists.
Bin norm_bin (const KernelImages& im, int dir, const Box& fbx, std::vector<int>& bins)
{
    std::vector<int> pts, offs, items;
    const int nimg = static_cast<int>(im.point.size());
    for (int m = 0; m < nimg; ++m) {
        if ((im.reach[m][dir] & fbx).isEmpty()) { continue; }
        if (pts.empty() || pts.back() != im.point[m]) {
            pts.push_back(im.point[m]);
            offs.push_back(static_cast<int>(items.size()));
        }
        items.push_back(m);
    }
    Bin b;
    b.start = static_cast<int>(bins.size());
    b.count = static_cast<int>(pts.size());
    if (b.count == 0) { return b; }
    offs.push_back(static_cast<int>(items.size()));
    const int item0 = b.start + 2 * b.count + 1;
    bins.insert(bins.end(), pts.begin(), pts.end());
    for (int o : offs) { bins.push_back(item0 + o); }
    bins.insert(bins.end(), items.begin(), items.end());
    return b;
}

// The second pass's bin of one box and face grid, appended to bins: the faces fbx cut into tiles of
// tile_size^3 faces from its low corner, and for each tile that some image reaches with a nonzero scale,
// the images that reach it, in image order. Laid out as the n tiles' low corners (3 each), then n + 1
// offsets (into bins) of their lists, then the lists.
Bin add_bin (const KernelImages& im, const std::vector<Real>& f_over_s, int dir, const Box& fbx, std::vector<int>& bins)
{
    const IntVect flo = fbx.smallEnd();
    const IntVect nt = (fbx.length() + IntVect(tile_size - 1)) / tile_size;
    const int ntiles = nt[0] * nt[1] * nt[2];
    std::vector<int> cnt(ntiles + 1, 0);
    std::vector<std::array<IntVect,2>> touched;   // per image: its tile range (empty when it reaches no tile)
    const int nimg = static_cast<int>(im.point.size());
    touched.resize(nimg, {{IntVect(0), IntVect(-1)}});
    for (int m = 0; m < nimg; ++m) {
        if (f_over_s[3*im.point[m]+dir] == Real(0.0)) { continue; }
        const Box c = im.reach[m][dir] & fbx;
        if (c.isEmpty()) { continue; }
        const IntVect t_lo = (c.smallEnd() - flo) / tile_size;
        const IntVect t_hi = (c.bigEnd() - flo) / tile_size;
        touched[m] = {{t_lo, t_hi}};
        for (int tz = t_lo[2]; tz <= t_hi[2]; ++tz) {
            for (int ty = t_lo[1]; ty <= t_hi[1]; ++ty) {
                for (int tx = t_lo[0]; tx <= t_hi[0]; ++tx) { ++cnt[1 + tx + nt[0] * (ty + nt[1] * tz)]; }
            }
        }
    }
    std::vector<int> used;    // the tiles some image reaches
    for (int t = 0; t < ntiles; ++t) { if (cnt[t+1] > 0) { used.push_back(t); } }
    Bin b;
    b.start = static_cast<int>(bins.size());
    b.count = static_cast<int>(used.size());
    if (b.count == 0) { return b; }
    // each used tile's list starts after the corners and offsets
    std::vector<int> next(ntiles, -1);
    int item = b.start + 4 * b.count + 1;
    for (int t : used) {
        const int tx = t % nt[0], ty = (t / nt[0]) % nt[1], tz = t / (nt[0] * nt[1]);
        bins.insert(bins.end(), {flo[0] + tile_size * tx, flo[1] + tile_size * ty, flo[2] + tile_size * tz});
    }
    for (int t : used) {
        bins.push_back(item);
        next[t] = item;
        item += cnt[t+1];
    }
    bins.push_back(item);
    bins.resize(item, -1);
    for (int m = 0; m < nimg; ++m) {
        const IntVect& t_lo = touched[m][0];
        const IntVect& t_hi = touched[m][1];
        for (int tz = t_lo[2]; tz <= t_hi[2]; ++tz) {
            for (int ty = t_lo[1]; ty <= t_hi[1]; ++ty) {
                for (int tx = t_lo[0]; tx <= t_hi[0]; ++tx) { bins[next[tx + nt[0] * (ty + nt[1] * tz)]++] = m; }
            }
        }
    }
    return b;
}

} // namespace

void
spread_partial_norms (const std::vector<Real>& pos, const std::vector<Real>& epsilon,
                      const MultiFab* z_phys_nd, const MultiFab* detJ_cc, const Geometry& geom,
                      const MultiFab& src_x, std::vector<Real>& norm, const ZBounds* zb_local,
                      ActuatorDeviceBuffers* buffers)
{
    BL_PROFILE("erf_actuator::spread_partial_norms()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(buffers == nullptr || !OpenMP::in_parallel(),
                                     "spread_partial_norms: kept buffers must not be shared by the threads of an OpenMP parallel region");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pos.size() == 3 * epsilon.size(), "spread_partial_norms: one kernel width per x,y,z point");
    const int npts = static_cast<int>(epsilon.size());
    norm.assign(3 * static_cast<std::size_t>(npts), Real(0.0));
    if (npts == 0) { return; }
    ActuatorDeviceBuffers own;
    ActuatorDeviceBuffers& buf = (buffers != nullptr) ? *buffers : own;
    KernelFaces kf(geom, z_phys_nd, detJ_cc, zb_local);
    const auto plo = kf.plo;
    const auto dx = kf.dx;
    const auto per = kf.per;
    const auto dlo = kf.dlo;
    const auto dhi = kf.dhi;
    const Real dxdydz = kf.dxdydz;
    const bool has_znd = kf.has_znd;
    const bool has_detj = kf.has_detj;

    // Work per point over the faces within its reach, so the cost is the sum over the points of
    // (6 eps_p / dx)^3 ((6 eps_p / dx)^2 x nz with z_phys_nd) rather than nfaces x npts. A point near a periodic
    // boundary is visited again as its periodic image(s), which is how the kernel wraps.
    const KernelImages im(kf, pos, epsilon, buf);
    buf.norms.host.assign(3 * static_cast<std::size_t>(npts), Real(0.0));
    buf.norms.upload();

    // the bins of every box and face grid, sent in one copy
    std::vector<std::array<Bin,3>> bin;
    buf.bins.host.clear();
    for (MFIter mfi(src_x, false); mfi.isValid(); ++mfi) {
        const Box vbx = amrex::enclosedCells(mfi.validbox());
        bin.push_back({{norm_bin(im, 0, owned_faces(0, vbx), buf.bins.host),
                        norm_bin(im, 1, owned_faces(1, vbx), buf.bins.host),
                        norm_bin(im, 2, owned_faces(2, vbx), buf.bins.host)}});
    }
    buf.bins.upload();
    const int* p_bins = buf.bins.data();
    const Real* p_q = buf.image_pos.data();
    const int* p_reach = buf.image_reach.data();
    const Real* p_eps = buf.point_eps.data();
    Real* p_norm = buf.norms.data();

    // S_p = sum over the faces of w_p dV, each physical face once (the owned faces of this rank's
    // boxes, boundary faces excluded). One launch per box and face grid, one thread per point (a
    // segment of the images): a thread adds its point's terms in the order of the images and, for an
    // image, of the faces (k, then j, then i), with no atomics, so that the sum does not depend on
    // the device's scheduling. A point whose kernel reaches two boxes is updated by both launches,
    // which read and write its slot of p_norm without atomics: the boxes therefore share one device
    // stream, so that their launches run one after the other, in box order
    int ib = 0;
    for (MFIter mfi(src_x, MFItInfo().SetNumStreams(1)); mfi.isValid(); ++mfi, ++ib) {
        const Box vbx = amrex::enclosedCells(mfi.validbox());
        Array4<Real const> znd = has_znd ? z_phys_nd->const_array(mfi) : Array4<Real const>{};
        Array4<Real const> detj = has_detj ? detJ_cc->const_array(mfi) : Array4<Real const>{};
        for (int dir = 0; dir < 3; ++dir) {
            const Bin b = bin[ib][dir];
            if (b.count == 0) { continue; }
            const Box fbx = owned_faces(dir, vbx);
            const GpuArray<int,3> flo{{fbx.smallEnd(0), fbx.smallEnd(1), fbx.smallEnd(2)}};
            const GpuArray<int,3> fhi{{fbx.bigEnd(0), fbx.bigEnd(1), fbx.bigEnd(2)}};
            const int* seg_pt = p_bins + b.start;
            const int* seg_off = p_bins + b.start + b.count;
            ++buf.kernel_launches;
            ParallelFor(b.count, [=] AMREX_GPU_DEVICE (int s) noexcept
            {
                const int p = seg_pt[s];
                const Real inv_eps2 = p_eps[2*p], reach2 = p_eps[2*p+1];
                Real acc = p_norm[3*p+dir];
                for (int u = seg_off[s]; u < seg_off[s+1]; ++u) {
                    const int m = p_bins[u];
                    const Real qx = p_q[3*m], qy = p_q[3*m+1], qz = p_q[3*m+2];
                    const int* r = p_reach + 18 * m + 6 * dir;
                    const int ilo = amrex::max(r[0], flo[0]), jlo = amrex::max(r[1], flo[1]), klo = amrex::max(r[2], flo[2]);
                    const int ihi = amrex::min(r[3], fhi[0]), jhi = amrex::min(r[4], fhi[1]), khi = amrex::min(r[5], fhi[2]);
                    for (int k = klo; k <= khi; ++k) {
                        for (int j = jlo; j <= jhi; ++j) {
                            for (int i = ilo; i <= ihi; ++i) {
                                if (boundary_face(dir, i, j, k, per, dlo, dhi)) { continue; }
                                const Real w = kernel_weight(dir, i, j, k, qx, qy, qz, plo, dx, znd, has_znd, inv_eps2, reach2);
                                // the product in a statement of its own: fused with the addition it would
                                // round differently on a target with a fused multiply-add
                                if (w > Real(0.0)) {
                                    const Real term = w * face_volume(dir, i, j, k, detj, has_detj, dxdydz);
                                    acc += term;
                                }
                            }
                        }
                    }
                }
                p_norm[3*p+dir] = acc;
            });
        }
    }
    buf.norms.download_async(norm.size());
    Gpu::streamSynchronize();
    std::copy(buf.norms.pinned.begin(), buf.norms.pinned.begin() + norm.size(), norm.begin());
}


void
add_scaled_kernels (const std::vector<Real>& pos, const std::vector<Real>& f_over_s, const std::vector<Real>& epsilon,
                    const MultiFab* z_phys_nd, const MultiFab* detJ_cc, const Geometry& geom,
                    MultiFab& src_x, MultiFab& src_y, MultiFab& src_z, const ZBounds* zb_local,
                    ActuatorDeviceBuffers* buffers)
{
    BL_PROFILE("erf_actuator::add_scaled_kernels()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(buffers == nullptr || !OpenMP::in_parallel(),
                                     "add_scaled_kernels: kept buffers must not be shared by the threads of an OpenMP parallel region");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pos.size() == 3 * epsilon.size() && f_over_s.size() == pos.size(),
                                     "add_scaled_kernels: one kernel width and one x,y,z scale per x,y,z point");
    const int npts = static_cast<int>(epsilon.size());
    if (npts == 0) { return; }
    ActuatorDeviceBuffers own;
    ActuatorDeviceBuffers& buf = (buffers != nullptr) ? *buffers : own;
    KernelFaces kf(geom, z_phys_nd, detJ_cc, zb_local);
    const auto plo = kf.plo;
    const auto dx = kf.dx;
    const auto per = kf.per;
    const auto dlo = kf.dlo;
    const auto dhi = kf.dhi;
    const bool has_znd = kf.has_znd;
    MultiFab* src[3] = {&src_x, &src_y, &src_z};

    const KernelImages im(kf, pos, epsilon, buf);
#ifndef AMREX_USE_GPU
    if (!buf.force_face_major && !OpenMP::in_parallel()) {
        // On the host (where ParallelFor runs its iterations in order, OpenMP or not, as this pass is not called
        // from a parallel region): each image adds into the faces of its own reach box, in image (point, then
        // image) order. A face then takes the same terms in the same order as from the face-major pass below,
        // whose threads (one per face, needed on a device, where images would race on a face) each walk every
        // image binned to their tile, half of which miss the face
        const std::vector<Real>& eps_host = buf.point_eps.host;
        for (MFIter mfi(src_x, false); mfi.isValid(); ++mfi) {
            Array4<Real const> znd = has_znd ? z_phys_nd->const_array(mfi) : Array4<Real const>{};
            for (int dir = 0; dir < 3; ++dir) {
                const Box fbx = mfi.validbox().convert(IntVect::TheDimensionVector(dir));
                Array4<Real> const& s = src[dir]->array(mfi);
                for (std::size_t m = 0; m < im.point.size(); ++m) {
                    const int p = im.point[m];
                    const Real scale = f_over_s[3*p+dir];
                    if (scale == Real(0.0)) { continue; }
                    const Box c = im.reach[m][dir] & fbx;
                    if (c.isEmpty()) { continue; }
                    const Real qx = im.pos[m][0], qy = im.pos[m][1], qz = im.pos[m][2];
                    const Real inv_eps2 = eps_host[2*p], reach2 = eps_host[2*p+1];
                    const auto lo = amrex::lbound(c);
                    const auto hi = amrex::ubound(c);
                    for (int k = lo.z; k <= hi.z; ++k) {
                        for (int j = lo.y; j <= hi.y; ++j) {
                            for (int i = lo.x; i <= hi.x; ++i) {
                                if (boundary_face(dir, i, j, k, per, dlo, dhi)) { continue; }
                                const Real w = kernel_weight(dir, i, j, k, qx, qy, qz, plo, dx, znd, has_znd, inv_eps2, reach2);
                                if (w > Real(0.0)) {
                                    const Real term = scale * w;   // a statement of its own, as in the first pass
                                    s(i,j,k) += term;
                                }
                            }
                        }
                    }
                }
            }
        }
        return;
    }
#endif
    buf.scales.host.assign(f_over_s.begin(), f_over_s.end());
    buf.scales.upload();
    std::vector<std::array<Bin,3>> bin;
    buf.bins.host.clear();
    for (MFIter mfi(src_x, false); mfi.isValid(); ++mfi) {
        std::array<Bin,3> bd;
        for (int dir = 0; dir < 3; ++dir) {
            bd[dir] = add_bin(im, f_over_s, dir, mfi.validbox().convert(IntVect::TheDimensionVector(dir)), buf.bins.host);
        }
        bin.push_back(bd);
    }
    buf.bins.upload();
    const int* p_bins = buf.bins.data();
    const Real* p_q = buf.image_pos.data();
    const int* p_pt = buf.image_point.data();
    const int* p_reach = buf.image_reach.data();
    const Real* p_eps = buf.point_eps.data();
    const Real* p_scale = buf.scales.data();

    // src(face) += (F_p / S_p) w_p on every face of every box, so that both copies of a shared or
    // periodic image face carry the same value. One launch per box and face grid over the faces of
    // the tiles some kernel reaches, one thread per face: a face adds the terms of the images binned
    // to its tile, in point then image order, with no atomics
    constexpr int t3 = tile_size * tile_size * tile_size;
    int ib = 0;
    for (MFIter mfi(src_x, false); mfi.isValid(); ++mfi, ++ib) {
        Array4<Real const> znd = has_znd ? z_phys_nd->const_array(mfi) : Array4<Real const>{};
        for (int dir = 0; dir < 3; ++dir) {
            const Bin b = bin[ib][dir];
            if (b.count == 0) { continue; }
            const Box fbx = mfi.validbox().convert(IntVect::TheDimensionVector(dir));
            const GpuArray<int,3> fhi{{fbx.bigEnd(0), fbx.bigEnd(1), fbx.bigEnd(2)}};
            Array4<Real> const& s = src[dir]->array(mfi);
            const int* corner = p_bins + b.start;
            const int* offs = p_bins + b.start + 3 * b.count;
            ++buf.kernel_launches;
            ParallelFor(b.count * t3, [=] AMREX_GPU_DEVICE (int idx) noexcept
            {
                const int t = idx / t3;
                const int l = idx - t * t3;
                const int i = corner[3*t]   + l % tile_size;
                const int j = corner[3*t+1] + (l / tile_size) % tile_size;
                const int k = corner[3*t+2] + l / (tile_size * tile_size);
                if (i > fhi[0] || j > fhi[1] || k > fhi[2]) { return; }
                if (boundary_face(dir, i, j, k, per, dlo, dhi)) { return; }
                Real acc = s(i,j,k);
                bool added = false;
                for (int u = offs[t]; u < offs[t+1]; ++u) {
                    const int m = p_bins[u];
                    const int* r = p_reach + 18 * m + 6 * dir;
                    if (i < r[0] || j < r[1] || k < r[2] || i > r[3] || j > r[4] || k > r[5]) { continue; }
                    const int p = p_pt[m];
                    const Real w = kernel_weight(dir, i, j, k, p_q[3*m], p_q[3*m+1], p_q[3*m+2], plo, dx, znd, has_znd,
                                                 p_eps[2*p], p_eps[2*p+1]);
                    if (w > Real(0.0)) {
                        const Real term = p_scale[3*p+dir] * w;   // a statement of its own, as in the first pass
                        acc += term;
                        added = true;
                    }
                }
                if (added) { s(i,j,k) = acc; }
            });
        }
    }
    Gpu::streamSynchronize();
}

void
add_spread_forces (const std::vector<Real>& pos, const std::vector<Real>& force, const std::vector<Real>& epsilon,
                   const MultiFab* z_phys_nd, const MultiFab* detJ_cc, const Geometry& geom,
                   MultiFab& src_x, MultiFab& src_y, MultiFab& src_z)
{
    check_spread_points(pos, force, epsilon);
    check_spread_layout(z_phys_nd, detJ_cc, src_x, src_y, src_z);
    const int npts = static_cast<int>(pos.size() / 3);
    check_same_points_everywhere(npts, "spread_forces");
    if (npts == 0) { return; }

    // pass 1: S_p, this rank's part, summed across ranks (this rank's bounds of a fitted mesh, once for both passes)
    ZBounds zb;
    if (z_phys_nd != nullptr) { zb = mesh_z_bounds(*z_phys_nd, geom, false); }
    ActuatorDeviceBuffers buf;
    std::vector<Real> norm;
    spread_partial_norms(pos, epsilon, z_phys_nd, detJ_cc, geom, src_x, norm, &zb, &buf);
    ParallelDescriptor::ReduceRealSum(norm.data(), 3 * npts);
    for (int p = 0; p < npts; ++p) {
        for (int dir = 0; dir < 3; ++dir) {
            if (!(norm[3*p+dir] > Real(0.0))) {
                Abort("spread_forces: point (" + std::to_string(pos[3*p]) + ", " + std::to_string(pos[3*p+1]) +
                      ", " + std::to_string(pos[3*p+2]) + ") reaches no face of the mesh within 3 epsilon = " +
                      std::to_string(3.0 * epsilon[p]) + " m; actuator points must lie inside the domain");
            }
        }
    }
    // pass 2: F_p w_p / S_p on the faces
    std::vector<Real> f_over_s(3 * static_cast<std::size_t>(npts));
    for (std::size_t q = 0; q < f_over_s.size(); ++q) { f_over_s[q] = force[q] / norm[q]; }
    add_scaled_kernels(pos, f_over_s, epsilon, z_phys_nd, detJ_cc, geom, src_x, src_y, src_z, &zb, &buf);
}

void
spread_forces (const std::vector<Real>& pos, const std::vector<Real>& force, const std::vector<Real>& epsilon,
               const MultiFab* z_phys_nd, const MultiFab* detJ_cc, const Geometry& geom,
               MultiFab& src_x, MultiFab& src_y, MultiFab& src_z)
{
    check_spread_points(pos, force, epsilon);
    check_spread_layout(z_phys_nd, detJ_cc, src_x, src_y, src_z);
    src_x.setVal(0.0);
    src_y.setVal(0.0);
    src_z.setVal(0.0);
    add_spread_forces(pos, force, epsilon, z_phys_nd, detJ_cc, geom, src_x, src_y, src_z);
}

Real
integrate_source (int dir, const MultiFab& src, const MultiFab* detJ_cc, const Geometry& geom)
{
    BL_PROFILE("erf_actuator::integrate_source()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dir >= 0 && dir <= 2 && src.ixType() == IndexType(IntVect::TheDimensionVector(dir)),
                                     "integrate_source: src must be face-centred in direction dir (0, 1 or 2)");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(detJ_cc == nullptr || detJ_cc->nGrow() >= 1, "integrate_source: detJ_cc needs one ghost cell");
    const auto dx = geom.CellSizeArray();
    const Real dxdydz = dx[0] * dx[1] * dx[2];
    const bool has_detj = (detJ_cc != nullptr);
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (MFIter mfi(src, false); mfi.isValid(); ++mfi) {
        // each physical face once: a box's cells' low faces
        const Box fbx = owned_faces(dir, amrex::enclosedCells(mfi.validbox()));
        Array4<Real const> const& s = src.const_array(mfi);
        Array4<Real const> detj = has_detj ? detJ_cc->const_array(mfi) : Array4<Real const>{};
        reduce_op.eval(fbx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            return { s(i,j,k) * face_volume(dir, i, j, k, detj, has_detj, dxdydz) };
        });
    }
    Real total = amrex::get<0>(reduce_data.value(reduce_op));
    ParallelDescriptor::ReduceRealSum(total);
    return total;
}

} // namespace erf_actuator
