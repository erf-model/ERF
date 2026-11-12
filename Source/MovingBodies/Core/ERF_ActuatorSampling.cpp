// Point sampling of face and cell fields and of the terrain height, and the coverage check.

#include "ERF_ActuatorSampling.H"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <array>
#include <limits>
#include <string>
#include <vector>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

#include "ERF_ActuatorGeometry.H"

using namespace amrex;

namespace erf_actuator {

namespace {

// The layout the samplers rely on: a field read through the MFIter of `ref` must share its
// distribution and cell boxes and have `ngrow` ghost cells (host-side check).
void check_layout (const MultiFab& ref, const MultiFab* other, int ngrow, const std::string& what)
{
    if (other == nullptr) { return; }
    if (!(other->nGrow() >= ngrow && other->DistributionMap() == ref.DistributionMap() && other->boxArray().CellEqual(ref.boxArray()))) {
        Abort(what + " must have " + std::to_string(ngrow) + " ghost cell(s) and the boxes and distribution of the field sampled");
    }
}

// Value of the face field f at height z in the column (i,j) of its staggered grid, linear in
// the physical height between the two faces that bracket z. The faces klo..khi are searched;
// outside them the nearest pair extrapolates, so a field linear in height stays exact down to
// the ground and up to the top without reading boundary ghost cells.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real column_value (Array4<Real const> const& f, int dir, int i, int j, int klo, int khi, Real z,
                   Array4<Real const> const& znd, bool has_znd, Real zlo, Real dz)
{
    int k0 = klo;
    for (int k = klo; k < khi; ++k) {
        if (z >= face_height(dir, i, j, k, znd, has_znd, zlo, dz)) { k0 = k; } else { break; }
    }
    const Real h0 = face_height(dir, i, j, k0,   znd, has_znd, zlo, dz);
    const Real h1 = face_height(dir, i, j, k0+1, znd, has_znd, zlo, dz);
    const Real w  = (z - h0) / (h1 - h0);
    return (Real(1.0) - w) * f(i,j,k0) + w * f(i,j,k0+1);
}

// One velocity component at (x,y,z): bilinear in the two horizontal directions between the
// four staggered columns around the point, each column evaluated at the physical height z.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real sample_component (int dir, Array4<Real const> const& f, Real x, Real y, Real z,
                       const GpuArray<Real,AMREX_SPACEDIM>& plo,
                       const GpuArray<Real,AMREX_SPACEDIM>& dxi,
                       int klo, int khi,
                       Array4<Real const> const& znd, bool has_znd, Real dz)
{
    // fractional index on the component's own grid: whole numbers on faces, half on centres
    const Real xi = (x - plo[0]) * dxi[0] - ((dir == 0) ? Real(0.0) : Real(0.5));
    const Real yi = (y - plo[1]) * dxi[1] - ((dir == 1) ? Real(0.0) : Real(0.5));
    const int i0 = static_cast<int>(std::floor(xi));
    const int j0 = static_cast<int>(std::floor(yi));
    const Real wx = xi - i0;
    const Real wy = yi - j0;
    const Real f00 = column_value(f, dir, i0,   j0,   klo, khi, z, znd, has_znd, plo[2], dz);
    const Real f10 = column_value(f, dir, i0+1, j0,   klo, khi, z, znd, has_znd, plo[2], dz);
    const Real f01 = column_value(f, dir, i0,   j0+1, klo, khi, z, znd, has_znd, plo[2], dz);
    const Real f11 = column_value(f, dir, i0+1, j0+1, klo, khi, z, znd, has_znd, plo[2], dz);
    return (Real(1.0) - wy) * ((Real(1.0) - wx) * f00 + wx * f10)
         +               wy  * ((Real(1.0) - wx) * f01 + wx * f11);
}

} // namespace

// The point list must be the same on every rank: the reductions below are collective, so a rank
// with another count (or none while the others have points) would hang or mis-sum. Checked
// with one collective, which every rank reaches whatever its count.
void check_same_points_everywhere (int npts, const char* who)
{
    int lo = npts, hi = npts;
    ParallelDescriptor::ReduceIntMin(lo);
    ParallelDescriptor::ReduceIntMax(hi);
    if (lo != hi) {
        Abort(std::string(who) + ": the ranks disagree on the number of points (" + std::to_string(lo) + " to " +
              std::to_string(hi) + "); every rank must pass the same point list");
    }
}

std::vector<Real>
wrap_periodic (const std::vector<Real>& pos, const Geometry& geom)
{
    std::vector<Real> p(pos);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (!geom.isPeriodic(d)) { continue; }
        const Real lo = static_cast<Real>(geom.ProbLo(d)), len = static_cast<Real>(geom.ProbLength(d));
        const Real hi = static_cast<Real>(geom.ProbHi(d));
        for (std::size_t i = static_cast<std::size_t>(d); i < p.size(); i += 3) {
            Real x = std::fmod(p[i] - lo, len);
            if (x < Real(0.0)) { x += len; }
            p[i] = lo + x;
            // a point a few ulps below lo can round onto hi or past it (with lo != 0), which no box owns
            if (p[i] >= hi) { p[i] = lo; }
        }
    }
    return p;
}

void
sample_velocity_local (const MultiFab& U, const MultiFab& V, const MultiFab& W,
                       const MultiFab* z_phys_nd, const Geometry& geom,
                       const std::vector<Real>& wpos, std::vector<Real>& vel, std::vector<int>& cnt,
                       ActuatorDeviceBuffers* buffers)
{
    BL_PROFILE("erf_actuator::sample_velocity_local()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(wpos.size() % 3 == 0, "sample_velocity: pos holds x,y,z triples");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(U.nGrow() >= 1 && V.nGrow() >= 1 && W.nGrow() >= 1,
                                     "sample_velocity needs one filled ghost cell on the velocities");
    check_layout(U, &V, 1, "sample_velocity: V");
    check_layout(U, &W, 1, "sample_velocity: W");
    check_layout(U, z_phys_nd, 1, "sample_velocity: z_phys_nd");
    const int npts = static_cast<int>(wpos.size() / 3);
    vel.assign(3 * static_cast<std::size_t>(npts), Real(0.0));
    cnt.assign(static_cast<std::size_t>(npts), 0);
    if (npts == 0) { return; }

    const Box& domain = geom.Domain();
    const auto plo = geom.ProbLoArray();
    const auto phi = geom.ProbHiArray();
    const auto dxi = geom.InvCellSizeArray();
    const Real dz = geom.CellSize(2);
    const int klo = domain.smallEnd(2);
    const int khi = domain.bigEnd(2);
    const int ihi = domain.bigEnd(0), jhi = domain.bigEnd(1);
    const bool has_znd = (z_phys_nd != nullptr);

    // the points, and zeros for the answers, to the device in one copy each (into buffers kept from call to call)
    ActuatorDeviceBuffers own;
    ActuatorDeviceBuffers& buf = (buffers != nullptr) ? *buffers : own;
    buf.sample_pos.host.assign(wpos.begin(), wpos.end());
    buf.sample_vel.host.assign(wpos.size(), Real(0.0));
    buf.sample_cnt.host.assign(static_cast<std::size_t>(npts), 0);
    buf.sample_pos.upload();
    buf.sample_vel.upload();
    buf.sample_cnt.upload();
    const Real* p_pos = buf.sample_pos.data();
    Real* p_vel = buf.sample_vel.data();
    int*  p_cnt = buf.sample_cnt.data();

    for (MFIter mfi(U, false); mfi.isValid(); ++mfi) {
        // the cells this rank owns; a point is sampled by the box holding its containing cell
        const Box vbx = amrex::enclosedCells(mfi.validbox());
        Array4<Real const> const& u = U.const_array(mfi);
        Array4<Real const> const& v = V.const_array(mfi);
        Array4<Real const> const& w = W.const_array(mfi);
        Array4<Real const> znd = has_znd ? z_phys_nd->const_array(mfi) : Array4<Real const>{};

        ++buf.kernel_launches;
        ParallelFor(npts, [=] AMREX_GPU_DEVICE (int p) noexcept
        {
            const Real x = p_pos[3*p], y = p_pos[3*p+1], z = p_pos[3*p+2];
            const Real xf = (x - plo[0]) * dxi[0], yf = (y - plo[1]) * dxi[1];
            int ic = static_cast<int>(std::floor(xf));
            int jc = static_cast<int>(std::floor(yf));
            // a point on the domain's upper face belongs to the last cell (decided in space: x/dx can
            // round one ulp above the cell count, as 1000 m on 88 cells does)
            if (ic > ihi && x <= phi[0]) { ic = ihi; }
            if (jc > jhi && y <= phi[1]) { jc = jhi; }
            // the columns are searched within this box and its first ghost layer (filled), so a box split in z,
            // or a fine level that stops below the domain's top, reads only its own data. A point the box owns
            // has its own column's bracketing pair in that range; a neighbouring column whose levels a slope
            // raises or lowers by more than about half a cell, at the box's top or bottom, extrapolates from
            // the range's end instead (a difference of the field's curvature, none for a field linear in height)
            const int kb = amrex::max(klo, vbx.smallEnd(2) - 1);
            const int kt = amrex::min(khi, vbx.bigEnd(2) + 1);
            if (ic < vbx.smallEnd(0) || ic > vbx.bigEnd(0) ||
                jc < vbx.smallEnd(1) || jc > vbx.bigEnd(1)) { return; }
            // the containing cell in z: between the z faces of column (ic,jc)
            int kc = -1;
            if (!has_znd) {
                int k = static_cast<int>(std::floor((z - plo[2]) * dxi[2]));
                // a point on the domain's top face belongs to the top cell (the cell, as on the upper x and y faces;
                // there the value is extrapolated from the top cells, as everywhere in z, not read from a ghost cell)
                if (k > khi && z <= phi[2]) { k = khi; }
                if (k >= klo && k <= kt) { kc = k; }
            } else {
                for (int k = kb; k <= kt; ++k) {
                    const Real zb = face_height(2, ic, jc, k,   znd, true, plo[2], dz);
                    const Real zt = face_height(2, ic, jc, k+1, znd, true, plo[2], dz);
                    // the top cell holds its top face too (the domain's top)
                    if (z >= zb && (z < zt || (k == khi && z <= zt))) { kc = k; break; }
                }
                // the bottom face is the mean of its four nodes, but the ground under (x, y) is
                // bilinear between them: on a slope a point just above the ground can lie below
                // the face, and belongs to the bottom cell (read by extrapolation) as long as it is
                // not below the lowest of the four nodes
                if (kc < 0 && kb == klo && z < face_height(2, ic, jc, klo, znd, true, plo[2], dz) &&
                    z >= amrex::min(amrex::min(znd(ic,jc,klo), znd(ic+1,jc,klo)), amrex::min(znd(ic,jc+1,klo), znd(ic+1,jc+1,klo)))) {
                    kc = klo;
                }
            }
            if (kc < vbx.smallEnd(2) || kc > vbx.bigEnd(2)) { return; }

            p_vel[3*p]   = sample_component(0, u, x, y, z, plo, dxi, kb, kt,   znd, has_znd, dz);
            p_vel[3*p+1] = sample_component(1, v, x, y, z, plo, dxi, kb, kt,   znd, has_znd, dz);
            p_vel[3*p+2] = sample_component(2, w, x, y, z, plo, dxi, kb, kt+1, znd, has_znd, dz);
            // a count, not a flag: boxes of one rank may run on different streams, and a point held by two of
            // them must still be refused by the caller
            Gpu::Atomic::AddNoRet(&p_cnt[p], 1);
        });
    }
    // the answers back with asynchronous copies and one synchronisation
    buf.sample_vel.download_async(vel.size());
    buf.sample_cnt.download_async(cnt.size());
    Gpu::streamSynchronize();
    std::copy(buf.sample_vel.pinned.begin(), buf.sample_vel.pinned.begin() + vel.size(), vel.begin());
    std::copy(buf.sample_cnt.pinned.begin(), buf.sample_cnt.pinned.begin() + cnt.size(), cnt.begin());
}

void
sample_velocity (const MultiFab& U, const MultiFab& V, const MultiFab& W,
                 const MultiFab* z_phys_nd, const Geometry& geom,
                 const std::vector<Real>& pos, std::vector<Real>& vel)
{
    BL_PROFILE("erf_actuator::sample_velocity()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pos.size() % 3 == 0, "sample_velocity: pos holds x,y,z triples");
    const int npts = static_cast<int>(pos.size() / 3);
    check_same_points_everywhere(npts, "sample_velocity");
    std::vector<int> cnt;
    sample_velocity_local(U, V, W, z_phys_nd, geom, wrap_periodic(pos, geom), vel, cnt);
    if (npts == 0) { return; }
    ParallelDescriptor::ReduceRealSum(vel.data(), static_cast<int>(vel.size()));
    ParallelDescriptor::ReduceIntSum(cnt.data(), npts);

    for (int p = 0; p < npts; ++p) {
        if (cnt[p] != 1) {
            Abort("sample_velocity: point (" + std::to_string(pos[3*p]) + ", " + std::to_string(pos[3*p+1]) +
                  ", " + std::to_string(pos[3*p+2]) + ") was sampled by " + std::to_string(cnt[p]) +
                  " boxes; actuator points must lie inside the domain");
        }
    }
}

namespace {

// height of the centre of cell (i,j,k): the mean of its bottom and top face heights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real centre_height (int i, int j, int k, Array4<Real const> const& znd, bool has_znd, Real zlo, Real dz)
{
    return Real(0.5) * (face_height(2, i, j, k, znd, has_znd, zlo, dz) + face_height(2, i, j, k+1, znd, has_znd, zlo, dz));
}

// the cell-centred field in column (i,j) at height z: linear between the two centres that
// bracket z, the nearest pair extrapolating beyond the first and last centre
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real centre_column_value (Array4<Real const> const& f, int comp, int i, int j, int klo, int khi, Real z,
                          Array4<Real const> const& znd, bool has_znd, Real zlo, Real dz)
{
    int k0 = klo;
    for (int k = klo; k < khi; ++k) {
        if (z >= centre_height(i, j, k, znd, has_znd, zlo, dz)) { k0 = k; } else { break; }
    }
    const Real h0 = centre_height(i, j, k0,   znd, has_znd, zlo, dz);
    const Real h1 = centre_height(i, j, k0+1, znd, has_znd, zlo, dz);
    const Real w  = (z - h0) / (h1 - h0);
    return (Real(1.0) - w) * f(i,j,k0,comp) + w * f(i,j,k0+1,comp);
}

} // namespace

void
sample_cell_scalar (const MultiFab& mf, int comp, const MultiFab* z_phys_nd, const Geometry& geom,
                    const std::vector<Real>& pos, std::vector<Real>& val)
{
    BL_PROFILE("erf_actuator::sample_cell_scalar()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pos.size() % 3 == 0, "sample_cell_scalar: pos holds x,y,z triples");
    const int npts = static_cast<int>(pos.size() / 3);
    val.assign(npts, 0.0);
    check_same_points_everywhere(npts, "sample_cell_scalar");
    if (npts == 0) { return; }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mf.nGrow() >= 1, "sample_cell_scalar: the field needs a filled ghost cell");
    check_layout(mf, z_phys_nd, 1, "sample_cell_scalar: z_phys_nd");
    const auto plo = geom.ProbLoArray();
    const auto phi = geom.ProbHiArray();
    const auto dxi = geom.InvCellSizeArray();
    const Real dz = geom.CellSize(2);
    const Box& domain = geom.Domain();
    const int klo = domain.smallEnd(2), khi = domain.bigEnd(2);
    const int ihi = domain.bigEnd(0), jhi = domain.bigEnd(1);
    const bool has_znd = (z_phys_nd != nullptr);
    const std::vector<Real> wpos = wrap_periodic(pos, geom);

    Gpu::DeviceVector<Real> d_pos(pos.size()), d_val(npts, 0.0);
    Gpu::DeviceVector<int> d_cnt(npts, 0);
    Gpu::copy(Gpu::hostToDevice, wpos.begin(), wpos.end(), d_pos.begin());
    Real* p_pos = d_pos.data();
    Real* p_val = d_val.data();
    int* p_cnt = d_cnt.data();

    for (MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const Box vbx = mfi.validbox();
        Array4<Real const> const& f = mf.const_array(mfi);
        Array4<Real const> znd = has_znd ? z_phys_nd->const_array(mfi) : Array4<Real const>{};
        ParallelFor(npts, [=] AMREX_GPU_DEVICE (int p) noexcept
        {
            const Real x = p_pos[3*p], y = p_pos[3*p+1], z = p_pos[3*p+2];
            const Real xf = (x - plo[0]) * dxi[0], yf = (y - plo[1]) * dxi[1];
            int ic = static_cast<int>(std::floor(xf));
            int jc = static_cast<int>(std::floor(yf));
            // a point on the domain's upper face belongs to the last cell (decided in space: x/dx can
            // round one ulp above the cell count, as 1000 m on 88 cells does)
            if (ic > ihi && x <= phi[0]) { ic = ihi; }
            if (jc > jhi && y <= phi[1]) { jc = jhi; }
            // the columns are searched within this box and its first ghost layer (filled), so a box split in z,
            // or a fine level that stops below the domain's top, reads only its own data. A point the box owns
            // has its own column's bracketing pair in that range; a neighbouring column whose levels a slope
            // raises or lowers by more than about half a cell, at the box's top or bottom, extrapolates from
            // the range's end instead (a difference of the field's curvature, none for a field linear in height)
            const int kb = amrex::max(klo, vbx.smallEnd(2) - 1);
            const int kt = amrex::min(khi, vbx.bigEnd(2) + 1);
            if (ic < vbx.smallEnd(0) || ic > vbx.bigEnd(0) ||
                jc < vbx.smallEnd(1) || jc > vbx.bigEnd(1)) { return; }
            int kc = -1;
            if (!has_znd) {
                int k = static_cast<int>(std::floor((z - plo[2]) * dxi[2]));
                // a point on the domain's top face belongs to the top cell (the cell, as on the upper x and y faces;
                // there the value is extrapolated from the top cells, as everywhere in z, not read from a ghost cell)
                if (k > khi && z <= phi[2]) { k = khi; }
                if (k >= klo && k <= kt) { kc = k; }
            } else {
                for (int k = kb; k <= kt; ++k) {
                    const Real zb = face_height(2, ic, jc, k,   znd, true, plo[2], dz);
                    const Real zt = face_height(2, ic, jc, k+1, znd, true, plo[2], dz);
                    // the top cell holds its top face too (the domain's top)
                    if (z >= zb && (z < zt || (k == khi && z <= zt))) { kc = k; break; }
                }
                // the bottom face is the mean of its four nodes, but the ground under (x, y) is
                // bilinear between them: on a slope a point just above the ground can lie below
                // the face, and belongs to the bottom cell (read by extrapolation) as long as it is
                // not below the lowest of the four nodes
                if (kc < 0 && kb == klo && z < face_height(2, ic, jc, klo, znd, true, plo[2], dz) &&
                    z >= amrex::min(amrex::min(znd(ic,jc,klo), znd(ic+1,jc,klo)), amrex::min(znd(ic,jc+1,klo), znd(ic+1,jc+1,klo)))) {
                    kc = klo;
                }
            }
            if (kc < vbx.smallEnd(2) || kc > vbx.bigEnd(2)) { return; }
            // bilinear between the cell centres around (x, y), each column at height z
            const Real xi = (x - plo[0]) * dxi[0] - Real(0.5);
            const Real yi = (y - plo[1]) * dxi[1] - Real(0.5);
            const int i0 = static_cast<int>(std::floor(xi));
            const int j0 = static_cast<int>(std::floor(yi));
            const Real wx = xi - i0, wy = yi - j0;
            const Real f00 = centre_column_value(f, comp, i0,   j0,   kb, kt, z, znd, has_znd, plo[2], dz);
            const Real f10 = centre_column_value(f, comp, i0+1, j0,   kb, kt, z, znd, has_znd, plo[2], dz);
            const Real f01 = centre_column_value(f, comp, i0,   j0+1, kb, kt, z, znd, has_znd, plo[2], dz);
            const Real f11 = centre_column_value(f, comp, i0+1, j0+1, kb, kt, z, znd, has_znd, plo[2], dz);
            p_val[p] = (Real(1.0) - wy) * ((Real(1.0) - wx) * f00 + wx * f10) + wy * ((Real(1.0) - wx) * f01 + wx * f11);
            // a count, not a flag: boxes of one rank may run on different streams, and a point held by two of
            // them must still be refused by the caller
            Gpu::Atomic::AddNoRet(&p_cnt[p], 1);
        });
    }
    Gpu::streamSynchronize();
    std::vector<int> cnt(npts, 0);
    Gpu::copy(Gpu::deviceToHost, d_val.begin(), d_val.end(), val.begin());
    Gpu::copy(Gpu::deviceToHost, d_cnt.begin(), d_cnt.end(), cnt.begin());
    ParallelDescriptor::ReduceRealSum(val.data(), npts);
    ParallelDescriptor::ReduceIntSum(cnt.data(), npts);
    for (int p = 0; p < npts; ++p) {
        if (cnt[p] != 1) {
            Abort("sample_cell_scalar: point (" + std::to_string(pos[3*p]) + ", " + std::to_string(pos[3*p+1]) +
                  ", " + std::to_string(pos[3*p+2]) + ") was sampled by " + std::to_string(cnt[p]) +
                  " boxes; the points must lie inside the domain");
        }
    }
}

void
terrain_heights (const MultiFab* z_phys_nd, const Geometry& geom, const std::vector<Real>& pos, std::vector<Real>& h)
{
    BL_PROFILE("erf_actuator::terrain_heights()");
    const int npts = static_cast<int>(pos.size() / 3);
    h.assign(npts, static_cast<Real>(geom.ProbLo(2)));
    check_same_points_everywhere(npts, "terrain_heights");
    if (npts == 0 || z_phys_nd == nullptr) { return; }
    const auto plo = geom.ProbLoArray();
    const auto phi = geom.ProbHiArray();
    const auto dxi = geom.InvCellSizeArray();
    const Box& domain = geom.Domain();
    const int ilo = domain.smallEnd(0), ihi = domain.bigEnd(0);
    const int jlo = domain.smallEnd(1), jhi = domain.bigEnd(1);
    const int klo = domain.smallEnd(2);
    const std::vector<Real> wpos = wrap_periodic(pos, geom);

    Gpu::DeviceVector<Real> d_pos(pos.size()), d_h(npts, 0.0);
    Gpu::DeviceVector<int> d_cnt(npts, 0);
    Gpu::copy(Gpu::hostToDevice, wpos.begin(), wpos.end(), d_pos.begin());
    Real* p_pos = d_pos.data();
    Real* p_h = d_h.data();
    int* p_cnt = d_cnt.data();

    for (MFIter mfi(*z_phys_nd, false); mfi.isValid(); ++mfi) {
        const Box nbx = mfi.validbox();
        if (nbx.smallEnd(2) > klo || nbx.bigEnd(2) < klo) { continue; }
        // the cells whose four lower nodes this nodal box holds; cell boxes are disjoint, so
        // every point has exactly one owner
        const Box cbx = amrex::enclosedCells(nbx);
        Array4<Real const> const& znd = z_phys_nd->const_array(mfi);
        ParallelFor(npts, [=] AMREX_GPU_DEVICE (int p) noexcept
        {
            const Real xi = (p_pos[3*p]   - plo[0]) * dxi[0];
            const Real yi = (p_pos[3*p+1] - plo[1]) * dxi[1];
            int ic = static_cast<int>(std::floor(xi));
            int jc = static_cast<int>(std::floor(yi));
            // a point on the domain's upper face belongs to the last cell
            if (ic > ihi && p_pos[3*p] <= phi[0]) { ic = ihi; }
            if (jc > jhi && p_pos[3*p+1] <= phi[1]) { jc = jhi; }
            if (ic < ilo || ic > ihi || jc < jlo || jc > jhi) { return; }
            if (ic < cbx.smallEnd(0) || ic > cbx.bigEnd(0) || jc < cbx.smallEnd(1) || jc > cbx.bigEnd(1)) { return; }
            const Real wx = xi - static_cast<Real>(ic), wy = yi - static_cast<Real>(jc);
            const Real z00 = znd(ic,   jc,   klo), z10 = znd(ic+1, jc,   klo);
            const Real z01 = znd(ic,   jc+1, klo), z11 = znd(ic+1, jc+1, klo);
            p_h[p] = (Real(1.0) - wy) * ((Real(1.0) - wx) * z00 + wx * z10) + wy * ((Real(1.0) - wx) * z01 + wx * z11);
            // a count, not a flag: boxes of one rank may run on different streams, and a point held by two of
            // them must still be refused by the caller
            Gpu::Atomic::AddNoRet(&p_cnt[p], 1);
        });
    }
    Gpu::streamSynchronize();
    std::vector<int> cnt(npts, 0);
    Gpu::copy(Gpu::deviceToHost, d_h.begin(), d_h.end(), h.begin());
    Gpu::copy(Gpu::deviceToHost, d_cnt.begin(), d_cnt.end(), cnt.begin());
    ParallelDescriptor::ReduceRealSum(h.data(), npts);
    ParallelDescriptor::ReduceIntSum(cnt.data(), npts);
    for (int p = 0; p < npts; ++p) {
        if (cnt[p] != 1) {
            Abort("terrain_heights: the point (" + std::to_string(pos[3*p]) + ", " + std::to_string(pos[3*p+1]) +
                  ") was found in " + std::to_string(cnt[p]) + " boxes: the point must lie inside the domain, over grids of "
                  "this level that reach the ground (a refined level aloft holds no ground)");
        }
    }
}

namespace {

// The table of this rank's boxes: each node level's lowest and highest height, and the cells the boxes
// hold. One kernel entry per node level of a box, which goes over that level's nodes and then makes one
// atomic min and one atomic max into the level's two values: the kernels of two boxes may run at the same
// time on GPU (MFIter spreads the boxes over the device's streams), and both may hold the level
void
node_level_table (const MultiFab& z_phys_nd, const Geometry& geom, ZBounds& b)
{
    const Real big = std::numeric_limits<Real>::max();
    const Box& dom = geom.Domain();
    b.k0 = dom.smallEnd(2);
    const int nnode = dom.length(2) + 1;
    const Box dom_nodes = amrex::surroundingNodes(dom);
    Gpu::DeviceVector<Real> lo(nnode, big), hi(nnode, -big);
    Real* plo = lo.data();
    Real* phi = hi.data();
    const int k0 = b.k0;
    b.has_cell.assign(nnode - 1, 0);
    for (MFIter mfi(z_phys_nd, false); mfi.isValid(); ++mfi) {
        const Box nbx = mfi.validbox();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dom_nodes.contains(nbx), "mesh_z_bounds: z_phys_nd's boxes must lie in geom's domain");
        Array4<Real const> const& znd = z_phys_nd.const_array(mfi);
        const int ilo = nbx.smallEnd(0), ihi = nbx.bigEnd(0);
        const int jlo = nbx.smallEnd(1), jhi = nbx.bigEnd(1);
        const int klo = nbx.smallEnd(2);
        ParallelFor(nbx.length(2), [=] AMREX_GPU_DEVICE (int n) noexcept
        {
            const int k = klo + n;
            Real mn = big;
            Real mx = -big;
            for (int j = jlo; j <= jhi; ++j) {
                for (int i = ilo; i <= ihi; ++i) {
                    mn = amrex::min(mn, znd(i,j,k));
                    mx = amrex::max(mx, znd(i,j,k));
                }
            }
            Gpu::Atomic::Min(plo + (k - k0), mn);
            Gpu::Atomic::Max(phi + (k - k0), mx);
        });
        // the box's cells: between its lowest and highest node levels
        for (int k = nbx.smallEnd(2); k < nbx.bigEnd(2); ++k) { b.has_cell[k - b.k0] = 1; }
    }
    b.node_lo.resize(nnode);
    b.node_hi.resize(nnode);
    Gpu::copyAsync(Gpu::deviceToHost, lo.begin(), lo.end(), b.node_lo.begin());
    Gpu::copyAsync(Gpu::deviceToHost, hi.begin(), hi.end(), b.node_hi.begin());
    Gpu::streamSynchronize();
}

} // namespace

ZBounds
mesh_z_bounds (const MultiFab& z_phys_nd, const Geometry& geom, bool over_all_ranks)
{
    BL_PROFILE("erf_actuator::mesh_z_bounds()");
    const int kground = geom.Domain().smallEnd(2);
    const Real big = std::numeric_limits<Real>::max();
    ReduceOps<ReduceOpMin, ReduceOpMax, ReduceOpMin, ReduceOpMax> ops;
    ReduceData<Real, Real, Real, Real> data(ops);
    using T = typename decltype(data)::Type;
    for (MFIter mfi(z_phys_nd, false); mfi.isValid(); ++mfi) {
        const Box nbx = mfi.validbox();
        const int ktop = nbx.bigEnd(2);
        Array4<Real const> const& znd = z_phys_nd.const_array(mfi);
        ops.eval(nbx, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> T
        {
            const bool ground = (k == kground);
            const bool above = (k < ktop);   // the spacing to the node above, in this box
            const Real dz = above ? znd(i,j,k+1) - znd(i,j,k) : Real(0.0);
            return {ground ? znd(i,j,k) : big, ground ? znd(i,j,k) : -big, above ? dz : big, above ? dz : -big};
        });
    }
    const T r = data.value(ops);
    ZBounds b;
    b.ground_lo = amrex::get<0>(r);
    b.ground_hi = amrex::get<1>(r);
    b.dz_lo = amrex::get<2>(r);
    b.dz_hi = amrex::get<3>(r);
    int all = 1;
    for (MFIter mfi(z_phys_nd, false); mfi.isValid(); ++mfi) { all = (mfi.validbox().smallEnd(2) == kground) ? all : 0; }
    if (over_all_ranks) {
        ParallelDescriptor::ReduceRealMin(b.ground_lo);
        ParallelDescriptor::ReduceRealMax(b.ground_hi);
        ParallelDescriptor::ReduceRealMin(b.dz_lo);
        ParallelDescriptor::ReduceRealMax(b.dz_hi);
        ParallelDescriptor::ReduceIntMin(all);
    }
    b.all_ground = (all != 0);
    if (!over_all_ranks && !b.all_ground) {
        // where some of this rank's boxes do not reach the ground
        node_level_table(z_phys_nd, geom, b);
    }
    return b;
}

ZBounds
mesh_z_levels (const MultiFab& z_phys_nd, const Geometry& geom)
{
    BL_PROFILE("erf_actuator::mesh_z_levels()");
    ZBounds b = mesh_z_bounds(z_phys_nd, geom, true);
    node_level_table(z_phys_nd, geom, b);
    // over every rank's boxes: a rank with no box at a node level holds the largest and the lowest Real there
    const int nnode = static_cast<int>(b.node_lo.size());
    ParallelDescriptor::ReduceRealMin(b.node_lo.data(), nnode);
    ParallelDescriptor::ReduceRealMax(b.node_hi.data(), nnode);
    std::vector<int> has(b.has_cell.begin(), b.has_cell.end());
    ParallelDescriptor::ReduceIntMax(has.data(), static_cast<int>(has.size()));
    for (std::size_t n = 0; n < has.size(); ++n) { b.has_cell[n] = static_cast<char>(has[n]); }
    return b;
}

const ZBounds&
LevelHeights::get (const MultiFab& z_phys_nd, const Geometry& geom, bool rebuild)
{
    if (rebuild || m_builds == 0 || m_domain != geom.Domain() || m_ba != z_phys_nd.boxArray()) {
        m_zb = mesh_z_levels(z_phys_nd, geom);
        m_ba = z_phys_nd.boxArray();
        m_domain = geom.Domain();
        ++m_builds;
    }
    return m_zb;
}

namespace {

// Per direction across, the node index ranges of the cells a region touches, from its low (b[0], b[1])
// and high (b[3], b[4]) ends in metres: one range, or two in a periodic direction whose range crosses
// the seam (the part beyond it wrapped to the other side)
void
horizontal_node_ranges (const Real* b, const Geometry& geom, std::array<std::vector<std::array<int,2>>,2>& ranges)
{
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const Box& domain = geom.Domain();
    for (int d = 0; d < 2; ++d) {
        ranges[d].clear();
        const int dlo = domain.smallEnd(d), dhi = domain.bigEnd(d), nc = dhi - dlo + 1;
        // a bound more than a domain length outside (a diverged point) is clamped there, so the
        // index stays an int; such a region already covers or misses the whole domain
        const Real len = Real(nc) / dxi[d];
        const Real blo = amrex::Clamp(b[d],   plo[d] - len, plo[d] + Real(2.0) * len);
        const Real bhi = amrex::Clamp(b[d+3], plo[d] - len, plo[d] + Real(2.0) * len);
        const int clo = static_cast<int>(std::floor((blo - plo[d]) * dxi[d]));
        const int chi = static_cast<int>(std::floor((bhi - plo[d]) * dxi[d]));
        if (geom.isPeriodic(d)) {
            const int nnodes = chi + 2 - clo;   // nodes clo .. chi + 1
            if (nnodes >= nc) {
                ranges[d].push_back({{dlo, dhi + 1}});
            } else {
                const int a = dlo + (((clo - dlo) % nc) + nc) % nc;
                const int e = a + nnodes - 1;
                if (e <= dhi + 1) {
                    ranges[d].push_back({{a, e}});
                } else {
                    ranges[d].push_back({{a, dhi + 1}});
                    ranges[d].push_back({{dlo, e - nc}});
                }
            }
        } else {
            const int a = std::max(clo, dlo), e = std::min(chi + 1, dhi + 1);
            if (a <= e) { ranges[d].push_back({{a, e}}); }
        }
    }
}

} // namespace

namespace {

// one piece of the work of column_cell_ranges: point p over the node columns i0..i1 x j0..j1 of box b
struct ColumnWork
{
    int p, b, i0, i1, j0, j1;
};

} // namespace

void
column_cell_ranges (const MultiFab& z_phys_nd, const Geometry& geom, const std::vector<Real>& pos,
                    const std::vector<Real>& reach, ColumnCells& cells)
{
    BL_PROFILE("erf_actuator::column_cell_ranges()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pos.size() == 3 * reach.size(), "column_cell_ranges: one reach per point");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!geom.isPeriodic(2), "column_cell_ranges: z must not be periodic");
    const int npts = static_cast<int>(reach.size());
    const Box& domain = geom.Domain();
    const int dlo = domain.smallEnd(2), dhi = domain.bigEnd(2);
    // nothing touched: an empty range (kb > kt), the domain's cells being dlo .. dhi
    cells.kb.assign(npts, dhi + 1);
    cells.kt.assign(npts, dlo - 1);
    check_same_points_everywhere(npts, "column_cell_ranges");
    if (npts == 0) { return; }
    // the pieces: for every point, the node columns its reach across touches in each of this rank's boxes
    std::vector<ColumnWork> work;
    std::vector<Box> boxes;
    std::vector<int> grid;   // the box's index in z_phys_nd
    for (MFIter mfi(z_phys_nd, false); mfi.isValid(); ++mfi) { boxes.push_back(mfi.validbox()); grid.push_back(mfi.index()); }
    std::array<std::vector<std::array<int,2>>,2> ranges;
    for (int p = 0; p < npts; ++p) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(pos[3*p]) && std::isfinite(pos[3*p+1]) && std::isfinite(pos[3*p+2]) &&
                                         std::isfinite(reach[p]) && reach[p] >= Real(0.0),
                                         "column_cell_ranges: the points and their reach must be finite, the reach not negative");
        const Real b[6] = {pos[3*p] - reach[p], pos[3*p+1] - reach[p], Real(0.0), pos[3*p] + reach[p], pos[3*p+1] + reach[p], Real(0.0)};
        horizontal_node_ranges(b, geom, ranges);
        for (int ib = 0; ib < static_cast<int>(boxes.size()); ++ib) {
            const Box& nbx = boxes[ib];
            for (const auto& ri : ranges[0]) {
                for (const auto& rj : ranges[1]) {
                    const int i0 = std::max(ri[0], nbx.smallEnd(0)), i1 = std::min(ri[1], nbx.bigEnd(0));
                    const int j0 = std::max(rj[0], nbx.smallEnd(1)), j1 = std::min(rj[1], nbx.bigEnd(1));
                    if (i0 <= i1 && j0 <= j1) { work.push_back({p, ib, i0, i1, j0, j1}); }
                }
            }
        }
    }
    if (!work.empty()) {
        // the pieces sorted by box (one launch per box), with each box's first piece
        std::stable_sort(work.begin(), work.end(), [] (const ColumnWork& a, const ColumnWork& c) { return a.b < c.b; });
        std::vector<int> first_of(boxes.size() + 1, static_cast<int>(work.size()));
        for (int w = static_cast<int>(work.size()) - 1; w >= 0; --w) { first_of[work[w].b] = w; }
        for (int ib = static_cast<int>(boxes.size()) - 1; ib >= 0; --ib) { first_of[ib] = std::min(first_of[ib], first_of[ib+1]); }
        // every device vector lives until the synchronisation after the last launch
        Gpu::DeviceVector<ColumnWork> d_work(work.size());
        Gpu::DeviceVector<Real> d_pos(pos.size()), d_reach(reach.size());
        Gpu::DeviceVector<int> d_kb(npts, dhi + 1), d_kt(npts, dlo - 1);
        Gpu::copy(Gpu::hostToDevice, work.begin(), work.end(), d_work.begin());
        Gpu::copy(Gpu::hostToDevice, pos.begin(), pos.end(), d_pos.begin());
        Gpu::copy(Gpu::hostToDevice, reach.begin(), reach.end(), d_reach.begin());
        const ColumnWork* p_work = d_work.data();
        const Real* p_pos = d_pos.data();
        const Real* p_reach = d_reach.data();
        int* p_kb = d_kb.data();
        int* p_kt = d_kt.data();
        // the pieces of one box in one launch: each over its node columns, every cell of the column
        for (int ib = 0; ib < static_cast<int>(boxes.size()); ++ib) {
            const int w0 = first_of[ib], nw = first_of[ib+1] - first_of[ib];
            if (nw == 0) { continue; }
            const Box nbx = boxes[ib];
            const int kbot = nbx.smallEnd(2), ktop = nbx.bigEnd(2);
            Array4<Real const> const& znd = z_phys_nd.const_array(grid[ib]);
            ParallelFor(nw, [=] AMREX_GPU_DEVICE (int m) noexcept
            {
                const ColumnWork w = p_work[w0 + m];
                const Real zl0 = p_pos[3*w.p+2] - p_reach[w.p], zh0 = p_pos[3*w.p+2] + p_reach[w.p];
                // the slack of ZBounds::k_range: the heights are rounded sums
                const Real slack = Real(1.0e-5) * amrex::max(amrex::Math::abs(zl0), amrex::max(amrex::Math::abs(zh0), zh0 - zl0));
                const Real zl = zl0 - slack, zh = zh0 + slack;
                // from the empty range, not from the point's slot: another piece of the same point (a range
                // split at a periodic seam) may be updating it
                int kb = dhi + 1, kt = dlo - 1;
                for (int j = w.j0; j <= w.j1; ++j) {
                    for (int i = w.i0; i <= w.i1; ++i) {
                        // the cells of this column whose spacing the heights overlap
                        for (int k = kbot; k < ktop; ++k) {
                            if (znd(i,j,k) <= zh && znd(i,j,k+1) >= zl) { kb = amrex::min(kb, k); kt = amrex::max(kt, k); }
                        }
                        // below the box's bottom face (not the floor): the cell under it; above its top (not the domain's): the cell over it
                        if (kbot > dlo && znd(i,j,kbot) > zl) { kb = amrex::min(kb, kbot - 1); }
                        if (ktop < dhi + 1 && znd(i,j,ktop) < zh) { kt = amrex::max(kt, ktop); }
                    }
                }
                Gpu::Atomic::Min(p_kb + w.p, kb);
                Gpu::Atomic::Max(p_kt + w.p, kt);
            });
        }
        Gpu::streamSynchronize();
        Gpu::copy(Gpu::deviceToHost, d_kb.begin(), d_kb.end(), cells.kb.begin());
        Gpu::copy(Gpu::deviceToHost, d_kt.begin(), d_kt.end(), cells.kt.begin());
    }
    ParallelDescriptor::ReduceIntMin(cells.kb.data(), npts);
    ParallelDescriptor::ReduceIntMax(cells.kt.data(), npts);
    // nothing overlapped (a reach wholly above the level's top or below its bottom, or beside it across): one
    // empty range, whatever the boxes (a box's bottom or top face the heights lie beyond leaves its own cell)
    for (int p = 0; p < npts; ++p) {
        if (cells.kb[p] > cells.kt[p]) { cells.kb[p] = dhi + 1; cells.kt[p] = dlo - 1; }
    }
}

CoverageGrids::CoverageGrids (const BoxArray& ba, const Geometry& geom)
    : cells(amrex::convert(ba, IntVect::TheZeroVector()))
{
    BL_PROFILE("erf_actuator::CoverageGrids()");
    const Box& domain = geom.Domain();
    BoxList flat;
    for (int b = 0; b < static_cast<int>(cells.size()); ++b) {
        Box f = cells[b];
        f.setSmall(2, domain.smallEnd(2));
        f.setBig(2, domain.smallEnd(2));
        flat.push_back(f);
    }
    footprints = BoxArray(flat);
}

bool
points_covered_by (const BoxArray& ba, const Geometry& geom, const std::vector<Real>& pos, Real reach, std::string& first_outside,
                   CoverZ z, const ZBounds* zb, const ColumnCells* cells)
{
    const CoverageGrids grids(ba, geom);
    return points_covered_by(grids, geom, pos, reach, first_outside, z, zb, cells);
}

bool
points_covered_by (const CoverageGrids& grids, const Geometry& geom, const std::vector<Real>& pos, Real reach,
                   std::string& first_outside, CoverZ z, const ZBounds* zb, const ColumnCells* cells)
{
    BL_PROFILE("erf_actuator::points_covered_by()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!geom.isPeriodic(2), "points_covered_by: z must not be periodic");
    // bounds of a level with no box at the ground (ground_lo > ground_hi) carry no spacing, and the whole
    // column is asked for; any other bounds must have a positive spacing
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(zb == nullptr || zb->ground_lo > zb->ground_hi ||
                                     (zb->dz_lo > Real(0.0) && zb->dz_hi >= zb->dz_lo),
                                     "points_covered_by: the levels' spacing must be positive");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(zb == nullptr || !zb->has_table(),
                                     "points_covered_by: zb must be the bounds over all ranks (mesh_z_bounds with over_all_ranks = true), "
                                     "not one rank's table");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(cells == nullptr || (cells->kb.size() == pos.size() / 3 && cells->kt.size() == pos.size() / 3),
                                     "points_covered_by: one cell range per point");
    first_outside.clear();
    const BoxArray& cc = grids.cells;
    const BoxArray& footprints = grids.footprints;
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const Box& domain = geom.Domain();
    for (std::size_t p = 0; p < pos.size() / 3; ++p) {
        // per direction, the index ranges of the cells needed: one range, or two in a periodic
        // direction whose range crosses the seam (the part beyond it wrapped to the other side)
        std::array<std::vector<std::array<int,2>>,3> ranges;
        bool none = false;   // the point needs cells the level has nowhere
        for (int d = 0; d < 3; ++d) {
            int lo = static_cast<int>(std::floor((pos[3*p+d] - reach - plo[d]) * dxi[d]));
            int hi = static_cast<int>(std::floor((pos[3*p+d] + reach - plo[d]) * dxi[d]));
            const int dlo = domain.smallEnd(d), dhi = domain.bigEnd(d), n = dhi - dlo + 1;
            if (geom.isPeriodic(d)) {
                // shift by whole periods so that lo lies in the domain
                const int s = static_cast<int>(std::floor(static_cast<double>(lo - dlo) / n));
                lo -= s * n;
                hi -= s * n;
                if (hi - lo + 1 >= n) {
                    ranges[d].push_back({{dlo, dhi}});
                } else if (hi <= dhi) {
                    ranges[d].push_back({{lo, hi}});
                } else {
                    ranges[d].push_back({{lo, dhi}});
                    ranges[d].push_back({{dlo, hi - n}});
                }
            } else if (d == 2 && z == CoverZ::Footprint) {
                ranges[d].push_back({{dlo, dlo}});
            } else if (d == 2 && z == CoverZ::Column) {
                if (cells != nullptr) {
                    // the cells of the point's own columns (column_cell_ranges); none: not covered
                    if (cells->kb[p] > cells->kt[p]) { none = true; }
                    ranges[d].push_back({{std::max(cells->kb[p], dlo), std::min(cells->kt[p], dhi)}});
                } else if (zb == nullptr) {
                    ranges[d].push_back({{dlo, dhi}});
                } else {
                    // every cell a height within reach may lie in, whatever column it is in (the global
                    // bounds: a rank's own table would make the answer depend on the rank)
                    int kb = dlo, kt = dhi;
                    zb->k_range(pos[3*p+2] - reach, pos[3*p+2] + reach, dlo, dhi, kb, kt);
                    ranges[d].push_back({{std::max(kb, dlo), std::min(kt, dhi)}});
                }
            } else {
                ranges[d].push_back({{std::max(lo, dlo), std::min(hi, dhi)}});
            }
        }
        bool covered = !none;
        for (const auto& rx : ranges[0]) {
            if (!covered) { break; }
            for (const auto& ry : ranges[1]) {
                for (const auto& rz : ranges[2]) {
                    const Box needed(IntVect(rx[0], ry[0], rz[0]), IntVect(rx[1], ry[1], rz[1]));
                    // a range outside a non-periodic domain gives an empty box, which is not covered
                    if ((z == CoverZ::Footprint) ? !footprints.contains(needed, false) : !cc.contains(needed, true)) { covered = false; }
                }
            }
        }
        if (!covered) {
            first_outside = "(" + std::to_string(pos[3*p]) + ", " + std::to_string(pos[3*p+1]) + ", " + std::to_string(pos[3*p+2]) + ") m";
            return false;
        }
    }
    return true;
}

bool
region_touches_grids (const BoxArray& ba, const Geometry& geom, const RealBox& region, const ZBounds* zb)
{
    BL_PROFILE("erf_actuator::region_touches_grids()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!geom.isPeriodic(2), "region_touches_grids: z must not be periodic");
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        // a NaN first: comparing one may trap (amrex.fpe_trap_invalid) before the message is printed
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!std::isnan(region.lo(d)) && !std::isnan(region.hi(d)) && region.lo(d) <= region.hi(d),
                                         "region_touches_grids: the region needs lo <= hi in every direction");
    }
    const Box& domain = geom.Domain();
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    IntVect lo, hi;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        // The region's cells, as points_covered_by finds a point's, kept as doubles until they are bounded: an int
        // holds no infinity and no index far outside the domain
        const int dlo = domain.smallEnd(d), dhi = domain.bigEnd(d), n = dhi - dlo + 1;
        double l = std::floor(static_cast<double>((region.lo(d) - plo[d]) * dxi[d]));
        double h = std::floor(static_cast<double>((region.hi(d) - plo[d]) * dxi[d]));
        if (geom.isPeriodic(d)) {
            // an infinite end is replaced by 0 before any difference is taken: inf - inf is invalid and may
            // trap (amrex.fpe_trap_invalid), also where a compiler evaluates it ahead of a test that skips it
            const bool unbounded = std::isinf(l) || std::isinf(h);
            const double fl = unbounded ? 0.0 : l;
            const double fh = unbounded ? 0.0 : h;
            if (unbounded || !(fh - fl + 1 < n)) {
                // a period or more, or an infinite extent: every cell
                l = dlo;
                h = dhi;
            } else {
                // by whole periods until lo lies in the domain; hi is then less than a period above it, and the
                // shifts below bring the part beyond the seam back in
                const double s = std::floor((fl - dlo) / n);
                l = fl - s * n;
                h = fh - s * n;
            }
        }
        // one cell outside the domain stays outside it (and bounds what rounding leaves of a far periodic region)
        lo[d] = static_cast<int>(std::min<double>(static_cast<double>(dhi) + n, std::max<double>(static_cast<double>(dlo) - 1, l)));
        hi[d] = static_cast<int>(std::min<double>(static_cast<double>(dhi) + n, std::max<double>(static_cast<double>(dlo) - 1, h)));
    }
    if (zb != nullptr) {
        // every cell a height of the region may lie in, whatever its column. With the table, from the level's own
        // node heights. With the bounds alone, from the ground and the spacing, which describe a column only where
        // every box stands on the ground; with a box that does not, the column is taken whole
        int kb = domain.smallEnd(2), kt = domain.bigEnd(2);
        if (zb->has_table() || zb->all_ground) { zb->k_range(region.lo(2), region.hi(2), domain.smallEnd(2), domain.bigEnd(2), kb, kt); }
        lo[2] = kb;
        hi[2] = kt;
    }
    if (lo[2] > hi[2]) { return false; }   // the table holds no cell at these heights
    const BoxArray cc = amrex::convert(ba, IntVect::TheZeroVector());
    const Box needed(lo, hi);
    for (const IntVect& s : geom.periodicity().shiftIntVect()) {
        const Box b = (needed + s) & domain;
        if (b.ok() && cc.intersects(b)) { return true; }
    }
    return false;
}

} // namespace erf_actuator
