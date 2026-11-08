// Point sampling of face and cell fields and of the terrain height, and the coverage check.

#include "ERF_ActuatorSampling.H"

#include <cmath>
#include <array>
#include <limits>
#include <string>
#include <vector>

#include <AMReX.H>
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
// the ground and up to the top without reading the ghost cells below the ground or above the top.
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
sample_velocity (const MultiFab& U, const MultiFab& V, const MultiFab& W,
                 const MultiFab* z_phys_nd, const Geometry& geom,
                 const std::vector<Real>& pos, std::vector<Real>& vel)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pos.size() % 3 == 0, "sample_velocity: pos holds x,y,z triples");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(U.nGrow() >= 1 && V.nGrow() >= 1 && W.nGrow() >= 1,
                                     "sample_velocity needs one filled ghost cell on the velocities");
    check_layout(U, &V, 1, "sample_velocity: V");
    check_layout(U, &W, 1, "sample_velocity: W");
    check_layout(U, z_phys_nd, 1, "sample_velocity: z_phys_nd");
    const int npts = static_cast<int>(pos.size() / 3);
    vel.assign(3 * static_cast<std::size_t>(npts), Real(0.0));
    check_same_points_everywhere(npts, "sample_velocity");
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
    const std::vector<Real> wpos = wrap_periodic(pos, geom);

    Gpu::DeviceVector<Real> d_pos(pos.size());
    Gpu::DeviceVector<Real> d_vel(pos.size(), Real(0.0));
    Gpu::DeviceVector<int>  d_cnt(npts, 0);
    Gpu::copyAsync(Gpu::hostToDevice, wpos.begin(), wpos.end(), d_pos.begin());
    Real* p_pos = d_pos.data();
    Real* p_vel = d_vel.data();
    int*  p_cnt = d_cnt.data();

    for (MFIter mfi(U, false); mfi.isValid(); ++mfi) {
        // the cells this rank owns; a point is sampled by the box holding its containing cell
        const Box vbx = amrex::enclosedCells(mfi.validbox());
        Array4<Real const> const& u = U.const_array(mfi);
        Array4<Real const> const& v = V.const_array(mfi);
        Array4<Real const> const& w = W.const_array(mfi);
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
            // the containing cell in z: between the z faces of column (ic,jc)
            int kc = -1;
            if (!has_znd) {
                int k = static_cast<int>(std::floor((z - plo[2]) * dxi[2]));
                // a point on the domain's top belongs to the top cell, as on the upper x and y faces
                if (k > khi && z <= phi[2]) { k = khi; }
                if (k >= klo && k <= kt) { kc = k; }
            } else {
                for (int k = kb; k <= kt; ++k) {
                    const Real zb = face_height(2, ic, jc, k,   znd, true, plo[2], dz);
                    const Real zt = face_height(2, ic, jc, k+1, znd, true, plo[2], dz);
                    // the top cell takes its top face too: the domain's top, as the upper x and y faces
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
            p_cnt[p] += 1;
        });
    }
    Gpu::streamSynchronize();

    std::vector<int> cnt(npts, 0);
    Gpu::copy(Gpu::deviceToHost, d_vel.begin(), d_vel.end(), vel.begin());
    Gpu::copy(Gpu::deviceToHost, d_cnt.begin(), d_cnt.end(), cnt.begin());
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
                // a point on the domain's top belongs to the top cell, as on the upper x and y faces
                if (k > khi && z <= phi[2]) { k = khi; }
                if (k >= klo && k <= kt) { kc = k; }
            } else {
                for (int k = kb; k <= kt; ++k) {
                    const Real zb = face_height(2, ic, jc, k,   znd, true, plo[2], dz);
                    const Real zt = face_height(2, ic, jc, k+1, znd, true, plo[2], dz);
                    // the top cell takes its top face too: the domain's top, as the upper x and y faces
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
            p_cnt[p] += 1;
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
            p_cnt[p] += 1;
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
                  ") was found in " + std::to_string(cnt[p]) + " boxes; the bodies' bases must lie inside the domain");
        }
    }
}

ZBounds
mesh_z_bounds (const MultiFab& z_phys_nd, const Geometry& geom)
{
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
    ParallelDescriptor::ReduceRealMin(b.ground_lo);
    ParallelDescriptor::ReduceRealMax(b.ground_hi);
    ParallelDescriptor::ReduceRealMin(b.dz_lo);
    ParallelDescriptor::ReduceRealMax(b.dz_hi);
    return b;
}

bool
points_covered_by (const BoxArray& ba, const Geometry& geom, const std::vector<Real>& pos, Real reach, std::string& first_outside,
                   CoverZ z, const ZBounds* zb)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!geom.isPeriodic(2), "points_covered_by: z must not be periodic");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(zb == nullptr || (zb->dz_lo > Real(0.0) && zb->dz_hi >= zb->dz_lo),
                                     "points_covered_by: the levels' spacing must be positive");
    first_outside.clear();
    const BoxArray cc = amrex::convert(ba, IntVect::TheZeroVector());
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const Box& domain = geom.Domain();
    // for Footprint, each box flattened onto the bottom layer (they may overlap)
    BoxList flat;
    if (z == CoverZ::Footprint) {
        for (int b = 0; b < static_cast<int>(cc.size()); ++b) {
            Box f = cc[b];
            f.setSmall(2, domain.smallEnd(2));
            f.setBig(2, domain.smallEnd(2));
            flat.push_back(f);
        }
    }
    const BoxArray footprints(flat);
    for (std::size_t p = 0; p < pos.size() / 3; ++p) {
        // per direction, the index ranges of the cells needed: one range, or two in a periodic
        // direction whose range crosses the seam (the part beyond it wrapped to the other side)
        std::array<std::vector<std::array<int,2>>,3> ranges;
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
                if (zb == nullptr) {
                    ranges[d].push_back({{dlo, dhi}});
                } else {
                    // every cell a height within reach may lie in, whatever column it is in
                    const int kb = static_cast<int>(std::floor((pos[3*p+2] - reach - zb->ground_hi) / zb->dz_hi));
                    const int kt = static_cast<int>(std::floor((pos[3*p+2] + reach - zb->ground_lo) / zb->dz_lo));
                    ranges[d].push_back({{std::max(kb, dlo), std::min(kt, dhi)}});
                }
            } else {
                ranges[d].push_back({{std::max(lo, dlo), std::min(hi, dhi)}});
            }
        }
        bool covered = true;
        for (const auto& rx : ranges[0]) {
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

} // namespace erf_actuator
