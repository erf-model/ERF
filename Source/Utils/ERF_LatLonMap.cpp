/**
 * \file ERF_LatLonMap.cpp
 */
#include "ERF_LatLonMap.H"
#include "ERF_NumericalConstants.H"

#include <algorithm>

#include <AMReX_Arena.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>

using namespace amrex;

namespace {

// The centred-difference stencil around (i,j), collapsed to one side at the
// edge of the domain
struct GradientStencil
{
    int im, ip, jm, jp;
    Real inv_di, inv_dj;
};

GradientStencil
gradient_stencil (const Box& dom, int i, int j)
{
    GradientStencil s{};
    s.im = std::max(i-1, dom.smallEnd(0));
    s.ip = std::min(i+1, dom.bigEnd(0));
    s.jm = std::max(j-1, dom.smallEnd(1));
    s.jp = std::min(j+1, dom.bigEnd(1));
    s.inv_di = (s.ip > s.im) ? Real(1.0)/Real(s.ip - s.im) : Real(0.0);
    s.inv_dj = (s.jp > s.jm) ? Real(1.0)/Real(s.jp - s.jm) : Real(0.0);
    return s;
}

} // namespace

//
// Copy the latitude and longitude into the two components of latlon.  A free
// function: nvcc does not allow an extended device lambda in a constructor.
//
void
pack_latlon (MultiFab& latlon, const MultiFab& lat, const MultiFab& lon)
{
    for (MFIter mfi(latlon, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const Array4<Real>&       ll  = latlon.array(mfi);
        const Array4<const Real>& lat_arr = lat.const_array(mfi);
        const Array4<const Real>& lon_arr = lon.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            ll(i,j,k,0) = lat_arr(i,j,0);
            ll(i,j,k,1) = lon_arr(i,j,0);
        });
    }
}

void
grid_rotation_from_latlon (const Array4<const Real>& ll, const Box& dom, int i, int j,
                           Real& cos_alpha, Real& sin_alpha)
{
    const GradientStencil s = gradient_stencil(dom, i, j);

    // Displacement on the sphere per cell in i: east and north components, in
    // degrees of latitude (the common factor, the Earth's radius, cancels)
    const Real lat   = ll(i,j,0,0);
    const Real north = (ll(s.ip,j,0,0) - ll(s.im,j,0,0)) * s.inv_di;
    const Real east  = wrap_longitude_difference(ll(s.ip,j,0,1) - ll(s.im,j,0,1)) * s.inv_di
                     * std::cos(lat * PI / Real(180.0));
    const Real norm  = std::sqrt(east*east + north*north);

    if (norm > Real(0.0)) {
        cos_alpha = east  / norm;
        sin_alpha = north / norm;
    } else {
        // A domain one cell wide in x has no i direction to measure
        cos_alpha = Real(1.0);
        sin_alpha = Real(0.0);
    }
}

LatLonStatus
locate_latlon_on_grid (const Array4<const Real>& ll, const Box& dom,
                       const GpuArray<Real, AMREX_SPACEDIM>& problo,
                       const GpuArray<Real, AMREX_SPACEDIM>& dx,
                       Real req_lat, Real req_lon, LatLonLocation& loc)
{
    const int ilo = dom.smallEnd(0), ihi = dom.bigEnd(0);
    const int jlo = dom.smallEnd(1), jhi = dom.bigEnd(1);

    const Real cosfac = std::cos(req_lat * PI / Real(180.0));
    auto dist2 = [&](int i, int j) {
        const Real dlat = ll(i,j,0,0) - req_lat;
        const Real dlon = wrap_longitude_difference(ll(i,j,0,1) - req_lon) * cosfac;
        return dlat*dlat + dlon*dlon;
    };

    int  bi = ilo, bj = jlo;
    Real best = dist2(ilo,jlo);
    for (int j = jlo; j <= jhi; ++j) {
        for (int i = ilo; i <= ihi; ++i) {
            const Real d = dist2(i,j);
            if (d < best) { best = d; bi = i; bj = j; }
        }
    }

    loc.near_lat = ll(bi,bj,0,0);
    loc.near_lon = ll(bi,bj,0,1);

    // One linear solve in index space inverts the map: over a cell the
    // projection is linear to well below a cell width.
    const GradientStencil s = gradient_stencil(dom, bi, bj);

    const Real dlat_di = (ll(s.ip,bj,0,0) - ll(s.im,bj,0,0)) * s.inv_di;
    const Real dlon_di = wrap_longitude_difference(ll(s.ip,bj,0,1) - ll(s.im,bj,0,1)) * s.inv_di;
    const Real dlat_dj = (ll(bi,s.jp,0,0) - ll(bi,s.jm,0,0)) * s.inv_dj;
    const Real dlon_dj = wrap_longitude_difference(ll(bi,s.jp,0,1) - ll(bi,s.jm,0,1)) * s.inv_dj;

    const Real rlat = req_lat - ll(bi,bj,0,0);
    const Real rlon = wrap_longitude_difference(req_lon - ll(bi,bj,0,1));
    Real di = Real(0.0), dj = Real(0.0);
    if (!invert_latlon_offset(rlat, rlon, dlat_di, dlon_di, dlat_dj, dlon_dj, di, dj)) {
        return LatLonStatus::Degenerate;
    }

    // The solve places the point relative to the nearest grid point, so more
    // than a cell away means the request cannot be trusted: either it is
    // outside the domain, or the map is too poorly conditioned near it for one
    // linear solve.
    if (std::abs(di) > Real(1.0) || std::abs(dj) > Real(1.0)) {
        return LatLonStatus::TooFar;
    }

    loc.x   = problo[0] + (Real(bi) + Real(0.5) + di) * dx[0];
    loc.y   = problo[1] + (Real(bj) + Real(0.5) + dj) * dx[1];
    loc.lat = ll(bi,bj,0,0) + dlat_di*di + dlat_dj*dj;
    // Resolving a point from a nearest grid point on the other side of the
    // antimeridian steps the longitude just past +-180; report the longitude
    // that names the place, not the one that runs off the end of the range
    loc.lon = wrap_longitude_difference(ll(bi,bj,0,1) + dlon_di*di + dlon_dj*dj);

    grid_rotation_from_latlon(ll, dom, bi, bj, loc.cos_alpha, loc.sin_alpha);

    return LatLonStatus::Ok;
}

LatLonMap::LatLonMap (const MultiFab& lat, const MultiFab& lon,
                      const BoxArray& ba2d, const DistributionMapping& dm,
                      const Geometry& geom0)
    : m_problo(geom0.ProbLoArray()),
      m_dx(geom0.CellSizeArray())
{
    m_dom = geom0.Domain();
    m_dom.setRange(2,0);

    // Gather the level-0 mass-point latitude and longitude onto the IO rank.
    // This is a setup-time, level-0, 2D array, so the gather is affordable and
    // keeps the search and the inverse map as plain host code.
    //
    // The ceiling is one rank holding 2 Reals per level-0 column -- 100 MB or
    // so for the largest domains ERF is run on -- and a search of that array
    // per point, once.
    MultiFab latlon(ba2d, dm, 2, 0);
    pack_latlon(latlon, lat, lon);

    BoxArray ba_one(m_dom);
    Vector<int> pmap(1, ParallelDescriptor::IOProcessorNumber());
    DistributionMapping dm_one(pmap);
    m_ll.define(ba_one, dm_one, 2, 0, MFInfo().SetArena(The_Pinned_Arena()));
    m_ll.ParallelCopy(latlon, 0, 0, 2);

    // The copy is device-side in a GPU build; the searches are host code
    Gpu::streamSynchronize();
}

LatLonStatus
LatLonMap::locate (Real req_lat, Real req_lon, LatLonLocation& loc) const
{
    Real out[9] = {Real(0.0), Real(0.0), Real(0.0), Real(0.0), Real(0.0),
                   Real(0.0), Real(0.0), Real(0.0), Real(0.0)};

    if (ParallelDescriptor::IOProcessor()) {
        LatLonLocation l;
        const LatLonStatus status = locate_latlon_on_grid(m_ll.const_array(0), m_dom, m_problo, m_dx,
                                                          req_lat, req_lon, l);
        out[0] = l.x;        out[1] = l.y;
        out[2] = l.lat;      out[3] = l.lon;
        out[4] = l.near_lat; out[5] = l.near_lon;
        out[6] = l.cos_alpha; out[7] = l.sin_alpha;
        out[8] = static_cast<Real>(static_cast<int>(status));
    }

    ParallelDescriptor::Bcast(out, 9, ParallelDescriptor::IOProcessorNumber());

    loc.x = out[0];        loc.y = out[1];
    loc.lat = out[2];      loc.lon = out[3];
    loc.near_lat = out[4]; loc.near_lon = out[5];
    loc.cos_alpha = out[6]; loc.sin_alpha = out[7];
    return static_cast<LatLonStatus>(static_cast<int>(std::lround(out[8])));
}

void
LatLonMap::rotation_at (Real x, Real y, Real& cos_alpha, Real& sin_alpha) const
{
    Real out[2] = {Real(1.0), Real(0.0)};

    if (ParallelDescriptor::IOProcessor()) {
        const int i = std::clamp(static_cast<int>(std::floor((x - m_problo[0]) / m_dx[0])),
                                 m_dom.smallEnd(0), m_dom.bigEnd(0));
        const int j = std::clamp(static_cast<int>(std::floor((y - m_problo[1]) / m_dx[1])),
                                 m_dom.smallEnd(1), m_dom.bigEnd(1));
        grid_rotation_from_latlon(m_ll.const_array(0), m_dom, i, j, out[0], out[1]);
    }

    ParallelDescriptor::Bcast(out, 2, ParallelDescriptor::IOProcessorNumber());
    cos_alpha = out[0];
    sin_alpha = out[1];
}
