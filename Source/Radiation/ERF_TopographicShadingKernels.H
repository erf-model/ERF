#ifndef ERF_TOPOGRAPHIC_SHADING_KERNELS_H_
#define ERF_TOPOGRAPHIC_SHADING_KERNELS_H_

#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <ERF_Constants.H>
#include <cmath>
#include <limits>

namespace erf_topography {
using amrex::Real;

struct Sun { Real x, y, z; };

// Geographic east/north/up, using the same orbital hour-angle convention as
// orbital_cos_zenith_instant. Return a vector rather than recovering an azimuth
// from the radiation solver's possibly interval-averaged zenith cosine.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Sun sun_direction (double calday, double declin, Real lat_deg, Real lon_deg)
{
    const double lat = double(lat_deg) * double(PI) / 180.0;
    const double h = 2.0 * double(PI) * (calday - std::floor(calday))
                   + double(lon_deg) * double(PI) / 180.0 - double(PI);
    return {Real(-std::cos(declin) * std::sin(h)),
            Real(std::cos(lat)*std::sin(declin) - std::sin(lat)*std::cos(declin)*std::cos(h)),
            Real(std::sin(lat)*std::sin(declin) + std::cos(lat)*std::cos(declin)*std::cos(h))};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real incidence (Real hx, Real hy, Sun s)
{
    return amrex::max(Real(0), (s.z - hx*s.x - hy*s.y) / std::sqrt(Real(1) + hx*hx + hy*hy));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int wrap (int i, int lo, int n) { return lo + ((i-lo) % n + n) % n; }

// Scan grid-line crossings toward the sun, interpolating the resolved terrain
// between mass points (WRF's toposhad construction). dx/dy are physical lengths
// at the target, including its map factors; the ray uses that local tangent plane.
// The search cutoff is horizontal distance, not a count of grid rows.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool shadowed (int i, int j, const amrex::Array4<const Real>& terrain,
               int ilo, int ihi, int jlo, int jhi, bool periodic_x, bool periodic_y,
               Real dx, Real dy, Sun s, Real search_length)
{
    const Real horizontal = std::sqrt(s.x*s.x + s.y*s.y);
    if (s.z <= Real(1.e-4) || horizontal <= Real(1.e-8)) { return false; }
    // Trigonometric roundoff at cardinal azimuths must not send a ray out
    // through a domain edge to which it is parallel (e.g. cos(90 degrees)).
    constexpr Real angular_eps = Real(8)*std::numeric_limits<Real>::epsilon();
    const Real di = std::abs(s.x) < angular_eps*horizontal ? Real(0) : s.x / (horizontal * dx);
    const Real dj = std::abs(s.y) < angular_eps*horizontal ? Real(0) : s.y / (horizontal * dy);
    const Real ds = Real(1) / amrex::max(std::abs(di), std::abs(dj));
    const Real z = terrain(i,j,0,0);
    for (Real r = ds; r <= search_length; r += ds) {
        const Real x = Real(i) + r*di, y = Real(j) + r*dj;
        if ((!periodic_x && (x < ilo || x > ihi)) ||
            (!periodic_y && (y < jlo || y > jhi))) { break; }
        int ix = int(std::floor(x)), jy = int(std::floor(y));
        const Real wx = x-Real(ix), wy = y-Real(jy);
        int ip = ix+1, jp = jy+1;
        if (periodic_x) { ix = wrap(ix,ilo,ihi-ilo+1); ip = wrap(ip,ilo,ihi-ilo+1); }
        else { ip = amrex::min(ip,ihi); }
        if (periodic_y) { jy = wrap(jy,jlo,jhi-jlo+1); jp = wrap(jp,jlo,jhi-jlo+1); }
        else { jp = amrex::min(jp,jhi); }
        AMREX_ASSERT(terrain.contains(ix,jy,0) && terrain.contains(ip,jp,0));
        const Real h = (Real(1)-wy)*((Real(1)-wx)*terrain(ix,jy,0,0) + wx*terrain(ip,jy,0,0))
                     + wy*((Real(1)-wx)*terrain(ix,jp,0,0) + wx*terrain(ip,jp,0,0));
        // A small vertical tolerance avoids self-shadowing a numerically flat
        // or exactly sun-parallel surface through interpolation roundoff.
        if (h > z + r*s.z/horizontal + Real(1.e-4)) { return true; }
    }
    return false;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real corrected_shortwave (Real direct, Real diffuse, Real factor)
{
    return diffuse + direct*factor;
}
} // namespace erf_topography
#endif
