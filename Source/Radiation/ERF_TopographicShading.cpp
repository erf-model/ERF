#include <ERF_TopographicShading.H>
#include <ERF_TopographicShadingKernels.H>
#include <ERF_LatLonMap.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

using namespace amrex;

namespace {
constexpr Real missing = Real(-1.e30);

// Fill only the part covered by this ancestor. Iterating ancestors from coarse
// to fine supplies terrain outside a nested patch without inventing mountains
// at its boundary. The same piecewise-constant interpolation is used for raw
// surface radiation, consistent with ERF's existing land-forcing interpolation.
void prolong (MultiFab& dst, const MultiFab& src, const IntVect& ratio, int nc)
{
    MultiFab coarse(amrex::coarsen(dst.boxArray(), ratio), dst.DistributionMap(), nc, 0);
    coarse.setVal(missing);
    coarse.ParallelCopy(src, 0, 0, nc);
    for (MFIter mfi(dst); mfi.isValid(); ++mfi) {
        auto d = dst.array(mfi);
        auto c = coarse.const_array(mfi);
        ParallelFor(mfi.validbox(), nc, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
            const int ic = amrex::coarsen(i, ratio[0]);
            const int jc = amrex::coarsen(j, ratio[1]);
            if (c(ic,jc,0,n) != missing) { d(i,j,k,n) = c(ic,jc,0,n); }
        });
    }
}

IntVect horizontal_ratio (const Geometry& coarse, const Geometry& fine)
{
    return IntVect(fine.Domain().length(0)/coarse.Domain().length(0),
                   fine.Domain().length(1)/coarse.Domain().length(1), 1);
}
}

void TopographicShadingChoice::read (const std::string& prefix)
{
    ParmParse pp(prefix);
    pp.queryAdd("slope_rad", slope_rad);
    pp.queryAdd("topo_shading", topo_shading);
    if (!active()) { return; }
    pp.queryAdd("shadlen", shadlen);
    pp.query("fixed_solar_zenith_angle", fixed_mu);
    const bool has_azimuth = pp.query("fixed_solar_azimuth_angle", fixed_azimuth);
    pp.query("rad_cons_lat", latitude);
    pp.query("rad_cons_lon", longitude);
    pp.query("rad_freq_in_steps", radiation_interval);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(shadlen) && shadlen > 0,
        "erf.shadlen must be a positive finite horizontal distance in metres");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(fixed_mu) && fixed_mu <= 1,
        "erf.fixed_solar_zenith_angle is a cosine and must not exceed 1");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!has_azimuth || fixed_mu > 0,
        "erf.fixed_solar_azimuth_angle requires a positive fixed_solar_zenith_angle");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fixed_mu <= 0 ||
        (has_azimuth && std::isfinite(fixed_azimuth) && fixed_azimuth >= 0 && fixed_azimuth < 360),
        "Topographic radiation with fixed_solar_zenith_angle requires "
        "fixed_solar_azimuth_angle in [0,360), clockwise from geographic north");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(radiation_interval > 0,
        "Topographic radiation requires erf.rad_freq_in_steps > 0");
}

void TopographicShading::define_level (int lev, const BoxArray& ba,
                                       const DistributionMapping& dm, const Geometry& geom)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ba.isDisjoint(),
        "Topographic radiation requires horizontally disjoint column boxes; "
        "set amr.no_box_split_dir = 2");
    if (lev >= int(m_levels.size())) { m_levels.resize(lev+1); }
    m_levels[lev] = std::make_unique<Level>();
    auto& l = *m_levels[lev];
    l.geom = geom;
    l.ground.define(ba, dm, NGeom, 0);
    l.raw.define(ba, dm, NRad, 0);
    l.diagnostic.define(ba, dm, NDiag, 0);
    l.raw.setVal(missing);
    l.diagnostic.setVal(-999);
    for (int n = 0; n < NRad; ++n) {
        l.aliases.push_back(std::make_unique<MultiFab>(l.raw, make_alias, n, 1));
    }
    // Descendant caches include this level's terrain outside their own grids.
    for (int f = lev+1; f < int(m_levels.size()); ++f) {
        if (m_levels[f]) { m_levels[f]->geometry_ready = false; }
    }
}

bool TopographicShading::geometry_ready (int lev) const
{
    return lev < int(m_levels.size()) && m_levels[lev] && m_levels[lev]->geometry_ready;
}

bool TopographicShading::forcing_ready (int lev) const
{
    return lev < int(m_levels.size()) && m_levels[lev] && m_levels[lev]->forcing_ready;
}

void TopographicShading::initialize_geometry (int lev, const TopographicShadingChoice& choice,
    const MultiFab& z_nd, const MultiFab& map_x, const MultiFab& map_y,
    const MultiFab* lat, const MultiFab* lon)
{
    auto& l = *m_levels.at(lev);
    const Box dom = makeSlab(l.geom.Domain(), 2, 0);
    const auto dx = l.geom.CellSizeArray();
    l.ground.setVal(missing);
    if (lev > 0) {
        AMREX_ALWAYS_ASSERT(geometry_ready(lev-1));
        prolong(l.ground, m_levels[lev-1]->ground,
                horizontal_ratio(m_levels[lev-1]->geom, l.geom), NGeom);
    }
    l.ground.ParallelCopy(map_x, 0, MapX, 1);
    l.ground.ParallelCopy(map_y, 0, MapY, 1);
    l.ground.setVal(choice.latitude, Latitude, 1);
    l.ground.setVal(choice.longitude, Longitude, 1);
    l.ground.setVal(1, CosAlpha, 1);
    l.ground.setVal(0, SinAlpha, 1);
    if (lat && lon) {
        l.ground.ParallelCopy(*lat, 0, Latitude, 1);
        l.ground.ParallelCopy(*lon, 0, Longitude, 1);
        MultiFab c(l.ground, make_alias, CosAlpha, 1);
        MultiFab s(l.ground, make_alias, SinAlpha, 1);
        fill_grid_rotation(c, s, *lat, *lon, l.geom);
    }
    const int kbot = l.geom.Domain().smallEnd(2);
    for (MFIter mfi(l.ground); mfi.isValid(); ++mfi) {
        // An elevated nested box has no ground of its own. Use the ancestor's.
        if (z_nd.boxArray()[mfi.index()].smallEnd(2) != kbot) { continue; }
        const auto z = z_nd.const_array(mfi);
        const auto g = l.ground.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            g(i,j,k,Height) = Real(0.25)*(z(i,j,kbot) + z(i+1,j,kbot)
                                        + z(i,j+1,kbot) + z(i+1,j+1,kbot));
        });
    }
    AMREX_ALWAYS_ASSERT(l.ground.min(Height) > missing);
    AMREX_ALWAYS_ASSERT(l.ground.min(MapX) > 0 && l.ground.min(MapY) > 0);

    // One search patch per owning rank, enlarged to include every possible ray
    // and the slope stencil. ParallelCopy brings remote terrain across MPI
    // boundaries once, at initialization/regrid, rather than at every sun update.
    const auto& ba = l.ground.boxArray();
    const auto& dm = l.ground.DistributionMap();
    BoxList boxes;
    Vector<int> owners;
    const Real max_map[] = {l.ground.max(MapX), l.ground.max(MapY)};
    for (int rank = 0; rank < ParallelDescriptor::NProcs(); ++rank) {
        Box patch;
        for (int b = 0; b < int(ba.size()); ++b) {
            if (dm[b] == rank) { patch.minBox(ba[b]); }
        }
        if (!patch.ok()) { continue; }
        for (int d = 0; d < 2; ++d) {
            const Real cells = choice.topo_shading ? choice.shadlen*max_map[d]/dx[d] : 0;
            const int ng = int(std::min(Real(dom.length(d)), std::ceil(cells))) + 2;
            patch.grow(d, ng);
            // Wrapped coordinates must be addressable in this one patch.
            if (l.geom.isPeriodic(d) &&
                (patch.smallEnd(d) < dom.smallEnd(d) || patch.bigEnd(d) > dom.bigEnd(d))) {
                patch.setSmall(d, dom.smallEnd(d));
                patch.setBig(d, dom.bigEnd(d));
            }
        }
        boxes.push_back(patch & dom);
        owners.push_back(rank);
    }
    l.patch.define(BoxArray(std::move(boxes)), DistributionMapping(std::move(owners)), 1, 0);
    l.patch.setVal(missing);
    for (int c = 0; c < lev; ++c) {
        AMREX_ALWAYS_ASSERT(geometry_ready(c));
        prolong(l.patch, m_levels[c]->ground, horizontal_ratio(m_levels[c]->geom, l.geom), 1);
    }
    l.patch.ParallelCopy(l.ground, Height, 0, 1);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(l.patch.min(0) > missing,
        "Topographic search patch contains terrain not supplied by any AMR level");
    l.geometry_ready = true;
}

Vector<MultiFab*> TopographicShading::raw_fields (int lev)
{
    Vector<MultiFab*> fields;
    for (auto& f : m_levels.at(lev)->aliases) { fields.push_back(f.get()); }
    return fields;
}

MultiFab* TopographicShading::raw (int lev) { return &m_levels.at(lev)->raw; }

const MultiFab* TopographicShading::diagnostics (int lev) const
{
    return lev < int(m_levels.size()) && m_levels[lev] ? &m_levels[lev]->diagnostic : nullptr;
}

void TopographicShading::set_sun (int lev, double calday, double declin)
{
    auto& l = *m_levels.at(lev);
    l.calday = calday;
    l.declin = declin;
    l.forcing_ready = true;
}

void TopographicShading::interpolate_from_coarse (int lev)
{
    AMREX_ALWAYS_ASSERT(lev > 0 && forcing_ready(lev-1));
    auto& l = *m_levels.at(lev);
    const auto& c = *m_levels[lev-1];
    prolong(l.raw, c.raw, horizontal_ratio(c.geom, l.geom), NRad);
    set_sun(lev, c.calday, c.declin);
}

void TopographicShading::apply (int lev, const TopographicShadingChoice& choice,
                               const Vector<MultiFab*>& destinations)
{
    AMREX_ALWAYS_ASSERT(geometry_ready(lev) && forcing_ready(lev));
    auto& l = *m_levels.at(lev);
    const Box dom = l.geom.Domain();
    const int ilo = dom.smallEnd(0), ihi = dom.bigEnd(0);
    const int jlo = dom.smallEnd(1), jhi = dom.bigEnd(1);
    const bool px = l.geom.isPeriodic(0), py = l.geom.isPeriodic(1);
    const auto dx = l.geom.CellSizeArray();
    const double calday = l.calday, declin = l.declin;
    Array4<const Real> terrain;
    for (MFIter mfi(l.patch); mfi.isValid(); ++mfi) { terrain = l.patch.const_array(mfi); }
    for (MFIter mfi(l.raw); mfi.isValid(); ++mfi) {
        auto g = l.ground.const_array(mfi);
        auto r = l.raw.const_array(mfi);
        auto diag = l.diagnostic.array(mfi);
        GpuArray<Array4<Real>, NRad> out;
        for (int n = 0; n < NRad; ++n) {
            out[n] = n < int(destinations.size()) && destinations[n]
                       ? destinations[n]->array(mfi) : Array4<Real>{};
        }
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            auto s = erf_topography::sun_direction(calday, declin, g(i,j,k,Latitude), g(i,j,k,Longitude));
            if (choice.fixed_mu > 0) {
                const Real az = choice.fixed_azimuth*PI/Real(180);
                const Real h = std::sqrt(Real(1)-choice.fixed_mu*choice.fixed_mu);
                s = {h*std::sin(az), h*std::cos(az), choice.fixed_mu};
            }
            const Real east = s.x, north = s.y;
            s.x = g(i,j,k,CosAlpha)*east + g(i,j,k,SinAlpha)*north;
            s.y = -g(i,j,k,SinAlpha)*east + g(i,j,k,CosAlpha)*north;
            const Real dxp = dx[0]/g(i,j,k,MapX), dyp = dx[1]/g(i,j,k,MapY);
            const int im = px ? erf_topography::wrap(i-1,ilo,ihi-ilo+1) : amrex::max(i-1,ilo);
            const int ip = px ? erf_topography::wrap(i+1,ilo,ihi-ilo+1) : amrex::min(i+1,ihi);
            const int jm = py ? erf_topography::wrap(j-1,jlo,jhi-jlo+1) : amrex::max(j-1,jlo);
            const int jp = py ? erf_topography::wrap(j+1,jlo,jhi-jlo+1) : amrex::min(j+1,jhi);
            const Real hx = (terrain(ip,j,0)-terrain(im,j,0))/(dxp*(px ? 2 : amrex::max(1,ip-im)));
            const Real hy = (terrain(i,jp,0)-terrain(i,jm,0))/(dyp*(py ? 2 : amrex::max(1,jp-jm)));
            const bool daytime = s.z > Real(1.e-4);
            const Real inc = daytime ? erf_topography::incidence(hx,hy,s) : Real(0);
            const bool shadow = daytime && choice.topo_shading &&
                erf_topography::shadowed(i,j,terrain,ilo,ihi,jlo,jhi,px,py,dxp,dyp,s,choice.shadlen);
            // Match WRF's near-horizon cutoff; avoid dividing by a tiny zenith cosine.
            const Real factor = shadow ? Real(0) :
                ((daytime && choice.slope_rad) ? inc/s.z : Real(1));
            const Real direct = r(i,j,k,DirectVis) + r(i,j,k,DirectNir);
            const Real diffuse = r(i,j,k,DiffuseVis) + r(i,j,k,DiffuseNir);
            // Preserve the exact raw broadband flux when no adjustment is needed.
            const Real sw = factor == Real(1) ? r(i,j,k,Shortwave) :
                erf_topography::corrected_shortwave(direct,diffuse,factor);
            diag(i,j,k,Shadow) = shadow ? Real(1) : Real(0);
            diag(i,j,k,DirectFactor) = factor;
            diag(i,j,k,Incidence) = inc;
            diag(i,j,k,UnshadedSW) = r(i,j,k,Shortwave);
            diag(i,j,k,TerrainSW) = sw;
            for (int n = 0; n < NRad; ++n) {
                if (out[n]) {
                    out[n](i,j,0) = n == Shortwave ? sw :
                        ((n == DirectVis || n == DirectNir) ? factor*r(i,j,k,n) : r(i,j,k,n));
                }
            }
        });
    }
}
