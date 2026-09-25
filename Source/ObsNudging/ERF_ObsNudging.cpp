/**
 * \file ERF_ObsNudging.cpp
 */
#include "ERF_ObsNudging.H"

#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <iomanip>
#include <sstream>

#include "ERF_Constants.H"
#include "ERF_IndexDefines.H"
#include "ERF_ProbCommon.H"
#include "ERF_TerrainSurfaceSlab.H"

using namespace amrex;

namespace obs_nudging {

//
// Where a nudged quantity lives and what it needs there: the position (x, y, z)
// of face Dir (0, 1, 2 for u, v, w) or of the cell centre (Dir = -1), the
// terrain height zg under it and the density on it.  HasZnd selects the
// terrain-fitted height (from z_phys_nd) over the uniform one.
//
// This is a device function rather than part of the kernels' lambdas because
// nvcc does not let an extended device lambda first-capture a variable inside
// an if-constexpr branch; the lambdas below only call it.
//
template <int Dir, bool HasZnd>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
nudged_point (int i, int j, int k,
              const Array4<const Real>& cons, const Array4<const Real>& znd,
              const Array4<const Real>& zs,
              const GpuArray<Real, AMREX_SPACEDIM>& problo,
              const GpuArray<Real, AMREX_SPACEDIM>& dx,
              Real& x, Real& y, Real& z, Real& zg, Real& rho) noexcept
{
    if constexpr (Dir == 0) {
        x   = problo[0] + Real(i) * dx[0];
        y   = problo[1] + (Real(j) + Real(0.5)) * dx[1];
        zg  = Real(0.5) * (zs(i,j,0) + zs(i,j+1,0));
        rho = Real(0.5) * (cons(i-1,j,k,Rho_comp) + cons(i,j,k,Rho_comp));
        if constexpr (HasZnd) {
            z = Real(0.25) * (znd(i,j,k) + znd(i,j+1,k) + znd(i,j,k+1) + znd(i,j+1,k+1));
        } else {
            z = problo[2] + (Real(k) + Real(0.5)) * dx[2];
        }
    } else if constexpr (Dir == 1) {
        x   = problo[0] + (Real(i) + Real(0.5)) * dx[0];
        y   = problo[1] + Real(j) * dx[1];
        zg  = Real(0.5) * (zs(i,j,0) + zs(i+1,j,0));
        rho = Real(0.5) * (cons(i,j-1,k,Rho_comp) + cons(i,j,k,Rho_comp));
        if constexpr (HasZnd) {
            z = Real(0.25) * (znd(i,j,k) + znd(i+1,j,k) + znd(i,j,k+1) + znd(i+1,j,k+1));
        } else {
            z = problo[2] + (Real(k) + Real(0.5)) * dx[2];
        }
    } else if constexpr (Dir == 2) {
        x   = problo[0] + (Real(i) + Real(0.5)) * dx[0];
        y   = problo[1] + (Real(j) + Real(0.5)) * dx[1];
        zg  = Real(0.25) * (zs(i,j,0) + zs(i+1,j,0) + zs(i,j+1,0) + zs(i+1,j+1,0));
        rho = Real(0.5) * (cons(i,j,k-1,Rho_comp) + cons(i,j,k,Rho_comp));
        if constexpr (HasZnd) {
            z = Real(0.25) * (znd(i,j,k) + znd(i+1,j,k) + znd(i,j+1,k) + znd(i+1,j+1,k));
        } else {
            z = problo[2] + Real(k) * dx[2];
        }
    } else {
        x   = problo[0] + (Real(i) + Real(0.5)) * dx[0];
        y   = problo[1] + (Real(j) + Real(0.5)) * dx[1];
        zg  = Real(0.25) * (zs(i,j,0) + zs(i+1,j,0) + zs(i,j+1,0) + zs(i+1,j+1,0));
        rho = cons(i,j,k,Rho_comp);
        if constexpr (HasZnd) {
            z = Real(0.125) * (znd(i,j,k  ) + znd(i+1,j,k  ) + znd(i,j+1,k  ) + znd(i+1,j+1,k  ) +
                               znd(i,j,k+1) + znd(i+1,j,k+1) + znd(i,j+1,k+1) + znd(i+1,j+1,k+1));
        } else {
            z = problo[2] + (Real(k) + Real(0.5)) * dx[2];
        }
    }
}

//
// The kernels.  Free functions in a named namespace, capturing nothing but
// values, so that the device lambdas are legal in every GPU build.
//
template <int Dir, bool HasZnd>
void
add_face_source (const Box& bx, const Array4<Real>& src, const Array4<const Real>& vel,
                 const Array4<const Real>& cons, const Array4<const Real>& znd,
                 const Array4<const Real>& zs, const Array4<const Real>& blank, bool has_blank,
                 const GpuArray<Real, AMREX_SPACEDIM> problo,
                 const GpuArray<Real, AMREX_SPACEDIM> dx,
                 int klo, int khi, const ObsNudgingView v, int comp)
{
    // w on the bottom and top of the domain is set by the boundary conditions
    const bool skip_walls = (Dir == 2);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (skip_walls && (k <= klo || k > khi)) { return; }

        Real x, y, z, zg, rho;
        nudged_point<Dir, HasZnd>(i, j, k, cons, znd, zs, problo, dx, x, y, z, zg, rho);

        Real tend = obs_nudging_tendency(v, comp, x, y, z, zg, vel(i,j,k));
        if (has_blank) { tend *= (Real(1.0) - blank(i,j,k)); }
        src(i,j,k) += rho * tend;
    });
}

template <bool HasZnd>
void
add_cell_theta_source (const Box& bx, const Array4<Real>& src,
                       const Array4<const Real>& cons, const Array4<const Real>& znd,
                       const Array4<const Real>& zs, const Array4<const Real>& blank, bool has_blank,
                       const GpuArray<Real, AMREX_SPACEDIM> problo,
                       const GpuArray<Real, AMREX_SPACEDIM> dx,
                       const ObsNudgingView v)
{
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x, y, z, zg, rho;
        nudged_point<-1, HasZnd>(i, j, k, cons, znd, zs, problo, dx, x, y, z, zg, rho);

        const Real theta = cons(i,j,k,RhoTheta_comp) / rho;

        Real tend = obs_nudging_tendency(v, Theta, x, y, z, zg, theta);
        if (has_blank) { tend *= (Real(1.0) - blank(i,j,k)); }
        src(i,j,k,RhoTheta_comp) += rho * tend;
    });
}

} // namespace obs_nudging

namespace {

void
read_station_file (ObsNudgingStation& st, Real missing_value)
{
    if (!FileExists(st.file)) {
        Abort("erf.obs_nudging." + st.name + ".file: cannot open '" + st.file + "'");
    }
    Vector<char> buffer;
    ParallelDescriptor::ReadAndBcastFile(st.file, buffer);
    std::istringstream is(std::string(buffer.dataPtr()), std::istringstream::in);

    std::string error;
    if (!obs_nudging::parse_station_series(is, st.file, missing_value, st.series, error)) {
        Abort("erf.obs_nudging." + st.name + ".file: " + error);
    }
}

} // namespace

ObsNudging::ObsNudging (TerrainType terrain_type, ProblemBase* prob,
                        bool use_datetime, double start_time)
    : m_terrain_type(terrain_type),
      m_prob(prob)
{
    // Also refused by SolverChoice as soon as the inputs are read
    if (terrain_type == TerrainType::EB) {
        Abort("erf.nudging_from_observations does not support erf.terrain_type = EB; "
              "use a terrain-fitted mesh or immersed forcing");
    }

    ParmParse pp("erf.obs_nudging");

    if (!pp.query("tau", m_tau)) {
        Abort("erf.nudging_from_observations needs erf.obs_nudging.tau, the relaxation time "
              "scale in seconds; there is no default");
    }
    if (!(m_tau > Real(0.0))) {
        Abort("erf.obs_nudging.tau must be positive");
    }

    pp.query("horizontal_radius", m_horizontal_radius);
    pp.query("vertical_radius", m_vertical_radius);
    pp.query("cutoff", m_cutoff);
    pp.query("sigma_factor", m_sigma_factor);
    pp.query("missing_value", m_missing_value);
    pp.query("nudge_wind", m_nudge_wind);
    pp.query("nudge_w", m_nudge_w);
    pp.query("nudge_theta", m_nudge_theta);

    if (!(m_horizontal_radius > Real(0.0))) {
        Abort("erf.obs_nudging.horizontal_radius must be positive");
    }
    if (!(m_vertical_radius > Real(0.0))) {
        Abort("erf.obs_nudging.vertical_radius must be positive");
    }
    if (!(m_cutoff > Real(0.0))) {
        Abort("erf.obs_nudging.cutoff must be positive");
    }
    if (!(m_sigma_factor >= Real(0.0))) {
        Abort("erf.obs_nudging.sigma_factor must not be negative");
    }
    if (!m_nudge_wind && !m_nudge_w && !m_nudge_theta) {
        Abort("erf.obs_nudging: nudge_wind, nudge_w and nudge_theta are all false, so "
              "nothing would be nudged");
    }

    std::string time_type = "elapsed";
    pp.query("time_type", time_type);
    if (time_type == "elapsed") {
        m_time_origin = 0.0;
    } else if (time_type == "epoch") {
        if (!use_datetime) {
            Abort("erf.obs_nudging.time_type = epoch needs the run to know its start date "
                  "(start_datetime, or a WRF or metgrid initialization)");
        }
        m_time_origin = start_time;
    } else {
        Abort("erf.obs_nudging.time_type must be elapsed or epoch, not '" + time_type + "'");
    }

    const int nst = pp.countval("stations");
    if (nst == 0) {
        Abort("erf.nudging_from_observations needs erf.obs_nudging.stations, the names of "
              "the stations");
    }
    Vector<std::string> names;
    pp.getarr("stations", names);

    for (const auto& name : names)
    {
        for (const auto& st : m_stations) {
            if (st.name == name) {
                Abort("erf.obs_nudging.stations names '" + name + "' twice");
            }
        }

        ObsNudgingStation st;
        st.name = name;

        ParmParse pps("erf.obs_nudging." + name);

        if (!pps.query("file", st.file)) {
            Abort("erf.obs_nudging." + name + ".file is required");
        }

        const bool has_lat = pps.contains("lat");
        const bool has_lon = pps.contains("long") || pps.contains("lon");
        const bool has_x   = pps.contains("x");
        const bool has_y   = pps.contains("y");
        if ((has_lat || has_lon) && (has_x || has_y)) {
            Abort("erf.obs_nudging." + name + ": give either lat/long or x/y, not both");
        }
        if (has_lat || has_lon) {
            if (!has_lat || !has_lon) {
                Abort("erf.obs_nudging." + name + ": lat and long must both be given");
            }
            st.use_latlon = true;
            pps.get("lat", st.req_lat);
            if (!pps.query("long", st.req_lon)) { pps.get("lon", st.req_lon); }
            if (std::abs(st.req_lat) > Real(90.0) || std::abs(st.req_lon) > Real(360.0)) {
                Abort("erf.obs_nudging." + name + ": lat must be in [-90, 90] and long in "
                      "[-360, 360] degrees");
            }
        } else {
            if (!has_x || !has_y) {
                Abort("erf.obs_nudging." + name + ": give lat and long, or x and y in domain "
                      "coordinates");
            }
            pps.get("x", st.x);
            pps.get("y", st.y);
        }

        std::string frame = "earth";
        pps.query("wind_frame", frame);
        if (frame == "earth") {
            st.earth_frame = true;
        } else if (frame == "grid") {
            st.earth_frame = false;
        } else {
            Abort("erf.obs_nudging." + name + ".wind_frame must be earth or grid, not '" + frame + "'");
        }

        std::string height_ref = "agl";
        pps.query("height_ref", height_ref);
        if (height_ref == "agl") {
            st.agl = true;
        } else if (height_ref == "msl") {
            st.agl = false;
        } else {
            Abort("erf.obs_nudging." + name + ".height_ref must be agl or msl, not '" +
                  height_ref + "'");
        }

        read_station_file(st, m_missing_value);

        const auto& has = st.series.has;
        const bool used = (m_nudge_wind && has[obs_nudging::U]) ||
                          (m_nudge_w && has[obs_nudging::W]) ||
                          (m_nudge_theta && has[obs_nudging::Theta]);
        if (!used) {
            // A lidar without temperature in a run that nudges only theta, say:
            // nothing to do at this station, which is worth saying but not fatal
            Warning("erf.obs_nudging." + name + ": the file " + st.file + " has none of the "
                    "quantities this run nudges (see nudge_wind, nudge_w, nudge_theta), so the "
                    "station is not used");
            continue;
        }

        m_stations.push_back(std::move(st));
    }

    if (m_stations.empty()) {
        Abort("erf.obs_nudging: no station measures any of the quantities this run nudges "
              "(see nudge_wind, nudge_w, nudge_theta)");
    }

    m_d_x.resize(m_stations.size());
    m_d_y.resize(m_stations.size());
    m_d_agl.resize(m_stations.size());
    m_d_offset.resize(obs_nudging::NComp * m_stations.size());
    m_d_count.resize(obs_nudging::NComp * m_stations.size());
}

bool
ObsNudging::wants_latlon () const
{
    for (const auto& st : m_stations) {
        if (st.use_latlon || st.earth_frame) { return true; }
    }
    return false;
}

void
ObsNudging::resolve_positions (const Geometry& geom0, const LatLonMap* latlon)
{
    const auto problo = geom0.ProbLoArray();
    const auto probhi = geom0.ProbHiArray();

    for (auto& st : m_stations)
    {
        if (st.use_latlon) {
            if (latlon == nullptr) {
                Abort("erf.obs_nudging." + st.name + ": lat/long needs a run with latitude/"
                      "longitude arrays (a WRF or metgrid initialization, or a restart from "
                      "one).  Place the station with x and y in domain coordinates instead");
            }
            LatLonLocation where;
            const LatLonStatus status = latlon->locate(st.req_lat, st.req_lon, where);
            if (status == LatLonStatus::Degenerate) {
                Abort("erf.obs_nudging." + st.name + ": the latitude/longitude arrays are "
                      "degenerate near the requested point, so it cannot be inverted");
            }
            if (status == LatLonStatus::TooFar) {
                Abort("erf.obs_nudging." + st.name + ": lat=" + std::to_string(st.req_lat) +
                      " long=" + std::to_string(st.req_lon) + " is more than one cell from the "
                      "nearest grid point (lat=" + std::to_string(where.near_lat) + " long=" +
                      std::to_string(where.near_lon) + "), so it is outside the domain");
            }
            st.x = where.x;
            st.y = where.y;
            st.got_lat = where.lat;
            st.got_lon = where.lon;
            if (st.earth_frame) {
                st.cos_alpha = where.cos_alpha;
                st.sin_alpha = where.sin_alpha;
            }
        } else if (st.earth_frame && latlon != nullptr) {
            latlon->rotation_at(st.x, st.y, st.cos_alpha, st.sin_alpha);
        }
        // Otherwise the grid is aligned with the frame of the file: no rotation

        if (st.x < problo[0] || st.x > probhi[0] || st.y < problo[1] || st.y > probhi[1]) {
            Abort("erf.obs_nudging." + st.name + ": x=" + std::to_string(st.x) + " y=" +
                  std::to_string(st.y) + " is outside the problem domain");
        }
    }
}

void
ObsNudging::print_summary () const
{
    Print() << "Observation nudging: " << m_stations.size() << " station(s), tau = " << m_tau
            << " s, R_h = " << m_horizontal_radius << " m, R_z = " << m_vertical_radius
            << " m, cutoff = " << m_cutoff << ", sigma_factor = " << m_sigma_factor << "\n";
    for (const auto& st : m_stations)
    {
        const auto& s = st.series;
        std::ostringstream os;
        os << "  " << st.name << ": x = " << st.x << ", y = " << st.y;
        if (st.use_latlon) {
            os << " (lat " << st.req_lat << ", long " << st.req_lon << " -> "
               << st.got_lat << ", " << st.got_lon << ")";
        }
        const Real alpha = std::atan2(st.sin_alpha, st.cos_alpha) * Real(180.0) / PI;
        os << ", " << s.heights.size() << " height(s) " << s.heights.front() << " to "
           << s.heights.back() << " m " << (st.agl ? "above the terrain" : "above z = 0")
           << ", " << s.records.size() << " time(s) " << s.records.front().time << " to "
           << s.records.back().time << " s, wind frame " << (st.earth_frame ? "earth" : "grid")
           << " (rotated " << alpha << " deg), nudging";
        if (m_nudge_wind  && s.has[obs_nudging::U])     { os << " u v"; }
        if (m_nudge_w     && s.has[obs_nudging::W])     { os << " w"; }
        if (m_nudge_theta && s.has[obs_nudging::Theta]) { os << " theta"; }
        Print() << os.str() << "\n";
    }
}

void
ObsNudging::update_targets (double time)
{
    if (m_have_targets && time == m_targets_time) { return; }

    const int nst = static_cast<int>(m_stations.size());
    const bool nudge[obs_nudging::NComp] = {m_nudge_wind, m_nudge_wind, m_nudge_w, m_nudge_theta};

    Vector<Real> h_x(nst), h_y(nst);
    Vector<int>  h_agl(nst);
    Vector<int>  h_offset(obs_nudging::NComp * nst, 0), h_count(obs_nudging::NComp * nst, 0);
    Vector<Real> h_z, h_mean, h_sigma;

    for (int s = 0; s < nst; ++s)
    {
        const auto& st = m_stations[s];
        h_x[s]   = st.x;
        h_y[s]   = st.y;
        h_agl[s] = st.agl ? 1 : 0;

        std::array<obs_nudging::Profile, obs_nudging::NComp> prof;
        obs_nudging::profiles_at_time(st.series, time + m_time_origin,
                                      st.cos_alpha, st.sin_alpha, prof);

        for (int c = 0; c < obs_nudging::NComp; ++c) {
            if (!nudge[c]) { continue; }
            const int idx = c*nst + s;
            h_offset[idx] = static_cast<int>(h_z.size());
            h_count[idx]  = static_cast<int>(prof[c].z.size());
            h_z.insert(h_z.end(), prof[c].z.begin(), prof[c].z.end());
            h_mean.insert(h_mean.end(), prof[c].mean.begin(), prof[c].mean.end());
            h_sigma.insert(h_sigma.end(), prof[c].sigma.begin(), prof[c].sigma.end());
        }
    }

    // A zero-length device array has no data pointer to hand out; keep one
    // element so that the view is always valid, and let the counts say it is unused
    if (h_z.empty()) { h_z.push_back(Real(0.0)); h_mean.push_back(Real(0.0)); h_sigma.push_back(Real(0.0)); }

    m_d_z.resize(h_z.size());
    m_d_mean.resize(h_mean.size());
    m_d_sigma.resize(h_sigma.size());

    Gpu::copy(Gpu::hostToDevice, h_x.begin(), h_x.end(), m_d_x.begin());
    Gpu::copy(Gpu::hostToDevice, h_y.begin(), h_y.end(), m_d_y.begin());
    Gpu::copy(Gpu::hostToDevice, h_agl.begin(), h_agl.end(), m_d_agl.begin());
    Gpu::copy(Gpu::hostToDevice, h_offset.begin(), h_offset.end(), m_d_offset.begin());
    Gpu::copy(Gpu::hostToDevice, h_count.begin(), h_count.end(), m_d_count.begin());
    Gpu::copy(Gpu::hostToDevice, h_z.begin(), h_z.end(), m_d_z.begin());
    Gpu::copy(Gpu::hostToDevice, h_mean.begin(), h_mean.end(), m_d_mean.begin());
    Gpu::copy(Gpu::hostToDevice, h_sigma.begin(), h_sigma.end(), m_d_sigma.begin());
    Gpu::streamSynchronize();

    m_have_targets = true;
    m_targets_time = time;
}

obs_nudging::ObsNudgingView
ObsNudging::view (Real dt, const Geometry& geom) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dt > Real(0.0), "ObsNudging: the time step must be positive");
    AMREX_ALWAYS_ASSERT(m_have_targets);

    obs_nudging::ObsNudgingView v;
    v.nstations    = static_cast<int>(m_stations.size());
    v.x            = m_d_x.data();
    v.y            = m_d_y.data();
    v.agl          = m_d_agl.data();
    v.offset       = m_d_offset.data();
    v.count        = m_d_count.data();
    v.z            = m_d_z.data();
    v.mean         = m_d_mean.data();
    v.sigma        = m_d_sigma.data();
    v.inv_rh2      = Real(1.0) / (m_horizontal_radius * m_horizontal_radius);
    v.inv_rz2      = Real(1.0) / (m_vertical_radius * m_vertical_radius);
    v.qmax         = m_cutoff * m_cutoff;
    v.sigma_factor = m_sigma_factor;
    v.inv_tau      = Real(1.0) / m_tau;
    v.max_rate     = Real(1.0) / dt;
    v.period_x     = geom.isPeriodic(0) ? static_cast<Real>(geom.ProbLength(0)) : Real(0.0);
    v.period_y     = geom.isPeriodic(1) ? static_cast<Real>(geom.ProbLength(1)) : Real(0.0);
    return v;
}

const MultiFab&
ObsNudging::surface (int lev, const Geometry& geom, const MultiFab& cons,
                     const MultiFab* z_phys_nd, double time)
{
    if (lev >= static_cast<int>(m_zsurf.size())) {
        m_zsurf.resize(lev+1);
        m_zsurf_ba.resize(lev+1);
        m_zsurf_dm.resize(lev+1);
        m_zsurf_time.resize(lev+1, 0.0);   // read only when m_zsurf[lev] is set
    }

    const bool fitted = (m_terrain_type == TerrainType::StaticFittedMesh ||
                         m_terrain_type == TerrainType::MovingFittedMesh);
    const bool moving = (m_terrain_type == TerrainType::MovingFittedMesh);

    // On a moving fitted mesh the terrain is rebuilt for each RK stage from the
    // stage's start time (update_terrain_stage), so it is a function of the
    // grids and of that time; on every other mesh it depends on the grids
    // alone.  Both source calls of a stage pass the same time, so the slab is
    // built once per stage rather than once per call: it costs a MultiFab, a
    // ParallelCopy and a global reduction that would otherwise be paid twice.
    if (m_zsurf[lev] &&
        m_zsurf_ba[lev] == cons.boxArray() && m_zsurf_dm[lev] == cons.DistributionMap() &&
        (!moving || m_zsurf_time[lev] == time)) {
        return *m_zsurf[lev];
    }

    // The nodes of the bottom slab under every grid of the level, laid out
    // like the level so that an MFIter over the state indexes it too
    m_zsurf[lev] = std::make_unique<MultiFab>(bottom_node_slab(cons.boxArray(), geom),
                                              cons.DistributionMap(), 1, 0);
    MultiFab& zs = *m_zsurf[lev];

    if (fitted)
    {
        // The bottom of the mesh is the terrain.  The surface nodes live only
        // in the boxes that touch the bottom of the domain, so gather them.
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(z_phys_nd != nullptr,
            "ObsNudging: a terrain-fitted mesh without z_phys_nd");
        zs.setVal(bogus_large_value);
        zs.ParallelCopy(*z_phys_nd, 0, 0, 1, IntVect(0), IntVect(0));
        const Real zmax = zs.max(0);
        if (zmax > Real(0.5) * bogus_large_value) {
            Abort("erf.nudging_from_observations: the grids of level " + std::to_string(lev) +
                  " do not reach the bottom of the domain everywhere, so the terrain under "
                  "them is not known on that level.  With a terrain-fitted mesh every level "
                  "must reach the ground where it is refined");
        }
    }
    else if (m_terrain_type == TerrainType::ImmersedForcing)
    {
        // The terrain surface the immersed boundary is built from, at this
        // level's resolution
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_prob != nullptr, "ObsNudging: no problem object");
        fill_terrain_surface_slab(zs, geom, *m_prob, time);
    }
    else
    {
        // Flat: the ground is the bottom of the domain
        zs.setVal(geom.ProbLo(2));
    }

    m_zsurf_ba[lev] = cons.boxArray();
    m_zsurf_dm[lev] = cons.DistributionMap();
    m_zsurf_time[lev] = time;
    return zs;
}

void
ObsNudging::prepare_level (int lev, const Geometry& geom, const MultiFab& cons,
                           const MultiFab* z_phys_nd, double time)
{
    surface(lev, geom, cons, z_phys_nd, time);
}

void
ObsNudging::add_theta_source (int lev, double time, Real dt, const Geometry& geom,
                              const MultiFab& cons, MultiFab& cc_src,
                              const MultiFab* z_phys_nd, const MultiFab* blank)
{
    if (!m_nudge_theta) { return; }

    update_targets(time);
    const obs_nudging::ObsNudgingView v = view(dt, geom);
    const MultiFab& zs = surface(lev, geom, cons, z_phys_nd, time);

    const auto problo = geom.ProbLoArray();
    const auto dx     = geom.CellSizeArray();
    const bool has_blank = (blank != nullptr);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Array4<Real>&       src = cc_src.array(mfi);
        const Array4<const Real>& c   = cons.const_array(mfi);
        const Array4<const Real>& s   = zs.const_array(mfi);
        const Array4<const Real>  b   = has_blank ? blank->const_array(mfi) : Array4<const Real>{};
        if (z_phys_nd) {
            obs_nudging::add_cell_theta_source<true>(bx, src, c, z_phys_nd->const_array(mfi), s, b,
                                                     has_blank, problo, dx, v);
        } else {
            obs_nudging::add_cell_theta_source<false>(bx, src, c, Array4<const Real>{}, s, b,
                                                      has_blank, problo, dx, v);
        }
    }
}

void
ObsNudging::add_momentum_sources (int lev, double time, Real dt, const Geometry& geom,
                                  const MultiFab& cons,
                                  const MultiFab& xvel, const MultiFab& yvel, const MultiFab& zvel,
                                  MultiFab& xmom_src, MultiFab& ymom_src, MultiFab& zmom_src,
                                  const MultiFab* z_phys_nd,
                                  const MultiFab* blank_x, const MultiFab* blank_y,
                                  const MultiFab* blank_z)
{
    if (!m_nudge_wind && !m_nudge_w) { return; }

    update_targets(time);
    const obs_nudging::ObsNudgingView v = view(dt, geom);
    const MultiFab& zs = surface(lev, geom, cons, z_phys_nd, time);

    const auto problo = geom.ProbLoArray();
    const auto dx     = geom.CellSizeArray();
    const int  klo    = geom.Domain().smallEnd(2);
    const int  khi    = geom.Domain().bigEnd(2);
    const bool has_blank = (blank_x != nullptr);
    AMREX_ALWAYS_ASSERT((blank_x != nullptr) == (blank_y != nullptr) &&
                        (blank_x != nullptr) == (blank_z != nullptr));

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Array4<const Real>& c   = cons.const_array(mfi);
        const Array4<const Real>& s   = zs.const_array(mfi);
        const Array4<const Real>  znd = z_phys_nd ? z_phys_nd->const_array(mfi) : Array4<const Real>{};
        const Array4<const Real>  bx_ = has_blank ? blank_x->const_array(mfi) : Array4<const Real>{};
        const Array4<const Real>  by_ = has_blank ? blank_y->const_array(mfi) : Array4<const Real>{};
        const Array4<const Real>  bz_ = has_blank ? blank_z->const_array(mfi) : Array4<const Real>{};

        if (m_nudge_wind) {
            const Box tbx = mfi.nodaltilebox(0);
            const Box tby = mfi.nodaltilebox(1);
            if (z_phys_nd) {
                obs_nudging::add_face_source<0,true>(tbx, xmom_src.array(mfi), xvel.const_array(mfi), c, znd,
                                                     s, bx_, has_blank, problo, dx, klo, khi, v, obs_nudging::U);
                obs_nudging::add_face_source<1,true>(tby, ymom_src.array(mfi), yvel.const_array(mfi), c, znd,
                                                     s, by_, has_blank, problo, dx, klo, khi, v, obs_nudging::V);
            } else {
                obs_nudging::add_face_source<0,false>(tbx, xmom_src.array(mfi), xvel.const_array(mfi), c, znd,
                                                      s, bx_, has_blank, problo, dx, klo, khi, v, obs_nudging::U);
                obs_nudging::add_face_source<1,false>(tby, ymom_src.array(mfi), yvel.const_array(mfi), c, znd,
                                                      s, by_, has_blank, problo, dx, klo, khi, v, obs_nudging::V);
            }
        }
        if (m_nudge_w) {
            const Box tbz = mfi.nodaltilebox(2);
            if (z_phys_nd) {
                obs_nudging::add_face_source<2,true>(tbz, zmom_src.array(mfi), zvel.const_array(mfi), c, znd,
                                                     s, bz_, has_blank, problo, dx, klo, khi, v, obs_nudging::W);
            } else {
                obs_nudging::add_face_source<2,false>(tbz, zmom_src.array(mfi), zvel.const_array(mfi), c, znd,
                                                      s, bz_, has_blank, problo, dx, klo, khi, v, obs_nudging::W);
            }
        }
    }
}
