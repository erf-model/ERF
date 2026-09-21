/**
 * \file ERF_StationSampler.cpp
 */

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <utility>

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

#include <ERF.H>
#include <ERF_StationSampler.H>
#include <ERF_Plotfile2DCatalog.H>
#include <ERF_Plotfile2DSampledLevel.H>
#include <ERF_Plotfile2DUtils.H>
#include <ERF_Plotfile2DWaterPath.H>
#include <ERF_Constants.H>
#include <ERF_EpochTime.H>
#include <ERF_Utils.H>

using namespace amrex;

namespace {

constexpr int datwidth      = 20;
constexpr int datprecision  = 10;
constexpr int timeprecision = 13;

std::string
station_key_error (const std::string& station, const std::string& key, const std::string& why)
{
    return "Station '" + station + "': erf." + station + "." + key + " " + why;
}

} // namespace

//
// Parse erf.station_names and the per-station keys.  Nothing here needs to know
// which variables the run can produce; that check is ERF::init_stations, which
// has the catalogs.
//
StationSampler::StationSampler (const std::string& pp_prefix)
{
    ParmParse pp(pp_prefix);

    Vector<std::string> names;
    const int nstation = pp.countval("station_names");
    if (nstation > 0) {
        pp.getarr("station_names", names, 0, nstation);
    }

    pp.queryAdd("station_buffer_steps", m_buffer_steps);
    if (m_buffer_steps < 1) { m_buffer_steps = 1; }

    pp.queryAdd("station_output_dir", m_dir);

    for (const auto& name : names)
    {
        Station station;
        station.name = name;

        ParmParse pps(pp_prefix + "." + name);

        const int nfield = pps.countval("field");
        if (nfield <= 0) {
            Abort(station_key_error(name, "field", "must name at least one variable"));
        }
        Vector<std::string> fields;
        pps.getarr("field", fields, 0, nfield);
        for (const auto& field : fields) {
            StationVar var;
            var.name = field;
            station.vars.push_back(var);
        }

        // A station is placed either in lat/lon or in domain coordinates, never both.
        const std::string lon_key = pps.contains("long") ? "long" : "lon";
        const bool has_lat = pps.contains("lat");
        const bool has_lon = pps.contains("long") || pps.contains("lon");
        const bool has_x   = pps.contains("x");
        const bool has_y   = pps.contains("y");

        if ((has_lat || has_lon) && (has_x || has_y)) {
            Abort("Station '" + name + "': specify either lat/long or x/y, not both");
        }

        if (has_lat || has_lon) {
            if (!has_lat || !has_lon) {
                Abort("Station '" + name + "': lat and long must both be given");
            }
            const int nlat = pps.countval("lat");
            const int nlon = pps.countval(lon_key.c_str());
            if (nlat != nlon) {
                Abort("Station '" + name + "': lat has " + std::to_string(nlat) +
                      " values but " + lon_key + " has " + std::to_string(nlon) +
                      "; they are paired, so the counts must match");
            }
            Vector<Real> lat, lon;
            pps.getarr("lat", lat, 0, nlat);
            pps.getarr(lon_key.c_str(), lon, 0, nlon);
            for (int i = 0; i < nlat; ++i) {
                StationLoc loc;
                loc.use_latlon = true;
                loc.req_lat = lat[i];
                loc.req_lon = lon[i];
                station.locs.push_back(loc);
            }
        } else {
            if (!has_x || !has_y) {
                Abort("Station '" + name + "': give either lat/long, or x and y in domain coordinates");
            }
            const int nx = pps.countval("x");
            const int ny = pps.countval("y");
            if (nx != ny) {
                Abort("Station '" + name + "': x has " + std::to_string(nx) +
                      " values but y has " + std::to_string(ny) +
                      "; they are paired, so the counts must match");
            }
            Vector<Real> xv, yv;
            pps.getarr("x", xv, 0, nx);
            pps.getarr("y", yv, 0, ny);
            for (int i = 0; i < nx; ++i) {
                StationLoc loc;
                loc.use_latlon = false;
                loc.x = xv[i];
                loc.y = yv[i];
                station.locs.push_back(loc);
            }
        }

        const int nheight = pps.countval("height");
        if (nheight > 0) {
            pps.getarr("height", station.heights, 0, nheight);
        }

        m_stations.push_back(std::move(station));
    }
}

//
// Column order within one station's file is location-major: for each location,
// the 2D variables first (they have no height), then, for each height in the
// order requested, the 3D variables in the order requested.
//
void
StationSampler::buildColumns ()
{
    m_ncolumns = 0;
    for (auto& station : m_stations) {
        station.columns.clear();
        for (int il = 0; il < static_cast<int>(station.locs.size()); ++il) {
            station.locs[il].col_begin = static_cast<int>(station.columns.size());
            for (int iv = 0; iv < static_cast<int>(station.vars.size()); ++iv) {
                if (station.vars[iv].is_2d) { station.columns.push_back(m_ncolumns++); }
            }
            for (int ih = 0; ih < static_cast<int>(station.heights.size()); ++ih) {
                for (int iv = 0; iv < static_cast<int>(station.vars.size()); ++iv) {
                    if (!station.vars[iv].is_2d) { station.columns.push_back(m_ncolumns++); }
                }
            }
        }
    }
}

void
StationSampler::appendRow (Real time, Real epoch_time, const Vector<Real>& values)
{
    AMREX_ALWAYS_ASSERT(static_cast<int>(values.size()) == m_ncolumns);

    if (ParallelDescriptor::IOProcessor()) {
        m_times.push_back(time);
        m_epoch_times.push_back(epoch_time);
        m_rows.insert(m_rows.end(), values.begin(), values.end());

        if (static_cast<int>(m_times.size()) >= m_buffer_steps) { flush(); }
    }
}

void
StationSampler::writeHeader (const Station& station, std::ostream& os) const
{
    // Coordinates are written with enough digits that a reader can tell two
    // neighbouring stations apart
    os << std::setprecision(datprecision);
    os << "# ERF station time series\n";
    os << "# station: " << station.name << "\n";

    int col = 1;
    os << "# column " << col++ << ": time [s]\n";
    if (m_write_timestamp) {
        os << "# column " << col++ << ": UTC timestamp (" << m_datetime_format
           << ", with 'T' where the format has a space, so the column is one field)\n";
    }

    auto describe = [&](const StationLoc& loc, const StationVar& var, Real height, bool with_height)
    {
        os << "# column " << col++ << ": " << var.name;
        if (!var.units.empty()) { os << " [" << var.units << "]"; }
        if (loc.use_latlon) {
            os << " at requested lat=" << loc.req_lat << " long=" << loc.req_lon
               << " (sampled lat=" << loc.got_lat << " long=" << loc.got_lon << ")";
        }
        os << " at x=" << loc.x << " y=" << loc.y;
        if (with_height) { os << ", height=" << height << " m above local terrain"; }
        if (var.is_missing) { os << "  [NOT AVAILABLE in this run: constant " << var.missing_value << "]"; }
        os << "\n";
    };

    for (const auto& loc : station.locs) {
        for (const auto& var : station.vars) {
            if (var.is_2d) { describe(loc, var, Real(0.0), false); }
        }
        for (const auto height : station.heights) {
            for (const auto& var : station.vars) {
                if (!var.is_2d) { describe(loc, var, height, true); }
            }
        }
    }
    os << "#\n";
    os << "# Values are bilinearly interpolated in the horizontal and linearly\n";
    os << "# interpolated in the vertical, from the finest level that covers the\n";
    os << "# interpolation stencil.\n";
}

void
StationSampler::flush ()
{
    if (!ParallelDescriptor::IOProcessor()) { return; }

    const int nrow = static_cast<int>(m_times.size());
    if (nrow == 0) { return; }

    if (!m_opened) {
        if (!UtilCreateDirectory(m_dir, 0755)) { CreateDirectoryFailed(m_dir); }
    }

    for (const auto& station : m_stations)
    {
        const std::string filename = m_dir + "/" + station.name + ".dat";

        // A run that is restarting appends to the file it wrote before, so the
        // series is continuous; the restart is marked in a comment rather than
        // by starting a new file.  A run that is not restarting starts afresh.
        const bool append = m_opened || (m_is_restart && FileExists(filename));

        std::ofstream os(filename, append ? (std::ios::out | std::ios::app)
                                          : (std::ios::out | std::ios::trunc));
        if (!os.good()) { FileOpenFailed(filename); }

        if (!m_opened) {
            if (append) {
                os << "# restarted run resumes at t = "
                   << std::setprecision(timeprecision) << m_times[0] << "\n";
            } else {
                writeHeader(station, os);
            }
        }

        for (int irow = 0; irow < nrow; ++irow) {
            os << std::setw(datwidth) << std::setprecision(timeprecision) << m_times[irow];
            if (m_write_timestamp) {
                // Every column has to be one whitespace-delimited field for the
                // file to be readable by column, and the usual datetime formats
                // put a space between the date and the time.
                std::string stamp = getTimestamp(static_cast<double>(m_epoch_times[irow]),
                                                 m_datetime_format, false);
                std::replace(stamp.begin(), stamp.end(), ' ', 'T');
                os << std::setw(datwidth) << stamp;
            }
            for (const int icol : station.columns) {
                os << std::setw(datwidth) << std::setprecision(datprecision)
                   << m_rows[static_cast<std::size_t>(irow)*m_ncolumns + icol];
            }
            os << "\n";
        }
        os.close();
    }

    m_opened = true;
    m_times.clear();
    m_epoch_times.clear();
    m_rows.clear();
}

//
// Resolve the requested variable names against the variables this run can
// produce, and the requested lat/lon against the model grid.  Aborts on a name
// or a location that cannot be honored: a station column that silently vanished
// or silently moved would be worse than a failed run.
//
void
ERF::init_stations ()
{
    if (!station_sampler) { return; }

    auto& stations = station_sampler->stations();

    // The set of 2D diagnostics this configuration can actually fill, selected
    // exactly as Write2DPlotFile selects them.
    const bool has_surface_layer =
        phys_bc_type[Orientation(Direction::z, Orientation::low)] == ERF_BC::surface_layer;
    const auto active_lsm_names = lsm.Get_DataNames();
    const auto available_2d = plotfile2d::available_diagnostic_names(solverChoice,
                                                                    has_surface_layer,
                                                                    active_lsm_names);

    Vector<std::string> requested_3d;
    Vector<std::string> requested_2d;

    for (auto& station : stations)
    {
        for (auto& var : station.vars)
        {
            // Is it a 3D plot variable that this configuration can produce?
            Vector<std::string> probe{var.name};
            canonicalizePlot3DVariables(probe);
            appendPlot3DVariables(Vector<std::string>{var.name}, probe);
            if (containerHasElement(probe, var.name)) {
                var.is_2d = false;
                if (!containerHasElement(requested_3d, var.name)) {
                    requested_3d.push_back(var.name);
                }
                continue;
            }

            // Is it a 2D diagnostic?
            const auto* descriptor = plotfile2d::find_diagnostic(var.name);
            if (descriptor == nullptr) {
                descriptor = plotfile2d::find_dynamic_soil_diagnostic(var.name);
            }
            if (descriptor != nullptr) {
                var.is_2d  = true;
                var.units  = descriptor->units;
                if (containerHasElement(available_2d, var.name)) {
                    if (!containerHasElement(requested_2d, var.name)) {
                        requested_2d.push_back(var.name);
                    }
                } else {
                    // Same convention as the 2D plotfile: a diagnostic that this
                    // run does not compute is written as its documented missing
                    // value rather than dropped.
                    var.is_missing = true;
                    var.missing_value =
                        (descriptor->missing_policy == plotfile2d::MissingPolicy::FillZeroWhenUnavailable)
                        ? Real(0.0) : Real(-999.0);
                }
                continue;
            }

            // Neither: say which of the two failures it is.
            if (containerHasElement(derived_names, var.name) ||
                containerHasElement(cons_names, var.name)    ||
                var.name == "x_velocity" || var.name == "y_velocity" || var.name == "z_velocity")
            {
                Abort("Station '" + station.name + "': '" + var.name +
                      "' is a 3D plotfile variable but is not available in this configuration");
            }
            Abort("Station '" + station.name + "': '" + var.name +
                  "' is not a 3D or 2D plotfile variable.  See erf.plot_vars_1 and "
                  "erf.plot2d_vars_1 in the documentation for the names that can be requested");
        }
    }

    // The fill lists, in the order the fill routines produce their components.
    station_vars_3d = requested_3d;
    canonicalizePlot3DVariables(station_vars_3d);
    appendPlot3DVariables(requested_3d, station_vars_3d);

    station_vars_2d = plotfile2d::select_requested_plot_variables(requested_2d, available_2d).accepted;

    for (auto& station : stations) {
        for (auto& var : station.vars) {
            if (var.is_missing) { continue; }
            const auto& list = var.is_2d ? station_vars_2d : station_vars_3d;
            for (int i = 0; i < static_cast<int>(list.size()); ++i) {
                if (list[i] == var.name) { var.comp = i; break; }
            }
            AMREX_ALWAYS_ASSERT(var.comp >= 0);
        }

        if (station.has_3d_vars() && station.heights.empty()) {
            Abort(station_key_error(station.name, "height",
                                    "must be given, since this station requests a 3D variable"));
        }
    }

    resolve_station_positions();

    station_sampler->buildColumns();
    station_sampler->setWriteTimestamp(use_datetime, datetime_format);

    if (verbose > 0) {
        Print() << "Station output: " << stations.size() << " station(s), "
                << station_sampler->numColumns() << " column(s) per output step" << std::endl;
    }
}

//
// Turn every requested lat/lon into a position in domain coordinates, and check
// that every station lies inside the domain.  Done once: the level-0 geometry
// and the lat/lon arrays do not change.
//
void
ERF::resolve_station_positions ()
{
    auto& stations = station_sampler->stations();

    bool any_latlon = false;
    for (const auto& station : stations) {
        for (const auto& loc : station.locs) { any_latlon = any_latlon || loc.use_latlon; }
    }

    const auto  problo = geom[0].ProbLoArray();
    const auto  probhi = geom[0].ProbHiArray();
    const auto  dx0    = geom[0].CellSizeArray();
    const Box&  dom0   = geom[0].Domain();

    if (any_latlon)
    {
        if (lat_m[0] == nullptr || lon_m[0] == nullptr) {
            Abort("Station output: .lat/.long need a run with latitude/longitude arrays "
                  "(a WRF or metgrid initialization, or a restart from one).  Place the "
                  "stations with .x and .y in domain coordinates instead");
        }

        // Gather the level-0 mass-point latitude and longitude onto the IO rank.
        // This is a setup-time, level-0, 2D array, so the gather is affordable and
        // keeps the search and the inverse map as plain host code.
        //
        // NOTE: init_from_wrfinput reads WRF's staggered XLAT_V / XLONG_U into
        //       lat_m / lon_m (see the staggering contract in ERF.H), so the mass
        //       point is the average of the two bracketing edges.  init_from_metgrid
        //       reads mass-point values already.
        const bool destagger = (solverChoice.init_type == InitType::WRFInput);

        MultiFab latlon(ba2d[0], dmap[0], 2, 0);
        for (MFIter mfi(latlon, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            const Array4<Real>&       ll  = latlon.array(mfi);
            const Array4<const Real>& lat = lat_m[0]->const_array(mfi);
            const Array4<const Real>& lon = lon_m[0]->const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                ll(i,j,k,0) = destagger ? Real(0.5)*(lat(i,j,0) + lat(i,j+1,0)) : lat(i,j,0);
                ll(i,j,k,1) = destagger ? Real(0.5)*(lon(i,j,0) + lon(i+1,j,0)) : lon(i,j,0);
            });
        }

        Box dom2d(dom0); dom2d.setRange(2,0);
        BoxArray ba_one(dom2d);
        Vector<int> pmap(1, ParallelDescriptor::IOProcessorNumber());
        DistributionMapping dm_one(pmap);
        MultiFab latlon_all(ba_one, dm_one, 2, 0, MFInfo().SetArena(The_Pinned_Arena()));
        latlon_all.ParallelCopy(latlon, 0, 0, 2);

        const int ilo = dom0.smallEnd(0), ihi = dom0.bigEnd(0);
        const int jlo = dom0.smallEnd(1), jhi = dom0.bigEnd(1);

        for (auto& station : stations)
        {
            for (auto& loc : station.locs)
            {
                if (!loc.use_latlon) { continue; }

                Real out[4] = {Real(0.0), Real(0.0), Real(0.0), Real(0.0)};

                if (ParallelDescriptor::IOProcessor())
                {
                    const Array4<const Real>& ll = latlon_all.const_array(0);

                    auto cosfac = std::cos(loc.req_lat * PI / Real(180.0));
                    auto dist2  = [&](int i, int j) {
                        const Real dlat = ll(i,j,0,0) - loc.req_lat;
                        const Real dlon = (ll(i,j,0,1) - loc.req_lon) * cosfac;
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

                    // One linear solve in index space inverts the map: over a cell
                    // the projection is linear to well below a cell width.
                    const int im = std::max(bi-1, ilo), ip = std::min(bi+1, ihi);
                    const int jm = std::max(bj-1, jlo), jp = std::min(bj+1, jhi);
                    const Real inv_di = (ip > im) ? Real(1.0)/Real(ip-im) : Real(0.0);
                    const Real inv_dj = (jp > jm) ? Real(1.0)/Real(jp-jm) : Real(0.0);

                    const Real dlat_di = (ll(ip,bj,0,0) - ll(im,bj,0,0)) * inv_di;
                    const Real dlon_di = (ll(ip,bj,0,1) - ll(im,bj,0,1)) * inv_di;
                    const Real dlat_dj = (ll(bi,jp,0,0) - ll(bi,jm,0,0)) * inv_dj;
                    const Real dlon_dj = (ll(bi,jp,0,1) - ll(bi,jm,0,1)) * inv_dj;

                    const Real rlat = loc.req_lat - ll(bi,bj,0,0);
                    const Real rlon = loc.req_lon - ll(bi,bj,0,1);
                    Real di = Real(0.0), dj = Real(0.0);
                    if (!invert_latlon_offset(rlat, rlon, dlat_di, dlon_di, dlat_dj, dlon_dj, di, dj)) {
                        Abort("Station '" + station.name + "': the latitude/longitude arrays are "
                              "degenerate near the requested point, so it cannot be inverted");
                    }

                    if (std::abs(di) > Real(1.0) || std::abs(dj) > Real(1.0)) {
                        Abort("Station '" + station.name + "': requested lat=" +
                              std::to_string(loc.req_lat) + " long=" + std::to_string(loc.req_lon) +
                              " is outside the domain (the nearest grid point is more than one "
                              "cell away from it)");
                    }

                    out[0] = problo[0] + (Real(bi) + Real(0.5) + di) * dx0[0];
                    out[1] = problo[1] + (Real(bj) + Real(0.5) + dj) * dx0[1];
                    out[2] = ll(bi,bj,0,0) + dlat_di*di + dlat_dj*dj;
                    out[3] = ll(bi,bj,0,1) + dlon_di*di + dlon_dj*dj;
                }

                ParallelDescriptor::Bcast(out, 4, ParallelDescriptor::IOProcessorNumber());
                loc.x = out[0]; loc.y = out[1]; loc.got_lat = out[2]; loc.got_lon = out[3];
            }
        }
    }

    for (const auto& station : stations) {
        for (const auto& loc : station.locs) {
            if (loc.x < problo[0] || loc.x > probhi[0] ||
                loc.y < problo[1] || loc.y > probhi[1]) {
                Abort("Station '" + station.name + "': x=" + std::to_string(loc.x) +
                      " y=" + std::to_string(loc.y) + " is outside the problem domain");
            }
        }
    }
}

//
// Choose, for every station location, the finest level whose valid region
// covers the whole interpolation stencil, and record the stencil there.  Must
// be redone after every regrid, so it is simply redone every output step.
//
void
ERF::resolve_station_stencils ()
{
    auto& stations = station_sampler->stations();

    for (auto& station : stations)
    {
        for (auto& loc : station.locs)
        {
            bool found = false;

            for (int lev = finest_level; lev >= 0 && !found; --lev)
            {
                const auto  problo = geom[lev].ProbLoArray();
                const auto  dx     = geom[lev].CellSizeArray();
                const Box&  dom    = geom[lev].Domain();

                const int ilo = dom.smallEnd(0), ihi = dom.bigEnd(0);
                const int jlo = dom.smallEnd(1), jhi = dom.bigEnd(1);

                int ic = ilo + static_cast<int>(std::floor((loc.x - problo[0]) / dx[0]));
                int jc = jlo + static_cast<int>(std::floor((loc.y - problo[1]) / dx[1]));
                ic = std::min(std::max(ic, ilo), ihi);
                jc = std::min(std::max(jc, jlo), jhi);

                const Real tx = (loc.x - (problo[0] + (Real(ic-ilo) + Real(0.5))*dx[0])) / dx[0];
                const Real ty = (loc.y - (problo[1] + (Real(jc-jlo) + Real(0.5))*dx[1])) / dx[1];

                int  i0 = (tx >= Real(0.0)) ? ic : ic-1;
                int  j0 = (ty >= Real(0.0)) ? jc : jc-1;
                Real wx = (tx >= Real(0.0)) ? tx : tx + Real(1.0);
                Real wy = (ty >= Real(0.0)) ? ty : ty + Real(1.0);
                int  i1 = i0 + 1;
                int  j1 = j0 + 1;

                // Within the outer half cell of a non-periodic boundary there is
                // no second cell to interpolate from, so the stencil collapses to
                // the edge cell.
                if (!geom[lev].isPeriodic(0)) {
                    if (i0 <  ilo) { i0 = i1 = ilo; wx = Real(0.0); }
                    if (i1 >  ihi) { i0 = i1 = ihi; wx = Real(0.0); }
                }
                if (!geom[lev].isPeriodic(1)) {
                    if (j0 <  jlo) { j0 = j1 = jlo; wy = Real(0.0); }
                    if (j1 >  jhi) { j0 = j1 = jhi; wy = Real(0.0); }
                }

                // Every cell of the stencil must be covered by this level's valid
                // region -- after wrapping, where the direction is periodic -- or
                // the interpolation would read data this level does not have.
                auto covered = [&](int i, int j)
                {
                    int ii = i, jj = j;
                    if (geom[lev].isPeriodic(0)) {
                        const int nx = dom.length(0);
                        ii = ilo + ((i - ilo) % nx + nx) % nx;
                    }
                    if (geom[lev].isPeriodic(1)) {
                        const int ny = dom.length(1);
                        jj = jlo + ((j - jlo) % ny + ny) % ny;
                    }
                    if (ii < ilo || ii > ihi || jj < jlo || jj > jhi) { return false; }
                    Box column(IntVect(ii,jj,dom.smallEnd(2)), IntVect(ii,jj,dom.bigEnd(2)));
                    return grids[lev].contains(column);
                };

                if (covered(i0,j0) && covered(i1,j0) && covered(i0,j1) && covered(i1,j1)) {
                    loc.lev = lev;
                    loc.ic = ic; loc.jc = jc;
                    loc.i0 = i0; loc.j0 = j0;
                    loc.i1 = i1; loc.j1 = j1;
                    loc.wx = wx; loc.wy = wy;
                    found = true;
                }
            }

            // Level 0 covers the whole domain, so a station inside the domain
            // always resolves.
            AMREX_ALWAYS_ASSERT(found);
        }
    }
}

//
// Sample every station column and buffer one row.
//
// The values are gathered onto the IO rank by copying the small stencil boxes
// there and interpolating on that rank, rather than by interpolating in place
// and reducing: which rank owns a station changes at every regrid, and a copy
// of a few 2x2 columns is cheaper than the bookkeeping that would avoid it.
//
void
ERF::sample_stations (Real time)
{
    if (!station_sampler || station_sampler->empty()) { return; }

    BL_PROFILE("ERF::sample_stations()");

    auto& ss       = *station_sampler;
    auto& stations = ss.stations();

    resolve_station_stencils();

    const int n3d = static_cast<int>(station_vars_3d.size());
    const int n2d = static_cast<int>(station_vars_2d.size());

    // ---- fill the plot variables, on the levels that host a station ----
    Vector<int> lev_has_station(finest_level+1, 0);
    for (const auto& station : stations) {
        for (const auto& loc : station.locs) { lev_has_station[loc.lev] = 1; }
    }

    Vector<MultiFab> mf3d(finest_level+1);
    Vector<MultiFab> mf2d(finest_level+1);

    if (n3d > 0) {
        // The scratch has to be built for every level: the vorticity path fills
        // level l from level l-1.
        Plot3DScratch scratch;
        BuildPlot3DScratch(station_vars_3d, scratch);
        for (int lev = 0; lev <= finest_level; ++lev) {
            if (!lev_has_station[lev]) { continue; }
            mf3d[lev].define(grids[lev], dmap[lev], n3d, 0);
            FillPlot3DVars(lev, station_vars_3d, scratch, mf3d[lev], n3d, time);
        }
    }

    if (n2d > 0) {
        const auto descriptors =
            plotfile2d::build_sampled_level_output_descriptors_from_definitions({},
                                                                               station_vars_2d,
                                                                               solverChoice);
        for (int lev = 0; lev <= finest_level; ++lev) {
            if (!lev_has_station[lev]) { continue; }
            mf2d[lev].define(ba2d[lev], dmap[lev], n2d, 0);
            FillPlot2DVars(lev, station_vars_2d, station_vars_2d, descriptors, mf2d[lev], n2d);
        }
    }

    // ---- gather the stencils onto the IO rank ----
    const int nstation = static_cast<int>(stations.size());

    // Flat list of (station, location) pairs, so a stencil box and the location
    // it belongs to share an index.
    Vector<std::pair<int,int>> all_locs;
    for (int is = 0; is < nstation; ++is) {
        for (int il = 0; il < static_cast<int>(stations[is].locs.size()); ++il) {
            all_locs.emplace_back(is, il);
        }
    }

    const int io_rank = ParallelDescriptor::IOProcessorNumber();

    // Per level: the stencil boxes, in all_locs order, of the locations on it
    Vector<Vector<int>> lev_locs(finest_level+1);
    for (int ip = 0; ip < static_cast<int>(all_locs.size()); ++ip) {
        const auto& loc = stations[all_locs[ip].first].locs[all_locs[ip].second];
        lev_locs[loc.lev].push_back(ip);
    }

    // Stencil data brought to the IO rank, indexed by the position within lev_locs
    Vector<MultiFab> st3d(finest_level+1), st2d(finest_level+1);
    Vector<MultiFab> stzc(finest_level+1), stznd(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        if (lev_locs[lev].empty()) { continue; }

        const Box& dom = geom[lev].Domain();
        const int  klo = dom.smallEnd(2);
        const int  khi = dom.bigEnd(2);

        // Periodicity restricted to the horizontal: the stencils are full columns
        // in z, so there is nothing to wrap there.
        const Periodicity period(IntVect(geom[lev].isPeriodic(0) ? dom.length(0) : 0,
                                         geom[lev].isPeriodic(1) ? dom.length(1) : 0,
                                         0));

        BoxList bl_cc, bl_2d, bl_nd;
        for (const int ip : lev_locs[lev]) {
            const auto& loc = stations[all_locs[ip].first].locs[all_locs[ip].second];
            Box cc(IntVect(loc.i0,loc.j0,klo), IntVect(loc.i1,loc.j1,khi));
            bl_cc.push_back(cc);
            Box b2(IntVect(loc.i0,loc.j0,0), IntVect(loc.i1,loc.j1,0));
            bl_2d.push_back(b2);
            bl_nd.push_back(Box(IntVect(loc.i0,loc.j0,klo), IntVect(loc.i1+1,loc.j1+1,klo),
                                IntVect::TheNodeVector()));
        }

        const int nbox = static_cast<int>(bl_cc.size());
        Vector<int> pmap(nbox, io_rank);
        DistributionMapping dm_io(pmap);
        MFInfo pinned = MFInfo().SetArena(The_Pinned_Arena());

        if (n3d > 0) {
            BoxArray ba(bl_cc);
            st3d[lev].define(ba, dm_io, n3d, 0, pinned);
            st3d[lev].setVal(Real(0.0));
            st3d[lev].ParallelCopy(mf3d[lev], 0, 0, n3d, IntVect(0), IntVect(0), period);
        }
        if (n2d > 0) {
            BoxArray ba(bl_2d);
            st2d[lev].define(ba, dm_io, n2d, 0, pinned);
            st2d[lev].setVal(Real(0.0));
            st2d[lev].ParallelCopy(mf2d[lev], 0, 0, n2d, IntVect(0), IntVect(0), period);
        }
        if (z_phys_cc[lev]) {
            BoxArray ba(bl_cc);
            stzc[lev].define(ba, dm_io, 1, 0, pinned);
            stzc[lev].setVal(Real(0.0));
            stzc[lev].ParallelCopy(*z_phys_cc[lev], 0, 0, 1, IntVect(0), IntVect(0), period);
        }
        if (z_phys_nd[lev]) {
            BoxArray ba(bl_nd);
            stznd[lev].define(ba, dm_io, 1, 0, pinned);
            stznd[lev].setVal(Real(0.0));
            stznd[lev].ParallelCopy(*z_phys_nd[lev], 0, 0, 1, IntVect(0), IntVect(0), period);
        }
    }

    // ---- interpolate, on the IO rank ----
    Vector<Real> row(ss.numColumns(), Real(0.0));

    if (ParallelDescriptor::IOProcessor())
    {
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            if (lev_locs[lev].empty()) { continue; }

            const auto  problo = geom[lev].ProbLoArray();
            const auto  dx     = geom[lev].CellSizeArray();
            const Box&  dom    = geom[lev].Domain();
            const int   klo    = dom.smallEnd(2);
            const int   khi    = dom.bigEnd(2);

            for (int ib = 0; ib < static_cast<int>(lev_locs[lev].size()); ++ib)
            {
                const int ip = lev_locs[lev][ib];
                const Station&    station = stations[all_locs[ip].first];
                const StationLoc& loc     = station.locs[all_locs[ip].second];

                const int i0 = loc.i0, j0 = loc.j0, i1 = loc.i1, j1 = loc.j1;
                const Real wx = loc.wx, wy = loc.wy;

                auto bilinear = [&](const Array4<const Real>& a, int comp, int k)
                {
                    return (Real(1.0)-wy) * ( (Real(1.0)-wx)*a(i0,j0,k,comp) + wx*a(i1,j0,k,comp) )
                         +            wy  * ( (Real(1.0)-wx)*a(i0,j1,k,comp) + wx*a(i1,j1,k,comp) );
                };

                // Height of the cell centres in the station's column
                Array4<const Real> zc_arr;
                if (z_phys_cc[lev]) { zc_arr = stzc[lev].const_array(ib); }
                auto z_of_k = [&](int k)
                {
                    if (zc_arr) { return bilinear(zc_arr, 0, k); }
                    return problo[2] + (Real(k-klo) + Real(0.5)) * dx[2];
                };

                // Terrain elevation under the station
                Real z_surf = problo[2];
                if (z_phys_nd[lev]) {
                    const Array4<const Real>& nd = stznd[lev].const_array(ib);
                    auto corner_avg = [&](int i, int j) {
                        return Real(0.25) * (nd(i,j,klo) + nd(i+1,j,klo) + nd(i,j+1,klo) + nd(i+1,j+1,klo));
                    };
                    z_surf = (Real(1.0)-wy) * ( (Real(1.0)-wx)*corner_avg(i0,j0) + wx*corner_avg(i1,j0) )
                           +            wy  * ( (Real(1.0)-wx)*corner_avg(i0,j1) + wx*corner_avg(i1,j1) );
                }

                Array4<const Real> a3, a2;
                if (n3d > 0) { a3 = st3d[lev].const_array(ib); }
                if (n2d > 0) { a2 = st2d[lev].const_array(ib); }

                // The columns of this location, in the order buildColumns assigned them
                int icol = loc.col_begin;

                for (const auto& var : station.vars) {
                    if (!var.is_2d) { continue; }
                    const int gcol = station.columns[icol++];
                    row[gcol] = var.is_missing ? var.missing_value : bilinear(a2, var.comp, 0);
                }

                for (const Real height : station.heights)
                {
                    const Real z_target = z_surf + height;

                    // Bracket the target height in the station's column
                    int kk = klo;
                    while (kk < khi && z_of_k(kk+1) < z_target) { ++kk; }
                    const Real zk  = z_of_k(kk);
                    const Real zk1 = z_of_k(std::min(kk+1, khi));
                    Real wz = (zk1 > zk) ? (z_target - zk) / (zk1 - zk) : Real(0.0);
                    wz = std::min(std::max(wz, Real(0.0)), Real(1.0));   // clamp below the first
                    const int k1 = std::min(kk+1, khi);                  // cell centre and above the top

                    for (const auto& var : station.vars) {
                        if (var.is_2d) { continue; }
                        const int gcol = station.columns[icol++];
                        row[gcol] = var.is_missing
                                  ? var.missing_value
                                  : (Real(1.0)-wz)*bilinear(a3, var.comp, kk) + wz*bilinear(a3, var.comp, k1);
                    }
                }
            }
        }
    }

    ss.appendRow(time, static_cast<Real>(start_time) + time, row);
}

void
ERF::flush_stations ()
{
    if (station_sampler) { station_sampler->flush(); }
}
