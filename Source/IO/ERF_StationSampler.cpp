/**
 * \file ERF_StationSampler.cpp
 */

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <utility>
#include <vector>

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
#include <ERF_LatLonMap.H>
#include <ERF_TerrainSurfaceSlab.H>
#include <ERF_Utils.H>

using namespace amrex;

namespace {

constexpr int datwidth      = 20;
constexpr int datprecision  = 10;
constexpr int timeprecision = 13;

// The tag of the machine-readable line a restart compares.  Bump it when the
// meaning of the signature changes, so that an old file is rejected as an old
// file rather than as a changed configuration.
const char* const station_format_tag = "erf-station-v1";

// Requested coordinates and heights go into the signature at seven significant
// digits: enough to tell two stations apart, few enough that the comparison
// does not depend on the last bit of a Real, which is seven digits wide in a
// single-precision build.
constexpr int sigprecision = 7;

// FNV-1a, 64 bit.  Any stable hash would do; this one is a few lines and needs
// no library, and the signature is a guard against a changed configuration, not
// against a forged file.
std::string
signature_hash (const std::string& s)
{
    std::uint64_t h = 14695981039346656037ULL;
    for (const unsigned char c : s) {
        h ^= static_cast<std::uint64_t>(c);
        h *= 1099511628211ULL;
    }
    std::ostringstream os;
    os << std::hex << std::setw(16) << std::setfill('0') << h;
    return os.str();
}

// The format line of a station file, if it has one
std::string
find_format_line (const std::vector<std::string>& header)
{
    const std::string prefix = "# format: ";
    for (const auto& line : header) {
        if (line.compare(0, prefix.size(), prefix) == 0) {
            return line.substr(prefix.size());
        }
    }
    return {};
}

std::string
station_key_error (const std::string& station, const std::string& key, const std::string& why)
{
    return "Station '" + station + "': erf." + station + "." + key + " " + why;
}

//
// The leading comment block of a station file, which is the header the run that
// created it wrote.  Comments that appear later (the restart markers) are not
// part of it.
//
std::vector<std::string>
read_station_header (const std::string& filename)
{
    std::vector<std::string> header;
    std::ifstream is(filename);
    std::string line;
    while (std::getline(is, line)) {
        if (line.empty() || line[0] != '#') { break; }
        header.push_back(line);
    }
    return header;
}

//
// Drop the rows at or after t_first.  A restart is from a checkpoint, and the
// run that wrote the checkpoint may have gone on past it and flushed rows the
// restart is about to write again; without this the series would run backwards
// in time at the seam.  Returns the number of rows dropped.
//
int
truncate_station_file (const std::string& filename, amrex::Real t_first)
{
    std::vector<std::string> keep;
    int dropped = 0;

    // Keep the rows strictly before t_first.  The times are printed with
    // timeprecision significant digits, so the comparison needs a slack of a few
    // units in the last place to separate "the row we are about to write" from
    // "a row from earlier in the series"; scaling it by the magnitude of t_first
    // says that in one expression that stays right at t_first = 0 and for a
    // negative time, rather than relying on the run never reaching either.
    const double t_scale = std::max(1.0, std::abs(static_cast<double>(t_first)));
    const double t_cut   = static_cast<double>(t_first) - 1.0e-12 * t_scale;

    {
        std::ifstream is(filename);
        std::string line;
        while (std::getline(is, line)) {
            // Comments before the first dropped row are the header and the
            // restart markers of the part of the series that survives; comments
            // after it belong to the part that does not.
            if (line.empty() || line[0] == '#') {
                if (dropped == 0) { keep.push_back(line); }
                continue;
            }
            std::istringstream iss(line);
            double t = 0.0;
            if (!(iss >> t) || t < t_cut) { keep.push_back(line); } else { ++dropped; }
        }
    }

    // Write the survivors beside the file and move them into place, rather than
    // truncating the file and refilling it.  The rewrite is the one moment the
    // whole series exists only in this process's memory, and a crash there would
    // take the run's history with it; where the move replaces the file in one
    // step the worst a crash can leave behind is the original file and a stray
    // .tmp.
    if (dropped > 0) {
        const std::string tmpname = filename + ".tmp";
        {
            std::ofstream os(tmpname, std::ios::out | std::ios::trunc);
            if (!os.good()) { amrex::FileOpenFailed(tmpname); }
            for (const auto& line : keep) { os << line << "\n"; }
            os.flush();
            if (!os.good()) { amrex::FileOpenFailed(tmpname); }
        }
        // std::rename replaces an existing destination on POSIX and fails on it on
        // Windows, where a restart that had rows to drop would otherwise abort here.
        // The one-step replace is tried first, so nothing is lost to a crash on the
        // platforms that offer it, and only the platforms that do not pay for the
        // moment between the remove and the rename in which the series is the .tmp.
        if (std::rename(tmpname.c_str(), filename.c_str()) != 0) {
            std::remove(filename.c_str());
            if (std::rename(tmpname.c_str(), filename.c_str()) != 0) {
                amrex::Abort("Station output: could not move " + tmpname + " onto " + filename +
                             " while dropping the rows written past the restart point");
            }
        }
    }
    return dropped;
}

//
// A station name is both a ParmParse prefix and the name of the file the series
// is written to, so it has to be a plain file name: anything that could reach
// outside the output directory, or that a shell or a reader would have to quote,
// is refused here rather than turned into a surprising path.
//
void
check_station_name (const std::string& name)
{
    if (name.empty()) {
        amrex::Abort("Station names in erf.station_names cannot be empty");
    }
    if (name == "." || name == "..") {
        amrex::Abort("Station name '" + name + "' is not usable as a file name");
    }
    if (!(std::isalpha(static_cast<unsigned char>(name[0])) || name[0] == '_')) {
        amrex::Abort("Station name '" + name + "' must start with a letter or an underscore");
    }
    for (const char c : name) {
        if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-' || c == '.')) {
            amrex::Abort("Station name '" + name + "' contains '" + std::string(1,c) +
                         "'; station names are used as file names, so they are limited to "
                         "letters, digits, '_', '-' and '.'");
        }
    }
}

//
// The horizontal interpolation stencil of one location on one level, and how
// far up that level's valid region covers it.
//
struct StationLevelStencil
{
    int ic = 0, jc = 0;                   // cell containing (x,y)
    int i0 = 0, j0 = 0;                   // lower-left corner of the stencil
    int i1 = 0, j1 = 0;                   // upper-right corner
    amrex::Real wx = amrex::Real(0.0);
    amrex::Real wy = amrex::Real(0.0);
    int kcov = -1;                        // top of the coverage that reaches the bottom
                                          // of the domain; below the domain's smallEnd(2)
                                          // when the level cannot supply this location
};

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
        check_station_name(name);
        for (const auto& earlier : m_stations) {
            if (earlier.name == name) {
                amrex::Abort("Station '" + name + "' is named twice in erf.station_names; the "
                             "two would write the same file");
            }
        }

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

        // Heights come in one of two ways, and the difference is only what the
        // zero of the column is: .height_agl measures from the local terrain,
        // .height_abs from the bottom of the domain, in the same z the geometry
        // is given in.  Both are read into the one list the rest of the sampler
        // uses, with a flag saying which zero it is, so nothing downstream has
        // to ask which key the user typed.
        const int n_agl = pps.countval("height_agl");
        const int n_abs = pps.countval("height_abs");

        if (n_agl > 0 && n_abs > 0) {
            Abort("Station '" + name + "': give heights as either erf." + name +
                  ".height_agl (above the local terrain) or erf." + name +
                  ".height_abs (in the model's z coordinate), not both");
        }
        if (pps.contains("height")) {
            Abort("Station '" + name + "': erf." + name + ".height is not a key; say "
                  "erf." + name + ".height_agl for metres above the local terrain, or erf." +
                  name + ".height_abs for metres in the model's z coordinate");
        }

        station.heights_are_agl = (n_abs == 0);
        const char* height_key  = station.heights_are_agl ? "height_agl" : "height_abs";
        const int   nheight     = std::max(n_agl, n_abs);
        if (nheight > 0) {
            pps.getarr(height_key, station.heights, 0, nheight);
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

//
// A canonical description of what the columns of this station's file are: the
// station's name, whether it carries a timestamp column, and one token per
// column in file order, giving the variable, its units, the location, the
// height and whether the variable is one this run can actually fill.
//
// This is deliberately not the header the reader sees.  It leaves out the prose,
// so that rewording the header does not make every existing series
// un-restartable, and it uses the coordinates the inputs file asked for rather
// than the ones the setup resolved: a resolved x is derived from the lat/lon
// arrays and can move with the build or with init_type without the
// configuration having changed at all.
//
std::string
StationSampler::columnSignature (const Station& station) const
{
    std::ostringstream os;
    os << std::setprecision(sigprecision);
    os << station_format_tag << '|' << station.name
       << "|ts=" << (m_write_timestamp ? 1 : 0);

    auto describe = [&](const StationVar& var, Real height, bool with_height)
    {
        os << '|' << var.name << ':' << var.units;
        if (with_height) {
            os << ':' << height << (station.heights_are_agl ? "agl" : "abs");
        } else {
            os << ":2d";
        }
        if (var.is_missing) { os << ":na"; }
    };

    for (const auto& loc : station.locs) {
        os << "|@";
        if (loc.use_latlon) {
            os << "ll=" << loc.req_lat << ',' << loc.req_lon;
        } else {
            os << "xy=" << loc.x << ',' << loc.y;
        }
        for (const auto& var : station.vars) {
            if (var.is_2d) { describe(var, Real(0.0), false); }
        }
        for (const auto height : station.heights) {
            for (const auto& var : station.vars) {
                if (!var.is_2d) { describe(var, height, true); }
            }
        }
    }
    return std::string(station_format_tag) + " " + signature_hash(os.str());
}

void
StationSampler::writeHeader (const Station& station, std::ostream& os) const
{
    // Coordinates are written with enough digits that a reader can tell two
    // neighbouring stations apart
    os << std::setprecision(datprecision);
    os << "# ERF station time series\n";
    // The one line a restart compares.  Everything below it is for the reader.
    os << "# format: " << columnSignature(station) << "\n";
    os << "# station: " << station.name << "\n";

    int col = 1;
    os << "# column " << col++ << ": time [s]\n";
    if (m_write_timestamp) {
        os << "# column " << col++ << ": UTC timestamp (" << m_datetime_format
           << ", with 'T' where the format has a space, so the column is one field)\n";
    }

    const bool station_heights_are_agl = station.heights_are_agl;

    auto describe = [&](const StationLoc& loc, const StationVar& var, Real height, bool with_height)
    {
        os << "# column " << col++ << ": " << var.name;
        if (!var.units.empty()) { os << " [" << var.units << "]"; }
        if (loc.use_latlon) {
            os << " at requested lat=" << loc.req_lat << " long=" << loc.req_lon
               << " (sampled lat=" << loc.got_lat << " long=" << loc.got_lon << ")";
        }
        os << " at x=" << loc.x << " y=" << loc.y;
        if (with_height) {
            os << ", height=" << height
               << (station_heights_are_agl ? " m above local terrain" : " m (absolute, model z)");
        }
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

//
// A restart appends to the file the earlier run wrote, so the two runs have to
// agree about what the columns are.  What is compared is the signature line, not
// the header as a whole: the rest of the header is prose and resolved
// coordinates, and neither should decide whether a series can be continued.
//
// The signature decides; the human header is read only to say what changed, and
// only once the signature has already failed.
//
void
StationSampler::checkHeaderMatches (const Station& station, const std::string& filename) const
{
    const std::vector<std::string> found = read_station_header(filename);
    const std::string found_format = find_format_line(found);
    const std::string want_format  = columnSignature(station);

    if (found_format == want_format) { return; }

    if (found_format.empty()) {
        Abort("Station '" + station.name + "': " + filename + " has no '# format:' line, so it "
              "was not written by this version of ERF's station output and cannot be appended "
              "to.  Move it aside and start a new series.");
    }

    // A tag mismatch is a different thing from a hash mismatch and deserves a
    // different message: nothing about the run is wrong, the file is just old.
    const std::string found_tag = found_format.substr(0, found_format.find(' '));
    if (found_tag != std::string(station_format_tag)) {
        Abort("Station '" + station.name + "': " + filename + " is in station file format '" +
              found_tag + "' and this ERF writes '" + station_format_tag + "'.  Move it aside "
              "and start a new series.");
    }

    // Same format, different columns.  Point at the first column line that
    // differs, which is the useful thing to say; the signature is what decided.
    std::ostringstream expected;
    writeHeader(station, expected);
    std::vector<std::string> want;
    {
        std::istringstream is(expected.str());
        std::string line;
        while (std::getline(is, line)) { want.push_back(line); }
    }

    // The signature lines differ by construction, so they say nothing a reader
    // does not already know; what is worth reporting is the first column that
    // differs.
    auto without_signature = [](const std::vector<std::string>& lines) {
        std::vector<std::string> out;
        for (const auto& line : lines) {
            if (line.compare(0, 10, "# format: ") != 0) { out.push_back(line); }
        }
        return out;
    };
    const std::vector<std::string> want_cmp  = without_signature(want);
    const std::vector<std::string> found_cmp = without_signature(found);

    std::string detail = "\n  the difference is not one the header spells out";
    for (std::size_t i = 0; i < want_cmp.size() || i < found_cmp.size(); ++i) {
        const bool have_want  = (i < want_cmp.size());
        const bool have_found = (i < found_cmp.size());
        if (have_want && have_found && want_cmp[i] == found_cmp[i]) { continue; }
        if (have_want && have_found) {
            detail = "\n  the file says:        " + found_cmp[i] +
                     "\n  this run would write: " + want_cmp[i];
        } else if (have_found) {
            detail = "\n  the file has a column this run does not: " + found_cmp[i];
        } else {
            detail = "\n  this run has a column the file does not: " + want_cmp[i];
        }
        break;
    }

    Abort("Station '" + station.name + "': " + filename + " describes a different set of "
          "columns than this restart would write -- a changed field, location, height or "
          "units list -- so appending to it would produce columns its header does not "
          "describe.  Change the configuration back, or move the file aside and start a "
          "new series." + detail);
}

//
// A station file this restart would append to has to describe the same columns.
// That is a setup error when it is wrong, so it is checked here, before the run
// does any work, rather than at the first flush -- which is the first checkpoint
// or station_buffer_steps rows in, by which time the run has spent real time on
// a series it cannot write.
//
void
StationSampler::checkRestartFiles () const
{
    if (!m_is_restart) { return; }
    if (!ParallelDescriptor::IOProcessor()) { return; }

    for (const auto& station : m_stations) {
        const std::string filename = m_dir + "/" + station.name + ".dat";
        if (FileExists(filename)) { checkHeaderMatches(station, filename); }
    }
}

//
// Whether the stencils stored in the locations were resolved against these grids.
// BoxArray equality is a pointer comparison when the two share a representation,
// which is the common case between steps that did not regrid, and a content
// comparison otherwise; either way it is far cheaper than resolving again.
//
bool
StationSampler::stencilsValidFor (const Vector<BoxArray>& ba, int finest_level) const
{
    if (static_cast<int>(m_stencil_grids.size()) != finest_level+1) { return false; }
    for (int lev = 0; lev <= finest_level; ++lev) {
        if (m_stencil_grids[lev] != ba[lev]) { return false; }
    }
    return true;
}

void
StationSampler::rememberStencilGrids (const Vector<BoxArray>& ba, int finest_level)
{
    m_stencil_grids.resize(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) { m_stencil_grids[lev] = ba[lev]; }
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

        // The first thing a restart writes into an existing file: drop anything
        // the earlier run wrote past the checkpoint this run restarted from.
        // That the file is the same series was settled at setup, by
        // checkRestartFiles.
        if (append && !m_opened) {
            const int dropped = truncate_station_file(filename, m_times[0]);
            if (dropped > 0) {
                Print() << "Station output: dropped " << dropped << " row(s) at or after t = "
                        << m_times[0] << " from " << filename << ", written by the run before "
                        << "the restart past the checkpoint it restarted from" << std::endl;
            }
        }

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
                  "' is neither a 3D plotfile variable nor a built-in 2D diagnostic, which "
                  "are the two catalogs a station can draw on.  Note that a 2D plotfile can "
                  "also carry sampled-level fields, named for the field and the level such "
                  "as 'theta_z100m', and those are not among them: ask for the 3D variable "
                  "itself and give the station a height instead.  See erf.plot_vars_1 and "
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
            Abort("Station '" + station.name + "': erf." + station.name + ".height_agl or erf." +
                  station.name + ".height_abs must be given, since this station requests a 3D "
                  "variable");
        }
    }

    resolve_station_positions();

    station_sampler->buildColumns();
    station_sampler->setWriteTimestamp(use_datetime, datetime_format);

    // Everything the signature describes is now resolved, so a restart that
    // would append to a file describing different columns stops here, before the
    // run does any work.
    station_sampler->checkRestartFiles();

    if (verbose > 0) {
        Print() << "Station output: " << stations.size() << " station(s), "
                << station_sampler->numColumns() << " column(s) per output step" << std::endl;
    }
}

//
// The terrain elevation at a station: bilinear in the elevations of the nodes
// of the bottom of the mesh (or of the immersed terrain surface) around it.
// The nodes gathered for a station are those of its cell-centred stencil,
// i0..i1+1 by j0..j1+1, which always contain the cell of nodes around (x, y).
//
// This is the elevation of the ground at the station itself.  Interpolating
// the cell averages of the elevation between cell centres instead flattens a
// curved hill, and put a mast on the flank of a 100 m hill 3.4 m too low at
// 50 m resolution.
//
static Real
terrain_at_station (const Array4<const Real>& nd, int klo, Real x, Real y,
                    const GpuArray<Real,AMREX_SPACEDIM>& problo,
                    const GpuArray<Real,AMREX_SPACEDIM>& dx,
                    int i0, int i1, int j0, int j1)
{
    const Real fx = (x - problo[0]) / dx[0];
    const Real fy = (y - problo[1]) / dx[1];
    const int  in = std::clamp(static_cast<int>(std::floor(fx)), i0, i1);
    const int  jn = std::clamp(static_cast<int>(std::floor(fy)), j0, j1);
    const Real tx = std::clamp(fx - Real(in), Real(0.0), Real(1.0));
    const Real ty = std::clamp(fy - Real(jn), Real(0.0), Real(1.0));
    return (Real(1.0)-ty) * ( (Real(1.0)-tx)*nd(in,jn  ,klo) + tx*nd(in+1,jn  ,klo) )
         +            ty  * ( (Real(1.0)-tx)*nd(in,jn+1,klo) + tx*nd(in+1,jn+1,klo) );
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

    if (any_latlon)
    {
        if (lat_m[0] == nullptr || lon_m[0] == nullptr) {
            Abort("Station output: .lat/.long need a run with latitude/longitude arrays "
                  "(a WRF or metgrid initialization, or a restart from one).  Place the "
                  "stations with .x and .y in domain coordinates instead");
        }

        // NOTE: every init path fills lat_m / lon_m with mass-point values --
        //       init_from_wrfinput reads WRF's XLAT / XLONG and init_from_metgrid
        //       reads XLAT_M / XLONG_M -- so no averaging is needed here.
        const LatLonMap latlon_map(*lat_m[0], *lon_m[0], ba2d[0], dmap[0], geom[0]);

        for (auto& station : stations)
        {
            for (auto& loc : station.locs)
            {
                if (!loc.use_latlon) { continue; }

                LatLonLocation where;
                const LatLonStatus status = latlon_map.locate(loc.req_lat, loc.req_lon, where);

                if (status == LatLonStatus::Degenerate) {
                    Abort("Station '" + station.name + "': the latitude/longitude arrays are "
                          "degenerate near the requested point, so it cannot be inverted");
                }
                if (status == LatLonStatus::TooFar) {
                    Abort("Station '" + station.name + "': requested lat=" +
                          std::to_string(loc.req_lat) + " long=" + std::to_string(loc.req_lon) +
                          " resolved to a point more than one cell from the nearest grid "
                          "point (lat=" + std::to_string(where.near_lat) + " long=" +
                          std::to_string(where.near_lon) + "), which means it is outside the "
                          "domain, or that the latitude/longitude arrays are too distorted "
                          "near it to invert");
                }

                loc.x = where.x; loc.y = where.y; loc.got_lat = where.lat; loc.got_lon = where.lon;
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
// Choose, for every station location, the finest level that can supply the
// interpolation, and record the stencil and the top of the column it reads.
// Must be redone after every regrid, so it is simply redone every output step.
//
// A level can supply a location when its valid region covers the four columns
// of the horizontal stencil from the bottom of the domain up through the cells
// the vertical interpolation reads.  Coverage is required from the bottom
// because a station height is measured from the local terrain, which is the
// bottom of the column, and only as far up as the requested heights reach: a
// level that refines the boundary layer and stops part way up the domain is
// exactly the level a station wants, and requiring it to cover its whole
// vertical extent -- which is the refinement of the coarse domain, not of the
// refined region -- would reject it for every station.
//
void
ERF::resolve_station_stencils ()
{
    // What this computes -- which level supplies each location, which four cells
    // its horizontal stencil reads, and how far up the column -- is a function of
    // the grids and of the terrain under the station, and it costs a ParallelCopy,
    // a stream synchronization and a broadcast to compute.  Neither input changes
    // between samples unless the grids do, so the answer already in the locations
    // stands and none of that work is needed.
    //
    // The exception is a mesh whose terrain moves: there the column heights change
    // under a fixed grid, so the level a height falls in can change with the grids
    // unchanged and the stencils have to be resolved every time.
    const bool terrain_moves = (solverChoice.terrain_type == TerrainType::MovingFittedMesh);
    if (!terrain_moves && station_sampler->stencilsValidFor(grids, finest_level)) { return; }

    auto& stations = station_sampler->stations();

    // Flat list of (station, location) pairs, so a stencil and the location it
    // belongs to share an index.
    Vector<std::pair<int,int>> all_locs;
    for (int is = 0; is < static_cast<int>(stations.size()); ++is) {
        for (int il = 0; il < static_cast<int>(stations[is].locs.size()); ++il) {
            all_locs.emplace_back(is, il);
        }
    }
    const int nloc = static_cast<int>(all_locs.size());

    // ---- the horizontal stencil, and the coverage above it, on every level ----
    Vector<Vector<StationLevelStencil>> cand(nloc, Vector<StationLevelStencil>(finest_level+1));

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        const auto  problo = geom[lev].ProbLoArray();
        const auto  dx     = geom[lev].CellSizeArray();
        const Box&  dom    = geom[lev].Domain();

        const int ilo = dom.smallEnd(0), ihi = dom.bigEnd(0);
        const int jlo = dom.smallEnd(1), jhi = dom.bigEnd(1);
        const int klo = dom.smallEnd(2), khi = dom.bigEnd(2);

        // How far up this level's valid region covers the column at (i,j) --
        // after wrapping, where the direction is periodic -- counting only
        // coverage that is contiguous from the bottom of the domain.  Returns
        // klo-1 when even the bottom cell is not covered.
        auto column_top = [&](int i, int j)
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
            if (ii < ilo || ii > ihi || jj < jlo || jj > jhi) { return klo-1; }

            std::vector<std::pair<int,Box>> isects;
            grids[lev].intersections(Box(IntVect(ii,jj,klo), IntVect(ii,jj,khi)), isects);

            // The pieces of the column come back in no particular order, so grow
            // the covered span from the bottom until nothing else abuts it.
            int ktop = klo-1;
            bool grew = true;
            while (grew) {
                grew = false;
                for (const auto& is : isects) {
                    if (is.second.smallEnd(2) <= ktop+1 && is.second.bigEnd(2) > ktop) {
                        ktop = is.second.bigEnd(2);
                        grew = true;
                    }
                }
            }
            return ktop;
        };

        for (int ip = 0; ip < nloc; ++ip)
        {
            const auto& loc = stations[all_locs[ip].first].locs[all_locs[ip].second];
            auto& c = cand[ip][lev];

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

            c.ic = ic; c.jc = jc;
            c.i0 = i0; c.j0 = j0;
            c.i1 = i1; c.j1 = j1;
            c.wx = wx; c.wy = wy;
            c.kcov = std::min(std::min(column_top(i0,j0), column_top(i1,j0)),
                              std::min(column_top(i0,j1), column_top(i1,j1)));
        }
    }

    // ---- the cell centre heights of the candidate columns, on the IO rank ----
    //
    // Which cells the vertical interpolation reads depends on the terrain under
    // the station, so choosing the level needs the heights themselves.  They are
    // one component of a 2x2 column, so they are gathered where the choice is
    // made, on the same rank that later interpolates.
    //
    const int io_rank = ParallelDescriptor::IOProcessorNumber();

    Vector<Vector<int>> lev_cand(finest_level+1);      // ip of each gathered box, in box order
    Vector<MultiFab> zc(finest_level+1), znd(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        const Box& dom = geom[lev].Domain();
        const int  klo = dom.smallEnd(2);

        for (int ip = 0; ip < nloc; ++ip) {
            if (cand[ip][lev].kcov >= klo) { lev_cand[lev].push_back(ip); }
        }
        // Terrain gives both arrays or neither: z_of_k below reads z_phys_cc and
        // falls back to a uniform dz without it, while z_surf reads z_phys_nd, so
        // a level with one and not the other would measure a height from the
        // terrain and then look it up on a flat column.  No init path produces
        // that, and this says so rather than leaving it to be discovered.
        AMREX_ALWAYS_ASSERT((z_phys_cc[lev] != nullptr) == (z_phys_nd[lev] != nullptr));
        if (lev_cand[lev].empty() || !z_phys_cc[lev]) { continue; }

        // Periodicity restricted to the horizontal: the columns start at the
        // bottom of the domain, so there is nothing to wrap in z.
        const Periodicity period(IntVect(geom[lev].isPeriodic(0) ? dom.length(0) : 0,
                                         geom[lev].isPeriodic(1) ? dom.length(1) : 0,
                                         0));

        BoxList bl_cc, bl_nd;
        for (const int ip : lev_cand[lev]) {
            const auto& c = cand[ip][lev];
            bl_cc.push_back(Box(IntVect(c.i0,c.j0,klo), IntVect(c.i1,c.j1,c.kcov)));
            bl_nd.push_back(Box(IntVect(c.i0,c.j0,klo), IntVect(c.i1+1,c.j1+1,klo),
                                IntVect::TheNodeVector()));
        }

        Vector<int> pmap(static_cast<int>(bl_cc.size()), io_rank);
        DistributionMapping dm_io(pmap);
        MFInfo pinned = MFInfo().SetArena(The_Pinned_Arena());

        if (z_phys_cc[lev]) {
            zc[lev].define(BoxArray(bl_cc), dm_io, 1, 0, pinned);
            zc[lev].setVal(Real(0.0));
            zc[lev].ParallelCopy(*z_phys_cc[lev], 0, 0, 1, IntVect(0), IntVect(0), period);
        }
        if (z_phys_nd[lev]) {
            znd[lev].define(BoxArray(bl_nd), dm_io, 1, 0, pinned);
            znd[lev].setVal(Real(0.0));
            znd[lev].ParallelCopy(station_ground(lev), 0, 0, 1, IntVect(0), IntVect(0), period);
        }
    }

    // The gathers above are device-side in a GPU build, and the choice below
    // reads the pinned destinations on the host.
    Gpu::streamSynchronize();

    // ---- choose the level, on the IO rank, and tell everyone ----
    //
    // Only the level and the top of the column have to be broadcast: every rank
    // computed the same stencils above.
    //
    Vector<int> chosen(2*nloc, -1);

    if (ParallelDescriptor::IOProcessor())
    {
        for (int ip = 0; ip < nloc; ++ip)
        {
            const Station& station = stations[all_locs[ip].first];

            // A station of 2D variables only never leaves the bottom of the column
            const bool reads_column = station.has_3d_vars() && !station.heights.empty();

            for (int lev = finest_level; lev >= 0 && chosen[2*ip] < 0; --lev)
            {
                const auto& c = cand[ip][lev];

                const auto problo = geom[lev].ProbLoArray();
                const auto dx     = geom[lev].CellSizeArray();
                const Box& dom    = geom[lev].Domain();
                const int  klo    = dom.smallEnd(2);
                const int  khi    = dom.bigEnd(2);

                if (c.kcov < klo) { continue; }

                if (!reads_column) {
                    chosen[2*ip] = lev; chosen[2*ip+1] = klo;
                    break;
                }

                const int ib = static_cast<int>(
                    std::find(lev_cand[lev].begin(), lev_cand[lev].end(), ip) - lev_cand[lev].begin());

                auto bilinear = [&](const Array4<const Real>& a, int k)
                {
                    return (Real(1.0)-c.wy) * ( (Real(1.0)-c.wx)*a(c.i0,c.j0,k) + c.wx*a(c.i1,c.j0,k) )
                         +            c.wy  * ( (Real(1.0)-c.wx)*a(c.i0,c.j1,k) + c.wx*a(c.i1,c.j1,k) );
                };

                Array4<const Real> zc_arr;
                if (z_phys_cc[lev]) { zc_arr = zc[lev].const_array(ib); }
                auto z_of_k = [&](int k)
                {
                    if (zc_arr) { return bilinear(zc_arr, k); }
                    return problo[2] + (Real(k-klo) + Real(0.5)) * dx[2];
                };

                // Terrain elevation under the station, which is where an AGL
                // height is measured from
                Real z_surf = problo[2];
                if (z_phys_nd[lev]) {
                    const StationLoc& sloc = station.locs[all_locs[ip].second];
                    z_surf = terrain_at_station(znd[lev].const_array(ib), klo, sloc.x, sloc.y,
                                                problo, dx, c.i0, c.i1, c.j0, c.j1);
                }

                // Every requested height must be bracketed within the covered
                // part of the column.  Running out of coverage is fatal to the
                // level unless the coverage reaches the top of the domain, where
                // the interpolation legitimately clamps.
                // The zero of the requested heights: the terrain under the
                // station, or the bottom of the domain.
                const Real z_zero = station.heights_are_agl ? z_surf : Real(0.0);

                int  ktop = klo;
                bool usable = true;
                bool below_first_cell = false;
                for (const Real height : station.heights)
                {
                    const Real z_target = z_zero + height;
                    if (z_of_k(c.kcov) < z_target && c.kcov < khi) { usable = false; break; }

                    if (z_target < z_of_k(klo)) { below_first_cell = true; }

                    int kk = klo;
                    while (kk < c.kcov && z_of_k(kk+1) < z_target) { ++kk; }
                    ktop = std::max(ktop, std::min(kk+1, c.kcov));
                }

                if (usable) {
                    chosen[2*ip] = lev;
                    chosen[2*ip+1] = ktop;

                    // Nothing in the model lives between the terrain and the first
                    // cell centre, so a height there is the first cell centre's
                    // value.  That is a trap for exactly the heights an observation
                    // comparison asks for, so say so, once.
                    if (below_first_cell && !station_sampler->lowHeightWarned()) {
                        station_sampler->setLowHeightWarned();
                        Warning("Station '" + station.name + "': a requested height is below the "
                                "first cell centre of the level it is sampled from, so it is "
                                "reported as the value at that cell centre -- there is no "
                                "similarity extrapolation to the requested height.  For 2 m or "
                                "10 m quantities, request the 2D diagnostics (temperature_2m, "
                                "water_vapor_mixing_ratio_2m, and the surface-layer diagnostics) "
                                "instead of a 3D variable at that height");
                    }
                }
            }
        }
    }

    ParallelDescriptor::Bcast(chosen.data(), chosen.size(), io_rank);

    for (int ip = 0; ip < nloc; ++ip)
    {
        // Level 0 covers the whole domain from the ground up, so a station
        // inside the domain always resolves.
        AMREX_ALWAYS_ASSERT(chosen[2*ip] >= 0);

        auto& loc     = stations[all_locs[ip].first].locs[all_locs[ip].second];
        const auto& c = cand[ip][chosen[2*ip]];

        loc.lev = chosen[2*ip];
        loc.ic = c.ic; loc.jc = c.jc;
        loc.i0 = c.i0; loc.j0 = c.j0;
        loc.i1 = c.i1; loc.j1 = c.j1;
        loc.wx = c.wx; loc.wy = c.wy;
        loc.ktop = chosen[2*ip+1];
    }

    station_sampler->rememberStencilGrids(grids, finest_level);

    // Which level a station is sampled from follows from the grids rather than
    // from anything in the inputs file, so report it once.
    if (verbose > 0 && !station_sampler->levelsReported())
    {
        station_sampler->setLevelsReported();
        for (int ip = 0; ip < nloc; ++ip) {
            const Station&    station = stations[all_locs[ip].first];
            const StationLoc& loc     = station.locs[all_locs[ip].second];
            Print() << "Station output: '" << station.name << "' at x=" << loc.x
                    << " y=" << loc.y << " is sampled from level " << loc.lev
                    << " (cell " << loc.ic << "," << loc.jc
                    << ", column through k=" << loc.ktop << ")" << std::endl;
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
    int top_station_lev = 0;
    for (const auto& station : stations) {
        for (const auto& loc : station.locs) {
            lev_has_station[loc.lev] = 1;
            top_station_lev = std::max(top_station_lev, loc.lev);
        }
    }

    Vector<MultiFab> mf3d(finest_level+1);
    Vector<MultiFab> mf2d(finest_level+1);

    if (n3d > 0) {
        // The scratch has to be built for every level up to the finest one a
        // station is on, not just the levels that host a station: the vorticity
        // path fills level l from level l-1.  Levels above that are not read, so
        // they are not built.  The sampler runs at every output step, so it also
        // asks BuildPlot3DScratch not to average the microphysics state down:
        // turning station output on must not change the answer.
        Plot3DScratch scratch;
        BuildPlot3DScratch(station_vars_3d, scratch, top_station_lev, false);
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

        // Periodicity restricted to the horizontal: the stencils start at the
        // bottom of the domain, so there is nothing to wrap in z.
        const Periodicity period(IntVect(geom[lev].isPeriodic(0) ? dom.length(0) : 0,
                                         geom[lev].isPeriodic(1) ? dom.length(1) : 0,
                                         0));

        // Each stencil reaches only as far up as resolve_station_stencils found
        // the vertical interpolation needs, which is as far up as this level is
        // guaranteed to cover it.
        // One box per location, in lev_locs order, so a gathered box and the
        // location it belongs to share an index.  Two locations inside the same
        // cell give identical boxes, which is allowed: the BoxArray is only ever
        // a list of what to copy and what to index, never a cover of a region,
        // and ParallelCopy fills each of a pair of identical boxes from the same
        // source.  The 2D boxes sit at k = 0 rather than at klo because ba2d is
        // built by compressing the 3D BoxArray with setRange(2,0), so 0 is where
        // its data is whatever the domain's smallEnd is (ERF_MakeNewArrays.cpp).
        BoxList bl_cc, bl_2d, bl_nd;
        for (const int ip : lev_locs[lev]) {
            const auto& loc = stations[all_locs[ip].first].locs[all_locs[ip].second];
            Box cc(IntVect(loc.i0,loc.j0,klo), IntVect(loc.i1,loc.j1,loc.ktop));
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
            stznd[lev].ParallelCopy(station_ground(lev), 0, 0, 1, IntVect(0), IntVect(0), period);
        }
    }

    // The gathers above are device-side in a GPU build, and what follows reads
    // the pinned destinations on the host.
    Gpu::streamSynchronize();

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
                    z_surf = terrain_at_station(stznd[lev].const_array(ib), klo, loc.x, loc.y,
                                                problo, dx, i0, i1, j0, j1);
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

                const Real z_zero = station.heights_are_agl ? z_surf : Real(0.0);

                for (const Real height : station.heights)
                {
                    const Real z_target = z_zero + height;

                    // Bracket the target height in the part of the column this
                    // level was chosen to cover
                    int kk = klo;
                    while (kk < loc.ktop && z_of_k(kk+1) < z_target) { ++kk; }
                    const Real zk  = z_of_k(kk);
                    const Real zk1 = z_of_k(std::min(kk+1, loc.ktop));
                    Real wz = (zk1 > zk) ? (z_target - zk) / (zk1 - zk) : Real(0.0);
                    wz = std::min(std::max(wz, Real(0.0)), Real(1.0));   // clamp below the first
                    const int k1 = std::min(kk+1, loc.ktop);             // cell centre and above the top

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

//
// The nodes a station's terrain elevation is read from.  On a terrain-fitted
// mesh (and on a flat one) that is the bottom of the mesh, z_phys_nd.  With
// immersed-forcing or embedded-boundary terrain the mesh is flat and the
// terrain is immersed in it, so a height above the local terrain has to be
// measured from the surface the immersed boundary is built from; that is built
// here once per set of grids.  Both build the surface from init_terrain_surface
// (see ERF::initializeEB and the immersed-boundary setup), so the same slab serves
// both.
//
// An EB run that specifies no terrain surface has a flat boundary, but its
// z_phys is shifted so that an eb2.geometry = plane boundary lies at zero (the
// z_offset of init_default_zphys); the ground is that zero, not z_phys_nd,
// which holds the shifted bottom of the index space.
//
const MultiFab&
ERF::station_ground (int lev)
{
    AMREX_ALWAYS_ASSERT(z_phys_nd[lev] != nullptr);
    const bool terrain_is_immersed = (solverChoice.terrain_type == TerrainType::ImmersedForcing ||
                                      solverChoice.terrain_type == TerrainType::EB);
    if (!terrain_is_immersed) {
        return *z_phys_nd[lev];
    }

    if (static_cast<int>(station_ib_ground.size()) <= lev) {
        station_ib_ground.resize(lev+1);
    }
    const BoxArray slab = bottom_node_slab(grids[lev], geom[lev]);
    if (!station_ib_ground[lev] ||
        !(station_ib_ground[lev]->boxArray() == slab) ||
        !(station_ib_ground[lev]->DistributionMap() == dmap[lev]))
    {
        station_ib_ground[lev] = std::make_unique<MultiFab>(slab, dmap[lev], 1, 0);
        if (solverChoice.terrain_type == TerrainType::EB && !prob->terrain_is_specified()) {
            station_ib_ground[lev]->setVal(Real(0.0));
        } else {
            fill_terrain_surface_slab(*station_ib_ground[lev], geom[lev], *prob, t_new[lev]);
        }
    }
    return *station_ib_ground[lev];
}

void
ERF::flush_stations () const
{
    if (station_sampler) { station_sampler->flush(); }
}
