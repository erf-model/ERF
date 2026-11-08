/**
 * \file ERF_StationSampler.H
 */
#ifndef ERF_STATIONSAMPLER_H_
#define ERF_STATIONSAMPLER_H_

#include <algorithm>
#include <cmath>
#include <string>

#include <AMReX_BoxArray.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

// invert_latlon_offset, used to place stations given by latitude/longitude
#include "ERF_LatLonMap.H"

//
// Station time-series output.
//
// A station is a named set of geographic points at which a user-chosen set of
// variables is written to an ASCII time series, for comparison against
// meteorological tower and surface-station observations.  The inputs follow the
// same shape as erf.refinement_indicators:
//
//     erf.station_names  = Lake1 Forests
//     erf.Lake1.field    = rain_accum
//     erf.Lake1.lat      = 45.13
//     erf.Lake1.long     = -122.34
//     erf.Forests.field  = magvel local_helicity
//     erf.Forests.lat    = 45.20 45.41
//     erf.Forests.long   = -122.10 -122.02
//     erf.Forests.height_agl = 10.0 80.0
//
// lat/long are paired positionally, so the example above is two locations, not
// four.  A station may instead be placed with .x/.y in domain coordinates, which
// is the only option for a run that has no lat/lon arrays.  The heights apply to
// every location of the station and are given one way or the other, never both:
// .height_agl is metres above the local terrain, .height_abs is metres in the
// model's own vertical coordinate, the same z as geometry.prob_lo/prob_hi.
//
// The variable names accepted are the 3D plotfile variables and the built-in 2D
// diagnostics, and the values come from the same fill routines
// (ERF::FillPlot3DVars / ERF::FillPlot2DVars), so a station column and the
// corresponding plotfile component cannot disagree.  A 2D plotfile can also
// carry sampled-level fields, named for the field and the level such as
// 'theta_z100m'; those are not accepted here, since a station asks for the 3D
// variable and a height directly.
//
// Four more names are accepted, which neither catalog holds as a 3D variable:
// the earth-relative winds u_east, v_north, wind_speed and wind_dir.  They are
// formed from the station's interpolated x_velocity and y_velocity, rotated by
// the grid orientation at the station (see ERF_EarthRelativeWind.H), so that a
// column compares directly with a tower on a rotated (e.g. Lambert) grid.
//

// Which earth-relative wind a station variable is; None for every other variable
enum class StationEarthWind { None, UEast, VNorth, WindSpeed, WindDir };

inline StationEarthWind
station_earth_wind_from_name (const std::string& name)
{
    if (name == "u_east")     { return StationEarthWind::UEast; }
    if (name == "v_north")    { return StationEarthWind::VNorth; }
    if (name == "wind_speed") { return StationEarthWind::WindSpeed; }
    if (name == "wind_dir")   { return StationEarthWind::WindDir; }
    return StationEarthWind::None;
}

struct StationVar
{
    std::string name;
    std::string units;
    bool is_2d      = false;              // comes from the 2D diagnostic catalog
    bool is_missing = false;              // valid name this run does not produce (an inactive
                                          // 2D diagnostic, or with
                                          // erf.station_allow_unavailable a 3D variable)
    bool cells_may_be_missing = false;    // active 2D field whose cells may hold missing_value
    amrex::Real missing_value = amrex::Real(0.0);
    int comp = -1;                        // component in the filled MultiFab
    // An earth-relative wind is formed from two 3D components: comp holds
    // x_velocity and comp_v y_velocity
    StationEarthWind earth_wind = StationEarthWind::None;
    int comp_v = -1;
};

//
// Bilinear interpolation of the four values around a station, with weights
// (1-wx)(1-wy), wx(1-wy), (1-wx)wy and wx*wy.  A value equal to missing_value, or
// not finite, is left out and the weights of the others are renormalised, so a
// missing neighbour cannot be blended in as a number; missing_value is returned
// only when every neighbour with a nonzero weight is missing.
//
inline amrex::Real
station_bilinear_skip_missing (amrex::Real v00, amrex::Real v10,
                               amrex::Real v01, amrex::Real v11,
                               amrex::Real wx, amrex::Real wy,
                               amrex::Real missing_value)
{
    const amrex::Real v[4] = {v00, v10, v01, v11};
    const amrex::Real w[4] = {(amrex::Real(1.0)-wx)*(amrex::Real(1.0)-wy),
                              wx*(amrex::Real(1.0)-wy),
                              (amrex::Real(1.0)-wx)*wy,
                              wx*wy};
    amrex::Real sum = amrex::Real(0.0);
    amrex::Real wsum = amrex::Real(0.0);
    for (int n = 0; n < 4; ++n) {
        if (w[n] > amrex::Real(0.0) && std::isfinite(v[n]) && v[n] != missing_value) {
            sum  += w[n] * v[n];
            wsum += w[n];
        }
    }
    return (wsum > amrex::Real(0.0)) ? sum / wsum : missing_value;
}

struct StationLoc
{
    bool use_latlon = false;

    // As requested in the inputs file
    amrex::Real req_lat = amrex::Real(0.0);
    amrex::Real req_lon = amrex::Real(0.0);

    // Domain coordinates of the station, resolved once at setup
    amrex::Real x = amrex::Real(0.0);
    amrex::Real y = amrex::Real(0.0);

    // lat/lon actually reached by (x,y); equals the request to within the
    // accuracy of the inverse map, and is reported in the file header
    amrex::Real got_lat = amrex::Real(0.0);
    amrex::Real got_lon = amrex::Real(0.0);

    // Orientation of the grid x axis at the station, counterclockwise from true
    // east, for the earth-relative winds.  Without lat/lon arrays the grid is
    // taken to be aligned with east.
    amrex::Real cos_alpha = amrex::Real(1.0);
    amrex::Real sin_alpha = amrex::Real(0.0);

    // Resolved against the current grids at every output step
    int lev = 0;
    int ic  = 0;                          // cell containing (x,y)
    int jc  = 0;
    int i0  = 0;                          // lower-left corner of the interpolation stencil
    int j0  = 0;
    int i1  = 0;                          // upper-right corner; equals i0 where the
    int j1  = 0;                          // stencil collapses at a physical boundary
    amrex::Real wx = amrex::Real(0.0);    // weight of i1 relative to i0
    amrex::Real wy = amrex::Real(0.0);
    int ktop = 0;                         // highest cell the vertical interpolation
                                          // reads, which is as high as the chosen
                                          // level is known to cover the stencil

    // Height of the chosen level's first cell centre above the zero the requested
    // heights are measured from; a requested height below it is sampled there
    amrex::Real first_cc_height = amrex::Real(0.0);

    // Index into the owning station's column list of this location's first column
    int col_begin = 0;
};

//
// One scalar written every output step: one variable, at one location, at one
// height.  Heights do not apply to 2D variables, whose height index is -1.
//
struct StationColumn
{
    int station = 0;
    int loc     = 0;
    int height  = -1;
    int var     = 0;
};

struct Station
{
    std::string name;
    amrex::Vector<StationVar>  vars;
    amrex::Vector<StationLoc>  locs;
    // The requested heights, in metres, and what they are measured from.  The
    // two input keys collapse to one list here: which one was given only decides
    // where the column's zero is, and that is a single flag rather than two
    // lists everything downstream would have to choose between.
    amrex::Vector<amrex::Real> heights;
    bool heights_are_agl = true;          // .height_agl, else .height_abs

    // Global column indices owned by this station, in file column order
    amrex::Vector<int> columns;

    bool has_3d_vars () const {
        for (const auto& v : vars) { if (!v.is_2d) { return true; } }
        return false;
    }
};

class StationSampler
{
public:
    // Parses erf.station_names and the per-station keys.  Name validation needs
    // the run's variable catalogs and is done by ERF::init_stations.
    explicit StationSampler (const std::string& pp_prefix);

    [[nodiscard]] bool empty () const { return m_stations.empty(); }

    amrex::Vector<Station>&       stations ()       { return m_stations; }
    [[nodiscard]] const amrex::Vector<Station>& stations () const { return m_stations; }

    [[nodiscard]] int numColumns () const { return m_ncolumns; }

    // Assigns the global column index of every (location, height, variable)
    // triple.  Called by ERF::init_stations once the variables are resolved.
    void buildColumns ();

    // Buffer one row of values.  Writes through to disk once the buffer is full.
    void appendRow (amrex::Real time, amrex::Real epoch_time, const amrex::Vector<amrex::Real>& values);

    // Write any buffered rows.  Safe to call from any rank; only the IO rank writes.
    void flush ();

    // A run that is not restarting starts its station files afresh; a restart
    // appends to what the earlier run wrote.
    void setRestart (bool is_restart) { m_is_restart = is_restart; }

    // On a restart, stop now if a file this run would append to describes a
    // different set of columns.  Called by ERF::init_stations once the columns
    // are resolved, so that a setup error costs a setup rather than a run.
    void checkRestartFiles () const;

    void setWriteTimestamp (bool write_timestamp, const std::string& datetime_format) {
        m_write_timestamp = write_timestamp;
        m_datetime_format = datetime_format;
    }

    [[nodiscard]] int bufferSteps () const { return m_buffer_steps; }

    // erf.station_allow_unavailable: write a 3D variable this configuration
    // cannot produce as a constant missing column instead of stopping the run
    [[nodiscard]] bool allowUnavailable () const { return m_allow_unavailable; }

    // A height below the first cell centre is honored by clamping, which is
    // worth saying once and not once per output step.
    [[nodiscard]] bool lowHeightWarned () const { return m_low_height_warned; }
    void setLowHeightWarned () { m_low_height_warned = true; }

    // Which level each location was resolved to is reported once, since it is
    // the thing about a station that is least obvious from the inputs file.
    [[nodiscard]] bool levelsReported () const { return m_levels_reported; }
    void setLevelsReported () { m_levels_reported = true; }

    // Resolving which level and which cells a station is read from costs a
    // ParallelCopy, a stream synchronization and a broadcast, and the answer is a
    // function of the grids and of the terrain under the station.  With terrain
    // that does not move, grids that have not changed give the answer already
    // stored in the locations, so the work can be skipped: these remember which
    // grids that was, and report whether it is still those.
    [[nodiscard]] bool stencilsValidFor (const amrex::Vector<amrex::BoxArray>& ba,
                                         int finest_level) const;
    void rememberStencilGrids (const amrex::Vector<amrex::BoxArray>& ba, int finest_level);

private:
    void writeHeader (const Station& station, std::ostream& os) const;

    // The format tag and a hash of what the columns are: the variables, their
    // units, the requested locations, the heights and their order.  This is what
    // a restart compares, rather than the header as a whole, so that rewording
    // the header or resolving a coordinate to a slightly different value does
    // not make an existing series un-restartable.
    [[nodiscard]] std::string columnSignature (const Station& station) const;

    // Stops a restart that would append columns the file's header does not
    // describe, which would otherwise pass unnoticed until the file was read.
    void checkHeaderMatches (const Station& station, const std::string& filename) const;

    amrex::Vector<Station> m_stations;
    int m_ncolumns    = 0;
    int m_buffer_steps = 100;
    bool m_allow_unavailable = false;

    // m_rows[irow] holds the values of every column at m_times[irow]
    amrex::Vector<amrex::Real> m_times;
    amrex::Vector<amrex::Real> m_epoch_times;
    amrex::Vector<amrex::Real> m_rows;

    // The grids the stencils currently stored in the locations were resolved
    // against; empty until the first resolve.
    amrex::Vector<amrex::BoxArray> m_stencil_grids;

    bool m_is_restart = false;
    bool m_low_height_warned = false;
    bool m_levels_reported = false;
    bool m_write_timestamp = false;
    std::string m_datetime_format;

    std::string m_dir = "Output_Stations";
    bool m_opened = false;
};

#endif
