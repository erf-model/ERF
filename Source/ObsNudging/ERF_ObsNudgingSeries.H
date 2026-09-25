/**
 * \file ERF_ObsNudgingSeries.H
 */
#ifndef ERF_OBSNUDGINGSERIES_H_
#define ERF_OBSNUDGINGSERIES_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <array>
#include <istream>
#include <string>

//
// The measurements of one observation-nudging station, as read from its file,
// and the profiles they give at a moment of the run.
//
// A station file is plain text.  Everything from a '#' to the end of a line is
// a comment.  The first line that is not blank names the columns, in any order:
//
//     time z u v w theta su sv sw stheta
//
// time and z are required; the wind is given either as u v (east and north, or
// grid x and y, see wind_frame) or as speed direction (metres per second, and
// the meteorological direction the wind blows from, in degrees clockwise from
// north).  w, theta and the standard deviations su sv sw stheta are optional,
// and a component with no column is not nudged at this station.  Each row that
// follows is one height at one time.  Rows are grouped by time in increasing
// order, and the heights within a time are in increasing order and the same at
// every time.  An entry that is not a finite number ("nan"), or that equals the
// run's missing value, is missing: that component is not nudged at that height
// and time.
//
namespace obs_nudging {

enum Comp : int { U = 0, V = 1, W = 2, Theta = 3, NComp = 4 };

// The name of a component, for messages
const char* comp_name (int comp);

// The measurements at one time.  mean[c][h] and sigma[c][h] are component c at
// height h; a missing mean is NaN, and a missing sigma is 0.
struct StationRecord
{
    double time = 0.0;
    std::array<amrex::Vector<amrex::Real>, NComp> mean;
    std::array<amrex::Vector<amrex::Real>, NComp> sigma;
};

struct StationSeries
{
    amrex::Vector<amrex::Real>    heights;   // strictly increasing
    amrex::Vector<StationRecord>  records;   // strictly increasing in time
    std::array<bool, NComp>       has = {false, false, false, false};
};

// Parse a station file from a stream.  Returns false and sets error (naming the
// line) when the file is malformed.  source is the name used in messages.
bool parse_station_series (std::istream& is, const std::string& source,
                           amrex::Real missing_value,
                           StationSeries& series, std::string& error);

// The measured profile of one component at one moment: the heights at which
// it is available, with its mean and standard deviation there.
struct Profile
{
    amrex::Vector<amrex::Real> z;
    amrex::Vector<amrex::Real> mean;
    amrex::Vector<amrex::Real> sigma;

    void clear () { z.clear(); mean.clear(); sigma.clear(); }
};

//
// The profiles of every component at time t.
//
// A file with one time holds for the whole run.  A file with several is
// interpolated linearly in time between the two records that bracket t, and a
// height is used only where both records have the component; before its first
// time and after its last the station is inactive and every profile is empty.
//
// The wind is then rotated from the frame of the file into the grid frame:
// (cos_alpha, sin_alpha) is the direction of the grid x axis measured
// counterclockwise from the file's x axis (true east for an earth-relative
// wind).  Unless the rotation is the identity it needs u and v together, so a
// height at which either is missing is dropped from both, and the standard
// deviations are rotated as those of independent components.
//
// Returns true when the station is active at t.
//
bool profiles_at_time (const StationSeries& series, double t,
                       amrex::Real cos_alpha, amrex::Real sin_alpha,
                       std::array<Profile, NComp>& profiles);

} // namespace obs_nudging

#endif
