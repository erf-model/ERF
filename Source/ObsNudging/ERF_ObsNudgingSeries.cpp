/**
 * \file ERF_ObsNudgingSeries.cpp
 */
#include "ERF_ObsNudgingSeries.H"
#include "ERF_NumericalConstants.H"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <sstream>

using namespace amrex;

namespace obs_nudging {

const char*
comp_name (int comp)
{
    switch (comp) {
        case U:     return "u";
        case V:     return "v";
        case W:     return "w";
        case Theta: return "theta";
        default:    return "?";
    }
}

namespace {

// The columns a station file may name
enum Column : int { ColTime = 0, ColZ, ColU, ColV, ColW, ColTheta,
                    ColSU, ColSV, ColSW, ColSTheta, ColSpeed, ColDirection, NColumns };

const char* column_names[NColumns] = {"time", "z", "u", "v", "w", "theta",
                                      "su", "sv", "sw", "stheta", "speed", "direction"};

// Strip a comment and report whether anything is left
bool
strip_comment (std::string& line)
{
    const auto pos = line.find('#');
    if (pos != std::string::npos) { line.erase(pos); }
    return line.find_first_not_of(" \t\r") != std::string::npos;
}

// Parse a whole token as a number.  "nan" parses (to NaN); anything with
// trailing characters does not.
bool
parse_number (const std::string& token, double& value)
{
    const char* begin = token.c_str();
    char* end = nullptr;
    value = std::strtod(begin, &end);
    return end != begin && *end == '\0';
}

std::string
where (const std::string& source, int line_no)
{
    return source + ", line " + std::to_string(line_no) + ": ";
}

} // namespace

bool
parse_station_series (std::istream& is, const std::string& source, Real missing_value,
                      StationSeries& series, std::string& error)
{
    series = StationSeries{};

    int col_of[NColumns];
    std::fill(col_of, col_of + NColumns, -1);
    int ncols = 0;

    const auto is_missing = [missing_value] (double v) {
        return !std::isfinite(v) ||
               std::abs(v - static_cast<double>(missing_value)) <=
               1.0e-6 * std::max(1.0, std::abs(static_cast<double>(missing_value)));
    };

    std::string line;
    int line_no = 0;
    bool have_header = false;
    int nheights = 0;           // heights per record, fixed by the first record
    int ih = 0;                 // height index within the current record

    while (std::getline(is, line))
    {
        ++line_no;
        if (!strip_comment(line)) { continue; }

        std::istringstream ls(line);
        std::vector<std::string> tokens;
        std::string tok;
        while (ls >> tok) { tokens.push_back(tok); }

        if (!have_header)
        {
            for (auto& t : tokens) {
                std::transform(t.begin(), t.end(), t.begin(),
                               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                int found = -1;
                for (int c = 0; c < NColumns; ++c) {
                    if (t == column_names[c]) { found = c; break; }
                }
                if (found < 0) {
                    error = where(source, line_no) + "unknown column '" + t + "'; the columns are "
                            "time z u v w theta su sv sw stheta speed direction";
                    return false;
                }
                if (col_of[found] >= 0) {
                    error = where(source, line_no) + "column '" + t + "' is named twice";
                    return false;
                }
                col_of[found] = ncols++;
            }

            const auto has = [&](int c) { return col_of[c] >= 0; };
            if (!has(ColTime) || !has(ColZ)) {
                error = where(source, line_no) + "the header must name the columns time and z";
                return false;
            }
            if (has(ColU) != has(ColV)) {
                error = where(source, line_no) + "u and v must be given together";
                return false;
            }
            if (has(ColSpeed) != has(ColDirection)) {
                error = where(source, line_no) + "speed and direction must be given together";
                return false;
            }
            if (has(ColU) && has(ColSpeed)) {
                error = where(source, line_no) + "give the wind as u v or as speed direction, not both";
                return false;
            }
            const bool has_wind = has(ColU) || has(ColSpeed);
            if ((has(ColSU) || has(ColSV)) && !has_wind) {
                error = where(source, line_no) + "su/sv need a wind column";
                return false;
            }
            if (has(ColSW) && !has(ColW)) {
                error = where(source, line_no) + "sw needs the column w";
                return false;
            }
            if (has(ColSTheta) && !has(ColTheta)) {
                error = where(source, line_no) + "stheta needs the column theta";
                return false;
            }
            if (!has_wind && !has(ColW) && !has(ColTheta)) {
                error = where(source, line_no) + "the file has none of the nudged quantities "
                        "(u v, speed direction, w, theta)";
                return false;
            }

            series.has[U]     = has_wind;
            series.has[V]     = has_wind;
            series.has[W]     = has(ColW);
            series.has[Theta] = has(ColTheta);
            have_header = true;
            continue;
        }

        if (static_cast<int>(tokens.size()) != ncols) {
            error = where(source, line_no) + "expected " + std::to_string(ncols) +
                    " entries, found " + std::to_string(tokens.size());
            return false;
        }

        double v[NColumns];
        std::fill(v, v + NColumns, std::numeric_limits<double>::quiet_NaN());
        for (int c = 0; c < NColumns; ++c) {
            if (col_of[c] < 0) { continue; }
            if (!parse_number(tokens[col_of[c]], v[c])) {
                error = where(source, line_no) + "'" + tokens[col_of[c]] + "' in column " +
                        column_names[c] + " is not a number";
                return false;
            }
        }

        if (is_missing(v[ColTime]) || is_missing(v[ColZ])) {
            error = where(source, line_no) + "time and z may not be missing";
            return false;
        }

        // A new time starts a new record
        const double time = v[ColTime];
        const Real   z    = static_cast<Real>(v[ColZ]);
        if (series.records.empty() || time != series.records.back().time)
        {
            if (!series.records.empty()) {
                if (time < series.records.back().time) {
                    error = where(source, line_no) + "the times must increase down the file";
                    return false;
                }
                if (ih != nheights) {
                    error = where(source, line_no) + "the record at time " +
                            std::to_string(series.records.back().time) + " has " +
                            std::to_string(ih) + " heights, but the first record has " +
                            std::to_string(nheights) + "; every time needs the same heights";
                    return false;
                }
            }
            StationRecord rec;
            rec.time = time;
            series.records.push_back(rec);
            ih = 0;
        }

        StationRecord& rec = series.records.back();

        if (series.records.size() == 1) {
            if (!series.heights.empty() && z <= series.heights.back()) {
                error = where(source, line_no) + "the heights within a time must increase";
                return false;
            }
            series.heights.push_back(z);
            nheights = static_cast<int>(series.heights.size());
        } else if (ih >= nheights || z != series.heights[ih]) {
            error = where(source, line_no) + "height " + std::to_string(z) + " at time " +
                    std::to_string(time) + " does not match the heights of the first time; "
                    "every time needs the same heights, in the same order";
            return false;
        }

        // The measured values of each component, NaN where missing
        const Real nan = std::numeric_limits<Real>::quiet_NaN();
        Real mean[NComp]  = {nan, nan, nan, nan};
        Real sigma[NComp] = {Real(0.0), Real(0.0), Real(0.0), Real(0.0)};

        if (col_of[ColU] >= 0 && !is_missing(v[ColU]) && !is_missing(v[ColV])) {
            mean[U] = static_cast<Real>(v[ColU]);
            mean[V] = static_cast<Real>(v[ColV]);
        }
        if (col_of[ColSpeed] >= 0 && !is_missing(v[ColSpeed]) && !is_missing(v[ColDirection])) {
            const double speed = v[ColSpeed];
            const double dir   = v[ColDirection];
            if (speed < 0.0) {
                error = where(source, line_no) + "the wind speed may not be negative";
                return false;
            }
            if (dir < 0.0 || dir > 360.0) {
                error = where(source, line_no) + "the wind direction must be in [0, 360] degrees";
                return false;
            }
            // Meteorological convention: the direction the wind blows from,
            // clockwise from north
            const double rad = dir * static_cast<double>(PI) / 180.0;
            mean[U] = static_cast<Real>(-speed * std::sin(rad));
            mean[V] = static_cast<Real>(-speed * std::cos(rad));
        }
        if (col_of[ColW] >= 0 && !is_missing(v[ColW])) {
            mean[W] = static_cast<Real>(v[ColW]);
        }
        if (col_of[ColTheta] >= 0 && !is_missing(v[ColTheta])) {
            if (v[ColTheta] <= 0.0) {
                error = where(source, line_no) + "theta must be positive (it is in kelvin)";
                return false;
            }
            mean[Theta] = static_cast<Real>(v[ColTheta]);
        }

        const int scol[NComp] = {ColSU, ColSV, ColSW, ColSTheta};
        for (int c = 0; c < NComp; ++c) {
            if (col_of[scol[c]] < 0 || is_missing(v[scol[c]])) { continue; }
            if (v[scol[c]] < 0.0) {
                error = where(source, line_no) + "the standard deviation " +
                        column_names[scol[c]] + " may not be negative";
                return false;
            }
            sigma[c] = static_cast<Real>(v[scol[c]]);
        }

        for (int c = 0; c < NComp; ++c) {
            rec.mean[c].push_back(mean[c]);
            rec.sigma[c].push_back(sigma[c]);
        }
        ++ih;
    }

    if (!have_header) {
        error = source + ": the file has no header line naming its columns";
        return false;
    }
    if (series.records.empty()) {
        error = source + ": the file has a header but no measurements";
        return false;
    }
    if (ih != nheights) {
        error = source + ": the last record has " + std::to_string(ih) +
                " heights, but the first record has " + std::to_string(nheights) +
                "; every time needs the same heights";
        return false;
    }
    return true;
}

bool
profiles_at_time (const StationSeries& series, double t,
                  Real cos_alpha, Real sin_alpha,
                  std::array<Profile, NComp>& profiles)
{
    for (auto& p : profiles) { p.clear(); }

    const auto& recs = series.records;
    const int nrec = static_cast<int>(recs.size());
    if (nrec == 0) { return false; }

    // The two records bracketing t and the weight of the later one
    int n0 = 0, n1 = 0;
    double f = 0.0;
    if (nrec > 1) {
        if (t < recs.front().time || t > recs.back().time) { return false; }
        while (n0 + 1 < nrec - 1 && recs[n0+1].time <= t) { ++n0; }
        n1 = n0 + 1;
        f = (t - recs[n0].time) / (recs[n1].time - recs[n0].time);
    }

    const int nh = static_cast<int>(series.heights.size());
    const Real nan = std::numeric_limits<Real>::quiet_NaN();

    // The value of component c at height h at time t, NaN when missing
    auto at_time = [&](int c, int h, Real& mean, Real& sigma) {
        const Real m0 = recs[n0].mean[c][h];
        const Real m1 = recs[n1].mean[c][h];
        const Real s0 = recs[n0].sigma[c][h];
        const Real s1 = recs[n1].sigma[c][h];
        if (f <= 0.0) {
            mean = m0; sigma = s0;
        } else if (f >= 1.0) {
            mean = m1; sigma = s1;
        } else if (std::isfinite(m0) && std::isfinite(m1)) {
            const Real w1 = static_cast<Real>(f);
            const Real w0 = Real(1.0) - w1;
            mean  = w0*m0 + w1*m1;
            sigma = w0*s0 + w1*s1;
        } else {
            mean = nan; sigma = Real(0.0);
        }
    };

    const bool rotate = (cos_alpha != Real(1.0) || sin_alpha != Real(0.0));

    for (int h = 0; h < nh; ++h)
    {
        Real mean[NComp], sigma[NComp];
        for (int c = 0; c < NComp; ++c) { at_time(c, h, mean[c], sigma[c]); }

        if (rotate) {
            if (std::isfinite(mean[U]) && std::isfinite(mean[V])) {
                const Real ug =  cos_alpha*mean[U] + sin_alpha*mean[V];
                const Real vg = -sin_alpha*mean[U] + cos_alpha*mean[V];
                const Real su = std::sqrt(cos_alpha*cos_alpha*sigma[U]*sigma[U] +
                                          sin_alpha*sin_alpha*sigma[V]*sigma[V]);
                const Real sv = std::sqrt(sin_alpha*sin_alpha*sigma[U]*sigma[U] +
                                          cos_alpha*cos_alpha*sigma[V]*sigma[V]);
                mean[U] = ug; mean[V] = vg; sigma[U] = su; sigma[V] = sv;
            } else {
                mean[U] = nan; mean[V] = nan;
            }
        }

        for (int c = 0; c < NComp; ++c) {
            if (!std::isfinite(mean[c])) { continue; }
            profiles[c].z.push_back(series.heights[h]);
            profiles[c].mean.push_back(mean[c]);
            profiles[c].sigma.push_back(sigma[c]);
        }
    }

    return true;
}

} // namespace obs_nudging
