#ifndef ERF_PRESCRIBED_SURFACE_FLUX_H_
#define ERF_PRESCRIBED_SURFACE_FLUX_H_

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace prescribed_surface_flux {

struct TimeBracket
{
    std::size_t lower{0};
    std::size_t upper{0};
    double weight{0.0};
};

inline std::string
validate_time_layout (const std::vector<std::string>& dimensions)
{
    return dimensions == std::vector<std::string>{"time"}
             ? std::string{}
             : "time variable must have dimensions (time)";
}

inline std::string
validate_field_layout (const std::string& name,
                       const std::vector<std::string>& dimensions,
                       const std::vector<std::size_t>& shape,
                       const std::vector<std::size_t>& expected_shape)
{
    if (dimensions != std::vector<std::string>{"time", "y", "x"}) {
        return name + " must have dimensions (time,y,x)";
    }
    if (shape != expected_shape) { return name + " has the wrong shape"; }
    return {};
}

inline std::string
validate_time_axis (const std::vector<double>& time)
{
    if (time.empty()) { return "time dimension must contain at least one value"; }
    for (std::size_t n = 0; n < time.size(); ++n) {
        if (!std::isfinite(time[n])) { return "time contains a nonfinite value"; }
        if (n > 0 && !(time[n] > time[n-1])) {
            return "time values must be strictly increasing";
        }
    }
    return {};
}

inline std::string
validate_field (const std::string& name,
                const std::vector<double>& values,
                const std::size_t expected_size,
                const bool require_nonnegative,
                const std::vector<double>& missing_values = {})
{
    if (values.size() != expected_size) {
        return name + " has the wrong number of values";
    }
    for (const double value : values) {
        if (!std::isfinite(value)) { return name + " contains a nonfinite value"; }
        for (const double missing : missing_values) {
            if (value == missing) { return name + " contains missing data"; }
        }
        if (require_nonnegative && value < 0.0) {
            return name + " contains a negative value";
        }
    }
    return {};
}

inline bool
find_time_bracket (const std::vector<double>& time,
                   const double query_time,
                   TimeBracket& bracket,
                   std::string& error)
{
    error = validate_time_axis(time);
    if (!error.empty()) { return false; }
    if (!std::isfinite(query_time)) {
        error = "requested forcing time is nonfinite";
        return false;
    }
    if (query_time < time.front() || query_time > time.back()) {
        error = "requested forcing time is outside file coverage";
        return false;
    }
    if (query_time == time.front()) {
        bracket = {0, 0, 0.0};
        return true;
    }
    if (query_time == time.back()) {
        const std::size_t last = time.size() - 1;
        bracket = {last, last, 0.0};
        return true;
    }
    std::size_t upper = 1;
    while (time[upper] < query_time) { ++upper; }
    if (time[upper] == query_time) {
        bracket = {upper, upper, 0.0};
    } else {
        const std::size_t lower = upper - 1;
        bracket = {lower, upper,
                   (query_time - time[lower]) / (time[upper] - time[lower])};
    }
    return true;
}

inline double
interpolate (const double lower, const double upper, const double weight)
{
    return lower + weight * (upper - lower);
}

} // namespace prescribed_surface_flux

#endif
