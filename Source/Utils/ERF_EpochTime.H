#ifndef ERF_EPOCH_TIME_H_
#define ERF_EPOCH_TIME_H_

#include <string>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <cctype>
#include <cmath>
#include <cstdio>

#ifdef _WIN32
#define timegm _mkgmtime
#endif

/**
 * Parse an integer of fixed width from a string.
 * @param[in,out] dt String to parse from; advanced by width on success.
 * @param[in] width Number of characters to parse.
 * @param[out] value Parsed integer value.
 * @return True if parsing was successful, false otherwise.
 */
AMREX_GPU_HOST AMREX_FORCE_INLINE
bool
parse_fixed_width_int (const char*& dt, int width, int& value)
{
    value = 0;

    for (int idx = 0; idx < width; ++idx) {
        unsigned char ch = static_cast<unsigned char>(dt[idx]);
        if (ch == '\0' || !std::isdigit(ch)) {
            return false;
        }
        value = 10 * value + (dt[idx] - '0');
    }

    dt += width;
    return true;
}

/**
 * Convert a formatted date-time string to epoch time.
 * @param[in] dateTime String containing the date and time.
 * @param[in] dateTimeFormat Format string describing the layout of dateTime.
 * @return Epoch time as std::time_t, or -1 on failure.
 */
AMREX_GPU_HOST AMREX_FORCE_INLINE
std::time_t
getEpochTime (const std::string& dateTime, const std::string& dateTimeFormat)
{
    std::tm tmTime{};
    const char* dt  = dateTime.c_str();
    const char* fmt = dateTimeFormat.c_str();

    // Walk both the format string and the date string together
    while (*fmt && *dt) {
        if (*fmt == '%') {
            ++fmt;
            int parsed_value = 0;
            switch (*fmt) {
            case 'Y': // 4-digit year
                if (!parse_fixed_width_int(dt, 4, parsed_value)) return -1;
                tmTime.tm_year = parsed_value;
                tmTime.tm_year -= 1900;
                break;
            case 'm': // 2-digit month
                if (!parse_fixed_width_int(dt, 2, parsed_value)) return -1;
                tmTime.tm_mon = parsed_value;
                tmTime.tm_mon -= 1; // struct tm expects 0-11
                break;
            case 'd': // 2-digit day
                if (!parse_fixed_width_int(dt, 2, parsed_value)) return -1;
                tmTime.tm_mday = parsed_value;
                break;
            case 'H': // 2-digit hour
                if (!parse_fixed_width_int(dt, 2, parsed_value)) return -1;
                tmTime.tm_hour = parsed_value;
                break;
            case 'M': // 2-digit minute
                if (!parse_fixed_width_int(dt, 2, parsed_value)) return -1;
                tmTime.tm_min = parsed_value;
                break;
            case 'S': // 2-digit second
                if (!parse_fixed_width_int(dt, 2, parsed_value)) return -1;
                tmTime.tm_sec = parsed_value;
                break;
            default:
                // Unsupported format specifier
                return -1;
            }
        } else {
            // Literal character, must match
            if (*fmt != *dt) {
                return -1; // mismatch
            }
            ++dt;
        }
        ++fmt;
    }

    if (*fmt != '\0' || *dt != '\0') {
        return -1;
    }

    if (tmTime.tm_mon < 0 || tmTime.tm_mon > 11 ||
        tmTime.tm_mday < 1 || tmTime.tm_mday > 31 ||
        tmTime.tm_hour < 0 || tmTime.tm_hour > 23 ||
        tmTime.tm_min < 0 || tmTime.tm_min > 59 ||
        tmTime.tm_sec < 0 || tmTime.tm_sec > 60) {
        return -1;
    }

    return timegm(&tmTime);
}

/**
 * Convert epoch time to a formatted timestamp string.
 * @param[in] epoch_real Epoch time in seconds.
 * @param[in] datetime_format Format string for the resulting timestamp.
 * @param[in] add_long_frac Whether to append fractional seconds to the output.
 * @return Formatted timestamp string.
 */
AMREX_FORCE_INLINE
std::string
getTimestamp (const double epoch_real, const std::string& datetime_format, bool add_long_frac=true)
{
    std::time_t epoch_nearest_sec;
    long long   frac_usec = 0;

    if (add_long_frac) {
        //
        // NOTE: when we print the fraction we must derive it and the whole seconds from
        //       the *same* rounded value.  If we instead truncated the seconds and then
        //       printed the leftover fraction to six places, a fraction just below one
        //       would print as "1.000000" next to the un-incremented second and the
        //       timestamp would come out a full second early.
        //
        // We round to the microsecond we are about to print, then split, rounding the
        //    seconds toward -infinity so that the fraction is never negative -- which
        //    also gives the right answer for times before the epoch.
        //
        constexpr long long usec_per_sec = 1000000;
        const long long total_usec = std::llround(epoch_real * 1.e6);
        long long whole_sec = total_usec / usec_per_sec;
        frac_usec           = total_usec % usec_per_sec;
        if (frac_usec < 0) { whole_sec -= 1; frac_usec += usec_per_sec; }
        epoch_nearest_sec = static_cast<std::time_t>(whole_sec);
    } else {
        epoch_nearest_sec = static_cast<std::time_t>(epoch_real);
    }

    std::tm *time_info = std::gmtime(&epoch_nearest_sec);

    char buffer[80];
    std::strftime(buffer, sizeof(buffer), datetime_format.c_str(), time_info);
    std::string str_nearest_sec(buffer);

    if (add_long_frac) {
        snprintf(buffer, 80, ".%06lld", frac_usec);
        return str_nearest_sec + std::string(buffer);
    } else {
        return str_nearest_sec;
    }
}


#endif
