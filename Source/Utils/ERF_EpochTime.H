#ifndef ERF_EPOCH_TIME_H_
#define ERF_EPOCH_TIME_H_

#include <string>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <cctype>

#ifdef _WIN32
#define timegm _mkgmtime
#endif

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

AMREX_FORCE_INLINE
std::string
getTimestamp (const amrex::Real epoch_real, const std::string& datetime_format, bool add_long_frac=true)
{
    auto epoch_nearest_sec = static_cast<std::time_t>(epoch_real);
    std::tm *time_info = std::gmtime(&epoch_nearest_sec);

    char buffer[80];
    std::strftime(buffer, sizeof(buffer), datetime_format.c_str(), time_info);
    std::string str_nearest_sec(buffer);

    if (add_long_frac) {
        // We use abs here because the diff would come out as negative
        //    if the times themselves were negative
        double frac_sec = std::abs(epoch_real - epoch_nearest_sec);
        snprintf(buffer, 80, "%.6f", frac_sec);
        AMREX_ASSERT(buffer[0] == '0');
        std::string str_frac_sec(buffer);
        return str_nearest_sec + str_frac_sec.substr(1);
    } else {
        return str_nearest_sec;
    }
}


#endif
