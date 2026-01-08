#ifndef ERF_EPOCH_TIME_H_
#define ERF_EPOCH_TIME_H_

#include <string>
#include <sstream>
#include <iomanip>
#include <ctime>

#ifdef _WIN32
#define timegm _mkgmtime
#endif

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
            switch (*fmt) {
            case 'Y': // 4-digit year
                tmTime.tm_year = (dt[0]-'0')*1000 + (dt[1]-'0')*100
                               + (dt[2]-'0')*10   + (dt[3]-'0');
                tmTime.tm_year -= 1900;
                dt += 4;
                break;
            case 'm': // 2-digit month
                tmTime.tm_mon = (dt[0]-'0')*10 + (dt[1]-'0');
                tmTime.tm_mon -= 1; // struct tm expects 0-11
                dt += 2;
                break;
            case 'd': // 2-digit day
                tmTime.tm_mday = (dt[0]-'0')*10 + (dt[1]-'0');
                dt += 2;
                break;
            case 'H': // 2-digit hour
                tmTime.tm_hour = (dt[0]-'0')*10 + (dt[1]-'0');
                dt += 2;
                break;
            case 'M': // 2-digit minute
                tmTime.tm_min = (dt[0]-'0')*10 + (dt[1]-'0');
                dt += 2;
                break;
            case 'S': // 2-digit second
                tmTime.tm_sec = (dt[0]-'0')*10 + (dt[1]-'0');
                dt += 2;
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
        double frac_sec = epoch_real - epoch_nearest_sec;
        snprintf(buffer, 80, "%.6f", frac_sec);
        AMREX_ASSERT(buffer[0] == '0');
        std::string str_frac_sec(buffer);
        return str_nearest_sec + str_frac_sec.substr(1);
    } else {
        return str_nearest_sec;
    }
}


#endif
