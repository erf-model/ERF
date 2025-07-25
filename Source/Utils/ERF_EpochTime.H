#ifndef ERF_EPOCH_TIME_H_
#define ERF_EPOCH_TIME_H_

#include <string>
#include <sstream>
#include <iomanip>
#include <ctime>

#ifdef _WIN32
#define timegm _mkgmtime
#endif

AMREX_FORCE_INLINE
std::time_t getEpochTime(const std::string& dateTime, const std::string& dateTimeFormat) {
    std::tm tmTime{};

    std::istringstream ss(dateTime);
    ss >> std::get_time(&tmTime, dateTimeFormat.c_str());

    if (ss.fail()) {
        amrex::Error("Failed to parse datetime");
    }

    // Convert to epoch time (UTC)
    return timegm(&tmTime);
}

AMREX_FORCE_INLINE
std::string
getTimestamp (const amrex::Real epoch_real, const std::string& datetime_format)
{
    auto epoch_nearest_sec = static_cast<std::time_t>(epoch_real);
    std::tm *time_info = std::gmtime(&epoch_nearest_sec);

    char buffer[80];
    std::strftime(buffer, sizeof(buffer), datetime_format.c_str(), time_info);
    std::string str_nearest_sec(buffer);

    double frac_sec = epoch_real - epoch_nearest_sec;
    snprintf(buffer, 80, "%.6f", frac_sec);
    AMREX_ASSERT(buffer[0] == '0');
    std::string str_frac_sec(buffer);

    return str_nearest_sec + str_frac_sec.substr(1);
}


#endif
