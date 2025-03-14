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
std::time_t
getEpochTime (const std::string& dateTime, const std::string& dateTimeFormat)
{
    // Create a tm object to store the parsed date and time.
    std::tm tmTime;

    // Create a stream which we will use to parse the string,
    // which we provide to constructor of stream to fill the buffer.
    std::istringstream ss{ dateTime };

    // Now we read from buffer using get_time manipulator
    // and formatting the input appropriately.
    ss >> std::get_time(&tmTime, dateTimeFormat.c_str());

    // Convert the tm structure to time_t value and return.
    // Here we use timegm since the output should be relative to UTC.
    auto epoch = timegm(&tmTime);
    // Print() << "Time Stamp: "<< std::put_time(&tmTime, "%c")
    //         << " , Epoch: " << epoch << std::endl;

    return epoch;
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
