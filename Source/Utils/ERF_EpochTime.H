#ifndef ERF_EPOCH_TIME_H_
#define ERF_EPOCH_TIME_H_

#include <string>
#include <ctime>

AMREX_FORCE_INLINE
std::time_t
getEpochTime (const std::string& dateTime, const std::string& dateTimeFormat)
{
    // Create a stream which we will use to parse the string,
    // which we provide to constructor of stream to fill the buffer.
    std::istringstream ss{ dateTime };

    // Create a tm object to store the parsed date and time.
    std::tm tmTime;
    memset(&tmTime, 0, sizeof(tmTime));

    // Now we read from buffer using get_time manipulator
    // and formatting the input appropriately.
    strptime(dateTime.c_str(), dateTimeFormat.c_str(), &tmTime);

    // Convert the tm structure to time_t value and return.
    // Here we use timegm since the output should be relative to UTC.
    auto epoch = timegm(&tmTime);
    // Print() << "Time Stamp: "<< std::put_time(&tmTime, "%c")
    //         << " , Epoch: " << epoch << std::endl;

    return epoch;
}

#endif
