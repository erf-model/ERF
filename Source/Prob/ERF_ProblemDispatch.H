#ifndef ERF_PROBLEM_DISPATCH_H_
#define ERF_PROBLEM_DISPATCH_H_

#include <string>

namespace erf_problem_dispatch {

enum class CustomVelocityInitializer : unsigned char {
    None,
    CloudChamber,
    UserDefined
};

inline CustomVelocityInitializer
custom_velocity_initializer (const std::string& problem_name) noexcept
{
    if (problem_name == "cloud chamber" || problem_name == "cloudchamber") {
        return CustomVelocityInitializer::CloudChamber;
    }
    if (problem_name == "userdefined") {
        return CustomVelocityInitializer::UserDefined;
    }
    return CustomVelocityInitializer::None;
}

} // namespace erf_problem_dispatch

#endif
