#ifndef ERF_GTEST_PRINTERS_H_
#define ERF_GTEST_PRINTERS_H_

#include <ostream>

namespace erf_gtest {

template <typename T>
inline void print_field (std::ostream* os,
                         const char* name,
                         const T& value,
                         const bool first = false)
{
    if (!first) {
        *os << ", ";
    }
    *os << name << "=" << value;
}

} // namespace erf_gtest

#endif
