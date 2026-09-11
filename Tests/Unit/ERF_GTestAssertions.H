#ifndef ERF_GTEST_ASSERTIONS_H_
#define ERF_GTEST_ASSERTIONS_H_

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include <gtest/gtest.h>

namespace erf_gtest {

inline const char* precision_mode ()
{
#ifdef AMREX_USE_FLOAT
    return "single";
#else
    return "double";
#endif
}

template <typename T>
::testing::AssertionResult near_with_diagnostics (const char* actual_expression,
                                                  const char* expected_expression,
                                                  const char* tolerance_expression,
                                                  const T& actual,
                                                  const T& expected,
                                                  const T& tolerance)
{
    const double actual_value = static_cast<double>(actual);
    const double expected_value = static_cast<double>(expected);
    const double tolerance_value = static_cast<double>(tolerance);
    const double absolute_error = std::abs(actual_value - expected_value);
    const double scale = std::max({1.0, std::abs(actual_value), std::abs(expected_value)});
    const double relative_error = absolute_error / scale;

    if (std::isfinite(actual_value) && absolute_error <= tolerance_value) {
        return ::testing::AssertionSuccess();
    }

    std::ostringstream message;
    message << std::setprecision(17)
            << "actual (" << actual_expression << ") = " << actual_value << '\n'
            << "expected (" << expected_expression << ") = " << expected_value << '\n'
            << "absolute error = " << absolute_error << '\n'
            << "relative error = " << relative_error << '\n'
            << "allowed tolerance (" << tolerance_expression << ") = " << tolerance_value << '\n'
            << "precision mode = " << precision_mode();
    return ::testing::AssertionFailure() << message.str();
}

template <typename T>
::testing::AssertionResult finite_value (const char* expression, const T& value)
{
    if (std::isfinite(static_cast<double>(value))) {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure()
           << expression << " = " << value << " is not finite (precision mode = "
           << precision_mode() << ")";
}

template <typename T>
::testing::AssertionResult nonnegative_with_tolerance (const char* expression,
                                                        const char* tolerance_expression,
                                                        const T& value,
                                                        const T& tolerance)
{
    if (value >= -tolerance) {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure()
           << expression << " = " << value << " is below zero by more than "
           << tolerance_expression << " = " << tolerance;
}

} // namespace erf_gtest

#define ERF_EXPECT_NEAR(actual, expected, tolerance) \
    EXPECT_PRED_FORMAT3(::erf_gtest::near_with_diagnostics, actual, expected, tolerance)

#define ERF_ASSERT_NEAR(actual, expected, tolerance) \
    ASSERT_PRED_FORMAT3(::erf_gtest::near_with_diagnostics, actual, expected, tolerance)

#define ERF_EXPECT_FINITE(value) \
    EXPECT_PRED_FORMAT1(::erf_gtest::finite_value, value)

#define ERF_EXPECT_NONNEGATIVE(value, tolerance) \
    EXPECT_PRED_FORMAT2(::erf_gtest::nonnegative_with_tolerance, value, tolerance)

#endif
