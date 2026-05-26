#include <array>

#include <gtest/gtest.h>

#include "ERF_GTestInterpolationCommon.H"

using namespace interpolation_test;

namespace {

struct ReferenceCase {
    AdvType adv_type;
    int sigma;
};

inline void expect_near_relative (const amrex::Real actual,
                                  const amrex::Real expected,
                                  const amrex::Real factor)
{
    EXPECT_NEAR(actual, expected, scaled_tol(actual, expected, factor));
}

std::array<StencilValues, 3> asymmetric_stencils ()
{
    return {{
        {amrex::Real(1.0),  amrex::Real(-2.5), amrex::Real(0.75),  amrex::Real(3.25),  amrex::Real(-1.125), amrex::Real(2.0)},
        {amrex::Real(-4.0), amrex::Real(0.5),  amrex::Real(7.0),   amrex::Real(-3.0),  amrex::Real(1.5),    amrex::Real(2.25)},
        {amrex::Real(2.0),  amrex::Real(2.5),  amrex::Real(-1.0),  amrex::Real(4.5),   amrex::Real(-3.75),  amrex::Real(0.125)}}};
}

std::string trace_label (const AdvType adv_type,
                         const int sigma,
                         const int sample_index,
                         const amrex::Real actual,
                         const amrex::Real expected,
                         const amrex::Real factor)
{
    return "scheme=" + scheme_name(adv_type) +
           " sigma=" + std::to_string(sigma) +
           " sample=" + std::to_string(sample_index) +
           " actual=" + std::to_string(static_cast<double>(actual)) +
           " expected=" + std::to_string(static_cast<double>(expected)) +
           " normalized_error=" + std::to_string(static_cast<double>(normalized_error(actual, expected, factor)));
}

} // namespace

// Motivation: The reference stencils must satisfy the independently derived
// finite-volume exactness orders before they are used in device-path comparisons.
TEST(InterpolationScalar, ReferencePolynomialExactnessFiniteVolume)
{
    const std::array<ReferenceCase, 7> cases = {{
        {AdvType::Centered_2nd, 0},
        {AdvType::Upwind_3rd, -1},
        {AdvType::Upwind_3rd, 1},
        {AdvType::Centered_4th, 0},
        {AdvType::Upwind_5th, -1},
        {AdvType::Upwind_5th, 1},
        {AdvType::Centered_6th, 0}}};

    for (const ReferenceCase& test_case : cases) {
        for (int degree = 0; degree <= exactness_degree(test_case.adv_type); ++degree) {
            const StencilValues stencil = make_monomial_stencil(degree);
            const amrex::Real actual = reference_reconstruction(stencil, test_case.adv_type, test_case.sigma);
            const amrex::Real expected = (degree == 0) ? amrex::Real(1.0) : amrex::Real(0.0);

            SCOPED_TRACE(trace_label(test_case.adv_type, test_case.sigma, degree,
                                     actual, expected, kPolynomialRelTol));
            expect_near_relative(actual, expected, kPolynomialRelTol);
        }
    }
}

// Motivation: Leading-error constants lock stencil orientation and distinguish
// 2nd/3rd/4th/5th/6th order behavior.
TEST(InterpolationScalar, ReferenceLeadingErrorCoefficients)
{
    const std::array<ReferenceCase, 7> cases = {{
        {AdvType::Centered_2nd, 0},
        {AdvType::Upwind_3rd, -1},
        {AdvType::Upwind_3rd, 1},
        {AdvType::Centered_4th, 0},
        {AdvType::Upwind_5th, -1},
        {AdvType::Upwind_5th, 1},
        {AdvType::Centered_6th, 0}}};

    for (const ReferenceCase& test_case : cases) {
        const int degree = leading_error_degree(test_case.adv_type);
        const StencilValues stencil = make_monomial_stencil(degree);
        const amrex::Real actual = reference_reconstruction(stencil, test_case.adv_type, test_case.sigma);
        const amrex::Real expected = leading_error_expected(test_case.adv_type, test_case.sigma);

        SCOPED_TRACE(trace_label(test_case.adv_type, test_case.sigma, degree,
                                 actual, expected, kPolynomialRelTol));
        expect_near_relative(actual, expected, kPolynomialRelTol);
    }
}

// Motivation: A zero velocity removes the odd upwind correction, so upwind 3/5
// reduce to centered 4/6.
TEST(InterpolationScalar, ZeroUpwindReferenceReducesToCentered)
{
    const auto stencils = asymmetric_stencils();

    for (int sample_index = 0; sample_index < static_cast<int>(stencils.size()); ++sample_index) {
        const StencilValues& stencil = stencils[sample_index];

        const amrex::Real upwind3_zero = reference_reconstruction(stencil, AdvType::Upwind_3rd, 0);
        const amrex::Real centered4 = reference_reconstruction(stencil, AdvType::Centered_4th, 0);
        SCOPED_TRACE(trace_label(AdvType::Upwind_3rd, 0, sample_index,
                                 upwind3_zero, centered4, kValueRelTol));
        expect_near_relative(upwind3_zero, centered4, kValueRelTol);

        const amrex::Real upwind5_zero = reference_reconstruction(stencil, AdvType::Upwind_5th, 0);
        const amrex::Real centered6 = reference_reconstruction(stencil, AdvType::Centered_6th, 0);
        SCOPED_TRACE(trace_label(AdvType::Upwind_5th, 0, sample_index,
                                 upwind5_zero, centered6, kValueRelTol));
        expect_near_relative(upwind5_zero, centered6, kValueRelTol);
    }
}

// Motivation: Flipping the upwind sign and mirroring the stencil must give the
// same face value.
TEST(InterpolationScalar, PositiveNegativeUpwindReferenceMirrorSymmetry)
{
    const auto stencils = asymmetric_stencils();
    const std::array<AdvType, 2> upwind_schemes = {{AdvType::Upwind_3rd, AdvType::Upwind_5th}};

    for (const AdvType adv_type : upwind_schemes) {
        for (int sample_index = 0; sample_index < static_cast<int>(stencils.size()); ++sample_index) {
            const StencilValues& stencil = stencils[sample_index];
            const amrex::Real positive = reference_reconstruction(stencil, adv_type, 1);
            const amrex::Real negative_mirror =
                reference_reconstruction(mirrored_stencil(stencil), adv_type, -1);

            SCOPED_TRACE(trace_label(adv_type, 1, sample_index,
                                     positive, negative_mirror, kValueRelTol));
            expect_near_relative(positive, negative_mirror, kValueRelTol);
        }
    }
}

// Motivation: interpolatedVal is a non-WENO centered/upwind helper; this
// documents the invalid-domain guard without relying on brittle AMReX death tests.
TEST(InterpolationScalar, SupportedAdvTypePredicateClassifiesOnlyNonWenoSchemes)
{
    for (const AdvType adv_type : kSupportedAdvTypes) {
        EXPECT_TRUE(isSupportedInterpolationAdvType(adv_type)) << scheme_name(adv_type);
    }

    for (const AdvType adv_type : kUnsupportedAdvTypes) {
        EXPECT_FALSE(isSupportedInterpolationAdvType(adv_type)) << scheme_name(adv_type);
    }
}