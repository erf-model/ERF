#include <array>

#include <AMReX_FArrayBox.H>
#include <AMReX_GpuContainers.H>

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

namespace {

std::array<CenteredFaceStencil7, 3> asymmetric_centered_face_stencils ()
{
    return {{
        {amrex::Real(-2.75), amrex::Real(1.5),  amrex::Real(-0.25), amrex::Real(3.0),
         amrex::Real(-1.25), amrex::Real(2.0), amrex::Real(0.75)},
        {amrex::Real(1.25),  amrex::Real(-3.5), amrex::Real(2.75),  amrex::Real(-0.5),
         amrex::Real(4.0),   amrex::Real(-1.5), amrex::Real(2.25)},
        {amrex::Real(0.5),   amrex::Real(2.5),  amrex::Real(-4.0),  amrex::Real(1.75),
         amrex::Real(-2.25), amrex::Real(3.5),  amrex::Real(-0.75)}}};
}

CenteredFaceStencil7 mirror_face_plus_half_stencil7 (const CenteredFaceStencil7& stencil)
{
    return CenteredFaceStencil7{
        amrex::Real(0.0),
        stencil.q3,
        stencil.q2,
        stencil.q1,
        stencil.q0,
        stencil.qm1,
        stencil.qm2};
}

inline amrex::Box upwind_direct_box ()
{
    return amrex::Box(amrex::IntVect(-2, 0, 0), amrex::IntVect(3, 0, 0));
}

void fill_upwind_direct_box (amrex::FArrayBox& qty,
                             const CenteredFaceStencil7& stencil)
{
    auto qty_arr = qty.array();
    qty_arr(-2, 0, 0, 0) = stencil.qm2;
    qty_arr(-1, 0, 0, 0) = stencil.qm1;
    qty_arr( 0, 0, 0, 0) = stencil.q0;
    qty_arr( 1, 0, 0, 0) = stencil.q1;
    qty_arr( 2, 0, 0, 0) = stencil.q2;
    qty_arr( 3, 0, 0, 0) = stencil.q3;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real normalized_upwind_alpha (const amrex::Real upw,
                                     const amrex::Real upw_frac) noexcept
{
    return amrex::Real(signum(upw)) * upw_frac;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real reference_upwind3_blended (const CenteredFaceStencil7& stencil,
                                       const amrex::Real upw,
                                       const amrex::Real upw_frac) noexcept
{
    const amrex::Real alpha = normalized_upwind_alpha(upw, upw_frac);
    return ((-amrex::Real(1.0) + alpha) / amrex::Real(12.0)) * stencil.q2 +
           (( amrex::Real(7.0) - amrex::Real(3.0) * alpha) / amrex::Real(12.0)) * stencil.q1 +
           (( amrex::Real(7.0) + amrex::Real(3.0) * alpha) / amrex::Real(12.0)) * stencil.q0 +
           ((-amrex::Real(1.0) - alpha) / amrex::Real(12.0)) * stencil.qm1;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real reference_upwind5_blended (const CenteredFaceStencil7& stencil,
                                       const amrex::Real upw,
                                       const amrex::Real upw_frac) noexcept
{
    const amrex::Real alpha = normalized_upwind_alpha(upw, upw_frac);
    return (( amrex::Real(1.0) - alpha) / amrex::Real(60.0)) * stencil.q3 +
           ((-amrex::Real(8.0) + amrex::Real(5.0) * alpha) / amrex::Real(60.0)) * stencil.q2 +
           (( amrex::Real(37.0) - amrex::Real(10.0) * alpha) / amrex::Real(60.0)) * stencil.q1 +
           (( amrex::Real(37.0) + amrex::Real(10.0) * alpha) / amrex::Real(60.0)) * stencil.q0 +
           ((-amrex::Real(8.0) - amrex::Real(5.0) * alpha) / amrex::Real(60.0)) * stencil.qm1 +
           (( amrex::Real(1.0) + alpha) / amrex::Real(60.0)) * stencil.qm2;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real reference_upwind_blended (const CenteredFaceStencil7& stencil,
                                      const AdvType adv_type,
                                      const amrex::Real upw,
                                      const amrex::Real upw_frac) noexcept
{
    return (adv_type == AdvType::Upwind_3rd)
        ? reference_upwind3_blended(stencil, upw, upw_frac)
        : reference_upwind5_blended(stencil, upw, upw_frac);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real evaluate_direct_upwind_policy_device (const amrex::Array4<const amrex::Real>& qty,
                                                  const AdvType adv_type,
                                                  const amrex::Real upw,
                                                  const amrex::Real upw_frac) noexcept
{
    amrex::Real value = amrex::Real(0.0);
    if (adv_type == AdvType::Upwind_3rd) {
        UPWIND3 policy(qty, upw_frac);
        policy.InterpolateInX(1, 0, 0, kQtyComp, value, upw);
    } else {
        UPWIND5 policy(qty, upw_frac);
        policy.InterpolateInX(1, 0, 0, kQtyComp, value, upw);
    }
    return value;
}

amrex::Real evaluate_direct_upwind_policy_on_device (const CenteredFaceStencil7& stencil,
                                                     const AdvType adv_type,
                                                     const amrex::Real upw,
                                                     const amrex::Real upw_frac)
{
    amrex::FArrayBox qty(upwind_direct_box(), 1);
    fill_upwind_direct_box(qty, stencil);
    const auto qty_arr = qty.const_array();
    amrex::Gpu::DeviceVector<amrex::Real> device_result(1, amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real> host_result(1, amrex::Real(0.0));
    amrex::Real* result_ptr = device_result.data();

    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int) noexcept {
        result_ptr[0] = evaluate_direct_upwind_policy_device(qty_arr, adv_type, upw, upw_frac);
    });
    gpu_sync();

    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_result.begin(), device_result.end(), host_result.begin());
    return host_result[0];
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real evaluate_upwindall_device (const amrex::Array4<const amrex::Real>& qty,
                                       const AdvType adv_type,
                                       const amrex::Real upw,
                                       const amrex::Real upw_frac) noexcept
{
    UPWINDALL policy(qty, upw_frac);
    amrex::Real value = amrex::Real(0.0);
    policy.InterpolateInZ(0, 0, 1, kQtyComp, value, upw, adv_type);
    return value;
}

amrex::Real evaluate_upwindall_with_fraction_on_device (const CenteredFaceStencil7& stencil,
                                                        const AdvType adv_type,
                                                        const amrex::Real upw,
                                                        const amrex::Real upw_frac)
{
    amrex::FArrayBox qty(amrex::Box(amrex::IntVect(0, 0, -2), amrex::IntVect(0, 0, 3)), 1);
    auto qty_arr = qty.array();
    qty_arr(0, 0, -2, 0) = stencil.qm2;
    qty_arr(0, 0, -1, 0) = stencil.qm1;
    qty_arr(0, 0,  0, 0) = stencil.q0;
    qty_arr(0, 0,  1, 0) = stencil.q1;
    qty_arr(0, 0,  2, 0) = stencil.q2;
    qty_arr(0, 0,  3, 0) = stencil.q3;

    const auto qty_const_arr = qty.const_array();
    amrex::Gpu::DeviceVector<amrex::Real> device_result(1, amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real> host_result(1, amrex::Real(0.0));
    amrex::Real* result_ptr = device_result.data();

    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int) noexcept {
        result_ptr[0] = evaluate_upwindall_device(qty_const_arr, adv_type, upw, upw_frac);
    });
    gpu_sync();

    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_result.begin(), device_result.end(), host_result.begin());
    return host_result[0];
}

std::array<amrex::Real, 2> evaluate_direct_policy_before_after_set_upwinding_on_device (
    const CenteredFaceStencil7& stencil,
    const AdvType adv_type,
    const amrex::Real upw,
    const amrex::Real initial_upw_frac,
    const amrex::Real updated_upw_frac)
{
    amrex::FArrayBox qty(upwind_direct_box(), 1);
    fill_upwind_direct_box(qty, stencil);
    const auto qty_arr = qty.const_array();
    amrex::Gpu::DeviceVector<amrex::Real> device_result(2, amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real> host_result(2, amrex::Real(0.0));
    amrex::Real* result_ptr = device_result.data();

    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int) noexcept {
        amrex::Real before = amrex::Real(0.0);
        amrex::Real after = amrex::Real(0.0);

        if (adv_type == AdvType::Upwind_3rd) {
            UPWIND3 policy(qty_arr, initial_upw_frac);
            policy.InterpolateInX(1, 0, 0, kQtyComp, before, upw);
            policy.SetUpwinding(updated_upw_frac);
            policy.InterpolateInX(1, 0, 0, kQtyComp, after, upw);
        } else {
            UPWIND5 policy(qty_arr, initial_upw_frac);
            policy.InterpolateInX(1, 0, 0, kQtyComp, before, upw);
            policy.SetUpwinding(updated_upw_frac);
            policy.InterpolateInX(1, 0, 0, kQtyComp, after, upw);
        }

        result_ptr[0] = before;
        result_ptr[1] = after;
    });
    gpu_sync();

    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_result.begin(), device_result.end(), host_result.begin());
    return {{host_result[0], host_result[1]}};
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real evaluate_direct_policy_device (const CenteredFaceStencil7& stencil,
                                           const AdvType adv_type,
                                           const int sigma,
                                           const amrex::Array4<const amrex::Real>& dummy_arr) noexcept
{
    if (adv_type == AdvType::Centered_2nd) {
        CENTERED2 policy(dummy_arr, amrex::Real(1.0));
        return policy.Evaluate(stencil.q1, stencil.q0);
    }

    if (adv_type == AdvType::Upwind_3rd) {
        UPWIND3 policy(dummy_arr, amrex::Real(1.0));
        return policy.Evaluate(stencil.q2, stencil.q1, stencil.q0, stencil.qm1, amrex::Real(sigma));
    }

    if (adv_type == AdvType::Centered_4th) {
        CENTERED4 policy(dummy_arr, amrex::Real(1.0));
        return policy.Evaluate(stencil.q2, stencil.q1, stencil.q0, stencil.qm1);
    }

    if (adv_type == AdvType::Upwind_5th) {
        UPWIND5 policy(dummy_arr, amrex::Real(1.0));
        return policy.Evaluate(stencil.q3, stencil.q2, stencil.q1, stencil.q0,
                               stencil.qm1, stencil.qm2, amrex::Real(sigma));
    }

    if (adv_type == AdvType::Centered_6th) {
        CENTERED6 policy(dummy_arr, amrex::Real(1.0));
        return policy.Evaluate(stencil.q3, stencil.q2, stencil.q1, stencil.q0,
                               stencil.qm1, stencil.qm2);
    }

    return amrex::Real(0.0);
}

amrex::Real evaluate_direct_policy_on_device (const CenteredFaceStencil7& stencil,
                                              const AdvType adv_type,
                                              const int sigma)
{
    amrex::FArrayBox dummy(stencil_box(), 1);
    const auto dummy_arr = dummy.const_array();
    amrex::Gpu::DeviceVector<amrex::Real> device_result(1, amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real> host_result(1, amrex::Real(0.0));
    amrex::Real* result_ptr = device_result.data();

    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int) noexcept {
        result_ptr[0] = evaluate_direct_policy_device(stencil, adv_type, sigma, dummy_arr);
    });
    gpu_sync();

    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_result.begin(), device_result.end(), host_result.begin());
    return host_result[0];
}

amrex::Real evaluate_upwindall_on_device (const CenteredFaceStencil7& stencil,
                                          const AdvType adv_type,
                                          const int sigma)
{
    amrex::FArrayBox dummy(stencil_box(), 1);
    const auto dummy_arr = dummy.const_array();
    amrex::Gpu::DeviceVector<amrex::Real> device_result(1, amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real> host_result(1, amrex::Real(0.0));
    amrex::Real* result_ptr = device_result.data();

    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int) noexcept {
        UPWINDALL policy(dummy_arr, amrex::Real(1.0));
        result_ptr[0] = policy.Evaluate(stencil.q3, stencil.q2, stencil.q1, stencil.q0,
                                        stencil.qm1, stencil.qm2, amrex::Real(sigma), adv_type);
    });
    gpu_sync();

    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_result.begin(), device_result.end(), host_result.begin());
    return host_result[0];
}

amrex::Real evaluate_upwind3sl_on_device (const amrex::Real s,
                                          const amrex::Real sm1,
                                          const amrex::Real sm2)
{
    amrex::FArrayBox dummy(stencil_box(), 1);
    const auto dummy_arr = dummy.const_array();
    amrex::Gpu::DeviceVector<amrex::Real> device_result(1, amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real> host_result(1, amrex::Real(0.0));
    amrex::Real* result_ptr = device_result.data();

    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int) noexcept {
        UPWIND3SL policy(dummy_arr, amrex::Real(1.0));
        result_ptr[0] = policy.Evaluate(s, sm1, sm2);
    });
    gpu_sync();

    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_result.begin(), device_result.end(), host_result.begin());
    return host_result[0];
}

std::string centered_face_trace (const AdvType adv_type,
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

std::string blended_trace (const AdvType adv_type,
                           const amrex::Real upw,
                           const amrex::Real upw_frac,
                           const int sample_index,
                           const amrex::Real actual,
                           const amrex::Real expected,
                           const amrex::Real factor)
{
    return "scheme=" + scheme_name(adv_type) +
           " upw=" + std::to_string(static_cast<double>(upw)) +
           " upw_frac=" + std::to_string(static_cast<double>(upw_frac)) +
           " sample=" + std::to_string(sample_index) +
           " actual=" + std::to_string(static_cast<double>(actual)) +
           " expected=" + std::to_string(static_cast<double>(expected)) +
           " normalized_error=" + std::to_string(static_cast<double>(normalized_error(actual, expected, factor)));
}

} // namespace

// Motivation: The direct centered/upwind policy helpers are the lowest-level
// finite-volume reconstructions. Monomial cell averages through each design
// degree must reconstruct the analytic face value at x = +1/2.
TEST(InterpolationScalar, DirectPolicyEvaluateExactnessMatchesFiniteVolumeMonomials)
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
            const CenteredFaceStencil7 stencil = make_centered_monomial_stencil7(degree);
            const amrex::Real actual = evaluate_direct_policy_on_device(stencil, test_case.adv_type, test_case.sigma);
            const amrex::Real expected = centered_face_value_monomial(degree, kFacePlusHalf);

            SCOPED_TRACE(centered_face_trace(test_case.adv_type, test_case.sigma, degree,
                                            actual, expected, kPolynomialRelTol));
            expect_near_relative(actual, expected, kPolynomialRelTol);
        }
    }
}

// Motivation: Zero upwind removes the odd correction in the direct policy
// helpers, so UPWIND3/5 must collapse to CENTERED4/6 on the same stencil.
TEST(InterpolationScalar, DirectUpwindPoliciesReduceToCenteredAtZeroUpwind)
{
    const auto stencils = asymmetric_centered_face_stencils();

    for (int sample_index = 0; sample_index < static_cast<int>(stencils.size()); ++sample_index) {
        const CenteredFaceStencil7& stencil = stencils[sample_index];

        const amrex::Real upwind3_zero = evaluate_direct_policy_on_device(stencil, AdvType::Upwind_3rd, 0);
        const amrex::Real centered4 = evaluate_direct_policy_on_device(stencil, AdvType::Centered_4th, 0);
        SCOPED_TRACE(centered_face_trace(AdvType::Upwind_3rd, 0, sample_index,
                                         upwind3_zero, centered4, kValueRelTol));
        expect_near_relative(upwind3_zero, centered4, kValueRelTol);

        const amrex::Real upwind5_zero = evaluate_direct_policy_on_device(stencil, AdvType::Upwind_5th, 0);
        const amrex::Real centered6 = evaluate_direct_policy_on_device(stencil, AdvType::Centered_6th, 0);
        SCOPED_TRACE(centered_face_trace(AdvType::Upwind_5th, 0, sample_index,
                                         upwind5_zero, centered6, kValueRelTol));
        expect_near_relative(upwind5_zero, centered6, kValueRelTol);
    }
}

// Motivation: Flipping the upwind sign and reflecting the centered face stencil
// must reconstruct the same face value.
TEST(InterpolationScalar, DirectUpwindPoliciesRespectMirrorSymmetry)
{
    const auto stencils = asymmetric_centered_face_stencils();
    const std::array<AdvType, 2> schemes = {{AdvType::Upwind_3rd, AdvType::Upwind_5th}};

    for (const AdvType adv_type : schemes) {
        for (int sample_index = 0; sample_index < static_cast<int>(stencils.size()); ++sample_index) {
            const CenteredFaceStencil7& stencil = stencils[sample_index];
            const CenteredFaceStencil7 mirrored = mirror_face_plus_half_stencil7(stencil);

            const amrex::Real positive = evaluate_direct_policy_on_device(stencil, adv_type, 1);
            const amrex::Real negative_mirror = evaluate_direct_policy_on_device(mirrored, adv_type, -1);

            SCOPED_TRACE(centered_face_trace(adv_type, 1, sample_index,
                                            positive, negative_mirror, kValueRelTol));
            expect_near_relative(positive, negative_mirror, kValueRelTol);
        }
    }
}

// Motivation: The slope-limited UPWIND3 helper is a four-branch limiter. Each
// representative ratio below protects a distinct limiter branch away from the
// branch boundaries.
TEST(InterpolationScalar, Upwind3SLLimiterBranchTable)
{
    const amrex::Real sm2 = amrex::Real(-1.0);
    const amrex::Real sm1 = amrex::Real(3.0);
    const std::array<std::pair<amrex::Real, amrex::Real>, 4> cases = {{
        {amrex::Real(1.0),  amrex::Real(0.0)},
        {amrex::Real(3.5),  amrex::Real(0.25)},
        {amrex::Real(7.0),  amrex::Real(1.0)},
        {amrex::Real(15.0), amrex::Real(2.0)}}};

    for (int sample_index = 0; sample_index < static_cast<int>(cases.size()); ++sample_index) {
        const amrex::Real s = cases[sample_index].first;
        const amrex::Real expected_phi = cases[sample_index].second;
        const amrex::Real actual = evaluate_upwind3sl_on_device(s, sm1, sm2);
        const amrex::Real expected = sm1 + myhalf * expected_phi * (sm1 - sm2);

        EXPECT_NEAR(actual, expected, scaled_tol(actual, expected, kValueRelTol))
            << "sample=" << sample_index;
    }
}

// Motivation: UPWINDALL is the consolidated centered/upwind helper. It must
// agree with the dedicated helper for every supported AdvType and upwind sign.
TEST(InterpolationScalar, UpwindAllMatchesDedicatedPolicyEvaluate)
{
    const auto stencils = asymmetric_centered_face_stencils();

    for (const AdvType adv_type : kSupportedAdvTypes) {
        for (const int sigma : kRepresentativeUpwindSigns) {
            for (int sample_index = 0; sample_index < static_cast<int>(stencils.size()); ++sample_index) {
                const CenteredFaceStencil7& stencil = stencils[sample_index];
                const amrex::Real actual = evaluate_upwindall_on_device(stencil, adv_type, sigma);
                const amrex::Real expected = evaluate_direct_policy_on_device(stencil, adv_type, sigma);

                SCOPED_TRACE(centered_face_trace(adv_type, sigma, sample_index,
                                                actual, expected, kIdentityRelTol));
                expect_near_relative(actual, expected, kIdentityRelTol);
            }
        }
    }
}

// Motivation: UPWIND3/5 blend the odd upwind correction by m_upw_frac after
// normalizing the velocity sign. Fractional blending is state carried by the
// policy object and is not covered by full-upwind exactness tests.
TEST(InterpolationScalar, DirectUpwindPoliciesHonorUpwindFraction)
{
    const auto stencils = asymmetric_centered_face_stencils();
    const std::array<AdvType, 2> schemes = {{AdvType::Upwind_3rd, AdvType::Upwind_5th}};
    const std::array<amrex::Real, 4> fractions = {{amrex::Real(0.0), amrex::Real(0.25), amrex::Real(0.5), amrex::Real(1.0)}};
    const std::array<amrex::Real, 2> upwinds = {{amrex::Real(1.0), -amrex::Real(1.0)}};

    for (const AdvType adv_type : schemes) {
        for (const amrex::Real upw_frac : fractions) {
            for (const amrex::Real upw : upwinds) {
                for (int sample_index = 0; sample_index < static_cast<int>(stencils.size()); ++sample_index) {
                    const CenteredFaceStencil7& stencil = stencils[sample_index];
                    const amrex::Real actual =
                        evaluate_direct_upwind_policy_on_device(stencil, adv_type, upw, upw_frac);
                    const amrex::Real expected =
                        reference_upwind_blended(stencil, adv_type, upw, upw_frac);

                    SCOPED_TRACE(blended_trace(adv_type, upw, upw_frac, sample_index,
                                               actual, expected, kValueRelTol));
                    expect_near_relative(actual, expected, kValueRelTol);
                }
            }
        }
    }
}

// Motivation: UPWINDALL stores the same fractional blend state as the dedicated
// UPWIND3/5 policies and must reproduce the same independent reference values.
TEST(InterpolationScalar, UpwindAllHonorsUpwindFraction)
{
    const auto stencils = asymmetric_centered_face_stencils();
    const std::array<AdvType, 2> schemes = {{AdvType::Upwind_3rd, AdvType::Upwind_5th}};
    const std::array<amrex::Real, 4> fractions = {{amrex::Real(0.0), amrex::Real(0.25), amrex::Real(0.5), amrex::Real(1.0)}};
    const std::array<amrex::Real, 2> upwinds = {{amrex::Real(1.0), -amrex::Real(1.0)}};

    for (const AdvType adv_type : schemes) {
        for (const amrex::Real upw_frac : fractions) {
            for (const amrex::Real upw : upwinds) {
                for (int sample_index = 0; sample_index < static_cast<int>(stencils.size()); ++sample_index) {
                    const CenteredFaceStencil7& stencil = stencils[sample_index];
                    const amrex::Real actual =
                        evaluate_upwindall_with_fraction_on_device(stencil, adv_type, upw, upw_frac);
                    const amrex::Real expected =
                        reference_upwind_blended(stencil, adv_type, upw, upw_frac);

                    SCOPED_TRACE(blended_trace(adv_type, upw, upw_frac, sample_index,
                                               actual, expected, kValueRelTol));
                    expect_near_relative(actual, expected, kValueRelTol);
                }
            }
        }
    }
}

// Motivation: The policy wrappers normalize every nonzero upwind input to its
// sign before reconstruction. This protects the branch for arbitrary velocity
// magnitudes instead of only testing -1, 0, and +1.
TEST(InterpolationScalar, DirectUpwindPoliciesNormalizeNonzeroUpwindMagnitude)
{
    const auto stencils = asymmetric_centered_face_stencils();
    const std::array<AdvType, 2> schemes = {{AdvType::Upwind_3rd, AdvType::Upwind_5th}};
    const std::array<std::pair<amrex::Real, amrex::Real>, 4> cases = {{
        {amrex::Real(7.3), amrex::Real(1.0)},
        {-amrex::Real(2.0), -amrex::Real(1.0)},
        {amrex::Real(1.0e-12), amrex::Real(1.0)},
        {-amrex::Real(1.0e-12), -amrex::Real(1.0)}}};

    for (const AdvType adv_type : schemes) {
        for (int sample_index = 0; sample_index < static_cast<int>(stencils.size()); ++sample_index) {
            const CenteredFaceStencil7& stencil = stencils[sample_index];
            for (const auto& test_case : cases) {
                const amrex::Real actual =
                    evaluate_direct_upwind_policy_on_device(stencil, adv_type, test_case.first, amrex::Real(1.0));
                const amrex::Real expected =
                    evaluate_direct_upwind_policy_on_device(stencil, adv_type, test_case.second, amrex::Real(1.0));

                SCOPED_TRACE(blended_trace(adv_type, test_case.first, amrex::Real(1.0), sample_index,
                                           actual, expected, kIdentityRelTol));
                expect_near_relative(actual, expected, kIdentityRelTol);
            }
        }
    }
}

// Motivation: UPWINDALL must normalize nonzero velocity magnitudes the same way
// as the dedicated UPWIND3/5 policies before applying the stored blend fraction.
TEST(InterpolationScalar, UpwindAllNormalizesNonzeroUpwindMagnitude)
{
    const auto stencils = asymmetric_centered_face_stencils();
    const std::array<AdvType, 2> schemes = {{AdvType::Upwind_3rd, AdvType::Upwind_5th}};
    const std::array<std::pair<amrex::Real, amrex::Real>, 4> cases = {{
        {amrex::Real(7.3), amrex::Real(1.0)},
        {-amrex::Real(2.0), -amrex::Real(1.0)},
        {amrex::Real(1.0e-12), amrex::Real(1.0)},
        {-amrex::Real(1.0e-12), -amrex::Real(1.0)}}};

    for (const AdvType adv_type : schemes) {
        for (int sample_index = 0; sample_index < static_cast<int>(stencils.size()); ++sample_index) {
            const CenteredFaceStencil7& stencil = stencils[sample_index];
            for (const auto& test_case : cases) {
                const amrex::Real actual =
                    evaluate_upwindall_with_fraction_on_device(stencil, adv_type, test_case.first, amrex::Real(1.0));
                const amrex::Real expected =
                    evaluate_upwindall_with_fraction_on_device(stencil, adv_type, test_case.second, amrex::Real(1.0));

                SCOPED_TRACE(blended_trace(adv_type, test_case.first, amrex::Real(1.0), sample_index,
                                           actual, expected, kIdentityRelTol));
                expect_near_relative(actual, expected, kIdentityRelTol);
            }
        }
    }
}

// Motivation: GetUpwindCellNumber documents the required one-sided ghost depth
// for each policy. This test ties the advertised stencil radius to the
// reconstruction family.
TEST(InterpolationScalar, PolicyStencilRadiusMatchesGetUpwindCellNumber)
{
    amrex::FArrayBox dummy(stencil_box(), 1);
    const auto dummy_arr = dummy.const_array();

    CENTERED2 centered2(dummy_arr, amrex::Real(1.0));
    UPWIND3 upwind3(dummy_arr, amrex::Real(1.0));
    UPWIND3SL upwind3sl(dummy_arr, amrex::Real(1.0));
    CENTERED4 centered4(dummy_arr, amrex::Real(1.0));
    UPWIND5 upwind5(dummy_arr, amrex::Real(1.0));
    CENTERED6 centered6(dummy_arr, amrex::Real(1.0));

    EXPECT_EQ(centered2.GetUpwindCellNumber(), 1);
    EXPECT_EQ(upwind3.GetUpwindCellNumber(), 2);
    EXPECT_EQ(upwind3sl.GetUpwindCellNumber(), 2);
    EXPECT_EQ(centered4.GetUpwindCellNumber(), 2);
    EXPECT_EQ(upwind5.GetUpwindCellNumber(), 3);
    EXPECT_EQ(centered6.GetUpwindCellNumber(), 3);
}

// Motivation: UPWIND3 and UPWIND5 carry the upwind fraction as mutable policy
// state. SetUpwinding must update that state without reconstructing the policy
// object.
TEST(InterpolationScalar, SetUpwindingUpdatesPolicyBlend)
{
    const auto stencils = asymmetric_centered_face_stencils();
    const std::array<AdvType, 2> schemes = {{AdvType::Upwind_3rd, AdvType::Upwind_5th}};

    for (const AdvType adv_type : schemes) {
        for (int sample_index = 0; sample_index < static_cast<int>(stencils.size()); ++sample_index) {
            const CenteredFaceStencil7& stencil = stencils[sample_index];
            const std::array<amrex::Real, 2> actual =
                evaluate_direct_policy_before_after_set_upwinding_on_device(
                    stencil, adv_type, amrex::Real(1.0), amrex::Real(1.0), amrex::Real(0.25));

            const amrex::Real expected_before = reference_upwind_blended(stencil, adv_type, amrex::Real(1.0), amrex::Real(1.0));
            const amrex::Real expected_after = reference_upwind_blended(stencil, adv_type, amrex::Real(1.0), amrex::Real(0.25));

            SCOPED_TRACE(blended_trace(adv_type, amrex::Real(1.0), amrex::Real(1.0), sample_index,
                                       actual[0], expected_before, kValueRelTol));
            expect_near_relative(actual[0], expected_before, kValueRelTol);

            SCOPED_TRACE(blended_trace(adv_type, amrex::Real(1.0), amrex::Real(0.25), sample_index,
                                       actual[1], expected_after, kValueRelTol));
            expect_near_relative(actual[1], expected_after, kValueRelTol);

            const amrex::Real centered_expected = reference_upwind_blended(stencil, adv_type, amrex::Real(1.0), amrex::Real(0.0));
            const amrex::Real zero_fraction_actual =
                evaluate_direct_upwind_policy_on_device(stencil, adv_type, amrex::Real(1.0), amrex::Real(0.0));
            SCOPED_TRACE(blended_trace(adv_type, amrex::Real(1.0), amrex::Real(0.0), sample_index,
                                       zero_fraction_actual, centered_expected, kValueRelTol));
            expect_near_relative(zero_fraction_actual, centered_expected, kValueRelTol);
        }
    }
}