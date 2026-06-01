#include <array>
#include <cmath>
#include <string>
#include <utility>

#include <AMReX_FArrayBox.H>
#include <AMReX_GpuContainers.H>

#include <gtest/gtest.h>

#include "ERF_GTestInterpolationCommon.H"

using namespace interpolation_test;

namespace {

inline amrex::Box centered_face_box ()
{
    return amrex::Box(amrex::IntVect(-4, 0, 0), amrex::IntVect(4, 0, 0));
}

void fill_centered_face_box (amrex::FArrayBox& qty,
                             const CenteredFaceStencil9& stencil)
{
    auto qty_arr = qty.array();
    amrex::single_task([=] AMREX_GPU_DEVICE () noexcept {
        qty_arr(-4, 0, 0, 0) = stencil.qm4;
        qty_arr(-3, 0, 0, 0) = stencil.qm3;
        qty_arr(-2, 0, 0, 0) = stencil.qm2;
        qty_arr(-1, 0, 0, 0) = stencil.qm1;
        qty_arr( 0, 0, 0, 0) = stencil.q0;
        qty_arr( 1, 0, 0, 0) = stencil.q1;
        qty_arr( 2, 0, 0, 0) = stencil.q2;
        qty_arr( 3, 0, 0, 0) = stencil.q3;
        qty_arr( 4, 0, 0, 0) = stencil.q4;
    });
    gpu_sync();
}

CenteredFaceStencil9 make_constant_stencil9 (const amrex::Real value)
{
    return CenteredFaceStencil9{value, value, value, value, value, value, value, value, value};
}

CenteredFaceStencil9 mirror_face_plus_half_stencil9 (const CenteredFaceStencil9& stencil)
{
    return CenteredFaceStencil9{
        amrex::Real(0.0),
        stencil.q4,
        stencil.q3,
        stencil.q2,
        stencil.q1,
        stencil.q0,
        stencil.qm1,
        stencil.qm2,
        stencil.qm3};
}

CenteredFaceStencil3 take_stencil3 (const CenteredFaceStencil9& stencil)
{
    return CenteredFaceStencil3{stencil.qm1, stencil.q0, stencil.q1};
}

CenteredFaceStencil5 take_stencil5 (const CenteredFaceStencil9& stencil)
{
    return CenteredFaceStencil5{stencil.qm2, stencil.qm1, stencil.q0, stencil.q1, stencil.q2};
}

CenteredFaceStencil7 take_stencil7 (const CenteredFaceStencil9& stencil)
{
    return CenteredFaceStencil7{stencil.qm3, stencil.qm2, stencil.qm1, stencil.q0,
                                stencil.q1, stencil.q2, stencil.q3};
}

CenteredFaceStencil3 take_right_stencil3 (const CenteredFaceStencil9& stencil)
{
    return CenteredFaceStencil3{stencil.q2, stencil.q1, stencil.q0};
}

CenteredFaceStencil5 take_right_stencil5 (const CenteredFaceStencil9& stencil)
{
    return CenteredFaceStencil5{stencil.q3, stencil.q2, stencil.q1, stencil.q0, stencil.qm1};
}

CenteredFaceStencil7 take_right_stencil7 (const CenteredFaceStencil9& stencil)
{
    return CenteredFaceStencil7{stencil.q4, stencil.q3, stencil.q2, stencil.q1,
                                stencil.q0, stencil.qm1, stencil.qm2};
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real evaluate_weno_x_device (const AdvType adv_type,
                                    const amrex::Array4<const amrex::Real>& qty,
                                    const int sigma) noexcept
{
    const amrex::Real upw = amrex::Real(sigma);
    amrex::Real value = amrex::Real(0.0);

    if (adv_type == AdvType::Weno_3) {
        WENO3 policy(qty, amrex::Real(1.0));
        policy.InterpolateInX(1, 0, 0, kQtyComp, value, upw);
    } else if (adv_type == AdvType::Weno_5) {
        WENO5 policy(qty, amrex::Real(1.0));
        policy.InterpolateInX(1, 0, 0, kQtyComp, value, upw);
    } else if (adv_type == AdvType::Weno_7) {
        WENO7 policy(qty, amrex::Real(1.0));
        policy.InterpolateInX(1, 0, 0, kQtyComp, value, upw);
    } else if (adv_type == AdvType::Weno_3Z) {
        WENO_Z3 policy(qty, amrex::Real(1.0));
        policy.InterpolateInX(1, 0, 0, kQtyComp, value, upw);
    } else if (adv_type == AdvType::Weno_3MZQ) {
        WENO_MZQ3 policy(qty, amrex::Real(1.0));
        policy.InterpolateInX(1, 0, 0, kQtyComp, value, upw);
    } else if (adv_type == AdvType::Weno_5Z) {
        WENO_Z5 policy(qty, amrex::Real(1.0));
        policy.InterpolateInX(1, 0, 0, kQtyComp, value, upw);
    } else {
        WENO_Z7 policy(qty, amrex::Real(1.0));
        policy.InterpolateInX(1, 0, 0, kQtyComp, value, upw);
    }

    return value;
}

amrex::Real evaluate_weno_x_on_device (const CenteredFaceStencil9& stencil,
                                       const AdvType adv_type,
                                       const int sigma)
{
    amrex::FArrayBox qty(centered_face_box(), 1);
    fill_centered_face_box(qty, stencil);
    const auto qty_arr = qty.const_array();
    amrex::Gpu::DeviceVector<amrex::Real> device_result(1, amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real> host_result(1, amrex::Real(0.0));
    amrex::Real* result_ptr = device_result.data();

    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int) noexcept {
        result_ptr[0] = evaluate_weno_x_device(adv_type, qty_arr, sigma);
    });
    gpu_sync();

    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_result.begin(), device_result.end(), host_result.begin());
    return host_result[0];
}

// MUST MATCH: these raw tau helpers intentionally encode the WENO-Z tau
// definitions so the symmetry tests can select cases that exercise both raw
// tau signs before the scheme takes abs(tau).
amrex::Real raw_tau3 (const CenteredFaceStencil3& stencil)
{
    return weno3_beta1(stencil) - weno3_beta0(stencil);
}

amrex::Real raw_tau5 (const CenteredFaceStencil5& stencil)
{
    return weno5_beta2(stencil) - weno5_beta0(stencil);
}

amrex::Real raw_tau7 (const CenteredFaceStencil7& stencil)
{
    return weno7_beta0(stencil) - weno7_beta1(stencil) -
           weno7_beta2(stencil) + weno7_beta3(stencil);
}

amrex::Real weno3_reference_value (const CenteredFaceStencil3& stencil)
{
    const amrex::Real w0 = weno_classic_weight(amrex::Real(1.0) / amrex::Real(3.0), weno3_beta0(stencil));
    const amrex::Real w1 = weno_classic_weight(amrex::Real(2.0) / amrex::Real(3.0), weno3_beta1(stencil));
    const amrex::Real wsum = w0 + w1;
    return (w0 * weno3_fv_candidate0(stencil) + w1 * weno3_fv_candidate1(stencil)) / wsum;
}

amrex::Real weno5_reference_value (const CenteredFaceStencil5& stencil)
{
    const amrex::Real w0 = weno_classic_weight(amrex::Real(1.0) / amrex::Real(10.0), weno5_beta0(stencil));
    const amrex::Real w1 = weno_classic_weight(amrex::Real(3.0) / amrex::Real(5.0), weno5_beta1(stencil));
    const amrex::Real w2 = weno_classic_weight(amrex::Real(3.0) / amrex::Real(10.0), weno5_beta2(stencil));
    const amrex::Real wsum = w0 + w1 + w2;
    return (w0 * weno5_fv_candidate0(stencil) +
            w1 * weno5_fv_candidate1(stencil) +
            w2 * weno5_fv_candidate2(stencil)) / wsum;
}

amrex::Real weno7_reference_value (const CenteredFaceStencil7& stencil)
{
    const amrex::Real w0 = weno_classic_weight(amrex::Real(1.0) / amrex::Real(35.0), weno7_beta0(stencil));
    const amrex::Real w1 = weno_classic_weight(amrex::Real(12.0) / amrex::Real(35.0), weno7_beta1(stencil));
    const amrex::Real w2 = weno_classic_weight(amrex::Real(18.0) / amrex::Real(35.0), weno7_beta2(stencil));
    const amrex::Real w3 = weno_classic_weight(amrex::Real(4.0) / amrex::Real(35.0), weno7_beta3(stencil));
    const amrex::Real wsum = w0 + w1 + w2 + w3;
    return (w0 * weno7_fv_candidate0(stencil) +
            w1 * weno7_fv_candidate1(stencil) +
            w2 * weno7_fv_candidate2(stencil) +
            w3 * weno7_fv_candidate3(stencil)) / wsum;
}

amrex::Real weno3z_reference_value (const CenteredFaceStencil3& stencil)
{
    const amrex::Real tau = weno3_tau(stencil);
    const amrex::Real w0 = weno_z_weight(amrex::Real(1.0) / amrex::Real(3.0), weno3_beta0(stencil), tau);
    const amrex::Real w1 = weno_z_weight(amrex::Real(2.0) / amrex::Real(3.0), weno3_beta1(stencil), tau);
    const amrex::Real wsum = w0 + w1;
    return (w0 * weno3_fv_candidate0(stencil) + w1 * weno3_fv_candidate1(stencil)) / wsum;
}

amrex::Real weno_mzq3_reference_value (const CenteredFaceStencil3& stencil)
{
    const amrex::Real beta0 = weno3_beta0(stencil);
    const amrex::Real beta1 = weno3_beta1(stencil);
    const amrex::Real beta_c = weno_mzq3_beta_c(stencil);
    // The MZQ tau scaling is method-specific; tests document and protect the
    // existing contract unless a derivation shows otherwise.
    const amrex::Real tau = (amrex::Math::abs(beta_c - beta0) +
                             amrex::Math::abs(beta_c - beta1)) / amrex::Real(32.0);
    const amrex::Real a0 = weno_z_weight(amrex::Real(1.0) / amrex::Real(3.0), beta0, tau);
    const amrex::Real a1 = weno_z_weight(amrex::Real(1.0) / amrex::Real(3.0), beta1, tau);
    const amrex::Real a2 = weno_z_weight(amrex::Real(1.0) / amrex::Real(3.0), beta_c, tau);
    const amrex::Real asum = a0 + a1 + a2;
    const amrex::Real w0 = a0 / asum;
    const amrex::Real w1 = a1 / asum;
    const amrex::Real w2 = a2 / asum;
    const amrex::Real v0 = weno3_fv_candidate0(stencil);
    const amrex::Real v1 = weno3_fv_candidate1(stencil);
    const amrex::Real v2 = (-stencil.qm1 + amrex::Real(5.0) * stencil.q0 + amrex::Real(2.0) * stencil.q1) /
                           amrex::Real(6.0);
    return (w2 / (amrex::Real(1.0) / amrex::Real(3.0))) *
               (v2 - (amrex::Real(1.0) / amrex::Real(3.0)) * v0 - (amrex::Real(1.0) / amrex::Real(3.0)) * v1) +
           w0 * v0 + w1 * v1;
}

amrex::Real weno5z_reference_value (const CenteredFaceStencil5& stencil)
{
    const amrex::Real tau = weno5_tau(stencil);
    const amrex::Real w0 = weno_z_weight(amrex::Real(1.0) / amrex::Real(10.0), weno5_beta0(stencil), tau);
    const amrex::Real w1 = weno_z_weight(amrex::Real(3.0) / amrex::Real(5.0), weno5_beta1(stencil), tau);
    const amrex::Real w2 = weno_z_weight(amrex::Real(3.0) / amrex::Real(10.0), weno5_beta2(stencil), tau);
    const amrex::Real wsum = w0 + w1 + w2;
    return (w0 * weno5_fv_candidate0(stencil) +
            w1 * weno5_fv_candidate1(stencil) +
            w2 * weno5_fv_candidate2(stencil)) / wsum;
}

amrex::Real weno7z_reference_value (const CenteredFaceStencil7& stencil)
{
    const amrex::Real tau = weno7_tau(stencil);
    const amrex::Real w0 = weno_z_weight(amrex::Real(1.0) / amrex::Real(35.0), weno7_beta0(stencil), tau);
    const amrex::Real w1 = weno_z_weight(amrex::Real(12.0) / amrex::Real(35.0), weno7_beta1(stencil), tau);
    const amrex::Real w2 = weno_z_weight(amrex::Real(18.0) / amrex::Real(35.0), weno7_beta2(stencil), tau);
    const amrex::Real w3 = weno_z_weight(amrex::Real(4.0) / amrex::Real(35.0), weno7_beta3(stencil), tau);
    const amrex::Real wsum = w0 + w1 + w2 + w3;
    return (w0 * weno7_fv_candidate0(stencil) +
            w1 * weno7_fv_candidate1(stencil) +
            w2 * weno7_fv_candidate2(stencil) +
            w3 * weno7_fv_candidate3(stencil)) / wsum;
}

amrex::Real reference_weno_face_value (const CenteredFaceStencil9& stencil,
                                       const AdvType adv_type,
                                       const int sigma)
{
    if (adv_type == AdvType::Weno_3) {
        return weno3_reference_value((sigma >= 0) ? take_stencil3(stencil) : take_right_stencil3(stencil));
    }
    if (adv_type == AdvType::Weno_5) {
        return weno5_reference_value((sigma >= 0) ? take_stencil5(stencil) : take_right_stencil5(stencil));
    }
    if (adv_type == AdvType::Weno_7) {
        return weno7_reference_value((sigma >= 0) ? take_stencil7(stencil) : take_right_stencil7(stencil));
    }
    if (adv_type == AdvType::Weno_3Z) {
        return weno3z_reference_value((sigma >= 0) ? take_stencil3(stencil) : take_right_stencil3(stencil));
    }
    if (adv_type == AdvType::Weno_3MZQ) {
        return weno_mzq3_reference_value((sigma >= 0) ? take_stencil3(stencil) : take_right_stencil3(stencil));
    }
    if (adv_type == AdvType::Weno_5Z) {
        return weno5z_reference_value((sigma >= 0) ? take_stencil5(stencil) : take_right_stencil5(stencil));
    }
    return weno7z_reference_value((sigma >= 0) ? take_stencil7(stencil) : take_right_stencil7(stencil));
}

std::string weno_trace (const AdvType adv_type,
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

// Motivation: The nonlinear WENO path must preserve constants and stay finite
// even when every beta collapses to zero.
TEST(InterpolationWENOScalar, ClassicWenoPreservesConstantsAndReturnsFiniteValues)
{
    const std::array<AdvType, 3> schemes = {{AdvType::Weno_3, AdvType::Weno_5, AdvType::Weno_7}};
    const CenteredFaceStencil9 stencil = make_constant_stencil9(amrex::Real(1.75));

    for (const AdvType adv_type : schemes) {
        for (const int sigma : std::array<int, 2>{{-1, 1}}) {
            const amrex::Real actual = evaluate_weno_x_on_device(stencil, adv_type, sigma);
            EXPECT_TRUE(std::isfinite(actual));
            EXPECT_NEAR(actual, stencil.q0, scaled_tol(actual, stencil.q0, kWenoIdentityRelTol))
                << weno_trace(adv_type, sigma, 0, actual, stencil.q0, kWenoIdentityRelTol);
        }
    }
}

// Motivation: A linear finite-volume field is the minimum smooth contract for
// every WENO family, regardless of the nonlinear weights.
TEST(InterpolationWENOScalar, ClassicWenoIsExactForLinearFiniteVolumeData)
{
    const std::array<AdvType, 3> schemes = {{AdvType::Weno_3, AdvType::Weno_5, AdvType::Weno_7}};
    const CenteredFaceStencil9 stencil = make_centered_monomial_stencil9(1);
    const amrex::Real expected = centered_face_value_monomial(1, kFacePlusHalf);

    for (const AdvType adv_type : schemes) {
        for (const int sigma : std::array<int, 2>{{-1, 1}}) {
            const amrex::Real actual = evaluate_weno_x_on_device(stencil, adv_type, sigma);
            EXPECT_NEAR(actual, expected, scaled_tol(actual, expected, kWenoValueRelTol))
                << weno_trace(adv_type, sigma, 1, actual, expected, kWenoValueRelTol);
        }
    }
}

// Motivation: The optimal linear combinations document the intended smooth
// finite-volume design order for the classic WENO families.
TEST(InterpolationWENOScalar, ClassicWenoOptimalFiniteVolumeReferencesMatchDesignDegree)
{
    for (int degree = 0; degree <= 2; ++degree) {
        const CenteredFaceStencil3 stencil = make_centered_monomial_stencil3(degree);
        const amrex::Real actual = weno3_fv_optimal_face_value(stencil);
        const amrex::Real expected = centered_face_value_monomial(degree, kFacePlusHalf);
        EXPECT_NEAR(actual, expected, scaled_tol(actual, expected, kPolynomialRelTol));
    }

    for (int degree = 0; degree <= 4; ++degree) {
        const CenteredFaceStencil5 stencil = make_centered_monomial_stencil5(degree);
        const amrex::Real actual = weno5_fv_optimal_face_value(stencil);
        const amrex::Real expected = centered_face_value_monomial(degree, kFacePlusHalf);
        EXPECT_NEAR(actual, expected, scaled_tol(actual, expected, kPolynomialRelTol));
    }

    for (int degree = 0; degree <= 6; ++degree) {
        const CenteredFaceStencil7 stencil = make_centered_monomial_stencil7(degree);
        const amrex::Real actual = weno7_fv_optimal_face_value(stencil);
        const amrex::Real expected = centered_face_value_monomial(degree, kFacePlusHalf);
        EXPECT_NEAR(actual, expected, scaled_tol(actual, expected, kPolynomialRelTol));
    }
}

// Motivation: Each candidate polynomial should exactly reconstruct the face
// value for the degree supported by its substencil.
TEST(InterpolationWENOScalar, CandidateFiniteVolumeReferencesMatchSubstencilDegree)
{
    for (int degree = 0; degree <= 1; ++degree) {
        const CenteredFaceStencil3 stencil = make_centered_monomial_stencil3(degree);
        const amrex::Real expected = centered_face_value_monomial(degree, kFacePlusHalf);
        EXPECT_NEAR(weno3_fv_candidate0(stencil), expected,
                    scaled_tol(weno3_fv_candidate0(stencil), expected, kPolynomialRelTol));
        EXPECT_NEAR(weno3_fv_candidate1(stencil), expected,
                    scaled_tol(weno3_fv_candidate1(stencil), expected, kPolynomialRelTol));
    }

    for (int degree = 0; degree <= 2; ++degree) {
        const CenteredFaceStencil5 stencil = make_centered_monomial_stencil5(degree);
        const amrex::Real expected = centered_face_value_monomial(degree, kFacePlusHalf);
        EXPECT_NEAR(weno5_fv_candidate0(stencil), expected,
                    scaled_tol(weno5_fv_candidate0(stencil), expected, kPolynomialRelTol));
        EXPECT_NEAR(weno5_fv_candidate1(stencil), expected,
                    scaled_tol(weno5_fv_candidate1(stencil), expected, kPolynomialRelTol));
        EXPECT_NEAR(weno5_fv_candidate2(stencil), expected,
                    scaled_tol(weno5_fv_candidate2(stencil), expected, kPolynomialRelTol));
    }

    for (int degree = 0; degree <= 3; ++degree) {
        const CenteredFaceStencil7 stencil = make_centered_monomial_stencil7(degree);
        const amrex::Real expected = centered_face_value_monomial(degree, kFacePlusHalf);
        EXPECT_NEAR(weno7_fv_candidate0(stencil), expected,
                    scaled_tol(weno7_fv_candidate0(stencil), expected, kPolynomialRelTol));
        EXPECT_NEAR(weno7_fv_candidate1(stencil), expected,
                    scaled_tol(weno7_fv_candidate1(stencil), expected, kPolynomialRelTol));
        EXPECT_NEAR(weno7_fv_candidate2(stencil), expected,
                    scaled_tol(weno7_fv_candidate2(stencil), expected, kPolynomialRelTol));
        EXPECT_NEAR(weno7_fv_candidate3(stencil), expected,
                    scaled_tol(weno7_fv_candidate3(stencil), expected, kPolynomialRelTol));
    }
}

// Motivation: This is the direct finite-volume conformance check: the source
// implementation must match an independent nonlinear reference, not just the
// smooth-limit optimal stencil.
TEST(InterpolationWENOScalar, NonlinearSchemesMatchIndependentFiniteVolumeReferences)
{
    const std::array<AdvType, 7> schemes = {{AdvType::Weno_3, AdvType::Weno_5, AdvType::Weno_7,
                                             AdvType::Weno_3Z, AdvType::Weno_3MZQ,
                                             AdvType::Weno_5Z, AdvType::Weno_7Z}};
    const std::array<CenteredFaceStencil9, 3> samples = {{
        CenteredFaceStencil9{amrex::Real(0.8), amrex::Real(-1.2), amrex::Real(0.45), amrex::Real(-0.6),
                            amrex::Real(1.1), amrex::Real(-0.35), amrex::Real(0.9), amrex::Real(-1.4), amrex::Real(0.5)},
        CenteredFaceStencil9{amrex::Real(-0.25), amrex::Real(0.4), amrex::Real(-0.9), amrex::Real(1.35),
                            amrex::Real(-0.2), amrex::Real(0.75), amrex::Real(-1.1), amrex::Real(0.6), amrex::Real(-0.3)},
        CenteredFaceStencil9{amrex::Real(1.0), amrex::Real(0.6), amrex::Real(-0.4), amrex::Real(0.2),
                            amrex::Real(-0.1), amrex::Real(0.3), amrex::Real(-0.7), amrex::Real(1.4), amrex::Real(-0.8)} }};

    for (const AdvType adv_type : schemes) {
        for (int sample_index = 0; sample_index < static_cast<int>(samples.size()); ++sample_index) {
            for (const int sigma : std::array<int, 2>{{-1, 1}}) {
                const amrex::Real actual = evaluate_weno_x_on_device(samples[sample_index], adv_type, sigma);
                const amrex::Real expected = reference_weno_face_value(samples[sample_index], adv_type, sigma);
                EXPECT_NEAR(actual, expected, scaled_tol(actual, expected, kWenoValueRelTol))
                    << weno_trace(adv_type, sigma, sample_index, actual, expected, kWenoValueRelTol);
            }
        }
    }
}

// Motivation: Mirroring the stencil and flipping the upwind sign must produce
// the same face state for each nonlinear WENO variant.
TEST(InterpolationWENOScalar, NonlinearSchemesRespectMirrorSymmetry)
{
    const std::array<AdvType, 7> schemes = {{AdvType::Weno_3, AdvType::Weno_5, AdvType::Weno_7,
                                             AdvType::Weno_3Z, AdvType::Weno_3MZQ,
                                             AdvType::Weno_5Z, AdvType::Weno_7Z}};
    const std::array<CenteredFaceStencil9, 2> samples = {{
        CenteredFaceStencil9{amrex::Real(-0.4), amrex::Real(0.3), amrex::Real(1.2), amrex::Real(-0.8),
                            amrex::Real(0.5), amrex::Real(-1.1), amrex::Real(0.9), amrex::Real(0.2), amrex::Real(-0.6)},
        CenteredFaceStencil9{amrex::Real(0.15), amrex::Real(-0.9), amrex::Real(0.7), amrex::Real(0.1),
                            amrex::Real(-0.2), amrex::Real(1.0), amrex::Real(-0.5), amrex::Real(0.8), amrex::Real(-1.3)} }};

    for (const AdvType adv_type : schemes) {
        for (int sample_index = 0; sample_index < static_cast<int>(samples.size()); ++sample_index) {
            const CenteredFaceStencil9 mirrored = mirror_face_plus_half_stencil9(samples[sample_index]);
            const amrex::Real positive = evaluate_weno_x_on_device(samples[sample_index], adv_type, 1);
            const amrex::Real negative = evaluate_weno_x_on_device(mirrored, adv_type, -1);
            EXPECT_NEAR(positive, negative, scaled_tol(positive, negative, kWenoValueRelTol))
                << weno_trace(adv_type, 1, sample_index, positive, negative, kWenoValueRelTol);
        }
    }
}

// Motivation: The WENO-Z families should inherit the same constant and linear
// finite-volume invariants as the classic WENO schemes.
TEST(InterpolationWENOScalar, WenoZVariantsPreserveConstantsAndLinearFields)
{
    const std::array<AdvType, 4> schemes = {{AdvType::Weno_3Z, AdvType::Weno_3MZQ,
                                             AdvType::Weno_5Z, AdvType::Weno_7Z}};
    const CenteredFaceStencil9 constant_stencil = make_constant_stencil9(amrex::Real(-0.875));
    const CenteredFaceStencil9 linear_stencil = make_centered_monomial_stencil9(1);
    const amrex::Real linear_expected = centered_face_value_monomial(1, kFacePlusHalf);

    for (const AdvType adv_type : schemes) {
        for (const int sigma : std::array<int, 2>{{-1, 1}}) {
            const amrex::Real constant_actual = evaluate_weno_x_on_device(constant_stencil, adv_type, sigma);
            EXPECT_TRUE(std::isfinite(constant_actual));
            EXPECT_NEAR(constant_actual, constant_stencil.q0,
                        scaled_tol(constant_actual, constant_stencil.q0, kWenoIdentityRelTol))
                << weno_trace(adv_type, sigma, 0, constant_actual, constant_stencil.q0, kWenoIdentityRelTol);

            const amrex::Real linear_actual = evaluate_weno_x_on_device(linear_stencil, adv_type, sigma);
            EXPECT_NEAR(linear_actual, linear_expected,
                        scaled_tol(linear_actual, linear_expected, kWenoValueRelTol))
                << weno_trace(adv_type, sigma, 1, linear_actual, linear_expected, kWenoValueRelTol);
        }
    }
}

// Motivation: The WENO-Z smooth limit must recover the independently derived
// finite-volume optimal face value through the advertised design degree.
TEST(InterpolationWENOScalar, WenoZOptimalFiniteVolumeReferencesMatchSmoothLimit)
{
    for (int degree = 0; degree <= 2; ++degree) {
        const CenteredFaceStencil3 stencil = make_centered_monomial_stencil3(degree);
        const amrex::Real actual = weno3_fv_optimal_face_value(stencil);
        const amrex::Real expected = centered_face_value_monomial(degree, kFacePlusHalf);
        EXPECT_NEAR(actual, expected, scaled_tol(actual, expected, kPolynomialRelTol));
    }

    for (int degree = 0; degree <= 4; ++degree) {
        const CenteredFaceStencil5 stencil = make_centered_monomial_stencil5(degree);
        const amrex::Real actual = weno5_fv_optimal_face_value(stencil);
        const amrex::Real expected = centered_face_value_monomial(degree, kFacePlusHalf);
        EXPECT_NEAR(actual, expected, scaled_tol(actual, expected, kPolynomialRelTol));
    }

    for (int degree = 0; degree <= 6; ++degree) {
        const CenteredFaceStencil7 stencil = make_centered_monomial_stencil7(degree);
        const amrex::Real actual = weno7_fv_optimal_face_value(stencil);
        const amrex::Real expected = centered_face_value_monomial(degree, kFacePlusHalf);
        EXPECT_NEAR(actual, expected, scaled_tol(actual, expected, kPolynomialRelTol));
    }
}

// Motivation: These cases force both raw tau signs so the tests protect the
// documented absolute-value tau contract used by the WENO-Z schemes.
TEST(InterpolationWENOScalar, WenoZTauAbsoluteValueMirrorSymmetry)
{
    const CenteredFaceStencil9 z3_positive{amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0),
                                           amrex::Real(0.0), amrex::Real(-0.5), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0)};
    const CenteredFaceStencil9 z3_negative{amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(-0.5),
                                           amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0)};
    const CenteredFaceStencil9 z5_positive{amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0),
                                           amrex::Real(0.0), amrex::Real(-0.5), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0)};
    const CenteredFaceStencil9 z5_negative{amrex::Real(0.0), amrex::Real(0.0), amrex::Real(-0.5), amrex::Real(0.0),
                                           amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0)};
    const CenteredFaceStencil9 z7_positive{amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0),
                                           amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(-0.5), amrex::Real(0.0)};
    const CenteredFaceStencil9 z7_negative{amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0),
                                           amrex::Real(0.0), amrex::Real(-0.5), amrex::Real(-0.5), amrex::Real(0.0), amrex::Real(0.0)};

    EXPECT_GT(raw_tau3(take_stencil3(z3_positive)), amrex::Real(0.0));
    EXPECT_LT(raw_tau3(take_stencil3(z3_negative)), amrex::Real(0.0));
    EXPECT_GT(raw_tau5(take_stencil5(z5_positive)), amrex::Real(0.0));
    EXPECT_LT(raw_tau5(take_stencil5(z5_negative)), amrex::Real(0.0));
    EXPECT_GT(raw_tau7(take_stencil7(z7_positive)), amrex::Real(0.0));
    EXPECT_LT(raw_tau7(take_stencil7(z7_negative)), amrex::Real(0.0));

    const std::array<std::pair<AdvType, CenteredFaceStencil9>, 6> cases = {{
        {AdvType::Weno_3Z, z3_positive},
        {AdvType::Weno_3Z, z3_negative},
        {AdvType::Weno_5Z, z5_positive},
        {AdvType::Weno_5Z, z5_negative},
        {AdvType::Weno_7Z, z7_positive},
        {AdvType::Weno_7Z, z7_negative}}};

    for (int sample_index = 0; sample_index < static_cast<int>(cases.size()); ++sample_index) {
        const AdvType adv_type = cases[sample_index].first;
        const CenteredFaceStencil9& stencil = cases[sample_index].second;
        const CenteredFaceStencil9 mirrored = mirror_face_plus_half_stencil9(stencil);
        const amrex::Real positive = evaluate_weno_x_on_device(stencil, adv_type, 1);
        const amrex::Real negative_mirror = evaluate_weno_x_on_device(mirrored, adv_type, -1);

        EXPECT_NEAR(positive, negative_mirror,
                    scaled_tol(positive, negative_mirror, kWenoValueRelTol))
            << weno_trace(adv_type, 1, sample_index, positive, negative_mirror, kWenoValueRelTol);
    }
}
