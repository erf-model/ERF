#ifndef ERF_GTEST_INTERPOLATION_COMMON_H_
#define ERF_GTEST_INTERPOLATION_COMMON_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <AMReX_IntVect.H>
#include <AMReX_Math.H>

#include <ERF_DataStruct.H>
#include <ERF_IndexDefines.H>
#include <ERF_Interpolation.H>

namespace interpolation_test {

// Interpolation unit-test contract
// --------------------------------
// These tests check finite-volume reconstruction contracts for ERF's
// non-WENO interpolation helpers, not implementation accidents.
//
// Core invariants:
//   1. Constant fields are preserved.
//   2. Centered/upwind schemes satisfy their finite-volume polynomial
//      exactness order.
//   3. Zero upwind reduces odd upwind corrections to the next centered
//      stencil.
//   4. Positive and negative upwind signs are mirror-symmetric.
//   5. InterpolatePertFromCell reconstructs qty - r0_arr.
//   6. Kernel tests exercise AMReX device callability over branch-dense
//      sample grids.
//
// Portability:
//   - Device kernels compute values or normalized errors only.
//   - GTest EXPECT/ASSERT calls stay on the host.
//   - No dynamic allocation is allowed inside device code.
//   - Device coverage uses representative case tables, not single-sample
//     smoke tests.

// Relative tolerance for direct value comparisons between independent host
// references and production device-path reconstructions.
// Provenance: machine epsilon, float-vs-double scaling, and a small safety
// factor for a handful of adds/multiplies with asymmetric stencils.
#ifdef AMREX_USE_FLOAT
constexpr amrex::Real kValueRelTol = amrex::Real(3.0e-5);
constexpr amrex::Real kIdentityRelTol = amrex::Real(2.0e-5);
constexpr amrex::Real kPolynomialRelTol = amrex::Real(6.0e-5);
constexpr amrex::Real kDeviceRelTol = amrex::Real(8.0e-5);
#else
constexpr amrex::Real kValueRelTol = amrex::Real(1.0e-12);
constexpr amrex::Real kIdentityRelTol = amrex::Real(1.0e-12);
constexpr amrex::Real kPolynomialRelTol = amrex::Real(2.0e-12);
constexpr amrex::Real kDeviceRelTol = amrex::Real(5.0e-12);
#endif

constexpr int kQtyComp = 0;
constexpr int kEvalLo = -1;
constexpr int kEvalHi = 1;
constexpr amrex::Real kPerturbationConstant = amrex::Real(10.0);

constexpr std::array<AdvType, 5> kSupportedAdvTypes = {{
    AdvType::Centered_2nd,
    AdvType::Upwind_3rd,
    AdvType::Centered_4th,
    AdvType::Upwind_5th,
    AdvType::Centered_6th}};

constexpr std::array<int, 3> kRepresentativeUpwindSigns = {{-1, 0, 1}};

constexpr std::array<AdvType, 9> kUnsupportedAdvTypes = {{
    AdvType::Upwind_3rd_SL,
    AdvType::Weno_3,
    AdvType::Weno_3Z,
    AdvType::Weno_5,
    AdvType::Weno_5Z,
    AdvType::Weno_3MZQ,
    AdvType::Weno_7,
    AdvType::Weno_7Z,
    AdvType::Unknown}};

struct StencilValues {
    amrex::Real qm3;
    amrex::Real qm2;
    amrex::Real qm1;
    amrex::Real q0;
    amrex::Real q1;
    amrex::Real q2;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int signum (const amrex::Real value) noexcept
{
    return (value > amrex::Real(0)) - (value < amrex::Real(0));
}

inline std::string scheme_name (const AdvType adv_type)
{
    switch (adv_type) {
    case AdvType::Centered_2nd: return "Centered_2nd";
    case AdvType::Upwind_3rd:   return "Upwind_3rd";
    case AdvType::Centered_4th: return "Centered_4th";
    case AdvType::Upwind_5th:   return "Upwind_5th";
    case AdvType::Centered_6th: return "Centered_6th";
    case AdvType::Upwind_3rd_SL:return "Upwind_3rd_SL";
    case AdvType::Weno_3:       return "Weno_3";
    case AdvType::Weno_3Z:      return "Weno_3Z";
    case AdvType::Weno_5:       return "Weno_5";
    case AdvType::Weno_5Z:      return "Weno_5Z";
    case AdvType::Weno_3MZQ:    return "Weno_3MZQ";
    case AdvType::Weno_7:       return "Weno_7";
    case AdvType::Weno_7Z:      return "Weno_7Z";
    default:                    return "Unknown";
    }
}

inline std::string direction_name (const Coord direction)
{
    switch (direction) {
    case Coord::x: return "x";
    case Coord::y: return "y";
    default:       return "z";
    }
}

inline amrex::Real scaled_tol (const amrex::Real actual,
                               const amrex::Real expected,
                               const amrex::Real factor)
{
    return factor * std::max({amrex::Real(1.0), std::abs(actual), std::abs(expected)});
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real portable_scaled_tol (const amrex::Real actual,
                                 const amrex::Real expected,
                                 const amrex::Real factor) noexcept
{
    return factor * amrex::max(amrex::Real(1.0),
                               amrex::max(amrex::Math::abs(actual), amrex::Math::abs(expected)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real normalized_error (const amrex::Real actual,
                              const amrex::Real expected,
                              const amrex::Real factor) noexcept
{
    return amrex::Math::abs(actual - expected) / portable_scaled_tol(actual, expected, factor);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real integer_power (const amrex::Real base, const int exponent) noexcept
{
    amrex::Real value = amrex::Real(1.0);
    for (int n = 0; n < exponent; ++n) {
        value *= base;
    }
    return value;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real finite_volume_cell_average_monomial (const int degree,
                                                 const int cell_lo) noexcept
{
    if (degree == 0) {
        return amrex::Real(1.0);
    }

    const int antiderivative_degree = degree + 1;
    return (integer_power(amrex::Real(cell_lo + 1), antiderivative_degree) -
            integer_power(amrex::Real(cell_lo), antiderivative_degree)) /
           amrex::Real(antiderivative_degree);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
StencilValues make_monomial_stencil (const int degree) noexcept
{
    return StencilValues{
        finite_volume_cell_average_monomial(degree, -3),
        finite_volume_cell_average_monomial(degree, -2),
        finite_volume_cell_average_monomial(degree, -1),
        finite_volume_cell_average_monomial(degree, 0),
        finite_volume_cell_average_monomial(degree, 1),
        finite_volume_cell_average_monomial(degree, 2)};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
StencilValues mirrored_stencil (const StencilValues& stencil) noexcept
{
    return StencilValues{stencil.q2, stencil.q1, stencil.q0,
                         stencil.qm1, stencil.qm2, stencil.qm3};
}

// Reference coefficients are derived from finite-volume polynomial
// moment conditions about the face x = 0. They intentionally avoid
// calling interpolatedVal or using production avg/diff assembly.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real reference_reconstruction (const StencilValues& stencil,
                                      const AdvType adv_type,
                                      const int sigma) noexcept
{
    if (adv_type == AdvType::Centered_2nd) {
        return myhalf * (stencil.qm1 + stencil.q0);
    }

    if (adv_type == AdvType::Upwind_3rd) {
        return ((amrex::Real(-1.0) + amrex::Real(sigma)) / amrex::Real(12.0)) * stencil.q1 +
               ((amrex::Real(7.0)  - amrex::Real(3.0 * sigma)) / amrex::Real(12.0)) * stencil.q0 +
               ((amrex::Real(7.0)  + amrex::Real(3.0 * sigma)) / amrex::Real(12.0)) * stencil.qm1 +
               ((amrex::Real(-1.0) - amrex::Real(sigma)) / amrex::Real(12.0)) * stencil.qm2;
    }

    if (adv_type == AdvType::Centered_4th) {
        return (-amrex::Real(1.0) / amrex::Real(12.0)) * stencil.q1 +
               ( amrex::Real(7.0) / amrex::Real(12.0)) * stencil.q0 +
               ( amrex::Real(7.0) / amrex::Real(12.0)) * stencil.qm1 +
               (-amrex::Real(1.0) / amrex::Real(12.0)) * stencil.qm2;
    }

    if (adv_type == AdvType::Upwind_5th) {
        return ((amrex::Real(1.0)  - amrex::Real(sigma)) / amrex::Real(60.0)) * stencil.q2 +
               ((-amrex::Real(8.0) + amrex::Real(5.0 * sigma)) / amrex::Real(60.0)) * stencil.q1 +
               (( amrex::Real(37.0) - amrex::Real(10.0 * sigma)) / amrex::Real(60.0)) * stencil.q0 +
               (( amrex::Real(37.0) + amrex::Real(10.0 * sigma)) / amrex::Real(60.0)) * stencil.qm1 +
               ((-amrex::Real(8.0) - amrex::Real(5.0 * sigma)) / amrex::Real(60.0)) * stencil.qm2 +
               (( amrex::Real(1.0) + amrex::Real(sigma)) / amrex::Real(60.0)) * stencil.qm3;
    }

    return ( amrex::Real(1.0) / amrex::Real(60.0)) * stencil.q2 +
           (-amrex::Real(8.0) / amrex::Real(60.0)) * stencil.q1 +
           ( amrex::Real(37.0) / amrex::Real(60.0)) * stencil.q0 +
           ( amrex::Real(37.0) / amrex::Real(60.0)) * stencil.qm1 +
           (-amrex::Real(8.0) / amrex::Real(60.0)) * stencil.qm2 +
           ( amrex::Real(1.0) / amrex::Real(60.0)) * stencil.qm3;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int exactness_degree (const AdvType adv_type) noexcept
{
    if (adv_type == AdvType::Centered_2nd) return 1;
    if (adv_type == AdvType::Upwind_3rd)   return 2;
    if (adv_type == AdvType::Centered_4th) return 3;
    if (adv_type == AdvType::Upwind_5th)   return 4;
    return 5;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int leading_error_degree (const AdvType adv_type) noexcept
{
    return exactness_degree(adv_type) + 1;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real leading_error_expected (const AdvType adv_type,
                                    const int sigma) noexcept
{
    if (adv_type == AdvType::Centered_2nd) return amrex::Real(1.0) / amrex::Real(3.0);
    if (adv_type == AdvType::Upwind_3rd)   return amrex::Real(sigma) / amrex::Real(2.0);
    if (adv_type == AdvType::Centered_4th) return -amrex::Real(4.0) / amrex::Real(5.0);
    if (adv_type == AdvType::Upwind_5th)   return -amrex::Real(2.0) * amrex::Real(sigma);
    return amrex::Real(36.0) / amrex::Real(7.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
AdvType centered_reduction (const AdvType adv_type) noexcept
{
    return (adv_type == AdvType::Upwind_3rd) ? AdvType::Centered_4th : AdvType::Centered_6th;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
StencilValues load_direction_stencil (const amrex::Array4<const amrex::Real>& qty,
                                      const int i, const int j, const int k,
                                      const Coord direction,
                                      const int comp = kQtyComp) noexcept
{
    if (direction == Coord::x) {
        return StencilValues{qty(i-3,j,k,comp), qty(i-2,j,k,comp), qty(i-1,j,k,comp),
                             qty(i  ,j,k,comp), qty(i+1,j,k,comp), qty(i+2,j,k,comp)};
    }

    if (direction == Coord::y) {
        return StencilValues{qty(i,j-3,k,comp), qty(i,j-2,k,comp), qty(i,j-1,k,comp),
                             qty(i,j  ,k,comp), qty(i,j+1,k,comp), qty(i,j+2,k,comp)};
    }

    return StencilValues{qty(i,j,k-3,comp), qty(i,j,k-2,comp), qty(i,j,k-1,comp),
                         qty(i,j,k  ,comp), qty(i,j,k+1,comp), qty(i,j,k+2,comp)};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void fill_avg_diff (const StencilValues& stencil,
                    amrex::Real& avg1,  amrex::Real& avg2,  amrex::Real& avg3,
                    amrex::Real& diff1, amrex::Real& diff2, amrex::Real& diff3) noexcept
{
    avg1  = stencil.q0 + stencil.qm1;
    diff1 = stencil.q0 - stencil.qm1;
    avg2  = stencil.q1 + stencil.qm2;
    diff2 = stencil.q1 - stencil.qm2;
    avg3  = stencil.q2 + stencil.qm3;
    diff3 = stencil.q2 - stencil.qm3;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real directional_field_value (const int i,
                                     const int j,
                                     const int k) noexcept
{
    return amrex::Real(0.5)   * finite_volume_cell_average_monomial(5, i) -
           amrex::Real(1.25)  * finite_volume_cell_average_monomial(4, i) +
           amrex::Real(0.75)  * finite_volume_cell_average_monomial(3, j) -
           amrex::Real(0.5)   * finite_volume_cell_average_monomial(2, j) +
           amrex::Real(1.5)   * finite_volume_cell_average_monomial(4, k) +
           amrex::Real(0.25)  * amrex::Real(i * j) -
           amrex::Real(0.125) * amrex::Real(j * k) +
           amrex::Real(0.0625)* amrex::Real(i * k) +
           amrex::Real(3.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real perturbation_background (const Coord direction,
                                     const int i,
                                     const int j,
                                     const int k) noexcept
{
    if (direction == Coord::x) {
        const amrex::Real x = amrex::Real(i);
        return integer_power(x, 3) + amrex::Real(0.25) * integer_power(x, 5);
    }

    if (direction == Coord::y) {
        const amrex::Real y = amrex::Real(j);
        return integer_power(y, 3) + amrex::Real(0.25) * integer_power(y, 5);
    }

    const amrex::Real z = amrex::Real(k);
    return integer_power(z, 3) + amrex::Real(0.25) * integer_power(z, 5);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real density_quantity_value (const int i,
                                    const int j,
                                    const int k) noexcept
{
    return amrex::Real(2.0) +
           amrex::Real(0.75) * finite_volume_cell_average_monomial(4, i) -
           amrex::Real(0.5)  * finite_volume_cell_average_monomial(3, j) +
           amrex::Real(0.25) * finite_volume_cell_average_monomial(2, k) +
           amrex::Real(0.125) * amrex::Real(i - 2 * j + 3 * k);
}

inline amrex::Box stencil_box ()
{
    return amrex::Box(amrex::IntVect(-4, -4, -4), amrex::IntVect(4, 4, 4));
}

inline amrex::Box eval_box ()
{
    return amrex::Box(amrex::IntVect(kEvalLo, kEvalLo, kEvalLo),
                      amrex::IntVect(kEvalHi, kEvalHi, kEvalHi));
}

inline void gpu_sync ()
{
    amrex::Gpu::streamSynchronize();
}

} // namespace interpolation_test

#endif