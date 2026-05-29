#include <array>
#include <cmath>

#include <AMReX_FArrayBox.H>

#include <gtest/gtest.h>

#include "ERF_GTestInterpolationCommon.H"

using namespace interpolation_test;

namespace {

#ifdef AMREX_USE_FLOAT
constexpr amrex::Real kObservedOrderErrorFloor = amrex::Real(1.0e-6);
#else
constexpr amrex::Real kObservedOrderErrorFloor = amrex::Real(1.0e-12);
#endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real evaluate_weno_in_x (const int i,
                                const amrex::Array4<const amrex::Real>& qty,
                                const AdvType adv_type,
                                const int sigma) noexcept
{
    const amrex::Real upw = amrex::Real(sigma);
    amrex::Real value = amrex::Real(0.0);

    if (adv_type == AdvType::Weno_3) {
        WENO3 policy(qty, amrex::Real(1.0));
        policy.InterpolateInX(i, 0, 0, kQtyComp, value, upw);
    } else if (adv_type == AdvType::Weno_5) {
        WENO5 policy(qty, amrex::Real(1.0));
        policy.InterpolateInX(i, 0, 0, kQtyComp, value, upw);
    } else if (adv_type == AdvType::Weno_7) {
        WENO7 policy(qty, amrex::Real(1.0));
        policy.InterpolateInX(i, 0, 0, kQtyComp, value, upw);
    } else if (adv_type == AdvType::Weno_3Z) {
        WENO_Z3 policy(qty, amrex::Real(1.0));
        policy.InterpolateInX(i, 0, 0, kQtyComp, value, upw);
    } else if (adv_type == AdvType::Weno_3MZQ) {
        WENO_MZQ3 policy(qty, amrex::Real(1.0));
        policy.InterpolateInX(i, 0, 0, kQtyComp, value, upw);
    } else if (adv_type == AdvType::Weno_5Z) {
        WENO_Z5 policy(qty, amrex::Real(1.0));
        policy.InterpolateInX(i, 0, 0, kQtyComp, value, upw);
    } else {
        WENO_Z7 policy(qty, amrex::Real(1.0));
        policy.InterpolateInX(i, 0, 0, kQtyComp, value, upw);
    }

    return value;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real evaluate_weno_in_y (const int j,
                                const amrex::Array4<const amrex::Real>& qty,
                                const AdvType adv_type,
                                const int sigma) noexcept
{
    const amrex::Real upw = amrex::Real(sigma);
    amrex::Real value = amrex::Real(0.0);

    if (adv_type == AdvType::Weno_3) {
        WENO3 policy(qty, amrex::Real(1.0));
        policy.InterpolateInY(0, j, 0, kQtyComp, value, upw);
    } else if (adv_type == AdvType::Weno_5) {
        WENO5 policy(qty, amrex::Real(1.0));
        policy.InterpolateInY(0, j, 0, kQtyComp, value, upw);
    } else if (adv_type == AdvType::Weno_7) {
        WENO7 policy(qty, amrex::Real(1.0));
        policy.InterpolateInY(0, j, 0, kQtyComp, value, upw);
    } else if (adv_type == AdvType::Weno_3Z) {
        WENO_Z3 policy(qty, amrex::Real(1.0));
        policy.InterpolateInY(0, j, 0, kQtyComp, value, upw);
    } else if (adv_type == AdvType::Weno_3MZQ) {
        WENO_MZQ3 policy(qty, amrex::Real(1.0));
        policy.InterpolateInY(0, j, 0, kQtyComp, value, upw);
    } else if (adv_type == AdvType::Weno_5Z) {
        WENO_Z5 policy(qty, amrex::Real(1.0));
        policy.InterpolateInY(0, j, 0, kQtyComp, value, upw);
    } else {
        WENO_Z7 policy(qty, amrex::Real(1.0));
        policy.InterpolateInY(0, j, 0, kQtyComp, value, upw);
    }

    return value;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real evaluate_weno_in_z (const int k,
                                const amrex::Array4<const amrex::Real>& qty,
                                const AdvType adv_type,
                                const int sigma) noexcept
{
    const amrex::Real upw = amrex::Real(sigma);
    amrex::Real value = amrex::Real(0.0);

    if (adv_type == AdvType::Weno_3) {
        WENO3 policy(qty, amrex::Real(1.0));
        policy.InterpolateInZ(0, 0, k, kQtyComp, value, upw);
    } else if (adv_type == AdvType::Weno_5) {
        WENO5 policy(qty, amrex::Real(1.0));
        policy.InterpolateInZ(0, 0, k, kQtyComp, value, upw);
    } else if (adv_type == AdvType::Weno_7) {
        WENO7 policy(qty, amrex::Real(1.0));
        policy.InterpolateInZ(0, 0, k, kQtyComp, value, upw);
    } else if (adv_type == AdvType::Weno_3Z) {
        WENO_Z3 policy(qty, amrex::Real(1.0));
        policy.InterpolateInZ(0, 0, k, kQtyComp, value, upw);
    } else if (adv_type == AdvType::Weno_3MZQ) {
        WENO_MZQ3 policy(qty, amrex::Real(1.0));
        policy.InterpolateInZ(0, 0, k, kQtyComp, value, upw);
    } else if (adv_type == AdvType::Weno_5Z) {
        WENO_Z5 policy(qty, amrex::Real(1.0));
        policy.InterpolateInZ(0, 0, k, kQtyComp, value, upw);
    } else {
        WENO_Z7 policy(qty, amrex::Real(1.0));
        policy.InterpolateInZ(0, 0, k, kQtyComp, value, upw);
    }

    return value;
}

amrex::Real periodic_l2_error (const AdvType adv_type,
                               const int num_cells)
{
    const amrex::Box qty_box(amrex::IntVect(-4, 0, 0), amrex::IntVect(num_cells + 4, 0, 0));
    const amrex::Box out_box(amrex::IntVect(0, 0, 0), amrex::IntVect(num_cells - 1, 0, 0));
    amrex::FArrayBox qty(qty_box, 1);
    amrex::FArrayBox output(out_box, 1);

    auto qty_arr = qty.array();
    amrex::ParallelFor(qty_box, [=] AMREX_GPU_DEVICE (int i, int, int) noexcept {
        int periodic_index = i % num_cells;
        if (periodic_index < 0) periodic_index += num_cells;
        qty_arr(i, 0, 0, 0) = smooth_periodic_cell_average(periodic_index, num_cells);
    });

    const auto qty_const_arr = qty.const_array();
    auto out_arr = output.array();
    amrex::ParallelFor(out_box, [=] AMREX_GPU_DEVICE (int i, int, int) noexcept {
        out_arr(i, 0, 0, 0) = evaluate_weno_in_x(i, qty_const_arr, adv_type, 0);
    });
    gpu_sync();

    // Use ReduceOps since we're operating on a single FArrayBox
    const auto out_const_arr = output.const_array();
    amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);

    reduce_op.eval(out_box, reduce_data,
    [=] AMREX_GPU_DEVICE (int i, int, int) -> amrex::GpuTuple<amrex::Real> {
        const amrex::Real diff = out_const_arr(i, 0, 0, 0)
            - smooth_periodic_face_value(i, num_cells);
        return diff * diff;
    });

    amrex::Real error_sum = amrex::get<0>(reduce_data.value());

    return std::sqrt(error_sum / amrex::Real(num_cells));
}

inline amrex::Box wrapper_branch_dense_box ()
{
    return amrex::Box(amrex::IntVect(-7, -7, -7), amrex::IntVect(7, 7, 7));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real branch_dense_profile (const int index)
{
    return smooth_periodic_cell_average(index + 9, 32) + amrex::Real(0.05) * amrex::Real(index * index - 3 * index);
}

void fill_branch_dense_profile_3d (const amrex::Box& box,
                                   const amrex::Array4<amrex::Real>& qty_x,
                                   const amrex::Array4<amrex::Real>& qty_y,
                                   const amrex::Array4<amrex::Real>& qty_z)
{
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        qty_x(i, j, k, 0) = branch_dense_profile(i);
        qty_y(i, j, k, 0) = branch_dense_profile(j);
        qty_z(i, j, k, 0) = branch_dense_profile(k);
    });
    gpu_sync();
}

void weno7_fill_qty_with_index (const amrex::Box& box,
                                const amrex::Array4<amrex::Real>& arr)
{
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int, int, int k) noexcept {
        arr(0, 0, k, 0) = amrex::Real(k);
    });
    gpu_sync();
}

void launch_weno7_z_indexing_case (const amrex::Array4<const amrex::Real>& qty,
                                   const amrex::Array4<amrex::Real>& output)
{
    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int) noexcept {
        output(0, 0, 0, 0) = evaluate_weno_in_z(0, qty, AdvType::Weno_7, 1);
        output(0, 0, 0, 1) = evaluate_weno_in_z(0, qty, AdvType::Weno_7, -1);
    });
    gpu_sync();
}

void launch_weno7_directional_wrapper_cases (const amrex::Box& out_box,
                                             const amrex::Array4<const amrex::Real>& qty_x,
                                             const amrex::Array4<const amrex::Real>& qty_y,
                                             const amrex::Array4<const amrex::Real>& qty_z,
                                             const amrex::Array4<amrex::Real>& output)
{
    amrex::ParallelFor(out_box, 6, [=] AMREX_GPU_DEVICE (int i, int, int, int n) noexcept {
        const int sigma = (n / 3 == 0) ? -1 : 1;
        const int axis = n % 3;

        if (axis == 0) {
            output(i, 0, 0, n) = evaluate_weno_in_x(i + 1, qty_x, AdvType::Weno_7, sigma);
        } else if (axis == 1) {
            output(i, 0, 0, n) = evaluate_weno_in_y(i + 1, qty_y, AdvType::Weno_7, sigma);
        } else {
            output(i, 0, 0, n) = evaluate_weno_in_z(i + 1, qty_z, AdvType::Weno_7, sigma);
        }
    });
    gpu_sync();
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
AdvType wenoz_directional_wrapper_scheme (const int group) noexcept
{
    return (group == 0) ? AdvType::Weno_3Z :
           (group == 1) ? AdvType::Weno_3MZQ :
           (group == 2) ? AdvType::Weno_5Z :
                          AdvType::Weno_7Z;
}

void launch_wenoz_directional_wrapper_cases (const amrex::Box& out_box,
                                             const int ncomp,
                                             const amrex::Array4<const amrex::Real>& qty_x,
                                             const amrex::Array4<const amrex::Real>& qty_y,
                                             const amrex::Array4<const amrex::Real>& qty_z,
                                             const amrex::Array4<amrex::Real>& output)
{
    amrex::ParallelFor(out_box, ncomp, [=] AMREX_GPU_DEVICE (int i, int, int, int n) noexcept {
        const int group = n / 3;
        const int axis = n % 3;
        const AdvType adv_type = wenoz_directional_wrapper_scheme(group / 2);
        const int sigma = (group % 2 == 0) ? -1 : 1;

        if (axis == 0) {
            output(i, 0, 0, n) = evaluate_weno_in_x(i + 1, qty_x, adv_type, sigma);
        } else if (axis == 1) {
            output(i, 0, 0, n) = evaluate_weno_in_y(i + 1, qty_y, adv_type, sigma);
        } else {
            output(i, 0, 0, n) = evaluate_weno_in_z(i + 1, qty_z, adv_type, sigma);
        }
    });
    gpu_sync();
}

} // namespace

// Motivation: Smooth periodic convergence is the kernel-side end-to-end check
// that the classic WENO wrappers recover their finite-volume design order.
TEST(InterpolationWENOKernel, ClassicWenoObservedConvergenceOnSmoothPeriodicField)
{
    const std::array<AdvType, 3> schemes = {{AdvType::Weno_3, AdvType::Weno_5, AdvType::Weno_7}};
    const std::array<int, 4> sizes = {{128, 256, 512, 1024}};

    for (const AdvType adv_type : schemes) {
        std::array<amrex::Real, 4> errors = {{amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0)}};
        for (int idx = 0; idx < static_cast<int>(sizes.size()); ++idx) {
            errors[idx] = periodic_l2_error(adv_type, sizes[idx]);
        }

        const amrex::Real threshold =
            (adv_type == AdvType::Weno_3) ? kWeno3ObservedOrderMin :
            (adv_type == AdvType::Weno_5) ? kWeno5ObservedOrderMin :
                                            kWeno7ObservedOrderMin;
        const int pair_count = (adv_type == AdvType::Weno_3) ? 2 : static_cast<int>(sizes.size()) - 1;

        for (int idx = 0; idx < pair_count; ++idx) {
            if ((errors[idx] < kObservedOrderErrorFloor) || (errors[idx + 1] < kObservedOrderErrorFloor)) {
                break;
            }
            EXPECT_GE(observed_order(errors[idx], errors[idx + 1]), threshold)
                << "scheme=" << scheme_name(adv_type)
                << " coarse_N=" << sizes[idx]
                << " coarse_error=" << static_cast<double>(errors[idx])
                << " fine_error=" << static_cast<double>(errors[idx + 1]);
        }
    }
}

// Motivation: This regression locks the WENO7 Z low-face stencil indexing so a
// repeated k+2 or k-3 load cannot silently return.
TEST(InterpolationWENOKernel, Weno7DirectionalLowFaceIndexingInZ)
{
    const amrex::Box qty_box(amrex::IntVect(0, 0, -4), amrex::IntVect(0, 0, 4));
    amrex::FArrayBox qty(qty_box, 1);
    weno7_fill_qty_with_index(qty_box, qty.array());

    const auto qty_const_arr = qty.const_array();
    amrex::FArrayBox output(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 0)), 2);
    launch_weno7_z_indexing_case(qty_const_arr, output.array());

    amrex::FArrayBox host_output(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 0)),
                                 2, amrex::The_Pinned_Arena());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     output.dataPtr(0), output.dataPtr(0) + output.size(),
                     host_output.dataPtr(0));

    const auto out_const_arr = host_output.const_array();
    EXPECT_NEAR(out_const_arr(0, 0, 0, 0), amrex::Real(-0.5),
                scaled_tol(out_const_arr(0, 0, 0, 0), amrex::Real(-0.5), kWenoDeviceRelTol));
    EXPECT_NEAR(out_const_arr(0, 0, 0, 1), amrex::Real(-0.5),
                scaled_tol(out_const_arr(0, 0, 0, 1), amrex::Real(-0.5), kWenoDeviceRelTol));
}

// Motivation: WENO7 should be directionally identical when the data profile is
// the same along x, y, and z.
TEST(InterpolationWENOKernel, Weno7DirectionalWrappersAreEquivalentOnBranchDenseProfiles)
{
    const amrex::Box qty_box = wrapper_branch_dense_box();
    const amrex::Box out_box(amrex::IntVect(0, 0, 0), amrex::IntVect(3, 0, 0));
    amrex::FArrayBox qty_x(qty_box, 1);
    amrex::FArrayBox qty_y(qty_box, 1);
    amrex::FArrayBox qty_z(qty_box, 1);

    auto qty_x_arr = qty_x.array();
    auto qty_y_arr = qty_y.array();
    auto qty_z_arr = qty_z.array();
    fill_branch_dense_profile_3d(qty_box, qty_x_arr, qty_y_arr, qty_z_arr);

    const auto qty_x_const_arr = qty_x.const_array();
    const auto qty_y_const_arr = qty_y.const_array();
    const auto qty_z_const_arr = qty_z.const_array();
    amrex::FArrayBox output(out_box, 6);
    launch_weno7_directional_wrapper_cases(out_box, qty_x_const_arr, qty_y_const_arr,
                                           qty_z_const_arr, output.array());

    amrex::FArrayBox host_output(out_box, 6, amrex::The_Pinned_Arena());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     output.dataPtr(0), output.dataPtr(0) + output.size(),
                     host_output.dataPtr(0));

    const auto out_const_arr = host_output.const_array();
    for (int sigma_index = 0; sigma_index < 2; ++sigma_index) {
        const int sigma = (sigma_index == 0) ? -1 : 1;
        for (int i = 0; i < 4; ++i) {
            const int base = 3 * sigma_index;
            const amrex::Real x_value = out_const_arr(i, 0, 0, base + 0);
            const amrex::Real y_value = out_const_arr(i, 0, 0, base + 1);
            const amrex::Real z_value = out_const_arr(i, 0, 0, base + 2);

            EXPECT_NEAR(x_value, y_value, scaled_tol(x_value, y_value, kWenoDeviceRelTol))
                << "scheme=" << scheme_name(AdvType::Weno_7) << " sigma=" << sigma << " face=" << i;
            EXPECT_NEAR(x_value, z_value, scaled_tol(x_value, z_value, kWenoDeviceRelTol))
                << "scheme=" << scheme_name(AdvType::Weno_7) << " sigma=" << sigma << " face=" << i;
        }
    }
}

// Motivation: The WENO-Z wrappers should share the same directional behavior on
// branch-dense inputs, independent of the axis-specific wrapper path.
TEST(InterpolationWENOKernel, WenoZDirectionalWrappersAreEquivalentOnBranchDenseProfiles)
{
    const std::array<AdvType, 4> schemes = {{AdvType::Weno_3Z, AdvType::Weno_3MZQ,
                                             AdvType::Weno_5Z, AdvType::Weno_7Z}};
    const amrex::Box qty_box = wrapper_branch_dense_box();
    const amrex::Box out_box(amrex::IntVect(0, 0, 0), amrex::IntVect(3, 0, 0));
    amrex::FArrayBox qty_x(qty_box, 1);
    amrex::FArrayBox qty_y(qty_box, 1);
    amrex::FArrayBox qty_z(qty_box, 1);

    auto qty_x_arr = qty_x.array();
    auto qty_y_arr = qty_y.array();
    auto qty_z_arr = qty_z.array();
    fill_branch_dense_profile_3d(qty_box, qty_x_arr, qty_y_arr, qty_z_arr);

    const auto qty_x_const_arr = qty_x.const_array();
    const auto qty_y_const_arr = qty_y.const_array();
    const auto qty_z_const_arr = qty_z.const_array();
    const int ncomp = 2 * static_cast<int>(schemes.size()) * 3;
    amrex::FArrayBox output(out_box, ncomp);
    launch_wenoz_directional_wrapper_cases(out_box, ncomp, qty_x_const_arr, qty_y_const_arr,
                                           qty_z_const_arr, output.array());

    amrex::FArrayBox host_output(out_box, ncomp, amrex::The_Pinned_Arena());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     output.dataPtr(0), output.dataPtr(0) + output.size(),
                     host_output.dataPtr(0));

    const auto out_const_arr = host_output.const_array();
    for (int n = 0; n < 2 * static_cast<int>(schemes.size()); ++n) {
        const AdvType adv_type = schemes[n / 2];
        const int sigma = (n % 2 == 0) ? -1 : 1;
        for (int i = 0; i < 4; ++i) {
            const int base = 3 * n;
            const amrex::Real x_value = out_const_arr(i, 0, 0, base + 0);
            const amrex::Real y_value = out_const_arr(i, 0, 0, base + 1);
            const amrex::Real z_value = out_const_arr(i, 0, 0, base + 2);

            EXPECT_NEAR(x_value, y_value, scaled_tol(x_value, y_value, kWenoDeviceRelTol))
                << "scheme=" << scheme_name(adv_type) << " sigma=" << sigma << " face=" << i;
            EXPECT_NEAR(x_value, z_value, scaled_tol(x_value, z_value, kWenoDeviceRelTol))
                << "scheme=" << scheme_name(adv_type) << " sigma=" << sigma << " face=" << i;
        }
    }
}

// Motivation: This is the kernel-side smooth-order check for the WENO-Z
// families after the tau and beta contracts are applied in production code.
TEST(InterpolationWENOKernel, WenoZObservedConvergenceOnSmoothPeriodicField)
{
    const std::array<AdvType, 4> schemes = {{AdvType::Weno_3Z, AdvType::Weno_3MZQ,
                                             AdvType::Weno_5Z, AdvType::Weno_7Z}};
    const std::array<int, 4> sizes = {{128, 256, 512, 1024}};

    for (const AdvType adv_type : schemes) {
        std::array<amrex::Real, 4> errors = {{amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0)}};
        for (int idx = 0; idx < static_cast<int>(sizes.size()); ++idx) {
            errors[idx] = periodic_l2_error(adv_type, sizes[idx]);
        }

        const amrex::Real threshold =
            (adv_type == AdvType::Weno_3Z)   ? kWeno3ObservedOrderMin :
            (adv_type == AdvType::Weno_3MZQ) ? kWenoMZQ3ObservedOrderMin :
            (adv_type == AdvType::Weno_5Z)   ? kWeno5ObservedOrderMin :
                                               kWeno7ObservedOrderMin;
        const int pair_count =
            ((adv_type == AdvType::Weno_3Z) || (adv_type == AdvType::Weno_3MZQ))
                ? 2
                : static_cast<int>(sizes.size()) - 1;

        for (int idx = 0; idx < pair_count; ++idx) {
            if ((errors[idx] < kObservedOrderErrorFloor) || (errors[idx + 1] < kObservedOrderErrorFloor)) {
                break;
            }
            EXPECT_GE(observed_order(errors[idx], errors[idx + 1]), threshold)
                << "scheme=" << scheme_name(adv_type)
                << " coarse_N=" << sizes[idx]
                << " coarse_error=" << static_cast<double>(errors[idx])
                << " fine_error=" << static_cast<double>(errors[idx + 1]);
        }
    }
}
