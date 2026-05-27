#include <array>
#include <string>
#include <vector>

#include <AMReX_FArrayBox.H>
#include <AMReX_GpuContainers.H>

#include <gtest/gtest.h>

#include "ERF_GTestInterpolationCommon.H"

using namespace interpolation_test;

namespace {

struct PolicyWrapperCase {
    AdvType adv_type;
    Coord direction;
    int upwind;
};

std::vector<PolicyWrapperCase> make_policy_wrapper_cases ()
{
    std::vector<PolicyWrapperCase> cases;
    const std::array<Coord, 3> directions = {{Coord::x, Coord::y, Coord::z}};
    const std::array<AdvType, 6> schemes = {{
        AdvType::Centered_2nd,
        AdvType::Upwind_3rd,
        AdvType::Upwind_3rd_SL,
        AdvType::Centered_4th,
        AdvType::Upwind_5th,
        AdvType::Centered_6th}};

    for (const AdvType adv_type : schemes) {
        for (const int upwind : kRepresentativeUpwindSigns) {
            for (const Coord direction : directions) {
                cases.push_back(PolicyWrapperCase{adv_type, direction, upwind});
            }
        }
    }

    return cases;
}

std::vector<PolicyWrapperCase> make_upwindall_cases ()
{
    std::vector<PolicyWrapperCase> cases;

    for (const AdvType adv_type : kSupportedAdvTypes) {
        for (const int upwind : kRepresentativeUpwindSigns) {
            cases.push_back(PolicyWrapperCase{adv_type, Coord::z, upwind});
        }
    }

    return cases;
}

void fill_directional_qty (amrex::FArrayBox& qty)
{
    auto qty_arr = qty.array();
    const amrex::Box box = qty.box();

    for (int k = box.smallEnd(2); k <= box.bigEnd(2); ++k) {
        for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
            for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
                qty_arr(i, j, k, 0) = directional_field_value(i, j, k);
            }
        }
    }
}

std::string sample_trace (const AdvType adv_type,
                          const Coord direction,
                          const int upwind,
                          const int sample_index,
                          const amrex::Real actual,
                          const amrex::Real expected,
                          const amrex::Real factor)
{
    return "scheme=" + scheme_name(adv_type) +
           " direction=" + direction_name(direction) +
           " upwind=" + std::to_string(upwind) +
           " sample=" + std::to_string(sample_index) +
           " actual=" + std::to_string(static_cast<double>(actual)) +
           " expected=" + std::to_string(static_cast<double>(expected)) +
           " normalized_error=" + std::to_string(static_cast<double>(normalized_error(actual, expected, factor)));
}

void launch_policy_struct_cases (const amrex::Box& box,
                                 const amrex::Array4<const amrex::Real>& qty,
                                 const PolicyWrapperCase* cases_ptr,
                                 const int ncases,
                                 const amrex::Array4<amrex::Real>& output)
{
    amrex::ParallelFor(box, ncases, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
        const PolicyWrapperCase test_case = cases_ptr[n];
        const amrex::Real upw = amrex::Real(test_case.upwind);
        amrex::Real value = amrex::Real(0.0);

        if (test_case.adv_type == AdvType::Centered_2nd) {
            CENTERED2 policy(qty, amrex::Real(1.0));
            if (test_case.direction == Coord::x) policy.InterpolateInX(i, j, k, kQtyComp, value, upw);
            else if (test_case.direction == Coord::y) policy.InterpolateInY(i, j, k, kQtyComp, value, upw);
            else policy.InterpolateInZ(i, j, k, kQtyComp, value, upw);
        } else if (test_case.adv_type == AdvType::Upwind_3rd) {
            UPWIND3 policy(qty, amrex::Real(1.0));
            if (test_case.direction == Coord::x) policy.InterpolateInX(i, j, k, kQtyComp, value, upw);
            else if (test_case.direction == Coord::y) policy.InterpolateInY(i, j, k, kQtyComp, value, upw);
            else policy.InterpolateInZ(i, j, k, kQtyComp, value, upw);
        } else if (test_case.adv_type == AdvType::Upwind_3rd_SL) {
            UPWIND3SL policy(qty, amrex::Real(1.0));
            if (test_case.direction == Coord::x) policy.InterpolateInX(i, j, k, kQtyComp, value, upw);
            else if (test_case.direction == Coord::y) policy.InterpolateInY(i, j, k, kQtyComp, value, upw);
            else policy.InterpolateInZ(i, j, k, kQtyComp, value, upw);
        } else if (test_case.adv_type == AdvType::Centered_4th) {
            CENTERED4 policy(qty, amrex::Real(1.0));
            if (test_case.direction == Coord::x) policy.InterpolateInX(i, j, k, kQtyComp, value, upw);
            else if (test_case.direction == Coord::y) policy.InterpolateInY(i, j, k, kQtyComp, value, upw);
            else policy.InterpolateInZ(i, j, k, kQtyComp, value, upw);
        } else if (test_case.adv_type == AdvType::Upwind_5th) {
            UPWIND5 policy(qty, amrex::Real(1.0));
            if (test_case.direction == Coord::x) policy.InterpolateInX(i, j, k, kQtyComp, value, upw);
            else if (test_case.direction == Coord::y) policy.InterpolateInY(i, j, k, kQtyComp, value, upw);
            else policy.InterpolateInZ(i, j, k, kQtyComp, value, upw);
        } else {
            CENTERED6 policy(qty, amrex::Real(1.0));
            if (test_case.direction == Coord::x) policy.InterpolateInX(i, j, k, kQtyComp, value, upw);
            else if (test_case.direction == Coord::y) policy.InterpolateInY(i, j, k, kQtyComp, value, upw);
            else policy.InterpolateInZ(i, j, k, kQtyComp, value, upw);
        }

        output(i, j, k, n) = value;
    });

    gpu_sync();
}

void launch_upwindall_cases (const amrex::Box& box,
                             const amrex::Array4<const amrex::Real>& qty,
                             const PolicyWrapperCase* cases_ptr,
                             const int ncases,
                             const amrex::Array4<amrex::Real>& output)
{
    amrex::ParallelFor(box, ncases, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
        const PolicyWrapperCase test_case = cases_ptr[n];
        UPWINDALL policy(qty, amrex::Real(1.0));
        amrex::Real value = amrex::Real(0.0);
        policy.InterpolateInZ(i, j, k, kQtyComp, value, amrex::Real(test_case.upwind), test_case.adv_type);
        output(i, j, k, n) = value;
    });

    gpu_sync();
}

} // namespace

// Motivation: The policy structs are the hot-path interpolation implementations.
// Each direct device wrapper must honor the same finite-volume contract as the
// independent host reference, including the slope-limited UPWIND3SL branch.
TEST(InterpolationPolicyKernel, PolicyStructWrappersMatchIndependentReference)
{
    const std::vector<PolicyWrapperCase> cases = make_policy_wrapper_cases();
    amrex::FArrayBox qty(stencil_box(), 1);
    amrex::FArrayBox output(eval_box(), static_cast<int>(cases.size()));
    fill_directional_qty(qty);

    amrex::Gpu::DeviceVector<PolicyWrapperCase> device_cases(cases.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, cases.begin(), cases.end(), device_cases.begin());
    launch_policy_struct_cases(eval_box(), qty.const_array(), device_cases.data(),
                               static_cast<int>(cases.size()), output.array());

    const auto qty_arr = qty.const_array();
    const auto output_arr = output.const_array();
    int sample_index = 0;

    for (int n = 0; n < static_cast<int>(cases.size()); ++n) {
        for (int k = eval_box().smallEnd(2); k <= eval_box().bigEnd(2); ++k) {
            for (int j = eval_box().smallEnd(1); j <= eval_box().bigEnd(1); ++j) {
                for (int i = eval_box().smallEnd(0); i <= eval_box().bigEnd(0); ++i) {
                    const PolicyWrapperCase& test_case = cases[n];
                    const StencilValues stencil = load_direction_stencil(qty_arr, i, j, k, test_case.direction);
                    const amrex::Real expected =
                        (test_case.adv_type == AdvType::Upwind_3rd_SL)
                        ? upwind3sl_reference_flux(stencil.q1, stencil.q0, stencil.qm1, stencil.qm2, test_case.upwind)
                        : reference_reconstruction(stencil, test_case.adv_type, test_case.upwind);
                    const amrex::Real actual = output_arr(i, j, k, n);

                    EXPECT_LE(normalized_error(actual, expected, kDeviceRelTol), amrex::Real(1.0))
                        << sample_trace(test_case.adv_type, test_case.direction, test_case.upwind,
                                        sample_index, actual, expected, kDeviceRelTol);
                    ++sample_index;
                }
            }
        }
    }
}

// Motivation: UPWINDALL consolidates the centered/upwind Z-face policies. The
// device path must remain identical to the dedicated helper for each AdvType.
TEST(InterpolationPolicyKernel, UpwindAllMatchesDedicatedZPolicy)
{
    const std::vector<PolicyWrapperCase> cases = make_upwindall_cases();
    amrex::FArrayBox qty(stencil_box(), 1);
    amrex::FArrayBox output(eval_box(), static_cast<int>(cases.size()));
    fill_directional_qty(qty);

    amrex::Gpu::DeviceVector<PolicyWrapperCase> device_cases(cases.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, cases.begin(), cases.end(), device_cases.begin());
    launch_upwindall_cases(eval_box(), qty.const_array(), device_cases.data(),
                           static_cast<int>(cases.size()), output.array());

    const auto qty_arr = qty.const_array();
    const auto output_arr = output.const_array();
    int sample_index = 0;

    for (int n = 0; n < static_cast<int>(cases.size()); ++n) {
        for (int k = eval_box().smallEnd(2); k <= eval_box().bigEnd(2); ++k) {
            for (int j = eval_box().smallEnd(1); j <= eval_box().bigEnd(1); ++j) {
                for (int i = eval_box().smallEnd(0); i <= eval_box().bigEnd(0); ++i) {
                    const PolicyWrapperCase& test_case = cases[n];
                    const StencilValues stencil = load_direction_stencil(qty_arr, i, j, k, Coord::z);
                    const amrex::Real expected = reference_reconstruction(stencil, test_case.adv_type, test_case.upwind);
                    const amrex::Real actual = output_arr(i, j, k, n);

                    EXPECT_LE(normalized_error(actual, expected, kDeviceRelTol), amrex::Real(1.0))
                        << sample_trace(test_case.adv_type, Coord::z, test_case.upwind,
                                        sample_index, actual, expected, kDeviceRelTol);
                    ++sample_index;
                }
            }
        }
    }
}