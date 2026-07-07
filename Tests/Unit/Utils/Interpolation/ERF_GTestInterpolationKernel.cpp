#include <array>
#include <string>
#include <vector>

#include <AMReX_FArrayBox.H>
#include <AMReX_GpuContainers.H>

#include <gtest/gtest.h>

#include "ERF_GTestInterpolationCommon.H"

using namespace interpolation_test;

namespace {

struct InterpolatedValCase {
    StencilValues stencil;
    AdvType adv_type;
    int upwind;
    int sample_index;
};

struct WrapperCase {
    AdvType adv_type;
    Coord direction;
    int upwind;
};

struct SentinelCase {
    AdvType adv_type;
    Coord direction;
    int upwind;
};

std::array<StencilValues, 3> asymmetric_stencils ()
{
    return {{
        {amrex::Real(1.0),  amrex::Real(-2.5), amrex::Real(0.75),  amrex::Real(3.25),  amrex::Real(-1.125), amrex::Real(2.0)},
        {amrex::Real(-4.0), amrex::Real(0.5),  amrex::Real(7.0),   amrex::Real(-3.0),  amrex::Real(1.5),    amrex::Real(2.25)},
        {amrex::Real(2.0),  amrex::Real(2.5),  amrex::Real(-1.0),  amrex::Real(4.5),   amrex::Real(-3.75),  amrex::Real(0.125)}}};
}

std::vector<InterpolatedValCase> make_interpolated_val_cases ()
{
    std::vector<InterpolatedValCase> cases;
    const auto stencils = asymmetric_stencils();

    for (int sample_index = 0; sample_index < static_cast<int>(stencils.size()); ++sample_index) {
        for (const AdvType adv_type : kSupportedAdvTypes) {
            for (const int upwind : kRepresentativeUpwindSigns) {
                cases.push_back(InterpolatedValCase{stencils[sample_index], adv_type, upwind, sample_index});
            }
        }
    }

    return cases;
}

std::vector<WrapperCase> make_wrapper_cases ()
{
    std::vector<WrapperCase> cases;
    const std::array<Coord, 3> directions = {{Coord::x, Coord::y, Coord::z}};

    for (const AdvType adv_type : kSupportedAdvTypes) {
        for (const int upwind : kRepresentativeUpwindSigns) {
            for (const Coord direction : directions) {
                cases.push_back(WrapperCase{adv_type, direction, upwind});
            }
        }
    }

    return cases;
}

std::vector<WrapperCase> make_zero_upwind_wrapper_cases ()
{
    std::vector<WrapperCase> cases;
    const std::array<Coord, 3> directions = {{Coord::x, Coord::y, Coord::z}};

    for (const Coord direction : directions) {
        cases.push_back(WrapperCase{AdvType::Upwind_3rd, direction, 0});
        cases.push_back(WrapperCase{AdvType::Upwind_5th, direction, 0});
    }

    return cases;
}

std::vector<SentinelCase> make_lower_order_sentinel_cases ()
{
    std::vector<SentinelCase> cases;
    const std::array<Coord, 3> directions = {{Coord::x, Coord::y, Coord::z}};
    const std::array<AdvType, 3> schemes = {{AdvType::Centered_2nd, AdvType::Upwind_3rd, AdvType::Centered_4th}};

    for (const AdvType adv_type : schemes) {
        for (const int upwind : kRepresentativeUpwindSigns) {
            for (const Coord direction : directions) {
                cases.push_back(SentinelCase{adv_type, direction, upwind});
            }
        }
    }

    return cases;
}

std::vector<SentinelCase> make_high_order_sentinel_cases ()
{
    std::vector<SentinelCase> cases;
    const std::array<Coord, 3> directions = {{Coord::x, Coord::y, Coord::z}};
    const std::array<AdvType, 2> schemes = {{AdvType::Upwind_5th, AdvType::Centered_6th}};

    for (const AdvType adv_type : schemes) {
        for (const int upwind : kRepresentativeUpwindSigns) {
            for (const Coord direction : directions) {
                cases.push_back(SentinelCase{adv_type, direction, upwind});
            }
        }
    }

    return cases;
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

void fill_directional_qty (amrex::FArrayBox& qty)
{
    auto qty_arr = qty.array();
    const amrex::Box box = qty.box();

    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        qty_arr(i, j, k, 0) = directional_field_value(i, j, k);
    });
    gpu_sync();
}

void fill_perturbation_boxes (amrex::FArrayBox& qty,
                              amrex::FArrayBox& r0,
                              const Coord direction)
{
    auto qty_arr = qty.array();
    auto r0_arr = r0.array();
    const amrex::Box box = qty.box();

    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        const amrex::Real background = perturbation_background(direction, i, j, k);
        r0_arr(i, j, k, 0) = background;
        qty_arr(i, j, k, 0) = background + kPerturbationConstant;
    });
    gpu_sync();
}

void fill_density_boxes (amrex::FArrayBox& cons_in,
                         amrex::FArrayBox& r0)
{
    auto cons_arr = cons_in.array();
    auto r0_arr = r0.array();
    const amrex::Box box = cons_in.box();

    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        const amrex::Real background = amrex::Real(1.0) +
            amrex::Real(0.5) * finite_volume_cell_average_monomial(4, i) -
            amrex::Real(0.25) * finite_volume_cell_average_monomial(3, j) +
            amrex::Real(0.125) * finite_volume_cell_average_monomial(2, k);

        r0_arr(i, j, k, 0) = background;
        cons_arr(i, j, k, 0) = background + density_quantity_value(i, j, k);
    });
    gpu_sync();
}

void apply_far_sentinels (amrex::FArrayBox& qty)
{
    auto qty_arr = qty.array();

    amrex::single_task([=] AMREX_GPU_DEVICE () noexcept {
        qty_arr( 2, 0, 0, 0) += amrex::Real(1000.0);
        qty_arr(-3, 0, 0, 0) -= amrex::Real(750.0);
        qty_arr( 0, 2, 0, 0) += amrex::Real(900.0);
        qty_arr( 0,-3, 0, 0) -= amrex::Real(650.0);
        qty_arr( 0, 0, 2, 0) += amrex::Real(800.0);
        qty_arr( 0, 0,-3, 0) -= amrex::Real(550.0);
    });
    gpu_sync();
}

void launch_interpolated_val_cases (const int ncases,
                                    const InterpolatedValCase* cases_ptr,
                                    amrex::Real* results_ptr)
{
    amrex::ParallelFor(ncases, [=] AMREX_GPU_DEVICE (int idx) noexcept {
        const InterpolatedValCase test_case = cases_ptr[idx];
        amrex::Real avg1 = amrex::Real(0.0);
        amrex::Real avg2 = amrex::Real(0.0);
        amrex::Real avg3 = amrex::Real(0.0);
        amrex::Real diff1 = amrex::Real(0.0);
        amrex::Real diff2 = amrex::Real(0.0);
        amrex::Real diff3 = amrex::Real(0.0);

        fill_avg_diff(test_case.stencil, avg1, avg2, avg3, diff1, diff2, diff3);
        results_ptr[idx] = interpolatedVal(avg1, avg2, avg3, diff1, diff2, diff3,
                                           amrex::Real(test_case.upwind), test_case.adv_type);
    });

    gpu_sync();
}

void launch_wrapper_cases (const amrex::Box& box,
                           const amrex::Array4<const amrex::Real>& qty,
                           const WrapperCase* cases_ptr,
                           const int ncases,
                           const amrex::Array4<amrex::Real>& output)
{
    amrex::ParallelFor(box, ncases, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
        const WrapperCase test_case = cases_ptr[n];

        if (test_case.direction == Coord::x) {
            output(i, j, k, n) = InterpolateInX(i, j, k, qty, kQtyComp,
                                                amrex::Real(test_case.upwind), test_case.adv_type);
        } else if (test_case.direction == Coord::y) {
            output(i, j, k, n) = InterpolateInY(i, j, k, qty, kQtyComp,
                                                amrex::Real(test_case.upwind), test_case.adv_type);
        } else {
            output(i, j, k, n) = InterpolateInZ(i, j, k, qty, kQtyComp,
                                                amrex::Real(test_case.upwind), test_case.adv_type);
        }
    });

    gpu_sync();
}

void launch_perturbation_cases (const amrex::Box& box,
                                const WrapperCase* cases_ptr,
                                const int ncases,
                                const amrex::Array4<const amrex::Real>& qty_x,
                                const amrex::Array4<const amrex::Real>& r0_x,
                                const amrex::Array4<const amrex::Real>& qty_y,
                                const amrex::Array4<const amrex::Real>& r0_y,
                                const amrex::Array4<const amrex::Real>& qty_z,
                                const amrex::Array4<const amrex::Real>& r0_z,
                                const amrex::Array4<amrex::Real>& output)
{
    amrex::ParallelFor(box, ncases, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
        const WrapperCase test_case = cases_ptr[n];

        if (test_case.direction == Coord::x) {
            output(i, j, k, n) = InterpolatePertFromCell(i, j, k, qty_x, kQtyComp,
                                                         amrex::Real(test_case.upwind), Coord::x,
                                                         test_case.adv_type, r0_x);
        } else if (test_case.direction == Coord::y) {
            output(i, j, k, n) = InterpolatePertFromCell(i, j, k, qty_y, kQtyComp,
                                                         amrex::Real(test_case.upwind), Coord::y,
                                                         test_case.adv_type, r0_y);
        } else {
            output(i, j, k, n) = InterpolatePertFromCell(i, j, k, qty_z, kQtyComp,
                                                         amrex::Real(test_case.upwind), Coord::z,
                                                         test_case.adv_type, r0_z);
        }
    });

    gpu_sync();
}

void launch_density_identity_cases (const amrex::Box& box,
                                    const WrapperCase* cases_ptr,
                                    const int ncases,
                                    const amrex::Array4<const amrex::Real>& cons_in,
                                    const amrex::Array4<const amrex::Real>& r0,
                                    const amrex::Array4<amrex::Real>& generic_output,
                                    const amrex::Array4<amrex::Real>& density_output)
{
    amrex::ParallelFor(box, ncases, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
        const WrapperCase test_case = cases_ptr[n];

        generic_output(i, j, k, n) = InterpolatePertFromCell(i, j, k, cons_in, Rho_comp,
                                                             amrex::Real(test_case.upwind), test_case.direction,
                                                             test_case.adv_type, r0);
        density_output(i, j, k, n) = InterpolateDensityPertFromCellToFace(i, j, k, cons_in,
                                                                          amrex::Real(test_case.upwind), test_case.direction,
                                                                          test_case.adv_type, r0);
    });

    gpu_sync();
}

void launch_sentinel_cases (const int ncases,
                            const SentinelCase* cases_ptr,
                            const amrex::Array4<const amrex::Real>& base_qty,
                            const amrex::Array4<const amrex::Real>& modified_qty,
                            amrex::Real* base_results,
                            amrex::Real* modified_results)
{
    amrex::ParallelFor(ncases, [=] AMREX_GPU_DEVICE (int idx) noexcept {
        const SentinelCase test_case = cases_ptr[idx];

        if (test_case.direction == Coord::x) {
            base_results[idx] = InterpolateInX(0, 0, 0, base_qty, kQtyComp,
                                               amrex::Real(test_case.upwind), test_case.adv_type);
            modified_results[idx] = InterpolateInX(0, 0, 0, modified_qty, kQtyComp,
                                                   amrex::Real(test_case.upwind), test_case.adv_type);
        } else if (test_case.direction == Coord::y) {
            base_results[idx] = InterpolateInY(0, 0, 0, base_qty, kQtyComp,
                                               amrex::Real(test_case.upwind), test_case.adv_type);
            modified_results[idx] = InterpolateInY(0, 0, 0, modified_qty, kQtyComp,
                                                   amrex::Real(test_case.upwind), test_case.adv_type);
        } else {
            base_results[idx] = InterpolateInZ(0, 0, 0, base_qty, kQtyComp,
                                               amrex::Real(test_case.upwind), test_case.adv_type);
            modified_results[idx] = InterpolateInZ(0, 0, 0, modified_qty, kQtyComp,
                                                   amrex::Real(test_case.upwind), test_case.adv_type);
        }
    });

    gpu_sync();
}

} // namespace

// Motivation: Production kernels assemble avg/diff terms and call
// interpolatedVal on device. This sweep checks every supported scheme and
// representative upwind sign against an independent finite-volume reference.
TEST(InterpolationKernel, InterpolatedValMatchesIndependentReference)
{
    const std::vector<InterpolatedValCase> cases = make_interpolated_val_cases();
    amrex::Gpu::DeviceVector<InterpolatedValCase> device_cases(cases.size());
    amrex::Gpu::DeviceVector<amrex::Real> device_results(cases.size(), amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real> host_results(cases.size(), amrex::Real(0.0));

    amrex::Gpu::copy(amrex::Gpu::hostToDevice, cases.begin(), cases.end(), device_cases.begin());
    launch_interpolated_val_cases(static_cast<int>(cases.size()), device_cases.data(), device_results.data());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_results.begin(), device_results.end(), host_results.begin());

    for (int sample_index = 0; sample_index < static_cast<int>(cases.size()); ++sample_index) {
        const InterpolatedValCase& test_case = cases[sample_index];
        const amrex::Real expected = reference_reconstruction(test_case.stencil, test_case.adv_type, test_case.upwind);
        const amrex::Real actual = host_results[sample_index];

        EXPECT_LE(normalized_error(actual, expected, kDeviceRelTol), amrex::Real(1.0))
            << sample_trace(test_case.adv_type, Coord::x, test_case.upwind,
                            test_case.sample_index, actual, expected, kDeviceRelTol);
    }
}

// Motivation: InterpolateInX/Y/Z are thin directional wrappers around the same
// finite-volume reconstruction contract. This checks each wrapper on device
// against the independent one-dimensional reference assembled from the same
// cell values.
TEST(InterpolationKernel, DirectionalWrappersMatchIndependentReference)
{
    const std::vector<WrapperCase> cases = make_wrapper_cases();
    amrex::FArrayBox qty(stencil_box(), 1);
    amrex::FArrayBox output(eval_box(), static_cast<int>(cases.size()));
    fill_directional_qty(qty);

    amrex::Gpu::DeviceVector<WrapperCase> device_cases(cases.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, cases.begin(), cases.end(), device_cases.begin());
    launch_wrapper_cases(eval_box(), qty.const_array(), device_cases.data(),
                         static_cast<int>(cases.size()), output.array());

    amrex::FArrayBox host_qty(qty.box(), qty.nComp(), amrex::The_Pinned_Arena());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     qty.dataPtr(0), qty.dataPtr(0) + qty.size(), host_qty.dataPtr(0));
    amrex::FArrayBox host_output(output.box(), output.nComp(), amrex::The_Pinned_Arena());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     output.dataPtr(0), output.dataPtr(0) + output.size(), host_output.dataPtr(0));

    const auto qty_arr = host_qty.const_array();
    const auto output_arr = host_output.const_array();
    int sample_index = 0;

    for (int n = 0; n < static_cast<int>(cases.size()); ++n) {
        for (int k = eval_box().smallEnd(2); k <= eval_box().bigEnd(2); ++k) {
            for (int j = eval_box().smallEnd(1); j <= eval_box().bigEnd(1); ++j) {
                for (int i = eval_box().smallEnd(0); i <= eval_box().bigEnd(0); ++i) {
                    const WrapperCase& test_case = cases[n];
                    const StencilValues stencil = load_direction_stencil(qty_arr, i, j, k, test_case.direction);
                    const amrex::Real expected = reference_reconstruction(stencil, test_case.adv_type, test_case.upwind);
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

// Motivation: A zero upwind sign removes the odd upwind correction, so the
// wrapper reductions must agree with centered fourth- and sixth-order references.
TEST(InterpolationKernel, ZeroUpwindWrappersReduceToCentered)
{
    const std::vector<WrapperCase> cases = make_zero_upwind_wrapper_cases();
    amrex::FArrayBox qty(stencil_box(), 1);
    amrex::FArrayBox output(eval_box(), static_cast<int>(cases.size()));
    fill_directional_qty(qty);

    amrex::Gpu::DeviceVector<WrapperCase> device_cases(cases.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, cases.begin(), cases.end(), device_cases.begin());
    launch_wrapper_cases(eval_box(), qty.const_array(), device_cases.data(),
                         static_cast<int>(cases.size()), output.array());

    amrex::FArrayBox host_qty(qty.box(), qty.nComp(), amrex::The_Pinned_Arena());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     qty.dataPtr(0), qty.dataPtr(0) + qty.size(), host_qty.dataPtr(0));
    amrex::FArrayBox host_output(output.box(), output.nComp(), amrex::The_Pinned_Arena());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     output.dataPtr(0), output.dataPtr(0) + output.size(), host_output.dataPtr(0));

    const auto qty_arr = host_qty.const_array();
    const auto output_arr = host_output.const_array();
    int sample_index = 0;

    for (int n = 0; n < static_cast<int>(cases.size()); ++n) {
        for (int k = eval_box().smallEnd(2); k <= eval_box().bigEnd(2); ++k) {
            for (int j = eval_box().smallEnd(1); j <= eval_box().bigEnd(1); ++j) {
                for (int i = eval_box().smallEnd(0); i <= eval_box().bigEnd(0); ++i) {
                    const WrapperCase& test_case = cases[n];
                    const StencilValues stencil = load_direction_stencil(qty_arr, i, j, k, test_case.direction);
                    const amrex::Real expected =
                        reference_reconstruction(stencil, centered_reduction(test_case.adv_type), 0);
                    const amrex::Real actual = output_arr(i, j, k, n);

                    EXPECT_LE(normalized_error(actual, expected, kDeviceRelTol), amrex::Real(1.0))
                        << sample_trace(test_case.adv_type, test_case.direction, 0,
                                        sample_index, actual, expected, kDeviceRelTol);
                    ++sample_index;
                }
            }
        }
    }
}

// Motivation: InterpolatePertFromCell reconstructs qty - r0_arr. With a
// nonlinear background and qty = r0_arr + C, every supported scheme and every
// direction must return exactly C. The pre-fix x/y upwind bug violated this,
// including the x/upwind-3/+1 reproducer at i=j=k=0.
TEST(InterpolationKernel, PerturbationConstantInvariantAllDirections)
{
    const std::vector<WrapperCase> cases = make_wrapper_cases();
    amrex::FArrayBox qty_x(stencil_box(), 1);
    amrex::FArrayBox r0_x(stencil_box(), 1);
    amrex::FArrayBox qty_y(stencil_box(), 1);
    amrex::FArrayBox r0_y(stencil_box(), 1);
    amrex::FArrayBox qty_z(stencil_box(), 1);
    amrex::FArrayBox r0_z(stencil_box(), 1);
    amrex::FArrayBox output(eval_box(), static_cast<int>(cases.size()));

    fill_perturbation_boxes(qty_x, r0_x, Coord::x);
    fill_perturbation_boxes(qty_y, r0_y, Coord::y);
    fill_perturbation_boxes(qty_z, r0_z, Coord::z);

    amrex::Gpu::DeviceVector<WrapperCase> device_cases(cases.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, cases.begin(), cases.end(), device_cases.begin());
    launch_perturbation_cases(eval_box(), device_cases.data(), static_cast<int>(cases.size()),
                              qty_x.const_array(), r0_x.const_array(),
                              qty_y.const_array(), r0_y.const_array(),
                              qty_z.const_array(), r0_z.const_array(),
                              output.array());

    amrex::FArrayBox host_output(output.box(), output.nComp(), amrex::The_Pinned_Arena());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     output.dataPtr(0), output.dataPtr(0) + output.size(), host_output.dataPtr(0));

    const auto output_arr = host_output.const_array();
    int sample_index = 0;

    for (int n = 0; n < static_cast<int>(cases.size()); ++n) {
        for (int k = eval_box().smallEnd(2); k <= eval_box().bigEnd(2); ++k) {
            for (int j = eval_box().smallEnd(1); j <= eval_box().bigEnd(1); ++j) {
                for (int i = eval_box().smallEnd(0); i <= eval_box().bigEnd(0); ++i) {
                    const WrapperCase& test_case = cases[n];
                    const amrex::Real actual = output_arr(i, j, k, n);
                    const amrex::Real expected = kPerturbationConstant;

                    EXPECT_LE(normalized_error(actual, expected, kIdentityRelTol), amrex::Real(1.0))
                        << sample_trace(test_case.adv_type, test_case.direction, test_case.upwind,
                                        sample_index, actual, expected, kIdentityRelTol);
                    ++sample_index;
                }
            }
        }
    }
}

// Motivation: The density wrapper is a named specialization of the generic
// perturbation helper at Rho_comp. This keeps that identity explicit on the
// AMReX device path.
TEST(InterpolationKernel, DensityWrapperMatchesGenericPerturbationHelper)
{
    const std::vector<WrapperCase> cases = make_wrapper_cases();
    amrex::FArrayBox cons_in(stencil_box(), 1);
    amrex::FArrayBox r0(stencil_box(), 1);
    amrex::FArrayBox generic_output(eval_box(), static_cast<int>(cases.size()));
    amrex::FArrayBox density_output(eval_box(), static_cast<int>(cases.size()));
    fill_density_boxes(cons_in, r0);

    amrex::Gpu::DeviceVector<WrapperCase> device_cases(cases.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, cases.begin(), cases.end(), device_cases.begin());
    launch_density_identity_cases(eval_box(), device_cases.data(), static_cast<int>(cases.size()),
                                  cons_in.const_array(), r0.const_array(),
                                  generic_output.array(), density_output.array());

    amrex::FArrayBox host_generic(generic_output.box(), generic_output.nComp(), amrex::The_Pinned_Arena());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     generic_output.dataPtr(0), generic_output.dataPtr(0) + generic_output.size(),
                     host_generic.dataPtr(0));
    amrex::FArrayBox host_density(density_output.box(), density_output.nComp(), amrex::The_Pinned_Arena());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     density_output.dataPtr(0), density_output.dataPtr(0) + density_output.size(),
                     host_density.dataPtr(0));

    const auto generic_arr = host_generic.const_array();
    const auto density_arr = host_density.const_array();
    int sample_index = 0;

    for (int n = 0; n < static_cast<int>(cases.size()); ++n) {
        for (int k = eval_box().smallEnd(2); k <= eval_box().bigEnd(2); ++k) {
            for (int j = eval_box().smallEnd(1); j <= eval_box().bigEnd(1); ++j) {
                for (int i = eval_box().smallEnd(0); i <= eval_box().bigEnd(0); ++i) {
                    const WrapperCase& test_case = cases[n];
                    const amrex::Real expected = generic_arr(i, j, k, n);
                    const amrex::Real actual = density_arr(i, j, k, n);

                    EXPECT_LE(normalized_error(actual, expected, kIdentityRelTol), amrex::Real(1.0))
                        << sample_trace(test_case.adv_type, test_case.direction, test_case.upwind,
                                        sample_index, actual, expected, kIdentityRelTol);
                    ++sample_index;
                }
            }
        }
    }
}

// Motivation: 2nd/3rd/4th-order branches must ignore the +2/-3 stencil cells.
// Changing only those far cells should leave the face value unchanged.
TEST(InterpolationKernel, LowerOrderSchemesIgnoreFarStencilSentinels)
{
    const std::vector<SentinelCase> cases = make_lower_order_sentinel_cases();
    amrex::FArrayBox base_qty(stencil_box(), 1);
    amrex::FArrayBox modified_qty(stencil_box(), 1);
    fill_directional_qty(base_qty);
    fill_directional_qty(modified_qty);
    apply_far_sentinels(modified_qty);

    amrex::Gpu::DeviceVector<SentinelCase> device_cases(cases.size());
    amrex::Gpu::DeviceVector<amrex::Real> device_base(cases.size(), amrex::Real(0.0));
    amrex::Gpu::DeviceVector<amrex::Real> device_modified(cases.size(), amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real> host_base(cases.size(), amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real> host_modified(cases.size(), amrex::Real(0.0));

    amrex::Gpu::copy(amrex::Gpu::hostToDevice, cases.begin(), cases.end(), device_cases.begin());
    launch_sentinel_cases(static_cast<int>(cases.size()), device_cases.data(),
                          base_qty.const_array(), modified_qty.const_array(),
                          device_base.data(), device_modified.data());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_base.begin(), device_base.end(), host_base.begin());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_modified.begin(), device_modified.end(), host_modified.begin());

    amrex::FArrayBox host_base_qty(base_qty.box(), base_qty.nComp(), amrex::The_Pinned_Arena());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     base_qty.dataPtr(0), base_qty.dataPtr(0) + base_qty.size(),
                     host_base_qty.dataPtr(0));
    amrex::FArrayBox host_modified_qty(modified_qty.box(), modified_qty.nComp(), amrex::The_Pinned_Arena());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     modified_qty.dataPtr(0), modified_qty.dataPtr(0) + modified_qty.size(),
                     host_modified_qty.dataPtr(0));

    const auto base_arr = host_base_qty.const_array();
    const auto modified_arr = host_modified_qty.const_array();

    for (int sample_index = 0; sample_index < static_cast<int>(cases.size()); ++sample_index) {
        const SentinelCase& test_case = cases[sample_index];
        const StencilValues base_stencil = load_direction_stencil(base_arr, 0, 0, 0, test_case.direction);
        const StencilValues modified_stencil = load_direction_stencil(modified_arr, 0, 0, 0, test_case.direction);
        const amrex::Real expected_base = reference_reconstruction(base_stencil, test_case.adv_type, test_case.upwind);
        const amrex::Real expected_modified = reference_reconstruction(modified_stencil, test_case.adv_type, test_case.upwind);

        EXPECT_LE(normalized_error(host_base[sample_index], expected_base, kDeviceRelTol), amrex::Real(1.0))
            << sample_trace(test_case.adv_type, test_case.direction, test_case.upwind,
                            sample_index, host_base[sample_index], expected_base, kDeviceRelTol);
        EXPECT_LE(normalized_error(host_modified[sample_index], expected_modified, kDeviceRelTol), amrex::Real(1.0))
            << sample_trace(test_case.adv_type, test_case.direction, test_case.upwind,
                            sample_index, host_modified[sample_index], expected_modified, kDeviceRelTol);
        EXPECT_LE(normalized_error(host_modified[sample_index], host_base[sample_index], kIdentityRelTol), amrex::Real(1.0))
            << sample_trace(test_case.adv_type, test_case.direction, test_case.upwind,
                            sample_index, host_modified[sample_index], host_base[sample_index], kIdentityRelTol);
    }
}

// Motivation: 5th/6th-order branches must use the +2/-3 cells. Changing only
// those far stencil entries should measurably change the face value.
TEST(InterpolationKernel, FifthAndSixthOrderSchemesUseFarStencil)
{
    const std::vector<SentinelCase> cases = make_high_order_sentinel_cases();
    amrex::FArrayBox base_qty(stencil_box(), 1);
    amrex::FArrayBox modified_qty(stencil_box(), 1);
    fill_directional_qty(base_qty);
    fill_directional_qty(modified_qty);
    apply_far_sentinels(modified_qty);

    amrex::Gpu::DeviceVector<SentinelCase> device_cases(cases.size());
    amrex::Gpu::DeviceVector<amrex::Real> device_base(cases.size(), amrex::Real(0.0));
    amrex::Gpu::DeviceVector<amrex::Real> device_modified(cases.size(), amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real> host_base(cases.size(), amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real> host_modified(cases.size(), amrex::Real(0.0));

    amrex::Gpu::copy(amrex::Gpu::hostToDevice, cases.begin(), cases.end(), device_cases.begin());
    launch_sentinel_cases(static_cast<int>(cases.size()), device_cases.data(),
                          base_qty.const_array(), modified_qty.const_array(),
                          device_base.data(), device_modified.data());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_base.begin(), device_base.end(), host_base.begin());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_modified.begin(), device_modified.end(), host_modified.begin());

    amrex::FArrayBox host_base_qty(base_qty.box(), base_qty.nComp(), amrex::The_Pinned_Arena());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     base_qty.dataPtr(0), base_qty.dataPtr(0) + base_qty.size(),
                     host_base_qty.dataPtr(0));
    amrex::FArrayBox host_modified_qty(modified_qty.box(), modified_qty.nComp(), amrex::The_Pinned_Arena());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     modified_qty.dataPtr(0), modified_qty.dataPtr(0) + modified_qty.size(),
                     host_modified_qty.dataPtr(0));

    const auto base_arr = host_base_qty.const_array();
    const auto modified_arr = host_modified_qty.const_array();

    for (int sample_index = 0; sample_index < static_cast<int>(cases.size()); ++sample_index) {
        const SentinelCase& test_case = cases[sample_index];
        const StencilValues base_stencil = load_direction_stencil(base_arr, 0, 0, 0, test_case.direction);
        const StencilValues modified_stencil = load_direction_stencil(modified_arr, 0, 0, 0, test_case.direction);
        const amrex::Real expected_base = reference_reconstruction(base_stencil, test_case.adv_type, test_case.upwind);
        const amrex::Real expected_modified = reference_reconstruction(modified_stencil, test_case.adv_type, test_case.upwind);
        const amrex::Real expected_delta = std::abs(expected_modified - expected_base);
        const amrex::Real actual_delta = std::abs(host_modified[sample_index] - host_base[sample_index]);

        EXPECT_LE(normalized_error(host_base[sample_index], expected_base, kDeviceRelTol), amrex::Real(1.0))
            << sample_trace(test_case.adv_type, test_case.direction, test_case.upwind,
                            sample_index, host_base[sample_index], expected_base, kDeviceRelTol);
        EXPECT_LE(normalized_error(host_modified[sample_index], expected_modified, kDeviceRelTol), amrex::Real(1.0))
            << sample_trace(test_case.adv_type, test_case.direction, test_case.upwind,
                            sample_index, host_modified[sample_index], expected_modified, kDeviceRelTol);
        EXPECT_GT(expected_delta, scaled_tol(expected_modified, expected_base, kIdentityRelTol))
            << sample_trace(test_case.adv_type, test_case.direction, test_case.upwind,
                            sample_index, expected_modified, expected_base, kIdentityRelTol);
        EXPECT_GT(actual_delta, scaled_tol(host_modified[sample_index], host_base[sample_index], kIdentityRelTol))
            << sample_trace(test_case.adv_type, test_case.direction, test_case.upwind,
                            sample_index, host_modified[sample_index], host_base[sample_index], kIdentityRelTol);
    }
}
