#include <algorithm>
#include <string>

#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>

#include <gtest/gtest.h>

#include "ERF_GTestEOSCommon.H"

using namespace eos_test;

namespace {

void launch_inverse_property_errors (amrex::Gpu::DeviceVector<amrex::Real>& errors)
{
    // Keep this device-path metric calculation explicit. The scalar tests
    // check the same physical invariants on the host, while this kernel test
    // independently proves that the EOS utilities are callable inside AMReX
    // ParallelFor and can report normalized errors back to host-side GTest checks.
    auto* error_ptr = errors.data();

    amrex::ParallelFor(kNumEOSStates, [=] AMREX_GPU_DEVICE (int index) noexcept {
        const EOSState state = get_test_state(index);
        const amrex::Real theta = getThgivenTandP(state.temperature, state.pressure, kRdOcp);
        const amrex::Real rho = getRhogivenTandPress(state.temperature, state.pressure, state.qv);
        const amrex::Real rhotheta = rho * theta;
        const amrex::Real rho_from_theta = getRhogivenThetaPress(theta, state.pressure, kRdOcp, state.qv);
        const amrex::Real pressure = getPgivenRTh(rhotheta, state.qv);
        const amrex::Real exner_from_pressure = getExnergivenP(state.pressure, kRdOcp);
        const amrex::Real exner_from_rhotheta =
            getExnergivenRTh(getRhoThetagivenP(state.pressure, state.qv), kRdOcp, state.qv);
        const amrex::Real dpdrho = getdPdRgivenConstantTheta(rho, theta, state.qv);
        const amrex::Real expected_pressure = rho * R_d * state.temperature * moisture_factor(state.qv);
        const amrex::Real ps = amrex::Real(2350.0) + amrex::Real(50.0) * index;

        const amrex::Real vapor_error = amrex::max(
            normalized_error(compute_vapor_pressure(ps, amrex::Real(0.0)), amrex::Real(0.0), kEOSRelTol),
            amrex::max(
                normalized_error(compute_vapor_pressure(ps, amrex::Real(1.0)), ps, kEOSRelTol),
                normalized_error(compute_vapor_pressure(ps, amrex::Real(0.75)), amrex::Real(0.75) * ps,
                                 kEOSRelTol)));

        const int offset = index * kKernelMetricsPerState;
        error_ptr[offset + 0] =
            normalized_error(getTgivenPandTh(state.pressure, theta, kRdOcp), state.temperature, kEOSRelTol);
        error_ptr[offset + 1] = amrex::max(
            normalized_error(pressure, state.pressure, kEOSRelTol),
            normalized_error(getRhoThetagivenP(state.pressure, state.qv), rhotheta, kEOSRelTol));
        error_ptr[offset + 2] = amrex::max(
            normalized_error(rho_from_theta, rho, kEOSRelTol),
            normalized_error(expected_pressure, state.pressure, kEOSRelTol));
        error_ptr[offset + 3] = normalized_error(exner_from_rhotheta, exner_from_pressure, kEOSRelTol);
        error_ptr[offset + 4] = normalized_error(dpdrho, Gamma * pressure / rho, kEOSRelTol);
        error_ptr[offset + 5] = vapor_error;
    });
}

std::string metric_label (const int metric)
{
    switch (metric) {
    case 0:
        return "T/theta inverse";
    case 1:
        return "pressure/rhotheta inverse";
    case 2:
        return "density closure";
    case 3:
        return "Exner closure";
    case 4:
        return "dPdrho closure";
    default:
        return "fractional RH vapor pressure";
    }
}

} // namespace

// Motivation: This test proves the EOS and its helpers are callable through the
// AMReX device path while keeping all GTest assertions on the host.
TEST(ERFEOSKernel, InversePropertiesMatchHostPortable)
{
    amrex::Gpu::DeviceVector<amrex::Real> device_errors(kNumEOSStates * kKernelMetricsPerState,
                                                        amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real> host_errors(device_errors.size(), amrex::Real(0.0));

    launch_inverse_property_errors(device_errors);
    amrex::Gpu::streamSynchronize();
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_errors.begin(), device_errors.end(), host_errors.begin());

    const amrex::Real max_error = *std::max_element(host_errors.begin(), host_errors.end());
    EXPECT_LE(max_error, kKernelErrorLimit);

    for (int index = 0; index < kNumEOSStates; ++index) {
        const EOSState state = get_test_state(index);

        for (int metric = 0; metric < kKernelMetricsPerState; ++metric) {
            const amrex::Real error = host_errors[index * kKernelMetricsPerState + metric];
            EXPECT_LE(error, kKernelErrorLimit)
                << state_label(index, state) << " metric=" << metric_label(metric);
        }
    }
}