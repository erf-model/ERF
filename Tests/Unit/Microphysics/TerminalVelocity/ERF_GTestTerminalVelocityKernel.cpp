#include <algorithm>
#include <string>

#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>

#include <gtest/gtest.h>

#include "ERF_GTestTerminalVelocityCommon.H"

using namespace terminal_velocity_test;

namespace {

// Evaluate every fall-speed property inside an AMReX ParallelFor so the
// terminal-velocity routines run on the device when ERF is built for GPU, then
// report the violations back to host-side GTest checks. The functor methods are
// AMREX_GPU_HOST_DEVICE, so the same code runs on host and device.
void launch_terminal_velocity_metrics (amrex::Gpu::DeviceVector<amrex::Real>& metrics)
{
    auto* metric_ptr = metrics.data();
    amrex::ParallelFor(kNumMetrics, [=] AMREX_GPU_DEVICE (int i) noexcept {
        metric_ptr[i] = terminal_velocity_metric(i);
    });
}

} // namespace

// Motivation: This proves the terminal-velocity routines are callable through
// the AMReX device path and satisfy the same physical properties as the host
// scalar test, while keeping all GTest assertions on the host.
TEST(TerminalVelocityKernel, FallSpeedPropertiesHoldOnDevice)
{
    amrex::Gpu::DeviceVector<amrex::Real> device_metrics(kNumMetrics, amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real>   host_metrics(kNumMetrics, amrex::Real(0.0));

    launch_terminal_velocity_metrics(device_metrics);
    amrex::Gpu::streamSynchronize();
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_metrics.begin(), device_metrics.end(), host_metrics.begin());

    for (int i = 0; i < kNumMetrics; ++i) {
        EXPECT_LE(host_metrics[i], kTol) << metric_label(i);
    }
}
