#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>

#include <gtest/gtest.h>

#include "ERF_GTestCoalescenceCommon.H"

using namespace coalescence_test;

namespace {

// Evaluate every kernel-property violation inside an AMReX ParallelFor so the
// collision kernels run on the device when ERF is built for GPU, then report the
// violations to host-side GTest checks. coalescence_metric is
// AMREX_GPU_HOST_DEVICE, so the same code runs on host and device.
void launch_coalescence_metrics (amrex::Gpu::DeviceVector<amrex::Real>& metrics)
{
    auto* metric_ptr = metrics.data();
    amrex::ParallelFor(kNumMetrics, [=] AMREX_GPU_DEVICE (int i) noexcept {
        metric_ptr[i] = coalescence_metric(i);
    });
}

} // namespace

// Motivation: prove the collision and riming-efficiency kernels are callable
// through the AMReX device path and satisfy the same physical properties as the
// host scalar test, with all assertions on the host.
TEST(CoalescenceKernel, KernelPropertiesHoldOnDevice)
{
    amrex::Gpu::DeviceVector<amrex::Real> device_metrics(kNumMetrics, amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real>   host_metrics(kNumMetrics, amrex::Real(0.0));

    launch_coalescence_metrics(device_metrics);
    amrex::Gpu::streamSynchronize();
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_metrics.begin(), device_metrics.end(), host_metrics.begin());

    for (int i = 0; i < kNumMetrics; ++i) {
        EXPECT_LE(host_metrics[i], kTol) << metric_label(i);
    }
}
