#include <algorithm>

#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>

#include <gtest/gtest.h>

#include "ERF_GTestIceDepositionCommon.H"

using namespace ice_deposition_test;

namespace {

// Grow the crystal at a few temperatures inside an AMReX ParallelFor so the
// deposition routines run on the device when ERF is built for GPU. We store the
// aspect ratio phi=c/a for each temperature and assert the habit on the host.
void launch_aspect_ratios (amrex::Gpu::DeviceVector<amrex::Real>& phi,
                           const amrex::Gpu::DeviceVector<amrex::Real>& T_C)
{
    auto* phi_ptr = phi.data();
    const auto* T_ptr = T_C.data();
    const int n = static_cast<int>(T_C.size());
    amrex::ParallelFor(n, [=] AMREX_GPU_DEVICE (int i) noexcept {
        const CrystalState x = grow_crystal(amrex::Real(273.15) + T_ptr[i]);
        phi_ptr[i] = aspect_ratio(x);
    });
}

} // namespace

// Motivation: prove the deposition/habit routines are callable through the AMReX
// device path and reproduce the TH91 habit (plate at -15 C, column at -6 C) on
// the device, with all assertions on the host.
TEST(IceDepositionKernel, HabitOnDeviceMatchesTakahashi)
{
    amrex::Gpu::DeviceVector<amrex::Real> T_C(3);
    amrex::Gpu::HostVector<amrex::Real>   T_host = {amrex::Real(-23), amrex::Real(-15), amrex::Real(-6)};
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, T_host.begin(), T_host.end(), T_C.begin());

    amrex::Gpu::DeviceVector<amrex::Real> phi(3, amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real>   phi_host(3, amrex::Real(0.0));

    launch_aspect_ratios(phi, T_C);
    amrex::Gpu::streamSynchronize();
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, phi.begin(), phi.end(), phi_host.begin());

    EXPECT_LT(phi_host[1], 1.0) << "phi(-15C)=" << phi_host[1];   // plate
    EXPECT_GT(phi_host[2], 1.0) << "phi(-6C)="  << phi_host[2];   // column
    EXPECT_LT(phi_host[1], phi_host[0]);                          // -15 flatter than -23
}
