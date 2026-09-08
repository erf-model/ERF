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

// Evaluate the ice mass rate at a supersaturated and a subsaturated state on the
// device, so the deposition/sublimation branch runs through the AMReX device path.
void launch_subl_rates (amrex::Gpu::DeviceVector<amrex::Real>& rate,
                        const amrex::Gpu::DeviceVector<amrex::Real>& S)
{
    auto* rate_ptr = rate.data();
    const auto* S_ptr = S.data();
    const int n = static_cast<int>(S.size());
    amrex::ParallelFor(n, [=] AMREX_GPU_DEVICE (int i) noexcept {
        rate_ptr[i] = deposition_mass_rate(amrex::Real(273.15 - 15.0), S_ptr[i]);
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

// Motivation: prove the deposition/sublimation mass rate (dMdt::rhs_func) runs on
// the device and keeps the right sign -- growth when supersaturated over ice,
// mass loss (sublimation) when subsaturated.
TEST(IceDepositionKernel, SublimationSignOnDevice)
{
    amrex::Gpu::DeviceVector<amrex::Real> S(2);
    amrex::Gpu::HostVector<amrex::Real>   S_host = {amrex::Real(1.05), amrex::Real(0.95)};
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, S_host.begin(), S_host.end(), S.begin());

    amrex::Gpu::DeviceVector<amrex::Real> rate(2, amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real>   rate_host(2, amrex::Real(0.0));

    launch_subl_rates(rate, S);
    amrex::Gpu::streamSynchronize();
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, rate.begin(), rate.end(), rate_host.begin());

    EXPECT_GT(rate_host[0], amrex::Real(0.0)) << "deposition rate at S=1.05";
    EXPECT_LT(rate_host[1], amrex::Real(0.0)) << "sublimation rate at S=0.95";
}
