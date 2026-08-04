#include <array>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <gtest/gtest.h>

#include "ERF_ReadFromERFBdy.H"

namespace
{

class ScopedTestDirectory
{
public:
    explicit ScopedTestDirectory (std::filesystem::path path)
        : m_path(std::move(path))
    {
        std::filesystem::remove_all(m_path);
        std::filesystem::create_directories(m_path / "Time_000000");
    }

    ~ScopedTestDirectory () { std::filesystem::remove_all(m_path); }

    const std::filesystem::path& path () const { return m_path; }

private:
    std::filesystem::path m_path;
};

void write_boundary_fab (const std::filesystem::path& filename, amrex::Real value)
{
    const amrex::Box bx(amrex::IntVect(-1, 0, 0), amrex::IntVect(0, 1, 1));
    amrex::FArrayBox fab(bx, 1, amrex::The_Pinned_Arena());
    fab.template setVal<amrex::RunOn::Host>(value);
    std::ofstream ofs(filename, std::ios::out | std::ios::binary);
    fab.writeOn(ofs);
}

} // namespace

TEST(ReadERFBdy, KeepsPersistentBoundaryDataDeviceAccessible)
{
    ScopedTestDirectory test_dir(
        std::filesystem::current_path() / "erf-erfbdy-read-device-test");
    const auto time_dir = test_dir.path() / "Time_000000";

    write_boundary_fab(time_dir / "BdyData_xlo_var0", amrex::Real(1.0));
    write_boundary_fab(time_dir / "BdyData_xhi_var0", amrex::Real(2.0));
    write_boundary_fab(time_dir / "BdyData_ylo_var0", amrex::Real(3.0));
    write_boundary_fab(time_dir / "BdyData_yhi_var0", amrex::Real(4.0));

    amrex::Vector<amrex::Vector<amrex::FArrayBox>> xlo(1), xhi(1), ylo(1), yhi(1);
    read_from_erfbdy(0, test_dir.path().string(), xlo, xhi, ylo, yhi, 1, 1);

    EXPECT_EQ(xlo[0][0].arena(), amrex::The_Arena());
    EXPECT_EQ(xhi[0][0].arena(), amrex::The_Arena());
    EXPECT_EQ(ylo[0][0].arena(), amrex::The_Arena());
    EXPECT_EQ(yhi[0][0].arena(), amrex::The_Arena());

    amrex::Gpu::DeviceVector<amrex::Real> observed_device(4);
    const auto xlo_arr = xlo[0][0].const_array();
    const auto xhi_arr = xhi[0][0].const_array();
    const auto ylo_arr = ylo[0][0].const_array();
    const auto yhi_arr = yhi[0][0].const_array();
    amrex::Real* observed = observed_device.data();
    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int) noexcept {
        observed[0] = xlo_arr(-1, 0, 0);
        observed[1] = xhi_arr(-1, 0, 0);
        observed[2] = ylo_arr(-1, 0, 0);
        observed[3] = yhi_arr(-1, 0, 0);
    });

    std::array<amrex::Real, 4> observed_host{};
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     observed_device.begin(), observed_device.end(),
                     observed_host.begin());
    EXPECT_EQ(observed_host[0], amrex::Real(1.0));
    EXPECT_EQ(observed_host[1], amrex::Real(2.0));
    EXPECT_EQ(observed_host[2], amrex::Real(3.0));
    EXPECT_EQ(observed_host[3], amrex::Real(4.0));
}
