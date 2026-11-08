#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <array>
#include <cmath>
#include <iostream>
#include <string>

int main (int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "usage: erf_immersed_cutoff_check plotfile active|coarse|active_native|coarse_native|coarse_ibseb\n";
        return 2;
    }

    amrex::Initialize(argc, argv, false);
    const std::string mode(argv[2]);
    amrex::PlotFileData plotfile(argv[1]);
    amrex::MultiFab mask = plotfile.get(0, "terrain_IB_mask");
    amrex::MultiFab xvel = plotfile.get(0, "x_velocity");

    const auto mask_min = mask.min(0);
    const auto mask_max = mask.max(0);
    const auto u_min = xvel.min(0);
    const auto u_max = xvel.max(0);
    bool passed = false;
    const bool active = mode == "active" || mode == "active_native";
    const bool coarse = mode == "coarse" || mode == "coarse_native" || mode == "coarse_ibseb";
    const bool native = mode == "active_native" || mode == "coarse_native" || mode == "coarse_ibseb";
    if (active) {
        passed = mask_min < 0.01 && mask_max > 0.99;
    } else if (coarse) {
        passed = mask_min == 0.0 && mask_max == 0.0 &&
                 std::abs(u_min - 3.0) < 1.e-6 && std::abs(u_max - 3.0) < 1.e-6;
    }
    if (!passed) {
        std::cerr << "Immersed building cutoff (" << mode << "): mask=["
                  << mask_min << ", " << mask_max << "], x_velocity=["
                  << u_min << ", " << u_max << "]\n";
    }
    if (native) {
        constexpr std::array<const char*, 3> directions = {"x", "y", "z"};
        for (const char* dir : directions) {
            const std::string name = std::string("terrain_blank_") + dir + "face";
            amrex::PlotFileData face_plotfile(std::string(argv[1]) + "_terrain_blank_" + dir + "face");
            amrex::MultiFab face_mask = face_plotfile.get(0, name);
            const auto face_min = face_mask.min(0);
            const auto face_max = face_mask.max(0);
            const bool face_ok = active ? (face_min < 0.01 && face_max > 0.99)
                                        : (face_min == 0.0 && face_max == 0.0);
            if (!face_ok) {
                passed = false;
                std::cerr << name << "=[" << face_min << ", " << face_max << "]\n";
            }
        }
    }
    amrex::Finalize();
    return passed ? 0 : 1;
}
