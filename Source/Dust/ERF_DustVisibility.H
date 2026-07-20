/**
 * @file ERF_DustVisibility.H
 * @brief Koschmieder horizontal visibility from near-surface PM10.
 *
 * V [m] = 3.912 / (k_ext * PM10 [kg/m³])
 * Capped at visibility_max_m (100,000 m) when PM10 < epsilon.
 *
 * References:
 *   Koschmieder, H. (1924). Beitr. Phys. Atmos., 12, 33-55.
 *   Seinfeld & Pandis (2006), Ch. 15.
 *   OSHA 29 CFR 1926.601: haul road closure < 300 m.
 */

#ifdef ERF_USE_DUST

#ifndef ERF_DUST_VISIBILITY_H
#define ERF_DUST_VISIBILITY_H

#include <AMReX_MultiFab.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <fstream>
#include <string>

namespace DustVisibilityConst {
    constexpr amrex::Real KOSCHMIEDER_COEFF = 3.912;  ///< Koschmieder constant [m]
    constexpr amrex::Real EPSILON            = 1.0e-12; ///< Cutoff for near-zero PM10 [kg/m³]
}

/**
 * @brief Compute Koschmieder visibility from near-surface PM10.
 *
 * Visibility [m] = 3.912 / (k_ext * PM10 [kg/m³])
 * Floor: 10 m when PM10 ~ 0.
 * Cap: visibility_max_m (typically 100,000 m).
 *
 * @param[out] vis_mf         Visibility MultiFab [m]
 * @param[in]  pm10_mf        PM10 concentration [µg/m³] on dust grid
 * @param[in]  k_ext          Mass extinction coefficient [m²/kg]
 * @param[in]  vis_max_m      Maximum visibility cap [m]
 */
inline void compute_visibility(
    amrex::MultiFab&       vis_mf,
    const amrex::MultiFab& pm10_mf,
    amrex::Real            k_ext,
    amrex::Real            vis_max_m)
{
    using namespace amrex;
    for (MFIter mfi(vis_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto vis = vis_mf.array(mfi);
        auto pm10 = pm10_mf.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Convert PM10 from µg/m³ to kg/m³: divide by 1e9
            Real pm10_kg_m3 = pm10(i,j,k) * 1.0e-9;

            Real visibility;
            if (pm10_kg_m3 < DustVisibilityConst::EPSILON) {
                // Near-zero PM10: visibility is maximum
                visibility = vis_max_m;
            } else {
                // Koschmieder: V = 3.912 / (k_ext * PM10)
                Real b_ext = k_ext * pm10_kg_m3;
                visibility = DustVisibilityConst::KOSCHMIEDER_COEFF / b_ext;
                // Floor at 10 m
                visibility = amrex::max(visibility, 10.0);
                // Cap at maximum
                visibility = amrex::min(visibility, vis_max_m);
            }
            vis(i,j,k) = visibility;
        });
    }
}

/**
 * @brief Compute haul road closure and reduced visibility warning flags.
 *
 * Sets:
 *   closure_flag(i,j) = 1 if visibility(i,j) < closure_thresh, else 0
 *   warning_flag(i,j) = 1 if visibility(i,j) < warning_thresh, else 0
 *
 * @param[out] closure_flag    Haul road closure flag [0/1]
 * @param[out] warning_flag    Reduced visibility warning flag [0/1]
 * @param[in]  vis_mf          Visibility MultiFab [m]
 * @param[in]  closure_thresh  Closure threshold [m]
 * @param[in]  warning_thresh  Warning threshold [m]
 */
inline void compute_visibility_flags(
    amrex::MultiFab&       closure_flag,
    amrex::MultiFab&       warning_flag,
    const amrex::MultiFab& vis_mf,
    amrex::Real            closure_thresh,
    amrex::Real            warning_thresh)
{
    using namespace amrex;
    for (MFIter mfi(closure_flag, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto cf = closure_flag.array(mfi);
        auto wf = warning_flag.array(mfi);
        auto vis = vis_mf.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            cf(i,j,k) = (vis(i,j,k) < closure_thresh) ? 1.0 : 0.0;
            wf(i,j,k) = (vis(i,j,k) < warning_thresh) ? 1.0 : 0.0;
        });
    }
}

/**
 * @brief Append visibility diagnostics to CSV file.
 *
 * Header format:
 *   step,time_s,visibility_min_m,visibility_mean_m,haul_road_closure_cells,reduced_visibility_cells
 *
 * Only IOProcessor writes; all ranks compute reductions.
 *
 * @param[in] step      Current step number
 * @param[in] time_s    Current time [s]
 * @param[in] vis_mf    Visibility MultiFab [m]
 * @param[in] closure   Closure flag MultiFab [0/1]
 * @param[in] warning   Warning flag MultiFab [0/1]
 * @param[in] filename  Output CSV filename
 */
inline void write_visibility_stats(int step, amrex::Real time_s,
                                   const amrex::MultiFab& vis_mf,
                                   const amrex::MultiFab& closure,
                                   const amrex::MultiFab& warning,
                                   const std::string& filename)
{
    using namespace amrex;

    // Compute domain reductions on all ranks
    Real vis_min = vis_mf.min(0);
    Real vis_mean = vis_mf.sum(0) / vis_mf.boxArray().numPts();
    Real n_closure = closure.sum(0);
    Real n_warning = warning.sum(0);

    if (!ParallelDescriptor::IOProcessor()) return;

    bool exists = std::ifstream(filename).good();
    std::ofstream out(filename, std::ios::app);
    if (!out.is_open()) {
        Print() << "[DUST VISIBILITY] WARNING: Cannot open file: " << filename << "\n";
        return;
    }

    if (!exists) {
        out << "# ERF-Dust visibility diagnostics (Koschmieder)\n"
            << "step,time_s,visibility_min_m,visibility_mean_m,haul_road_closure_cells,reduced_visibility_cells\n";
    }

    out << step << "," << time_s << "," << vis_min << "," << vis_mean << ","
        << (long)n_closure << "," << (long)n_warning << "\n";
}

#endif // ERF_DUST_VISIBILITY_H
#endif // ERF_USE_DUST
