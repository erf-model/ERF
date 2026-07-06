#ifndef ERF_FIRE_STATS_OUTPUT_H
#define ERF_FIRE_STATS_OUTPUT_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReEx_ParallelDescriptor.H>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

/**
 * @file ERF_FireStatsOutput.H
 * @brief CSV output functions for fire statistics (burned area, perimeter, etc.)
 *
 * Writes fire statistics including burned area [ha], perimeter [km], active cells,
 * head rate of spread [m/s], and fire ellipse dimensions to a CSV file.
 *
 * Based on wildfire_levelset write_xy_data.H implementation.
 */

/**
 * @brief Compute fire ellipse axes from area and perimeter
 *
 * Approximates the fire as an ellipse given its perimeter and burned area.
 * Uses Ramanujan's approximation for ellipse perimeter.
 *
 * @param[in] area_m2       Burned area [m²]
 * @param[in] perimeter_m   Fire perimeter [m]
 *
 * @return {major_axis_m, minor_axis_m} (full axes, not semi-axes)
 */
inline std::pair<amrex::Real, amrex::Real>
compute_fire_ellipse_axes(amrex::Real area_m2, amrex::Real perimeter_m)
{
    if (perimeter_m <= 0.0 || area_m2 <= 0.0) {
        return {0.0, 0.0};
    }

    constexpr amrex::Real PI = M_PI;
    
    // From Area = π * a * b and Perimeter ≈ π * (a + b)
    // we get: a + b ≈ Perimeter / π and a * b = Area / π
    Real sum_ab  = perimeter_m / PI;
    Real prod_ab = area_m2 / PI;
    
    if (sum_ab <= 0.0 || prod_ab <= 0.0) {
        return {0.0, 0.0};
    }
    
    // Solve for semi-major and semi-minor axes using quadratic formula
    Real discriminant = sum_ab * sum_ab - 4.0 * prod_ab;
    if (discriminant < 0.0) {
        // No real solution; return a circle approximation
        Real radius = std::sqrt(area_m2 / PI);
        return {2.0 * radius, 2.0 * radius};
    }
    
    Real sqrt_disc = std::sqrt(discriminant);
    Real semi_a = (sum_ab + sqrt_disc) / 2.0;  // semi-major axis
    Real semi_b = (sum_ab - sqrt_disc) / 2.0;  // semi-minor axis
    
    // Return full axes (2 * semi-axes)
    return {2.0 * semi_a, 2.0 * semi_b};
}

/**
 * @brief Write CSV header for fire statistics
 *
 * Writes the column headers for the fire statistics CSV file.
 * Should be called once at simulation startup before first append.
 *
 * @param[in] filename  Path to CSV file to create
 */
inline void
write_fire_stats_header(const std::string& filename)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) return;
    
    std::ofstream out(filename, std::ios::out | std::ios::trunc);
    if (!out.is_open()) {
        amrex::Print() << "[FIRE] WARNING: Cannot open fire_stats file: " << filename << "\n";
        return;
    }
    
    out << "step,time_s,burned_area_ha,perimeter_km,"
           "active_front_cells,head_ros_ms,major_axis_m,minor_axis_m,heat_flux_max_Wm2\n";
    out.close();
}

/**
 * @brief Append one row to fire statistics CSV file
 *
 * Computes burned area [ha], perimeter length [km], active front cell count,
 * head rate of spread [m/s], and fire ellipse axes, then appends one row to CSV.
 *
 * Burned area and perimeter are based on arrival_time >= 0 (every cell that
 * has ever been ignited, including the burned interior).
 *
 * @param[in] phi               Fire front phi field (current timestep)
 * @param[in] arrival_time      Cumulative arrival time field (burned region tracker)
 * @param[in] geom              Fire grid geometry
 * @param[in] step              Simulation step number
 * @param[in] time_s            Simulation time [s]
 * @param[in] filename          CSV file path
 * @param[in] R_mf              Optional rate-of-spread field for head ROS tracking
 * @param[in] Q_mf              Optional heat flux field for max heat flux tracking
 */
inline void
append_fire_stats(const amrex::MultiFab& phi,
                  const amrex::MultiFab& arrival_time,
                  const amrex::Geometry& geom,
                  int step,
                  amrex::Real time_s,
                  const std::string& filename,
                  const amrex::MultiFab* R_mf = nullptr,
                  const amrex::MultiFab* Q_mf = nullptr)
{
    using namespace amrex;
    
    const auto dx = geom.CellSize();
    Real cell_area = dx[0] * dx[1];

    Long n_burned = 0, n_active = 0;
    Real perim_m = Real(0.0);
    Real head_ros = Real(0.0);

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p  = phi.const_array(mfi);
        auto const at = arrival_time.const_array(mfi);
        const IntVect lo = bx.smallEnd();
        const IntVect hi = bx.bigEnd();

        Array4<const Real> R_arr;
        if (R_mf) R_arr = R_mf->const_array(mfi);

        amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
            bool burned = (at(i, j, k, 0) >= Real(0.0));
            Real pv     = p(i, j, k, 0);
            
            if (burned) {
                ++n_burned;
                
                // Track head ROS for perimeter cells
                bool on_perim = false;
                if (i > lo[0]) on_perim |= (at(i-1,j,k,0) < Real(0.0));
                if (i < hi[0]) on_perim |= (at(i+1,j,k,0) < Real(0.0));
                if (j > lo[1]) on_perim |= (at(i,j-1,k,0) < Real(0.0));
                if (j < hi[1]) on_perim |= (at(i,j+1,k,0) < Real(0.0));
                
                if (on_perim && R_mf) {
                    Real r = R_arr(i, j, k, 0);
                    if (r > head_ros) head_ros = r;
                }
            } else if (pv < Real(0.1) && pv >= Real(0.0)) {
                ++n_active;  // unburned cells within threshold of front
            }
            
            // Edge-crossing perimeter (arrival_time boundary)
            if (i < hi[0]) {
                bool burned_n = (at(i+1, j, k, 0) >= Real(0.0));
                if (burned != burned_n) perim_m += dx[1];
            }
            if (j < hi[1]) {
                bool burned_n = (at(i, j+1, k, 0) >= Real(0.0));
                if (burned != burned_n) perim_m += dx[0];
            }
        });
    }

    ParallelDescriptor::ReduceLongSum(n_burned);
    ParallelDescriptor::ReduceLongSum(n_active);
    ParallelDescriptor::ReduceRealSum(perim_m);
    if (R_mf) {
        ParallelDescriptor::ReduceRealMax(head_ros);
    }

    Real Q_max = Q_mf ? Q_mf->max(0) : 0.0_rt;
    if (Q_mf) {
        ParallelDescriptor::ReduceRealMax(Q_max);
    }

    Real burned_ha    = static_cast<Real>(n_burned) * cell_area / Real(1.0e4);
    Real perim_km     = perim_m / Real(1000.0);
     
    // Calculate fire ellipse major/minor axes from area and perimeter
    Real area_m2 = burned_ha * Real(1.0e4);  // convert ha to m²
    auto [major_axis_m, minor_axis_m] = compute_fire_ellipse_axes(area_m2, perim_m);

    if (!ParallelDescriptor::IOProcessor()) return;

    std::ofstream out(filename, std::ios::out | std::ios::app);
    if (!out.is_open()) {
        amrex::Print() << "[FIRE] WARNING: Cannot append to fire_stats file: " << filename << "\n";
        return;
    }
     
    out << step       << ","
        << time_s     << ","
        << burned_ha  << ","
        << perim_km   << ","
        << n_active   << ","
        << head_ros   << ","
        << major_axis_m << ","
        << minor_axis_m << ","
        << Q_max << "\n";
    out.close();
}

#endif
