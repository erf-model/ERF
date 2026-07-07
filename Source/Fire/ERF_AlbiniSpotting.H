#ifndef ERF_ALBINI_SPOTTING_H
#define ERF_ALBINI_SPOTTING_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <cmath>
#include <vector>
#include <random>
#include <string>

#include <ERF_FireParams.H>
#include <ERF_Rothermel.H>
#include <ERF_ScottSpottingTable.H>

/**
 * @file ERF_AlbiniSpotting.H
 * @brief Albini (1983) stochastic firebrand spotting model.
 *
 * Implements the lofting height, stochastic launch, forward-Euler trajectory
 * integration, and phi-negative stamping for ember spotting.
 *
 * Key references:
 *   Albini, F.A. (1983). "Potential spotting distance from wind-driven surface fires."
 *   USDA Forest Service Research Paper INT-309.
 *
 * Ported from: wildfire_levelset/src/albini_spotting.H
 * Changes from source:
 *   - InputParameters::AlbiniSpottingParams → FireParams::SpottingParams
 *   - RothermelParams replaced by RothermelComputed (from ERF_Rothermel.H)
 *   - plt_wind 3-D wind removed; fire_wind_eff is passed as the 2-D wind source
 *   - weise_data/WeiseBigingParams removed
 *   - Fuel system string passed for scott_max_spot_dist lookup
 *   - Operates on CPU host copies for GPU safety (stochastic trajectory is serial)
 *
 * Diagnostic MultiFab layout (4 components):
 *   comp 0: lofting height H_z [m] at source fire-front cells
 *   comp 1: firebrand count launched per source cell
 *   comp 2: maximum landing distance [m] from source cell
 *   comp 3: active flag (1.0) at cells that received a spot ignition
 */

// ─── unit conversion constants ─────────────────────────────────────────────
// 1 BTU/ft·s = 1.055 kW/m  (Byram fireline intensity conversion)
static constexpr amrex::Real ALBINI_BTU_FTS_TO_KW_M = 1.055;
// 1 ft/min = 0.00508 m/s
static constexpr amrex::Real ALBINI_FT_PER_MIN_TO_M_PER_S = 0.00508;

// ─── Albini lofting height ──────────────────────────────────────────────────

/**
 * @brief Compute firebrand lofting height [m].
 *
 * H_z = 12.2 × I_B^(1/3)   [Albini 1983 INT-309, eq. 1]
 *
 * Returns 0 when I_B_kW is below I_B_min_kW.
 *
 * @param[in] I_B_kW       Byram fireline intensity [kW/m]
 * @param[in] I_B_min_kW   Minimum intensity threshold [kW/m]
 * @return Lofting height H_z [m], or 0.0 if below threshold.
 */
inline amrex::Real
albini_lofting_height(amrex::Real I_B_kW, amrex::Real I_B_min_kW = 100.0)
{
    if (I_B_kW < I_B_min_kW) { return 0.0; }
    return 12.2 * std::cbrt(I_B_kW);
}

// ─── Byram fireline intensity ───────────────────────────────────────────────

/**
 * @brief Compute Byram fireline intensity [kW/m] from Rothermel quantities.
 *
 * I_B = R × w0 × h_heat / 60 × BTU_FTS_TO_KW_M
 *
 * where R is ROS in ft/min, w0 is total fuel load [lb/ft²], and h_heat is
 * heat content [BTU/lb].  The factor 1/60 converts lb·BTU/(ft·min) to
 * lb·BTU/(ft·s).
 *
 * @param[in] ros_ms  Rate of spread [m/s]
 * @param[in] fp      Fuel model parameters (load fractions, heat content)
 * @param[in] rc      Precomputed Rothermel coefficients (ros_conv for unit conversion)
 * @return Byram fireline intensity [kW/m].
 */
inline amrex::Real
albini_byram_intensity(amrex::Real ros_ms,
                       const FuelModelParams& fp,
                       const RothermelComputed& rc)
{
    // Total fuel load [lb/ft²]
    amrex::Real w0 = fp.w_d1 + fp.w_d10 + fp.w_d100 + fp.w_lh + fp.w_lw;
    // Convert ROS to ft/min using stored conversion factor (1/ros_conv = m/s → ft/min)
    amrex::Real R_ftmin = (rc.ros_conv > 0.0) ? ros_ms / rc.ros_conv : 0.0;
    // Byram intensity [BTU/ft·s] then convert to [kW/m]
    amrex::Real I_B = R_ftmin * w0 * fp.heat_content / 60.0 * ALBINI_BTU_FTS_TO_KW_M;
    return (I_B > 0.0) ? I_B : 0.0;
}

// ─── Bilinear wind interpolation ────────────────────────────────────────────

/**
 * @brief Bilinear interpolation of 2-component wind at continuous coordinates.
 *
 * @param[in]  wind_u     u-component array (comp 0 of wind MultiFab)
 * @param[in]  wind_v     v-component array (comp 1 of wind MultiFab)
 * @param[in]  nx         Number of cells in x
 * @param[in]  ny         Number of cells in y
 * @param[in]  dx         Cell size in x [m]
 * @param[in]  dy         Cell size in y [m]
 * @param[in]  prob_lo_x  Domain lower x [m]
 * @param[in]  prob_lo_y  Domain lower y [m]
 * @param[in]  x          Query x position [m]
 * @param[in]  y          Query y position [m]
 * @param[out] u_out      Interpolated u-wind [m/s]
 * @param[out] v_out      Interpolated v-wind [m/s]
 * @return true if position is within domain, false if outside (exits trajectory).
 */
inline bool
albini_interp_wind(const std::vector<amrex::Real>& wind_u,
                   const std::vector<amrex::Real>& wind_v,
                   int nx, int ny,
                   amrex::Real dx, amrex::Real dy,
                   amrex::Real prob_lo_x, amrex::Real prob_lo_y,
                   amrex::Real x, amrex::Real y,
                   amrex::Real& u_out, amrex::Real& v_out)
{
    // Cell-centre coordinates: cell (i,j) centre at prob_lo + (i+0.5)*dx
    amrex::Real fx = (x - prob_lo_x) / dx - 0.5;
    amrex::Real fy = (y - prob_lo_y) / dy - 0.5;

    int i0 = static_cast<int>(std::floor(fx));
    int j0 = static_cast<int>(std::floor(fy));

    // Domain exit check
    if (i0 < 0 || i0 >= nx - 1 || j0 < 0 || j0 >= ny - 1) {
        u_out = 0.0;
        v_out = 0.0;
        return false;
    }

    amrex::Real ax = fx - i0;
    amrex::Real ay = fy - j0;

    auto idx = [&](int ii, int jj) { return jj * nx + ii; };

    u_out = (1.0 - ax) * (1.0 - ay) * wind_u[idx(i0,   j0  )]
          + (      ax) * (1.0 - ay) * wind_u[idx(i0+1, j0  )]
          + (1.0 - ax) * (      ay) * wind_u[idx(i0,   j0+1)]
          + (      ax) * (      ay) * wind_u[idx(i0+1, j0+1)];

    v_out = (1.0 - ax) * (1.0 - ay) * wind_v[idx(i0,   j0  )]
          + (      ax) * (1.0 - ay) * wind_v[idx(i0+1, j0  )]
          + (1.0 - ax) * (      ay) * wind_v[idx(i0,   j0+1)]
          + (      ax) * (      ay) * wind_v[idx(i0+1, j0+1)];

    return true;
}

// ─── Helper: build RothermelComputed subset needed by Albini ────────────────

/**
 * @brief Build the RothermelComputed subset needed by albini_byram_intensity().
 *
 * Only ros_conv is required from RothermelComputed for the unit conversion in
 * albini_byram_intensity(); all other fields are filled to zero.
 *
 * @param[in] fp  Fuel model parameters.
 * @return RothermelComputed with ros_conv = 0.00508 m/s per ft/min.
 */
inline RothermelComputed
make_rothermel_params_for_albini(const FuelModelParams& /*fp*/)
{
    RothermelComputed rc{};
    rc.ros_conv = ALBINI_FT_PER_MIN_TO_M_PER_S;  // 1 ft/min = 0.00508 m/s
    return rc;
}

// ─── Main spotting function ─────────────────────────────────────────────────

/**
 * @brief Apply Albini (1983) stochastic firebrand spotting to the fire level-set.
 *
 * Algorithm:
 *  1. Build CPU host copies of phi, wind, and ROS (GPU-sync first).
 *  2. Identify fire-front cells: phi <= phi_thresh and ROS > 0.
 *  3. Compute Byram intensity I_B for each front cell.
 *  4. Compute lofting height H_z = 12.2 × I_B^(1/3).
 *  5. Stochastic launch: draw uniform random; launch if draw < P_launch.
 *  6. Forward-Euler trajectory over n_traj_steps sub-steps.
 *  7. Apply Scott/Albini per-fuel-model maximum distance cap.
 *  8. Apply P_catch catching probability filter.
 *  9. Skip if landing cell residual fuel < reentry_fuel_thresh.
 * 10. Stamp phi negative within spot_radius of landing location.
 * 11. Record diagnostics in albini_data.
 * 12. MPI: Allgatherv front cells so all ranks apply identical spots.
 *
 * @param[in,out] phi              Fire level-set MultiFab (modified at landing sites)
 * @param[in,out] albini_data      4-component diagnostic MultiFab (comp 0-3)
 * @param[in]     wind_eff         2-component effective wind on fire grid [m/s]
 * @param[in]     ros              Rate-of-spread MultiFab [m/s]
 * @param[in]     geom             Fire grid geometry
 * @param[in]     fp               Fuel model parameters
 * @param[in]     sp               Spotting parameters (FireParams::SpottingParams)
 * @param[in]     step             Current simulation step (for deterministic seeding)
 * @param[in]     fuel_load        Optional residual fuel MultiFab (for re-entry filter)
 * @param[in]     fuel_system      Fuel system string ("13" or "40") for distance table
 * @param[in]     fuel_model_id    Anderson/Scott-Burgan model ID for distance table lookup
 */
inline void
compute_albini_spotting(amrex::MultiFab&                  phi,
                        amrex::MultiFab&                  albini_data,
                        const amrex::MultiFab&            wind_eff,
                        const amrex::MultiFab&            ros,
                        const amrex::Geometry&            geom,
                        const FuelModelParams&            fp,
                        const FireParams::SpottingParams& sp,
                        int                               step,
                        const amrex::MultiFab*            fuel_load,
                        const std::string*                fuel_system,
                        int                               fuel_model_id = 1)
{
    using namespace amrex;

    const auto& dx      = geom.CellSize();
    const auto& prob_lo = geom.ProbLo();

    // Fire grid dimensions
    const Box domain = geom.Domain();
    int nx = domain.length(0);
    int ny = domain.length(1);

    // Build flat CPU host copies of phi, wind (u,v), ROS, and fuel_load.
    // This is necessary because stochastic trajectory integration is serial.
    Gpu::streamSynchronize();

    int ncells = nx * ny;
    std::vector<Real> h_phi(ncells, 0.0);
    std::vector<Real> h_u(ncells, 0.0);
    std::vector<Real> h_v(ncells, 0.0);
    std::vector<Real> h_ros(ncells, 0.0);
    std::vector<Real> h_fuel(ncells, 1.0);

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx    = mfi.validbox();
        auto phi_arr     = phi.const_array(mfi);
        auto wind_arr    = wind_eff.const_array(mfi);
        auto ros_arr     = ros.const_array(mfi);
        amrex::Array4<const amrex::Real> fuel_arr;
        bool has_fuel = false;
        if (fuel_load) {
            fuel_arr = fuel_load->const_array(mfi);
            has_fuel = true;
        }

        LoopOnCpu(bx, [&](int i, int j, int /*k*/) {
            int idx = j * nx + i;
            h_phi[idx]  = phi_arr(i, j, 0, 0);
            h_u[idx]    = wind_arr(i, j, 0, 0);
            h_v[idx]    = wind_arr(i, j, 0, 1);
            h_ros[idx]  = ros_arr(i, j, 0, 0);
            if (has_fuel) { h_fuel[idx] = fuel_arr(i, j, 0, 0); }
        });
    }

    // MPI: all ranks build the same global view of phi/wind/ROS.
    // Each rank owns disjoint boxes; ReduceRealSum gives each rank the complete field.
    ParallelDescriptor::ReduceRealSum(h_phi.data(),  ncells);
    ParallelDescriptor::ReduceRealSum(h_u.data(),    ncells);
    ParallelDescriptor::ReduceRealSum(h_v.data(),    ncells);
    ParallelDescriptor::ReduceRealSum(h_ros.data(),  ncells);
    if (fuel_load) {
        // h_fuel is initialised to 1.0; ReduceRealSum of disjoint patches would double-count
        // the default 1.0 values.  Initialise to 0.0, then reduce, then set unowned cells to 1.0.
        // Reset h_fuel to 0 and re-fill only owned cells:
        std::fill(h_fuel.begin(), h_fuel.end(), 0.0);
        for (MFIter mfi(*fuel_load); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto fuel_arr = fuel_load->const_array(mfi);
            LoopOnCpu(bx, [&](int i, int j, int /*k*/) {
                h_fuel[j * nx + i] = fuel_arr(i, j, 0, 0);
            });
        }
        ParallelDescriptor::ReduceRealSum(h_fuel.data(), ncells);
    }

    RothermelComputed rc = make_rothermel_params_for_albini(fp);

    // Identify fire-front cells using the global (reduced) fields.
    // All ranks see identical data after the reduce, so they generate identical spots.
    struct FrontCell {
        int  i, j;
        Real I_B;
        Real H_z;
        Real src_x, src_y;
    };

    std::vector<FrontCell> front_cells;
    constexpr Real phi_thresh = 0.1;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx  = j * nx + i;
            Real pv  = h_phi[idx];
            Real rsp = h_ros[idx];
            if (pv <= phi_thresh && rsp > 0.0) {
                Real I_B = albini_byram_intensity(rsp, fp, rc);
                Real H_z = albini_lofting_height(I_B, sp.I_B_min);
                if (H_z > 0.0) {
                    Real cx = prob_lo[0] + (i + 0.5) * dx[0];
                    Real cy = prob_lo[1] + (j + 0.5) * dx[1];
                    front_cells.push_back({i, j, I_B, H_z, cx, cy});
                }
            }
        }
    }

    int n_front = static_cast<int>(front_cells.size());
    if (n_front == 0) { return; }

    // RNG: fixed seed when sp.random_seed > 0, otherwise time+step based.
    unsigned seed = (sp.random_seed > 0)
                    ? static_cast<unsigned>(sp.random_seed + step)
                    : static_cast<unsigned>(std::time(nullptr)) + static_cast<unsigned>(step);
    // Broadcast seed so all ranks use identical RNG state.
    int iseed = static_cast<int>(seed);
    ParallelDescriptor::Bcast(&iseed, 1, ParallelDescriptor::IOProcessorNumber());
    seed = static_cast<unsigned>(iseed);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<Real> uniform(0.0, 1.0);

    // Maximum spotting distance cap from Scott/Albini table.
    std::string fuel_sys = fuel_system ? *fuel_system : "13";
    Real max_dist_cap = scott_max_spot_dist(fuel_model_id, fuel_sys);

    // Collect landing sites to stamp (avoid double-stamping during iteration).
    struct LandingSite {
        Real x, y;
        int  src_i, src_j;
        Real dist;
    };
    std::vector<LandingSite> landings;

    for (int fc = 0; fc < n_front; ++fc) {
        int  si  = front_cells[fc].i;
        int  sj  = front_cells[fc].j;
        Real I_B = front_cells[fc].I_B;
        Real H_z = front_cells[fc].H_z;
        Real cx  = front_cells[fc].src_x;
        Real cy  = front_cells[fc].src_y;

        // Stochastic launch: P_launch = P_base × min(1, I_B / I_B_min)
        Real P_launch = sp.P_base * std::min(Real(1.0), I_B / sp.I_B_min);
        if (uniform(rng) >= P_launch) { continue; }

        // Forward-Euler trajectory
        Real t_flight = (sp.terminal_velocity > 0.0) ? H_z / sp.terminal_velocity : 0.0;
        Real dt_traj  = (sp.n_traj_steps > 0) ? t_flight / sp.n_traj_steps : 0.0;

        Real px = cx;
        Real py = cy;
        bool in_domain = true;

        for (int ns = 0; ns < sp.n_traj_steps && in_domain; ++ns) {
            Real u_w = 0.0, v_w = 0.0;
            in_domain = albini_interp_wind(h_u, h_v, nx, ny, dx[0], dx[1],
                                           prob_lo[0], prob_lo[1],
                                           px, py, u_w, v_w);
            if (!in_domain) { break; }
            px += u_w * dt_traj;
            py += v_w * dt_traj;
        }

        if (!in_domain) { continue; }

        // Maximum distance cap
        Real dist = std::sqrt((px - cx) * (px - cx) + (py - cy) * (py - cy));
        if (dist > max_dist_cap) {
            // Clamp landing to max_dist_cap along trajectory direction
            if (dist > 0.0) {
                px = cx + (px - cx) / dist * max_dist_cap;
                py = cy + (py - cy) / dist * max_dist_cap;
                dist = max_dist_cap;
            }
        }

        // Catching probability filter
        if (sp.P_catch < 1.0 && uniform(rng) >= sp.P_catch) { continue; }

        // Re-entry fuel filter: skip if residual fuel below threshold
        int li = static_cast<int>((px - prob_lo[0]) / dx[0]);
        int lj = static_cast<int>((py - prob_lo[1]) / dx[1]);
        li = std::max(0, std::min(nx - 1, li));
        lj = std::max(0, std::min(ny - 1, lj));
        if (fuel_load) {
            Real rf = h_fuel[lj * nx + li];
            if (rf < sp.reentry_fuel_thresh) { continue; }
        }

        landings.push_back({px, py, si, sj, dist});

        // Update diagnostics comp 0 (H_z) and comp 1 (brand count) at source cell
        // These are done after the loop to avoid MultiFab conflicts mid-loop.
    }

    // Stamp phi at landing sites and update diagnostics
    for (const auto& ls : landings) {
        Real lx = ls.x;
        Real ly = ls.y;
        int  si = ls.src_i;
        int  sj = ls.src_j;

        // Stamp phi negative within spot_radius (SDF min-stamp, same as ignition)
        for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
            const Box& bx  = mfi.validbox();
            auto phi_arr   = phi.array(mfi);
            auto alb_arr   = albini_data.array(mfi);

            LoopOnCpu(bx, [&](int i, int j, int /*k*/) {
                Real ccx = prob_lo[0] + (i + 0.5) * dx[0];
                Real ccy = prob_lo[1] + (j + 0.5) * dx[1];
                Real d   = std::sqrt((ccx - lx) * (ccx - lx) + (ccy - ly) * (ccy - ly));
                if (d <= sp.spot_radius) {
                    // SDF stamp: phi = min(phi, d - spot_radius)
                    phi_arr(i, j, 0, 0) = std::min(phi_arr(i, j, 0, 0), d - sp.spot_radius);
                    alb_arr(i, j, 0, 3) = 1.0;  // active spot flag
                }

                // Source cell diagnostics
                if (i == si && j == sj) {
                    // Update max lofting height (comp 0) and brand count (comp 1)
                    // at source. These accumulate across multiple landing events.
                }
            });
        }
    }

    // Update source-cell diagnostics after stamping loop
    for (const auto& ls : landings) {
        int si = ls.src_i;
        int sj = ls.src_j;
        // Recompute H_z for this cell
        int idx_s = sj * nx + si;
        Real rsp  = h_ros[idx_s];
        Real I_B  = albini_byram_intensity(rsp, fp, rc);
        Real H_z  = albini_lofting_height(I_B, sp.I_B_min);

        for (MFIter mfi(albini_data); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            if (!bx.contains(IntVect(si, sj, 0))) { continue; }
            auto alb_arr = albini_data.array(mfi);
            alb_arr(si, sj, 0, 0) = H_z;            // lofting height
            alb_arr(si, sj, 0, 1) += 1.0;           // brand count
            alb_arr(si, sj, 0, 2) = std::max(alb_arr(si, sj, 0, 2), ls.dist);  // max landing dist
        }
    }

    phi.FillBoundary(geom.periodicity());
}

#endif
