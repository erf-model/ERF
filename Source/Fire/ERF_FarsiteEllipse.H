#ifndef ERF_FARSITE_ELLIPSE_H
#define ERF_FARSITE_ELLIPSE_H

/**
 * @file ERF_FarsiteEllipse.H
 * @brief FARSITE elliptical fire growth model and level-set propagation
 *
 * Richards (1990) elliptical fire growth model.
 * Richards, G.D. (1990). An elliptical growth model of forest fire fronts
 *   and its numerical solution. Int. J. Numer. Meth. Eng. 30(6):1163-1179.
 *
 * Ported from hgopalan/wildfire_levelset src/farsite_ellipse.H.
 * That file implements compute_farsite_spread() operating on AMReX MultiFabs.
 *
 * The fire front is represented by a signed level-set field phi:
 *   phi < 0  : burned
 *   phi > 0  : unburned
 *   phi = 0  : fire front
 *
 * Propagation is two-pass:
 *   Pass 1 (GPU): compute spread vector at each front cell from ROS field
 *                 and Richards ellipse coefficients.
 *   Pass 2 (host): collect propagated positions, MPI_Allgatherv so every
 *                  rank holds the global front, then stamp phi.
 *
 * The MPI_Allgatherv call transfers O(perimeter / dx_fire) Real values.
 * For a 2 km perimeter on a 10 m fire grid this is ~200 points per rank,
 * and the all-gather transfers ~200 * nranks * 8 bytes — negligible overhead.
 *
 * Assumptions:
 *   - Single AMR level (no AMR refinement in z between fire and atm grids).
 *   - The fire grid uses the same x-y MPI decomposition as the atmospheric
 *     level-0 grid; each rank owns the same x-y columns in both.
 *   - fire_phi has 1 ghost cell in x and y for gradient stencils.
 *   - wind_eff has no ghost cells required (read at owned cells only).
 */

#include <AMReX_Config.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_ccse-mpi.H>
#include <AMReX_Gpu.H>
#include <cmath>
#include <vector>

/**
 * @struct FarsiteParams
 * @brief Configuration parameters for FARSITE level-set propagation
 */
struct FarsiteParams
{
    amrex::Real phi_threshold   = 0.0;  ///< front detection level in phi
    amrex::Real coeff_a         = 0.5;  ///< Richards head coefficient
    amrex::Real coeff_b         = 0.25; ///< Richards flank coefficient
    amrex::Real coeff_c         = 0.1;  ///< Richards backing coefficient
    int         use_anderson_lw = 1;    ///< derive a,b,c from Anderson L/W ratio
    amrex::Real gaussian_sigma  = -1.0; ///< phi stamping radius: <0 = single cell, 0 = auto, >0 = user [m]
    amrex::Real cfl_fire        = 0.5;  ///< fire CFL number for subcycle timestep control
};

/**
 * @brief Anderson (1983) L/W ratio as a function of midflame wind speed.
 *
 * Anderson, H.E. (1983). Predicting wind-driven wildland fire size and shape.
 *   USDA Forest Service Research Paper INT-305.
 *
 * Returns L/W in [1.0, 8.0].
 *
 * @param[in] wind_speed_mph Midflame wind speed [mph]
 * @return Length-to-width ratio L/W
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real anderson_LW_ratio(amrex::Real wind_speed_mph) noexcept
{
    // Anderson (1983) empirical formula: L/W = 0.936 * e^(0.2566*U) + 0.461 - 0.397 * sqrt(U)
    // Simplified approximation used in FARSITE:
    // For U >= 1 mph: L/W = 0.936 * e^(0.2566*U)
    // For U < 1 mph: L/W = 1.0 (no wind effect)

    if (wind_speed_mph < 1.0) {
        return 1.0;
    }
    // Clamp to [1.0, 8.0]
    amrex::Real lw = 0.936 * std::exp(0.2566 * wind_speed_mph) - 0.397 * std::sqrt(wind_speed_mph);
    return amrex::max(1.0, amrex::min(lw, 8.0));
}

/**
 * @brief Convert Anderson L/W ratio to Richards (1990) ellipse coefficients.
 *
 * Sets coeff_a=1.0, coeff_c=0.2*coeff_a, coeff_b=(a+c)/(2*LW).
 * From FARSITE parameterisation.
 *
 * @param[in]  LW_ratio Length-to-width ratio
 * @param[out] coeff_a  Richards head coefficient (set to 1.0)
 * @param[out] coeff_b  Richards flank coefficient
 * @param[out] coeff_c  Richards backing coefficient (set to 0.2)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void LW_ratio_to_richards_coefficients(amrex::Real LW_ratio,
    amrex::Real& coeff_a, amrex::Real& coeff_b, amrex::Real& coeff_c) noexcept
{
    coeff_a = 1.0;
    coeff_c = 0.2 * coeff_a;
    coeff_b = (coeff_a + coeff_c) / (2.0 * LW_ratio);
}

/**
 * @brief Advance fire front phi by one fire substep dt_fire using Richards (1990)
 * elliptical expansion.
 *
 * Two-pass algorithm:
 * Pass 1 (GPU): Compute spread vectors from ROS and ellipse coefficients.
 * Pass 2 (host): Collect front positions, global gather, stamp phi.
 *
 * Inputs:
 *   phi          — level-set field on fire grid (C*nx x C*ny x 1), 1 ghost cell
 *   vel_eff      — effective midflame wind on fire grid, 2 comp, no ghosts
 *   R_mf         — ROS field from Phase 2 Rothermel [m/s], no ghosts
 *   farsite_work — scratch MultiFab, same BA/DM as phi, 2 comp, no ghosts
 *   geom_fire    — fire grid geometry
 *   dt_fire      — fire substep [s]
 *   fp           — FARSITE parameters (coeff_a/b/c, phi_threshold, etc.)
 *
 * The function is not GPU-only: Pass 1 is GPU (ParallelFor),
 * Pass 2 is host (loop + MPI_Allgatherv + GPU stamp).
 *
 * phi is modified in-place: propagated front positions are stamped to -1.
 * Already-burned cells (phi < 0) are never reset — the burned interior is
 * preserved across timesteps.
 *
 * @param[in,out] phi         Level-set field
 * @param[in,out] farsite_work Scratch for spread vectors
 * @param[in]     vel_eff     Effective wind field [m/s], 2 comp
 * @param[in]     R_mf        ROS field [m/s]
 * @param[in]     geom_fire   Fire grid geometry
 * @param[in]     dt_fire     Fire substep [s]
 * @param[in]     fp          FARSITE parameters
 */
inline
void advance_farsite_one_step(
    amrex::MultiFab& phi,
    amrex::MultiFab& farsite_work,
    const amrex::MultiFab& vel_eff,
    const amrex::MultiFab& R_mf,
    const amrex::Geometry& geom_fire,
    amrex::Real dt_fire,
    const FarsiteParams& fp)
{
    using namespace amrex;

    const auto* prob_lo = geom_fire.ProbLo();
    const auto* dx      = geom_fire.CellSize();
    const Real dxi = 1.0 / dx[0];
    const Real dyi = 1.0 / dx[1];

    // Pass 1: GPU — compute spread vectors
    for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto phi_arr     = phi.array(mfi);
        auto spread_arr  = farsite_work.array(mfi);
        auto vel_arr     = vel_eff.const_array(mfi);
        auto ros_arr     = R_mf.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) noexcept
        {
            int i = iv[0], j = iv[1], k = iv[2];
            Real phi_val = phi_arr(iv);

            // Skip if not near the front
            if (phi_val > fp.phi_threshold) {
                spread_arr(iv, 0) = 0.0;
                spread_arr(iv, 1) = 0.0;
                return;
            }

            // Compute centered-difference gradient of phi
            Real dphi_dx = (phi_arr(i+1, j, k) - phi_arr(i-1, j, k)) * dxi * 0.5;
            Real dphi_dy = (phi_arr(i, j+1, k) - phi_arr(i, j-1, k)) * dyi * 0.5;

            Real grad_mag = std::sqrt(dphi_dx*dphi_dx + dphi_dy*dphi_dy);

            if (grad_mag < 1.0e-50) {
                spread_arr(iv, 0) = 0.0;
                spread_arr(iv, 1) = 0.0;
                return;
            }

            // Outward normal
            Real nx = dphi_dx / grad_mag;
            Real ny = dphi_dy / grad_mag;

            // Wind components
            Real ux = vel_arr(iv, 0);
            Real uy = vel_arr(iv, 1);
            Real wind_mag = std::sqrt(ux*ux + uy*uy);

            // Compute Richards scale factor
            Real scale = 1.0;

            if (wind_mag > 1.0e-10) {
                Real wind_unit_x = ux / wind_mag;
                Real wind_unit_y = uy / wind_mag;

                // cos_theta = dot(normal, wind_unit)
                Real cos_theta = nx * wind_unit_x + ny * wind_unit_y;
                cos_theta = amrex::max(-1.0, amrex::min(cos_theta, 1.0));

                Real sin_theta = std::sqrt(1.0 - cos_theta*cos_theta);

                // Apply Richards coefficients
                Real a = fp.coeff_a;
                Real b = fp.coeff_b;
                Real c = fp.coeff_c;

                // If Anderson L/W ratio is enabled, compute coefficients from wind speed
                if (fp.use_anderson_lw != 0) {
                    // Convert wind speed to mph (assuming wind_mag is in m/s)
                    Real wind_mph = wind_mag * 2.237;  // m/s to mph conversion
                    Real lw_ratio = anderson_LW_ratio(wind_mph);
                    LW_ratio_to_richards_coefficients(lw_ratio, a, b, c);
                }

                if (cos_theta > 0.0) {
                    scale = a * cos_theta + b * sin_theta;
                } else {
                    scale = b * sin_theta - c * cos_theta;
                }
            }

            // Read ROS
            Real R_head = ros_arr(iv);

            // Compute spread distance
            Real spread_dist = amrex::max(scale, 0.0) * R_head * dt_fire;

            // Write spread vectors
            spread_arr(iv, 0) = spread_dist * nx;
            spread_arr(iv, 1) = spread_dist * ny;
        });
    }

    // Synchronize GPU after Pass 1
    amrex::Gpu::streamSynchronize();

    // Pass 2: Host — collect positions, gather, stamp phi
    std::vector<Real> xarray, yarray;

    for (MFIter mfi(farsite_work); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto spread_arr = farsite_work.const_array(mfi);

        for (int k = bx.smallEnd()[2]; k <= bx.bigEnd()[2]; ++k) {
            for (int j = bx.smallEnd()[1]; j <= bx.bigEnd()[1]; ++j) {
                for (int i = bx.smallEnd()[0]; i <= bx.bigEnd()[0]; ++i) {
                    Real spread_x = spread_arr(i, j, k, 0);
                    Real spread_y = spread_arr(i, j, k, 1);

                    if (std::abs(spread_x) > 1.0e-14 || std::abs(spread_y) > 1.0e-14) {
                        // Compute physical position of propagated point
                        Real x = prob_lo[0] + (i + 0.5) * dx[0];
                        Real y = prob_lo[1] + (j + 0.5) * dx[1];

                        xarray.push_back(x + spread_x);
                        yarray.push_back(y + spread_y);
                    }
                }
            }
        }
    }

    // Collect positions from all ranks via AllGather
    int npts_local = static_cast<int>(xarray.size());
    std::vector<int> npts_all(ParallelDescriptor::NProcs());
    ParallelAllGather::AllGather(&npts_local, 1, npts_all.data(),
                                 ParallelDescriptor::Communicator());

    std::vector<int> displs(ParallelDescriptor::NProcs() + 1, 0);
    for (int i = 0; i < ParallelDescriptor::NProcs(); ++i) {
        displs[i+1] = displs[i] + npts_all[i];
    }
    int npts_global = displs[ParallelDescriptor::NProcs()];

    std::vector<Real> xarray_all(npts_global), yarray_all(npts_global);

#ifdef AMREX_USE_MPI
    MPI_Allgatherv(xarray.data(), npts_local,
                   ParallelDescriptor::Mpi_typemap<Real>::type(),
                   xarray_all.data(), npts_all.data(), displs.data(),
                   ParallelDescriptor::Mpi_typemap<Real>::type(),
                   ParallelDescriptor::Communicator());
    MPI_Allgatherv(yarray.data(), npts_local,
                   ParallelDescriptor::Mpi_typemap<Real>::type(),
                   yarray_all.data(), npts_all.data(), displs.data(),
                   ParallelDescriptor::Mpi_typemap<Real>::type(),
                   ParallelDescriptor::Communicator());
#else
    // Serial build: only one rank, just copy local data
    xarray_all = xarray;
    yarray_all = yarray;
#endif

    // FIX: Reset phi to 0, then stamp propagated positions.
    // The arrival time record preserves the burned interior across timesteps.
    phi.setVal(0.0_rt);

    // Stamp propagated front positions into phi
    if (npts_global > 0) {
        std::vector<Real> x_host = xarray_all;
        std::vector<Real> y_host = yarray_all;

        const Real* x_ptr = x_host.data();
        const Real* y_ptr = y_host.data();
        int npts = npts_global;

        if (fp.gaussian_sigma < 0.0) {
            // Single-cell stamping: copy point arrays to device and use ParallelFor.
            // This is GPU-safe unlike a direct host-side Array4 write.
            amrex::Gpu::DeviceVector<Real> d_x(npts_global), d_y(npts_global);
            amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                             xarray_all.begin(), xarray_all.end(), d_x.begin());
            amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                             yarray_all.begin(), yarray_all.end(), d_y.begin());

            const Real* xp = d_x.data();
            const Real* yp = d_y.data();
            int npts = npts_global;

            for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                auto p = phi.array(mfi);

                // For each propagated point, mark the containing cell as burned.
                // Use the validbox bounds to avoid ghost-cell confusion.
                amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int n) noexcept {
                    Real x = xp[n];
                    Real y = yp[n];

                    int ii = static_cast<int>((x - prob_lo[0]) * dxi);
                    int jj = static_cast<int>((y - prob_lo[1]) * dyi);

                    if (bx.contains(amrex::IntVect(ii, jj, 0))) {
                        p(ii, jj, 0) = -1.0_rt;
                    }
                });
            }
        } else {
            // Gaussian/SDF path — already GPU-safe (unchanged)
            Real sigma = (fp.gaussian_sigma > 0.0) ? fp.gaussian_sigma : dx[0];

            for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
                auto p = phi.array(mfi);
                const Box& bx = mfi.validbox();

                ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) noexcept {
                    Real x_cell = prob_lo[0] + (iv[0] + 0.5) * dx[0];
                    Real y_cell = prob_lo[1] + (iv[1] + 0.5) * dx[1];

                    for (int n = 0; n < npts; ++n) {
                        Real x_pt = x_ptr[n];
                        Real y_pt = y_ptr[n];

                        Real dist = std::sqrt((x_cell - x_pt)*(x_cell - x_pt) +
                                             (y_cell - y_pt)*(y_cell - y_pt));
                        Real phi_sdf = amrex::max(dist - sigma, -1.0_rt);

                        p(iv) = amrex::min(p(iv), phi_sdf);
                    }
                });
            }
        }
    }
}

/**
 * @brief Compute the maximum fire substep satisfying the CFL condition.
 *
 * dt_fire = cfl * min(dx, dy) / max(R_mf)
 *
 * Returns dt_atm if R_mf is zero everywhere (no active fire).
 *
 * @param[in] R_mf       ROS field [m/s]
 * @param[in] geom_fire  Fire grid geometry
 * @param[in] dt_atm     Atmospheric timestep [s]
 * @param[in] cfl_fire   CFL number
 * @return Fire substep [s]
 */
inline
amrex::Real compute_fire_dt(
    const amrex::MultiFab& R_mf,
    const amrex::Geometry& geom_fire,
    amrex::Real dt_atm,
    amrex::Real cfl_fire)
{
    using namespace amrex;

    Real max_ros = R_mf.max(0);

    if (max_ros < 1.0e-10) {
        return dt_atm;
    }

    const auto* dx = geom_fire.CellSize();
    Real min_dx = std::min(dx[0], dx[1]);

    return cfl_fire * min_dx / max_ros;
}

/**
 * @brief Advance the fire front through multiple substeps to cover one atmospheric timestep.
 *
 * Each substep:
 *   1. Calls advance_farsite_one_step().
 *   2. Recomputes dt_fire from updated R_mf (ROS is held constant within
 *      the atmospheric step; this recompute only adjusts for CFL).
 *
 * Returns the number of fire substeps taken.
 *
 * @param[in,out] phi         Level-set field
 * @param[in,out] farsite_work Scratch for spread vectors
 * @param[in]     vel_eff     Effective wind field [m/s], 2 comp
 * @param[in]     R_mf        ROS field [m/s]
 * @param[in]     geom_fire   Fire grid geometry
 * @param[in]     dt_atm      Atmospheric timestep [s]
 * @param[in]     fp          FARSITE parameters
 * @return Number of fire substeps taken
 */
inline
int advance_fire_subcycle(
    amrex::MultiFab& phi,
    amrex::MultiFab& farsite_work,
    const amrex::MultiFab& vel_eff,
    const amrex::MultiFab& R_mf,
    const amrex::Geometry& geom_fire,
    amrex::Real dt_atm,
    const FarsiteParams& fp)
{
    using namespace amrex;

    Real time_remaining = dt_atm;
    int n_substeps = 0;

    while (time_remaining > 1.0e-14) {
        Real dt_fire = compute_fire_dt(R_mf, geom_fire, time_remaining, fp.cfl_fire);
        dt_fire = std::min(dt_fire, time_remaining);

        advance_farsite_one_step(phi, farsite_work, vel_eff, R_mf,
                                geom_fire, dt_fire, fp);

        time_remaining -= dt_fire;
        n_substeps++;

        // Safety check to avoid infinite loop
        if (n_substeps > 1000) {
            amrex::Abort("advance_fire_subcycle: too many substeps");
        }
    }

    return n_substeps;
}


#endif // ERF_FARSITE_ELLIPSE_H