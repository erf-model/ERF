#ifndef ERF_FARSITE_ELLIPSE_H
#define ERF_FARSITE_ELLIPSE_H

#include <AMReX_Config.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_ccse-mpi.H>
#include <AMReX_Gpu.H>
#include <cmath>
#include <vector>

struct FarsiteParams
{
    amrex::Real phi_threshold   = 0.1;
    amrex::Real coeff_a         = 0.5;
    amrex::Real coeff_b         = 0.25;
    amrex::Real coeff_c         = 0.1;
    int         use_anderson_lw = 1;
    amrex::Real gaussian_sigma  = -1.0;
    amrex::Real cfl_fire        = 0.5;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real anderson_LW_ratio(amrex::Real wind_speed_mph) noexcept
{
    if (wind_speed_mph < 1.0) return 1.0;
    amrex::Real lw = 0.936 * std::exp(0.2566 * wind_speed_mph)
                   - 0.397 * std::sqrt(wind_speed_mph);
    return amrex::max(1.0, amrex::min(lw, 8.0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void LW_ratio_to_richards_coefficients(amrex::Real LW_ratio,
    amrex::Real& coeff_a, amrex::Real& coeff_b, amrex::Real& coeff_c) noexcept
{
    coeff_a = 1.0;
    coeff_c = 0.2 * coeff_a;
    coeff_b = (coeff_a + coeff_c) / (2.0 * LW_ratio);
}

/**
 * @brief Advance fire front by one fire substep.
 *
 * Spread vectors are accumulated into disp_accum across atmospheric steps.
 * A cell is stamped burned only when its accumulated displacement magnitude
 * reaches stamp_thresh = 0.5 * min(dx, dy), i.e. when the front has
 * genuinely crossed into a neighbouring cell.
 */
inline
void advance_farsite_one_step(
    amrex::MultiFab& phi,
    amrex::MultiFab& farsite_work,
    amrex::MultiFab& disp_accum,
    amrex::MultiFab& arrival_time,       // now NON-CONST: updated after stamping
    const amrex::MultiFab& vel_eff,
    const amrex::MultiFab& R_mf,
    const amrex::Geometry& geom_fire,
    amrex::Real dt_fire,
    amrex::Real current_time,            // simulation time for arrival_time stamp
    const FarsiteParams& fp)
{
    using namespace amrex;

    const auto* prob_lo = geom_fire.ProbLo();
    const auto* dx      = geom_fire.CellSize();
    const Real dxi = 1.0_rt / dx[0];
    const Real dyi = 1.0_rt / dx[1];
    const Real stamp_thresh = 0.5_rt * std::min(dx[0], dx[1]);

    // Stamp threshold: fire front stamps a new cell once accumulated
    // displacement reaches one cell width in either direction.
    // Using min(dx,dy) ensures we stamp before skipping a cell.
    const Real stamp_threshold = amrex::min(dx[0], dx[1]);

    farsite_work.setVal(0.0_rt);

    // -------------------------------------------------------------------------
    // Pass 1 (GPU): compute spread vector and accumulate into disp_accum
    // -------------------------------------------------------------------------
    for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
        const Box& bx   = mfi.tilebox();
        auto phi_arr    = phi.array(mfi);
        auto spread_arr = farsite_work.array(mfi);
        auto accum_arr  = disp_accum.array(mfi);
        auto vel_arr    = vel_eff.const_array(mfi);
        auto ros_arr    = R_mf.const_array(mfi);
        auto at_arr     = arrival_time.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) noexcept
        {
            int i = iv[0], j = iv[1], k = iv[2];

            // Only process cells on or inside the fire front
            if (phi_val > fp.phi_threshold) {
                spread_arr(iv, 0) = 0.0_rt;
                spread_arr(iv, 1) = 0.0_rt;
                return;
            }

            // Centered-difference gradient of phi
            Real dphi_dx = (phi_arr(i+1,j,k) - phi_arr(i-1,j,k)) * dxi * 0.5_rt;
            Real dphi_dy = (phi_arr(i,j+1,k) - phi_arr(i,j-1,k)) * dyi * 0.5_rt;
            Real grad_mag = std::sqrt(dphi_dx*dphi_dx + dphi_dy*dphi_dy);

            if (grad_mag < 1.0e-50) {
                spread_arr(iv, 0) = 0.0_rt;
                spread_arr(iv, 1) = 0.0_rt;
                return;
            }

            Real nx = dphi_dx / grad_mag;
            Real ny = dphi_dy / grad_mag;

            Real ux = vel_arr(iv, 0);
            Real uy = vel_arr(iv, 1);
            Real wind_mag = std::sqrt(ux*ux + uy*uy);

            Real scale = 1.0_rt;
            if (wind_mag > 1.0e-10) {
                Real wx = ux / wind_mag;
                Real wy = uy / wind_mag;
                Real cos_theta = amrex::max(-1.0_rt, amrex::min(nx*wx + ny*wy, 1.0_rt));
                Real sin_theta = std::sqrt(1.0_rt - cos_theta*cos_theta);

                Real a = fp.coeff_a, b = fp.coeff_b, c = fp.coeff_c;
                if (fp.use_anderson_lw != 0) {
                    LW_ratio_to_richards_coefficients(
                        anderson_LW_ratio(wind_mag * 2.237_rt), a, b, c);
                }
                scale = (cos_theta > 0.0_rt)
                        ? a * cos_theta + b * sin_theta
                        : b * sin_theta - c * cos_theta;
            }

            Real spread_dist = amrex::max(scale, 0.0_rt) * ros_arr(iv) * dt_fire;
            spread_arr(iv, 0) = spread_dist * nx;
            spread_arr(iv, 1) = spread_dist * ny;

            // Accumulate displacement
            accum_arr(iv, 0) += spread_dist * nx;
            accum_arr(iv, 1) += spread_dist * ny;
        });
    }

    amrex::Gpu::streamSynchronize();

    // -------------------------------------------------------------------------
    // Pass 2 (host): collect cells whose accumulated displacement >= stamp_thresh
    // -------------------------------------------------------------------------
    std::vector<Real> xarray, yarray;

    for (MFIter mfi(disp_accum); mfi.isValid(); ++mfi) {
        const Box& bx  = mfi.tilebox();
        auto accum_arr = disp_accum.array(mfi);

        for (int k = bx.smallEnd()[2]; k <= bx.bigEnd()[2]; ++k) {
            for (int j = bx.smallEnd()[1]; j <= bx.bigEnd()[1]; ++j) {
                for (int i = bx.smallEnd()[0]; i <= bx.bigEnd()[0]; ++i) {
                    Real ax = accum_arr(i, j, k, 0);
                    Real ay = accum_arr(i, j, k, 1);

                    if (std::sqrt(ax*ax + ay*ay) >= stamp_thresh) {
                        Real x = prob_lo[0] + (i + 0.5_rt) * dx[0];
                        Real y = prob_lo[1] + (j + 0.5_rt) * dx[1];
                        xarray.push_back(x + ax);
                        yarray.push_back(y + ay);
                        // Reset accumulator — displacement has been "spent"
                        accum_arr(i, j, k, 0) = 0.0_rt;
                        accum_arr(i, j, k, 1) = 0.0_rt;
                    }
                }
    }

    // MPI gather
    int npts_local = static_cast<int>(xarray.size());
    std::vector<int> npts_all(ParallelDescriptor::NProcs());
    ParallelAllGather::AllGather(&npts_local, 1, npts_all.data(),
                                 ParallelDescriptor::Communicator());
    std::vector<int> displs(ParallelDescriptor::NProcs() + 1, 0);
    for (int r = 0; r < ParallelDescriptor::NProcs(); ++r)
        displs[r+1] = displs[r] + npts_all[r];
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
    xarray_all = xarray;
    yarray_all = yarray;
#endif

    // Reset phi to unburned, then stamp propagated positions
    phi.setVal(0.0_rt);

    if (npts_global > 0) {
        Gpu::DeviceVector<Real> d_x(npts_global), d_y(npts_global);
        Gpu::copy(Gpu::hostToDevice, xarray_all.begin(), xarray_all.end(), d_x.begin());
        Gpu::copy(Gpu::hostToDevice, yarray_all.begin(), yarray_all.end(), d_y.begin());

        const Real* xp = d_x.data();
        const Real* yp = d_y.data();
        int npts = npts_global;

        if (fp.gaussian_sigma < 0.0) {
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
                amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int n) noexcept {
                    int ii = static_cast<int>(
                        std::round((xp[n] - prob_lo[0]) / dx[0] - 0.5_rt));
                    int jj = static_cast<int>(
                        std::round((yp[n] - prob_lo[1]) / dx[1] - 0.5_rt));
                    if (bx.contains(amrex::IntVect(ii, jj, 0))) {
                        p(ii, jj, 0) = amrex::min(p(ii, jj, 0), -1.0_rt);
                    }
                });
            }
        } else {
            Real sigma = (fp.gaussian_sigma > 0.0) ? fp.gaussian_sigma : dx[0];
            const Real* x_ptr = xarray_all.data();
            const Real* y_ptr = yarray_all.data();
            int npts = npts_global;
            for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
                auto p = phi.array(mfi);
                const Box& bx = mfi.validbox();
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) noexcept {
                    Real xc = prob_lo[0] + (iv[0]+0.5)*dx[0];
                    Real yc = prob_lo[1] + (iv[1]+0.5)*dx[1];
                    for (int n = 0; n < npts; ++n) {
                        Real dist = std::sqrt((xc-x_ptr[n])*(xc-x_ptr[n])
                                           + (yc-y_ptr[n])*(yc-y_ptr[n]));
                        p(iv) = amrex::min(p(iv), amrex::max(dist-sigma,-1.0_rt));
                    }
                });
            }
        }
    }

    // Step D: KEY FIX — record arrival time for newly burned cells immediately
    // so that the next substep's Step B can restore them.
    // Without this, cells stamped this substep have arrival_time < 0 and are
    // reset to phi=0 (unburned) at the start of the next substep.
    {
        const Real t_now = current_time;
        for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
            auto p  = phi.const_array(mfi);
            auto at = arrival_time.array(mfi);
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (const IntVect& iv) noexcept {
                if (p(iv) < 0.0_rt && at(iv) < 0.0_rt) {
                    at(iv) = t_now;
                }
            });
        }
    }
}

inline
amrex::Real compute_fire_dt(
    const amrex::MultiFab& R_mf,
    const amrex::Geometry& geom_fire,
    amrex::Real dt_atm,
    amrex::Real cfl_fire)
{
    amrex::Real max_ros = R_mf.max(0);
    if (max_ros < 1.0e-10) return dt_atm;
    const auto* dx = geom_fire.CellSize();
    return cfl_fire * std::min(dx[0], dx[1]) / max_ros;
}

inline
int advance_fire_subcycle(
    amrex::MultiFab& phi,
    amrex::MultiFab& farsite_work,
    amrex::MultiFab& disp_accum,
    amrex::MultiFab& arrival_time,       // non-const: updated each substep
    const amrex::MultiFab& vel_eff,
    const amrex::MultiFab& R_mf,
    const amrex::Geometry& geom_fire,
    amrex::Real dt_atm,
    amrex::Real current_time,
    const FarsiteParams& fp)
{
    amrex::Real time_remaining = dt_atm;
    int n_substeps = 0;

    while (time_remaining > 1.0e-14) {
        amrex::Real dt_fire = std::min(
            compute_fire_dt(R_mf, geom_fire, time_remaining, fp.cfl_fire),
            time_remaining);

        advance_farsite_one_step(phi, farsite_work, disp_accum, vel_eff, R_mf,
                                 geom_fire, dt_fire, fp);
        time_remaining -= dt_fire;
        ++n_substeps;

        if (n_substeps > 1000)
            amrex::Abort("advance_fire_subcycle: too many substeps");
    }
    return n_substeps;
}

#endif // ERF_FARSITE_ELLIPSE_H