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
    if (wind_speed_mph < 1.0) { return 1.0; }
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

inline
void advance_farsite_one_step(
    amrex::MultiFab& phi,
    amrex::MultiFab& farsite_work,
    amrex::MultiFab& disp_accum,
    amrex::MultiFab& arrival_time,
    const amrex::MultiFab& vel_eff,
    const amrex::MultiFab& R_mf,
    const amrex::Geometry& geom_fire,
    amrex::Real dt_fire,
    amrex::Real current_time,
    const FarsiteParams& fp)
{
    using namespace amrex;

    const auto* prob_lo = geom_fire.ProbLo();
    const auto* dx      = geom_fire.CellSize();
    const Real dxi = 1.0 / dx[0];
    const Real dyi = 1.0 / dx[1];
    const Real stamp_threshold = amrex::min(dx[0], dx[1]);

    farsite_work.setVal(0.0_rt);

    // -------------------------------------------------------------------------
    // Pass 1: accumulate displacement; stamp absolute target when threshold met
    // -------------------------------------------------------------------------
    for (MFIter mfi(R_mf, false); mfi.isValid(); ++mfi) {
        const Box& bx   = mfi.validbox();
        auto phi_arr    = phi.array(mfi);
        auto spread_arr = farsite_work.array(mfi);
        auto disp_arr   = disp_accum.array(mfi);
        auto vel_arr    = vel_eff.const_array(mfi);
        auto ros_arr    = R_mf.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) noexcept
        {
            int i = iv[0], j = iv[1], k = iv[2];

            if (phi_arr(iv) > fp.phi_threshold) { return; }

            Real dphi_dx = (phi_arr(i+1,j,k) - phi_arr(i-1,j,k)) * dxi * 0.5;
            Real dphi_dy = (phi_arr(i,j+1,k) - phi_arr(i,j-1,k)) * dyi * 0.5;
            Real grad_mag = std::sqrt(dphi_dx*dphi_dx + dphi_dy*dphi_dy);
            if (grad_mag < 1.0e-10) { return; }

            Real nx = dphi_dx / grad_mag;
            Real ny = dphi_dy / grad_mag;

            Real ux = vel_arr(iv, 0);
            Real uy = vel_arr(iv, 1);
            Real wind_mag = std::sqrt(ux*ux + uy*uy);

            Real scale = 1.0;
            if (wind_mag > 1.0e-10) {
                Real wx = ux / wind_mag;
                Real wy = uy / wind_mag;
                Real cos_t = amrex::max(-1.0, amrex::min(nx*wx + ny*wy, 1.0));
                Real sin_t = std::sqrt(1.0 - cos_t*cos_t);
                Real a = fp.coeff_a, b = fp.coeff_b, c = fp.coeff_c;
                if (fp.use_anderson_lw != 0) {
                    LW_ratio_to_richards_coefficients(
                        anderson_LW_ratio(wind_mag * 2.237), a, b, c);
                }
                scale = (cos_t > 0.0)
                    ? a * cos_t + b * sin_t
                    : b * sin_t - c * cos_t;
            }

            Real ds = amrex::max(scale, 0.0) * ros_arr(iv) * dt_fire;
            disp_arr(iv, 0) += ds * nx;
            disp_arr(iv, 1) += ds * ny;

            Real accum_mag = std::sqrt(disp_arr(iv,0)*disp_arr(iv,0)
                                     + disp_arr(iv,1)*disp_arr(iv,1));
            if (accum_mag >= stamp_threshold) {
                Real x_src = prob_lo[0] + (i + 0.5) * dx[0];
                Real y_src = prob_lo[1] + (j + 0.5) * dx[1];
                spread_arr(iv, 0) = x_src + disp_arr(iv, 0);
                spread_arr(iv, 1) = y_src + disp_arr(iv, 1);
                disp_arr(iv, 0) = 0.0;
                disp_arr(iv, 1) = 0.0;
            }
        });
    }

    amrex::Gpu::streamSynchronize();

    // -------------------------------------------------------------------------
    // Pass 2: collect stamped positions on host
    // -------------------------------------------------------------------------
    std::vector<Real> xarray, yarray;

    for (MFIter mfi(farsite_work, false); mfi.isValid(); ++mfi) {
        const Box& bx   = mfi.validbox();
        auto spread_arr = farsite_work.const_array(mfi);
        for (int k = bx.smallEnd()[2]; k <= bx.bigEnd()[2]; ++k)
            for (int j = bx.smallEnd()[1]; j <= bx.bigEnd()[1]; ++j)
                for (int i = bx.smallEnd()[0]; i <= bx.bigEnd()[0]; ++i) {
                    Real sx = spread_arr(i,j,k,0);
                    Real sy = spread_arr(i,j,k,1);
                    if (sx != 0.0_rt || sy != 0.0_rt) {
                        xarray.push_back(sx);
                        yarray.push_back(sy);
                    }
                }
    }

    // -------------------------------------------------------------------------
    // MPI all-gather: use MPI_Allgather directly to avoid AMReX wrapper issues.
    // Falls back to a simple copy on non-MPI builds.
    // -------------------------------------------------------------------------
    int npts_local = static_cast<int>(xarray.size());
    int nranks     = ParallelDescriptor::NProcs();

    std::vector<int> npts_all(nranks, 0);
    std::vector<int> displs(nranks + 1, 0);

#ifdef AMREX_USE_MPI
    // Gather per-rank counts
    MPI_Allgather(&npts_local, 1, MPI_INT,
                  npts_all.data(), 1, MPI_INT,
                  ParallelDescriptor::Communicator());
#else
    npts_all[0] = npts_local;
#endif

    for (int r = 0; r < nranks; ++r)
        displs[r+1] = displs[r] + npts_all[r];
    int npts_global = displs[nranks];

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

    amrex::Print() << "[FIRE] Stamping " << npts_global << " new perimeter points\n";

    // -------------------------------------------------------------------------
    // Step A: reset phi to 0
    // Step B: restore burned interior from arrival_time
    // Step C: stamp new perimeter positions
    // Step D: update arrival_time for newly burned cells
    // -------------------------------------------------------------------------
    phi.setVal(0.0_rt);

    for (MFIter mfi(phi, false); mfi.isValid(); ++mfi) {
        auto p  = phi.array(mfi);
        auto at = arrival_time.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (const IntVect& iv) noexcept {
            if (at(iv) >= 0.0_rt) p(iv) = -1.0_rt;
        });
    }
    amrex::Gpu::streamSynchronize();

    if (npts_global > 0) {
        if (fp.gaussian_sigma < 0.0) {
            amrex::Gpu::DeviceVector<Real> d_x(npts_global), d_y(npts_global);
            amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                             xarray_all.begin(), xarray_all.end(), d_x.begin());
            amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                             yarray_all.begin(), yarray_all.end(), d_y.begin());
            const Real* xp = d_x.data();
            const Real* yp = d_y.data();
            int npts = npts_global;

            for (MFIter mfi(phi, false); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
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
            for (MFIter mfi(phi, false); mfi.isValid(); ++mfi) {
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

    // Step D: record arrival time for newly burned cells
    {
        const Real t_now = current_time;
        for (MFIter mfi(phi, false); mfi.isValid(); ++mfi) {
            auto p  = phi.const_array(mfi);
            auto at = arrival_time.array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (const IntVect& iv) noexcept {
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
    using namespace amrex;
    Real max_ros = R_mf.max(0);
    if (max_ros < 1.0e-10) { return dt_atm; }
    const auto* dx = geom_fire.CellSize();
    return cfl_fire * amrex::min(dx[0], dx[1]) / max_ros;
}

inline
int advance_fire_subcycle(
    amrex::MultiFab& phi,
    amrex::MultiFab& farsite_work,
    amrex::MultiFab& disp_accum,
    amrex::MultiFab& arrival_time,
    const amrex::MultiFab& vel_eff,
    const amrex::MultiFab& R_mf,
    const amrex::Geometry& geom_fire,
    amrex::Real dt_atm,
    amrex::Real current_time,
    const FarsiteParams& fp)
{
    using namespace amrex;

    Real time_remaining = dt_atm;
    int n_substeps = 0;

    while (time_remaining > 1.0e-14) {
        Real dt_fire = amrex::min(
            compute_fire_dt(R_mf, geom_fire, time_remaining, fp.cfl_fire),
            time_remaining);

        advance_farsite_one_step(phi, farsite_work, disp_accum,
                                 arrival_time,
                                 vel_eff, R_mf,
                                 geom_fire, dt_fire,
                                 current_time + (dt_atm - time_remaining),
                                 fp);

        time_remaining -= dt_fire;
        ++n_substeps;

        if (n_substeps > 1000) {
            amrex::Abort("advance_fire_subcycle: too many substeps");
        }
    }

    return n_substeps;
}

#endif // ERF_FARSITE_ELLIPSE_H