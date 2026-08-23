#include <ERF_Utils.H>
#include <ERF_IndexDefines.H>

using namespace amrex;

/*
 * Accumulate time averaged velocity fields
 */
void
Time_Avg_Vel_atCC (double dt_d,
                   double& t_avg_cnt,
                   MultiFab* vel_t_avg,
                   MultiFab& xvel,
                   MultiFab& yvel,
                   MultiFab& zvel)
{
    // Augment the counter
    t_avg_cnt += dt_d;

    Real dt = static_cast<Real>(dt_d);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*(vel_t_avg),TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // CC tilebox
        Box tbx = mfi.tilebox();

        // Velocity on faces
        const Array4<Real>& velx = xvel.array(mfi);
        const Array4<Real>& vely = yvel.array(mfi);
        const Array4<Real>& velz = zvel.array(mfi);

        // Time average at CC
        Array4<Real> vel_t_avg_arr = vel_t_avg->array(mfi);

        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real u_cc = myhalf * ( velx(i,j,k) + velx(i+1,j  ,k  ) );
            Real v_cc = myhalf * ( vely(i,j,k) + vely(i  ,j+1,k  ) );
            Real w_cc = myhalf * ( velz(i,j,k) + velz(i  ,j  ,k+1) );
            Real umag_cc = std::sqrt(u_cc*u_cc + v_cc*v_cc + w_cc*w_cc);
            vel_t_avg_arr(i,j,k,0) += u_cc * dt;
            vel_t_avg_arr(i,j,k,1) += v_cc * dt;
            vel_t_avg_arr(i,j,k,2) += w_cc * dt;
            vel_t_avg_arr(i,j,k,3) += umag_cc * dt;
        });
    }
}

void
Accumulate_Interval_Means (double dt_d,
                           double& t_mean_cnt,
                           MultiFab* interval_means,
                           MultiFab& xvel,
                           MultiFab& yvel,
                           MultiFab& zvel,
                           MultiFab& cons)
{
    AMREX_ALWAYS_ASSERT(interval_means != nullptr);

    t_mean_cnt += dt_d;
    const Real dt = static_cast<Real>(dt_d);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*interval_means, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Array4<const Real>& u = xvel.const_array(mfi);
        const Array4<const Real>& v = yvel.const_array(mfi);
        const Array4<const Real>& w = zvel.const_array(mfi);
        const Array4<const Real>& state = cons.const_array(mfi);
        const Array4<Real>& mean = interval_means->array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real u_cc = Real(0.5) * (u(i,j,k) + u(i+1,j,k));
            const Real v_cc = Real(0.5) * (v(i,j,k) + v(i,j+1,k));
            const Real w_cc = Real(0.5) * (w(i,j,k) + w(i,j,k+1));
            const Real theta = state(i,j,k,RhoTheta_comp) / state(i,j,k,Rho_comp);

            mean(i,j,k,0) += u_cc * dt;
            mean(i,j,k,1) += v_cc * dt;
            mean(i,j,k,2) += w_cc * dt;
            mean(i,j,k,3) += theta * dt;
            mean(i,j,k,4) += u_cc * u_cc * dt;
            mean(i,j,k,5) += v_cc * v_cc * dt;
            mean(i,j,k,6) += w_cc * w_cc * dt;
            mean(i,j,k,7) += u_cc * w_cc * dt;
            mean(i,j,k,8) += v_cc * w_cc * dt;
            mean(i,j,k,9) += w_cc * theta * dt;
        });
    }
}
