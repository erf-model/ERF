#include <Utils.H>

using namespace amrex;

/*
 * Accumulate time averaged velocity fields
 */
void
Time_Avg_Vel_atCC (const Real& dt,
                   Real& t_avg_cnt,
                   MultiFab* vel_t_avg,
                   MultiFab& xvel,
                   MultiFab& yvel,
                   MultiFab& zvel)
{
    // Augment the counter
    t_avg_cnt += dt;

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
            Real u_cc = 0.5 * ( velx(i,j,k) + velx(i+1,j  ,k  ) );
            Real v_cc = 0.5 * ( vely(i,j,k) + vely(i  ,j+1,k  ) );
            Real w_cc = 0.5 * ( velz(i,j,k) + velz(i  ,j  ,k+1) );
            Real umag_cc = std::sqrt(u_cc*u_cc + v_cc*v_cc + w_cc*w_cc);
            vel_t_avg_arr(i,j,k,0) += u_cc;
            vel_t_avg_arr(i,j,k,1) += v_cc;
            vel_t_avg_arr(i,j,k,2) += w_cc;
            vel_t_avg_arr(i,j,k,3) += umag_cc;
        });
    }
}
