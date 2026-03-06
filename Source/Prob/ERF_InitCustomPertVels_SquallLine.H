
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        Real z_cc = 0.125 * (
          z_nd(i,j,k  ) + z_nd(i+1,j,k  ) + z_nd(i,j+1,k  ) + z_nd(i+1,j+1,k  ) +
          z_nd(i,j,k+1) + z_nd(i+1,j,k+1) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1) );

        x_vel_pert(i,j,k) = -12.0*std::max(0.0, (2.5e3 - z_cc)/2.5e3);
    });

