  ParmParse pp_for_pert_vels("prob");
  Real U_0 = 0.0; pp_for_pert_vels.query("U_0",U_0);

  const int klo = geomdata.Domain().smallEnd()[2];
  const int khi = geomdata.Domain().bigEnd()[2];

  // Set the x-velocity
  ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      Real ztop = z_nd(i,j,khi+1);
      Real zht  = z_nd(i,j,klo);
      x_vel_pert(i, j, k) = U_0 * ztop / (ztop - zht);
  });

  const auto dx = geomdata.CellSize();
  amrex::GpuArray<Real, AMREX_SPACEDIM> dxInv;
  dxInv[0] = 1. / dx[0];
  dxInv[1] = 1. / dx[1];
  dxInv[2] = 1. / dx[2];

  // Set the z-velocity from impenetrable condition
  ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      z_vel_pert(i, j, k) = WFromOmega(i, j, k, 0.0,
                                       x_vel_pert, y_vel_pert,
                                       mf_u, mf_v, z_nd, dxInv);
  });
