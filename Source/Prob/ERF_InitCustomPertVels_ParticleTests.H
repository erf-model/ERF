  ParmParse pp_for_pert_vels("prob");
  Real U_0 = zero; pp_for_pert_vels.query("U_0",U_0);

  const Real ztop = geomdata.ProbHi()[2];

  // Set the x-velocity
  // NOTE: z_nd(i,j,0) is the terrain surface height. This works when the
  // refined patch abuts the terrain (klo=0) but will fail if L1 is elevated
  // above the surface (k=0 would be out of bounds for the local z_nd data).
  ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      Real zbot = z_nd(i,j,0);
      x_vel_pert(i, j, k) = U_0 * ztop / (ztop - zbot);
  });

  const auto dx = geomdata.CellSize();
  amrex::GpuArray<Real, AMREX_SPACEDIM> dxInv;
  dxInv[0] = one / dx[0];
  dxInv[1] = one / dx[1];
  dxInv[2] = one / dx[2];

  // Set the z-velocity from impenetrable condition
  ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      z_vel_pert(i, j, k) = WFromOmega(i, j, k, zero,
                                       x_vel_pert, y_vel_pert,
                                       mf_u, mf_v, z_nd, dxInv);
  });
