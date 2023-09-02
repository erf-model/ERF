void
init_isentropic_hse_terrain(int i, int j,
                            const Real& r_sfc, const Real& theta,
                                  Real* r,           Real* p,
                            const Array4<Real const> z_cc,
                            const int& khi)
{
  // r_sfc / p_0 are the density / pressure at the surface
  int k0  = 0;
  Real half_dz = z_cc(i,j,k0);

  // Initial guess
  r[k0] = r_sfc;
  p[k0] = p_0 - half_dz * r[k0] * CONST_GRAV;

  int MAX_ITER = 10;
  Real TOL = 1.e-8;

  {
      // We do a Newton iteration to satisfy the EOS & HSE (with constant theta)
      bool converged_hse = false;
      Real p_hse;
      Real p_eos;

      for (int iter = 0; iter < MAX_ITER && !converged_hse; iter++)
      {
          p_hse = p_0 - half_dz * r[k0] * CONST_GRAV;
          p_eos = getPgivenRTh(r[k0]*theta);

          Real A = p_hse - p_eos;

          Real dpdr = getdPdRgivenConstantTheta(r[k0],theta);

          Real drho = A / (dpdr + half_dz * CONST_GRAV);

          r[k0] = r[k0] + drho;
          p[k0] = getPgivenRTh(r[k0]*theta);

          if (std::abs(drho) < TOL)
          {
              converged_hse = true;
              break;
          }
      }

      //if (!converged_hse) amrex::Print() << "DOING ITERATIONS AT K = " << k0 << std::endl;
      //if (!converged_hse) amrex::Error("Didn't converge the iterations in init");
  }

  // To get values at k > 0 we do a Newton iteration to satisfy the EOS (with constant theta) and
  for (int k = 1; k <= khi; k++)
  {
      // To get values at k > 0 we do a Newton iteration to satisfy the EOS (with constant theta) and
      // to discretely satisfy HSE -- here we assume spatial_order = 2 -- we can generalize this later if needed
      bool converged_hse = false;

      Real dz_loc = (z_cc(i,j,k) - z_cc(i,j,k-1));

      r[k] = r[k-1];

      Real p_eos = getPgivenRTh(r[k]*theta);
      Real p_hse;

      for (int iter = 0; iter < MAX_ITER && !converged_hse; iter++)
      {

          p_hse = p[k-1] -  dz_loc * 0.5 * (r[k-1]+r[k]) * CONST_GRAV;
          p_eos = getPgivenRTh(r[k]*theta);

          Real A = p_hse - p_eos;

          Real dpdr = getdPdRgivenConstantTheta(r[k],theta);

          Real drho = A / (dpdr + dz_loc * CONST_GRAV);

          r[k] = r[k] + drho;
          p[k] = getPgivenRTh(r[k]*theta);

          if (std::abs(drho) < TOL * r[k-1])
          {
              converged_hse = true;
              //amrex::Print() << " converged " << std::endl;
              break;
          }
      }

      //if (!converged_hse) amrex::Print() << "DOING ITERATIONS AT K = " << k << std::endl;
      //if (!converged_hse) amrex::Error("Didn't converge the iterations in init");
  }
}

