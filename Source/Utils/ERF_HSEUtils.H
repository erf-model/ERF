#ifndef ERF_HSEUTIL_H_
#define ERF_HSEUTIL_H_

#include "ERF_MicrophysicsUtils.H"
#include "ERF_Constants.H"
#include "ERF_EOS.H"

/**
 * Utility functions for calculating a hydrostatic equilibrium (HSE) base state
*/

namespace HSEutils
{
    using namespace amrex;

    // to control Newton iterations in init_isentropic_hse*
    const int MAX_ITER = 10;
    const amrex::Real TOL = Real(1.e-8);

    /**
     * Function to calculate the hydrostatic density and pressure
     * from Newton-Raphson iterations
     *
     * @param[in]  m_tol     iteration tolerance
     * @param[in]  RdOCp     Rd/Cp
     * @param[out] dz        change in vertical height
     * @param[out] g         magnitude of gravity
     * @param[in]  C         sum of known terms in HSE balance
     * @param[in]  T         theta at the current cell center
     * @param[in] qt         total moisture (non-precip and precip)
     * @param[in] qv         water vapor
     * @param[in] P          pressure at cell center
     * @param[in] rd         dry density at cell center
     * @param[in] F          starting residual of non-linear eq
    */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    void
    Newton_Raphson_hse (const Real& m_tol,
                        const Real& RdOCp,
                        const Real& dz,
                        const Real& g,
                        const Real& C,
                        const Real& Th,
                        const Real& qt,
                        const Real& qv,
                        Real& P,
                        Real& rd,
                        Real& F)
    {
        int iter=0;
        int max_iter=20;
        Real eps  = Real(1.e-6);
        Real ieps = Real(1.e6);
        while (std::abs(F)>m_tol && iter<max_iter) {
            // Compute change in pressure
            Real hi   = getRhogivenThetaPress(Th, P + eps, RdOCp, qv);
            Real lo   = getRhogivenThetaPress(Th, P - eps, RdOCp, qv);
            Real dFdp = one + fourth*ieps*(hi - lo)*g*dz;
            P -= F/dFdp;
            if (P < Real(1.e3)) amrex::Warning("P < 1000 [Pa]; Domain height may be too large...");
            AMREX_ALWAYS_ASSERT(P > amrex::Real(0));

            // Diagnose density and residual
            rd = getRhogivenThetaPress(Th, P, RdOCp, qv);
            AMREX_ALWAYS_ASSERT(rd > amrex::Real(0));
            Real r_tot = rd * (one + qt);
            F = P + myhalf*r_tot*g*dz + C;
            ++iter;
        }
        if (iter>=max_iter) {
            AMREX_DEVICE_PRINTF("WARNING: HSE Newton iterations did not converge to tolerance!\n%s", "");
            AMREX_DEVICE_PRINTF("HSE Newton tol: %e %e\n",F,m_tol);
        }
    }


    /**
     * Function to calculate the hydrostatic density and pressure with constant dz
     *
     * @param[in]  r_sfc     surface density
     * @param[in]  theta     surface potential temperature
     * @param[out] r         hydrostatically balanced density profile
     * @param[out] p         hydrostatically balanced pressure profile
     * @param[in]  dz        vertical grid spacing (constant)
     * @param[in]  klo       z-index corresponding to the small end of the domain
     * @param[in]  khi       z-index corresponding to the big   end of the domain
    */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    void
    init_isentropic_hse_constant_dz (const amrex::Real& r_sfc,
                                     const amrex::Real& theta,
                                     amrex::Real* r,
                                     amrex::Real* p,
                                     const amrex::Real& dz,
                                     const int klo, const int khi)
    {
      int kstart;

      // r_sfc / p_0 are the density / pressure at klo

      // Initial guess
      Real myhalf_dz = myhalf*dz;
      r[klo] = r_sfc;
      p[klo] = p_0 - myhalf_dz * r[klo] * CONST_GRAV;

      if (klo == 0)
      {
          kstart = 1;

          // We do a Newton iteration to satisfy the EOS & HSE (with constant theta)
          bool converged_hse = false;
          Real p_hse;
          Real p_eos;

          for (int iter = 0; iter < MAX_ITER && !converged_hse; iter++)
          {
              p_hse = p_0 - myhalf_dz * r[klo] * CONST_GRAV;
              p_eos = getPgivenRTh(r[klo]*theta);

              Real A = p_hse - p_eos;

              Real dpdr = getdPdRgivenConstantTheta(r[klo],theta);

              Real drho = A / (dpdr + myhalf_dz * CONST_GRAV);

              r[klo] = r[klo] + drho;
              p[klo] = getPgivenRTh(r[klo]*theta);

              if (std::abs(drho) < TOL)
              {
                  converged_hse = true;
                  break;
              }
          }

          // if (!converged_hse) amrex::Print() << "DOING ITERATIONS AT K = " << klo << std::endl;
          // if (!converged_hse) amrex::Error("Didn't converge the iterations in init");
      } else {
          kstart = klo; // because we need to solve for r[klo] here; r[klo-1] is what was passed in r_sfc
          r[klo-1] = r_sfc; // r_sfc is really the interpolated density in the ghost cell in this case
          p[klo-1] = getPgivenRTh(r[klo-1]*theta);
      }

      // To get values at k > 0 we do a Newton iteration to satisfy the EOS (with constant theta) and
      for (int k = kstart; k <= khi; k++)
      {
          // To get values at k > 0 we do a Newton iteration to satisfy the EOS (with constant theta) and
          // to discretely satisfy HSE -- here we assume spatial_order = 2 -- we can generalize this later if needed
          bool converged_hse = false;

          r[k] = r[k-1];

          Real p_eos = getPgivenRTh(r[k]*theta);
          Real p_hse;

          for (int iter = 0; iter < MAX_ITER && !converged_hse; iter++)
          {
              Real r_avg = myhalf * (r[k-1]+r[k]);
              p_hse = p[k-1] -  dz * r_avg * CONST_GRAV;
              p_eos = getPgivenRTh(r[k]*theta);

              Real A = p_hse - p_eos;

              Real dpdr = getdPdRgivenConstantTheta(r[k],theta);
              // Gamma * p_0 * std::pow( (R_d * theta / p_0), Gamma) * std::pow(r[k], Gamma-one) ;

              Real drho = A / (dpdr + dz * CONST_GRAV);

              r[k] = r[k] + drho;
              p[k] = getPgivenRTh(r[k]*theta);

              if (std::abs(drho) < TOL * r[k-1])
              {
                  converged_hse = true;
                  // amrex::Print() << " converged " << std::endl;
                  break;
              }
          }

          // if (!converged_hse) amrex::Print() << "DOING ITERATIONS AT K = " << k << std::endl;
          // if (!converged_hse) amrex::Error("Didn't converge the iterations in init");
      }
      r[khi+1] = r[khi];
    }

    /**
     * Function to calculate the hydrostatic density and pressure with stretched dz
     *
     * @param[in]  r_sfc     surface density
     * @param[in]  theta     surface potential temperature
     * @param[out] r         hydrostatically balanced density profile
     * @param[out] p         hydrostatically balanced pressure profile
     * @param[in]  stretched_dz vertical grid spacing (stretched)
     * @param[in]  klo       z-index corresponding to the small end of the domain
     * @param[in]  khi       z-index corresponding to the big   end of the domain
    */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    void
    init_isentropic_hse_stretched_dz (const amrex::Real& r_sfc,
                                      const amrex::Real& theta,
                                      amrex::Real* r,
                                      amrex::Real* p,
                                      const Gpu::DeviceVector<Real>& stretched_dz_d,
                                      const int klo, const int khi)
    {
      int kstart;

      // r_sfc / p_0 are the density / pressure at klo

      // Initial guess
      r[klo] = r_sfc;
      p[klo] = p_0 - myhalf*stretched_dz_d[klo] * r[klo] * CONST_GRAV;

      if (klo == 0)
      {
          kstart = 1;

          // We do a Newton iteration to satisfy the EOS & HSE (with constant theta)
          bool converged_hse = false;
          Real p_hse;
          Real p_eos;

          for (int iter = 0; iter < MAX_ITER && !converged_hse; iter++)
          {
              p_hse = p_0 - myhalf*stretched_dz_d[klo] * r[klo] * CONST_GRAV;
              p_eos = getPgivenRTh(r[klo]*theta);

              Real A = p_hse - p_eos;

              Real dpdr = getdPdRgivenConstantTheta(r[klo],theta);

              Real drho = A / (dpdr + myhalf*stretched_dz_d[klo] * CONST_GRAV);

              r[klo] = r[klo] + drho;
              p[klo] = getPgivenRTh(r[klo]*theta);

              if (std::abs(drho) < TOL)
              {
                  converged_hse = true;
                  break;
              }
          }

          // if (!converged_hse) amrex::Print() << "DOING ITERATIONS AT K = " << klo << std::endl;
          // if (!converged_hse) amrex::Error("Didn't converge the iterations in init");
      } else {
          kstart = klo; // because we need to solve for r[klo] here; r[klo-1] is what was passed in r_sfc
          r[klo-1] = r_sfc; // r_sfc is really the interpolated density in the ghost cell in this case
          p[klo-1] = getPgivenRTh(r[klo-1]*theta);
      }

      // To get values at k > 0 we do a Newton iteration to satisfy the EOS (with constant theta) and
      for (int k = kstart; k <= khi; k++)
      {
          // To get values at k > 0 we do a Newton iteration to satisfy the EOS (with constant theta) and
          // to discretely satisfy HSE -- here we assume spatial_order = 2 -- we can generalize this later if needed
          bool converged_hse = false;

          r[k] = r[k-1];

          Real p_eos = getPgivenRTh(r[k]*theta);
          Real p_hse;

          for (int iter = 0; iter < MAX_ITER && !converged_hse; iter++)
          {
              Real  r_avg = myhalf * (r[k-1]+r[k]);
              Real dz_avg = myhalf * (stretched_dz_d[k-1] + stretched_dz_d[k]);
              p_hse = p[k-1] -  dz_avg * r_avg * CONST_GRAV;
              p_eos = getPgivenRTh(r[k]*theta);

              Real A = p_hse - p_eos;

              Real dpdr = getdPdRgivenConstantTheta(r[k],theta);
              // Gamma * p_0 * std::pow( (R_d * theta / p_0), Gamma) * std::pow(r[k], Gamma-one) ;

              Real drho = A / (dpdr + dz_avg * CONST_GRAV);

              r[k] = r[k] + drho;
              p[k] = getPgivenRTh(r[k]*theta);

              if (std::abs(drho) < TOL * r[k-1])
              {
                  converged_hse = true;
                  // amrex::Print() << " converged " << std::endl;
                  break;
              }
          }

          // if (!converged_hse) amrex::Print() << "DOING ITERATIONS AT K = " << k << std::endl;
          // if (!converged_hse) amrex::Error("Didn't converge the iterations in init");
      }
      r[khi+1] = r[khi];
    }


    /**
     * Function to calculate the hydrostatic density and pressure over terrain
     *
     * @param[in]  i         x-index
     * @param[in]  j         y-index
     * @param[in]  r_sfc     surface density
     * @param[in]  theta     surface potential temperature
     * @param[out] r         hydrostatically balanced density profile
     * @param[out] p         hydrostatically balanced pressure profile
     * @param[in]  z_cc      cell-center heights
     * @param[in]  khi       z-index corresponding to the big end of the domain
    */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    void
    init_isentropic_hse_terrain (int i,
                                 int j,
                                 const amrex::Real& r_sfc,
                                 const amrex::Real& theta,
                                 amrex::Real* r,
                                 amrex::Real* p,
                                 const amrex::Array4<amrex::Real const> z_cc,
                                 const int& klo, const int& khi)
    {
        int kstart;

        if (klo == 0)  {
            //
            // r_sfc / p_0 are the density / pressure at the surface
            //
            // Initial guess
            int k0 = 0;

            // Where we start the lower iteration
            kstart = 1;

            Real myhalf_dz = z_cc(i,j,k0);
            r[k0] = r_sfc;
            p[k0] = p_0 - myhalf_dz * r[k0] * CONST_GRAV;
            {
                // We do a Newton iteration to satisfy the EOS & HSE (with constant theta)
                bool converged_hse = false;
                Real p_hse;
                Real p_eos;

                for (int iter = 0; iter < MAX_ITER && !converged_hse; iter++)
                {
                    p_hse = p_0 - myhalf_dz * r[k0] * CONST_GRAV;
                    p_eos = getPgivenRTh(r[k0]*theta);

                    Real A = p_hse - p_eos;

                    Real dpdr = getdPdRgivenConstantTheta(r[k0],theta);

                    Real drho = A / (dpdr + myhalf_dz * CONST_GRAV);

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
        } else {
            kstart = klo; // because we need to solve for r[klo] here; r[klo-1] is what was passed in r_sfc
            r[klo-1] = r_sfc; // r_sfc is really the interpolated density in the ghost cell in this case
            p[klo-1] = getPgivenRTh(r[klo-1]*theta);
        }

         // To get values at k > 0 we do a Newton iteration to satisfy the EOS (with constant theta) and
         for (int k = kstart; k <= khi; k++)
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
                 p_hse = p[k-1] -  dz_loc * myhalf * (r[k-1]+r[k]) * CONST_GRAV;
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

             // if (!converged_hse) amrex::Print() << "DOING ITERATIONS AT K = " << k << std::endl;
             // if (!converged_hse) amrex::Error("Didn't converge the iterations in init");
        }
    }

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_saturation_pressure (const Real T_b, const bool use_empirical)
{
    Real p_s = erf_esatw(T_b,use_empirical);
    return p_s * Real(100.0);
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_relative_humidity (const Real p_b, const Real T_b, const bool use_empirical,
                                const int which_zone, const Real scaled_height)
{
    if (which_zone > 0) {
        if (which_zone == 1) { // z <= height
            Real p_s = compute_saturation_pressure(T_b, use_empirical);
            Real q_s = Rd_on_Rv*p_s/(p_b - p_s);
            return Real(0.014)/q_s;
        } else if (which_zone == 2) { // z > height and z <= z_tr; scaled_height = z/z_tr
            return one - Real(0.75)*std::pow(scaled_height,Real(1.25));
        } else { // z > z_tr
            return fourth;
        }
    } else {
        return one;
    }
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real vapor_mixing_ratio (const Real p_b, const Real T_b, const Real RH, const bool use_empirical,
                         int which_zone)
{
    Real p_s = compute_saturation_pressure(T_b,use_empirical);
    Real p_v = compute_vapor_pressure(p_s, RH);
    Real q_v = Rd_on_Rv*p_v/(p_b - p_v);

    if (which_zone == 1) { // z <= height
        return Real(0.014);
    } else {
        return q_v;
    }
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_F_for_temp_in_zone(const Real T_b, const Real p_b, const Real q_t, const Real eq_pot_temp,
                                const bool use_empirical, const int which_zone, const Real scaled_height)
{
    Real fac = Cp_d + Cp_l*q_t;
    Real RH  = compute_relative_humidity(p_b, T_b, use_empirical, which_zone, scaled_height);
    Real q_v = vapor_mixing_ratio(p_b, T_b, RH, use_empirical, which_zone);
    Real p_s = compute_saturation_pressure(T_b, use_empirical);
    Real p_v = compute_vapor_pressure(p_s, RH);
    return eq_pot_temp - T_b*std::pow((p_b - p_v)/p_0, -R_d/fac)*std::exp(L_v*q_v/(fac*T_b));
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_temperature (const Real p_b, const Real q_t, const Real eq_pot_temp, const bool use_empirical,
                          const int which_zone, const Real scaled_height)
{
    Real T_b = Real(200.0), delta_T; // Initial guess
    for (int iter=0; iter<20; iter++)
    {
        Real F         = compute_F_for_temp_in_zone(T_b      , p_b, q_t, eq_pot_temp, use_empirical, which_zone, scaled_height);
        Real F_plus_dF = compute_F_for_temp_in_zone(T_b+Real(1e-10), p_b, q_t, eq_pot_temp, use_empirical, which_zone, scaled_height);
        Real F_prime   = (F_plus_dF - F)/Real(1e-10);
        delta_T = -F/F_prime;
        T_b = T_b + delta_T;
    }

    if(std::fabs(delta_T) > 1e-8){
         amrex::Abort("Newton Raphson for temperature could not converge");
    }

    return T_b;
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_dewpoint_temperature (const Real T_b, const Real RH)
{
    Real T_dp, gamma, T;
    T = T_b - Real(273.15);

    Real b = Real(18.678), c = Real(257.14), d = Real(234.5);
    gamma = std::log(RH*std::exp((b - T/d)*T/(c + T)));

    T_dp = c*gamma/(b - gamma);

    return T_dp;
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_theta (const Real scaled_height, const Real theta_0,
                    const Real theta_tr, const Real z_tr, const Real T_tr)
{
    if(scaled_height <= one) {
        return theta_0 + (theta_tr - theta_0)*std::pow(scaled_height,Real(1.25));
    } else {
        return theta_tr*std::exp(CONST_GRAV/(Cp_d*T_tr)*z_tr*(scaled_height - one));
    }
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
void compute_rho (const Real& pressure, Real& theta, Real& rho, Real& q_v, Real& T_dp, Real& T_b,
                  const Real q_t, const Real eq_pot_temp, const bool use_empirical, const int which_zone,
                  const Real scaled_height, const bool T_from_theta,
                  const Real theta_0, const Real theta_tr, const Real z_tr, const Real T_tr)
{

    if (T_from_theta) {
        theta   = compute_theta(scaled_height, theta_0, theta_tr, z_tr, T_tr);
        T_b     = getTgivenPandTh(pressure, theta, (R_d/Cp_d));
    } else {
        T_b     = compute_temperature(pressure, q_t, eq_pot_temp, use_empirical, which_zone, scaled_height);
        theta   = getThgivenTandP(T_b, pressure, (R_d/Cp_d));
    }

    Real RH = compute_relative_humidity(pressure, T_b, use_empirical, which_zone, scaled_height);

    q_v     = vapor_mixing_ratio(pressure, T_b, RH, use_empirical, which_zone);

    rho     = getRhogivenTandPress(T_b, pressure, q_v);

    if (T_from_theta) {
        rho *= (one + q_v);
    } else {
        rho *= (one + q_t);
    }

    T_dp    = compute_dewpoint_temperature(T_b, RH);
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_F (const Real& p_k, const Real& p_k_minus_1, Real &theta_k, Real& rho_k, Real& q_v_k,
                Real& T_dp, Real& T_b, const Real& dz, const Real& rho_k_minus_1, const Real q_t,
                const Real eq_pot_temp, const bool use_empirical,
                const int which_zone, const Real scaled_height, const bool T_from_theta,
                const Real theta_0, const Real theta_tr, const Real z_tr, const Real T_tr)
{
    Real F;

    compute_rho(p_k, theta_k, rho_k, q_v_k, T_dp, T_b, q_t, eq_pot_temp, use_empirical, which_zone, scaled_height,
                T_from_theta, theta_0, theta_tr, z_tr, T_tr);

    if(rho_k_minus_1 == amrex::Real(0)) // This loop is for the first point above the ground
    {
        F = p_k - p_k_minus_1 + rho_k*CONST_GRAV*dz/two;
    }
    else
    {
        F = p_k - p_k_minus_1 + myhalf * (rho_k + rho_k_minus_1)*CONST_GRAV*dz;
    }

    return F;
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_p_k_in_zone (Real &p_k, const Real p_k_minus_1, Real &theta_k, Real& rho_k, Real& q_v_k,
                          Real& T_dp, Real& T_b, const Real dz, const Real rho_k_minus_1, const Real q_t,
                          const Real eq_pot_temp, const bool use_empirical,
                          const int which_zone, const Real scaled_height, const bool T_from_theta,
                          const Real theta_0, const Real theta_tr, const Real z_tr, const Real T_tr)
{
    Real delta_p_k;

#ifdef AMREX_USE_FLOAT
    Real eps = Real(1e-6);
    Real tol = Real(1e-8);
#else
    Real eps = Real(1e-10);
    Real tol = Real(1e-8);
#endif

    for(int iter=0; iter<20; iter++)
    {
        Real F         = compute_F(p_k      , p_k_minus_1, theta_k, rho_k, q_v_k, T_dp, T_b, dz, rho_k_minus_1,
                                   q_t, eq_pot_temp, use_empirical, which_zone, scaled_height, T_from_theta,
                                   theta_0, theta_tr, z_tr, T_tr);
        Real F_plus_dF = compute_F(p_k+eps, p_k_minus_1, theta_k, rho_k, q_v_k, T_dp, T_b, dz, rho_k_minus_1,
                                   q_t, eq_pot_temp, use_empirical, which_zone, scaled_height, T_from_theta,
                                   theta_0, theta_tr, z_tr, T_tr);
        Real F_prime   = (F_plus_dF - F)/eps;
        delta_p_k = -F/F_prime;
        p_k            = p_k + delta_p_k;
    }

    if (std::fabs(delta_p_k) > tol) {
        amrex::Abort("Newton Raphson for pressure could not converge");
    }

    return p_k;
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
void
init_isentropic_hse_no_terrain(Real *theta, Real* r, Real* p, Real *q_v,
                               const Real& dz, const int& khi, const Real q_t,
                               const Real eq_pot_temp, const bool use_empirical,
                               const bool T_from_theta = false,
                               const Real z_tr_1 = -one, const Real z_tr_2 = -one,
                               const Real theta_0 = amrex::Real(0), const Real theta_tr = amrex::Real(0), const Real T_tr = amrex::Real(0))
{
    Real T_b, T_dp;

    int which_zone = -1;
    Real scaled_height = amrex::Real(0);

    // Compute the quantities at z = myhalf*dz (first cell center)
    p[0] = p_0;
    if (z_tr_1 > amrex::Real(0)) {
        Real z = myhalf * dz;
        if (z <= z_tr_1) {
            which_zone = 1;
        } else if (z <= z_tr_2) {
            which_zone = 2;
        } else {
            which_zone = 3;
        }
        scaled_height = z/z_tr_2;
    }

    compute_p_k_in_zone(p[0], p_0, theta[0], r[0], q_v[0], T_dp, T_b, dz, amrex::Real(0), q_t, eq_pot_temp,
                        use_empirical, which_zone, scaled_height, T_from_theta,
                        theta_0, theta_tr, z_tr_2, T_tr);

    for (int k=1; k<=khi; k++)
    {
        p[k] = p[k-1];

        if (z_tr_1 > amrex::Real(0)) {
            Real z = (k+myhalf) * dz;
            if (z <= z_tr_1) {
                which_zone = 1;
            } else if (z <= z_tr_2) {
                which_zone = 2;
            } else {
                which_zone = 3;
            }
            scaled_height = z/z_tr_2;
        }

        compute_p_k_in_zone(p[k], p[k-1], theta[k], r[k], q_v[k], T_dp, T_b, dz, r[k-1], q_t, eq_pot_temp,
                            use_empirical, which_zone, scaled_height, T_from_theta,
                            theta_0, theta_tr, z_tr_2, T_tr);
    }

    r[khi+1] = r[khi];
}

} // namespace

#endif
