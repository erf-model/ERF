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
#ifdef AMREX_USE_FLOAT
    const amrex::Real TOL = Real(1.e-4);
    const amrex::Real TEMP_STEP_REL_TOL = Real(1.e-4);
    const amrex::Real TEMP_STEP_ABS_TOL = Real(1.e-3);
    const amrex::Real TEMP_RES_REL_TOL  = Real(1.e-6);
    const amrex::Real TEMP_RES_ABS_TOL  = Real(1.e-3);
    const amrex::Real PRES_STEP_REL_TOL = Real(1.e-4);
    const amrex::Real PRES_STEP_ABS_TOL = Real(1.e-2);
    const amrex::Real PRES_RES_REL_TOL  = Real(1.e-6);
    const amrex::Real PRES_RES_ABS_TOL  = Real(1.e-1);
#else
    const amrex::Real TOL = Real(1.e-8);
    const amrex::Real TEMP_STEP_REL_TOL = Real(1.e-8);
    const amrex::Real TEMP_STEP_ABS_TOL = Real(1.e-8);
    const amrex::Real TEMP_RES_REL_TOL  = Real(1.e-10);
    const amrex::Real TEMP_RES_ABS_TOL  = Real(1.e-8);
    const amrex::Real PRES_STEP_REL_TOL = Real(1.e-8);
    const amrex::Real PRES_STEP_ABS_TOL = Real(1.e-6);
    const amrex::Real PRES_RES_REL_TOL  = Real(1.e-10);
    const amrex::Real PRES_RES_ABS_TOL  = Real(1.e-6);
#endif

    /**
     * Function to calculate the hydrostatic density and pressure
     * from Newton-Raphson iterations at constant Th & Qv
     *
     * @param[in]  m_tol      iteration tolerance
     * @param[in]  RdOCp      Rd/Cp
     * @param[out] dz         change in vertical height
     * @param[out] g          magnitude of gravity
     * @param[in]  C          sum of known terms in HSE balance
     * @param[in]  Th         theta at the current cell center
     * @param[in] qt          total moisture (non-precip and precip)
     * @param[in] qv          water vapor
     * @param[in] P           pressure at cell center
     * @param[in] rd          dry density at cell center
     * @param[in] F           starting residual of non-linear eq
     * @param[in] maintain_Th maintain theta or temperature?
    */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    void
    Newton_Raphson_hse (const Real& m_tol,
                        const Real& RdoCp,
                        const Real& dz,
                        const Real& g,
                        const Real& C,
                        const Real& Th,
                        const Real& T,
                        const Real& qt,
                        const Real& qv,
                        Real& P,
                        Real& rd,
                        Real& F,
                        const bool& maintain_Th)
    {
        int iter=0;
        int max_iter=20;
        while (std::abs(F)>m_tol && iter<max_iter) {
            // Compute change in pressure
            Real dRdP = (maintain_Th) ? iGamma * rd / P : rd / P;
            Real dFdp = one + myhalf*dRdP*g*dz;
            P -= F/dFdp;
            if (P < Real(1.e3)) amrex::Warning("P < 1000 [Pa]; Domain height may be too large...");
            AMREX_ALWAYS_ASSERT(P > amrex::Real(0));

            // Diagnose density and residual
            rd = (maintain_Th) ? getRhogivenThetaPress(Th, P, RdoCp, qv) :
                                 getRhogivenTandPress(T , P, qv);
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
                                      const amrex::Real* stretched_dz,
                                      const int klo, const int khi)
    {
      int kstart;

      // r_sfc / p_0 are the density / pressure at klo

      // Initial guess
      r[klo] = r_sfc;
      p[klo] = p_0 - myhalf*stretched_dz[klo] * r[klo] * CONST_GRAV;

      if (klo == 0)
      {
          kstart = 1;

          // We do a Newton iteration to satisfy the EOS & HSE (with constant theta)
          bool converged_hse = false;
          Real p_hse;
          Real p_eos;

          for (int iter = 0; iter < MAX_ITER && !converged_hse; iter++)
          {
              p_hse = p_0 - myhalf*stretched_dz[klo] * r[klo] * CONST_GRAV;
              p_eos = getPgivenRTh(r[klo]*theta);

              Real A = p_hse - p_eos;

              Real dpdr = getdPdRgivenConstantTheta(r[klo],theta);

              Real drho = A / (dpdr + myhalf*stretched_dz[klo] * CONST_GRAV);

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
              Real dz_avg = myhalf * (stretched_dz[k-1] + stretched_dz[k]);
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
            Real q_s = RdoRv*p_s/(p_b - p_s);
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
    Real q_v = RdoRv*p_v/(p_b - p_v);

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

struct TempNewtonEval {
    Real F;
    Real dF;
    Real dFdp;
};

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
TempNewtonEval compute_F_and_dF_for_temp_in_zone(const Real T_b, const Real p_b, const Real q_t,
                                                 const Real eq_pot_temp, const bool use_empirical,
                                                 const int which_zone, const Real scaled_height)
{
    Real fac = Cp_d + Cp_l*q_t;
    Real RH  = compute_relative_humidity(p_b, T_b, use_empirical, which_zone, scaled_height);
    Real q_v = vapor_mixing_ratio(p_b, T_b, RH, use_empirical, which_zone);
    Real p_s = compute_saturation_pressure(T_b, use_empirical);
    Real p_v = compute_vapor_pressure(p_s, RH);
    Real G   = std::pow((p_b - p_v)/p_0, -R_d/fac) * std::exp(L_v*q_v/(fac*T_b));
    Real F   = eq_pot_temp - T_b*G;

    Real dp_s_dT = Real(100.0) * erf_dtesatw(T_b, use_empirical);
    Real dp_v_dT = Real(0.0);
    Real dq_v_dT = Real(0.0);
    Real dp_v_dp = Real(0.0);
    Real dq_v_dp = Real(0.0);

    if (which_zone == 1) {
        Real q_s;
        Real dq_s_dT;
        erf_qsatw(T_b, p_b, q_s);
        erf_dtqsatw(T_b, p_b, dq_s_dT);

        Real q_ref = Real(0.014);
        Real RH_ref = q_ref / q_s;
        Real dRH_dT = -q_ref * dq_s_dT / (q_s * q_s);
        dp_v_dT = dp_s_dT * RH_ref + p_s * dRH_dT;
        dp_v_dp = q_ref / RdoRv;
    } else {
        dp_v_dT = RH * dp_s_dT;
        dq_v_dT = RdoRv * p_b * dp_v_dT / ((p_b - p_v) * (p_b - p_v));
        dq_v_dp = -RdoRv * p_v / ((p_b - p_v) * (p_b - p_v));
    }

    Real dlnA_dT = (R_d / fac) * (dp_v_dT / (p_b - p_v));
    Real dB_dT   = (L_v / fac) * (dq_v_dT / T_b - q_v / (T_b * T_b));
    Real dF      = -G * (Real(1.0) + T_b * (dlnA_dT + dB_dT));

    Real dlnA_dp = - (R_d / fac) * (Real(1.0) - dp_v_dp) / (p_b - p_v);
    Real dB_dp   = (L_v / fac) * (dq_v_dp / T_b);
    Real dFdp    = -T_b * G * (dlnA_dp + dB_dp);

    return {F, dF, dFdp};
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_temperature (const Real p_b, const Real q_t, const Real eq_pot_temp, const bool use_empirical,
                          const int which_zone, const Real scaled_height, const Real T_guess = Real(200.0))
{
#ifdef AMREX_USE_FLOAT
    Real T_b = (T_guess > Real(0.0)) ? T_guess : Real(200.0); // Initial guess
    Real delta_T = Real(0.0);
    Real residual_T = Real(0.0);
    Real max_rel_step = Real(0.1);
    bool converged = false;

    for (int iter=0; iter<20; iter++)
    {
        TempNewtonEval newton_eval = compute_F_and_dF_for_temp_in_zone(T_b, p_b, q_t, eq_pot_temp,
                                                                       use_empirical, which_zone, scaled_height);
        residual_T = newton_eval.F;
        delta_T = -residual_T / newton_eval.dF;

        Real step_tol = TEMP_STEP_ABS_TOL + TEMP_STEP_REL_TOL * std::abs(T_b);
        Real resid_tol = TEMP_RES_ABS_TOL + TEMP_RES_REL_TOL * std::abs(eq_pot_temp);
        if (std::abs(delta_T) <= step_tol && std::abs(residual_T) <= resid_tol) {
            converged = true;
            break;
        }

        Real max_step = max_rel_step * std::abs(T_b);
        if (std::abs(delta_T) > max_step) {
            delta_T = std::copysign(max_step, delta_T);
        }

        T_b += delta_T;
    }

    if (!converged) {
         amrex::Print() << "Newton Raphson for temperature could not converge; last residual = "
                        << residual_T << ", last delta_T = " << delta_T << std::endl;
         amrex::Abort("Newton Raphson for temperature could not converge");
    }

    return T_b;
#else
    ignore_unused(T_guess);
    Real T_b = Real(200.0), delta_T; // Initial guess

    Real eps = Real(1.e-10) * T_b;

    for (int iter=0; iter<20; iter++)
    {
        Real F         = compute_F_for_temp_in_zone(T_b    , p_b, q_t, eq_pot_temp, use_empirical, which_zone, scaled_height);
        Real F_plus_dF = compute_F_for_temp_in_zone(T_b+eps, p_b, q_t, eq_pot_temp, use_empirical, which_zone, scaled_height);
        Real F_prime   = (F_plus_dF - F)/eps;
        delta_T = -F/F_prime;
        T_b = T_b + delta_T;
    }

    if (std::fabs(delta_T) > TOL * T_b) {
        amrex::Abort("Newton Raphson for temperature could not converge");
    }

    return T_b;
#endif
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
        T_b     = getTgivenPandTh(pressure, theta, RdoCp);
    } else {
        if (!(T_b > Real(0.0))) {
            T_b = eq_pot_temp;
        }
        T_b     = compute_temperature(pressure, q_t, eq_pot_temp, use_empirical, which_zone, scaled_height, T_b);
        theta   = getThgivenTandP(T_b, pressure, RdoCp);
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
#ifdef AMREX_USE_FLOAT
    Real delta_p_k = Real(0.0);
    Real residual_p = Real(0.0);
    Real max_rel_step = Real(0.1);
    Real rel_eps = Real(1e-6);
    bool converged = false;

    for(int iter=0; iter<20; iter++)
    {
        residual_p = compute_F(p_k, p_k_minus_1, theta_k, rho_k, q_v_k, T_dp, T_b, dz, rho_k_minus_1,
                               q_t, eq_pot_temp, use_empirical, which_zone, scaled_height, T_from_theta,
                               theta_0, theta_tr, z_tr, T_tr);

        Real dp = rel_eps * std::max(Real(1.0), std::abs(p_k));
        Real p_plus = p_k + dp;
        Real p_minus = p_k - dp;

        Real theta_tmp = theta_k;
        Real rho_tmp = rho_k;
        Real qv_tmp = q_v_k;
        Real Tdp_tmp = T_dp;
        Real Tb_tmp  = T_b;
        Real F_plus = compute_F(p_plus, p_k_minus_1, theta_tmp, rho_tmp, qv_tmp, Tdp_tmp, Tb_tmp, dz,
                                rho_k_minus_1, q_t, eq_pot_temp, use_empirical, which_zone,
                                scaled_height, T_from_theta, theta_0, theta_tr, z_tr, T_tr);

        theta_tmp = theta_k;
        rho_tmp = rho_k;
        qv_tmp = q_v_k;
        Tdp_tmp = T_dp;
        Tb_tmp  = T_b;
        Real F_minus = compute_F(p_minus, p_k_minus_1, theta_tmp, rho_tmp, qv_tmp, Tdp_tmp, Tb_tmp, dz,
                                 rho_k_minus_1, q_t, eq_pot_temp, use_empirical, which_zone,
                                 scaled_height, T_from_theta, theta_0, theta_tr, z_tr, T_tr);

        Real F_prime = (F_plus - F_minus) / (Real(2.0) * dp);
        if (std::abs(F_prime) < Real(1e-14)) {
            F_prime = std::copysign(Real(1e-14), F_prime == Real(0.0) ? Real(1.0) : F_prime);
        }

        delta_p_k = -residual_p / F_prime;
        Real step_tol = PRES_STEP_ABS_TOL + PRES_STEP_REL_TOL * std::abs(p_k);
        Real resid_scale = std::max(std::abs(p_k), std::max(std::abs(p_k_minus_1), Real(1.0)));
        Real resid_tol = PRES_RES_ABS_TOL + PRES_RES_REL_TOL * resid_scale;
        if (std::abs(delta_p_k) <= step_tol && std::abs(residual_p) <= resid_tol) {
            converged = true;
            break;
        }

        Real max_step = max_rel_step * std::abs(p_k);
        if (std::abs(delta_p_k) > max_step) {
            delta_p_k = std::copysign(max_step, delta_p_k);
        }
        p_k += delta_p_k;
    }

    if (!converged) {
        amrex::Print() << "Newton Raphson for pressure could not converge; last residual = "
                       << residual_p << ", last delta_p_k = " << delta_p_k << std::endl;
        amrex::Abort("Newton Raphson for pressure could not converge");
    }
    return p_k;
#else
    Real delta_p_k;

    Real eps = Real(1e-10);

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

    if (std::fabs(delta_p_k) > TOL) {
        amrex::Abort("Newton Raphson for pressure could not converge");
    }

    return p_k;
#endif
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
    Real T_b = eq_pot_temp, T_dp;

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
