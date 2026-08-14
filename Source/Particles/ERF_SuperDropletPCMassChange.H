#ifndef SUPERDROPLET_PC_MASSCHANGE_H_
#define SUPERDROPLET_PC_MASSCHANGE_H_

#include <cmath>
#include "ERF_Constants.H"

#ifdef ERF_USE_PARTICLES

/*! \brief Namespace with classes and functions for condensation/evaporation
 *  (liquid <--> vapour) */
namespace SDMassChangeUtils_LV
{
    /*! \brief Ventilation factor as a function of droplet radius (Bayley et al., 2025, Eq. 3) */
    template <typename RT /*!< real-type */ >
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RT ventilationFactor ( const RT a_R,       /*!< radius */
                           const RT a_alpha1,  /*!< fit coefficient alpha_1 */
                           const RT a_beta1,   /*!< fit exponent beta_1 */
                           const RT a_alpha2,  /*!< fit coefficient alpha_2 */
                           const RT a_beta2,   /*!< fit exponent beta_2 */
                           const RT a_fcap     /*!< maximum ventilation factor */ ) noexcept
    {
        RT x1 = a_alpha1 * std::pow(a_R, a_beta1);
        RT x2 = a_alpha2 * std::pow(a_R, a_beta2);
        RT fv = RT(1.0) + (x1*x2) / (x1 + x2);
        return std::min(fv, a_fcap);
    }

#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
    /*! \brief Phase change equation (in terms of R) */
    template <typename RT /*!< real-type */ >
    struct dRdt
    {
        RT L;         /*!< latent heat of vaporization (condensate) */
        RT K;         /*!< thermal conductivity (condensate) */
        RT Rv;        /*!< gas constant of air with vapour */
        RT rho_l;     /*!< density of condensate */

        /*! \brief Right-hand-side of the phase change ODE */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT operator () ( const RT a_R,   /*!< radius */
                         const RT a_S,   /*!< saturation ratio */
                         const RT a_T,   /*!< temperature */
                         const RT a_e_s, /*!< saturation pressure */
                         const RT a_D,   /*!< mol. diffusion coeff */
                         const RT a_a,   /*!< curvature coeff */
                         const RT a_b,   /*!< solute coeff */
                         const RT a_N_s, /*!< total solute moles */
                         const RT a_f_v  /*!< ventilation factor */ ) const noexcept
        {
            RT lambda_v = 2*a_D/std::sqrt(8*a_T*Rv/PI);
            RT Kn = lambda_v/a_R;
            RT dcf = (1+Kn) / (1+2*Kn*(1+Kn));

            RT F_k = ( L/(Rv*a_T) - RT(1.0)) * ((L*rho_l) / (K*a_T));
            RT F_d = (rho_l*Rv*a_T) / (dcf*a_D*a_e_s);

            RT R_inv = RT(1.0)/a_R;
            RT R_inv_cubed = R_inv*R_inv*R_inv;

            RT alpha = (a_S-RT(1.0)) / (F_k + F_d);
            RT retval = alpha*R_inv;

            RT beta = -(a_a/a_T) / (F_k + F_d);
            retval += beta*R_inv*R_inv;

            RT gamma = a_b * a_N_s / (F_k + F_d);
            retval += gamma*R_inv_cubed*R_inv;

            return a_f_v * retval;
        }

    };
#endif

    /*! \brief Phase change equation (in terms of R^2) */
    template <typename RT /*!< real-type */ >
    struct dRsqdt
    {
        RT L;         /*!< latent heat of vaporization (condensate) */
        RT K;         /*!< thermal conductivity (condensate) */
        RT Rv;        /*!< gas constant of air with vapour */
        RT rho_l;     /*!< density of condensate */

        /*! \brief Right-hand-side of the phase change ODE */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT rhs_func ( const RT a_R_sq, /*!< radius squared */
                      const RT a_S,    /*!< saturation ratio */
                      const RT a_T,    /*!< temperature */
                      const RT a_e_s,  /*!< saturation pressure */
                      const RT a_D,    /*!< mol. diffusion coeff */
                      const RT a_a,    /*!< curvature coeff */
                      const RT a_b,    /*!< solute coeff */
                      const RT a_N_s,  /*!< total solute moles */
                      const RT a_f_v   /*!< ventilation factor */ ) const noexcept
        {
            RT lambda_v = 2*a_D/std::sqrt(8*a_T*Rv/PI);
            RT Kn = lambda_v/std::sqrt(a_R_sq);
            RT dcf = (1+Kn) / (1+2*Kn*(1+Kn));

            RT F_k = ( L/(Rv*a_T) - RT(1.0)) * ((L*rho_l) / (K*a_T));
            RT F_d = (rho_l*Rv*a_T) / (dcf*a_D*a_e_s);

            RT R_inv = RT(1.0)/std::sqrt(a_R_sq);
            RT R_inv_cubed = R_inv*R_inv*R_inv;

            RT alpha = RT(2.0) * (a_S-RT(1.0)) / (F_k + F_d);
            RT retval = alpha;

            RT beta = -RT(2.0) * (a_a/a_T) / (F_k + F_d);
            retval += beta*R_inv;

            RT gamma = RT(2.0) * a_b * a_N_s / (F_k + F_d);
            retval += gamma*R_inv_cubed;

            return a_f_v * retval;
        }

        /*! \brief Jacobian of right-hand-side of the phase change ODE */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT rhs_jac ( const RT a_R_sq, /*!< radius squared */
                     const RT a_T,    /*!< temperature */
                     const RT a_e_s,  /*!< saturation pressure */
                     const RT a_D,    /*!< mol. diffusion coeff */
                     const RT a_a,    /*!< curvature coeff */
                     const RT a_b,    /*!< solute coeff */
                     const RT a_N_s,  /*!< total solute moles */
                     const RT a_f_v   /*!< ventilation factor */ ) const noexcept
        {
            RT lambda_v = 2*a_D/std::sqrt(8*a_T*Rv/PI);
            RT Kn = lambda_v/std::sqrt(a_R_sq);
            RT dcf = (1+Kn) / (1+2*Kn*(1+Kn));

            RT F_k = ( L/(Rv*a_T) - RT(1.0)) * ((L*rho_l) / (K*a_T));
            RT F_d = (rho_l*Rv*a_T) / (dcf*a_D*a_e_s);

            RT R_inv = RT(1.0)/std::sqrt(a_R_sq);
            RT R_inv_3 = R_inv*R_inv*R_inv;
            RT R_inv_5 = R_inv_3*R_inv*R_inv;

            RT retval = RT(0.0);

            RT beta = -RT(2.0) * (a_a/a_T) / (F_k + F_d);
            retval -= myhalf * beta*R_inv_3;

            RT gamma = RT(2.0) * a_b * a_N_s / (F_k + F_d);
            retval -= myhalf * three*gamma*R_inv_5;

            return a_f_v * retval;
        }

    };

    /*! \brief Scalar Newton solver for phase change equation
     *
     * Solves the following nonlinear equation:
     * mu * u - F(u) - R = 0,
     * where:
     *   u: solution variable
     *   mu: constant
     *   R: right-hand-side (constant)
     *   F(u): function
    */
    template<typename NE /*!< Nonlinear equation */, typename RT /*!< real-type */>
    struct NewtonSolver
    {
        const NE  m_ne;      /*!< nonlinear equation */

        RT  m_rtol;     /*!< relative tolerance */
        RT  m_atol;     /*!< absolute tolerance */
        RT  m_stol;     /*!< step size tolerance */
        int m_maxits;   /*!< max number of iterations */

        /*! \brief solve the nonlinear equation */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator()  (  RT&             a_u,            /*!< solution variable */
                            RT&             a_r,            /*!< right-hand-side */
                            const RT&       a_mu,           /*!< mu */
                            const RT&       a_S,            /*!< saturation ratio */
                            const RT&       a_T,            /*!< temperature */
                            const RT&       a_e_s,          /*!< saturation pressure */
                            const RT&       a_D,            /*!< mol. diff. coeff */
                            const RT&       a_a,            /*!< curvature coeff */
                            const RT&       a_b,            /*!< solute coeff */
                            const RT&       a_N_s,          /*!< total solute moles */
                            const RT&       a_f_v,          /*!< ventilation factor */
                            RT&             a_res_norm_a,   /*!< absolute norm at exit */
                            RT&             a_res_norm_r    /*!< relative norm at exit */,
                            bool&           a_converged     /*!< convergence status at exit */
                         ) const
        {
            a_converged = false;
            RT res_norm0 = RT(0.0);

            for (int k = 0; k < m_maxits; k++) {
                RT residual = a_mu * a_u
                              - (   a_r
                                  + m_ne.rhs_func( a_u, a_S, a_T, a_e_s, a_D, a_a, a_b, a_N_s, a_f_v ) );
                a_res_norm_a = std::sqrt(residual*residual);

                if (k == 0) {
                    if (a_res_norm_a > 0) {
                        res_norm0 = a_res_norm_a;
                    } else {
                        res_norm0 = RT(1.0);
                    }
                }
                a_res_norm_r = a_res_norm_a / res_norm0;

                if (a_res_norm_a <= m_atol) {
                    a_converged = true;
                    break;
                }
                if (a_res_norm_r <= m_rtol) {
                    a_converged = true;
                    break;
                }
                if (!std::isfinite(a_res_norm_a)) {
                    a_converged = false;
                    break;
                }

                RT slope = a_mu - m_ne.rhs_jac( a_u, a_T, a_e_s, a_D, a_a, a_b, a_N_s, a_f_v );
                RT du = RT(0.0);
                du = - residual / slope;

                RT du_norm = std::sqrt(du*du);
                RT u_norm = std::sqrt(a_u*a_u);
                if (du_norm/u_norm <= m_stol) {
                    a_converged = true;
                    break;
                }

                a_u += du;
                if (a_u <= RT(0.0)) {
                    a_converged = false;
                    break;
                }
            }
        }
    };

    /*! \brief Implicit and explicit time integrators for the phase change equation */
    template<   typename ODE /*!< ODE */,
                typename NewtonSolver /*!< Newton solver */,
                /**
                 * @brief Implicit and explicit time integrators for the phase change equation.
                 *
                 * @tparam ODE ODE type.
                 * @tparam NewtonSolver Newton solver type.
                 * @tparam RT Real-type.
                 */
                typename RT /*!< real-type */ >
    struct TI
    {
        const ODE m_ode; /*!< ODE */
        const NewtonSolver m_newton; /*!< Newton solver */

        RT m_t_final;   /*!< final time */
        RT m_max_steps; /*!< max number of timesteps */
        RT m_S;         /*!< saturation ratio */
        RT m_T;         /*!< temperature */
        RT m_e_s;       /*!< saturation pressure */
        RT m_D;         /*!< mol. diff. coeff */
        RT m_a;         /*!< coefficient of curvature */
        RT m_b;         /*!< coefficient of solute */
        RT m_N_s;       /*!< total solute moles */

        RT m_cfl;       /*!< CFL */
        RT m_atol;      /*!< absolute tolerance (for adaptive dt) */
        RT m_rtol;      /*!< absolute tolerance (for adaptive dt) */
        RT m_stol;      /*!< solution update tolerance for exit due to steady state */

        bool m_adapt_dt;  /*!< use error-based adaptive dt? */
        bool m_verbose;   /*!< verbosity */

        bool m_ventilation;   /*!< include ventilation factor? */
        RT m_vent_alpha1;     /*!< ventilation factor fit coefficient alpha_1 */
        RT m_vent_beta1;      /*!< ventilation factor fit exponent beta_1 */
        RT m_vent_alpha2;     /*!< ventilation factor fit coefficient alpha_2 */
        RT m_vent_beta2;      /*!< ventilation factor fit exponent beta_2 */
        RT m_vent_fcap;       /*!< maximum ventilation factor */

        /*! \brief Ventilation factor at the given radius-squared, frozen over a
         *  step to keep the implicit Jacobian exact (unity if disabled) */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT ventFactor (const RT& a_u) const
        {
            if (!m_ventilation) { return RT(1.0); }
            return ventilationFactor( std::sqrt(a_u),
                                      m_vent_alpha1, m_vent_beta1,
                                      m_vent_alpha2, m_vent_beta2,
                                      m_vent_fcap );
        }

        /*! \brief Compute timestep from stiffness estimate */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT computeTimestep (const RT& a_u) const
        {
            RT tau = m_ode.rhs_jac(a_u, m_T, m_e_s, m_D, m_a, m_b, m_N_s, ventFactor(a_u));
            return m_cfl / std::sqrt(tau*tau);
        }

        /*! \brief Compute stiffness estimate (tau) */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT computeTau (const RT& a_u) const
        {
            return m_ode.rhs_jac(a_u, m_T, m_e_s, m_D, m_a, m_b, m_N_s, ventFactor(a_u));
        }

        /*! \brief Limit timestep to not exceed final time */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT limitTimestep (RT dt, RT cur_time) const
        {
            if ((cur_time + dt) > m_t_final) {
                dt = m_t_final - cur_time;
            }
            return dt;
        }

        /*! \brief Check if timestep is too small to continue */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool isTimestepTooSmall (RT dt, RT tau) const
        {
            return (dt < (RT(1.0e-12) * m_cfl / std::sqrt(tau*tau)))
                && (dt < (RT(1.0e-12) * m_t_final));
        }

        /*! \brief Evaluate ODE right-hand side */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT evalRHS (const RT& a_u) const
        {
            return m_ode.rhs_func(a_u, m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s, ventFactor(a_u));
        }

        /*! \brief Print verbose step information */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void printStepInfo (RT cur_time, RT dt, RT tau, RT radius, RT snorm) const
        {
            if (m_verbose) {
                AMREX_DEVICE_PRINTF("Time %1.2e, dt = %1.2e, cfl = %1.1e, radius = %1.4e, snorm = %1.1e\n",
                                    cur_time, dt, dt * std::sqrt(tau*tau), radius, snorm);
            }
        }

        /*! \brief Print verbose step info with Newton solver details */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void printStepInfoNewton (RT cur_time, RT dt, RT tau, RT radius, RT snorm,
                                  RT res_norm_a, RT res_norm_r, bool converged) const
        {
            if (m_verbose) {
                AMREX_DEVICE_PRINTF("Time %1.2e, dt = %1.2e, cfl = %1.1e, radius = %1.4e, snorm = %1.1e\n",
                                    cur_time, dt, dt * std::sqrt(tau*tau), radius, snorm);
                AMREX_DEVICE_PRINTF("    norms = %1.3e (abs), %1.3e (rel), converged = %s\n",
                                    res_norm_a, res_norm_r, (converged ? "yes" : "no"));
            }
        }

        /*! \brief 3rd-order, 4-stage Bogacki-Shampine explicit RK method
         *  with 2nd order embedded method */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void rk3bs ( RT& a_u,  /*!< solution */
                     bool& a_success  /*!< success/failure flag */ ) const
        {
            RT cur_time = RT(0.0);
            a_success = true;

            RT tau = m_ode.rhs_jac(a_u, m_T, m_e_s, m_D, m_a, m_b, m_N_s, ventFactor(a_u));
            RT dt = m_cfl / std::sqrt(tau*tau);

            RT dt_new = dt;
            RT a_u_old = a_u;

            int n_step = 0;
            while ((cur_time < m_t_final) && (n_step < m_max_steps)) {

                if (!m_adapt_dt) {
                    tau = m_ode.rhs_jac(a_u, m_T, m_e_s, m_D, m_a, m_b, m_N_s, ventFactor(a_u));
                    dt = m_cfl / std::sqrt(tau*tau);
                } else {
                    dt = dt_new;
                }

                if ((cur_time + dt) > m_t_final) {
                    dt = m_t_final - cur_time;
                }
                if (!std::isfinite(dt)) {
                    a_success = false;
                    break;
                }

                RT u_new = RT(0.0);
                bool step_success = false;
                while (!step_success) {

                    if (    (dt < (RT(1.0e-12)*m_cfl/std::sqrt(tau*tau)))
                         && (dt < (RT(1.0e-12)*m_t_final)) ) {
                        break;
                    }

                    RT u1 = a_u;
                    if (u1 <= 0) {
                        dt *= myhalf;
                        continue;
                    }
                    RT f1 = m_ode.rhs_func(u1, m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s, ventFactor(a_u));

                    RT u2 = a_u + myhalf*dt*f1;
                    if (u2 <= 0) {
                        dt *= myhalf;
                        continue;
                    }
                    RT f2 = m_ode.rhs_func(u2, m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s, ventFactor(a_u));

                    RT u3 = a_u + RT(0.75)*dt*f2;
                    if (u3 <= 0) {
                        dt *= myhalf;
                        continue;
                    }
                    RT f3 = m_ode.rhs_func(u3, m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s, ventFactor(a_u));

                    RT u4 = a_u + (RT(1.0)/RT(9.0))*dt * (RT(2.0)*f1 + three*f2 + RT(4.0)*f3);
                    if (u4 <= 0) {
                        dt *= myhalf;
                        continue;
                    }
                    RT f4 = m_ode.rhs_func(u4, m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s, ventFactor(a_u));

                    u_new = u4;

                    if (m_adapt_dt) {
                        RT u_embed = a_u + (RT(1.0)/RT(24.0))*dt * (RT(7.0)*f1 + RT(6.0)*f2 + RT(8.0)*f3 + three*f4);
                        RT err = std::sqrt((u_new-u_embed)*(u_new-u_embed));
                        RT tol = m_atol + m_rtol * std::max(a_u, a_u_old);
                        RT E = err / tol;
                        dt_new = dt / std::cbrt(E);
                    }

                    if (std::isfinite(u_new)) {
                        if (u_new > 0) {
                            step_success = true;
                            break;
                        }
                    }
                    dt *= myhalf;
                }

                if (step_success) {

                    RT snorm = std::sqrt((a_u-u_new)*(a_u-u_new)/(a_u*a_u));
                    a_u_old = a_u;
                    a_u = u_new;
                    cur_time += dt;

                    if (m_verbose) {
                        AMREX_DEVICE_PRINTF( "Time %1.2e, dt = %1.2e, cfl = %1.1e, radius = %1.4e, snorm = %1.1e\n",
                                cur_time, dt, dt * std::sqrt(tau*tau), std::sqrt(a_u), snorm);
                    }
                    if (snorm < m_stol) {
                        break;
                    }

                } else {

                    a_success = false;
                    break;

                }

                n_step++;
            }

            return;
        }

        /*! \brief 4th-order, 4-stage explicit Runge-Kutta method */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void rk4 ( RT& a_u,  /*!< solution */
                   bool& a_success  /*!< success/failure flag */ ) const
        {
            RT cur_time = RT(0.0);
            a_success = true;

            int n_step = 0;
            while ((cur_time < m_t_final) && (n_step < m_max_steps)) {

                RT tau = computeTau(a_u);
                RT dt = m_cfl / std::sqrt(tau*tau);
                dt = limitTimestep(dt, cur_time);
                if (!std::isfinite(dt)) {
                    a_success = false;
                    break;
                }

                RT u_new = RT(0.0);
                bool step_success = false;
                while (!step_success) {

                    if (isTimestepTooSmall(dt, tau)) { break; }

                    RT u1 = a_u;
                    if (u1 <= 0) { dt *= myhalf; continue; }
                    RT f1 = evalRHS(u1);

                    RT u2 = a_u + myhalf*dt*f1;
                    if (u2 <= 0) { dt *= myhalf; continue; }
                    RT f2 = evalRHS(u2);

                    RT u3 = a_u + myhalf*dt*f2;
                    if (u3 <= 0) { dt *= myhalf; continue; }
                    RT f3 = evalRHS(u3);

                    RT u4 = a_u + dt*f3;
                    if (u4 <= 0) { dt *= myhalf; continue; }
                    RT f4 = evalRHS(u4);

                    u_new = a_u + dt*(f1 + RT(2.0)*f2 + RT(2.0)*f3 + f4)/RT(6.0);

                    if (std::isfinite(u_new) && (u_new > 0)) {
                        step_success = true;
                        break;
                    }
                    dt *= myhalf;
                }

                if (step_success) {
                    RT snorm = std::sqrt((a_u-u_new)*(a_u-u_new)/(a_u*a_u));
                    a_u = u_new;
                    cur_time += dt;
                    printStepInfo(cur_time, dt, tau, std::sqrt(a_u), snorm);
                    if (snorm < m_stol) { break; }
                } else {
                    a_success = false;
                    break;
                }

                n_step++;
            }

            return;
        }

        /*! \brief 1st order implicit backward Euler method */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void be ( RT& a_u,  /*!< solution */
                  bool& a_success  /*!< success/failure flag */ ) const
        {
            RT cur_time = RT(0.0);
            a_success = true;

            int n_step = 0;
            while ((cur_time < m_t_final) && (n_step < m_max_steps)) {

                RT tau = computeTau(a_u);
                RT dt = m_cfl / std::sqrt(tau*tau);
                dt = limitTimestep(dt, cur_time);
                if (!std::isfinite(dt)) {
                    a_success = false;
                    break;
                }

#ifdef AMREX_USE_FLOAT
                RT res_norm_a = std::numeric_limits<float>::max();
                RT res_norm_r = std::numeric_limits<float>::max();
#else
                RT res_norm_a = std::numeric_limits<double>::max();
                RT res_norm_r = std::numeric_limits<double>::max();
#endif
                bool converged = false;

                RT u_new = RT(0.0);
                bool step_success = false;
                while (!step_success) {

                    if (isTimestepTooSmall(dt, tau)) { break; }

                    RT mu = RT(1.0) / dt;
                    RT rhs = mu * a_u;
                    u_new = a_u;
                    m_newton(u_new, rhs, mu,
                             m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s, ventFactor(a_u),
                             res_norm_a, res_norm_r, converged);

                    if (converged && std::isfinite(u_new) && (u_new > 0)) {
                        step_success = true;
                        break;
                    }
                    dt *= myhalf;
                }

                if (step_success) {
                    RT snorm = std::sqrt((a_u-u_new)*(a_u-u_new)/(a_u*a_u));
                    a_u = u_new;
                    cur_time += dt;
                    printStepInfoNewton(cur_time, dt, tau, std::sqrt(a_u), snorm,
                                        res_norm_a, res_norm_r, converged);
                    if (snorm < m_stol) { break; }
                } else {
                    a_success = false;
                    break;
                }

                n_step++;
            }

            return;
        }

        /*! \brief 2nd-order implicit Crank-Nicolson method */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void cn ( RT& a_u,  /*!< solution */
                  bool& a_success  /*!< success/failure flag */ ) const
        {
            RT cur_time = RT(0.0);
            a_success = true;

            int n_step = 0;
            while ((cur_time < m_t_final) && (n_step < m_max_steps)) {

                RT tau = computeTau(a_u);
                RT dt = m_cfl / std::sqrt(tau*tau);
                dt = limitTimestep(dt, cur_time);
                if (!std::isfinite(dt)) {
                    a_success = false;
                    break;
                }

                RT res_norm_a = std::numeric_limits<RT>::max();
                RT res_norm_r = std::numeric_limits<RT>::max();
                bool converged = false;

                RT u_new = RT(0.0);
                bool step_success = false;
                while (!step_success) {

                    if (isTimestepTooSmall(dt, tau)) { break; }

                    RT mu = RT(1.0) / (myhalf*dt);
                    RT u1 = a_u;
                    RT f1 = evalRHS(u1);

                    RT u2 = u1;
                    RT rhs = mu * (a_u + myhalf*dt*f1);
                    m_newton(u2, rhs, mu,
                             m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s, ventFactor(a_u),
                             res_norm_a, res_norm_r, converged);
                    if (u2 <= 0) { dt *= myhalf; continue; }
                    RT f2 = evalRHS(u2);

                    u_new = a_u + myhalf * dt * (f1 + f2);

                    if (converged && std::isfinite(u_new) && (u_new > 0)) {
                        step_success = true;
                        break;
                    }
                    dt *= myhalf;
                }

                if (step_success) {
                    RT snorm = std::sqrt((a_u-u_new)*(a_u-u_new)/(a_u*a_u));
                    a_u = u_new;
                    cur_time += dt;
                    printStepInfoNewton(cur_time, dt, tau, std::sqrt(a_u), snorm,
                                        res_norm_a, res_norm_r, converged);
                    if (snorm < m_stol) { break; }
                } else {
                    a_success = false;
                    break;
                }

                n_step++;
            }

            return;
        }

        /*! \brief 2nd-order, 2-stage diagonally-implicit RK method */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void dirk212 ( RT& a_u,  /*!< solution */
                       bool& a_success  /*!< success/failure flag */ ) const
        {
            RT cur_time = RT(0.0);
            a_success = true;

            int n_step = 0;
            while ((cur_time < m_t_final) && (n_step < m_max_steps)) {

                RT tau = computeTau(a_u);
                RT dt = m_cfl / std::sqrt(tau*tau);
                dt = limitTimestep(dt, cur_time);
                if (!std::isfinite(dt)) {
                    a_success = false;
                    break;
                }

                bool converged = false;
                RT res_norm_a = std::numeric_limits<RT>::max();
                RT res_norm_r = std::numeric_limits<RT>::max();

                RT u_new = RT(0.0);
                bool step_success = false;
                while (!step_success) {

                    if (isTimestepTooSmall(dt, tau)) {
                        break;
                    }
                    RT mu = RT(1.0) / dt;

                    converged = true;
                    res_norm_a = RT(0.0);
                    res_norm_r = RT(0.0);

                    RT u1 = a_u;
                    {
                        RT rhs = mu * a_u;
                        RT res_norm_a_i = std::numeric_limits<RT>::max();
                        RT res_norm_r_i = std::numeric_limits<RT>::max();
                        bool converged_i = false;
                        m_newton( u1, rhs, mu,
                                  m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s, ventFactor(a_u),
                                  res_norm_a_i, res_norm_r_i, converged_i );
                        converged = converged && converged_i;
                        res_norm_a = std::max(res_norm_a, res_norm_a_i);
                        res_norm_r = std::max(res_norm_r, res_norm_r_i);
                    }
                    if (u1 <= 0) {
                        dt *= myhalf;
                        continue;
                    }
                    RT f1 = evalRHS(u1);

                    RT u2 = u1;
                    {
                        RT rhs = mu * (a_u - dt*f1);
                        RT res_norm_a_i = std::numeric_limits<RT>::max();
                        RT res_norm_r_i = std::numeric_limits<RT>::max();
                        bool converged_i = false;
                        m_newton( u2, rhs, mu,
                                  m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s, ventFactor(a_u),
                                  res_norm_a_i, res_norm_r_i, converged_i );
                        converged = converged && converged_i;
                        res_norm_a = std::max(res_norm_a, res_norm_a_i);
                        res_norm_r = std::max(res_norm_r, res_norm_r_i);
                    }
                    if (u2 <= 0) {
                        dt *= myhalf;
                        continue;
                    }
                    RT f2 = evalRHS(u2);

                    u_new = a_u + myhalf * dt * (f1 + f2);

                    if (converged) {
                        if (std::isfinite(u_new)) {
                            if (u_new > 0) {
                                step_success = true;
                                break;
                            }
                        }
                    }
                    dt *= myhalf;
                }

                if (step_success) {

                    RT snorm = std::sqrt((a_u-u_new)*(a_u-u_new)/(a_u*a_u));
                    a_u = u_new;
                    cur_time += dt;

                    printStepInfoNewton(cur_time, dt, tau, std::sqrt(a_u), snorm,
                                        res_norm_a, res_norm_r, converged);
                    if (snorm < m_stol) {
                        break;
                    }

                } else {

                    a_success = false;
                    break;

                }

                n_step++;
            }

            return;
        }

    };
}

/*! \brief Namespace with classes and functions for deposition/sublimation
 *  (solid <--> vapour) */
namespace SDMassChangeUtils_SV
{
    /*! \brief Phase change equation (in terms of mass) */
    template <typename RT /*!< real-type */ >
    struct dMdt
    {
        RT L;       /*!< latent heat of vaporization (condensate) */
        RT K;       /*!< thermal conductivity */
        RT Rv;      /*!< gas constant of air with vapour */
        RT rho_ice; /*!< True density of ice */

        /*! \brief Compute growth ratio */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT growthRatio ( const RT a_T /*!< temperature */ ) const noexcept
        {
            // Table of the inherent growth ratio (Chen and Lamb 1994)
            // T(n) = (1.375-n)/4. Here, T is the temperature in [degreeC], and n is the index.
            static constexpr amrex::Real tb_habit[121] = {
                amrex::Real(1.000000e+00), amrex::Real(9.799412e-01), amrex::Real(9.600636e-01), amrex::Real(9.399397e-01), amrex::Real(9.200258e-01),
                amrex::Real(9.001191e-01), amrex::Real(8.800351e-01), amrex::Real(8.570378e-01), amrex::Real(8.340652e-01), amrex::Real(8.109611e-01),
                amrex::Real(7.830689e-01), amrex::Real(7.550922e-01), amrex::Real(7.030723e-01), amrex::Real(5.370318e-01), amrex::Real(4.677351e-01),
                amrex::Real(5.248075e-01), amrex::Real(6.309573e-01), amrex::Real(8.128305e-01), amrex::Real(1.096478e+00), amrex::Real(1.479108e+00),
                amrex::Real(1.905461e+00), amrex::Real(2.089296e+00), amrex::Real(2.290868e+00), amrex::Real(2.398833e+00), amrex::Real(2.454709e+00),
                amrex::Real(2.426610e+00), amrex::Real(2.371374e+00), amrex::Real(2.290868e+00), amrex::Real(2.137962e+00), amrex::Real(1.995262e+00),
                amrex::Real(1.862087e+00), amrex::Real(1.737801e+00), amrex::Real(1.621810e+00), amrex::Real(1.513561e+00), amrex::Real(1.396368e+00),
                amrex::Real(1.288250e+00), amrex::Real(1.188502e+00), amrex::Real(1.096478e+00), amrex::Real(1.000000e+00), amrex::Real(9.225714e-01),
                amrex::Real(8.511380e-01), amrex::Real(7.852356e-01), amrex::Real(7.244360e-01), amrex::Real(6.683439e-01), amrex::Real(6.165950e-01),
                amrex::Real(5.754399e-01), amrex::Real(5.370318e-01), amrex::Real(5.011872e-01), amrex::Real(4.677351e-01), amrex::Real(4.365158e-01),
                amrex::Real(4.073803e-01), amrex::Real(3.801894e-01), amrex::Real(3.548134e-01), amrex::Real(3.311311e-01), amrex::Real(3.162278e-01),
                amrex::Real(3.019952e-01), amrex::Real(2.917427e-01), amrex::Real(2.851018e-01), amrex::Real(2.818383e-01), amrex::Real(2.786121e-01),
                amrex::Real(2.754229e-01), amrex::Real(2.786121e-01), amrex::Real(2.818383e-01), amrex::Real(2.851018e-01), amrex::Real(2.917427e-01),
                amrex::Real(2.985383e-01), amrex::Real(3.090295e-01), amrex::Real(3.198895e-01), amrex::Real(3.311311e-01), amrex::Real(3.467369e-01),
                amrex::Real(3.672823e-01), amrex::Real(3.935501e-01), amrex::Real(4.265795e-01), amrex::Real(4.570882e-01), amrex::Real(4.897788e-01),
                amrex::Real(5.248075e-01), amrex::Real(5.623413e-01), amrex::Real(6.095369e-01), amrex::Real(6.606934e-01), amrex::Real(7.161434e-01),
                amrex::Real(7.852356e-01), amrex::Real(8.609938e-01), amrex::Real(9.549926e-01), amrex::Real(1.047129e+00), amrex::Real(1.148154e+00),
                amrex::Real(1.258925e+00), amrex::Real(1.380384e+00), amrex::Real(1.496236e+00), amrex::Real(1.603245e+00), amrex::Real(1.698244e+00),
                amrex::Real(1.778279e+00), amrex::Real(1.840772e+00), amrex::Real(1.883649e+00), amrex::Real(1.905461e+00), amrex::Real(1.905461e+00),
                amrex::Real(1.883649e+00), amrex::Real(1.862087e+00), amrex::Real(1.840772e+00), amrex::Real(1.798871e+00), amrex::Real(1.737801e+00),
                amrex::Real(1.698244e+00), amrex::Real(1.640590e+00), amrex::Real(1.584893e+00), amrex::Real(1.549173e+00), amrex::Real(1.513910e+00),
                amrex::Real(1.476046e+00), amrex::Real(1.452112e+00), amrex::Real(1.428894e+00), amrex::Real(1.412863e+00), amrex::Real(1.393157e+00),
                amrex::Real(1.376892e+00), amrex::Real(1.361131e+00), amrex::Real(1.348963e+00), amrex::Real(1.336903e+00), amrex::Real(1.327089e+00),
                amrex::Real(1.317953e+00), amrex::Real(1.308881e+00), amrex::Real(1.302867e+00), amrex::Real(1.293898e+00), amrex::Real(1.287953e+00),
                amrex::Real(1.279087e+00)
            };
            auto T_degC = a_T - RT(273.15); // [K] => [degC]
            int n = std::min(std::max(static_cast<int>(std::round(RT(1.375)-RT(4.0)*T_degC)),1),121);
            return static_cast<RT>(tb_habit[n-1]);

        }

        /*! \brief Compute F_k + F_d */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT Fk_plus_Fd ( const RT a_T,    /*!< temperature */
                        const RT a_e_s,  /*!< saturation pressure */
                        const RT a_D     /*!< mol. diffusion coeff */ ) const noexcept
        {
            auto F_k = (L/(Rv*a_T) - RT(1.0)) * (L / (K*a_T));
            auto F_d = (Rv*a_T) / (a_D*a_e_s);
            return (F_k + F_d);
        }

        /*! \brief Compute dynamic viscosity */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT viscCoeff ( const RT a_T /*!< temperature */ ) const noexcept
        {
            auto T_degC = a_T - RT(273.15); // [K] => [degC]
            RT visc_coeff = RT(0.0);
            if( T_degC >= RT(0.0) ) {
                visc_coeff = ( RT(1.7180) + RT(4.9E-3)*T_degC ) * RT(1.E-5);
            } else {
                visc_coeff = ( RT(1.7180) + RT(4.9E-3)*T_degC -RT(1.2E-5)*T_degC*T_degC ) * RT(1.E-5);
            }
            return visc_coeff;
        }

        /*! \brief Compute capacitance */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT capacitance( const RT a_a, /*!< equatorial radius */
                        const RT a_c  /*!< polar radius */ ) const noexcept
        {
            auto phi = a_c / a_a;
            RT capaci = RT(0.0);
            if(phi >= RT(1.001)) {
                auto ecent = std::sqrt(RT(1.0) - RT(1.0)/(phi*phi));
                capaci = ecent*a_c / std::log((RT(1.0)+ecent)*phi);
            } else if (phi <= RT(0.999)) {
                auto ecent = std::sqrt(RT(1.0) - phi*phi);
                capaci = ecent*a_a / std::asin(ecent);
            } else {
                capaci = (RT(2.0)*a_a + a_c)/three;
            }
            return capaci;
        }

        /*! \brief Compute ventilation effect */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT ventX ( const RT a_a,    /*!< equatorial radius */
                   const RT a_c,    /*!< polar radius */
                   const RT a_vterm,/*!< terminal velocity */
                   const RT a_T,    /*!< temperature */
                   const RT a_rhom, /*!< moist density */
                   const RT a_D     /*!< mol. diffusion coeff */ ) const noexcept
        {
            auto visc_coeff = viscCoeff(a_T);
            auto maxdim = RT(2.0) * std::max(a_a, a_c);
            auto n_re = a_rhom * maxdim * a_vterm / visc_coeff;
            auto n_sc = visc_coeff / (a_rhom*a_D);
            auto X_vent = std::cbrt(n_sc) * std::sqrt(n_re);
            return X_vent;
        }

        /*! \brief Radius change */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT growthRatioStar ( const RT a_dmass,   /*!< delta mass */
                             const RT a_a,       /*!< equatorial radius */
                             const RT a_c,       /*!< polar radius */
                             const RT a_vterm,   /*!< terminal velocity */
                             const RT a_rhom,    /*!< moist air density */
                             const RT a_T,       /*!< temperature */
                             const RT a_D        /*!< mol. diffusion coeff */ ) const noexcept
        {
            auto capaci = capacitance(a_a, a_c);
            auto X_vent = ventX(a_a, a_c, a_vterm, a_T, a_rhom, a_D);
            RT fratio = RT(0.0);
            if(X_vent <= RT(1.0)) {
               fratio =    (RT(1.0) + RT(0.14)*X_vent*X_vent*std::sqrt(a_c/capaci))
                         / (RT(1.0) + RT(0.14)*X_vent*X_vent*std::sqrt(a_a/capaci));
            } else {
               fratio =   (RT(0.86) + RT(0.28)*X_vent*std::sqrt(a_c/capaci))
                        / (RT(0.86) + RT(0.28)*X_vent*std::sqrt(a_a/capaci));
            }

            auto gr = growthRatio(a_T);
            if ((RT(2.0)*std::max(a_a,a_c) < RT(1.0e-5)) || (a_dmass <= RT(0.0))) { gr = RT(1.0); }
            auto gr_star = gr * fratio;
            if((a_dmass > RT(0.0)) && (a_c/a_a > RT(40.0))) { gr_star = std::min(gr_star,RT(1.0)); }

            return gr_star;
        }

        /*! \brief Volume change */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT dVolume ( const RT a_dmass,        /*!< delta mass */
                     const RT a_a,            /*!< equatorial radius */
                     const RT a_c,            /*!< polar radius */
                     const RT a_rho,          /*!< ice density of particle */
                     const RT a_S,            /*!< saturation ratio */
                     const RT a_T,            /*!< temperature */
                     const RT a_e_s,          /*!< saturation pressure */
                     const RT a_es_ratio_wi,  /*!< ratio of saturation pressures of water to ice */
                     const RT a_D             /*!< mol. diffusion coeff */ ) const noexcept
        {
            RT rho_i_ds = a_rho;
            auto gr = growthRatio(a_T);
            if (a_dmass > RT(0.0)) {
                if ((gr < RT(1.0)) && (a_a < RT(1.0e-4))) {
                    rho_i_ds = rho_ice;
                } else {
                    if (RT(2.0)*std::max(a_a,a_c) < RT(1.0e-5)) { gr = RT(1.0); }
                    auto ex_vap_dens = RT(1000.0) * (std::min(a_S,a_es_ratio_wi) - RT(1.0)) / (Fk_plus_Fd(a_T,a_e_s,a_D) * a_D);
                    rho_i_ds = rho_ice * std::exp(-three*std::max(ex_vap_dens-RT(0.05),RT(0.0))/gr);
                }
            }

            auto dvol = a_dmass / rho_i_ds;
            return dvol;
        }

        /*! \brief Radius change */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT dLogRadius ( const RT a_grstar,  /*!< growth ratio */
                        const RT a_dlogvol  /*!< delta log volume */ ) const noexcept
        {
            return a_dlogvol / (a_grstar+RT(2.0));
        }

        /*! \brief Right-hand-side of the phase change ODE */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT rhs_func ( const RT a_a,    /*!< equatorial radius */
                      const RT a_c,    /*!< polar radius */
                      const RT a_vterm,/*!< terminal velocity */
                      const RT a_S,    /*!< saturation ratio */
                      const RT a_rhom, /*!< moist density */
                      const RT a_T,    /*!< temperature */
                      const RT a_e_s,  /*!< saturation pressure */
                      const RT a_D     /*!< mol. diffusion coeff */ ) const noexcept
        {
            auto capaci = capacitance(a_a, a_c);
            auto X_vent = ventX(a_a, a_c, a_vterm, a_T, a_rhom, a_D);

            RT f_vent = RT(0.0);
            if (X_vent <= RT(1.0)) { f_vent = RT(1.0) + RT(0.14)*X_vent*X_vent; }
            else               { f_vent = RT(0.86) + RT(0.28)*X_vent; }

            auto retval = RT(4.0)*PI * capaci * (a_S-RT(1.0)) * f_vent / Fk_plus_Fd(a_T, a_e_s, a_D);
            return retval;
        }

        /*! \brief Ice melt rate dm_ice/dt (Pruppacher-Klett; Seifert-Beheng 2006 Eq. 72).
         *  Negative when melting. Heat conducted from warm air plus latent heat of the
         *  vapour flux to the 0 degC wet surface, divided by the latent heat of fusion.
         *  Reuses capacitance and ventilation. For this functor built with the water
         *  material, L is the latent heat of vaporization, K the air thermal
         *  conductivity, and Rv the vapour gas constant. */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT meltRate ( const RT a_a,         /*!< equatorial radius */
                      const RT a_c,         /*!< polar radius */
                      const RT a_vterm,     /*!< terminal velocity */
                      const RT a_T,         /*!< temperature */
                      const RT a_e_inf,     /*!< ambient vapour pressure */
                      const RT a_e_sat_wT0, /*!< sat. pressure over water at melting point */
                      const RT a_rhom,      /*!< moist density */
                      const RT a_D,         /*!< vapour diffusion coeff */
                      const RT a_L_f        /*!< latent heat of fusion */ ) const noexcept
        {
            auto capaci = capacitance(a_a, a_c);
            auto X_vent = ventX(a_a, a_c, a_vterm, a_T, a_rhom, a_D);

            RT f_vent = RT(0.0);
            if (X_vent <= RT(1.0)) { f_vent = RT(1.0) + RT(0.14)*X_vent*X_vent; }
            else               { f_vent = RT(0.86) + RT(0.28)*X_vent; }

            RT T0 = RT(tmelt);
            RT heat =   K * f_vent * (a_T - T0)
                      + (L*a_D/Rv) * (a_e_inf/a_T - a_e_sat_wT0/T0) * f_vent;
            return -(RT(4.0)*PI*capaci/a_L_f) * heat;
        }

    };

    /*! \brief Implicit and explicit time integrators for the phase change equation */
    template<   typename ODE /*!< ODE */,
                typename RT /*!< real-type */ >
    struct TI
    {
        const ODE m_ode; /*!< ODE */

        RT m_t_final;   /*!< final time */
        mutable RT m_dt;/*!< dt */
        RT m_max_steps; /*!< max number of timesteps */

        RT m_a;         /*!< equatorial radius */
        RT m_c;         /*!< polar radius */
        RT m_vt;        /*!< terminal velocity */
        RT m_S;         /*!< saturation ratio */
        RT m_rhom;      /*!< moist density */
        RT m_T;         /*!< temperature */
        RT m_e_s;       /*!< saturation pressure */
        RT m_D;         /*!< mol. diff. coeff */

        bool m_verbose; /*!< verbosity */

        /*! \brief 1st-order, forward Euler method */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void fe ( RT& a_u,  /*!< solution */
                  bool& a_success  /*!< success/failure flag */ ) const
        {
            RT cur_time = RT(0.0);
            a_success = true;

            int n_step = 0;
            while ((cur_time < m_t_final) && (n_step < m_max_steps)) {

                if ((cur_time + m_dt) > m_t_final) {
                    m_dt = m_t_final - cur_time;
                }

                RT rhs = m_ode.rhs_func(m_a, m_c, m_vt, m_S, m_rhom, m_T, m_e_s, m_D);
                RT u_new = a_u + m_dt*rhs;
                if (!std::isfinite(u_new)) {
                    a_success = false;
                    break;
                }

                RT snorm = std::sqrt((a_u-u_new)*(a_u-u_new)/(a_u*a_u));
                a_u = u_new;
                cur_time += m_dt;

                if (m_verbose) {
                    printf( "Time %1.2e, m_dt = %1.2e, mass = %1.4e, snorm = %1.1e\n",
                            cur_time, m_dt, a_u, snorm);
                }

                n_step++;
            }

            return;
        }

    };
}

#endif
#endif
