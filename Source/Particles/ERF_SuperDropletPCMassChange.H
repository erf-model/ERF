#ifndef SUPERDROPLET_PC_MASSCHANGE_H_
#define SUPERDROPLET_PC_MASSCHANGE_H_

#include <cmath>
#include "ERF_Constants.H"

#ifdef ERF_USE_PARTICLES

/*! \brief Namespace with classes and functions for condensation/evaporation
 *  (liquid <--> vapour) */
namespace SDMassChangeUtils_LV
{
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
                         const RT a_N_s  /*!< total solute moles */ ) const noexcept
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

            return retval;
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
                      const RT a_N_s   /*!< total solute moles */ ) const noexcept
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

            return retval;
        }

        /*! \brief Jacobian of right-hand-side of the phase change ODE */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT rhs_jac ( const RT a_R_sq, /*!< radius squared */
                     const RT a_T,    /*!< temperature */
                     const RT a_e_s,  /*!< saturation pressure */
                     const RT a_D,    /*!< mol. diffusion coeff */
                     const RT a_a,    /*!< curvature coeff */
                     const RT a_b,    /*!< solute coeff */
                     const RT a_N_s   /*!< total solute moles */ ) const noexcept
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

            return retval;
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
                                  + m_ne.rhs_func( a_u, a_S, a_T, a_e_s, a_D, a_a, a_b, a_N_s ) );
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

                RT slope = a_mu - m_ne.rhs_jac( a_u, a_T, a_e_s, a_D, a_a, a_b, a_N_s );
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

        /*! \brief Compute timestep from stiffness estimate */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT computeTimestep (const RT& a_u) const
        {
            RT tau = m_ode.rhs_jac(a_u, m_T, m_e_s, m_D, m_a, m_b, m_N_s);
            return m_cfl / std::sqrt(tau*tau);
        }

        /*! \brief Compute stiffness estimate (tau) */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RT computeTau (const RT& a_u) const
        {
            return m_ode.rhs_jac(a_u, m_T, m_e_s, m_D, m_a, m_b, m_N_s);
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
            return m_ode.rhs_func(a_u, m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s);
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

            RT tau = m_ode.rhs_jac(a_u, m_T, m_e_s, m_D, m_a, m_b, m_N_s);
            RT dt = m_cfl / std::sqrt(tau*tau);

            RT dt_new = dt;
            RT a_u_old = a_u;

            int n_step = 0;
            while ((cur_time < m_t_final) && (n_step < m_max_steps)) {

                if (!m_adapt_dt) {
                    tau = m_ode.rhs_jac(a_u, m_T, m_e_s, m_D, m_a, m_b, m_N_s);
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

                RT u_new = zero;
                bool step_success = false;
                while (!step_success) {

                    if (    (dt < (amrex::Real(1.0e-12)*m_cfl/std::sqrt(tau*tau)))
                         && (dt < (amrex::Real(1.0e-12)*m_t_final)) ) {
                        break;
                    }

                    RT u1 = a_u;
                    if (u1 <= 0) {
                        dt *= myhalf;
                        continue;
                    }
                    RT f1 = m_ode.rhs_func(u1, m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s);

                    RT u2 = a_u + myhalf*dt*f1;
                    if (u2 <= 0) {
                        dt *= myhalf;
                        continue;
                    }
                    RT f2 = m_ode.rhs_func(u2, m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s);

                    RT u3 = a_u + amrex::Real(0.75)*dt*f2;
                    if (u3 <= 0) {
                        dt *= myhalf;
                        continue;
                    }
                    RT f3 = m_ode.rhs_func(u3, m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s);

                    RT u4 = a_u + (one/amrex::Real(9.0))*dt * (two*f1 + three*f2 + amrex::Real(4.0)*f3);
                    if (u4 <= 0) {
                        dt *= myhalf;
                        continue;
                    }
                    RT f4 = m_ode.rhs_func(u4, m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s);

                    u_new = u4;

                    if (m_adapt_dt) {
                        RT u_embed = a_u + (one/amrex::Real(24.0))*dt * (amrex::Real(7.0)*f1 + amrex::Real(6.0)*f2 + amrex::Real(8.0)*f3 + three*f4);
                        RT err = std::sqrt((u_new-u_embed)*(u_new-u_embed));
                        RT tol = m_atol + m_rtol * std::max(a_u, a_u_old);
                        RT E = err / tol;
                        dt_new = dt * std::exp((one/3)*std::log(one/E));
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
            RT cur_time = zero;
            a_success = true;

            int n_step = 0;
            while ((cur_time < m_t_final) && (n_step < m_max_steps)) {

                RT tau = m_ode.rhs_jac(a_u, m_T, m_e_s, m_D, m_a, m_b, m_N_s);
                RT dt = m_cfl / std::sqrt(tau*tau);
                if ((cur_time + dt) > m_t_final) {
                    dt = m_t_final - cur_time;
                }
                if (!std::isfinite(dt)) {
                    a_success = false;
                    break;
                }

                RT u_new = zero;
                bool step_success = false;
                while (!step_success) {

                    if (    (dt < (amrex::Real(1.0e-12)*m_cfl/std::sqrt(tau*tau)))
                         && (dt < (amrex::Real(1.0e-12)*m_t_final)) ) {
                        break;
                    }

                    RT u1 = a_u;
                    if (u1 <= 0) {
                        dt *= myhalf;
                        continue;
                    }
                    RT f1 = m_ode.rhs_func(u1, m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s);

                    RT u2 = a_u + myhalf*dt*f1;
                    if (u2 <= 0) {
                        dt *= myhalf;
                        continue;
                    }
                    RT f2 = m_ode.rhs_func(u2, m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s);

                    RT u3 = a_u + myhalf*dt*f2;
                    if (u3 <= 0) {
                        dt *= myhalf;
                        continue;
                    }
                    RT f3 = m_ode.rhs_func(u3, m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s);

                    RT u4 = a_u + one*dt*f3;
                    if (u4 <= 0) {
                        dt *= myhalf;
                        continue;
                    }
                    RT f4 = m_ode.rhs_func(u4, m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s);

                    u_new = a_u + dt*(f1+two*f2+two*f3+f4)/amrex::Real(6.0);

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

        /*! \brief 1st order implicit backward Euler method */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void be ( RT& a_u,  /*!< solution */
                  bool& a_success  /*!< success/failure flag */ ) const
        {
            RT cur_time = zero;
            a_success = true;

            int n_step = 0;
            while ((cur_time < m_t_final) && (n_step < m_max_steps)) {

                RT tau = m_ode.rhs_jac(a_u, m_T, m_e_s, m_D, m_a, m_b, m_N_s);
                RT dt = m_cfl / std::sqrt(tau*tau);
                if ((cur_time + dt) > m_t_final) {
                    dt = m_t_final - cur_time;
                }
                if (!std::isfinite(dt)) {
                    a_success = false;
                    break;
                }

                RT res_norm_a = DBL_MAX;
                RT res_norm_r = DBL_MAX;
                bool converged = false;

                RT u_new = zero;
                bool step_success = false;
                while (!step_success) {

                    if (    (dt < (amrex::Real(1.0e-12)*m_cfl/std::sqrt(tau*tau)))
                         && (dt < (amrex::Real(1.0e-12)*m_t_final)) ) {
                        break;
                    }

                    RT mu = one / dt;
                    RT rhs = mu * a_u;
                    u_new = a_u;
                    m_newton( u_new, rhs, mu,
                              m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s,
                              res_norm_a, res_norm_r, converged );

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

                    if (m_verbose) {
                        AMREX_DEVICE_PRINTF( "Time %1.2e, dt = %1.2e, cfl = %1.1e, radius = %1.4e, snorm = %1.1e\n",
                                cur_time, dt, dt * std::sqrt(tau*tau), std::sqrt(a_u), snorm);
                        AMREX_DEVICE_PRINTF( "    norms = %1.3e (abs), %1.3e (rel), converged = %s\n",
                                res_norm_a, res_norm_r,
                                (converged ? "yes" : "no") );
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

        /*! \brief 2nd-order implicit Crank-Nicolson method */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void cn ( RT& a_u,  /*!< solution */
                  bool& a_success  /*!< success/failure flag */ ) const
        {
            RT cur_time = zero;
            a_success = true;

            int n_step = 0;
            while ((cur_time < m_t_final) && (n_step < m_max_steps)) {

                RT tau = m_ode.rhs_jac(a_u, m_T, m_e_s, m_D, m_a, m_b, m_N_s);
                RT dt = m_cfl / std::sqrt(tau*tau);
                if ((cur_time + dt) > m_t_final) {
                    dt = m_t_final - cur_time;
                }
                if (!std::isfinite(dt)) {
                    a_success = false;
                    break;
                }

                RT res_norm_a = DBL_MAX;
                RT res_norm_r = DBL_MAX;
                bool converged = false;

                RT u_new = zero;
                bool step_success = false;
                while (!step_success) {

                    if (    (dt < (amrex::Real(1.0e-12)*m_cfl/std::sqrt(tau*tau)))
                         && (dt < (amrex::Real(1.0e-12)*m_t_final)) ) {
                        break;
                    }
                    RT mu = one / (myhalf*dt);

                    RT u1 = a_u;
                    RT f1 = m_ode.rhs_func(u1, m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s);

                    RT u2 = u1;
                    RT rhs = mu * (a_u + myhalf*dt*f1);
                    m_newton( u2, rhs, mu,
                              m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s,
                              res_norm_a, res_norm_r, converged );
                    if (u2 <= 0) {
                        dt *= myhalf;
                        continue;
                    }
                    RT f2 = m_ode.rhs_func(u2, m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s);

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

                    if (m_verbose) {
                        AMREX_DEVICE_PRINTF( "Time %1.2e, dt = %1.2e, cfl = %1.1e, radius = %1.4e, snorm = %1.1e\n",
                                cur_time, dt, dt * std::sqrt(tau*tau), std::sqrt(a_u), snorm);
                        AMREX_DEVICE_PRINTF( "    norms = %1.3e (abs), %1.3e (rel), converged = %s\n",
                                res_norm_a, res_norm_r,
                                (converged ? "yes" : "no") );
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

        /*! \brief 2nd-order, 2-stage diagonally-implicit RK method */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void dirk212 ( RT& a_u,  /*!< solution */
                       bool& a_success  /*!< success/failure flag */ ) const
        {
            RT cur_time = zero;
            a_success = true;

            int n_step = 0;
            while ((cur_time < m_t_final) && (n_step < m_max_steps)) {

                RT tau = m_ode.rhs_jac(a_u, m_T, m_e_s, m_D, m_a, m_b, m_N_s);
                RT dt = m_cfl / std::sqrt(tau*tau);
                if ((cur_time + dt) > m_t_final) {
                    dt = m_t_final - cur_time;
                }
                if (!std::isfinite(dt)) {
                    a_success = false;
                    break;
                }

                bool converged = false;
                RT res_norm_a = DBL_MAX;
                RT res_norm_r = DBL_MAX;

                RT u_new = zero;
                bool step_success = false;
                while (!step_success) {

                    if (    (dt < (amrex::Real(1.0e-12)*m_cfl/std::sqrt(tau*tau)))
                         && (dt < (amrex::Real(1.0e-12)*m_t_final)) ) {
                        break;
                    }
                    RT mu = one / dt;

                    converged = true;
                    res_norm_a = zero;
                    res_norm_r = zero;

                    RT u1 = a_u;
                    {
                        RT rhs = mu * a_u;
                        RT res_norm_a_i = DBL_MAX;
                        RT res_norm_r_i = DBL_MAX;
                        bool converged_i = false;
                        m_newton( u1, rhs, mu,
                                  m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s,
                                  res_norm_a_i, res_norm_r_i, converged_i );
                        converged = converged && converged_i;
                        res_norm_a = std::max(res_norm_a, res_norm_a_i);
                        res_norm_r = std::max(res_norm_r, res_norm_r_i);
                    }
                    if (u1 <= 0) {
                        dt *= myhalf;
                        continue;
                    }
                    RT f1 = m_ode.rhs_func(u1, m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s);

                    RT u2 = u1;
                    {
                        RT rhs = mu * (a_u - dt*f1);
                        RT res_norm_a_i = DBL_MAX;
                        RT res_norm_r_i = DBL_MAX;
                        bool converged_i = false;
                        m_newton( u2, rhs, mu,
                                  m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s,
                                  res_norm_a_i, res_norm_r_i, converged_i );
                        converged = converged && converged_i;
                        res_norm_a = std::max(res_norm_a, res_norm_a_i);
                        res_norm_r = std::max(res_norm_r, res_norm_r_i);
                    }
                    if (u2 <= 0) {
                        dt *= myhalf;
                        continue;
                    }
                    RT f2 = m_ode.rhs_func(u2, m_S, m_T, m_e_s, m_D, m_a, m_b, m_N_s);

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

                    if (m_verbose) {
                        AMREX_DEVICE_PRINTF( "Time %1.2e, dt = %1.2e, cfl = %1.1e, radius = %1.4e, snorm = %1.1e\n",
                                cur_time, dt, dt * std::sqrt(tau*tau), std::sqrt(a_u), snorm);
                        AMREX_DEVICE_PRINTF( "    norms = %1.3e (abs), %1.3e (rel), converged = %s\n",
                                res_norm_a, res_norm_r,
                                (converged ? "yes" : "no") );
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

    };
}

#endif
#endif
