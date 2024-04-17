#ifndef SUPERDROPLET_PC_COALESCENCE_H_
#define SUPERDROPLET_PC_COALESCENCE_H_

#include <cmath>

#ifdef ERF_USE_PARTICLES

namespace SDMassChangeUtils
{
    /*! \brief RHS function for super-droplet radius ODE */
    template<typename RT>
    struct dRsqdt_RHSFunc
    {
        RT a;
        RT b;
        RT L;
        RT K;
        RT Rv;
        RT rho_l;
        RT D;

        AMREX_GPU_HOST_DEVICE
        RT operator() ( const RT a_R_sq,
                        const RT a_S,
                        const RT a_T,
                        const RT a_e_s,
                        const RT a_M_s ) const noexcept
        {
            RT F_k = ( L/(Rv*a_T) - 1.0) * ((L*rho_l) / (K*a_T));
            RT F_d = (rho_l*Rv*a_T) / (D*a_e_s);

            RT alpha = 2.0 * (a_S-1.0) / (F_k + F_d);
            RT beta = -2.0 * (a/a_T) / (F_k + F_d);
            RT gamma = 2.0 * b * a_M_s / (F_k + F_d);

            RT R_inv = 1.0/std::sqrt(a_R_sq);
            RT R_inv_cubed = R_inv*R_inv*R_inv;

            RT retval = alpha + beta*R_inv + gamma*R_inv_cubed;
            return retval;
        }

    };

    /*! \brief Jacobian for super-droplet radius ODE */
    template<typename RT>
    struct dRsqdt_RHSJac
    {
        RT a;
        RT b;
        RT L;
        RT K;
        RT Rv;
        RT rho_l;
        RT D;

        AMREX_GPU_HOST_DEVICE
        RT operator() ( const RT a_R_sq,
                        const RT a_T,
                        const RT a_e_s,
                        const RT a_M_s ) const noexcept
        {
            RT F_k = ( L/(Rv*a_T) - 1.0) * ((L*rho_l) / (K*a_T));
            RT F_d = (rho_l*Rv*a_T) / (D*a_e_s);

            RT beta = -2.0 * (a/a_T) / (F_k + F_d);
            RT gamma = 2.0 * b * a_M_s / (F_k + F_d);

            RT R_inv = 1.0/std::sqrt(a_R_sq);
            RT R_inv_3 = R_inv*R_inv*R_inv;
            RT R_inv_5 = R_inv_3*R_inv*R_inv;

            RT retval = -0.5 * (beta*R_inv_3 + 3.0*gamma*R_inv_5);
            return retval;
        }

    };

    /*! \brief Newton solver for super-droplet radius ODE */
    template<typename RHSFunc, typename JacFunc, typename RT>
    struct NewtonSolver
    {
        RHSFunc m_rhs;
        JacFunc m_jac;

        RT  m_rtol;
        RT  m_atol;
        RT  m_stol;
        int m_maxits;

        AMREX_GPU_HOST_DEVICE
        void operator()  (  RT&             a_u,
                            RT&             a_r,
                            const RT&       a_dt,
                            const RT&       a_S,
                            const RT&       a_T,
                            const RT&       a_e_s,
                            const RT&       a_M_s,
                            RT&             a_res_norm_a,
                            RT&             a_res_norm_r,
                            bool&           a_converged ) const
        {
            a_converged = false;
            RT res_norm0 = 0.0;

            for (int k = 0; k < m_maxits; k++) {
                RT residual = a_u - (a_r + a_dt * m_rhs( a_u, a_S, a_T, a_e_s, a_M_s ));
                a_res_norm_a = std::sqrt(residual*residual);

                if (k == 0) { res_norm0 = a_res_norm_a; }
                a_res_norm_r = a_res_norm_a / res_norm0;

                if (a_res_norm_a <= m_atol) {
                    a_converged = true;
                    break;
                }
                if (a_res_norm_r <= m_rtol) {
                    a_converged = true;
                    break;
                }

                RT slope = 1.0 - a_dt * m_jac( a_u, a_T, a_e_s, a_M_s );
                RT du = 0.0;
                if ((slope*slope) < (1.0e-12*residual*residual)) {
                    du = 0.0;
                } else {
                    du = - residual / slope;
                }

                RT du_norm = std::sqrt(du*du);
                if (du_norm <= m_stol) {
                    a_converged = true;
                    break;
                }

                if (a_u + du < 0) {
                    a_u *= 0.99;
                } else {
                    a_u += du;
                }
            }
        }
    };
}

#endif
#endif
