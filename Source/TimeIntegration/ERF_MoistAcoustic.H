#ifndef ERF_MOIST_ACOUSTIC_H_
#define ERF_MOIST_ACOUSTIC_H_

#include <AMReX.H>

#include <ERF_Constants.H>

namespace ERFMoistAcoustic {

enum {
    FastRho_comp = 0,
    FastRhoTheta_comp,
    FastRhoQv_comp,
    NumFastComps
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real theta_m (const amrex::Real theta_d,
                     const amrex::Real qv,
                     const bool use_moisture)
{
    return use_moisture ? theta_d * (amrex::Real(1) + (R_v / R_d) * qv)
                        : theta_d;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real rho_theta_m (const amrex::Real rho_theta_d,
                         const amrex::Real qv,
                         const bool use_moisture)
{
    return use_moisture ? rho_theta_d * (amrex::Real(1) + (R_v / R_d) * qv)
                        : rho_theta_d;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real rho_theta_m_pert (const amrex::Real rho_pert,
                              const amrex::Real rho_theta_d_pert,
                              const amrex::Real rho_qv_pert,
                              const amrex::Real theta_d_stage,
                              const amrex::Real qv_stage,
                              const bool use_moisture)
{
    if (!use_moisture) {
        return rho_theta_d_pert;
    }

    const amrex::Real rvr = R_v / R_d;
    return (amrex::Real(1) + rvr * qv_stage) * rho_theta_d_pert
         + rvr * theta_d_stage * rho_qv_pert
         - rvr * theta_d_stage * qv_stage * rho_pert;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real rho_theta_m_rhs (const amrex::Real rho_rhs,
                             const amrex::Real rho_theta_d_rhs,
                             const amrex::Real rho_qv_rhs,
                             const amrex::Real theta_d_stage,
                             const amrex::Real qv_stage,
                             const bool use_moisture)
{
    return rho_theta_m_pert(rho_rhs, rho_theta_d_rhs, rho_qv_rhs,
                            theta_d_stage, qv_stage, use_moisture);
}

} // namespace ERFMoistAcoustic

#endif
