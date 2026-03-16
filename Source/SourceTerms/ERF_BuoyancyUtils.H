#ifndef ERF_BUOYANCY_UTILS_H_
#define ERF_BUOYANCY_UTILS_H_

#include <ERF_EOS.H>
#include <ERF_Constants.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_dry_anelastic (int& i, int& j, int& k,
                        amrex::Real const& grav_gpu,
                        const amrex::Array4<const amrex::Real>& r0_arr,
                        const amrex::Array4<const amrex::Real>& th0_arr,
                        const amrex::Array4<const amrex::Real>& cell_data)
{
    // Note: this is the same term as the moist anelastic buoyancy when qv = qc = qt = 0
    amrex::Real theta_d_hi = cell_data(i,j,k  ,RhoTheta_comp)/cell_data(i,j,k  ,Rho_comp);
    amrex::Real theta_d_lo = cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp);

    amrex::Real  theta_d_wface = amrex::Real(0.5) * (theta_d_lo + theta_d_hi);
    amrex::Real theta_d0_wface = amrex::Real(0.5) * (th0_arr(i,j,k) + th0_arr(i,j,k-1));
    amrex::Real     rho0_wface = amrex::Real(0.5) * (r0_arr(i,j,k) + r0_arr(i,j,k-1));

    return (-rho0_wface * grav_gpu * (theta_d_wface - theta_d0_wface) / theta_d0_wface);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_dry_anelastic_eb (int& i, int& j, int& k,
                        amrex::Real const& grav_gpu,
                        const amrex::Array4<const amrex::Real>& r0_arr,
                        const amrex::Array4<const amrex::Real>& th0_arr,
                        const amrex::Array4<const amrex::Real>& cell_data,
                        amrex::Array4<amrex::EBCellFlag const> const& flag)
{
    amrex::Real buoyancy = 0.0;
    if (flag(i,j,k).isRegular()) {
        amrex::Real theta_d_hi = cell_data(i,j,k  ,RhoTheta_comp)/cell_data(i,j,k  ,Rho_comp);
        amrex::Real theta_d_lo = cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp);
        amrex::Real  theta_d_wface = amrex::Real(0.5) * (theta_d_lo + theta_d_hi);
        amrex::Real theta_d0_wface = amrex::Real(0.5) * (th0_arr(i,j,k) + th0_arr(i,j,k-1));
        amrex::Real     rho0_wface = amrex::Real(0.5) * (r0_arr(i,j,k) + r0_arr(i,j,k-1));
        buoyancy = -rho0_wface * grav_gpu * (theta_d_wface - theta_d0_wface) / theta_d0_wface;
    } else if (flag(i,j,k).isSingleValued()) {
        if (flag(i,j,k-1).isCovered()) {
            amrex::Real theta_d_hihi = cell_data(i,j,k+1,RhoTheta_comp)/cell_data(i,j,k+1,Rho_comp);
            amrex::Real theta_d_hi = cell_data(i,j,k  ,RhoTheta_comp)/cell_data(i,j,k  ,Rho_comp);
            amrex::Real theta_d_wface  = amrex::Real(0.5) * (3.0*theta_d_hi - theta_d_hihi);
            amrex::Real theta_d0_wface = amrex::Real(0.5) * (3.0*th0_arr(i,j,k) - th0_arr(i,j,k+1));
            amrex::Real     rho0_wface = amrex::Real(0.5) * (3.0*r0_arr(i,j,k) - r0_arr(i,j,k+1));
            buoyancy = -rho0_wface * grav_gpu * (theta_d_wface - theta_d0_wface) / theta_d0_wface;
        } else if (flag(i,j,k).isCovered()) {
            amrex::Real theta_d_lo = cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp);
            amrex::Real theta_d_lolo = cell_data(i,j,k-2,RhoTheta_comp)/cell_data(i,j,k-2,Rho_comp);
            amrex::Real  theta_d_wface = amrex::Real(0.5) * (3.0*theta_d_lo - theta_d_lolo);
            amrex::Real theta_d0_wface = amrex::Real(0.5) * (3.0*th0_arr(i,j,k-1) - th0_arr(i,j,k-2));
            amrex::Real     rho0_wface = amrex::Real(0.5) * (3.0*r0_arr(i,j,k-1) - r0_arr(i,j,k-2));
            buoyancy = -rho0_wface * grav_gpu * (theta_d_wface - theta_d0_wface) / theta_d0_wface;
        } else {
            amrex::Real theta_d_hi = cell_data(i,j,k  ,RhoTheta_comp)/cell_data(i,j,k  ,Rho_comp);
            amrex::Real theta_d_lo = cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp);
            amrex::Real  theta_d_wface = amrex::Real(0.5) * (theta_d_lo + theta_d_hi);
            amrex::Real theta_d0_wface = amrex::Real(0.5) * (th0_arr(i,j,k) + th0_arr(i,j,k-1));
            amrex::Real     rho0_wface = amrex::Real(0.5) * (r0_arr(i,j,k) + r0_arr(i,j,k-1));
            buoyancy = -rho0_wface * grav_gpu * (theta_d_wface - theta_d0_wface) / theta_d0_wface;
        }
    }
    return buoyancy;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_moist_anelastic (int& i, int& j, int& k,
                          amrex::Real const& grav_gpu,
                          amrex::Real const& rv_over_rd,
                          const amrex::Array4<const amrex::Real>&  r0_arr,
                          const amrex::Array4<const amrex::Real>& th0_arr,
                          const amrex::Array4<const amrex::Real>& qv0_arr,
                          const amrex::Array4<const amrex::Real>& cell_data,
                          const amrex::Array4<const amrex::Real>& qt_arr)
{
    amrex::Real Fact = (1.0 - rv_over_rd); // = -0.61

    amrex::Real th_d_hi = cell_data(i,j,k  ,RhoTheta_comp)/cell_data(i,j,k  ,Rho_comp);
    amrex::Real th_d_lo = cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp);

    amrex::Real qv_hi   = cell_data(i,j,k  ,RhoQ1_comp)   /cell_data(i,j,k  ,Rho_comp);
    amrex::Real qv_lo   = cell_data(i,j,k-1,RhoQ1_comp)   /cell_data(i,j,k-1,Rho_comp);

    amrex::Real qt_hi   = qt_arr(i,j,k  ) - qv_hi;
    amrex::Real th_v_hi = th_d_hi * (1.0 - Fact*qt_hi - rv_over_rd*qt_hi);

    amrex::Real qt_lo   = qt_arr(i,j,k-1) - qv_lo;
    amrex::Real th_v_lo = th_d_lo * (1.0 - Fact*qt_lo - rv_over_rd*qt_lo);

    amrex::Real th_v0_hi = th0_arr(i,j,k  ) * (1.0 - Fact*qv0_arr(i,j,k  ));
    amrex::Real th_v0_lo = th0_arr(i,j,k-1) * (1.0 - Fact*qv0_arr(i,j,k-1));

    amrex::Real qv0_hi      = qv0_arr(i,j,k  );
    amrex::Real qv0_lo      = qv0_arr(i,j,k-1);

    amrex::Real q_hi = (qv_hi-qv0_hi) - qt_hi + (th_v_hi - th_v0_hi) / th_v0_hi;
    amrex::Real q_lo = (qv_lo-qv0_lo) - qt_lo + (th_v_lo - th_v0_lo) / th_v0_lo;

    amrex::Real qavg = amrex::Real(0.5) * (q_hi + q_lo);

    amrex::Real r0avg = amrex::Real(0.5) * (r0_arr(i,j,k) + r0_arr(i,j,k-1));

    return (-r0avg * grav_gpu * qavg);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_rhopert (int& i, int& j, int& k,
                  amrex::Real const& grav_gpu,
                  const amrex::Array4<const amrex::Real>& r0_arr,
                  const amrex::Array4<const amrex::Real>& cell_data,
                  const amrex::Array4<const amrex::Real>& qt_arr)
{
    amrex::Real rhop_hi = cell_data(i,j,k  ,Rho_comp) * (1.0 + qt_arr(i,j,k  ))- r0_arr(i,j,k  );
    amrex::Real rhop_lo = cell_data(i,j,k-1,Rho_comp) * (1.0 + qt_arr(i,j,k-1))- r0_arr(i,j,k-1);
    return( grav_gpu * amrex::Real(0.5) * ( rhop_hi + rhop_lo ) );
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_rhopert_eb (int& i, int& j, int& k,
                    amrex::Real const& grav_gpu,
                    const amrex::Array4<const amrex::Real>& r0_arr,
                    const amrex::Array4<const amrex::Real>& cell_data,
                    const amrex::Array4<const amrex::Real>& qt_arr,
                    amrex::Array4<amrex::EBCellFlag const> const& flag)
{
    amrex::Real buoyancy = 0.0;
    if (flag(i,j,k).isRegular()) {
        amrex::Real rhop_hi = cell_data(i,j,k  ,Rho_comp) * (1.0 + qt_arr(i,j,k  ))- r0_arr(i,j,k  );
        amrex::Real rhop_lo = cell_data(i,j,k-1,Rho_comp) * (1.0 + qt_arr(i,j,k-1))- r0_arr(i,j,k-1);
        buoyancy = grav_gpu * amrex::Real(0.5) * ( rhop_hi + rhop_lo );
    } else if (flag(i,j,k).isSingleValued()) {
        if (flag(i,j,k-1).isCovered()) {
            amrex::Real rhop_hihi = cell_data(i,j,k+1,Rho_comp) * (1.0 + qt_arr(i,j,k+1))- r0_arr(i,j,k+1);
            amrex::Real rhop_hi   = cell_data(i,j,k  ,Rho_comp) * (1.0 + qt_arr(i,j,k  ))- r0_arr(i,j,k  );
            buoyancy = grav_gpu * 0.5 * ( 3.0 * rhop_hi - rhop_hihi );
            // buoyancy = grav_gpu * rhop_hi;
        } else if (flag(i,j,k).isCovered()) {
            amrex::Real rhop_lo   = cell_data(i,j,k-1,Rho_comp) * (1.0 + qt_arr(i,j,k-1))- r0_arr(i,j,k-1);
            amrex::Real rhop_lolo = cell_data(i,j,k-2,Rho_comp) * (1.0 + qt_arr(i,j,k-2))- r0_arr(i,j,k-2);
            buoyancy = grav_gpu * 0.5 * ( 3.0 * rhop_lo - rhop_lolo );
            // buoyancy = grav_gpu * rhop_lo;
        } else {
            amrex::Real rhop_hi = cell_data(i,j,k  ,Rho_comp) * (1.0 + qt_arr(i,j,k  ))- r0_arr(i,j,k  );
            amrex::Real rhop_lo = cell_data(i,j,k-1,Rho_comp) * (1.0 + qt_arr(i,j,k-1))- r0_arr(i,j,k-1);
            buoyancy = grav_gpu * amrex::Real(0.5) * ( rhop_hi + rhop_lo );
        }
    }
    return buoyancy;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_dry_Tpert (int& i, int& j, int& k,
                    amrex::Real const& grav_gpu,
                    amrex::Real const& rd_over_cp,
                    const amrex::Array4<const amrex::Real>& r0_arr,
                    const amrex::Array4<const amrex::Real>& p0_arr,
                    const amrex::Array4<const amrex::Real>& th0_arr,
                    const amrex::Array4<const amrex::Real>& cell_data)
{
    amrex::Real t0_hi = getTgivenPandTh(p0_arr(i,j,k), th0_arr(i,j,k), rd_over_cp);
    amrex::Real  t_hi = getTgivenRandRTh(cell_data(i,j,k  ,Rho_comp), cell_data(i,j,k  ,RhoTheta_comp));

    amrex::Real t0_lo = getTgivenPandTh(p0_arr(i,j,k-1), th0_arr(i,j,k-1), rd_over_cp);
    amrex::Real  t_lo = getTgivenRandRTh(cell_data(i,j,k-1,Rho_comp), cell_data(i,j,k-1,RhoTheta_comp));

    amrex::Real tprime_hi = (t_hi-t0_hi)/t0_hi;
    amrex::Real tprime_lo = (t_lo-t0_lo)/t0_lo;

    amrex::Real tp_avg = amrex::Real(0.5) * (tprime_hi + tprime_lo);
    amrex::Real r0_avg = amrex::Real(0.5) * (r0_arr(i,j,k) + r0_arr(i,j,k-1));

    return ( -r0_avg * grav_gpu * tp_avg);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_dry_Thpert (int& i, int& j, int& k,
                     amrex::Real const& grav_gpu,
                     const amrex::Array4<const amrex::Real>& r0_arr,
                     const amrex::Array4<const amrex::Real>& th0_arr,
                     const amrex::Array4<const amrex::Real>& cell_prim)
{
    amrex::Real thetaprime_hi = (cell_prim(i,j,k  ,PrimTheta_comp) - th0_arr(i,j,k  )) / th0_arr(i,j,k  );
    amrex::Real thetaprime_lo = (cell_prim(i,j,k-1,PrimTheta_comp) - th0_arr(i,j,k-1)) / th0_arr(i,j,k-1);

    amrex::Real thp_avg = amrex::Real(0.5) * (thetaprime_hi + thetaprime_lo);
    amrex::Real  r0_avg = amrex::Real(0.5) * (r0_arr(i,j,k) + r0_arr(i,j,k-1));

    return ( -r0_avg * grav_gpu * thp_avg);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_moist_Tpert (int& i, int& j, int& k,
                      const int& n_qstate,
                      amrex::Real const& grav_gpu,
                      amrex::Real const& rd_over_cp,
                      const amrex::Array4<const amrex::Real>& r0_arr,
                      const amrex::Array4<const amrex::Real>& th0_arr,
                      const amrex::Array4<const amrex::Real>& qv0_arr,
                      const amrex::Array4<const amrex::Real>& p0_arr,
                      const amrex::Array4<const amrex::Real>& cell_prim,
                      const amrex::Array4<const amrex::Real>& cell_data,
                      const amrex::Array4<const amrex::Real>& qt_arr)
{
    amrex::Real qv_hi = (n_qstate >= 1) ? cell_prim(i,j,k  ,PrimQ1_comp) : 0.0;
    amrex::Real qv_lo = (n_qstate >= 1) ? cell_prim(i,j,k-1,PrimQ1_comp) : 0.0;

    amrex::Real t0_hi = getTgivenPandTh(p0_arr(i,j,k  ), th0_arr(i,j,k  ), rd_over_cp);
    amrex::Real t0_lo = getTgivenPandTh(p0_arr(i,j,k-1), th0_arr(i,j,k-1), rd_over_cp);

    amrex::Real dt_hi  = getTgivenRandRTh(cell_data(i,j,k  ,Rho_comp), cell_data(i,j,k  ,RhoTheta_comp), qv_hi) - t0_hi;
    amrex::Real dt_lo  = getTgivenRandRTh(cell_data(i,j,k-1,Rho_comp), cell_data(i,j,k-1,RhoTheta_comp), qv_lo) - t0_lo;

    amrex::Real qv0_hi = qv0_arr(i,j,k  );
    amrex::Real qv0_lo = qv0_arr(i,j,k-1);

    amrex::Real q_hi = amrex::Real(0.61) * (qv_hi-qv0_hi) - (qt_arr(i,j,k  )-qv_hi) + dt_hi/t0_hi;
    amrex::Real q_lo = amrex::Real(0.61) * (qv_lo-qv0_lo) - (qt_arr(i,j,k-1)-qv_lo) + dt_lo/t0_lo;

    amrex::Real qavg = amrex::Real(0.5)  * (q_hi + q_lo);

    amrex::Real r0avg = amrex::Real(0.5) * (r0_arr(i,j,k) + r0_arr(i,j,k-1));

    return ( -r0avg * grav_gpu * qavg);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_moist_Thpert (int& i, int& j, int& k,
                       const int& n_qstate,
                       amrex::Real const& grav_gpu,
                       const amrex::Array4<const amrex::Real>& r0_arr,
                       const amrex::Array4<const amrex::Real>& th0_arr,
                       const amrex::Array4<const amrex::Real>& qv0_arr,
                       const amrex::Array4<const amrex::Real>& cell_prim,
                       const amrex::Array4<const amrex::Real>& qt_arr)
{
    amrex::Real qv_hi = (n_qstate >= 1) ? cell_prim(i,j,k  ,PrimQ1_comp) : 0.0;
    amrex::Real qv_lo = (n_qstate >= 1) ? cell_prim(i,j,k-1,PrimQ1_comp) : 0.0;

    amrex::Real qt_hi = qt_arr(i,j,k  ) - qv_hi;
    amrex::Real qt_lo = qt_arr(i,j,k-1) - qv_lo;

    amrex::Real dth_hi = cell_prim(i,j,k  ,PrimTheta_comp) - th0_arr(i,j,k  );
    amrex::Real dth_lo = cell_prim(i,j,k-1,PrimTheta_comp) - th0_arr(i,j,k-1);

    amrex::Real qv0_hi = qv0_arr(i,j,k  );
    amrex::Real qv0_lo = qv0_arr(i,j,k-1);

    amrex::Real q_hi = amrex::Real(0.61) * (qv_hi-qv0_hi) - qt_hi + dth_hi/th0_arr(i,j,k  );
    amrex::Real q_lo = amrex::Real(0.61) * (qv_lo-qv0_lo) - qt_lo + dth_lo/th0_arr(i,j,k-1);

    amrex::Real qavg  = amrex::Real(0.5) * (q_hi + q_lo);

    amrex::Real r0avg = amrex::Real(0.5) * (r0_arr(i,j,k) + r0_arr(i,j,k-1));

    return ( -r0avg * grav_gpu * qavg);
}

// **************************************************************************************
// Routines below this line are not currently used
// **************************************************************************************

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_dry_anelastic_T (int& i, int& j, int& k,
                          amrex::Real const& grav_gpu,
                          amrex::Real const& rd_over_cp,
                          const amrex::Array4<const amrex::Real>& r0_arr,
                          const amrex::Array4<const amrex::Real>& p0_arr,
                          const amrex::Array4<const amrex::Real>& cell_data)
{
    amrex::Real rt0_hi = getRhoThetagivenP(p0_arr(i,j,k));
    amrex::Real  t0_hi = getTgivenPandTh(p0_arr(i,j,k), rt0_hi/r0_arr(i,j,k), rd_over_cp);
    amrex::Real   t_hi = getTgivenPandTh(p0_arr(i,j,k), cell_data(i,j,k,RhoTheta_comp)/r0_arr(i,j,k), rd_over_cp);
    amrex::Real   q_hi  = (t_hi-t0_hi)/t0_hi;

    amrex::Real rt0_lo = getRhoThetagivenP(p0_arr(i,j,k-1));
    amrex::Real  t0_lo = getTgivenPandTh(p0_arr(i,j,k-1), rt0_lo/r0_arr(i,j,k-1), rd_over_cp);
    amrex::Real   t_lo = getTgivenPandTh(p0_arr(i,j,k-1), cell_data(i,j,k-1,RhoTheta_comp)/r0_arr(i,j,k-1), rd_over_cp);
    amrex::Real   q_lo = (t_lo-t0_lo)/t0_lo;

    amrex::Real r0_q_avg = amrex::Real(0.5) * (r0_arr(i,j,k) * q_hi + r0_arr(i,j,k-1) * q_lo);
    return (-r0_q_avg * grav_gpu);
}

#endif
