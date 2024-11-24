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
    amrex::Real theta_d_lo = cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp);
    amrex::Real theta_d_hi = cell_data(i,j,k  ,RhoTheta_comp)/cell_data(i,j,k  ,Rho_comp);

    amrex::Real theta_d_wface = amrex::Real(0.5) * (theta_d_lo + theta_d_hi);

    amrex::Real theta_d_0_lo = th0_arr(i,j,k-1);
    amrex::Real theta_d_0_hi = th0_arr(i,j,k  );

    amrex::Real theta_d_0_wface = amrex::Real(0.5) * (theta_d_0_lo + theta_d_0_hi);

    amrex::Real rho0_wface = amrex::Real(0.5) * (r0_arr(i,j,k) + r0_arr(i,j,k-1));

    return (-rho0_wface * grav_gpu * (theta_d_wface - theta_d_0_wface) / theta_d_0_wface);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_moist_anelastic (int& i, int& j, int& k,
                          amrex::Real const& grav_gpu,
                          amrex::Real const& rv_over_rd,
                          const amrex::Array4<const amrex::Real>&  r0_arr,
                          const amrex::Array4<const amrex::Real>& th0_arr,
                          const amrex::Array4<const amrex::Real>& cell_data)
{
    amrex::Real theta_d_lo = cell_data(i,j,k-1,RhoTheta_comp)/r0_arr(i,j,k-1);
    amrex::Real qv_lo      = cell_data(i,j,k-1,RhoQ1_comp)   /r0_arr(i,j,k-1);
    amrex::Real qc_lo      = cell_data(i,j,k-1,RhoQ2_comp)   /r0_arr(i,j,k-1);
    amrex::Real qt_lo      = qv_lo + qc_lo;
    amrex::Real theta_v_lo = theta_d_lo * (1.0 - (1.0 - rv_over_rd)*qt_lo - rv_over_rd*qc_lo);

    amrex::Real theta_d_hi = cell_data(i,j,k,RhoTheta_comp)/r0_arr(i,j,k);
    amrex::Real qv_hi      = cell_data(i,j,k,RhoQ1_comp)   /r0_arr(i,j,k);
    amrex::Real qc_hi      = cell_data(i,j,k,RhoQ2_comp)   /r0_arr(i,j,k);
    amrex::Real qt_hi      = qv_hi + qc_hi;
    amrex::Real theta_v_hi = theta_d_hi * (1.0 - (1.0 - rv_over_rd)*qt_hi - rv_over_rd*qc_hi);

    amrex::Real theta_v_wface = amrex::Real(0.5) * (theta_v_lo + theta_v_hi);

    amrex::Real theta_v_0_lo = th0_arr(i,j,k-1) * (1.0 - (1.0 - rv_over_rd)*qt_lo - rv_over_rd*qc_lo);
    amrex::Real theta_v_0_hi = th0_arr(i,j,k  ) * (1.0 - (1.0 - rv_over_rd)*qt_hi - rv_over_rd*qc_hi);

    amrex::Real theta_v_0_wface = amrex::Real(0.5) * (theta_v_0_lo + theta_v_0_hi);

    amrex::Real rho0_wface = amrex::Real(0.5) * (r0_arr(i,j,k) + r0_arr(i,j,k-1));

    return (-rho0_wface * grav_gpu * (theta_v_wface - theta_v_0_wface) / theta_v_0_wface);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_dry_default (int& i, int& j, int& k,
                      amrex::Real const& grav_gpu,
                      amrex::Real const& rd_over_cp,
                      const amrex::Array4<const amrex::Real>& r0_arr,
                      const amrex::Array4<const amrex::Real>& p0_arr,
                      const amrex::Array4<const amrex::Real>& th0_arr,
                      const amrex::Array4<const amrex::Real>& cell_data)
{
    amrex::Real  t0_hi = getTgivenPandTh(p0_arr(i,j,k), th0_arr(i,j,k), rd_over_cp);
    amrex::Real   t_hi = getTgivenRandRTh(cell_data(i,j,k  ,Rho_comp), cell_data(i,j,k  ,RhoTheta_comp));
    amrex::Real qplus  = (t_hi-t0_hi)/t0_hi;

    amrex::Real  t0_lo = getTgivenPandTh(p0_arr(i,j,k-1), th0_arr(i,j,k-1), rd_over_cp);
    amrex::Real   t_lo = getTgivenRandRTh(cell_data(i,j,k-1,Rho_comp), cell_data(i,j,k-1,RhoTheta_comp));
    amrex::Real qminus = (t_lo-t0_lo)/t0_lo;

    amrex::Real r0_q_avg = amrex::Real(0.5) * (r0_arr(i,j,k) * qplus + r0_arr(i,j,k-1) * qminus);
    return (-r0_q_avg * grav_gpu);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_rhopert (int& i, int& j, int& k,
                  const int& n_qstate,
                  amrex::Real const& grav_gpu,
                  const amrex::Array4<const amrex::Real>& r0_arr,
                  const amrex::Array4<const amrex::Real>& cell_data)
{
    amrex::Real rhop_hi = cell_data(i,j,k  ,Rho_comp) - r0_arr(i,j,k  );
    amrex::Real rhop_lo = cell_data(i,j,k-1,Rho_comp) - r0_arr(i,j,k-1);
    for (int q_offset(0); q_offset<n_qstate; ++q_offset) {
        rhop_hi += cell_data(i,j,k  ,RhoQ1_comp+q_offset);
        rhop_lo += cell_data(i,j,k-1,RhoQ1_comp+q_offset);
    }
    return( grav_gpu * amrex::Real(0.5) * ( rhop_hi + rhop_lo ) );
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_type2 (int& i, int& j, int& k,
                const int& n_qstate,
                amrex::Real const& grav_gpu,
                amrex::Real* rho_d_ptr,
                amrex::Real* theta_d_ptr,
                amrex::Real* qv_d_ptr,
                amrex::Real* qc_d_ptr,
                amrex::Real* qp_d_ptr,
                const amrex::Array4<const amrex::Real>& cell_prim,
                const amrex::Array4<const amrex::Real>& cell_data)
{
    amrex::Real tempp1d = getTgivenRandRTh(rho_d_ptr[k  ], rho_d_ptr[k  ]*theta_d_ptr[k  ], qv_d_ptr[k  ]);
    amrex::Real tempm1d = getTgivenRandRTh(rho_d_ptr[k-1], rho_d_ptr[k-1]*theta_d_ptr[k-1], qv_d_ptr[k-1]);

    amrex::Real tempp3d  = getTgivenRandRTh(cell_data(i,j,k  ,Rho_comp),
                                            cell_data(i,j,k  ,RhoTheta_comp),
                                            cell_data(i,j,k  ,RhoQ1_comp)/cell_data(i,j,k  ,Rho_comp));
    amrex::Real tempm3d  = getTgivenRandRTh(cell_data(i,j,k-1,Rho_comp),
                                            cell_data(i,j,k-1,RhoTheta_comp),
                                            cell_data(i,j,k-1,RhoQ1_comp)/cell_data(i,j,k-1,Rho_comp));

    amrex::Real qv_plus  = (n_qstate >= 1) ? cell_prim(i,j,k  ,PrimQ1_comp) : 0.0;
    amrex::Real qv_minus = (n_qstate >= 1) ? cell_prim(i,j,k-1,PrimQ1_comp) : 0.0;

    amrex::Real qc_plus  = (n_qstate >= 2) ? cell_prim(i,j,k  ,PrimQ2_comp) : 0.0;
    amrex::Real qc_minus = (n_qstate >= 2) ? cell_prim(i,j,k-1,PrimQ2_comp) : 0.0;

    amrex::Real qp_plus  = (n_qstate >= 3) ? cell_prim(i,j,k  ,PrimQ3_comp) : 0.0;
    amrex::Real qp_minus = (n_qstate >= 3) ? cell_prim(i,j,k-1,PrimQ3_comp) : 0.0;

    amrex::Real qplus  = 0.61 * ( qv_plus - qv_d_ptr[k] ) -
                                ( qc_plus - qc_d_ptr[k]   +
                                  qp_plus - qp_d_ptr[k] )
                       + (tempp3d-tempp1d)/tempp1d*(amrex::Real(1.0) + amrex::Real(0.61)*qv_d_ptr[k]-qc_d_ptr[k]-qp_d_ptr[k]);

    amrex::Real qminus = 0.61 * ( qv_minus - qv_d_ptr[k-1] ) -
                                ( qc_minus - qc_d_ptr[k-1]   +
                                  qp_minus - qp_d_ptr[k-1] )
                       + (tempm3d-tempm1d)/tempm1d*(amrex::Real(1.0) + amrex::Real(0.61)*qv_d_ptr[k-1]-qc_d_ptr[k-1]-qp_d_ptr[k-1]);

    amrex::Real qavg  = amrex::Real(0.5) * (qplus + qminus);
    amrex::Real r0avg = amrex::Real(0.5) * (rho_d_ptr[k] + rho_d_ptr[k-1]);

    return (-qavg * r0avg * grav_gpu);
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_type3 (int& i, int& j, int& k,
                const int& n_qstate,
                amrex::Real const& grav_gpu,
                amrex::Real* rho_d_ptr,
                amrex::Real* theta_d_ptr,
                amrex::Real* qv_d_ptr,
                const amrex::Array4<const amrex::Real>& cell_prim,
                const amrex::Array4<const amrex::Real>& cell_data)
{
    amrex::Real tempp1d = getTgivenRandRTh(rho_d_ptr[k  ], rho_d_ptr[k  ]*theta_d_ptr[k  ], qv_d_ptr[k  ]);
    amrex::Real tempm1d = getTgivenRandRTh(rho_d_ptr[k-1], rho_d_ptr[k-1]*theta_d_ptr[k-1], qv_d_ptr[k-1]);

    amrex::Real tempp3d  = getTgivenRandRTh(cell_data(i,j,k  ,Rho_comp),
                                            cell_data(i,j,k  ,RhoTheta_comp),
                                            cell_data(i,j,k  ,RhoQ1_comp)/cell_data(i,j,k  ,Rho_comp));
    amrex::Real tempm3d  = getTgivenRandRTh(cell_data(i,j,k-1,Rho_comp),
                                            cell_data(i,j,k-1,RhoTheta_comp),
                                            cell_data(i,j,k-1,RhoQ1_comp)/cell_data(i,j,k-1,Rho_comp));

    amrex::Real qv_plus  = (n_qstate >= 1) ? cell_prim(i,j,k  ,PrimQ1_comp) : 0.0;
    amrex::Real qv_minus = (n_qstate >= 1) ? cell_prim(i,j,k-1,PrimQ1_comp) : 0.0;

    amrex::Real qc_plus  = (n_qstate >= 2) ? cell_prim(i,j,k  ,PrimQ2_comp) : 0.0;
    amrex::Real qc_minus = (n_qstate >= 2) ? cell_prim(i,j,k-1,PrimQ2_comp) : 0.0;

    amrex::Real qp_plus  = (n_qstate >= 3) ? cell_prim(i,j,k  ,PrimQ3_comp) : 0.0;
    amrex::Real qp_minus = (n_qstate >= 3) ? cell_prim(i,j,k-1,PrimQ3_comp) : 0.0;

    amrex::Real qplus  = 0.61 * qv_plus  - (qc_plus  + qp_plus)  + (tempp3d-tempp1d)/tempp1d;

    amrex::Real qminus = 0.61 * qv_minus - (qc_minus + qp_minus) + (tempm3d-tempm1d)/tempm1d;

    amrex::Real qavg  = amrex::Real(0.5) * (qplus + qminus);
    amrex::Real r0avg = amrex::Real(0.5) * (rho_d_ptr[k] + rho_d_ptr[k-1]);

    return ( -qavg * r0avg * grav_gpu );
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_type4 (int& i, int& j, int& k,
                const int& n_qstate,
                amrex::Real const& grav_gpu,
                amrex::Real* rho_d_ptr,
                amrex::Real* theta_d_ptr,
                amrex::Real* qv_d_ptr,
                amrex::Real* qc_d_ptr,
                amrex::Real* qp_d_ptr,
                const amrex::Array4<const amrex::Real>& cell_prim,
                const amrex::Array4<const amrex::Real>& cell_data)
{
    amrex::Real qv_plus  = (n_qstate >= 1) ? cell_prim(i,j,k  ,PrimQ1_comp) : 0.0;
    amrex::Real qv_minus = (n_qstate >= 1) ? cell_prim(i,j,k-1,PrimQ1_comp) : 0.0;

    amrex::Real qc_plus  = (n_qstate >= 2) ? cell_prim(i,j,k  ,PrimQ2_comp) : 0.0;
    amrex::Real qc_minus = (n_qstate >= 2) ? cell_prim(i,j,k-1,PrimQ2_comp) : 0.0;

    amrex::Real qp_plus  = (n_qstate >= 3) ? cell_prim(i,j,k  ,PrimQ3_comp) : 0.0;
    amrex::Real qp_minus = (n_qstate >= 3) ? cell_prim(i,j,k-1,PrimQ3_comp) : 0.0;

    amrex::Real qplus  = amrex::Real(0.61) * ( qv_plus - qv_d_ptr[k] ) -
                                             ( qc_plus - qc_d_ptr[k]   +
                                               qp_plus - qp_d_ptr[k] )
                       + (cell_data(i,j,k  ,RhoTheta_comp)/cell_data(i,j,k  ,Rho_comp) - theta_d_ptr[k  ])/theta_d_ptr[k  ];

    amrex::Real qminus = amrex::Real(0.61) * ( qv_minus - qv_d_ptr[k-1] ) -
                                             ( qc_minus - qc_d_ptr[k-1]   +
                                               qp_minus - qp_d_ptr[k-1] )
                       + (cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp) - theta_d_ptr[k-1])/theta_d_ptr[k-1];

    amrex::Real qavg  = amrex::Real(0.5) * (qplus + qminus);
    amrex::Real r0avg = amrex::Real(0.5) * (rho_d_ptr[k] + rho_d_ptr[k-1]);

    return (-qavg * r0avg * grav_gpu);
}

// **************************************************************************************
// Routines below this line are not currently used
// **************************************************************************************

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
    amrex::Real  t0_hi = getTgivenPandTh(p0_arr(i,j,k), th0_arr(i,j,k), rd_over_cp);
    amrex::Real   t_hi = getTgivenRandRTh(cell_data(i,j,k  ,Rho_comp), cell_data(i,j,k  ,RhoTheta_comp));
    amrex::Real qplus  = (t_hi-t0_hi)/t0_hi;

    amrex::Real  t0_lo = getTgivenPandTh(p0_arr(i,j,k-1), th0_arr(i,j,k-1), rd_over_cp);
    amrex::Real   t_lo = getTgivenRandRTh(cell_data(i,j,k-1,Rho_comp), cell_data(i,j,k-1,RhoTheta_comp));
    amrex::Real qminus = (t_lo-t0_lo)/t0_lo;

    amrex::Real r0_q_avg = amrex::Real(0.5) * (r0_arr(i,j,k) * qplus + r0_arr(i,j,k-1) * qminus);
    return (-r0_q_avg * grav_gpu);
}

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
    amrex::Real qplus  = (t_hi-t0_hi)/t0_hi;

    amrex::Real rt0_lo = getRhoThetagivenP(p0_arr(i,j,k-1));
    amrex::Real  t0_lo = getTgivenPandTh(p0_arr(i,j,k-1), rt0_lo/r0_arr(i,j,k-1), rd_over_cp);
    amrex::Real   t_lo = getTgivenPandTh(p0_arr(i,j,k-1), cell_data(i,j,k-1,RhoTheta_comp)/r0_arr(i,j,k-1), rd_over_cp);
    amrex::Real qminus = (t_lo-t0_lo)/t0_lo;

    amrex::Real r0_q_avg = amrex::Real(0.5) * (r0_arr(i,j,k) * qplus + r0_arr(i,j,k-1) * qminus);
    return (-r0_q_avg * grav_gpu);
}

#endif
