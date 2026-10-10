#ifndef ERF_BUOYANCY_UTILS_H_
#define ERF_BUOYANCY_UTILS_H_

#include <ERF_EOS.H>
#include <ERF_Constants.H>

/**
 * Compute dry anelastic vertical perturbational gravitational force density [N/m^3].
 *
 * @param[in] i x-index
 * @param[in] j y-index
 * @param[in] k z-index
 * @param[in] grav_gpu Gravity
 * @param[in] r0_arr Dry reference density field
 * @param[in] th0_arr Reference potential temperature field
 * @param[in] cell_data State data FAB
 * @return Dry anelastic vertical force density [N/m^3]
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_dry_anelastic (int& i, int& j, int& k,
                        const amrex::Real& grav_gpu,
                        const amrex::Array4<const amrex::Real>& r0_arr,
                        const amrex::Array4<const amrex::Real>& th0_arr,
                        const amrex::Array4<const amrex::Real>& cell_data)
{
    // The dry anelastic branch averages theta and theta0 before forming their ratio.
    amrex::Real theta_d_hi = cell_data(i,j,k  ,RhoTheta_comp)/cell_data(i,j,k  ,Rho_comp);
    amrex::Real theta_d_lo = cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp);

    amrex::Real  theta_d_wface = myhalf * (theta_d_lo + theta_d_hi);
    amrex::Real theta_d0_wface = myhalf * (th0_arr(i,j,k) + th0_arr(i,j,k-1));
    amrex::Real     rho0_wface = myhalf * (r0_arr(i,j,k) + r0_arr(i,j,k-1));

    return (-rho0_wface * grav_gpu * (theta_d_wface - theta_d0_wface) / theta_d0_wface);
}

/**
 * Compute dry anelastic force density with embedded boundary handling.
 *
 * @param[in] i x-index
 * @param[in] j y-index
 * @param[in] k z-index
 * @param[in] grav_gpu Gravity
 * @param[in] r0_arr Reference density field
 * @param[in] th0_arr Reference potential temperature field
 * @param[in] cell_data State data FAB
 * @param[in] flag EB cell flags
 * @return Dry anelastic vertical force density [N/m^3]
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_dry_anelastic_eb (int& i, int& j, int& k,
                           const amrex::Real& grav_gpu,
                           const amrex::Array4<const amrex::Real>& r0_arr,
                           const amrex::Array4<const amrex::Real>& th0_arr,
                           const amrex::Array4<const amrex::Real>& cell_data,
                           amrex::Array4<amrex::EBCellFlag const> const& flag)
{
    amrex::Real buoyancy = zero;
    if (flag(i,j,k).isRegular()) {
        amrex::Real theta_d_hi = cell_data(i,j,k  ,RhoTheta_comp)/cell_data(i,j,k  ,Rho_comp);
        amrex::Real theta_d_lo = cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp);
        amrex::Real  theta_d_wface = myhalf * (theta_d_lo + theta_d_hi);
        amrex::Real theta_d0_wface = myhalf * (th0_arr(i,j,k) + th0_arr(i,j,k-1));
        amrex::Real     rho0_wface = myhalf * (r0_arr(i,j,k) + r0_arr(i,j,k-1));
        buoyancy = -rho0_wface * grav_gpu * (theta_d_wface - theta_d0_wface) / theta_d0_wface;
    } else if (flag(i,j,k).isSingleValued()) {
        if (flag(i,j,k-1).isCovered()) {
            amrex::Real theta_d_hihi = cell_data(i,j,k+1,RhoTheta_comp)/cell_data(i,j,k+1,Rho_comp);
            amrex::Real theta_d_hi = cell_data(i,j,k  ,RhoTheta_comp)/cell_data(i,j,k  ,Rho_comp);
            amrex::Real theta_d_wface  = myhalf * (three*theta_d_hi - theta_d_hihi);
            amrex::Real theta_d0_wface = myhalf * (three*th0_arr(i,j,k) - th0_arr(i,j,k+1));
            amrex::Real     rho0_wface = myhalf * (three*r0_arr(i,j,k) - r0_arr(i,j,k+1));
            buoyancy = -rho0_wface * grav_gpu * (theta_d_wface - theta_d0_wface) / theta_d0_wface;
        } else if (flag(i,j,k).isCovered()) {
            amrex::Real theta_d_lo = cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp);
            amrex::Real theta_d_lolo = cell_data(i,j,k-2,RhoTheta_comp)/cell_data(i,j,k-2,Rho_comp);
            amrex::Real  theta_d_wface = myhalf * (three*theta_d_lo - theta_d_lolo);
            amrex::Real theta_d0_wface = myhalf * (three*th0_arr(i,j,k-1) - th0_arr(i,j,k-2));
            amrex::Real     rho0_wface = myhalf * (three*r0_arr(i,j,k-1) - r0_arr(i,j,k-2));
            buoyancy = -rho0_wface * grav_gpu * (theta_d_wface - theta_d0_wface) / theta_d0_wface;
        } else {
            amrex::Real theta_d_hi = cell_data(i,j,k  ,RhoTheta_comp)/cell_data(i,j,k  ,Rho_comp);
            amrex::Real theta_d_lo = cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp);
            amrex::Real  theta_d_wface = myhalf * (theta_d_lo + theta_d_hi);
            amrex::Real theta_d0_wface = myhalf * (th0_arr(i,j,k) + th0_arr(i,j,k-1));
            amrex::Real     rho0_wface = myhalf * (r0_arr(i,j,k) + r0_arr(i,j,k-1));
            buoyancy = -rho0_wface * grav_gpu * (theta_d_wface - theta_d0_wface) / theta_d0_wface;
        }
    }
    return buoyancy;
}

/**
 * Compute the exact total-density-perturbation gravitational force density [N/m^3].
 *
 * @param[in] i x-index
 * @param[in] j y-index
 * @param[in] k z-index
 * @param[in] grav_gpu Gravity
 * @param[in] r0_arr Dry reference density field
 * @param[in] qv0_arr Reference water vapor mixing ratio field
 * @param[in] cell_data Conserved state FAB
 * @param[in] qt_arr Total water mixing ratio field
 * @return Vertical perturbational gravitational force density [N/m^3]
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_rhopert (int& i, int& j, int& k,
                  const amrex::Real& grav_gpu,
                  const amrex::Array4<const amrex::Real>& r0_arr,
                  const amrex::Array4<const amrex::Real>& qv0_arr,
                  const amrex::Array4<const amrex::Real>& cell_data,
                  const amrex::Array4<const amrex::Real>& qt_arr)
{
    amrex::Real rhop_hi = cell_data(i,j,k  ,Rho_comp) * (one +  qt_arr(i,j,k  )) -
                             r0_arr(i,j,k           ) * (one + qv0_arr(i,j,k  ));
    amrex::Real rhop_lo = cell_data(i,j,k-1,Rho_comp) * (one +  qt_arr(i,j,k-1)) -
                             r0_arr(i,j,k-1         ) * (one + qv0_arr(i,j,k-1));
    return( grav_gpu * myhalf * ( rhop_hi + rhop_lo ) );
}

/**
 * Compute type-1 total-density-perturbation force density with EB handling.
 *
 * @param[in] i x-index
 * @param[in] j y-index
 * @param[in] k z-index
 * @param[in] grav_gpu Gravity
 * @param[in] r0_arr Reference density field
 * @param[in] qv0_arr Reference water vapor mixing ratio field
 * @param[in] cell_data State data FAB
 * @param[in] qt_arr Total water mixing ratio field
 * @param[in] flag EB cell flags
 * @return Vertical perturbational gravitational force density [N/m^3]
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_rhopert_eb (int& i, int& j, int& k,
                     const amrex::Real& grav_gpu,
                     const amrex::Array4<const amrex::Real>& r0_arr,
                     const amrex::Array4<const amrex::Real>& qv0_arr,
                     const amrex::Array4<const amrex::Real>& cell_data,
                     const amrex::Array4<const amrex::Real>& qt_arr,
                     amrex::Array4<amrex::EBCellFlag const> const& flag)
{
    amrex::Real buoyancy = zero;
    if (flag(i,j,k).isRegular()) {
        amrex::Real rhop_hi = cell_data(i,j,k  ,Rho_comp) * (one +  qt_arr(i,j,k  )) -
                                 r0_arr(i,j,k           ) * (one + qv0_arr(i,j,k  ));
        amrex::Real rhop_lo = cell_data(i,j,k-1,Rho_comp) * (one +  qt_arr(i,j,k-1)) -
                                 r0_arr(i,j,k-1         ) * (one + qv0_arr(i,j,k-1));
        buoyancy = grav_gpu * myhalf * ( rhop_hi + rhop_lo );
    } else if (flag(i,j,k).isSingleValued()) {
        if (flag(i,j,k-1).isCovered()) {
            amrex::Real rhop_hihi = cell_data(i,j,k+1,Rho_comp) * (one +  qt_arr(i,j,k+1)) -
                                       r0_arr(i,j,k+1         ) * (one + qv0_arr(i,j,k+1));
            amrex::Real rhop_hi   = cell_data(i,j,k  ,Rho_comp) * (one +  qt_arr(i,j,k  )) -
                                       r0_arr(i,j,k           ) * (one + qv0_arr(i,j,k  ));
            buoyancy = grav_gpu * myhalf * ( three * rhop_hi - rhop_hihi );
            // buoyancy = grav_gpu * rhop_hi;
        } else if (flag(i,j,k).isCovered()) {
            amrex::Real rhop_lo   = cell_data(i,j,k-1,Rho_comp) * (one +  qt_arr(i,j,k-1)) -
                                       r0_arr(i,j,k-1         ) * (one + qv0_arr(i,j,k-1));
            amrex::Real rhop_lolo = cell_data(i,j,k-2,Rho_comp) * (one +  qt_arr(i,j,k-2)) -
                                       r0_arr(i,j,k-2         ) * (one + qv0_arr(i,j,k-2));
            buoyancy = grav_gpu * myhalf * ( three * rhop_lo - rhop_lolo );
            // buoyancy = grav_gpu * rhop_lo;
        } else {
            amrex::Real rhop_hi = cell_data(i,j,k  ,Rho_comp) * (one +  qt_arr(i,j,k  )) -
                                     r0_arr(i,j,k           ) * (one + qv0_arr(i,j,k  ));
            amrex::Real rhop_lo = cell_data(i,j,k-1,Rho_comp) * (one +  qt_arr(i,j,k-1)) -
                                     r0_arr(i,j,k-1         ) * (one + qv0_arr(i,j,k-1));
            buoyancy = grav_gpu * myhalf * ( rhop_hi + rhop_lo );
        }
    }
    return buoyancy;
}

/**
 * Compute approximate dry temperature-perturbation force density [N/m^3].
 * Compressible selectors 2 and 3 both dispatch to this implementation.
 *
 * @param[in] i x-index
 * @param[in] j y-index
 * @param[in] k z-index
 * @param[in] grav_gpu Gravity
 * @param[in] r0_arr Reference density field
 * @param[in] p0_arr Reference pressure field
 * @param[in] th0_arr Reference potential temperature field
 * @param[in] cell_data State data FAB
 * @return Approximate dry temperature force density [N/m^3]
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_dry_Tpert (int& i, int& j, int& k,
                    const amrex::Real& grav_gpu,
                    const amrex::Array4<const amrex::Real>& r0_arr,
                    const amrex::Array4<const amrex::Real>& p0_arr,
                    const amrex::Array4<const amrex::Real>& th0_arr,
                    const amrex::Array4<const amrex::Real>& cell_data)
{
    amrex::Real t0_hi = getTgivenPandTh(p0_arr(i,j,k), th0_arr(i,j,k), RdoCp);
    amrex::Real  t_hi = getTgivenRandRTh(cell_data(i,j,k  ,Rho_comp), cell_data(i,j,k  ,RhoTheta_comp));

    amrex::Real t0_lo = getTgivenPandTh(p0_arr(i,j,k-1), th0_arr(i,j,k-1), RdoCp);
    amrex::Real  t_lo = getTgivenRandRTh(cell_data(i,j,k-1,Rho_comp), cell_data(i,j,k-1,RhoTheta_comp));

    amrex::Real tprime_hi = (t_hi-t0_hi)/t0_hi;
    amrex::Real tprime_lo = (t_lo-t0_lo)/t0_lo;

    amrex::Real tp_avg = myhalf * (tprime_hi + tprime_lo);
    amrex::Real r0_avg = myhalf * (r0_arr(i,j,k) + r0_arr(i,j,k-1));

    return ( -r0_avg * grav_gpu * tp_avg);
}

/**
 * Compute approximate dry potential-temperature force density [N/m^3].
 *
 * @param[in] i x-index
 * @param[in] j y-index
 * @param[in] k z-index
 * @param[in] grav_gpu Gravity
 * @param[in] r0_arr Reference density field
 * @param[in] th0_arr Reference potential temperature field
 * @param[in] cell_prim Primitive variable FAB
 * @return Approximate dry potential-temperature force density [N/m^3]
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_dry_Thpert (int& i, int& j, int& k,
                     const amrex::Real& grav_gpu,
                     const amrex::Array4<const amrex::Real>& r0_arr,
                     const amrex::Array4<const amrex::Real>& th0_arr,
                     const amrex::Array4<const amrex::Real>& cell_prim)
{
    amrex::Real thetaprime_hi = (cell_prim(i,j,k  ,PrimTheta_comp) - th0_arr(i,j,k  )) / th0_arr(i,j,k  );
    amrex::Real thetaprime_lo = (cell_prim(i,j,k-1,PrimTheta_comp) - th0_arr(i,j,k-1)) / th0_arr(i,j,k-1);

    amrex::Real thp_avg = myhalf * (thetaprime_hi + thetaprime_lo);
    amrex::Real  r0_avg = myhalf * (r0_arr(i,j,k) + r0_arr(i,j,k-1));

    return ( -r0_avg * grav_gpu * thp_avg);
}

/**
 * Compute approximate moist temperature-perturbation force density [N/m^3].
 * Compressible selectors 2 and 3 both dispatch to this implementation.
 *
 * @param[in] i x-index
 * @param[in] j y-index
 * @param[in] k z-index
 * @param[in] n_qstate Number of moisture species
 * @param[in] grav_gpu Gravity
 * @param[in] r0_arr Reference density field
 * @param[in] th0_arr Reference potential temperature field
 * @param[in] qv0_arr Reference water vapor mixing ratio field
 * @param[in] p0_arr Reference pressure field
 * @param[in] cell_prim Primitive variable FAB
 * @param[in] cell_data State data FAB
 * @param[in] qt_arr Total water mixing ratio field
 * @return Approximate moist temperature force density [N/m^3]
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_moist_Tpert (int& i, int& j, int& k,
                      const int& n_qstate,
                      const amrex::Real& grav_gpu,
                      const amrex::Array4<const amrex::Real>& r0_arr,
                      const amrex::Array4<const amrex::Real>& th0_arr,
                      const amrex::Array4<const amrex::Real>& qv0_arr,
                      const amrex::Array4<const amrex::Real>& p0_arr,
                      const amrex::Array4<const amrex::Real>& cell_prim,
                      const amrex::Array4<const amrex::Real>& cell_data,
                      const amrex::Array4<const amrex::Real>& qt_arr)
{
    amrex::Real qv_hi = (n_qstate >= 1) ? cell_prim(i,j,k  ,PrimQ1_comp) : zero;
    amrex::Real qv_lo = (n_qstate >= 1) ? cell_prim(i,j,k-1,PrimQ1_comp) : zero;

    amrex::Real t0_hi = getTgivenPandTh(p0_arr(i,j,k  ), th0_arr(i,j,k  ), RdoCp);
    amrex::Real t0_lo = getTgivenPandTh(p0_arr(i,j,k-1), th0_arr(i,j,k-1), RdoCp);

    amrex::Real dt_hi  = getTgivenRandRTh(cell_data(i,j,k  ,Rho_comp), cell_data(i,j,k  ,RhoTheta_comp), qv_hi) - t0_hi;
    amrex::Real dt_lo  = getTgivenRandRTh(cell_data(i,j,k-1,Rho_comp), cell_data(i,j,k-1,RhoTheta_comp), qv_lo) - t0_lo;

    amrex::Real qv0_hi = qv0_arr(i,j,k  );
    amrex::Real qv0_lo = qv0_arr(i,j,k-1);

    amrex::Real q_hi = epsv * (qv_hi-qv0_hi) - (qt_arr(i,j,k  )-qv_hi) + dt_hi/t0_hi;
    amrex::Real q_lo = epsv * (qv_lo-qv0_lo) - (qt_arr(i,j,k-1)-qv_lo) + dt_lo/t0_lo;

    amrex::Real qavg = myhalf  * (q_hi + q_lo);

    amrex::Real r0avg = myhalf * (r0_arr(i,j,k) + r0_arr(i,j,k-1));

    return ( -r0avg * grav_gpu * qavg);
}

/**
 * Compute approximate moist potential-temperature force density [N/m^3].
 *
 * @param[in] i x-index
 * @param[in] j y-index
 * @param[in] k z-index
 * @param[in] n_qstate Number of moisture species
 * @param[in] grav_gpu Gravity
 * @param[in] r0_arr Reference density field
 * @param[in] th0_arr Reference potential temperature field
 * @param[in] qv0_arr Reference water vapor mixing ratio field
 * @param[in] cell_prim Primitive variable FAB
 * @param[in] qt_arr Total water mixing ratio field
 * @return Approximate moist potential-temperature force density [N/m^3]
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
buoyancy_moist_Thpert (int& i, int& j, int& k,
                       const int& n_qstate,
                       const amrex::Real& grav_gpu,
                       const amrex::Array4<const amrex::Real>& r0_arr,
                       const amrex::Array4<const amrex::Real>& th0_arr,
                       const amrex::Array4<const amrex::Real>& qv0_arr,
                       const amrex::Array4<const amrex::Real>& cell_prim,
                       const amrex::Array4<const amrex::Real>& qt_arr)
{
    amrex::Real qv_hi = (n_qstate >= 1) ? cell_prim(i,j,k  ,PrimQ1_comp) : zero;
    amrex::Real qv_lo = (n_qstate >= 1) ? cell_prim(i,j,k-1,PrimQ1_comp) : zero;

    amrex::Real qt_hi = qt_arr(i,j,k  ) - qv_hi;
    amrex::Real qt_lo = qt_arr(i,j,k-1) - qv_lo;

    amrex::Real dth_hi = cell_prim(i,j,k  ,PrimTheta_comp) - th0_arr(i,j,k  );
    amrex::Real dth_lo = cell_prim(i,j,k-1,PrimTheta_comp) - th0_arr(i,j,k-1);

    amrex::Real qv0_hi = qv0_arr(i,j,k  );
    amrex::Real qv0_lo = qv0_arr(i,j,k-1);

    amrex::Real q_hi = epsv * (qv_hi-qv0_hi) - qt_hi + dth_hi/th0_arr(i,j,k  );
    amrex::Real q_lo = epsv * (qv_lo-qv0_lo) - qt_lo + dth_lo/th0_arr(i,j,k-1);

    amrex::Real qavg  = myhalf * (q_hi + q_lo);

    amrex::Real r0avg = myhalf * (r0_arr(i,j,k) + r0_arr(i,j,k-1));

    return ( -r0avg * grav_gpu * qavg);
}


#endif
