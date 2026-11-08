#ifndef ERF_SLOW_RHS_PRE_UTILS_H_
#define ERF_SLOW_RHS_PRE_UTILS_H_

#include <AMReX_Array4.H>

#include <ERF_Constants.H>

/**
 * Return the pressure and buoyancy contribution to the vertical momentum RHS.
 *
 * The moisture inertia is face-centered from the two adjacent cell values,
 * matching the z-momentum staggering in erf_slow_rhs_pre. In a dry run the
 * moisture array is not read.
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
slow_rhs_pre_z_pressure_buoyancy (const amrex::Array4<const amrex::Real>& qt_arr,
                                  int i, int j, int k, bool use_moisture,
                                  amrex::Real gpz, amrex::Real abl_pressure_grad_z,
                                  amrex::Real buoyancy)
{
    amrex::Real q = use_moisture
        ? myhalf * (qt_arr(i,j,k) + qt_arr(i,j,k-1))
        : zero;
    return (-gpz - abl_pressure_grad_z + buoyancy) / (one + q);
}

#endif
