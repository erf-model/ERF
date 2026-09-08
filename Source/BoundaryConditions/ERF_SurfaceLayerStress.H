#ifndef ERF_SURFACE_LAYER_STRESS_H_
#define ERF_SURFACE_LAYER_STRESS_H_

#include "AMReX_GpuQualifiers.H"
#include "AMReX_REAL.H"

namespace surface_layer_stress {

/**
 * LSM cache value from a conservative face stress.
 *
 * @param[in] conservative_stress conservative face stress
 * @param[in] rho density used to convert to kinematic stress
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
conservative_to_kinematic_stress(
  amrex::Real conservative_stress, amrex::Real rho)
{
  return conservative_stress / rho;
}

/**
 * Test one LSM side without coupling validity to another component's handle.
 *
 * @param[in] has_flux whether the LSM flux array exists
 * @param[in] is_land whether the side is classified as land
 * @param[in] flux flux value to test against the undefined sentinel
 * @param[in] undefined sentinel value marking undefined LSM fluxes
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
lsm_flux_is_valid(
  bool has_flux, bool is_land, amrex::Real flux, amrex::Real undefined)
{
  return has_flux && is_land && flux < undefined;
}

struct FaceStressResult
{
  amrex::Real face_stress;
  amrex::Real low_kinematic_cache;
  amrex::Real high_kinematic_cache;
};

/**
 * Combine cell-centered kinematic LSM stresses with a conservative MOST
 * fallback stress at one face. The returned cache values are always
 * kinematic, while face_stress is always conservative.
 *
 * @param[in] rho_low density on the low-side cell
 * @param[in] rho_high density on the high-side cell
 * @param[in] kinematic_low low-side LSM kinematic stress
 * @param[in] kinematic_high high-side LSM kinematic stress
 * @param[in] low_valid whether the low-side LSM stress is valid
 * @param[in] high_valid whether the high-side LSM stress is valid
 * @param[in] most_face_stress conservative MOST fallback stress at the face
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
FaceStressResult
combine_lsm_and_most_stress(
  amrex::Real rho_low,
  amrex::Real rho_high,
  amrex::Real kinematic_low,
  amrex::Real kinematic_high,
  bool low_valid,
  bool high_valid,
  amrex::Real most_face_stress)
{
  const amrex::Real half = amrex::Real(0.5);
  FaceStressResult result{
    amrex::Real(0.0),
    low_valid ? kinematic_low
              : conservative_to_kinematic_stress(most_face_stress, rho_low),
    high_valid ? kinematic_high
               : conservative_to_kinematic_stress(most_face_stress, rho_high)};

  if (low_valid && high_valid) {
    result.face_stress =
      half * rho_low * kinematic_low + half * rho_high * kinematic_high;
  } else if (low_valid) {
    result.face_stress =
      half * rho_low * kinematic_low + half * most_face_stress;
  } else if (high_valid) {
    result.face_stress =
      half * rho_high * kinematic_high + half * most_face_stress;
  } else {
    result.face_stress = most_face_stress;
  }

  return result;
}

} // namespace surface_layer_stress

#endif
