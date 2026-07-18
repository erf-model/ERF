#ifndef ERF_NOAHMP_RESULT_POLICY_H_
#define ERF_NOAHMP_RESULT_POLICY_H_

#include <cmath>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include "ERF_Constants.H"
#include "ERF_NOAHMP_Fields.H"

namespace noahmp_result_policy
{

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool result_is_valid (amrex::Real value) noexcept
{
    // Noah-MP uses approximately -9999 for unprocessed cells. ERF also
    // rejects nonfinite values before storing provider data.
    return std::isfinite(value) && value > amrex::Real(-9990.0) &&
           value < amrex::Real(0.5) * lsm_undefined;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool result_is_valid_for_output (amrex::Real value, int output_comp) noexcept
{
    if (!result_is_valid(value)) {
        return false;
    }
    if (output_comp == NoahmpOutputComp::o_noahmp_q2m_veg ||
        output_comp == NoahmpOutputComp::o_noahmp_q2m_bare) {
        return value >= amrex::Real(0.0);
    }
    if (output_comp == NoahmpOutputComp::o_noahmp_fveg) {
        return value >= amrex::Real(0.0) && value <= amrex::Real(1.0);
    }
    return true;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real select_result (bool cell_processed,
                           amrex::Real value,
                           amrex::Real undefined_value,
                           int output_comp) noexcept
{
    return cell_processed && result_is_valid_for_output(value, output_comp)
        ? value : undefined_value;
}

} // namespace noahmp_result_policy

#endif
