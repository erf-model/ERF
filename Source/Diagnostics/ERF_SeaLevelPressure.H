#ifndef ERF_SEA_LEVEL_PRESSURE_H_
#define ERF_SEA_LEVEL_PRESSURE_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>

#include "ERF_Constants.H"
#include "ERF_EOS.H"

/**
 * NMC/NGM sea-level-pressure reduction with the Shuell correction.
 *
 * The algorithm and empirical parameters are adapted from NOAA-EMC UPP,
 * sorc/ncep_post.fd/NGMSLP.f, commit
 * 1ee7b175c340a566287a987870d925e4708fe19b.
 *
 * This ERF implementation uses ERF's native thermodynamic constants and
 * reconstructs ground-interface pressure from the lowest cell-centered
 * atmospheric state. It is therefore an adaptation of the UPP algorithm,
 * not a bit-for-bit reproduction of UPP output.
 */

namespace sea_level_pressure_diagnostics
{

// Parameters inherited from NOAA-EMC UPP's NMC/NGM reduction. They are not
// ERF thermodynamic constants and intentionally remain local to this module.
constexpr amrex::Real lapse_rate = amrex::Real(0.0065); // [K m^-1]
constexpr amrex::Real shuell_critical_temperature = amrex::Real(290.66); // [K]
constexpr amrex::Real shuell_quadratic_coefficient = amrex::Real(0.005); // [K^-1]
constexpr amrex::Real sea_level_height = amrex::Real(0.0); // [m]
constexpr amrex::Real surface_height_deadband = amrex::Real(1.0); // [m]

/** Convert ERF dry-air water-vapor mixing ratio to specific humidity.
 *
 * @param[in] vapor_mixing_ratio Water vapor mass per unit dry-air mass [kg/kg].
 * @return Specific humidity used by the NMC/NGM reduction [kg/kg].
 * @pre The mixing ratio is finite and nonnegative.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real specific_humidity_from_dry_mixing_ratio (
    amrex::Real vapor_mixing_ratio) noexcept
{
    return vapor_mixing_ratio / (amrex::Real(1.0) + vapor_mixing_ratio);
}

/** Compute virtual temperature using ERF's native gas constants.
 *
 * @param[in] temperature Physical temperature [K].
 * @param[in] specific_humidity Specific humidity [kg/kg].
 * @return Virtual temperature [K].
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real virtual_temperature (amrex::Real temperature,
                                 amrex::Real specific_humidity) noexcept
{
    // UPP uses rounded 0.608; ERF uses the equivalent coefficient from its
    // shared R_d and R_v so the diagnostic remains EOS-consistent.
    return temperature * (amrex::Real(1.0) + epsv * specific_humidity);
}

/** Apply the NMC/NGM Shuell correction to the sea-level virtual temperature.
 *
 * The strict/non-strict comparisons reproduce NGMSLP.f exactly. The quadratic
 * temperature form is algebraically equivalent to UPP's tau-space form.
 *
 * @param[in] virtual_temperature_surface Ground-interface virtual temperature [K].
 * @param[in] virtual_temperature_sea_level Uncorrected sea-level virtual temperature [K].
 * @return Corrected sea-level virtual temperature [K].
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real apply_shuell_correction (amrex::Real virtual_temperature_surface,
                                     amrex::Real virtual_temperature_sea_level) noexcept
{
    if (virtual_temperature_sea_level > shuell_critical_temperature &&
        virtual_temperature_surface <= shuell_critical_temperature) {
        return shuell_critical_temperature;
    } else if (virtual_temperature_sea_level > shuell_critical_temperature &&
               virtual_temperature_surface > shuell_critical_temperature) {
        return shuell_critical_temperature -
               shuell_quadratic_coefficient *
               (virtual_temperature_surface - shuell_critical_temperature) *
               (virtual_temperature_surface - shuell_critical_temperature);
    }
    return virtual_temperature_sea_level;
}

/** Reconstruct pressure at the lower terrain interface from the lowest cell.
 *
 * @param[in] pressure_lowest_level Lowest cell-center pressure [Pa].
 * @param[in] virtual_temperature_lowest_level Lowest cell-center virtual temperature [K].
 * @param[in] virtual_temperature_surface Ground-interface virtual temperature [K].
 * @param[in] height_lowest_level Lowest cell-center physical height [m].
 * @param[in] height_surface Ground-interface physical height [m].
 * @return Ground-interface pressure [Pa].
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real reconstruct_surface_pressure (
    amrex::Real pressure_lowest_level,
    amrex::Real virtual_temperature_lowest_level,
    amrex::Real virtual_temperature_surface,
    amrex::Real height_lowest_level,
    amrex::Real height_surface) noexcept
{
    const amrex::Real mean_tau = (R_d / CONST_GRAV) *
        (virtual_temperature_lowest_level + virtual_temperature_surface) /
        amrex::Real(2.0);
    return pressure_lowest_level * std::exp(
        (height_lowest_level - height_surface) / mean_tau);
}

/** Reconstruct ground-interface pressure and reduce it to physical z=0.
 *
 * @param[in] surface_pressure Ground-interface pressure [Pa].
 * @param[in] temperature_lowest_level Lowest cell-center temperature [K].
 * @param[in] specific_humidity Lowest-level specific humidity [kg/kg].
 * @param[in] height_lowest_level Lowest cell-center physical height [m].
 * @param[in] height_surface Ground-interface physical height [m].
 * @return Pressure reduced to the physical z=0 datum [Pa].
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real reduce_surface_pressure_to_sea_level (
    amrex::Real surface_pressure,
    amrex::Real temperature_lowest_level,
    amrex::Real specific_humidity,
    amrex::Real height_lowest_level,
    amrex::Real height_surface) noexcept
{
    const amrex::Real virtual_temperature_lowest =
        virtual_temperature(temperature_lowest_level, specific_humidity);
    const amrex::Real virtual_temperature_surface =
        virtual_temperature_lowest + lapse_rate *
        (height_lowest_level - height_surface);
    const amrex::Real virtual_temperature_sea_level =
        virtual_temperature_lowest + lapse_rate *
        (height_lowest_level - sea_level_height);
    const amrex::Real corrected_virtual_temperature_sea_level =
        apply_shuell_correction(virtual_temperature_surface,
                                virtual_temperature_sea_level);
    const amrex::Real mean_tau_surface_sea_level =
        (R_d / CONST_GRAV) *
        (virtual_temperature_surface + corrected_virtual_temperature_sea_level) /
        amrex::Real(2.0);

    if (amrex::Math::abs(height_surface - sea_level_height) <=
        surface_height_deadband) {
        return surface_pressure;
    }
    return surface_pressure * std::exp(
        (height_surface - sea_level_height) / mean_tau_surface_sea_level);
}

/** Compute the complete NMC/NGM reduction from an ERF conserved state.
 *
 * ERF stores dry density, dry density times potential temperature, and (when
 * active) dry-density times vapor mixing ratio. The pressure and temperature
 * are diagnosed through ERF's EOS, then pressure is reconstructed at the
 * ground interface before the method-specific sea-level reduction.
 *
 * @param[in] dry_density Lowest-level dry-air density [kg/m^3].
 * @param[in] dry_density_times_potential_temperature Lowest-level rho_d theta [kg K/m^3].
 * @param[in] vapor_mixing_ratio Lowest-level vapor mixing ratio [kg/kg dry air].
 * @param[in] height_lowest_level Lowest cell-center physical height [m].
 * @param[in] height_surface Ground-interface physical height [m].
 * @return Pressure reduced to physical z=0 [Pa].
 * @pre State and heights are finite and physically valid; use zero mixing ratio for dry runs.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_from_erf_lowest_level_state (
    amrex::Real dry_density,
    amrex::Real dry_density_times_potential_temperature,
    amrex::Real vapor_mixing_ratio,
    amrex::Real height_lowest_level,
    amrex::Real height_surface) noexcept
{
    const amrex::Real specific_humidity =
        specific_humidity_from_dry_mixing_ratio(vapor_mixing_ratio);
    const amrex::Real pressure_lowest =
        getPgivenRTh(dry_density_times_potential_temperature,
                     vapor_mixing_ratio);
    const amrex::Real temperature_lowest =
        getTgivenRandRTh(dry_density,
                         dry_density_times_potential_temperature,
                         vapor_mixing_ratio);
    const amrex::Real virtual_temperature_lowest =
        virtual_temperature(temperature_lowest, specific_humidity);
    const amrex::Real virtual_temperature_surface =
        virtual_temperature_lowest + lapse_rate *
        (height_lowest_level - height_surface);
    const amrex::Real pressure_surface = reconstruct_surface_pressure(
        pressure_lowest, virtual_temperature_lowest, virtual_temperature_surface,
        height_lowest_level, height_surface);

    return reduce_surface_pressure_to_sea_level(
        pressure_surface, temperature_lowest, specific_humidity,
        height_lowest_level, height_surface);
}

/** Fill one selected 2D field from each column's lowest ERF state.
 *
 * @param[out] destination 2D destination MultiFab; only destination_component is written.
 * @param[in] destination_component Component to fill in destination.
 * @param[in] conserved_state Cell-centered ERF conserved state MultiFab.
 * @param[in] physical_nodal_height Absolute physical nodal heights [m].
 * @param[in] vapor_component Conserved vapor component, or -1 for a dry run.
 * @param[in] klo Lowest valid source-state vertical index for this AMR level.
 */
void fill_sea_level_pressure (amrex::MultiFab& destination,
                              int destination_component,
                              const amrex::MultiFab& conserved_state,
                              const amrex::MultiFab& physical_nodal_height,
                              int vapor_component,
                              int klo);

} // namespace sea_level_pressure_diagnostics

#endif
