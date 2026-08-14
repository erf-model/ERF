/**
 * \file ERF_Plotfile2DInterpolator.cpp
 */
#include "ERF_Plotfile2DInterpolator.H"
#include "ERF_Plotfile2DFill.H"

#include <cmath>
#include <sstream>

using namespace amrex;

namespace plotfile2d
{

namespace
{

struct SampledEarthWind
{
    Real u_east = Real(0.0);
    Real v_north = Real(0.0);
    Real w = Real(0.0);
};

struct SampledRotation
{
    Real cos_alpha = Real(1.0);
    Real sin_alpha = Real(0.0);
};

constexpr Real sampled_wind_dir_pi = Real(3.141592653589793238462643383279502884L);

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
SampledRotation rotation_at_mass_point (int i, int j,
                                        const Array4<const Real>& cos_alpha_arr,
                                        const Array4<const Real>& sin_alpha_arr,
                                        bool have_rotation) noexcept
{
    if (!have_rotation) {
        return {};
    }

    return {cos_alpha_arr(i, j, 0), sin_alpha_arr(i, j, 0)};
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
SampledEarthWind earth_wind_at_cell_center (int i, int j, int k,
                                            const Array4<const Real>& xvel_arr,
                                            const Array4<const Real>& yvel_arr,
                                            const Array4<const Real>& zvel_arr,
                                            const Array4<const Real>& cos_alpha_arr,
                                            const Array4<const Real>& sin_alpha_arr,
                                            bool have_rotation) noexcept
{
    const Real u_cc = Real(0.5) * (xvel_arr(i, j, k) + xvel_arr(i + 1, j, k));
    const Real v_cc = Real(0.5) * (yvel_arr(i, j, k) + yvel_arr(i, j + 1, k));
    const Real w_cc = Real(0.5) * (zvel_arr(i, j, k) + zvel_arr(i, j, k + 1));

    const SampledRotation rot = rotation_at_mass_point(i, j, cos_alpha_arr, sin_alpha_arr,
                                                       have_rotation);

    SampledEarthWind wind;
    wind.u_east = u_cc * rot.cos_alpha - v_cc * rot.sin_alpha;
    wind.v_north = u_cc * rot.sin_alpha + v_cc * rot.cos_alpha;
    wind.w = w_cc;
    return wind;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real horizontal_wind_speed (const SampledEarthWind& wind) noexcept
{
    return std::sqrt(wind.u_east * wind.u_east + wind.v_north * wind.v_north);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real meteorological_wind_direction (const SampledEarthWind& wind,
                                    Real missing_value) noexcept
{
    const Real speed = horizontal_wind_speed(wind);
    if (speed <= Real(1.0e-12)) {
        return missing_value;
    }

    Real direction = Real(270.0)
        - std::atan2(wind.v_north, wind.u_east) * Real(180.0) / sampled_wind_dir_pi;
    while (direction < Real(0.0)) {
        direction += Real(360.0);
    }
    while (direction >= Real(360.0)) {
        direction -= Real(360.0);
    }
    return direction;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real wind_field_from_earth_wind (SampledFieldID field_id,
                                 const SampledEarthWind& wind,
                                 Real missing_value) noexcept
{
    switch (field_id) {
    case SampledFieldID::UEast:     return wind.u_east;
    case SampledFieldID::VNorth:    return wind.v_north;
    case SampledFieldID::W:         return wind.w;
    case SampledFieldID::WindSpeed:  return horizontal_wind_speed(wind);
    case SampledFieldID::WindDir:    return meteorological_wind_direction(wind, missing_value);
    default:                        break;
    }
    return Real(0.0);
}

// A target is bracketed when it lies between adjacent coordinate values. This
// test works for increasing height and decreasing pressure.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
SampledBracket find_sampled_bracket (SampledCoordinate coordinate,
                                     Real target,
                                     int klo,
                                     int khi,
                                     int i, int j,
                                     const Array4<const Real>& cons_arr,
                                     const Array4<const Real>& z_phys_cc_arr,
                                     const Array4<const Real>& z_phys_nd_arr,
                                     const Array4<const Real>& p_hse_arr,
                                     bool have_z_phys_cc,
                                     bool use_hse_pressure,
                                     const MoistureComponentIndices& moisture_indices) noexcept
{
    SampledBracket bracket;
    if (khi < klo) {
        return bracket;
    }

    auto coordinate_value = [=] AMREX_GPU_DEVICE (int kk) noexcept -> Real {
        if (coordinate == SampledCoordinate::ModelIndex) {
            return static_cast<Real>(kk);
        }

        const auto field_id = (coordinate == SampledCoordinate::HeightMSL)
            ? SampledFieldID::HeightMSL
            : (coordinate == SampledCoordinate::HeightAGL)
                ? SampledFieldID::HeightAGL
                : SampledFieldID::Pressure;
        return sampled_field_value(field_id, cons_arr, z_phys_cc_arr, z_phys_nd_arr, p_hse_arr,
                                   have_z_phys_cc, use_hse_pressure, i, j, kk, moisture_indices);
    };

    const Real c0 = coordinate_value(klo);
    if (target == c0) {
        bracket.klo = klo;
        bracket.khi = klo;
        bracket.found = true;
        return bracket;
    }

    Real lo = c0;
    for (int k = klo; k < khi; ++k) {
        const Real hi = coordinate_value(k + 1);
        if (sampled_target_is_bracketed(target, lo, hi)) {
            bracket.klo = k;
            bracket.khi = k + 1;
            bracket.found = true;
            return bracket;
        }
        lo = hi;
    }

    return bracket;
}

} // namespace

void
fill_sampled_level_component (MultiFab& dst,
                              int dst_comp,
                              const Plotfile2DOutputDescriptor& descriptor,
                              const MultiFab& cons,
                              const MultiFab* z_phys_cc,
                              const MultiFab& z_phys_nd,
                              bool have_z_phys_cc,
                              const MoistureComponentIndices& moisture_indices,
                              int klo,
                              int khi,
                              const SampledWindSources& wind_sources,
                              const MultiFab* p_hse)
{
    const auto* sampled = descriptor.sampled_level ? &*descriptor.sampled_level : nullptr;
    if (sampled == nullptr) {
        fill_component_with_value(dst, dst_comp, descriptor.missing_value);
        return;
    }

    const auto* field_descriptor = find_sampled_field(sampled->source_field);
    if (field_descriptor == nullptr) {
        Abort("Unknown sampled-level source field '" + sampled->source_field + "'");
    }

    const auto coordinate = sampled_coordinate_from_string(sampled->vertical_coordinate.type);
    const auto field_id = field_descriptor->id;
    const Real target = sampled->vertical_coordinate.canonical_value;
    const int target_k = static_cast<int>(std::nearbyint(target));
    const Real missing_value = descriptor.missing_value;
    const bool is_wind_field = sampled_field_is_wind(field_id);
    const bool have_rotation = (wind_sources.cos_alpha != nullptr) || (wind_sources.sin_alpha != nullptr);
    // The caller passes the base state pressure only for an anelastic run, in which case
    //    that -- not the compressible EOS -- is the pressure of the system
    const bool use_hse_pressure = (p_hse != nullptr);

    if (is_wind_field) {
        if (wind_sources.xvel == nullptr || wind_sources.yvel == nullptr || wind_sources.zvel == nullptr) {
            Abort("Sampled-level wind field '" + sampled->source_field + "' requires xvel, yvel, and zvel sources");
        }
        if ((wind_sources.cos_alpha == nullptr) != (wind_sources.sin_alpha == nullptr)) {
            Abort("Sampled-level wind field '" + sampled->source_field + "' requires both cos_alpha and sin_alpha when rotation is supplied");
        }
    }

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const auto& dst_arr = dst.array(mfi);
        const auto& cons_arr = cons.const_array(mfi);
        const auto& z_phys_nd_arr = z_phys_nd.const_array(mfi);
        const auto& z_phys_cc_arr = have_z_phys_cc && z_phys_cc != nullptr
            ? z_phys_cc->const_array(mfi)
            : cons.const_array(mfi);
        const auto& p_hse_arr = use_hse_pressure ? p_hse->const_array(mfi) : cons.const_array(mfi);
        const auto& xvel_arr = is_wind_field ? wind_sources.xvel->const_array(mfi) : cons.const_array(mfi);
        const auto& yvel_arr = is_wind_field ? wind_sources.yvel->const_array(mfi) : cons.const_array(mfi);
        const auto& zvel_arr = is_wind_field ? wind_sources.zvel->const_array(mfi) : cons.const_array(mfi);
        const auto& cos_alpha_arr = have_rotation ? wind_sources.cos_alpha->const_array(mfi) : cons.const_array(mfi);
        const auto& sin_alpha_arr = have_rotation ? wind_sources.sin_alpha->const_array(mfi) : cons.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (coordinate == SampledCoordinate::ModelIndex) {
                if (target_k < klo || target_k > khi) {
                    dst_arr(i, j, k, dst_comp) = missing_value;
                    return;
                }

                if (is_wind_field) {
                    const SampledEarthWind wind = earth_wind_at_cell_center(
                        i, j, target_k, xvel_arr, yvel_arr, zvel_arr, cos_alpha_arr, sin_alpha_arr,
                        have_rotation);
                    dst_arr(i, j, k, dst_comp) = wind_field_from_earth_wind(field_id, wind, missing_value);
                } else {
                    dst_arr(i, j, k, dst_comp) =
                        sampled_field_value(field_id, cons_arr, z_phys_cc_arr, z_phys_nd_arr, p_hse_arr,
                                            have_z_phys_cc, use_hse_pressure, i, j, target_k, moisture_indices);
                }
                return;
            }

            const auto bracket = find_sampled_bracket(coordinate, target, klo, khi, i, j,
                                                      cons_arr, z_phys_cc_arr, z_phys_nd_arr, p_hse_arr,
                                                      have_z_phys_cc, use_hse_pressure, moisture_indices);
            if (!bracket.found) {
                dst_arr(i, j, k, dst_comp) = missing_value;
                return;
            }

            const Real coord_lo =
                sampled_field_value((coordinate == SampledCoordinate::HeightMSL)
                                        ? SampledFieldID::HeightMSL
                                        : (coordinate == SampledCoordinate::HeightAGL)
                                            ? SampledFieldID::HeightAGL
                                            : SampledFieldID::Pressure,
                                    cons_arr, z_phys_cc_arr, z_phys_nd_arr, p_hse_arr,
                                    have_z_phys_cc, use_hse_pressure, i, j, bracket.klo, moisture_indices);
            const Real coord_hi =
                sampled_field_value((coordinate == SampledCoordinate::HeightMSL)
                                        ? SampledFieldID::HeightMSL
                                        : (coordinate == SampledCoordinate::HeightAGL)
                                            ? SampledFieldID::HeightAGL
                                            : SampledFieldID::Pressure,
                                    cons_arr, z_phys_cc_arr, z_phys_nd_arr, p_hse_arr,
                                    have_z_phys_cc, use_hse_pressure, i, j, bracket.khi, moisture_indices);
            if (is_wind_field) {
                const SampledEarthWind wind_lo = earth_wind_at_cell_center(
                    i, j, bracket.klo, xvel_arr, yvel_arr, zvel_arr, cos_alpha_arr, sin_alpha_arr,
                    have_rotation);
                if (bracket.klo == bracket.khi) {
                    dst_arr(i, j, k, dst_comp) = wind_field_from_earth_wind(field_id, wind_lo, missing_value);
                    return;
                }

                const SampledEarthWind wind_hi = earth_wind_at_cell_center(
                    i, j, bracket.khi, xvel_arr, yvel_arr, zvel_arr, cos_alpha_arr, sin_alpha_arr,
                    have_rotation);
                SampledEarthWind wind_interp;
                wind_interp.u_east = linear_interpolate(wind_lo.u_east, wind_hi.u_east, coord_lo, coord_hi, target);
                wind_interp.v_north = linear_interpolate(wind_lo.v_north, wind_hi.v_north, coord_lo, coord_hi, target);
                wind_interp.w = linear_interpolate(wind_lo.w, wind_hi.w, coord_lo, coord_hi, target);
                dst_arr(i, j, k, dst_comp) = wind_field_from_earth_wind(field_id, wind_interp, missing_value);
            } else {
                const Real field_lo =
                    sampled_field_value(field_id, cons_arr, z_phys_cc_arr, z_phys_nd_arr, p_hse_arr,
                                        have_z_phys_cc, use_hse_pressure, i, j, bracket.klo, moisture_indices);
                if (bracket.klo == bracket.khi) {
                    dst_arr(i, j, k, dst_comp) = field_lo;
                    return;
                }

                const Real field_hi =
                    sampled_field_value(field_id, cons_arr, z_phys_cc_arr, z_phys_nd_arr, p_hse_arr,
                                        have_z_phys_cc, use_hse_pressure, i, j, bracket.khi, moisture_indices);
                dst_arr(i, j, k, dst_comp) =
                    linear_interpolate(field_lo, field_hi, coord_lo, coord_hi, target);
            }
        });
    }
}

} // namespace plotfile2d
