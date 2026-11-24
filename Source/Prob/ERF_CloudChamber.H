#ifndef ERF_CLOUD_CHAMBER_H_
#define ERF_CLOUD_CHAMBER_H_

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include "ERF_Constants.H"
#include "ERF_CloudChamberTemperatureBands.H"
#include "ERF_DataStruct.H"
#include "ERF_EOS.H"
#include "ERF_MicrophysicsUtils.H"
#include "BoundaryConditions/ERF_WallThermodynamics.H"
#include "Diffusion/ERF_WallGradientCoeffs.H"

#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace erf_cloud_chamber {

/**
 * Cloud Chamber wall metadata.  This is deliberately a small, trivially
 * copyable description: it can be captured by an AMReX device kernel without
 * virtual dispatch.  The per-channel contract supports resolved molecular,
 * fixed bulk-aerodynamic, and neutral roughness-log transfer, with MOST
 * coefficient generation on the supported horizontal faces.
 */
enum class WallFace : unsigned char { XLo, XHi, YLo, YHi, ZLo, ZHi };

/**
 * @brief Profile data used by the device initializer; no wall transport is exposed here.
 */
struct Initialization {
    bool active = false;
    bool has_moisture = false;
    bool physical_initialization = false;
    amrex::Real initial_temperature_bottom = amrex::Real(0.0);
    amrex::Real initial_temperature_top = amrex::Real(0.0);
    amrex::Real initial_relative_humidity = amrex::Real(0.0);
    amrex::Real temperature_perturbation_amplitude = amrex::Real(0.0);
    amrex::Real theta_bottom = amrex::Real(0.0); // legacy_theta_qv only
    amrex::Real theta_top = amrex::Real(0.0);    // legacy_theta_qv only
    amrex::Real theta_perturbation_amplitude = amrex::Real(0.0); // legacy
    amrex::Real qv_bottom = amrex::Real(0.0);    // legacy_theta_qv only
    amrex::Real qv_top = amrex::Real(0.0);       // legacy_theta_qv only
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo = {};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_hi = {};
};

struct Config : Initialization {
    // Only the full solver-side parser can publish walls. The default (inactive)
    // config has no boundary, rather than a plausible-looking array of zero speeds.
    [[nodiscard]] const erf_wall_thermodynamics::Boundary& transport_boundary () const &
    {
        return m_boundary.value();
    }
    const erf_wall_thermodynamics::Boundary& transport_boundary () const && = delete;

    [[nodiscard]] erf_wall_thermodynamics::SidewallBandViews temperature_band_views () const & noexcept
    {
        return {m_bands[0].view(), m_bands[1].view(), m_bands[2].view(), m_bands[3].view()};
    }
    erf_wall_thermodynamics::SidewallBandViews temperature_band_views () const && = delete;

    static bool has_bulk_scalar_wall (const erf_wall_thermodynamics::Boundary& walls) noexcept
    {
        using erf_wall_thermodynamics::ScalarModel;
        for (const auto& wall : walls) {
            if ((wall.heat_is_conductive() && wall.heat.model == ScalarModel::BulkAero) ||
                (wall.vapor_is_wet() && wall.vapor.model == ScalarModel::BulkAero)) {
                return true;
            }
        }
        return false;
    }

    bool has_wall_rate_channel () const
    {
        for (const auto& wall : transport_boundary()) {
            if (wall.has_wall_rate_channel()) { return true; }
        }
        return false;
    }

private:
    friend Config parse_config (const amrex::Geometry&, MoistureType, const DiffChoice&);
    std::optional<erf_wall_thermodynamics::Boundary> m_boundary;
    std::array<SidewallTemperatureBands, 4> m_bands;
};

enum class InitializationMode : unsigned char {
    PhysicalTemperatureRH,
    LegacyThetaQv
};

/**
 * @brief Contract used to validate initialization parameter consistency.
 */
struct InitializationContract {
    InitializationMode mode = InitializationMode::PhysicalTemperatureRH;
    // True for any active moisture model, including noncondensing transport.
    bool has_moisture = false;
    bool has_initial_relative_humidity = false;
    bool has_legacy_profile_keys = false;
    bool has_physical_profile_keys = false;
};

/**
 * @brief Validate the initialization contract and return an error message if violated.
 * @param[in] contract The contract to validate.
 * @return Error message string, empty if valid.
 */
inline std::string
initialization_contract_error (const InitializationContract& contract)
{
    switch (contract.mode) {
    case InitializationMode::PhysicalTemperatureRH:
        if (contract.has_moisture) {
            if (!contract.has_initial_relative_humidity) {
                return "Cloud Chamber: physical moist initialization requires prob.initial_relative_humidity";
            }
        } else if (contract.has_initial_relative_humidity) {
            return "Cloud Chamber: prob.initial_relative_humidity requires an active moisture model";
        }
        if (contract.has_legacy_profile_keys) {
            return "Cloud Chamber: physical_temperature_rh cannot be combined with legacy theta/qv profile keys";
        }
        break;
    case InitializationMode::LegacyThetaQv:
        if (contract.has_physical_profile_keys) {
            return "Cloud Chamber: legacy_theta_qv cannot be combined with physical temperature/RH profile keys";
        }
        break;
    }

    return {};
}

/**
 * @brief Return whether a moisture type is supported by the Cloud Chamber.
 *
 * This is intentionally an allowlist.  In particular, spectral SBM remains
 * outside the chamber's supported moisture configurations even though it is
 * a valid ERF moisture type.
 */
inline bool
supports_moisture_model (const MoistureType moisture_type) noexcept
{
    switch (moisture_type) {
    case MoistureType::None:
    case MoistureType::SatAdj:
    case MoistureType::MoistNoCondensation:
    case MoistureType::Kessler_NoRain:
    case MoistureType::Kessler:
    case MoistureType::SAM_NoPrecip_NoIce:
    case MoistureType::SAM_NoIce:
    case MoistureType::SAM:
    case MoistureType::WSM6:
    case MoistureType::WDM6:
    case MoistureType::Morrison_NoIce:
    case MoistureType::Morrison:
    case MoistureType::SuperDroplets:
        return true;
    case MoistureType::SBM:
        return false;
    }
    return false;
}

/**
 * @brief Compute a linear profile between two values.
 * @param[in] bottom Value at the lower boundary.
 * @param[in] top Value at the upper boundary.
 * @param[in] coordinate Current coordinate.
 * @param[in] lower Lower boundary coordinate.
 * @param[in] length Length of the domain.
 * @return Interpolated value.
 */
AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
amrex::Real
linear_profile (amrex::Real bottom, amrex::Real top,
                amrex::Real coordinate, amrex::Real lower,
                amrex::Real length) noexcept
{
    return bottom + (top - bottom) * (coordinate - lower) / length;
}

/**
 * @brief Compute a deterministic 3D sine perturbation.
 * @param[in] x x-coordinate.
 * @param[in] y y-coordinate.
 * @param[in] z z-coordinate.
 * @param[in] lo Lower boundary coordinates.
 * @param[in] length Domain lengths.
 * @param[in] amplitude Perturbation amplitude.
 * @return Perturbation value at the given coordinates.
 */
AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
amrex::Real
deterministic_perturbation (amrex::Real x, amrex::Real y, amrex::Real z,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& lo,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& length,
                            amrex::Real amplitude) noexcept
{
    constexpr amrex::Real pi = amrex::Real(3.1415926535897932384626433832795);
    const amrex::Real sx = std::sin(amrex::Real(2.0) * pi * (x - lo[0]) / length[0]);
    const amrex::Real sy = std::sin(amrex::Real(2.0) * pi * (y - lo[1]) / length[1]);
    const amrex::Real sz = std::sin(pi * (z - lo[2]) / length[2]);
    return amplitude * sx * sy * sz;
}

/**
 * @brief Discrete derivative ratio for the staggered initial velocity field.
 *
 * The ratio makes the finite-volume x/y divergence cancel on rectangular
 * grids.  A one-cell horizontal direction has an identically zero sine
 * pattern, so the perturbation is zero rather than dividing by a zero
 * discrete derivative factor.
 */
AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
amrex::Real
staggered_velocity_v_scale (amrex::Real lx, amrex::Real ly,
                            amrex::Real dx, amrex::Real dy,
                            int nx, int ny) noexcept
{
    if (nx <= 1 || ny <= 1) { return amrex::Real(0.0); }
    if (lx == ly && dx == dy) { return amrex::Real(1.0); }
    constexpr amrex::Real pi = amrex::Real(3.1415926535897932384626433832795);
    const amrex::Real kx_d = std::sin(pi * dx / lx) / dx;
    const amrex::Real ky_d = std::sin(pi * dy / ly) / dy;
    return ky_d == amrex::Real(0.0) ? amrex::Real(0.0) : kx_d / ky_d;
}

/**
 * @brief Evaluate the x-face component of the Cloud Chamber initial velocity.
 */
AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
amrex::Real
initial_u_velocity_perturbation (
    amrex::Real x, amrex::Real y, amrex::Real z,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& lo,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& length,
    amrex::Real U_0, int nx, int ny) noexcept
{
    if (nx <= 1 || ny <= 1) { return amrex::Real(0.0); }
    constexpr amrex::Real pi = amrex::Real(3.1415926535897932384626433832795);
    const amrex::Real sx = std::sin(pi * (x - lo[0]) / length[0]);
    const amrex::Real sz = std::sin(pi * (z - lo[2]) / length[2]);
    return U_0 * sx * sx *
        std::sin(amrex::Real(2.0) * pi * (y - lo[1]) / length[1]) * sz * sz;
}

/**
 * @brief Evaluate the y-face component of the Cloud Chamber initial velocity.
 */
AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
amrex::Real
initial_v_velocity_perturbation (
    amrex::Real x, amrex::Real y, amrex::Real z,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& lo,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& length,
    amrex::Real U_0, amrex::Real v_scale, int nx, int ny) noexcept
{
    if (nx <= 1 || ny <= 1) { return amrex::Real(0.0); }
    constexpr amrex::Real pi = amrex::Real(3.1415926535897932384626433832795);
    const amrex::Real sy = std::sin(pi * (y - lo[1]) / length[1]);
    const amrex::Real sz = std::sin(pi * (z - lo[2]) / length[2]);
    return -U_0 * v_scale *
        std::sin(amrex::Real(2.0) * pi * (x - lo[0]) / length[0]) *
        sy * sy * sz * sz;
}

/**
 * @brief The Cloud Chamber initial vertical velocity perturbation.
 */
AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
amrex::Real
initial_w_velocity_perturbation () noexcept
{
    return amrex::Real(0.0);
}

/**
 * @brief Compute the potential temperature at a given location.
 * @param[in] config Configuration parameters.
 * @param[in] x x-coordinate.
 * @param[in] y y-coordinate.
 * @param[in] z z-coordinate.
 * @return Potential temperature at the specified location.
 */
AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
amrex::Real
theta_at (const Initialization& config, amrex::Real x, amrex::Real y,
          amrex::Real z) noexcept
{
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> length = {
        config.prob_hi[0] - config.prob_lo[0],
        config.prob_hi[1] - config.prob_lo[1],
        config.prob_hi[2] - config.prob_lo[2]};
    const amrex::Real bottom = config.physical_initialization ?
        config.initial_temperature_bottom : config.theta_bottom;
    const amrex::Real top = config.physical_initialization ?
        config.initial_temperature_top : config.theta_top;
    const amrex::Real amplitude = config.physical_initialization ?
        config.temperature_perturbation_amplitude : config.theta_perturbation_amplitude;
    return linear_profile(bottom, top, z,
                          config.prob_lo[2], length[2]) +
           deterministic_perturbation(x, y, z, config.prob_lo, length,
                                      amplitude);
}

/**
 * @brief Compute the water vapor mixing ratio at a given height.
 * @param[in] config Configuration parameters.
 * @param[in] z Vertical coordinate.
 * @return Vapor mixing ratio at the specified height.
 */
AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
amrex::Real
qv_at (const Initialization& config, amrex::Real z) noexcept
{
    return linear_profile(config.qv_bottom, config.qv_top, z,
                          config.prob_lo[2],
                          config.prob_hi[2] - config.prob_lo[2]);
}

/**
 * @brief Convert relative humidity to water vapor mixing ratio.
 * @param[in] temperature_K Temperature in Kelvin.
 * @param[in] pressure_Pa Pressure in Pascals.
 * @param[in] relative_humidity Relative humidity [0,1].
 * @return Water vapor mixing ratio.
 */
AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
amrex::Real
vapor_mixing_ratio_from_relative_humidity (amrex::Real temperature_K,
                                           amrex::Real pressure_Pa,
                                           amrex::Real relative_humidity) noexcept
{
    // erf_esatw returns hPa; convert saturation and vapor pressure to Pa before
    // applying the dry-air vapor mixing-ratio definition.
    const amrex::Real es_Pa = erf_esatw(temperature_K) * amrex::Real(100.0);
    const amrex::Real ev_Pa = relative_humidity * es_Pa;
    return RdoRv * ev_Pa / (pressure_Pa - ev_Pa);
}

inline bool finite (amrex::Real value) noexcept
{
    return std::isfinite(static_cast<double>(value));
}

inline void require_finite (const std::string& key, amrex::Real value)
{
    if (!finite(value)) {
        amrex::Error("Cloud Chamber: " + key + " must be finite");
    }
}

inline std::string face_name (WallFace face)
{
    switch (face) {
    case WallFace::XLo: return "xlo";
    case WallFace::XHi: return "xhi";
    case WallFace::YLo: return "ylo";
    case WallFace::YHi: return "yhi";
    case WallFace::ZLo: return "zlo";
    case WallFace::ZHi: return "zhi";
    }
    return "unknown";
}

inline bool valid_wall_wetness (amrex::Real rh) noexcept
{
    return finite(rh) && rh >= amrex::Real(0.0) && rh <= amrex::Real(1.0);
}

inline std::string
parse_wall_wetness (amrex::ParmParse& pp, const std::string& face,
                    erf_wall_thermodynamics::FaceWall& wall)
{
    if (!pp.contains("wetness")) { return {}; }
    if (wall.moisture != erf_wall_thermodynamics::MoistureMode::WetEquilibrium) {
        return "Cloud Chamber: " + face + ".wetness is only valid on physical walls with moisture = wet";
    }
    amrex::Real rh;
    pp.get("wetness", rh);
    if (!valid_wall_wetness(rh)) {
        return "Cloud Chamber: " + face + ".wetness must be a finite relative humidity fraction in [0,1]";
    }
    wall.wetness = rh;
    return {};
}

struct WallTransferContract {
    bool legacy_aggregate_specified = false;
    bool momentum_model_specified = false;
    bool heat_model_specified = false;
    bool vapor_model_specified = false;
    bool coefficient_source_specified = false;
    bool momentum_coefficient_specified = false;
    bool heat_coefficient_specified = false;
    bool vapor_coefficient_specified = false;
    bool z0_m_specified = false;
    bool z0_h_specified = false;
    bool z0_q_specified = false;
    std::string momentum_model;
    std::string heat_model;
    std::string vapor_model;
    std::string coefficient_source;
    amrex::Real momentum_coefficient = amrex::Real(0.0);
    amrex::Real heat_coefficient = amrex::Real(0.0);
    amrex::Real vapor_coefficient = amrex::Real(0.0);
    amrex::Real z0_m = amrex::Real(0.0);
    amrex::Real z0_h = amrex::Real(0.0);
    amrex::Real z0_q = amrex::Real(0.0);
};

inline WallTransferContract read_wall_transfer_contract (const amrex::ParmParse& pp,
                                                         const std::string& face)
{
    WallTransferContract contract;
    auto selector = [&](const char* key, std::string& value) {
        if (!pp.query(key, value)) { return false; }
        value = amrex::toLower(value);
        return true;
    };
    std::string aggregate;
    contract.legacy_aggregate_specified = selector("wall_transfer_model", aggregate);
    if (contract.legacy_aggregate_specified && aggregate != "resolved_molecular") {
        amrex::Error("Cloud Chamber: " + face + ".wall_transfer_model only accepts "
                     "resolved_molecular; select other models with the per-channel keys");
    }
    contract.momentum_model_specified = selector("momentum_transfer_model", contract.momentum_model);
    contract.heat_model_specified = selector("heat_transfer_model", contract.heat_model);
    contract.vapor_model_specified = selector("vapor_transfer_model", contract.vapor_model);
    contract.coefficient_source_specified = selector("coefficient_source", contract.coefficient_source);
    auto coefficient = [&](const char* key, amrex::Real& value) {
        if (!pp.query(key, value)) { return false; }
        require_finite(face + "." + key, value);
        return true;
    };
    contract.momentum_coefficient_specified = coefficient("C_D", contract.momentum_coefficient);
    contract.heat_coefficient_specified = coefficient("C_H", contract.heat_coefficient);
    contract.vapor_coefficient_specified = coefficient("C_E", contract.vapor_coefficient);
    contract.z0_m_specified = coefficient("z0_m", contract.z0_m);
    contract.z0_h_specified = coefficient("z0_h", contract.z0_h);
    contract.z0_q_specified = coefficient("z0_q", contract.z0_q);
    return contract;
}

inline bool is_bulk_model (const std::string& model)
{
    return model == "bulk_aero" || model == "bulk_aero_plus_molecular";
}

inline std::string
wall_transfer_contract_error (const WallTransferContract& contract,
                              const std::string& face)
{
    const bool momentum_neutral = contract.momentum_model_specified &&
        (contract.momentum_model == "neutral_roughness_log" ||
         contract.momentum_model == "law_of_wall_momentum");
    const bool momentum_bulk = contract.momentum_model_specified &&
        is_bulk_model(contract.momentum_model);
    const bool heat_neutral = contract.heat_model_specified &&
        contract.heat_model == "neutral_roughness_log";
    const bool vapor_neutral = contract.vapor_model_specified &&
        contract.vapor_model == "neutral_roughness_log";
    const bool any_neutral =
        momentum_neutral || heat_neutral || vapor_neutral;

    const bool any_channel_key =
        contract.momentum_model_specified ||
        contract.heat_model_specified ||
        contract.vapor_model_specified ||
        contract.coefficient_source_specified ||
        contract.heat_coefficient_specified ||
        contract.vapor_coefficient_specified ||
        contract.momentum_coefficient_specified ||
        contract.z0_m_specified ||
        contract.z0_h_specified ||
        contract.z0_q_specified;
    if (contract.legacy_aggregate_specified && any_channel_key) {
        return "Cloud Chamber: " + face +
            " cannot combine wall_transfer_model with per-channel wall model keys";
    }

    if (contract.momentum_model_specified &&
        contract.momentum_model != "resolved_noslip" &&
        !momentum_bulk && !momentum_neutral) {
        return "Cloud Chamber: " + face +
            ".momentum_transfer_model must be resolved_noslip, bulk_aero, bulk_aero_plus_molecular, "
            "neutral_roughness_log, or law_of_wall_momentum";
    }

    const bool heat_bulk = contract.heat_model_specified &&
        is_bulk_model(contract.heat_model);
    const bool vapor_bulk = contract.vapor_model_specified &&
        is_bulk_model(contract.vapor_model);
    if (contract.heat_model_specified &&
        contract.heat_model != "resolved_molecular" && !heat_bulk && !heat_neutral) {
        return "Cloud Chamber: " + face +
            ".heat_transfer_model must be resolved_molecular, bulk_aero, bulk_aero_plus_molecular, or neutral_roughness_log";
    }
    if (contract.vapor_model_specified &&
        contract.vapor_model != "resolved_molecular" && !vapor_bulk && !vapor_neutral) {
        return "Cloud Chamber: " + face +
            ".vapor_transfer_model must be resolved_molecular, bulk_aero, bulk_aero_plus_molecular, or neutral_roughness_log";
    }

    const bool bulk_active = heat_bulk || vapor_bulk || momentum_bulk;
    if ((contract.momentum_model == "bulk_aero_plus_molecular" ||
         contract.heat_model == "bulk_aero_plus_molecular" ||
         contract.vapor_model == "bulk_aero_plus_molecular") && contract.coefficient_source != "fixed") {
        return "Cloud Chamber: " + face + " bulk_aero_plus_molecular requires coefficient_source = fixed";
    }
    const bool most_provider = contract.coefficient_source_specified &&
        contract.coefficient_source == "most";
    if (most_provider && !(face == "zlo" || face == "zhi")) {
        return "Cloud Chamber: " + face +
            " supports coefficient_source = most only on zlo or zhi";
    }
    if ((any_neutral || (bulk_active && most_provider)) &&
        !contract.z0_m_specified) {
        return "Cloud Chamber: " + face +
            " requires z0_m for the active neutral/MOST wall-transfer channels";
    }
    if ((heat_neutral || (heat_bulk && most_provider)) &&
        !contract.z0_h_specified) {
        return "Cloud Chamber: " + face +
            " requires z0_h for neutral roughness-log heat or MOST heat";
    }
    if ((vapor_neutral || (vapor_bulk && most_provider)) &&
        !contract.z0_q_specified) {
        return "Cloud Chamber: " + face +
            " requires z0_q for neutral roughness-log vapor or MOST vapor";
    }
    if (contract.z0_m_specified &&
        !(any_neutral || (bulk_active && most_provider))) {
        return "Cloud Chamber: " + face +
            ".z0_m requires a neutral_roughness_log or MOST channel";
    }
    if (contract.z0_h_specified &&
        !(heat_neutral || (heat_bulk && most_provider))) {
        return "Cloud Chamber: " + face +
            ".z0_h requires neutral_roughness_log heat or MOST heat";
    }
    if (contract.z0_q_specified &&
        !(vapor_neutral || (vapor_bulk && most_provider))) {
        return "Cloud Chamber: " + face +
            ".z0_q requires neutral_roughness_log vapor or MOST vapor";
    }
    if ((contract.z0_m_specified &&
         (!finite(contract.z0_m) || contract.z0_m <= amrex::Real(0.0))) ||
        (contract.z0_h_specified &&
         (!finite(contract.z0_h) || contract.z0_h <= amrex::Real(0.0))) ||
        (contract.z0_q_specified &&
         (!finite(contract.z0_q) || contract.z0_q <= amrex::Real(0.0)))) {
        return "Cloud Chamber: " + face +
            " active roughness lengths z0_m, z0_h, and z0_q must be finite and positive";
    }

    if (bulk_active) {
        if (!contract.coefficient_source_specified ||
            (contract.coefficient_source != "fixed" &&
             contract.coefficient_source != "most")) {
            return "Cloud Chamber: " + face +
                " active bulk_aero channels require coefficient_source = fixed or most";
        }
    } else if (contract.coefficient_source_specified) {
        return "Cloud Chamber: " + face +
            ".coefficient_source requires a bulk_aero or bulk_aero_plus_molecular channel";
    }

    if (momentum_bulk && !most_provider &&
        !contract.momentum_coefficient_specified) {
        return "Cloud Chamber: " + face +
            " requires C_D when momentum_transfer_model = " + contract.momentum_model +
            " with coefficient_source = fixed";
    }
    if (momentum_bulk && most_provider && contract.momentum_coefficient_specified) {
        return "Cloud Chamber: " + face +
            ".C_D is not valid with coefficient_source = most";
    }
    if (heat_bulk && !most_provider && !contract.heat_coefficient_specified) {
        return "Cloud Chamber: " + face +
            " requires C_H when heat_transfer_model = " + contract.heat_model +
            " with coefficient_source = fixed";
    }
    if (heat_bulk && most_provider && contract.heat_coefficient_specified) {
        return "Cloud Chamber: " + face +
            ".C_H is not valid with coefficient_source = most";
    }
    if (vapor_bulk && !most_provider && !contract.vapor_coefficient_specified) {
        return "Cloud Chamber: " + face +
            " requires C_E when vapor_transfer_model = " + contract.vapor_model +
            " with coefficient_source = fixed";
    }
    if (vapor_bulk && most_provider && contract.vapor_coefficient_specified) {
        return "Cloud Chamber: " + face +
            ".C_E is not valid with coefficient_source = most";
    }
    if (contract.momentum_coefficient_specified && !momentum_bulk) {
        return "Cloud Chamber: " + face +
            ".C_D requires momentum_transfer_model = bulk_aero or bulk_aero_plus_molecular";
    }
    if (contract.heat_coefficient_specified && !heat_bulk) {
        return "Cloud Chamber: " + face +
            ".C_H requires heat_transfer_model = bulk_aero or bulk_aero_plus_molecular";
    }
    if (contract.vapor_coefficient_specified && !vapor_bulk) {
        return "Cloud Chamber: " + face +
            ".C_E requires vapor_transfer_model = bulk_aero or bulk_aero_plus_molecular";
    }
    if ((contract.momentum_coefficient_specified &&
         (!finite(contract.momentum_coefficient) ||
          contract.momentum_coefficient < amrex::Real(0.0))) ||
        (contract.heat_coefficient_specified &&
         (!finite(contract.heat_coefficient) ||
          contract.heat_coefficient < amrex::Real(0.0))) ||
        (contract.vapor_coefficient_specified &&
         (!finite(contract.vapor_coefficient) ||
          contract.vapor_coefficient < amrex::Real(0.0)))) {
        return "Cloud Chamber: " + face +
            " bulk coefficients C_D, C_H, and C_E must be finite and nonnegative";
    }
    return {};
}

// Wet-wall saturation needs T_wall independently of sensible heat transfer.
// A dry adiabatic wall does not need T_wall, even for MOST drag.
inline std::string
parse_wall_heat (amrex::ParmParse& pp, const std::string& face, bool physical,
                 erf_wall_thermodynamics::FaceWall& wall, const amrex::Geometry& geom,
                 SidewallTemperatureBands* sidewall_bands)
{
    using namespace erf_wall_thermodynamics;
    const bool edges = pp.contains("temperature_height_edges"), values = pp.contains("temperature_values");
    if (edges || values) {
        const auto prefix = "Cloud Chamber: " + face + ": ";
        if (!physical || !sidewall_bands) { return prefix + "temperature bands require a physical x or y sidewall"; }
        if (!(edges && values)) { return prefix + "temperature_height_edges and temperature_values must be supplied together"; }
        if (pp.contains("temperature")) { return prefix + "temperature and temperature bands are mutually exclusive"; }
    }
    if (!physical) {
        if (pp.contains("heat") || pp.contains("thermal")) {
            return "Cloud Chamber: " + face + ".heat and its .thermal alias require physical_temperature_rh";
        }
        return {};
    }

    std::string heat = "conductive";
    const bool explicit_heat = pp.query("heat", heat);
    heat = amrex::toLower(heat);
    if (heat != "conductive" && heat != "adiabatic") {
        return "Cloud Chamber: " + face + ".heat must be conductive or adiabatic";
    }
    std::string alias;
    if (pp.query("thermal", alias)) {
        alias = amrex::toLower(alias);
        if (alias == "fixed_temperature") { alias = "conductive"; }
        else if (alias != "adiabatic") {
            return "Cloud Chamber: " + face + ".thermal must be fixed_temperature or adiabatic";
        }
        if (explicit_heat && heat != alias) {
            return "Cloud Chamber: conflicting " + face + ".heat and " + face + ".thermal";
        }
        heat = alias;
    }

    const bool conductive = heat == "conductive";
    ThermalBoundary thermal;
    thermal.mode = conductive ? ThermalMode::FixedPhysicalTemperature : ThermalMode::Adiabatic;
    if (edges) {
        amrex::Vector<amrex::Real> z, temperatures;
        pp.getarr("temperature_height_edges", z);
        pp.getarr("temperature_values", temperatures);
        const auto error = sidewall_bands->define(std::move(z), temperatures, geom);
        if (!error.empty()) { return "Cloud Chamber: " + face + ": " + error; }
        // Uniform bands use the scalar path. For nonuniform bands each active
        // scalar flux uses its panel temperature; MOST cannot run on this face.
        thermal.temperature_K = temperatures.front();
    } else if (pp.query("temperature", thermal.temperature_K)) {
        if (!finite(thermal.temperature_K) || thermal.temperature_K <= 0) {
            return "Cloud Chamber: " + face + ".temperature must be finite and positive Kelvin";
        }
    } else if (conductive || wall.vapor_is_wet()) {
        return "Cloud Chamber: missing required physical " + face +
            ".temperature (needed for conductive heat or wet-wall saturation)";
    }
    wall.thermal = thermal;
    return {};
}

inline std::string
wall_coupled_most_roughness_error (const erf_wall_thermodynamics::FaceWall& wall,
                                   const std::string& face)
{
    if (wall.heat_uses_most() && wall.vapor_uses_most() && wall.heat.z0 != wall.vapor.z0) {
        return "Cloud Chamber: " + face +
            " coupled heat/vapor MOST requires z0_h == z0_q because the "
            "prescribed-temperature bulk-Richardson closure uses one scalar "
            "similarity resistance; z0_h=" +
            std::to_string(static_cast<double>(wall.heat.z0)) +
            ", z0_q=" +
            std::to_string(static_cast<double>(wall.vapor.z0));
    }
    return {};
}

inline void require_zero_periodicity ()
{
    amrex::ParmParse pp_geom("geometry");
    std::vector<int> periodic;
    pp_geom.getarr("is_periodic", periodic);
    if (periodic.size() != AMREX_SPACEDIM) {
        amrex::Error("Cloud Chamber: geometry.is_periodic must provide three values");
    }
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (periodic[d] != 0) {
            amrex::Error("Cloud Chamber: geometry.is_periodic must be 0 0 0");
        }
    }
}

inline void require_no_unsupported_geometry ()
{
    amrex::ParmParse pp_erf("erf");
    std::string value;
    if (pp_erf.query("terrain_type", value) &&
        amrex::toLower(value) != "none") {
        amrex::Error("Cloud Chamber: erf.terrain_type must be None; terrain and embedded boundaries are out of scope");
    }
    value.clear();
    if (pp_erf.query("buildings_type", value) &&
        amrex::toLower(value) != "none") {
        amrex::Error("Cloud Chamber: erf.buildings_type must be None; immersed buildings are out of scope");
    }
    value.clear();
    if (pp_erf.query("mesh_type", value) &&
        amrex::toLower(value) != "constantdz") {
        amrex::Error("Cloud Chamber: erf.mesh_type must be ConstantDz");
    }
    // A stretched or listed vertical grid makes the mesh StretchedDz without erf.mesh_type
    if (SolverChoice::mesh_type != MeshType::ConstantDz) {
        amrex::Error("Cloud Chamber: the vertical grid must be uniform (no erf.grid_stretching_ratio, "
                     "erf.initial_dz or erf.terrain_z_levels)");
    }
}

inline void reject_unsupported_face_keys (amrex::ParmParse& pp,
                                          const std::string& face)
{
    if (pp.contains("physical_temperature") || pp.contains("relative_humidity") ||
        pp.contains("rh") || pp.contains("roughness") ||
        pp.contains("roughness_length") || pp.contains("z0") ||
        pp.contains("most") || pp.contains("MOST")) {
        amrex::Error("Cloud Chamber: " + face +
                     " must use temperature and moisture for physical walls; "
                     "roughness aliases and legacy aliases are unsupported");
    }
}

inline std::string
wall_roughness_geometry_error (
    const erf_wall_thermodynamics::Boundary& walls,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx)
{
    for (int n = 0; n < 2 * AMREX_SPACEDIM; ++n) {
        const auto& wall = walls[n];
        const int dir = n / 2;
        const amrex::Real z_ref = amrex::Real(0.5) * dx[dir];
        const std::string face = face_name(static_cast<WallFace>(n));
        const auto check = [&](bool active, const char* key,
                               amrex::Real z0) -> std::string {
            if (!active) { return {}; }
            if (!finite(z_ref) || !finite(z0) ||
                !(z0 > amrex::Real(0.0)) || !(z_ref > z0)) {
                return "Cloud Chamber: " + face + " " + key +
                    " requires 0 < " + key + " < z_ref; " + key + "=" +
                    std::to_string(static_cast<double>(z0)) +
                    ", z_ref=" +
                    std::to_string(static_cast<double>(z_ref));
            }
            return {};
        };
        // An adiabatic face never evaluates its heat resistance. Keep the
        // existing geometry contract for configured vapor models, even if dry.
        using namespace erf_wall_thermodynamics;
        const bool heat_rough = wall.heat_uses_most() ||
            (wall.heat_is_conductive() && wall.heat.model == ScalarModel::NeutralRoughnessLog);
        const bool vapor_rough = wall.vapor.model == ScalarModel::NeutralRoughnessLog ||
            (wall.vapor.model == ScalarModel::BulkAero &&
             wall.vapor.provider == CoefficientProvider::MOST);
        const auto momentum_error = check(
            wall.momentum.model == MomentumModel::NeutralRoughnessLog ||
            wall.momentum_uses_most() || heat_rough || vapor_rough,
            "z0_m", wall.momentum.z0_m);
        if (!momentum_error.empty()) { return momentum_error; }
        const auto heat_error = check(
            heat_rough,
            "z0_h", wall.heat.z0);
        if (!heat_error.empty()) { return heat_error; }
        const auto vapor_error = check(
            vapor_rough,
            "z0_q", wall.vapor.z0);
        if (!vapor_error.empty()) { return vapor_error; }
    }
    return {};
}

inline void
validate_wall_roughness_geometry (
    const erf_wall_thermodynamics::Boundary& walls,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx)
{
    const auto error = wall_roughness_geometry_error(walls, dx);
    if (!error.empty()) { amrex::Error(error); }
}

inline Initialization parse_initialization (const amrex::Real* problo,
                                            const amrex::Real* probhi,
                                            const MoistureType selected_moisture_model)
{
    Initialization config;
    config.active = true;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        config.prob_lo[d] = problo[d];
        config.prob_hi[d] = probhi[d];
        require_finite("geometry extent", config.prob_lo[d]);
        require_finite("geometry extent", config.prob_hi[d]);
        if (!(config.prob_hi[d] > config.prob_lo[d])) {
            amrex::Error("Cloud Chamber: geometry.prob_hi must exceed geometry.prob_lo in every direction");
        }
    }

    amrex::ParmParse pp_erf("erf");
    CompressibilityType compressibility_type = CompressibilityType::FullyCompressible;
    int anelastic = 0;
    if (pp_erf.query("anelastic", anelastic)) {
        if (anelastic == 1) { compressibility_type = CompressibilityType::Anelastic; }
    } else {
        pp_erf.query_enum_case_insensitive("compressibility_type", compressibility_type);
    }
    // The chamber initializer corrects only rho theta and the moisture at the fixed base-state
    // density, so it needs the anelastic equations; PseudoIncompressible is not supported
    if (compressibility_type != CompressibilityType::Anelastic) {
        amrex::Error("Cloud Chamber: erf.compressibility_type = Anelastic is required "
                     "(FullyCompressible and PseudoIncompressible are not supported)");
    }
    bool use_gravity = false;
    if (!pp_erf.query("use_gravity", use_gravity) || !use_gravity) {
        amrex::Error("Cloud Chamber: erf.use_gravity = true is required");
    }

    amrex::ParmParse pp_amr("amr");
    int max_level = 0;
    pp_amr.query("max_level", max_level);
    if (max_level != 0) {
        amrex::Error("Cloud Chamber: amr.max_level must be 0");
    }
    require_zero_periodicity();
    require_no_unsupported_geometry();

    amrex::ParmParse pp_prob("prob");
    amrex::Real p_inf = p_0;
    if (pp_prob.query("p_inf", p_inf)) {
        require_finite("prob.p_inf", p_inf);
        const amrex::Real pressure_tolerance = amrex::Real(1.0e-6) * p_0;
        if (std::abs(p_inf - p_0) > pressure_tolerance) {
            amrex::Error("Cloud Chamber: prob.p_inf must match ERF's p_0 = 100000 Pa within 1e-6 relative tolerance");
        }
    }

    if (!supports_moisture_model(selected_moisture_model)) {
        amrex::Error("Cloud Chamber: selected moisture model is unsupported; use an Eulerian model from the supported moisture allowlist");
    }
    config.has_moisture = selected_moisture_model != MoistureType::None;
    if (!config.has_moisture && pp_erf.contains("moisture_model")) {
        amrex::Error("Cloud Chamber: omit erf.moisture_model for dry mode instead of explicitly selecting None");
    }
    if (!config.has_moisture &&
        (pp_prob.contains("qv_bottom") || pp_prob.contains("qv_top"))) {
        amrex::Error("Cloud Chamber: qv_bottom/qv_top require an active moisture model");
    }

    std::string init_mode;
    pp_prob.get("thermodynamic_initialization", init_mode);
    init_mode = amrex::toLower(init_mode);
    if (init_mode == "physical_temperature_rh") {
        config.physical_initialization = true;
        std::string init_type;
        pp_erf.get("init_type", init_type);
        if (amrex::toLower(init_type) != "constantdensity") {
            amrex::Error("Cloud Chamber: physical initialization requires the existing HSE-backed erf.init_type = ConstantDensity path");
        }
        bool vert_implicit = false;
        pp_erf.query("vert_implicit", vert_implicit);
        if (vert_implicit) {
            amrex::Error("Cloud Chamber: physical wall transfer requires erf.vert_implicit = false");
        }
        pp_prob.get("initial_temperature_bottom", config.initial_temperature_bottom);
        pp_prob.get("initial_temperature_top", config.initial_temperature_top);
        const bool has_initial_relative_humidity = pp_prob.contains("initial_relative_humidity");
        auto physical_contract_error = [&](bool has_legacy_profile_keys) {
            const InitializationContract contract {
                InitializationMode::PhysicalTemperatureRH,
                config.has_moisture,
                has_initial_relative_humidity,
                has_legacy_profile_keys,
                false};
            return initialization_contract_error(contract);
        };
        std::string contract_error = physical_contract_error(false);
        if (!contract_error.empty()) {
            amrex::Error(contract_error);
        }
        if (config.has_moisture) {
            pp_prob.get("initial_relative_humidity", config.initial_relative_humidity);
        }
        pp_prob.query("temperature_perturbation_amplitude",
                      config.temperature_perturbation_amplitude);
        const bool has_legacy_profile_keys = pp_prob.contains("theta_bottom") ||
            pp_prob.contains("theta_top") ||
            pp_prob.contains("theta_perturbation_amplitude") ||
            pp_prob.contains("qv_bottom") ||
            pp_prob.contains("qv_top");
        contract_error = physical_contract_error(has_legacy_profile_keys);
        if (!contract_error.empty()) {
            amrex::Error(contract_error);
        }
        require_finite("prob.initial_temperature_bottom", config.initial_temperature_bottom);
        require_finite("prob.initial_temperature_top", config.initial_temperature_top);
        require_finite("prob.temperature_perturbation_amplitude",
                       config.temperature_perturbation_amplitude);
        if (config.initial_temperature_bottom <= amrex::Real(0.0) ||
            config.initial_temperature_top <= amrex::Real(0.0)) {
            amrex::Error("Cloud Chamber: physical initial temperatures must be positive");
        }
        if (config.has_moisture) {
            require_finite("prob.initial_relative_humidity", config.initial_relative_humidity);
            if (config.initial_relative_humidity < amrex::Real(0.0) ||
                config.initial_relative_humidity > amrex::Real(1.0)) {
                amrex::Error("Cloud Chamber: initial_relative_humidity must lie in [0,1]");
            }
        }
    } else if (init_mode == "legacy_theta_qv") {
        config.physical_initialization = false;
        const InitializationContract contract {
            InitializationMode::LegacyThetaQv,
            false,
            false,
            false,
            pp_prob.contains("initial_temperature_bottom") ||
                pp_prob.contains("initial_temperature_top") ||
                pp_prob.contains("initial_relative_humidity") ||
                pp_prob.contains("temperature_perturbation_amplitude")};
        const std::string contract_error = initialization_contract_error(contract);
        if (!contract_error.empty()) {
            amrex::Error(contract_error);
        }
        pp_prob.get("theta_bottom", config.theta_bottom);
        pp_prob.get("theta_top", config.theta_top);
        pp_prob.query("theta_perturbation_amplitude",
                      config.theta_perturbation_amplitude);
    } else {
        amrex::Error("Cloud Chamber: prob.thermodynamic_initialization must be physical_temperature_rh or legacy_theta_qv");
    }
    std::string perturbation_mode = "deterministic_sine";
    pp_prob.query("perturbation_mode", perturbation_mode);
    if (amrex::toLower(perturbation_mode) != "deterministic_sine") {
        amrex::Error("Cloud Chamber: prob.perturbation_mode must be deterministic_sine");
    }
    if (!config.physical_initialization) {
        require_finite("prob.theta_bottom", config.theta_bottom);
        require_finite("prob.theta_top", config.theta_top);
        require_finite("prob.theta_perturbation_amplitude",
                       config.theta_perturbation_amplitude);
        if (config.theta_bottom <= amrex::Real(0.0) ||
            config.theta_top <= amrex::Real(0.0)) {
            amrex::Error("Cloud Chamber: theta profile values must be positive");
        }
    }

    if (config.has_moisture) {
        if (!config.physical_initialization) {
            pp_prob.get("qv_bottom", config.qv_bottom);
            pp_prob.get("qv_top", config.qv_top);
            require_finite("prob.qv_bottom", config.qv_bottom);
            require_finite("prob.qv_top", config.qv_top);
            if (config.qv_bottom < amrex::Real(0.0) ||
                config.qv_top < amrex::Real(0.0)) {
                amrex::Error("Cloud Chamber: qv profile values must be nonnegative");
            }
        }
    }

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const amrex::Real length = config.prob_hi[d] - config.prob_lo[d];
        if (!(length > amrex::Real(0.0)) || !finite(length)) {
            amrex::Error("Cloud Chamber: every rectangular dimension must be finite and positive");
        }
    }
    return config;
}

// Set up one physical face using the already parsed ConstantAlpha solver choice (no SGS).
// Validate a selected combined channel even on an impermeable wall, just as for C_E;
// only an active channel gets a speed. Ordinary bulk/resolved inputs retain their defaults.
inline std::string set_wall_transport_speeds (
    erf_wall_thermodynamics::FaceWall& wall, const WallTransferContract& contract,
    const DiffChoice& dc, amrex::Real dx_inv, const std::string& face)
{
    using namespace erf_wall_thermodynamics;
    const auto positive = [](amrex::Real x) { return finite(x) && x > 0; };
    const auto error = [&](const char* channel, const char* key) {
        return "Cloud Chamber: " + face + "." + channel + "_transfer_model = "
            "bulk_aero_plus_molecular requires finite positive " + key;
    };
    const auto inv_distance = amrex::Real(2) * dx_inv;
    if (contract.momentum_model == "bulk_aero_plus_molecular") {
        if (!positive(dc.dynamic_viscosity)) { return error("momentum", "dynamic_viscosity"); }
        if (!positive(dc.rho0_trans)) { return error("momentum", "rho0_trans"); }
        const auto speed = (dc.dynamic_viscosity / dc.rho0_trans) * inv_distance;
        if (!positive(speed)) {
            return error("momentum", "(dynamic_viscosity/rho0_trans) / half-cell distance");
        }
        wall.momentum.molecular_speed = speed;
    } else if (wall.momentum.model == MomentumModel::ResolvedNoSlip &&
               positive(dc.dynamic_viscosity) && positive(dc.rho0_trans)) {
        // ComputeStrain_N uses a prescribed wall value and two interior velocities.
        // Only the latter vary in the wall-stress row; the normal velocity is prescribed too.
        const auto w = WallGradientCoeffsUniform();
        wall.momentum.resolved_row_rate = (std::abs(w.c_0) + std::abs(w.c_1)) *
            (dc.dynamic_viscosity / dc.rho0_trans) * dx_inv * dx_inv;
    }
    const auto scalar = [&](ScalarTransfer& transfer, const std::string& model, bool active,
                            amrex::Real alpha, const char* channel, const char* key) -> std::string {
        if (model == "bulk_aero_plus_molecular") {
            if (!positive(alpha)) { return error(channel, key); }
            if (active) {
                const auto speed = alpha * inv_distance;
                if (!positive(speed)) { return error(channel, "molecular transfer speed"); }
                transfer.molecular_speed = speed;
            }
        } else if (active && transfer.model == ScalarModel::ResolvedMolecular) {
            // Guard metadata only: resolved flux already includes this molecular exchange.
            transfer.resolved_speed = alpha * inv_distance;
        }
        return {};
    };
    const auto heat_error = scalar(wall.heat, contract.heat_model,
        wall.heat_is_conductive(), dc.alpha_T, "heat", "alpha_T");
    if (!heat_error.empty()) { return heat_error; }
    return scalar(wall.vapor, contract.vapor_model,
        wall.vapor_is_wet(), dc.alpha_C, "vapor", "alpha_C");
}

inline Config parse_config (const amrex::Geometry& geom, MoistureType selected_moisture_model,
                            const DiffChoice& dc)
{
    Config config;
    static_cast<Initialization&>(config) = parse_initialization(
        geom.ProbLo(), geom.ProbHi(), selected_moisture_model);
    if (config.physical_initialization && dc.molec_diff_type != MolecDiffType::ConstantAlpha) {
        amrex::Error("Cloud Chamber: physical wall transfer requires erf.molec_diff_type = ConstantAlpha");
    }
    // max_level == 0 and ConstantDz were checked above. This is the only geometry
    // from which the RHS, timestep and roughness validation get their face metadata.
    erf_wall_thermodynamics::Boundary walls{};
    const auto dx_inv = geom.InvCellSizeArray();
    const WallFace faces[] = {WallFace::XLo, WallFace::XHi, WallFace::YLo,
                              WallFace::YHi, WallFace::ZLo, WallFace::ZHi};
    for (int n = 0; n < 2 * AMREX_SPACEDIM; ++n) {
        const std::string face = face_name(faces[n]);
        amrex::ParmParse pp_face(face.c_str());
        reject_unsupported_face_keys(pp_face, face);
        std::string type;
        if (!pp_face.query("type", type) ||
            amrex::toLower(type) != "noslipwall") {
            amrex::Error("Cloud Chamber: " + face + ".type must be NoSlipWall");
        }
        auto& wall = walls[n];
        auto transfer_contract = read_wall_transfer_contract(pp_face, face);
        const std::string transfer_error =
            wall_transfer_contract_error(transfer_contract, face);
        if (!transfer_error.empty()) {
            amrex::Error(transfer_error);
        }
        const bool any_channel_key =
            transfer_contract.momentum_model_specified ||
            transfer_contract.heat_model_specified ||
            transfer_contract.vapor_model_specified ||
            transfer_contract.coefficient_source_specified ||
            transfer_contract.momentum_coefficient_specified ||
            transfer_contract.heat_coefficient_specified ||
            transfer_contract.vapor_coefficient_specified ||
            transfer_contract.z0_m_specified ||
            transfer_contract.z0_h_specified ||
            transfer_contract.z0_q_specified;
        if (!config.physical_initialization && any_channel_key) {
            amrex::Error("Cloud Chamber: per-channel wall transfer keys require physical temperature initialization");
        }
        if (config.physical_initialization) {
            if (pp_face.contains("theta") || pp_face.contains("qv") ||
                pp_face.contains("theta_grad") || pp_face.contains("qv_grad")) {
                amrex::Error("Cloud Chamber: physical walls cannot combine temperature/moisture with legacy theta/qv wall keys");
            }
            std::string moisture;
            pp_face.get("moisture", moisture);
            moisture = amrex::toLower(moisture);
            if (moisture == "dry") {
                wall.moisture = erf_wall_thermodynamics::MoistureMode::DryImpermeable;
            } else if (moisture == "wet") {
                wall.moisture = erf_wall_thermodynamics::MoistureMode::WetEquilibrium;
                if (!config.has_moisture) {
                    amrex::Error("Cloud Chamber: wet walls require an active moisture model");
                }
            } else {
                amrex::Error("Cloud Chamber: " + face + ".moisture must be dry or wet");
            }
            if (transfer_contract.momentum_model == "law_of_wall_momentum") {
                transfer_contract.momentum_model = "neutral_roughness_log";
            }
            wall.momentum.model =
                transfer_contract.momentum_model == "neutral_roughness_log" ?
                erf_wall_thermodynamics::MomentumModel::NeutralRoughnessLog :
                (is_bulk_model(transfer_contract.momentum_model) ?
                 erf_wall_thermodynamics::MomentumModel::BulkAero :
                 erf_wall_thermodynamics::MomentumModel::ResolvedNoSlip);
            wall.heat.model =
                is_bulk_model(transfer_contract.heat_model) ?
                erf_wall_thermodynamics::ScalarModel::BulkAero :
                (transfer_contract.heat_model == "neutral_roughness_log" ?
                 erf_wall_thermodynamics::ScalarModel::NeutralRoughnessLog :
                 erf_wall_thermodynamics::ScalarModel::ResolvedMolecular);
            wall.vapor.model =
                is_bulk_model(transfer_contract.vapor_model) ?
                erf_wall_thermodynamics::ScalarModel::BulkAero :
                (transfer_contract.vapor_model == "neutral_roughness_log" ?
                 erf_wall_thermodynamics::ScalarModel::NeutralRoughnessLog :
                 erf_wall_thermodynamics::ScalarModel::ResolvedMolecular);
            wall.momentum.z0_m = transfer_contract.z0_m;
            wall.momentum.C_D = transfer_contract.momentum_coefficient;
            wall.heat.z0 = transfer_contract.z0_h;
            wall.vapor.z0 = transfer_contract.z0_q;
            if (wall.momentum.model == erf_wall_thermodynamics::MomentumModel::BulkAero) {
                wall.momentum.provider = transfer_contract.coefficient_source == "most" ?
                    erf_wall_thermodynamics::CoefficientProvider::MOST :
                    erf_wall_thermodynamics::CoefficientProvider::Fixed;
            }
            if (wall.heat.model == erf_wall_thermodynamics::ScalarModel::BulkAero) {
                wall.heat.provider = transfer_contract.coefficient_source == "most" ?
                    erf_wall_thermodynamics::CoefficientProvider::MOST :
                    erf_wall_thermodynamics::CoefficientProvider::Fixed;
                wall.heat.coefficient = transfer_contract.heat_coefficient;
            }
            if (wall.vapor.model == erf_wall_thermodynamics::ScalarModel::BulkAero) {
                wall.vapor.provider = transfer_contract.coefficient_source == "most" ?
                    erf_wall_thermodynamics::CoefficientProvider::MOST :
                    erf_wall_thermodynamics::CoefficientProvider::Fixed;
                wall.vapor.coefficient = transfer_contract.vapor_coefficient;
            }
        } else {
            if (!pp_face.query("theta", wall.thermal.temperature_K)) {
                amrex::Error("Cloud Chamber: missing required " + face + ".theta");
            }
            wall.thermal.mode = erf_wall_thermodynamics::ThermalMode::LegacyNumerical;
            wall.moisture = config.has_moisture ?
                erf_wall_thermodynamics::MoistureMode::LegacyNumerical :
                erf_wall_thermodynamics::MoistureMode::DryImpermeable;
            if (config.has_moisture && !pp_face.contains("qv")) {
                amrex::Error("Cloud Chamber: legacy cloudy walls require " + face + ".qv");
            }
        }
        const auto heat_error = parse_wall_heat(pp_face, face, config.physical_initialization, wall,
                            geom, n < 4 ? &config.m_bands[n] : nullptr);
        if (!heat_error.empty()) { amrex::Error(heat_error); }
        const auto coupled_error = wall_coupled_most_roughness_error(wall, face);
        if (!coupled_error.empty()) { amrex::Error(coupled_error); }
        if (config.physical_initialization) {
            // After parse_wall_heat: only conductive heat gets a molecular speed.
            const auto speed_error = set_wall_transport_speeds(
                wall, transfer_contract, dc, dx_inv[n/2], face);
            if (!speed_error.empty()) { amrex::Error(speed_error); }
        }
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (pp_face.contains("velocity")) {
                amrex::Error("Cloud Chamber: wall velocity metadata is not active; use NoSlipWall with zero velocity");
            }
            wall.wall_velocity[d] = amrex::Real(0.0);
        }
        const auto wetness_error = parse_wall_wetness(pp_face, face, wall);
        if (!wetness_error.empty()) { amrex::Error(wetness_error); }
    }
    validate_wall_roughness_geometry(walls, geom.CellSizeArray());
    config.m_boundary = walls;
    return config;
}

} // namespace erf_cloud_chamber

#endif
