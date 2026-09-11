#ifndef ERF_CLOUD_CHAMBER_H_
#define ERF_CLOUD_CHAMBER_H_

#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include "ERF_Constants.H"
#include "ERF_EOS.H"
#include "ERF_MicrophysicsUtils.H"
#include "BoundaryConditions/ERF_WallThermodynamics.H"

#include <cmath>
#include <limits>
#include <string>
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
 * @brief Metadata describing a chamber wall.
 */
struct WallSpec {
    WallFace face = WallFace::XLo;
    erf_wall_thermodynamics::FaceWall wall;
};

/**
 * @brief Configuration parameters for cloud chamber initialization.
 */
struct Config {
    bool active = false;
    bool cloudy = false;
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
    amrex::GpuArray<WallSpec, 2 * AMREX_SPACEDIM> walls = {};

    /**
     * @brief Collect thermodynamic boundary conditions for all chamber walls.
     * @return Boundary array containing thermodynamics for each face.
     */
    erf_wall_thermodynamics::Boundary wall_boundary () const noexcept
    {
        erf_wall_thermodynamics::Boundary boundary{};
        for (int n = 0; n < 2 * AMREX_SPACEDIM; ++n) {
            boundary[n] = walls[n].wall;
        }
        return boundary;
    }

    bool has_bulk_scalar_wall () const noexcept
    {
        for (const auto& wall : walls) {
            if (wall.wall.heat.model == erf_wall_thermodynamics::ScalarModel::BulkAero ||
                (wall.wall.vapor.model == erf_wall_thermodynamics::ScalarModel::BulkAero &&
                 wall.wall.moisture ==
                     erf_wall_thermodynamics::MoistureMode::WetEquilibrium)) {
                return true;
            }
        }
        return false;
    }

    bool has_wall_rate_channel () const noexcept
    {
        for (const auto& wall : walls) {
            if (wall.wall.momentum.model ==
                    erf_wall_thermodynamics::MomentumModel::BulkAero ||
                wall.wall.momentum.model ==
                    erf_wall_thermodynamics::MomentumModel::NeutralRoughnessLog ||
                wall.wall.heat.model ==
                    erf_wall_thermodynamics::ScalarModel::BulkAero ||
                wall.wall.heat.model ==
                    erf_wall_thermodynamics::ScalarModel::NeutralRoughnessLog ||
                (wall.wall.vapor.model ==
                     erf_wall_thermodynamics::ScalarModel::BulkAero &&
                 wall.wall.moisture ==
                     erf_wall_thermodynamics::MoistureMode::WetEquilibrium) ||
                (wall.wall.vapor.model ==
                     erf_wall_thermodynamics::ScalarModel::NeutralRoughnessLog &&
                 wall.wall.moisture ==
                     erf_wall_thermodynamics::MoistureMode::WetEquilibrium)) {
                return true;
            }
        }
        return false;
    }
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
    bool cloudy = false;
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
        if (contract.cloudy) {
            if (!contract.has_initial_relative_humidity) {
                return "Cloud Chamber: physical SatAdj initialization requires prob.initial_relative_humidity";
            }
        } else if (contract.has_initial_relative_humidity) {
            return "Cloud Chamber: prob.initial_relative_humidity is only used with erf.moisture_model = SatAdj";
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
theta_at (const Config& config, amrex::Real x, amrex::Real y,
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
qv_at (const Config& config, amrex::Real z) noexcept
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

inline std::string
wall_transfer_contract_error (const WallTransferContract& contract,
                              const std::string& face)
{
    const bool momentum_neutral = contract.momentum_model_specified &&
        (contract.momentum_model == "neutral_roughness_log" ||
         contract.momentum_model == "law_of_wall_momentum");
    const bool momentum_bulk = contract.momentum_model_specified &&
        contract.momentum_model == "bulk_aero";
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
            ".momentum_transfer_model must be resolved_noslip, bulk_aero, "
            "neutral_roughness_log, or law_of_wall_momentum";
    }

    const bool heat_bulk = contract.heat_model_specified &&
        contract.heat_model == "bulk_aero";
    const bool vapor_bulk = contract.vapor_model_specified &&
        contract.vapor_model == "bulk_aero";
    if (contract.heat_model_specified &&
        contract.heat_model != "resolved_molecular" && !heat_bulk && !heat_neutral) {
        return "Cloud Chamber: " + face +
            ".heat_transfer_model must be resolved_molecular, bulk_aero, or neutral_roughness_log";
    }
    if (contract.vapor_model_specified &&
        contract.vapor_model != "resolved_molecular" && !vapor_bulk && !vapor_neutral) {
        return "Cloud Chamber: " + face +
            ".vapor_transfer_model must be resolved_molecular, bulk_aero, or neutral_roughness_log";
    }

    const bool bulk_active = heat_bulk || vapor_bulk || momentum_bulk;
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
            ".z0_m requires a neutral_roughness_log channel";
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
            ".coefficient_source is only valid when bulk_aero is active";
    }

    if (momentum_bulk && !most_provider &&
        !contract.momentum_coefficient_specified) {
        return "Cloud Chamber: " + face +
            " requires C_D when momentum_transfer_model = bulk_aero with coefficient_source = fixed";
    }
    if (momentum_bulk && most_provider && contract.momentum_coefficient_specified) {
        return "Cloud Chamber: " + face +
            ".C_D is not valid with coefficient_source = most";
    }
    if (heat_bulk && !most_provider && !contract.heat_coefficient_specified) {
        return "Cloud Chamber: " + face +
            " requires C_H when heat_transfer_model = bulk_aero";
    }
    if (heat_bulk && most_provider && contract.heat_coefficient_specified) {
        return "Cloud Chamber: " + face +
            ".C_H is not valid with coefficient_source = most";
    }
    if (vapor_bulk && !most_provider && !contract.vapor_coefficient_specified) {
        return "Cloud Chamber: " + face +
            " requires C_E when vapor_transfer_model = bulk_aero";
    }
    if (vapor_bulk && most_provider && contract.vapor_coefficient_specified) {
        return "Cloud Chamber: " + face +
            ".C_E is not valid with coefficient_source = most";
    }
    if (contract.momentum_coefficient_specified && !momentum_bulk) {
        return "Cloud Chamber: " + face +
            ".C_D requires momentum_transfer_model = bulk_aero";
    }
    if (contract.heat_coefficient_specified && !heat_bulk) {
        return "Cloud Chamber: " + face +
            ".C_H requires heat_transfer_model = bulk_aero";
    }
    if (contract.vapor_coefficient_specified && !vapor_bulk) {
        return "Cloud Chamber: " + face +
            ".C_E requires vapor_transfer_model = bulk_aero";
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
        const bool momentum_neutral =
            wall.momentum.model ==
                erf_wall_thermodynamics::MomentumModel::NeutralRoughnessLog;
        const bool heat_neutral =
            wall.heat.model ==
                erf_wall_thermodynamics::ScalarModel::NeutralRoughnessLog;
        const bool vapor_neutral =
            wall.vapor.model ==
                erf_wall_thermodynamics::ScalarModel::NeutralRoughnessLog;
        const bool any_neutral =
            momentum_neutral || heat_neutral || vapor_neutral;
        const bool most_momentum =
            wall.momentum.model == erf_wall_thermodynamics::MomentumModel::BulkAero &&
            wall.momentum.provider == erf_wall_thermodynamics::CoefficientProvider::MOST;
        const bool most_heat =
            wall.heat.model == erf_wall_thermodynamics::ScalarModel::BulkAero &&
            wall.heat.provider == erf_wall_thermodynamics::CoefficientProvider::MOST;
        const bool most_vapor =
            wall.vapor.model == erf_wall_thermodynamics::ScalarModel::BulkAero &&
            wall.vapor.provider == erf_wall_thermodynamics::CoefficientProvider::MOST;
        const auto momentum_error = check(
            any_neutral || most_momentum || most_heat || most_vapor,
            "z0_m", wall.momentum.z0_m);
        if (!momentum_error.empty()) { return momentum_error; }
        const auto heat_error = check(
            heat_neutral || most_heat,
            "z0_h", wall.heat.z0);
        if (!heat_error.empty()) { return heat_error; }
        const auto vapor_error = check(
            vapor_neutral || most_vapor,
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

inline Config parse_config (const amrex::Real* problo,
                            const amrex::Real* probhi)
{
    Config config;
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
    int anelastic = 0;
    if (!pp_erf.query("anelastic", anelastic) || anelastic != 1) {
        amrex::Error("Cloud Chamber: erf.anelastic = 1 is required");
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

    std::string moisture_model;
    const bool moisture_specified = pp_erf.query("moisture_model", moisture_model);
    if (moisture_specified) {
        moisture_model = amrex::toLower(moisture_model);
        if (moisture_model != "satadj") {
            amrex::Error("Cloud Chamber: moisture_model must be omitted for dry mode or set to SatAdj");
        }
        config.cloudy = true;
    } else if (pp_prob.contains("qv_bottom") || pp_prob.contains("qv_top")) {
        amrex::Error("Cloud Chamber: qv_bottom/qv_top require erf.moisture_model = SatAdj");
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
        std::string molec_diff_type;
        pp_erf.get("molec_diff_type", molec_diff_type);
        if (amrex::toLower(molec_diff_type) != "constantalpha") {
            amrex::Error("Cloud Chamber: physical wall transfer requires erf.molec_diff_type = ConstantAlpha");
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
                config.cloudy,
                has_initial_relative_humidity,
                has_legacy_profile_keys,
                false};
            return initialization_contract_error(contract);
        };
        std::string contract_error = physical_contract_error(false);
        if (!contract_error.empty()) {
            amrex::Error(contract_error);
        }
        if (config.cloudy) {
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
        if (config.cloudy) {
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

    if (moisture_specified) {
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
        WallSpec& wall = config.walls[n];
        wall.face = faces[n];
        WallTransferContract transfer_contract;
        std::string transfer = "resolved_molecular";
        if (pp_face.query("wall_transfer_model", transfer)) {
            transfer_contract.legacy_aggregate_specified = true;
            transfer = amrex::toLower(transfer);
            if (transfer != "resolved_molecular") {
                amrex::Error("Cloud Chamber: only resolved_molecular wall transfer is supported");
            }
        }
        if (pp_face.query("momentum_transfer_model", transfer_contract.momentum_model)) {
            transfer_contract.momentum_model_specified = true;
            transfer_contract.momentum_model = amrex::toLower(transfer_contract.momentum_model);
        }
        if (pp_face.query("heat_transfer_model", transfer_contract.heat_model)) {
            transfer_contract.heat_model_specified = true;
            transfer_contract.heat_model = amrex::toLower(transfer_contract.heat_model);
        }
        if (pp_face.query("vapor_transfer_model", transfer_contract.vapor_model)) {
            transfer_contract.vapor_model_specified = true;
            transfer_contract.vapor_model = amrex::toLower(transfer_contract.vapor_model);
        }
        if (pp_face.query("coefficient_source", transfer_contract.coefficient_source)) {
            transfer_contract.coefficient_source_specified = true;
            transfer_contract.coefficient_source = amrex::toLower(transfer_contract.coefficient_source);
        }
        if (pp_face.query("C_D", transfer_contract.momentum_coefficient)) {
            transfer_contract.momentum_coefficient_specified = true;
            require_finite(face + ".C_D", transfer_contract.momentum_coefficient);
        }
        if (pp_face.query("C_H", transfer_contract.heat_coefficient)) {
            transfer_contract.heat_coefficient_specified = true;
            require_finite(face + ".C_H", transfer_contract.heat_coefficient);
        }
        if (pp_face.query("C_E", transfer_contract.vapor_coefficient)) {
            transfer_contract.vapor_coefficient_specified = true;
            require_finite(face + ".C_E", transfer_contract.vapor_coefficient);
        }
        if (pp_face.query("z0_m", transfer_contract.z0_m)) {
            transfer_contract.z0_m_specified = true;
            require_finite(face + ".z0_m", transfer_contract.z0_m);
        }
        if (pp_face.query("z0_h", transfer_contract.z0_h)) {
            transfer_contract.z0_h_specified = true;
            require_finite(face + ".z0_h", transfer_contract.z0_h);
        }
        if (pp_face.query("z0_q", transfer_contract.z0_q)) {
            transfer_contract.z0_q_specified = true;
            require_finite(face + ".z0_q", transfer_contract.z0_q);
        }
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
            if (!pp_face.query("temperature", wall.wall.thermal.temperature_K)) {
                amrex::Error("Cloud Chamber: missing required physical " + face + ".temperature");
            }
            require_finite(face + ".temperature", wall.wall.thermal.temperature_K);
            if (wall.wall.thermal.temperature_K <= amrex::Real(0.0)) {
                amrex::Error("Cloud Chamber: " + face + ".temperature must be positive Kelvin");
            }
            std::string moisture;
            pp_face.get("moisture", moisture);
            moisture = amrex::toLower(moisture);
            if (moisture == "dry") {
                wall.wall.moisture = erf_wall_thermodynamics::MoistureMode::DryImpermeable;
            } else if (moisture == "wet") {
                wall.wall.moisture = erf_wall_thermodynamics::MoistureMode::WetEquilibrium;
                if (!config.cloudy) {
                    amrex::Error("Cloud Chamber: wet walls require SatAdj moisture");
                }
            } else {
                amrex::Error("Cloud Chamber: " + face + ".moisture must be dry or wet");
            }
            wall.wall.thermal.mode = erf_wall_thermodynamics::ThermalMode::FixedPhysicalTemperature;
            if (transfer_contract.momentum_model == "law_of_wall_momentum") {
                transfer_contract.momentum_model = "neutral_roughness_log";
            }
            wall.wall.momentum.model =
                transfer_contract.momentum_model == "neutral_roughness_log" ?
                erf_wall_thermodynamics::MomentumModel::NeutralRoughnessLog :
                (transfer_contract.momentum_model == "bulk_aero" ?
                 erf_wall_thermodynamics::MomentumModel::BulkAero :
                 erf_wall_thermodynamics::MomentumModel::ResolvedNoSlip);
            wall.wall.heat.model =
                transfer_contract.heat_model == "bulk_aero" ?
                erf_wall_thermodynamics::ScalarModel::BulkAero :
                (transfer_contract.heat_model == "neutral_roughness_log" ?
                 erf_wall_thermodynamics::ScalarModel::NeutralRoughnessLog :
                 erf_wall_thermodynamics::ScalarModel::ResolvedMolecular);
            wall.wall.vapor.model =
                transfer_contract.vapor_model == "bulk_aero" ?
                erf_wall_thermodynamics::ScalarModel::BulkAero :
                (transfer_contract.vapor_model == "neutral_roughness_log" ?
                 erf_wall_thermodynamics::ScalarModel::NeutralRoughnessLog :
                 erf_wall_thermodynamics::ScalarModel::ResolvedMolecular);
            wall.wall.momentum.z0_m = transfer_contract.z0_m;
            wall.wall.momentum.C_D = transfer_contract.momentum_coefficient;
            wall.wall.heat.z0 = transfer_contract.z0_h;
            wall.wall.vapor.z0 = transfer_contract.z0_q;
            if (wall.wall.momentum.model == erf_wall_thermodynamics::MomentumModel::BulkAero) {
                wall.wall.momentum.provider = transfer_contract.coefficient_source == "most" ?
                    erf_wall_thermodynamics::CoefficientProvider::MOST :
                    erf_wall_thermodynamics::CoefficientProvider::Fixed;
            }
            if (wall.wall.heat.model == erf_wall_thermodynamics::ScalarModel::BulkAero) {
                wall.wall.heat.provider = transfer_contract.coefficient_source == "most" ?
                    erf_wall_thermodynamics::CoefficientProvider::MOST :
                    erf_wall_thermodynamics::CoefficientProvider::Fixed;
                wall.wall.heat.coefficient = transfer_contract.heat_coefficient;
            }
            if (wall.wall.vapor.model == erf_wall_thermodynamics::ScalarModel::BulkAero) {
                wall.wall.vapor.provider = transfer_contract.coefficient_source == "most" ?
                    erf_wall_thermodynamics::CoefficientProvider::MOST :
                    erf_wall_thermodynamics::CoefficientProvider::Fixed;
                wall.wall.vapor.coefficient = transfer_contract.vapor_coefficient;
            }
        } else {
            if (!pp_face.query("theta", wall.wall.thermal.temperature_K)) {
                amrex::Error("Cloud Chamber: missing required " + face + ".theta");
            }
            wall.wall.thermal.mode = erf_wall_thermodynamics::ThermalMode::LegacyNumerical;
            wall.wall.moisture = config.cloudy ?
                erf_wall_thermodynamics::MoistureMode::LegacyNumerical :
                erf_wall_thermodynamics::MoistureMode::DryImpermeable;
            if (config.cloudy && !pp_face.contains("qv")) {
                amrex::Error("Cloud Chamber: legacy cloudy walls require " + face + ".qv");
            }
        }
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (pp_face.contains("velocity")) {
                amrex::Error("Cloud Chamber: wall velocity metadata is not active; use NoSlipWall with zero velocity");
            }
            wall.wall.wall_velocity[d] = amrex::Real(0.0);
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

} // namespace erf_cloud_chamber

#endif
