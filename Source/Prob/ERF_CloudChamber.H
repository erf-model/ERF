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
 * Stage 1 wall metadata.  This is deliberately a small, trivially copyable
 * description: it can be captured by an AMReX device kernel without virtual
 * dispatch.  PrescribedState is the only active treatment.  A future wall
 * factory can consume the face identity and return momentum, thermal-scalar,
 * and vapor fluxes without changing the chamber initializer.
 */
enum class WallFace : unsigned char { XLo, XHi, YLo, YHi, ZLo, ZHi };

/**
 * @brief Metadata describing a chamber wall.
 */
struct WallSpec {
    WallFace face = WallFace::XLo;
    erf_wall_thermodynamics::Face thermodynamics;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> velocity = {};
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
            boundary[n] = walls[n].thermodynamics;
        }
        return boundary;
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
        pp.contains("rh") || pp.contains("roughness")) {
        amrex::Error("Cloud Chamber: " + face +
                     " must use temperature and moisture for physical walls; "
                     "roughness and legacy aliases are unsupported");
    }
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
        if (config.physical_initialization) {
            if (pp_face.contains("theta") || pp_face.contains("qv") ||
                pp_face.contains("theta_grad") || pp_face.contains("qv_grad")) {
                amrex::Error("Cloud Chamber: physical walls cannot combine temperature/moisture with legacy theta/qv wall keys");
            }
            if (!pp_face.query("temperature", wall.thermodynamics.physical_temperature_K)) {
                amrex::Error("Cloud Chamber: missing required physical " + face + ".temperature");
            }
            require_finite(face + ".temperature", wall.thermodynamics.physical_temperature_K);
            if (wall.thermodynamics.physical_temperature_K <= amrex::Real(0.0)) {
                amrex::Error("Cloud Chamber: " + face + ".temperature must be positive Kelvin");
            }
            std::string moisture;
            pp_face.get("moisture", moisture);
            moisture = amrex::toLower(moisture);
            if (moisture == "dry") {
                wall.thermodynamics.moisture_mode = erf_wall_thermodynamics::MoistureMode::DryImpermeable;
            } else if (moisture == "wet") {
                wall.thermodynamics.moisture_mode = erf_wall_thermodynamics::MoistureMode::WetEquilibrium;
                if (!config.cloudy) {
                    amrex::Error("Cloud Chamber: wet walls require SatAdj moisture");
                }
            } else {
                amrex::Error("Cloud Chamber: " + face + ".moisture must be dry or wet");
            }
            std::string transfer = "resolved_molecular";
            pp_face.query("wall_transfer_model", transfer);
            transfer = amrex::toLower(transfer);
            if (transfer != "resolved_molecular") {
                amrex::Error("Cloud Chamber: only resolved_molecular wall transfer is supported");
            }
            wall.thermodynamics.thermal_mode = erf_wall_thermodynamics::ThermalMode::FixedPhysicalTemperature;
            wall.thermodynamics.transfer_model = erf_wall_thermodynamics::TransferModel::ResolvedMolecular;
        } else {
            if (!pp_face.query("theta", config.walls[n].thermodynamics.physical_temperature_K)) {
                amrex::Error("Cloud Chamber: missing required " + face + ".theta");
            }
            wall.thermodynamics.thermal_mode = erf_wall_thermodynamics::ThermalMode::LegacyNumerical;
            wall.thermodynamics.moisture_mode = config.cloudy ?
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
            wall.velocity[d] = amrex::Real(0.0);
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
