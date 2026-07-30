#ifndef ERF_CLOUD_CHAMBER_H_
#define ERF_CLOUD_CHAMBER_H_

#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

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

enum class WallTreatment : unsigned char { PrescribedState };

struct WallSpec {
    WallFace face = WallFace::XLo;
    WallTreatment treatment = WallTreatment::PrescribedState;
    amrex::Real theta = amrex::Real(0.0);
    amrex::Real qv = amrex::Real(0.0);
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> velocity = {};
};

struct Config {
    bool active = false;
    bool cloudy = false;
    amrex::Real theta_bottom = amrex::Real(0.0);
    amrex::Real theta_top = amrex::Real(0.0);
    amrex::Real theta_perturbation_amplitude = amrex::Real(0.0);
    amrex::Real qv_bottom = amrex::Real(0.0);
    amrex::Real qv_top = amrex::Real(0.0);
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo = {};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_hi = {};
    amrex::GpuArray<WallSpec, 2 * AMREX_SPACEDIM> walls = {};
};

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
amrex::Real
linear_profile (amrex::Real bottom, amrex::Real top,
                amrex::Real coordinate, amrex::Real lower,
                amrex::Real length) noexcept
{
    return bottom + (top - bottom) * (coordinate - lower) / length;
}

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
    return linear_profile(config.theta_bottom, config.theta_top, z,
                          config.prob_lo[2], length[2]) +
           deterministic_perturbation(x, y, z, config.prob_lo, length,
                                      config.theta_perturbation_amplitude);
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
amrex::Real
qv_at (const Config& config, amrex::Real z) noexcept
{
    return linear_profile(config.qv_bottom, config.qv_top, z,
                          config.prob_lo[2],
                          config.prob_hi[2] - config.prob_lo[2]);
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
    for (const char* key : {"temperature", "physical_temperature",
                            "relative_humidity", "rh", "roughness"}) {
        if (pp.contains(key)) {
            amrex::Error("Cloud Chamber: " + face + "." + key +
                         " is unsupported; use prescribed numerical "
                         "theta and qv instead");
        }
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
    pp_prob.get("theta_bottom", config.theta_bottom);
    pp_prob.get("theta_top", config.theta_top);
    pp_prob.query("theta_perturbation_amplitude",
                  config.theta_perturbation_amplitude);
    std::string perturbation_mode = "deterministic_sine";
    pp_prob.query("perturbation_mode", perturbation_mode);
    if (amrex::toLower(perturbation_mode) != "deterministic_sine") {
        amrex::Error("Cloud Chamber: prob.perturbation_mode must be deterministic_sine");
    }
    require_finite("prob.theta_bottom", config.theta_bottom);
    require_finite("prob.theta_top", config.theta_top);
    require_finite("prob.theta_perturbation_amplitude",
                   config.theta_perturbation_amplitude);
    if (config.theta_bottom <= amrex::Real(0.0) ||
        config.theta_top <= amrex::Real(0.0)) {
        amrex::Error("Cloud Chamber: theta profile values must be positive");
    }

    std::string moisture_model;
    const bool moisture_specified = pp_erf.query("moisture_model", moisture_model);
    if (moisture_specified) {
        moisture_model = amrex::toLower(moisture_model);
        if (moisture_model != "satadj") {
            amrex::Error("Cloud Chamber: moisture_model must be omitted for dry mode or set to SatAdj");
        }
        config.cloudy = true;
        pp_prob.get("qv_bottom", config.qv_bottom);
        pp_prob.get("qv_top", config.qv_top);
        require_finite("prob.qv_bottom", config.qv_bottom);
        require_finite("prob.qv_top", config.qv_top);
        if (config.qv_bottom < amrex::Real(0.0) ||
            config.qv_top < amrex::Real(0.0)) {
            amrex::Error("Cloud Chamber: qv profile values must be nonnegative");
        }
    } else if (pp_prob.contains("qv_bottom") || pp_prob.contains("qv_top")) {
        amrex::Error("Cloud Chamber: qv_bottom/qv_top require erf.moisture_model = SatAdj");
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
        wall.treatment = WallTreatment::PrescribedState;
        if (!pp_face.query("theta", wall.theta)) {
            amrex::Error("Cloud Chamber: missing required " + face + ".theta");
        }
        require_finite(face + ".theta", wall.theta);
        if (wall.theta <= amrex::Real(0.0)) {
            amrex::Error("Cloud Chamber: " + face + ".theta must be positive");
        }
        if (config.cloudy) {
            if (!pp_face.query("qv", wall.qv)) {
                amrex::Error("Cloud Chamber: missing required " + face + ".qv for SatAdj");
            }
            require_finite(face + ".qv", wall.qv);
            if (wall.qv < amrex::Real(0.0)) {
                amrex::Error("Cloud Chamber: " + face + ".qv must be nonnegative");
            }
        } else if (pp_face.contains("qv")) {
            amrex::Error("Cloud Chamber: " + face + ".qv requires erf.moisture_model = SatAdj");
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
