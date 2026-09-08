/**
 * \file ERF_Plotfile2DSampledLevel.H
 */
#ifndef ERF_PLOTFILE2DSAMPLED_LEVEL_H_
#define ERF_PLOTFILE2DSAMPLED_LEVEL_H_

#include <optional>
#include <string>

#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include "ERF_Plotfile2DCatalog.H"
#include "ERF_Plotfile2DSampledField.H"

// User-facing sampled-level definitions for ERF 2D plotfiles. A level set
// names a vertical coordinate, target values, and fields to sample. This module
// parses and validates those definitions. It does not interpolate fields.

namespace plotfile2d
{

enum class SampledCoordinate
{
    ModelIndex,
    HeightMSL,
    HeightAGL,
    Pressure
};

enum class SampledInterpolation
{
    None,
    Linear
};

struct SampledVerticalCoordinateMetadata
{
    std::string type;
    amrex::Real value = amrex::Real(0.0);
    std::string units;
    amrex::Real canonical_value = amrex::Real(0.0);
    std::string canonical_units;
    std::string interpolation;
};

struct SampledLevelMetadata
{
    std::string level_set_name;
    std::string source_field;
    SampledVerticalCoordinateMetadata vertical_coordinate;
};

struct Plotfile2DOutputDescriptor
{
    std::string name;
    std::string long_name;
    std::string units;
    DiagnosticCategory category = DiagnosticCategory::SampledLevel;
    MissingPolicy missing_policy = MissingPolicy::FillMinus999WhenUnavailable;
    amrex::Real missing_value = amrex::Real(-999.0);
    const DiagnosticDescriptor* static_diagnostic = nullptr;
    std::optional<SampledLevelMetadata> sampled_level;
};

struct SampledLevelDefinition
{
    std::string name;
    SampledCoordinate coordinate = SampledCoordinate::ModelIndex;
    std::string units;
    amrex::Vector<amrex::Real> values;
    amrex::Vector<std::string> fields;
    SampledInterpolation interpolation = SampledInterpolation::None;
    amrex::Real missing_value = amrex::Real(-999.0);
};

const char* sampled_coordinate_to_string (SampledCoordinate coordinate) noexcept;
const char* sampled_coordinate_tag (SampledCoordinate coordinate) noexcept;
const char* sampled_coordinate_default_units (SampledCoordinate coordinate) noexcept;
const char* sampled_coordinate_default_interpolation (SampledCoordinate coordinate) noexcept;

const char* sampled_interpolation_to_string (SampledInterpolation interpolation) noexcept;
SampledCoordinate sampled_coordinate_from_string (const std::string& value);
SampledInterpolation sampled_interpolation_from_string (const std::string& value);
std::string validate_sampled_coordinate_string (const std::string& level_set_name,
                                                const std::string& value);

std::string sampled_level_error_prefix (const std::string& level_set_name,
                                       const std::string& parameter_name);

std::string sampled_level_value_tag (amrex::Real value, const std::string& units);
std::string sampled_output_name (const std::string& field_name,
                                SampledCoordinate coordinate,
                                amrex::Real value,
                                const std::string& units);

std::string validate_sampled_level_definition (const SampledLevelDefinition& level_set);

SampledLevelDefinition parse_sampled_level_definition (const std::string& level_set_name,
                                                      const std::string& pp_prefix);

amrex::Vector<std::string> parse_requested_sampled_level_sets (const std::string& pp_prefix,
                                                               int which);

amrex::Vector<Plotfile2DOutputDescriptor> build_sampled_level_output_descriptors_from_definitions (
    const amrex::Vector<SampledLevelDefinition>& level_sets,
    const amrex::Vector<std::string>& static_plot_vars,
    const SolverChoice& solver_choice);

amrex::Vector<Plotfile2DOutputDescriptor> build_sampled_level_output_descriptors (
    const std::string& pp_prefix,
    int which,
    const amrex::Vector<std::string>& static_plot_vars,
    const SolverChoice& solver_choice);

} // namespace plotfile2d

#endif
