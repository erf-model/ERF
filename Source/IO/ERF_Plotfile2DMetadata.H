/**
 * \file ERF_Plotfile2DMetadata.H
 */
#ifndef ERF_PLOTFILE2DMETADATA_H_
#define ERF_PLOTFILE2DMETADATA_H_

#include <string>

#include <AMReX_Vector.H>

#include "ERF_Plotfile2DCatalog.H"
#include "ERF_Plotfile2DSampledLevel.H"

// Utilities for writing machine-readable metadata for native AMReX 2D
// plotfiles. This layer formats the existing 2D diagnostic catalog. It must not
// compute diagnostics or duplicate diagnostic metadata strings.

namespace plotfile2d
{

const char* diagnostic_category_to_string (DiagnosticCategory category) noexcept;
const char* missing_policy_to_string (MissingPolicy policy) noexcept;

std::string missing_value_json (MissingPolicy policy);

std::string escape_json_string (const std::string& value);

std::string format_2d_metadata_json (const amrex::Vector<Plotfile2DOutputDescriptor>& descriptors);
std::string format_2d_metadata_json (const amrex::Vector<std::string>& varnames);

std::string metadata_json_filename (const std::string& plotfilename);

void write_2d_metadata_json (const std::string& plotfilename,
                             const amrex::Vector<Plotfile2DOutputDescriptor>& descriptors);

void write_2d_metadata_json (const std::string& plotfilename,
                             const amrex::Vector<std::string>& varnames);

} // namespace plotfile2d

#endif
