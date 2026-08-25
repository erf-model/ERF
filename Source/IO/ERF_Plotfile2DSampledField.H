/**
 * \file ERF_Plotfile2DSampledField.H
 */
#ifndef ERF_PLOTFILE2DSAMPLED_FIELD_H_
#define ERF_PLOTFILE2DSAMPLED_FIELD_H_

#include <string>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Vector.H>

#include "ERF_DataStruct.H"
#include "ERF_EOS.H"
#include "ERF_TerrainMetrics.H"

// Field registry for sampled-level 2D diagnostics. Each field defines its
// public name, units, availability, and cell-centered value. Derived fields can
// be added here without changing the sampled-level parser.

namespace plotfile2d
{

enum class SampledFieldID
{
    Rho,
    Theta,
    Temp,
    Pressure,
    HeightMSL,
    HeightAGL,
    Qv,
    Qc,
    Qi,
    Qr,
    Qs,
    Qg,
    UEast,
    VNorth,
    W,
    WindSpeed,
    WindDir
};

struct SampledFieldDescriptor
{
    SampledFieldID id;
    const char* name;
    const char* long_name;
    const char* units;
};

struct SampledFieldSelection
{
    amrex::Vector<SampledFieldDescriptor> accepted;
    amrex::Vector<std::string> unavailable;
};

const amrex::Vector<SampledFieldDescriptor>& sampled_field_catalog ();

const SampledFieldDescriptor* find_sampled_field (const std::string& name);

amrex::Vector<std::string> available_sampled_field_names (const SolverChoice& solver_choice);

SampledFieldSelection select_requested_sampled_fields (const amrex::Vector<std::string>& requested,
                                                       const SolverChoice& solver_choice);

std::string sampled_field_id_to_string (SampledFieldID field_id);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool sampled_field_is_wind (SampledFieldID field_id) noexcept
{
    switch (field_id) {
    case SampledFieldID::UEast:
    case SampledFieldID::VNorth:
    case SampledFieldID::W:
    case SampledFieldID::WindSpeed:
    case SampledFieldID::WindDir:
        return true;
    default:
        return false;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool sampled_field_is_scalar_state (SampledFieldID field_id) noexcept
{
    return !sampled_field_is_wind(field_id);
}

// NOTE: p_hse_arr holds the base state pressure and is used for the pressure field only
//       when use_hse_pressure is true, i.e. when the run is anelastic -- see the comment
//       at SampledFieldID::Pressure below.  When use_hse_pressure is false the array is
//       not read, so the caller may pass any valid array.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real sampled_field_value (SampledFieldID field_id,
                                 const amrex::Array4<const amrex::Real>& cons_arr,
                                 const amrex::Array4<const amrex::Real>& z_phys_cc_arr,
                                 const amrex::Array4<const amrex::Real>& z_phys_nd_arr,
                                 const amrex::Array4<const amrex::Real>& p_hse_arr,
                                 bool have_z_phys_cc,
                                 bool use_hse_pressure,
                                 int i, int j, int k,
                                 const MoistureComponentIndices& moisture_indices) noexcept
{
    const amrex::Real rho = cons_arr(i, j, k, Rho_comp);
    const bool has_moisture = (moisture_indices.qv >= 0);
    // Microphysical fields are sampled as mixing ratios. The conserved state
    // stores rho*q, so field evaluation divides by rho before interpolation.
    const amrex::Real qv_for_eos = has_moisture
        ? cons_arr(i, j, k, moisture_indices.qv) / rho
        : amrex::Real(0.0);

    switch (field_id) {
    case SampledFieldID::Rho:
        return rho;
    case SampledFieldID::Theta:
        return cons_arr(i, j, k, RhoTheta_comp) / rho;
    case SampledFieldID::Temp:
        return getTgivenRandRTh(rho, cons_arr(i, j, k, RhoTheta_comp), qv_for_eos);
    case SampledFieldID::Pressure:
        // In an anelastic run the compressible EOS does not give the pressure of the
        //    system -- rho is the (frozen) base state density, so getPgivenRTh would
        //    respond to theta perturbations.  We use the base state pressure instead,
        //    just as the 3-D plotfile path does.
        return use_hse_pressure ? p_hse_arr(i, j, k)
                                : getPgivenRTh(cons_arr(i, j, k, RhoTheta_comp), qv_for_eos);
    case SampledFieldID::HeightMSL:
        if (have_z_phys_cc) {
            return z_phys_cc_arr(i, j, k);
        }
        return Compute_Z_AtCellCenter(i, j, k, z_phys_nd_arr);
    case SampledFieldID::HeightAGL:
        return Compute_Zrel_AtCellCenter(i, j, k, z_phys_nd_arr);
    case SampledFieldID::Qv:
        return has_moisture ? cons_arr(i, j, k, moisture_indices.qv) / rho : amrex::Real(0.0);
    case SampledFieldID::Qc:
        return (moisture_indices.qc >= 0) ? cons_arr(i, j, k, moisture_indices.qc) / rho : amrex::Real(0.0);
    case SampledFieldID::Qi:
        return (moisture_indices.qi >= 0) ? cons_arr(i, j, k, moisture_indices.qi) / rho : amrex::Real(0.0);
    case SampledFieldID::Qr:
        return (moisture_indices.qr >= 0) ? cons_arr(i, j, k, moisture_indices.qr) / rho : amrex::Real(0.0);
    case SampledFieldID::Qs:
        return (moisture_indices.qs >= 0) ? cons_arr(i, j, k, moisture_indices.qs) / rho : amrex::Real(0.0);
    case SampledFieldID::Qg:
        return (moisture_indices.qg >= 0) ? cons_arr(i, j, k, moisture_indices.qg) / rho : amrex::Real(0.0);
    case SampledFieldID::UEast:
    case SampledFieldID::VNorth:
    case SampledFieldID::W:
    case SampledFieldID::WindSpeed:
    case SampledFieldID::WindDir:
        return amrex::Real(0.0);
    }

    return amrex::Real(0.0);
}

inline const char*
sampled_field_name (SampledFieldID field_id) noexcept
{
    switch (field_id) {
    case SampledFieldID::Rho:       return "rho";
    case SampledFieldID::Theta:     return "theta";
    case SampledFieldID::Temp:      return "temp";
    case SampledFieldID::Pressure:  return "pressure";
    case SampledFieldID::HeightMSL: return "height_msl";
    case SampledFieldID::HeightAGL: return "height_agl";
    case SampledFieldID::Qv:        return "qv";
    case SampledFieldID::Qc:        return "qc";
    case SampledFieldID::Qi:        return "qi";
    case SampledFieldID::Qr:        return "qr";
    case SampledFieldID::Qs:        return "qs";
    case SampledFieldID::Qg:        return "qg";
    case SampledFieldID::UEast:     return "u_east";
    case SampledFieldID::VNorth:    return "v_north";
    case SampledFieldID::W:         return "w";
    case SampledFieldID::WindSpeed: return "wind_speed";
    case SampledFieldID::WindDir:    return "wind_dir";
    }

    return "unknown";
}

} // namespace plotfile2d

#endif
