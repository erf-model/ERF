#ifndef ERF_WALL_THERMODYNAMICS_H_
#define ERF_WALL_THERMODYNAMICS_H_

#include <AMReX_Array.H>
#include <AMReX_REAL.H>

namespace erf_wall_thermodynamics {

// These compact enums and the face record are copied into GPU kernels.  Keep
// parsing strings and all other host-only state outside this representation.
enum class ThermalMode : unsigned char {
    LegacyNumerical,
    FixedPhysicalTemperature
};

enum class MoistureMode : unsigned char {
    LegacyNumerical,
    DryImpermeable,
    WetEquilibrium
};

enum class TransferModel : unsigned char {
    None,
    ResolvedMolecular
};

struct Face {
    ThermalMode thermal_mode = ThermalMode::LegacyNumerical;
    MoistureMode moisture_mode = MoistureMode::LegacyNumerical;
    TransferModel transfer_model = TransferModel::None;
    amrex::Real physical_temperature_K = amrex::Real(0.0);
};

using Boundary = amrex::GpuArray<Face, 2 * AMREX_SPACEDIM>;

} // namespace erf_wall_thermodynamics

#endif
