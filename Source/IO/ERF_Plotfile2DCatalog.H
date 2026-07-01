#ifndef ERF_PLOTFILE2DCATALOG_H_
#define ERF_PLOTFILE2DCATALOG_H_

#include <string>

#include <AMReX_Vector.H>

namespace plotfile2d
{

enum class DiagnosticID
{
    ZSurf,
    LandMask,
    MapFac,
    LatM,
    LonM,
    UStar,
    WStar,
    TStar,
    QStar,
    Olen,
    Pblh,
    TSurf,
    QSurf,
    Z0,
    OLR,
    SensFlux,
    LatenFlux,
    SurfPres,
    IntegratedQv,
    SurfaceDiagnosticSource,
    SensibleHeatFlux,
    LatentHeatFlux
};

enum class DiagnosticCategory
{
    Geometry,
    SurfaceLayer,
    Radiation,
    SurfaceFlux,
    SurfaceState,
    ColumnIntegral
};

enum class MissingPolicy
{
    AlwaysAvailable,
    FillZeroWhenUnavailable,
    FillMinus999WhenUnavailable
};

struct DiagnosticDescriptor
{
    DiagnosticID id;
    const char* name;
    const char* long_name;
    const char* units;
    DiagnosticCategory category;
    MissingPolicy missing_policy;
};

const amrex::Vector<DiagnosticDescriptor>& diagnostic_catalog ();

amrex::Vector<std::string> diagnostic_names ();

const DiagnosticDescriptor* find_diagnostic (const std::string& name);

} // namespace plotfile2d

#endif
