/**
 * \file ERF_Plotfile2DCatalog.H
 */
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
    SeaLevelPressure,
    PrecipTotalAccum,
    PrecipRainAccum,
    PrecipSnowAccum,
    PrecipGraupelAccum,
    PrecipHailAccum,
    PrecipFrozenAccum,
    IntegratedQv,
    IntegratedQc,
    IntegratedQi,
    IntegratedQr,
    IntegratedQs,
    IntegratedQg,
    SurfaceDiagnosticSource,
    SensibleHeatFlux,
    LatentHeatFlux,
    ShocUStar,
    ShocOlen,
    ShocWthvSfc,
    LandSurfaceTsfC,
    LandSurfaceEmissivity,
    LandSurfaceAlbDirVis,
    LandSurfaceAlbDirNir,
    LandSurfaceAlbDifVis,
    LandSurfaceAlbDifNir,
    LandSurfaceCosZenith,
    LandSurfaceSwFluxDn,
    LandSurfaceSwFluxDnDirVis,
    LandSurfaceSwFluxDnDirNir,
    LandSurfaceSwFluxDnDifVis,
    LandSurfaceSwFluxDnDifNir,
    LandSurfaceLwFluxDn,
    LandSurfaceGrdflx,
    LandSurfaceFira,
    LandSurfaceSav,
    LandSurfaceSag,
    LandSurfaceAlbedo,
    LandSurfaceSfcrunoff,
    LandSurfaceUdrunoff,
    NoahmpTemperature2mVegetated,
    NoahmpTemperature2mBare,
    NoahmpWaterVaporMixingRatio2mVegetated,
    NoahmpWaterVaporMixingRatio2mBare,
    NoahmpVegetationFraction,
    Temperature2m,
    WaterVaporMixingRatio2m,
    NearSurfaceDiagnosticSource,
    DynamicSoilBase
};

enum class DiagnosticCategory
{
    Geometry,
    SurfaceLayer,
    Radiation,
    SurfaceFlux,
    PBL,
    SurfaceState,
    Precipitation,
    ColumnIntegral,
    LandSurface,
    SampledLevel
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

// Runtime Noah-MP soil descriptors use the provider's actual configured layer
// count. Their IDs are unique within the returned runtime inventory.
amrex::Vector<std::string> dynamic_soil_diagnostic_names (int nsoil);
amrex::Vector<std::string>
dynamic_soil_diagnostic_names (const amrex::Vector<std::string>& active_lsm_names);
const DiagnosticDescriptor* find_dynamic_soil_diagnostic (const std::string& name);
bool is_dynamic_soil_diagnostic_name (const std::string& name);

} // namespace plotfile2d

#endif
