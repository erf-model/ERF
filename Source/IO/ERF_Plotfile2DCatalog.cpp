/**
 * \file ERF_Plotfile2DCatalog.cpp
 */
#include "ERF_Plotfile2DCatalog.H"

#include <deque>
#include <string>

namespace plotfile2d
{

// This catalog is the source of truth for built-in 2D plotfile names and
// metadata. Its order defines the canonical component order used for user
// selection. Keep this order synchronized with the fill blocks in
// ERF_Plotfile2D.cpp until those fill blocks move to dedicated diagnostic
// modules.
// The catalog is metadata only. It must not become a home for diagnostic
// science.

namespace
{

const amrex::Vector<DiagnosticDescriptor>& catalog_storage ()
{
    static const amrex::Vector<DiagnosticDescriptor> catalog{
        {DiagnosticID::ZSurf,         "z_surf",       "Surface elevation",                               "m",         DiagnosticCategory::Geometry,       MissingPolicy::AlwaysAvailable},
        {DiagnosticID::LandMask,       "landmask",     "Land-sea mask",                                   "1",         DiagnosticCategory::Geometry,       MissingPolicy::AlwaysAvailable},
        {DiagnosticID::MapFac,         "mapfac",       "Map factor at mass points",                       "1",         DiagnosticCategory::Geometry,       MissingPolicy::AlwaysAvailable},
        {DiagnosticID::LatM,           "lat_m",        "Latitude at unstaggered mass points",             "deg",       DiagnosticCategory::Geometry,       MissingPolicy::FillZeroWhenUnavailable},
        {DiagnosticID::LonM,           "lon_m",        "Longitude at unstaggered mass points",            "deg",       DiagnosticCategory::Geometry,       MissingPolicy::FillZeroWhenUnavailable},
        {DiagnosticID::UStar,          "u_star",       "Friction velocity from the surface layer",        "m/s",       DiagnosticCategory::SurfaceLayer,    MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::WStar,          "w_star",       "Convective velocity scale from the surface layer","m/s",       DiagnosticCategory::SurfaceLayer,    MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::TStar,          "t_star",       "Temperature scale from the surface layer",        "K",         DiagnosticCategory::SurfaceLayer,    MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::QStar,          "q_star",       "Humidity scale from the surface layer",           "kg/kg",     DiagnosticCategory::SurfaceLayer,    MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::Olen,           "Olen",         "Obukhov length from the surface layer",          "m",         DiagnosticCategory::SurfaceLayer,    MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::Pblh,           "pblh",         "Planetary boundary layer height from the active PBL diagnostic provider", "m",         DiagnosticCategory::SurfaceLayer,    MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::TSurf,          "t_surf",       "Surface temperature from the surface layer",      "K",         DiagnosticCategory::SurfaceLayer,    MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::QSurf,          "q_surf",       "Surface humidity from the surface layer",         "kg/kg",     DiagnosticCategory::SurfaceLayer,    MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::Z0,             "z0",           "Roughness height from the surface layer",         "m",         DiagnosticCategory::SurfaceLayer,    MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::OLR,            "OLR",          "Outgoing longwave radiation at the model top",    "W/m^2",     DiagnosticCategory::Radiation,       MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::SensFlux,       "sens_flux",    "Surface sensible heat flux",                      "kg K m^-2 s^-1", DiagnosticCategory::SurfaceFlux, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LatenFlux,      "laten_flux",   "Surface moisture flux (legacy output name)",      "kg m^-2 s^-1", DiagnosticCategory::SurfaceFlux,    MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::SurfPres,       "surf_pres",    "Surface pressure",                                "Pa",        DiagnosticCategory::SurfaceState,   MissingPolicy::AlwaysAvailable},
        {DiagnosticID::SeaLevelPressure, "sea_level_pressure", "Sea-level pressure from the NMC/NGM reduction with Shuell correction", "Pa", DiagnosticCategory::SurfaceState, MissingPolicy::AlwaysAvailable},
        {DiagnosticID::PrecipTotalAccum,   "precip_total_accum",   "Accumulated surface precipitation, liquid-water equivalent", "kg/m^2", DiagnosticCategory::Precipitation, MissingPolicy::FillZeroWhenUnavailable},
        {DiagnosticID::PrecipRainAccum,    "precip_rain_accum",    "Accumulated surface rain precipitation, liquid-water equivalent", "kg/m^2", DiagnosticCategory::Precipitation, MissingPolicy::FillZeroWhenUnavailable},
        {DiagnosticID::PrecipSnowAccum,    "precip_snow_accum",    "Accumulated surface snow precipitation, liquid-water equivalent", "kg/m^2", DiagnosticCategory::Precipitation, MissingPolicy::FillZeroWhenUnavailable},
        {DiagnosticID::PrecipGraupelAccum, "precip_graupel_accum", "Accumulated surface graupel precipitation, liquid-water equivalent", "kg/m^2", DiagnosticCategory::Precipitation, MissingPolicy::FillZeroWhenUnavailable},
        {DiagnosticID::PrecipHailAccum,    "precip_hail_accum",    "Accumulated surface hail precipitation, liquid-water equivalent", "kg/m^2", DiagnosticCategory::Precipitation, MissingPolicy::FillZeroWhenUnavailable},
        {DiagnosticID::PrecipFrozenAccum,  "precip_frozen_accum",  "Accumulated frozen surface precipitation, liquid-water equivalent", "kg/m^2", DiagnosticCategory::Precipitation, MissingPolicy::FillZeroWhenUnavailable},
        {DiagnosticID::IntegratedQv,   "integrated_qv","Column-integrated water vapor",                   "kg/m^2",    DiagnosticCategory::ColumnIntegral, MissingPolicy::FillZeroWhenUnavailable},
        {DiagnosticID::IntegratedQc,   "integrated_qc","Column-integrated cloud liquid water",            "kg/m^2",    DiagnosticCategory::ColumnIntegral, MissingPolicy::FillZeroWhenUnavailable},
        {DiagnosticID::IntegratedQi,   "integrated_qi","Column-integrated cloud ice",                     "kg/m^2",    DiagnosticCategory::ColumnIntegral, MissingPolicy::FillZeroWhenUnavailable},
        {DiagnosticID::IntegratedQr,   "integrated_qr","Column-integrated rain water",                    "kg/m^2",    DiagnosticCategory::ColumnIntegral, MissingPolicy::FillZeroWhenUnavailable},
        {DiagnosticID::IntegratedQs,   "integrated_qs","Column-integrated snow",                          "kg/m^2",    DiagnosticCategory::ColumnIntegral, MissingPolicy::FillZeroWhenUnavailable},
        {DiagnosticID::IntegratedQg,   "integrated_qg","Column-integrated graupel",                       "kg/m^2",    DiagnosticCategory::ColumnIntegral, MissingPolicy::FillZeroWhenUnavailable},
        {DiagnosticID::SurfaceDiagnosticSource,
                                         "surface_diagnostic_source",
                                                        "Surface diagnostic source code",
                                                        "1",         DiagnosticCategory::SurfaceLayer,    MissingPolicy::AlwaysAvailable},
        {DiagnosticID::SensibleHeatFlux,"sensible_heat_flux","Surface sensible heat flux",                  "W m^-2",    DiagnosticCategory::SurfaceFlux,     MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LatentHeatFlux,  "latent_heat_flux","Surface latent heat flux",                    "W m^-2",    DiagnosticCategory::SurfaceFlux,     MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::ShocUStar,      "shoc_u_star",  "Native SHOC friction velocity diagnostic",        "m/s",       DiagnosticCategory::PBL,             MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::ShocOlen,       "shoc_Olen",    "Native SHOC Obukhov length diagnostic",            "m",         DiagnosticCategory::PBL,             MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::ShocWthvSfc,    "shoc_wthv_sfc","Native SHOC surface virtual potential temperature flux", "K m s^-1", DiagnosticCategory::PBL,             MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LandSurfaceTsfC, "t_sfc", "Noah-MP radiative surface temperature", "K", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LandSurfaceEmissivity, "sfc_emis", "Surface bulk emissivity", "1", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LandSurfaceAlbDirVis, "sfc_alb_dir_vis", "Direct visible surface albedo", "1", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LandSurfaceAlbDirNir, "sfc_alb_dir_nir", "Direct near-infrared surface albedo", "1", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LandSurfaceAlbDifVis, "sfc_alb_dif_vis", "Diffuse visible surface albedo", "1", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LandSurfaceAlbDifNir, "sfc_alb_dif_nir", "Diffuse near-infrared surface albedo", "1", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LandSurfaceCosZenith, "cos_zenith_angle", "Cosine of the solar zenith angle supplied to Noah-MP", "1", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LandSurfaceSwFluxDn, "sw_flux_dn", "Downwelling shortwave flux supplied to Noah-MP", "W m^-2", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LandSurfaceSwFluxDnDirVis, "sw_flux_dn_dir_vis", "Direct visible downwelling shortwave flux", "W m^-2", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LandSurfaceSwFluxDnDirNir, "sw_flux_dn_dir_nir", "Direct near-infrared downwelling shortwave flux", "W m^-2", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LandSurfaceSwFluxDnDifVis, "sw_flux_dn_dif_vis", "Diffuse visible downwelling shortwave flux", "W m^-2", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LandSurfaceSwFluxDnDifNir, "sw_flux_dn_dif_nir", "Diffuse near-infrared downwelling shortwave flux", "W m^-2", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LandSurfaceLwFluxDn, "lw_flux_dn", "Downwelling longwave flux supplied to Noah-MP", "W m^-2", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LandSurfaceGrdflx, "grdflx", "Ground or snow heat flux", "W m^-2", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LandSurfaceFira, "fira", "Total net longwave flux, positive toward the atmosphere", "W m^-2", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LandSurfaceSav, "sav", "Solar radiation absorbed by vegetation", "W m^-2", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LandSurfaceSag, "sag", "Solar radiation absorbed by the ground", "W m^-2", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LandSurfaceAlbedo, "albedo", "Broadband surface albedo", "1", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LandSurfaceSfcrunoff, "sfcrunoff", "Accumulated surface runoff", "m", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LandSurfaceUdrunoff, "udrunoff", "Accumulated subsurface runoff", "m", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::NoahmpTemperature2mVegetated, "noahmp_temperature_2m_vegetated", "Noah-MP 2-m temperature over the vegetated fraction", "K", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::NoahmpTemperature2mBare, "noahmp_temperature_2m_bare", "Noah-MP 2-m temperature over the bare fraction", "K", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::NoahmpWaterVaporMixingRatio2mVegetated, "noahmp_water_vapor_mixing_ratio_2m_vegetated", "Noah-MP 2-m water-vapor mixing ratio over the vegetated fraction", "kg kg^-1 dry air", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::NoahmpWaterVaporMixingRatio2mBare, "noahmp_water_vapor_mixing_ratio_2m_bare", "Noah-MP 2-m water-vapor mixing ratio over the bare fraction", "kg kg^-1 dry air", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::NoahmpVegetationFraction, "noahmp_vegetation_fraction", "Noah-MP vegetation fraction", "1", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::Temperature2m, "temperature_2m", "Physical air temperature 2 m above the local surface", "K", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::WaterVaporMixingRatio2m, "water_vapor_mixing_ratio_2m", "Water-vapor mixing ratio 2 m above the local surface per unit dry-air mass", "kg kg^-1 dry air", DiagnosticCategory::LandSurface, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::NearSurfaceDiagnosticSource, "near_surface_diagnostic_source", "Source code for the unified near-surface diagnostic bundle", "1", DiagnosticCategory::LandSurface, MissingPolicy::AlwaysAvailable},
    };

    return catalog;
}

struct DynamicDescriptorStorage
{
    std::string name;
    std::string long_name;
    std::string units;
    DiagnosticDescriptor descriptor;
};

std::deque<DynamicDescriptorStorage>& dynamic_storage ()
{
    static std::deque<DynamicDescriptorStorage> storage;
    return storage;
}

const DiagnosticDescriptor* ensure_dynamic_soil_descriptor (const std::string& name)
{
    for (const auto& item : dynamic_storage()) {
        if (item.name == name) {
            return &item.descriptor;
        }
    }

    const auto underscore = name.find('_');
    if (underscore == std::string::npos) {
        return nullptr;
    }

    const std::string group = name.substr(0, underscore);
    const int layer = std::stoi(name.substr(underscore + 1));
    const char* units = (group == "tslb") ? "K" : "m^3 m^-3";
    const char* meaning = (group == "smois") ? "Volumetric total soil moisture"
                         : (group == "sh2o") ? "Volumetric liquid soil water"
                                             : "Soil temperature";
    dynamic_storage().emplace_back();
    auto& item = dynamic_storage().back();
    item.name = name;
    item.long_name = std::string(meaning) + " in layer " + std::to_string(layer);
    item.units = units;
    const int group_id = (group == "smois") ? 0 : (group == "sh2o") ? 1 : 2;
    const int dynamic_id = static_cast<int>(DiagnosticID::DynamicSoilBase) +
                           group_id * 10000 + layer;
    item.descriptor = {static_cast<DiagnosticID>(dynamic_id), item.name.c_str(),
                       item.long_name.c_str(), item.units.c_str(),
                       DiagnosticCategory::LandSurface,
                       MissingPolicy::FillMinus999WhenUnavailable};
    return &item.descriptor;
}

} // namespace

const amrex::Vector<DiagnosticDescriptor>&
diagnostic_catalog ()
{
    return catalog_storage();
}

amrex::Vector<std::string>
diagnostic_names ()
{
    amrex::Vector<std::string> names;
    names.reserve(diagnostic_catalog().size());

    for (const auto& descriptor : diagnostic_catalog()) {
        names.push_back(descriptor.name);
    }

    return names;
}

const DiagnosticDescriptor*
find_diagnostic (const std::string& name)
{
    for (const auto& descriptor : diagnostic_catalog()) {
        if (name == descriptor.name) {
            return &descriptor;
        }
    }

    return nullptr;
}

amrex::Vector<std::string>
dynamic_soil_diagnostic_names (int nsoil)
{
    amrex::Vector<std::string> names;
    if (nsoil < 1) {
        return names;
    }
    names.reserve(3*nsoil);
    for (const char* group : {"smois", "sh2o", "tslb"}) {
        for (int layer = 1; layer <= nsoil; ++layer) {
            const std::string name = std::string(group) + "_" + std::to_string(layer);
            ensure_dynamic_soil_descriptor(name);
            names.push_back(name);
        }
    }
    return names;
}

amrex::Vector<std::string>
dynamic_soil_diagnostic_names (const amrex::Vector<std::string>& active_lsm_names)
{
    amrex::Vector<std::string> names;
    for (const auto& name : active_lsm_names) {
        if (is_dynamic_soil_diagnostic_name(name)) {
            ensure_dynamic_soil_descriptor(name);
            names.push_back(name);
        }
    }
    return names;
}

const DiagnosticDescriptor*
find_dynamic_soil_diagnostic (const std::string& name)
{
    if (!is_dynamic_soil_diagnostic_name(name)) {
        return nullptr;
    }
    return ensure_dynamic_soil_descriptor(name);
}

bool
is_dynamic_soil_diagnostic_name (const std::string& name)
{
    const auto underscore = name.find('_');
    if (underscore == std::string::npos) {
        return false;
    }
    const std::string group = name.substr(0, underscore);
    if (group != "smois" && group != "sh2o" && group != "tslb") {
        return false;
    }
    try {
        std::size_t parsed = 0;
        const int layer = std::stoi(name.substr(underscore + 1), &parsed);
        return parsed == name.size() - underscore - 1 && layer > 0;
    } catch (...) {
        return false;
    }
}

} // namespace plotfile2d
