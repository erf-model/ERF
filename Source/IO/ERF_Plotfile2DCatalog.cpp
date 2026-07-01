#include "ERF_Plotfile2DCatalog.H"

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
        {DiagnosticID::Pblh,           "pblh",         "Planetary boundary layer height",                "m",         DiagnosticCategory::SurfaceLayer,    MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::TSurf,          "t_surf",       "Surface temperature from the surface layer",      "K",         DiagnosticCategory::SurfaceLayer,    MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::QSurf,          "q_surf",       "Surface humidity from the surface layer",         "kg/kg",     DiagnosticCategory::SurfaceLayer,    MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::Z0,             "z0",           "Roughness height from the surface layer",         "m",         DiagnosticCategory::SurfaceLayer,    MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::OLR,            "OLR",          "Outgoing longwave radiation at the model top",    "W/m^2",     DiagnosticCategory::Radiation,       MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::SensFlux,       "sens_flux",    "Surface sensible heat flux",                      "kg K m^-2 s^-1", DiagnosticCategory::SurfaceFlux, MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LatenFlux,      "laten_flux",   "Surface moisture flux (legacy output name)",      "kg m^-2 s^-1", DiagnosticCategory::SurfaceFlux,    MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::SurfPres,       "surf_pres",    "Surface pressure",                                "Pa",        DiagnosticCategory::SurfaceState,   MissingPolicy::AlwaysAvailable},
        {DiagnosticID::IntegratedQv,   "integrated_qv","Column-integrated water vapor",                   "kg/m^2",    DiagnosticCategory::ColumnIntegral, MissingPolicy::FillZeroWhenUnavailable},
        {DiagnosticID::SurfaceDiagnosticSource,
                                         "surface_diagnostic_source",
                                                        "Surface diagnostic source code",
                                                        "1",         DiagnosticCategory::SurfaceLayer,    MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::SensibleHeatFlux,"sensible_heat_flux","Surface sensible heat flux",                  "W m^-2",    DiagnosticCategory::SurfaceFlux,     MissingPolicy::FillMinus999WhenUnavailable},
        {DiagnosticID::LatentHeatFlux,  "latent_heat_flux","Surface latent heat flux",                    "W m^-2",    DiagnosticCategory::SurfaceFlux,     MissingPolicy::FillMinus999WhenUnavailable},
    };

    return catalog;
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

} // namespace plotfile2d
