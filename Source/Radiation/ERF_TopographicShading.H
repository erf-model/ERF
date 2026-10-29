#ifndef ERF_TOPOGRAPHIC_SHADING_H_
#define ERF_TOPOGRAPHIC_SHADING_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <memory>

struct TopographicShadingChoice
{
    bool slope_rad = false;
    bool topo_shading = false;
    amrex::Real shadlen = 25000;
    amrex::Real fixed_mu = -1;
    amrex::Real fixed_azimuth = -1; // degrees clockwise from geographic north
    amrex::Real latitude = 39.809860;
    amrex::Real longitude = -98.555183;
    int radiation_interval = 1;

    bool active () const { return slope_rad || topo_shading; }
    void read (const std::string& prefix);
};

/** Surface-only correction of the radiation supplied to land models.
 * The atmospheric fluxes and heating remain the column solver's outputs.
 * Raw forcing is retained separately, including across radiation intervals,
 * so a refined level never interpolates an already corrected parent flux.
 */
class TopographicShading
{
public:
    // Same order as SurfaceModel's canonical radiation outputs.
    enum RadiationField { CosZen, Shortwave, DirectVis, DirectNir, DiffuseVis, DiffuseNir, Longwave, NRad };
    enum Diagnostic { Shadow, DirectFactor, Incidence, UnshadedSW, TerrainSW, NDiag };
    enum GeometryField { Height, MapX, MapY, Latitude, Longitude, CosAlpha, SinAlpha, NGeom };

    void define_level (int lev, const amrex::BoxArray& ba2d,
                       const amrex::DistributionMapping& dm, const amrex::Geometry& geom);
    bool geometry_ready (int lev) const;
    bool forcing_ready (int lev) const;
    void initialize_geometry (int lev, const TopographicShadingChoice& choice,
                              const amrex::MultiFab& z_nd,
                              const amrex::MultiFab& map_x, const amrex::MultiFab& map_y,
                              const amrex::MultiFab* lat, const amrex::MultiFab* lon);
    amrex::Vector<amrex::MultiFab*> raw_fields (int lev);
    amrex::MultiFab* raw (int lev);
    const amrex::MultiFab* diagnostics (int lev) const;
    void set_sun (int lev, double calday, double declin);
    void interpolate_from_coarse (int lev);
    void apply (int lev, const TopographicShadingChoice& choice,
                const amrex::Vector<amrex::MultiFab*>& destinations);

private:
    struct Level {
        amrex::Geometry geom;
        amrex::MultiFab ground, patch, raw, diagnostic;
        amrex::Vector<std::unique_ptr<amrex::MultiFab>> aliases;
        bool geometry_ready = false;
        bool forcing_ready = false;
        double calday = 0, declin = 0;
    };
    amrex::Vector<std::unique_ptr<Level>> m_levels;
};
#endif
