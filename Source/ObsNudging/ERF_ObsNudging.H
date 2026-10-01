/**
 * \file ERF_ObsNudging.H
 */
#ifndef ERF_OBSNUDGING_H_
#define ERF_OBSNUDGING_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <memory>
#include <string>

#include "ERF_DataStruct.H"
#include "ERF_LatLonMap.H"
#include "ERF_ObsNudgingKernel.H"
#include "ERF_ObsNudgingSeries.H"

class ProblemBase;

//
// Nudging towards observations at stations (met masts, lidars, profilers).
//
// Each station is a horizontal position and a time series of measured
// profiles (see ERF_ObsNudgingSeries.H for the file).  Near a station the
// velocity components and the potential temperature relax towards the
// measurements, with a weight that falls off with the distance from the
// station horizontally and from its measured height range vertically, and
// only the part of the local value outside the measured variability, mean +-
// alpha sigma, is relaxed (see ERF_ObsNudgingKernel.H for the formula).
//
// Heights are measured above the local terrain (or, per station, above the
// model's z = 0): on a terrain-fitted mesh the terrain is the bottom of the
// mesh, with immersed-forcing terrain it is the terrain surface the immersed
// boundary is built from, and cells inside the immersed terrain or buildings
// are not nudged.  A station can be placed by x/y in domain coordinates or, in
// a run with latitude/longitude arrays, by lat/long; an earth-relative wind is
// rotated into the grid frame with the orientation of the grid at the station.
//
// Inputs are under erf.obs_nudging; see the Inputs section of the docs.
//
struct ObsNudgingStation
{
    std::string name;
    std::string file;

    bool use_latlon = false;
    amrex::Real req_lat = amrex::Real(0.0);
    amrex::Real req_lon = amrex::Real(0.0);

    // Domain coordinates, given or resolved from lat/long
    amrex::Real x = amrex::Real(0.0);
    amrex::Real y = amrex::Real(0.0);

    // The lat/lon actually reached, for a station placed by lat/long
    amrex::Real got_lat = amrex::Real(0.0);
    amrex::Real got_lon = amrex::Real(0.0);

    bool earth_frame = true;   // the wind in the file is east/north
    bool agl = true;           // heights are above the local terrain

    // Rotation from the frame of the file to the grid frame
    amrex::Real cos_alpha = amrex::Real(1.0);
    amrex::Real sin_alpha = amrex::Real(0.0);

    obs_nudging::StationSeries series;
};

class ObsNudging
{
public:
    // Reads erf.obs_nudging.* and every station file.  prob supplies the
    // terrain surface of immersed-forcing terrain; start_time is the epoch time
    // of t = 0, used when use_datetime is true and the files give epoch times.
    ObsNudging (TerrainType terrain_type, ProblemBase* prob,
                bool use_datetime, double start_time);

    // Does any station need the latitude/longitude arrays (placed by lat/long,
    // or an earth-relative wind to rotate)?
    [[nodiscard]] bool wants_latlon () const;

    // Place the stations and set the wind rotation.  latlon may be null only
    // when the run has no latitude/longitude arrays, in which case a station
    // placed by lat/long aborts and the grid is taken to be aligned with east
    // and north.  Collective.
    void resolve_positions (const amrex::Geometry& geom0, const LatLonMap* latlon);

    void print_summary () const;

    // Find the terrain under level lev now rather than at its first step, so
    // that a level the nudging cannot handle is refused at start-up
    void prepare_level (int lev, const amrex::Geometry& geom, const amrex::MultiFab& cons,
                        const amrex::MultiFab* z_phys_nd, double time);

    // Add rho * d(theta)/dt to the RhoTheta component of cc_src at level lev.
    // z_phys_nd may be null on a flat mesh; blank is the cell-centred terrain
    // blanking of immersed forcing, or null.
    void add_theta_source (int lev, double time, amrex::Real dt, const amrex::Geometry& geom,
                           const amrex::MultiFab& cons, amrex::MultiFab& cc_src,
                           const amrex::MultiFab* z_phys_nd, const amrex::MultiFab* blank);

    // Add rho * d(u,v,w)/dt to the momentum sources at level lev.
    void add_momentum_sources (int lev, double time, amrex::Real dt, const amrex::Geometry& geom,
                               const amrex::MultiFab& cons,
                               const amrex::MultiFab& xvel, const amrex::MultiFab& yvel,
                               const amrex::MultiFab& zvel,
                               amrex::MultiFab& xmom_src, amrex::MultiFab& ymom_src,
                               amrex::MultiFab& zmom_src,
                               const amrex::MultiFab* z_phys_nd,
                               const amrex::MultiFab* blank_x, const amrex::MultiFab* blank_y,
                               const amrex::MultiFab* blank_z);

private:
    // Interpolate every station in time and copy the profiles to the device
    void update_targets (double time);

    // The terrain height under level lev, on the nodes of a k = 0 slab of its
    // grids; rebuilt when the grids change or the terrain moves
    const amrex::MultiFab& surface (int lev, const amrex::Geometry& geom,
                                    const amrex::MultiFab& cons,
                                    const amrex::MultiFab* z_phys_nd, double time);

    [[nodiscard]] obs_nudging::ObsNudgingView view (amrex::Real dt, const amrex::Geometry& geom) const;

    amrex::Vector<ObsNudgingStation> m_stations;

    TerrainType  m_terrain_type = TerrainType::None;
    ProblemBase* m_prob = nullptr;

    // Parameters
    amrex::Real m_tau = amrex::Real(-1.0);
    amrex::Real m_horizontal_radius = amrex::Real(500.0);
    amrex::Real m_vertical_radius = amrex::Real(25.0);
    amrex::Real m_cutoff = amrex::Real(6.0);
    amrex::Real m_sigma_factor = amrex::Real(1.0);
    amrex::Real m_missing_value = amrex::Real(-9999.0);
    bool m_nudge_wind = true;
    bool m_nudge_w = true;
    bool m_nudge_theta = true;
    double m_time_origin = 0.0;   // added to the run time to get the time of the files

    // The profiles at the time last asked for, on the device
    bool   m_have_targets = false;
    double m_targets_time = 0.0;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_x, m_d_y, m_d_z, m_d_mean, m_d_sigma;
    amrex::Gpu::DeviceVector<int> m_d_agl, m_d_offset, m_d_count;

    // Terrain height per level, with the grids it was built for.  On a moving
    // fitted mesh the terrain is also a function of the stage time, so that is
    // part of the key too; on every other mesh the terrain does not move and
    // the time is ignored.
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_zsurf;
    amrex::Vector<amrex::BoxArray> m_zsurf_ba;
    amrex::Vector<amrex::DistributionMapping> m_zsurf_dm;
    amrex::Vector<double> m_zsurf_time;
};

#endif
