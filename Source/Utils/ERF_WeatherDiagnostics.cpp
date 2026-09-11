#ifndef ERF_WEATHER_DIAGNOSTICS_H_
#define ERF_WEATHER_DIAGNOSTICS_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelReduce.H>
#include <limits>

#include "ERF_DataStruct.H"
#include "ERF.H"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace amrex;

namespace fs = std::filesystem;

/**
 * Routines to compute hurricane diagnostics
 */

#ifndef M_PI
#define M_PI Real(3.14159265358979323846)
#endif

namespace {

/**
 * Linearize a 2D (i,j) index relative to the domain so that an arg-min can be
 * carried out with a single atomic. Using one packed index rather than storing
 * i and j separately keeps the recorded location consistent with the recorded
 * minimum, and makes the choice among tied cells deterministic.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Long pack_ij (const int i, const int j, const int nx, const Dim3& dlo) noexcept
{
    return static_cast<Long>(j - dlo.y) * static_cast<Long>(nx)
         + static_cast<Long>(i - dlo.x);
}

/**
 * Invert pack_ij. A packed index that was never set (still the initial
 * sentinel) decodes to the (-1,-1) "no local candidate" marker.
 */
void unpack_ij (const Long idx, const int nx, const Dim3& dlo, int& i, int& j) noexcept
{
    if (idx == std::numeric_limits<Long>::max()) {
        i = -1;
        j = -1;
    } else {
        i = static_cast<int>(idx % static_cast<Long>(nx)) + dlo.x;
        j = static_cast<int>(idx / static_cast<Long>(nx)) + dlo.y;
    }
}

} // anonymous namespace


/**
 * Compute the global minimum and its location across all ranks.
 *
 * @param[in] sc Solver choices
 * @param[in] lev_geom Geometry of the current level
 * @param[in] S_data Conservative state data
 * @param[in] d_val_min_ptr Device pointer to local minimum value
 * @param[in] d_i_min_ptr Device pointer to local minimum i-index
 * @param[in] d_j_min_ptr Device pointer to local minimum j-index
 * @param[out] global_val_min Global minimum value
 * @param[out] station_loc_i Global minimum i-index
 * @param[out] station_loc_j Global minimum j-index
 */
std::pair<int, int> 
ComputeLocation (const SolverChoice& sc,
                 const Geometry& lev_geom,
                 const Vector<MultiFab>& S_data,
                 Real* d_val_min_ptr,
                 int* d_i_min_ptr,
                 int* d_j_min_ptr)
{
    Real h_val_min;
    int h_i_min, h_j_min;

    Gpu::copy(Gpu::deviceToHost, d_val_min_ptr, d_val_min_ptr + 1, &h_val_min);
    Gpu::copy(Gpu::deviceToHost, d_i_min_ptr, d_i_min_ptr + 1, &h_i_min);
    Gpu::copy(Gpu::deviceToHost, d_j_min_ptr, d_j_min_ptr + 1, &h_j_min);
    Gpu::synchronize();

    Real local_val_min = h_val_min;
    int local_i_min = h_i_min;
    int local_j_min = h_j_min;

    int rank = ParallelDescriptor::MyProc();

    // NOTE: reduce through the amrex wrappers rather than a hard-coded
    //       MPI_DOUBLE_INT MINLOC. The latter is a type mismatch whenever
    //       amrex::Real is float (ERF_PRECISION=SINGLE), in which case MPI
    //       reads and writes 8 bytes of a 4-byte member. Reducing the value
    //       and then taking the smallest rank that attains it reproduces
    //       MPI_MINLOC's lowest-rank tie-break without naming an MPI type.
    Real global_val_min = local_val_min;
    ParallelDescriptor::ReduceRealMin(global_val_min);

    int owner_rank = (local_val_min == global_val_min) ? rank : ParallelDescriptor::NProcs();
    ParallelDescriptor::ReduceIntMin(owner_rank);
    AMREX_ALWAYS_ASSERT(owner_rank < ParallelDescriptor::NProcs());

    // Broadcast the indices from the rank that owns the minimum
    int station_loc_i = local_i_min;
    int station_loc_j = local_j_min;

    ParallelDescriptor::Bcast(&station_loc_i, 1, owner_rank);
    ParallelDescriptor::Bcast(&station_loc_j, 1, owner_rank);

    if (rank == 0) {
        Print() << "Global minimum distance to station location (k=0): "
                       << global_val_min << " at (i,j) = ("
                       << station_loc_i << ", " << station_loc_j << ")\n";
    }
     return {station_loc_i, station_loc_j};
}

/**
 * Initialize the hurricane eye tracker using a given latitude and longitude.
 *
 * @param[in] sc Solver choices
 * @param[in] lev_geom Geometry of the current level
 * @param[in] S_data Conservative state data
 * @param[in] station_latitude Target latitude for the eye
 * @param[in] station_longitude Target longitude for the eye
 */
std::pair<int, int>
ERF::ComputeStationLocationIJ (const SolverChoice& sc,
                               const Geometry& lev_geom,
                               const Vector<MultiFab>& S_data,
                               const Real& station_latitude,
                               const Real& station_longitude)
{
    int levc = finest_level;
    Gpu::DeviceScalar<Real> d_val_min(1e10);
    Gpu::DeviceScalar<int> d_i_min(-1), d_j_min(-1);

    Real* d_val_min_ptr = d_val_min.dataPtr();
    int* d_i_min_ptr = d_i_min.dataPtr();
    int* d_j_min_ptr = d_j_min.dataPtr();

    // NOTE: the arg-min is done in two passes. A single pass that takes an
    //       atomic min of the distance and then plainly stores i and j is
    //       racy: a thread with a worse distance can still win the store and
    //       leave an eye location that does not belong to the recorded
    //       minimum. The first pass reduces the distance only; the second
    //       records the location of the cells that attain it, through one
    //       atomic on a packed index so value and location stay consistent.
    const Dim3 dlo = lbound(lev_geom.Domain());
    const int nx = lev_geom.Domain().length(0);

    Gpu::DeviceScalar<Long> d_idx_min(std::numeric_limits<Long>::max());
    Long* d_idx_min_ptr = d_idx_min.dataPtr();

    if(sc.init_type == InitType::WRFInput){
        for (MFIter mfi(S_data[IntVars::cons]); mfi.isValid(); ++mfi) {
            const Box& box = mfi.validbox();
            FArrayBox& fab_lat = (*(lat_m[levc]))[mfi];
            FArrayBox& fab_lon = (*(lon_m[levc]))[mfi];
            const Array4<Real>& lat_arr = fab_lat.array();
            const Array4<Real>& lon_arr = fab_lon.array();

            ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (k==0) {
                    Real dlat = lat_arr(i,j,0) - station_latitude;
                    Real dlon = lon_arr(i,j,0) - station_longitude;
                    Real dist = std::sqrt(dlat*dlat + dlon*dlon);
                    Gpu::Atomic::Min(&d_val_min_ptr[0], dist);
                }
            });
        }

        // The minimum over every box on this rank must be known before the
        // locating pass below can test against it.
        Gpu::synchronize();

        for (MFIter mfi(S_data[IntVars::cons]); mfi.isValid(); ++mfi) {
            const Box& box = mfi.validbox();
            FArrayBox& fab_lat = (*(lat_m[levc]))[mfi];
            FArrayBox& fab_lon = (*(lon_m[levc]))[mfi];
            const Array4<Real>& lat_arr = fab_lat.array();
            const Array4<Real>& lon_arr = fab_lon.array();

            ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (k==0) {

                    Real dlat = lat_arr(i,j,0) - station_latitude;
                    Real dlon = lon_arr(i,j,0) - station_longitude;
                    Real dist = std::sqrt(dlat*dlat + dlon*dlon);
                    if (dist == d_val_min_ptr[0]) {
                        Gpu::Atomic::Min(d_idx_min_ptr, pack_ij(i,j,nx,dlo));
                    }
                }
            });
        }
    }

    if(sc.init_type == InitType::HindCast){
        for (amrex::MFIter mfi(S_data[IntVars::cons]); mfi.isValid(); ++mfi) {
            const amrex::Box& box = mfi.validbox();
            const auto& mf_latlon = forecast_state_interp[levc][4];
            const auto latlon_arr = mf_latlon.array(mfi);

            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (k==0) {

                    amrex::Real dlat = latlon_arr(i,j,k,0) - station_latitude;
                    amrex::Real dlon = latlon_arr(i,j,k,1) - station_longitude;
                    amrex::Real dist = std::sqrt(dlat*dlat + dlon*dlon);
                    amrex::Gpu::Atomic::Min(&d_val_min_ptr[0], dist);
                }
            });
        }

        // The minimum over every box on this rank must be known before the
        // locating pass below can test against it.
        Gpu::synchronize();

        for (amrex::MFIter mfi(S_data[IntVars::cons]); mfi.isValid(); ++mfi) {
            const amrex::Box& box = mfi.validbox();
            const auto& mf_latlon = forecast_state_interp[levc][4];
            const auto latlon_arr = mf_latlon.array(mfi);

            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (k==0) {

                    amrex::Real dlat = latlon_arr(i,j,k,0) - station_latitude;
                    amrex::Real dlon = latlon_arr(i,j,k,1) - station_longitude;
                    amrex::Real dist = std::sqrt(dlat*dlat + dlon*dlon);
                    if (dist == d_val_min_ptr[0]) {
                        amrex::Gpu::Atomic::Min(d_idx_min_ptr, pack_ij(i,j,nx,dlo));
                    }
                }
            });
        }
    }

    Gpu::synchronize();

    // Unpack the located index back into the (i,j) device scalars that
    // ComputeLocation reads.
    {
        int h_i_min, h_j_min;
        unpack_ij(d_idx_min.dataValue(), nx, dlo, h_i_min, h_j_min);
        Gpu::copy(Gpu::hostToDevice, &h_i_min, &h_i_min + 1, d_i_min_ptr);
        Gpu::copy(Gpu::hostToDevice, &h_j_min, &h_j_min + 1, d_j_min_ptr);
    }

    std::pair<int, int> station_loc = ComputeLocation (sc, 
                                                       lev_geom, 
                                                       S_data,
                                                       d_val_min_ptr, 
                                                       d_i_min_ptr, 
                                                       d_j_min_ptr);
    return station_loc;
}

/**
 * Compute and track the maximum wind velocity near the hurricane eye.
 *
 * @param[in] lev_geom Geometry of the current level
 * @param[in] mf_cc_vel MultiFab containing cell-centered velocity
 * @param[in] time Current simulation time
 */
void
ERF::TrackerAtStation_RainAccumulation(const SolverChoice& sc,
                                       const int levc, 
                                       const std::pair<int,int>& station_loc, 
                                       const Real& time)
{
    const MoistureComponentIndices& mi = sc.moisture_indices;
    const int idx = mi.qmoist_index_for_var("rain_accum");

    AMREX_ALWAYS_ASSERT(idx >= 0);
    AMREX_ALWAYS_ASSERT(idx < static_cast<int>(qmoist[levc].size()));
    AMREX_ALWAYS_ASSERT(qmoist[levc][idx] != nullptr);

    MultiFab& rain_accum = *(qmoist[levc][idx]);

    auto [station_loc_i, station_loc_j] = station_loc;

    Gpu::DeviceVector<Real> d_rain_accum(1, -bogus_large_value);
    Real* d_rain_accum_ptr = d_rain_accum.data();

    const IntVect station_iv(station_loc_i, station_loc_j, 0);

    for (MFIter mfi(rain_accum); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();

        if (box.contains(station_iv)) {
            const auto& rain_arr = rain_accum.const_array(mfi);

            ParallelFor(
                Box(station_iv, station_iv),
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    d_rain_accum_ptr[0] = rain_arr(i, j, k, 0);
                });

            break;
        }
    }

    Gpu::synchronize();

    Real h_rain_accum_local = -bogus_large_value;
    Gpu::copy(Gpu::deviceToHost, d_rain_accum.begin(), d_rain_accum.end(), &h_rain_accum_local);

    Real h_rain_accum_global = -bogus_large_value;
    #ifdef AMREX_USE_MPI
        h_rain_accum_global = h_rain_accum_local;
        amrex::ParallelDescriptor::ReduceRealMax(h_rain_accum_global);
    #else
        h_rain_accum_global = h_rain_accum_local;
    #endif

    double time_in_hrs = time / 3600.0;
    station_rain_accum_vs_time.push_back({static_cast<Real>(time_in_hrs), h_rain_accum_global});
}

/**
 * Wrapper to track the hurricane eye position over time.
 *
 * @param[in] sc Solver choices
 */
void
ERF::WeatherDiagnosticsTracker (const SolverChoice& sc)
{
    static bool is_start = true;
    int levc=finest_level;

    const MoistureType moisture_type = sc.moisture_type;
    const Real station_latitude  = sc.station_latitude;
    const Real station_longitude = sc.station_longitude;

    std::pair<int, int> station_loc;

    if(is_start) {
        station_loc = ComputeStationLocationIJ(sc, 
                                               geom[levc],
                                               vars_new[levc],
                                               station_latitude,
                                               station_longitude);
    }

    TrackerAtStation_RainAccumulation(sc,
                                      levc, 
                                      station_loc,
                                      t_new[0]); 
}
#endif
