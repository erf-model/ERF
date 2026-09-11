#ifndef ERF_HURRICANE_DIAGNOSTICS_H_
#define ERF_HURRICANE_DIAGNOSTICS_H_

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
 * @param[out] global_i_min Global minimum i-index
 * @param[out] global_j_min Global minimum j-index
 */
void
ERF::ComputeGlobalMinLocation (const SolverChoice& sc,
                               const Geometry& lev_geom,
                               const Vector<MultiFab>& S_data,
                               Real* d_val_min_ptr,
                               int* d_i_min_ptr,
                               int* d_j_min_ptr,
                               Real& global_val_min,
                               int& global_i_min,
                               int& global_j_min)
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
    global_val_min = local_val_min;
    ParallelDescriptor::ReduceRealMin(global_val_min);

    int owner_rank = (local_val_min == global_val_min) ? rank : ParallelDescriptor::NProcs();
    ParallelDescriptor::ReduceIntMin(owner_rank);
    AMREX_ALWAYS_ASSERT(owner_rank < ParallelDescriptor::NProcs());

    // Broadcast the indices from the rank that owns the minimum
    global_i_min = local_i_min;
    global_j_min = local_j_min;

    ParallelDescriptor::Bcast(&global_i_min, 1, owner_rank);
    ParallelDescriptor::Bcast(&global_j_min, 1, owner_rank);

    if (rank == 0) {
        Print() << "Global minimum distance to hurricane eye (k=0): "
                       << global_val_min << " at (i,j) = ("
                       << global_i_min << ", " << global_j_min << ")\n";
    }

    Gpu::DeviceScalar<Real> d_eye_lat(zero), d_eye_lon(zero);

    Real* d_eye_lat_ptr = d_eye_lat.dataPtr();
    Real* d_eye_lon_ptr = d_eye_lon.dataPtr();

    int levc = finest_level;
    // On owner_rank, compute eye_lat and eye_lon
    if (sc.init_type == InitType::WRFInput and rank == owner_rank) {
        for (MFIter mfi(S_data[IntVars::cons]); mfi.isValid(); ++mfi) {
            const Box& box = mfi.validbox();
            FArrayBox& fab_lat = (*(lat_m[levc]))[mfi];
            FArrayBox& fab_lon = (*(lon_m[levc]))[mfi];
            const Array4<Real>& lat_arr = fab_lat.array();
            const Array4<Real>& lon_arr = fab_lon.array();

            if (box.smallEnd()[2] == 0) {
                Box bx2d = makeSlab(box,2,0);
                ParallelFor(bx2d, [=] AMREX_GPU_DEVICE(int i, int j, int ) {
                    if (i == global_i_min && j == global_j_min) {
                        *d_eye_lat_ptr = lat_arr(i,j,0);
                        *d_eye_lon_ptr = lon_arr(i,j,0);
                    }
                });
            }
        }
    }

    if (sc.init_type == InitType::HindCast and rank == owner_rank) {
        // On owner_rank, compute eye_lat and eye_lon
        if (rank == owner_rank) {
            for (amrex::MFIter mfi(S_data[IntVars::cons]); mfi.isValid(); ++mfi) {
                const amrex::Box& box = mfi.validbox();
                const auto& mf_latlon = forecast_state_interp[levc][4];
                const auto latlon_arr = mf_latlon.array(mfi);
                amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    if (i == global_i_min && j == global_j_min && k == 0) {
                        *d_eye_lat_ptr = latlon_arr(i,j,k,0);
                        *d_eye_lon_ptr = latlon_arr(i,j,k,1);
                    }
                });
            }
        }
    }

    Real eye_lat = d_eye_lat.dataValue();
    Real eye_lon = d_eye_lon.dataValue();

    // Synchronize to ensure the owner has computed values
    Gpu::synchronize();

    ParallelDescriptor::Bcast(&eye_lat, 1, owner_rank);
    ParallelDescriptor::Bcast(&eye_lon, 1, owner_rank);

    const auto dx = lev_geom.CellSizeArray();
    const auto prob_lo = lev_geom.ProbLoArray();

    Real eye_x =  prob_lo[0] + (global_i_min+myhalf)*dx[0];
    Real eye_y =  prob_lo[1] + (global_j_min+myhalf)*dx[1];

    hurricane_eye_track_xy.push_back({eye_x, eye_y});
    hurricane_eye_track_latlon.push_back({eye_lon, eye_lat});
}

/**
 * Generate a circular set of points around the last known eye position.
 */
void
ERF::HurricaneTrackerCircle ()
{
    // Check that there is at least one eye position
    if (hurricane_eye_track_xy.empty()) return;

    // Get the last known (x, y) position of the eye
    const auto [x_last, y_last] = hurricane_eye_track_xy.back();

    // Define circle properties
    const int n_points = 100;        // number of points on the circle
    const Real radius = 200e3; // radius in meters (example: 50 km)

    // Clear previous points and reserve space
    hurricane_tracker_circle.clear();
    hurricane_tracker_circle.reserve(n_points);

    // Fill the circle points
    for (int i = 0; i < n_points; ++i) {
        Real theta = two * static_cast<Real>(M_PI) * static_cast<Real>(i) / static_cast<Real>(n_points);
        Real x = x_last + radius * std::cos(theta);
        Real y = y_last + radius * std::sin(theta);
        hurricane_tracker_circle.push_back({x, y});
    }
}

/**
 * Initialize the hurricane eye tracker using a given latitude and longitude.
 *
 * @param[in] sc Solver choices
 * @param[in] lev_geom Geometry of the current level
 * @param[in] S_data Conservative state data
 * @param[in] hurricane_eye_latitude Target latitude for the eye
 * @param[in] hurricane_eye_longitude Target longitude for the eye
 */
void
ERF::HurricaneEyeTrackerInitial (const SolverChoice& sc,
                                 const Geometry& lev_geom,
                                 const Vector<MultiFab>& S_data,
                                 const Real& hurricane_eye_latitude,
                                 const Real& hurricane_eye_longitude)
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

                    Real dlat = lat_arr(i,j,0) - hurricane_eye_latitude;
                    Real dlon = lon_arr(i,j,0) - hurricane_eye_longitude;
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

                    Real dlat = lat_arr(i,j,0) - hurricane_eye_latitude;
                    Real dlon = lon_arr(i,j,0) - hurricane_eye_longitude;
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

                    amrex::Real dlat = latlon_arr(i,j,k,0) - hurricane_eye_latitude;
                    amrex::Real dlon = latlon_arr(i,j,k,1) - hurricane_eye_longitude;
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

                    amrex::Real dlat = latlon_arr(i,j,k,0) - hurricane_eye_latitude;
                    amrex::Real dlon = latlon_arr(i,j,k,1) - hurricane_eye_longitude;
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
    // ComputeGlobalMinLocation reads.
    {
        int h_i_min, h_j_min;
        unpack_ij(d_idx_min.dataValue(), nx, dlo, h_i_min, h_j_min);
        Gpu::copy(Gpu::hostToDevice, &h_i_min, &h_i_min + 1, d_i_min_ptr);
        Gpu::copy(Gpu::hostToDevice, &h_j_min, &h_j_min + 1, d_j_min_ptr);
    }

    Real global_val_min;
    int global_i_min, global_j_min;

    ComputeGlobalMinLocation(sc, lev_geom, S_data,
                             d_val_min_ptr, d_i_min_ptr, d_j_min_ptr,
                             global_val_min, global_i_min, global_j_min);
}

/**
 * Track the hurricane eye by searching for minimum pressure near the previous position.
 *
 * @param[in] sc Solver choices
 * @param[in] lev_geom Geometry of the current level
 * @param[in] S_data Conservative state data
 * @param[in] moisture_type Moisture model type
 */
void
ERF::HurricaneEyeTrackerNotInitial (const SolverChoice& sc,
                                    const Geometry& lev_geom,
                                    const Vector<MultiFab>& S_data,
                                    MoistureType moisture_type)
{

    if (hurricane_eye_track_xy.empty()) {
        Print() << "Error: hurricane_eye_track_xy is empty!\n";
        Abort("Attempted to access hurricane_eye_track_xy[0]");
    }

    Real tmp_x_eye = hurricane_eye_track_xy.back()[0];
    Real tmp_y_eye = hurricane_eye_track_xy.back()[1];

    if(ParallelDescriptor::IOProcessor()){
        std::cout << "The value of x y are " << tmp_x_eye << " " << tmp_y_eye << std::endl;
    }

    Gpu::DeviceScalar<Real> d_val_min(1e10);
    Gpu::DeviceScalar<int> d_i_min(-1), d_j_min(-1);

    Real* d_val_min_ptr = d_val_min.dataPtr();
    int* d_i_min_ptr = d_i_min.dataPtr();
    int* d_j_min_ptr = d_j_min.dataPtr();

    bool use_moisture = (moisture_type != MoistureType::None);
    const int ncomp = S_data[IntVars::cons].nComp();

    const auto dx = lev_geom.CellSizeArray();
    const auto prob_lo = lev_geom.ProbLoArray();

    // NOTE: see HurricaneEyeTrackerInitial -- the arg-min is split into a
    //       reducing pass and a locating pass so that the recorded (i,j)
    //       always belongs to the cell holding the recorded minimum.
    const Dim3 dlo = lbound(lev_geom.Domain());
    const int nx = lev_geom.Domain().length(0);

    Gpu::DeviceScalar<Long> d_idx_min(std::numeric_limits<Long>::max());
    Long* d_idx_min_ptr = d_idx_min.dataPtr();

    for (MFIter mfi(S_data[IntVars::cons]); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Array4<Real const>& S_arr = S_data[IntVars::cons].const_array(mfi);

        ParallelFor(box,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if(k==0) {
                Real x =  prob_lo[0] + (i+myhalf)*dx[0];
                Real y =  prob_lo[1] + (j+myhalf)*dx[1];
                Real dist = std::sqrt((x-tmp_x_eye)*(x-tmp_x_eye) + (y-tmp_y_eye)*(y-tmp_y_eye));
                if(dist < 200e3) {
                    Real qv_for_p = (use_moisture && (ncomp > RhoQ1_comp)) ? S_arr(i,j,k,RhoQ1_comp)/S_arr(i,j,k,Rho_comp) : 0;
                    const Real rhotheta = S_arr(i,j,k,RhoTheta_comp);
                    Real pressure = getPgivenRTh(rhotheta,qv_for_p);
                    Gpu::Atomic::Min(&d_val_min_ptr[0], pressure);
                }
            }
        });
    }

    // The minimum over every box on this rank must be known before the
    // locating pass below can test against it.
    Gpu::synchronize();

    for (MFIter mfi(S_data[IntVars::cons]); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Array4<Real const>& S_arr = S_data[IntVars::cons].const_array(mfi);

        ParallelFor(box,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if(k==0) {
                Real x =  prob_lo[0] + (i+myhalf)*dx[0];
                Real y =  prob_lo[1] + (j+myhalf)*dx[1];
                Real dist = std::sqrt((x-tmp_x_eye)*(x-tmp_x_eye) + (y-tmp_y_eye)*(y-tmp_y_eye));
                if(dist < 200e3) {
                    Real qv_for_p = (use_moisture && (ncomp > RhoQ1_comp)) ? S_arr(i,j,k,RhoQ1_comp)/S_arr(i,j,k,Rho_comp) : 0;
                    const Real rhotheta = S_arr(i,j,k,RhoTheta_comp);
                    Real pressure = getPgivenRTh(rhotheta,qv_for_p);
                    if (pressure == d_val_min_ptr[0]) {
                        Gpu::Atomic::Min(d_idx_min_ptr, pack_ij(i,j,nx,dlo));
                    }
                }
            }
        });
    }

    Gpu::synchronize();

    // Unpack the located index back into the (i,j) device scalars that
    // ComputeGlobalMinLocation reads.
    {
        int h_i_min, h_j_min;
        unpack_ij(d_idx_min.dataValue(), nx, dlo, h_i_min, h_j_min);
        Gpu::copy(Gpu::hostToDevice, &h_i_min, &h_i_min + 1, d_i_min_ptr);
        Gpu::copy(Gpu::hostToDevice, &h_j_min, &h_j_min + 1, d_j_min_ptr);
    }

    Real global_val_min;
    int global_i_min, global_j_min;

    ComputeGlobalMinLocation(sc, lev_geom, S_data,
                             d_val_min_ptr, d_i_min_ptr, d_j_min_ptr,
                             global_val_min, global_i_min, global_j_min);
}

/**
 * Read hurricane tracking history from restart files.
 */
void
ERF::ReadStormTrackerRestart ()
{
    hurricane_eye_track_xy.clear();
    hurricane_eye_track_latlon.clear();
    hurricane_maxvel_vs_time.clear();
    hurricane_minpressure_vs_time.clear();

    const fs::path base_dir("Output_StormTracker");

    // Nothing to do for a fresh run.
    if (!fs::exists(base_dir)) {
        return;
    }
    //
    // Return the alphabetically last file in a directory.
    // Since the filenames are zero-padded, this is also the newest output.
    //
    auto last_file = [](const fs::path& dir) -> fs::path
    {
        std::vector<fs::path> files;

        for (const auto& entry : fs::directory_iterator(dir)) {
            if (entry.is_regular_file()) {
                files.push_back(entry.path());
            }
        }

        if (files.empty()) {
            return fs::path{};
        }

        std::sort(files.begin(), files.end());

        return files.back();
    };

    //==========================================================
    // Read lat/lon file
    //==========================================================

    {
        fs::path file = last_file(base_dir / "latlon");

        if (!file.empty())
        {
            std::ifstream ifs(file);

            if (!ifs.is_open()) {
                Abort("Could not open " + file.string());
            }

            std::string line;

            // Skip the header line.
            std::getline(ifs, line);

            Real lat, lon;

            while (ifs >> lat >> lon)
            {
                hurricane_eye_track_latlon.push_back({lat, lon});
            }
        }
    }

    //==========================================================
    // Read maxvel tracker file
    //==========================================================

    {
        fs::path file = last_file(base_dir / "maxvel");

        if (!file.empty())
        {
            std::ifstream ifs(file);

            if (!ifs.is_open()) {
                Abort("Could not open " + file.string());
            }

            std::string line;

            // Skip the header line.
            std::getline(ifs, line);

            amrex::Real val1, val2;

            while (ifs >> val1 >> val2)
            {
                hurricane_maxvel_vs_time.push_back({val1, val2});
            }
        }
    }

    //==========================================================
    // Read minpressure tracker file
    //==========================================================

    {
        fs::path file = last_file(base_dir / "minpressure");

        if (!file.empty())
        {
            std::ifstream ifs(file);

            if (!ifs.is_open()) {
                Abort("Could not open " + file.string());
            }

            std::string line;

            // Skip the header line.
            std::getline(ifs, line);

            amrex::Real val1, val2;

            while (ifs >> val1 >> val2)
            {
                hurricane_minpressure_vs_time.push_back({val1, val2});
            }
        }
    }

    //==========================================================
    // Read XY VTK file
    //==========================================================

    {
        fs::path file = last_file(base_dir / "xy");

        if (!file.empty())
        {
            std::ifstream ifs(file);
            std::string line;

            // Skip the first four header lines.
            for (int i = 0; i < 4; ++i) {
                std::getline(ifs, line);
            }


            std::getline(ifs, line);
            std::istringstream iss(line);
            std::string keyword, datatype;
            int npoints;

            iss >> keyword >> npoints >> datatype;
            hurricane_eye_track_xy.reserve(npoints);

            for (int i = 0; i < npoints; ++i)
            {
                Real x, y, z;
                ifs >> x >> y >> z;
                hurricane_eye_track_xy.push_back({x, y});
            }
        }
    }
}

/**
 * Wrapper to track the hurricane eye position over time.
 *
 * @param[in] sc Solver choices
 */
void
ERF::HurricaneEyeTracker (const SolverChoice& sc)
{
    static bool is_start = true;
    int levc=finest_level;

    const MoistureType moisture_type   = sc.moisture_type;
    const Real hurricane_eye_latitude  = sc.hurricane_eye_latitude;
    const Real hurricane_eye_longitude = sc.hurricane_eye_longitude;

    if(is_start and restart_chkfile.empty()){
        HurricaneEyeTrackerInitial(sc, geom[levc],
                                   vars_new[levc],
                                   hurricane_eye_latitude,
                                   hurricane_eye_longitude);
        is_start = false;
    } else {
         if(!restart_chkfile.empty()) {
            ReadStormTrackerRestart();
        }
        HurricaneEyeTrackerNotInitial(sc, geom[levc], vars_new[levc],
                                      moisture_type);
    }
    HurricaneTrackerCircle();
}

/**
 * Compute and track the maximum wind velocity near the hurricane eye.
 *
 * @param[in] lev_geom Geometry of the current level
 * @param[in] mf_cc_vel MultiFab containing cell-centered velocity
 * @param[in] time Current simulation time
 */
void
ERF::HurricaneMaxVelTracker(const Geometry& lev_geom,
                            const MultiFab& mf_cc_vel,
                            const double& time)
{
    const int ncomp = AMREX_SPACEDIM;

    Real* d_val_max_ptr;
    Gpu::DeviceVector<Real> d_val_max(1, -bogus_large_value);
    d_val_max_ptr = d_val_max.data();

    const auto [x_last, y_last] = hurricane_eye_track_xy.back();
    const auto dx = lev_geom.CellSizeArray();
    const auto prob_lo = lev_geom.ProbLoArray();

    Real x_eye = x_last;
    Real y_eye = y_last;

    for (MFIter mfi(mf_cc_vel); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const auto& vel_arr = mf_cc_vel.const_array(mfi);

        ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            Real x = prob_lo[0] + (i+myhalf)*dx[0];
            Real y = prob_lo[1] + (j+myhalf)*dx[1];
            Real dist = std::sqrt((x-x_eye)*(x-x_eye) +
                                         (y-y_eye)*(y-y_eye));
            if(k==1 && dist < 200e3) {
                Real velmag = zero;
                for (int comp = 0; comp < ncomp; ++comp) {
                    Real vel = vel_arr(i, j, k, comp);
                    velmag += vel * vel;
                }
                velmag = std::sqrt(velmag)*Real(3.6); // km/hr
                Gpu::Atomic::Max(&d_val_max_ptr[0], velmag);
            }
        });
    }

    Gpu::synchronize();

    Real h_val_max_local = -bogus_large_value;
    Gpu::copy(Gpu::deviceToHost, d_val_max.begin(), d_val_max.end(), &h_val_max_local);

    Real h_val_max_global = -bogus_large_value;
    #ifdef AMREX_USE_MPI
        h_val_max_global = h_val_max_local;
        amrex::ParallelDescriptor::ReduceRealMax(h_val_max_global);
    #else
        h_val_max_global = h_val_max_local;
    #endif

    double time_in_hrs = time / 3600.0;
    hurricane_maxvel_vs_time.push_back({static_cast<Real>(time_in_hrs), h_val_max_global});
}

/**
 * Compute and track the minimum pressure near the hurricane eye.
 *
 * @param[in] moisture_type Moisture model type
 * @param[in] lev_geom Geometry of the current level
 * @param[in] mf_cons_var MultiFab containing conservative variables
 * @param[in] time Current simulation time
 */
void
ERF::HurricaneMinPressureTracker (MoistureType moisture_type,
                                  const Geometry& lev_geom,
                                  const MultiFab& mf_cons_var,
                                  const double& time)
{


    Real* d_val_min_ptr;
    Gpu::DeviceVector<Real> d_val_min(1, bogus_large_value);
    d_val_min_ptr = d_val_min.data();

    const Real x_last = hurricane_eye_track_xy.back()[0];
    const Real y_last = hurricane_eye_track_xy.back()[1];
    const auto dx = lev_geom.CellSizeArray();
    const auto prob_lo = lev_geom.ProbLoArray();

    const int ncomp = mf_cons_var.nComp();
    bool use_moisture = (moisture_type != MoistureType::None);

    for (MFIter mfi(mf_cons_var); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const auto& S_arr = mf_cons_var.const_array(mfi);

        ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            Real x = prob_lo[0] + (i+myhalf)*dx[0];
            Real y = prob_lo[1] + (j+myhalf)*dx[1];
            Real dist2 = (x-x_last)*(x-x_last) +
                                (y-y_last)*(y-y_last);
            if(k==1 && dist2 < 200e3*200e3) {
                const Real rhotheta = S_arr(i,j,k,RhoTheta_comp);
                const Real qv_for_p = (use_moisture && (ncomp > RhoQ1_comp)) ? S_arr(i,j,k,RhoQ1_comp)/S_arr(i,j,k,Rho_comp) : 0;
                const Real pressure = getPgivenRTh(rhotheta,qv_for_p);
                Gpu::Atomic::Min(&d_val_min_ptr[0], pressure);
            }
        });
    }

    Gpu::synchronize();

    Real h_val_min_local = bogus_large_value;
    Gpu::copy(Gpu::deviceToHost, d_val_min.begin(), d_val_min.end(), &h_val_min_local);

    // NOTE: use the amrex wrapper rather than a hard-coded MPI_DOUBLE, which is
    //       a type mismatch when amrex::Real is float (ERF_PRECISION=SINGLE).
    Real h_val_min_global = h_val_min_local;
    ParallelDescriptor::ReduceRealMin(h_val_min_global);

    double time_in_hrs = time / 3600.0;
    hurricane_minpressure_vs_time.push_back({static_cast<Real>(time_in_hrs), h_val_min_global});
}
#endif
