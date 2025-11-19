#ifndef ERF_HURRICANE_DIAGNOSTICS_H_
#define ERF_HURRICANE_DIAGNOSTICS_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelReduce.H>
#include <limits>

#include "ERF_DataStruct.H"

/**
 * Routines to compute hurricane diagnostics
 */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct {
    amrex::Real value;
    int rank;
    } in, out;

void ComputeGlobalMinLocation(const amrex::Geometry& geom,
                              const amrex::Vector<amrex::MultiFab>& S_data,
                              const amrex::Vector<amrex::MultiFab>* forecast_state_at_lev,
                              amrex::Real* d_val_min_ptr,
                              int* d_i_min_ptr,
                              int* d_j_min_ptr,
                              amrex::Real& global_val_min,
                              int& global_i_min,
                              int& global_j_min,
                              amrex::Vector<std::array<amrex::Real, 2>>& hurricane_eye_track_xy,
                              amrex::Vector<std::array<amrex::Real, 2>>& hurricane_eye_track_latlon)
{
    amrex::Real h_val_min;
    int h_i_min, h_j_min;

    amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_val_min_ptr, d_val_min_ptr + 1, &h_val_min);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_i_min_ptr, d_i_min_ptr + 1, &h_i_min);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_j_min_ptr, d_j_min_ptr + 1, &h_j_min);
    amrex::Gpu::synchronize();

    amrex::Real local_val_min = h_val_min;
    int local_i_min = h_i_min;
    int local_j_min = h_j_min;

    int rank = amrex::ParallelDescriptor::MyProc();

    in.value = local_val_min;
    in.rank  = rank;

    #ifdef AMREX_USE_MPI
        MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    #else
        out = in;
    #endif

    global_val_min = out.value;
    int owner_rank = out.rank;

    // Broadcast the indices from the rank that owns the minimum
    global_i_min = local_i_min;
    global_j_min = local_j_min;

    amrex::ParallelDescriptor::Bcast(&global_i_min, 1, owner_rank);
    amrex::ParallelDescriptor::Bcast(&global_j_min, 1, owner_rank);

    if (rank == 0) {
        amrex::Print() << "Global minimum distance to hurricane eye (k=0): "
                       << global_val_min << " at (i,j) = ("
                       << global_i_min << ", " << global_j_min << ")\n";
    }

    amrex::Gpu::DeviceScalar<amrex::Real> d_eye_lat(0.0), d_eye_lon(0.0);

    amrex::Real* d_eye_lat_ptr = d_eye_lat.dataPtr();
    amrex::Real* d_eye_lon_ptr = d_eye_lon.dataPtr();

    // On owner_rank, compute eye_lat and eye_lon
    if (rank == owner_rank) {
        for (amrex::MFIter mfi(S_data[IntVars::cons]); mfi.isValid(); ++mfi) {
            const amrex::Box& box = mfi.validbox();
            const auto latlon_arr = (*forecast_state_at_lev)[4].array(mfi);
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (i == global_i_min && j == global_j_min && k == 0) {
                    *d_eye_lat_ptr = latlon_arr(i,j,k,0);
                    *d_eye_lon_ptr = latlon_arr(i,j,k,1);
                }
            });
        }
    }

    amrex::Real eye_lat = d_eye_lat.dataValue();
    amrex::Real eye_lon = d_eye_lon.dataValue();

    // Synchronize to ensure the owner has computed values
    amrex::Gpu::synchronize();

    amrex::ParallelDescriptor::Bcast(&eye_lat, 1, owner_rank);
    amrex::ParallelDescriptor::Bcast(&eye_lon, 1, owner_rank);

    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();

    amrex::Real eye_x =  prob_lo[0] + (global_i_min+0.5)*dx[0];
    amrex::Real eye_y =  prob_lo[1] + (global_j_min+0.5)*dx[1];

    hurricane_eye_track_xy.push_back({eye_x, eye_y});
    hurricane_eye_track_latlon.push_back({eye_lon, eye_lat});
}

void
HurricaneTrackerCircle (const amrex::Vector<std::array<amrex::Real, 2>>& hurricane_eye_track_xy,
                        amrex::Vector<std::array<amrex::Real, 2>>& hurricane_tracker_circle)
{
    // Check that there is at least one eye position
    if (hurricane_eye_track_xy.empty()) return;

    // Get the last known (x, y) position of the eye
    const auto [x_last, y_last] = hurricane_eye_track_xy.back();

    // Define circle properties
    const int n_points = 100;        // number of points on the circle
    const amrex::Real radius = 200e3; // radius in meters (example: 50 km)

    // Clear previous points and reserve space
    hurricane_tracker_circle.clear();
    hurricane_tracker_circle.reserve(n_points);

    // Fill the circle points
    for (int i = 0; i < n_points; ++i) {
        amrex::Real theta = 2.0 * M_PI * i / n_points;
        amrex::Real x = x_last + radius * std::cos(theta);
        amrex::Real y = y_last + radius * std::sin(theta);
        hurricane_tracker_circle.push_back({x, y});
    }
}

void HurricaneEyeTrackerInitial (const amrex::Geometry& geom,
                                 const amrex::Vector<amrex::MultiFab>& S_data,
                                 const amrex::Vector<amrex::MultiFab>* forecast_state_at_lev,
                                 amrex::Vector<std::array<amrex::Real, 2>>& hurricane_eye_track_xy,
                                 amrex::Vector<std::array<amrex::Real, 2>>& hurricane_eye_track_latlon,
                                 const amrex::Real& hurricane_eye_latitude,
                                 const amrex::Real& hurricane_eye_longitude)
{
    amrex::Gpu::DeviceScalar<amrex::Real> d_val_min(1e10);
    amrex::Gpu::DeviceScalar<int> d_i_min(-1), d_j_min(-1);

    amrex::Real* d_val_min_ptr = d_val_min.dataPtr();
    int* d_i_min_ptr = d_i_min.dataPtr();
    int* d_j_min_ptr = d_j_min.dataPtr();

    for (amrex::MFIter mfi(S_data[IntVars::cons]); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.validbox();
        const auto latlon_arr = (*forecast_state_at_lev)[4].array(mfi);

        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if (k==0) {

                amrex::Real dlat = latlon_arr(i,j,k,0) - hurricane_eye_latitude;
                amrex::Real dlon = latlon_arr(i,j,k,1) - hurricane_eye_longitude;
                amrex::Real dist = std::sqrt(dlat*dlat + dlon*dlon);
                // Atomic min using device pointer from DeviceVector
                amrex::Real old = amrex::Gpu::Atomic::Min(&d_val_min_ptr[0], dist);
                //amrex::Gpu::Atomic::Min(&d_val_min_ptr[0], dist);
                if (dist < old) {
                    // We are the new minimum; record indices
                    d_i_min_ptr[0] = i;
                    d_j_min_ptr[0] = j;
                }
            }
        });
    }

    amrex::Real global_val_min;
    int global_i_min, global_j_min;

    ComputeGlobalMinLocation(geom, S_data, forecast_state_at_lev,
                             d_val_min_ptr, d_i_min_ptr, d_j_min_ptr,
                             global_val_min, global_i_min, global_j_min,
                             hurricane_eye_track_xy,
                             hurricane_eye_track_latlon);
}

void HurricaneEyeTrackerNotInitial (const amrex::Geometry& geom,
                                    const amrex::Vector<amrex::MultiFab>& S_data,
                                    MoistureType moisture_type,
                                    const amrex::Vector<amrex::MultiFab>* forecast_state_at_lev,
                                    amrex::Vector<std::array<amrex::Real, 2>>& hurricane_eye_track_xy,
                                    amrex::Vector<std::array<amrex::Real, 2>>& hurricane_eye_track_latlon)
{

    if (hurricane_eye_track_xy.empty()) {
        amrex::Print() << "Error: hurricane_eye_track_xy is empty!\n";
        amrex::Abort("Attempted to access hurricane_eye_track_xy[0]");
    }

    amrex::Real tmp_x_eye = hurricane_eye_track_xy.back()[0];
    amrex::Real tmp_y_eye = hurricane_eye_track_xy.back()[1];

    if(amrex::ParallelDescriptor::IOProcessor()){
        std::cout << "The value of x y are " << tmp_x_eye << " " << tmp_y_eye << std::endl;
    }

    amrex::Gpu::DeviceScalar<amrex::Real> d_val_min(1e10);
    amrex::Gpu::DeviceScalar<int> d_i_min(-1), d_j_min(-1);

    amrex::Real* d_val_min_ptr = d_val_min.dataPtr();
    int* d_i_min_ptr = d_i_min.dataPtr();
    int* d_j_min_ptr = d_j_min.dataPtr();

    bool use_moisture = (moisture_type != MoistureType::None);
    const int ncomp = S_data[IntVars::cons].nComp();

    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();

    for (amrex::MFIter mfi(S_data[IntVars::cons]); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.validbox();
        const amrex::Array4<amrex::Real const>& S_arr = S_data[IntVars::cons].const_array(mfi);

        amrex::ParallelFor(box,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if(k==0) {
                amrex::Real x =  prob_lo[0] + (i+0.5)*dx[0];
                amrex::Real y =  prob_lo[1] + (j+0.5)*dx[1];
                amrex::Real dist = std::sqrt((x-tmp_x_eye)*(x-tmp_x_eye) + (y-tmp_y_eye)*(y-tmp_y_eye));
                if(dist < 200e3) {
                    amrex::Real qv_for_p = (use_moisture && (ncomp > RhoQ1_comp)) ? S_arr(i,j,k,RhoQ1_comp)/S_arr(i,j,k,Rho_comp) : 0;
                    const amrex::Real rhotheta = S_arr(i,j,k,RhoTheta_comp);
                    amrex::Real pressure = getPgivenRTh(rhotheta,qv_for_p);
                    amrex::Real old = amrex::Gpu::Atomic::Min(&d_val_min_ptr[0], pressure);
                    //amrex::Gpu::Atomic::Min(&d_val_min_ptr[0], dist);
                    if (old > pressure) {
                        // We are the new minimum; record indices
                        d_i_min_ptr[0] = i;
                        d_j_min_ptr[0] = j;
                    }
                }
            }
        });
    }

    amrex::Real global_val_min;
    int global_i_min, global_j_min;

    ComputeGlobalMinLocation(geom, S_data, forecast_state_at_lev,
                             d_val_min_ptr, d_i_min_ptr, d_j_min_ptr,
                             global_val_min, global_i_min, global_j_min,
                             hurricane_eye_track_xy,
                             hurricane_eye_track_latlon);
}

void HurricaneEyeTracker (const amrex::Geometry& geom,
                          const amrex::Vector<amrex::MultiFab>& S_data,
                          MoistureType moisture_type,
                          const amrex::Vector<amrex::MultiFab>* forecast_state_at_lev,
                          const amrex::Real& hurricane_eye_latitude,
                          const amrex::Real& hurricane_eye_longitude,
                          amrex::Vector<std::array<amrex::Real, 2>>& hurricane_eye_track_xy,
                          amrex::Vector<std::array<amrex::Real, 2>>& hurricane_eye_track_latlon,
                          amrex::Vector<std::array<amrex::Real, 2>>& hurricane_tracker_circle)
{
    static bool is_start = true;
    if(is_start){
        HurricaneEyeTrackerInitial(geom, S_data,
                                   forecast_state_at_lev, hurricane_eye_track_xy,
                                   hurricane_eye_track_latlon,
                                   hurricane_eye_latitude, hurricane_eye_longitude);
        is_start = false;
    } else {
        HurricaneEyeTrackerNotInitial(geom, S_data, moisture_type,
                                      forecast_state_at_lev, hurricane_eye_track_xy,
                                      hurricane_eye_track_latlon);
    }
    HurricaneTrackerCircle(hurricane_eye_track_xy, hurricane_tracker_circle);
}

void
HurricaneMaxVelTracker(const amrex::Geometry& geom,
                       const amrex::MultiFab& mf_cc_vel,
                       const amrex::Real& time,
                       const amrex::Vector<std::array<amrex::Real, 2>>& hurricane_eye_track_xy,
                       amrex::Vector<std::array<amrex::Real, 2>>& hurricane_maxvel_vs_time)
{
    const int ncomp = AMREX_SPACEDIM;

    amrex::Real* d_val_max_ptr;
    amrex::Gpu::DeviceVector<amrex::Real> d_val_max(1, -1.0e30);
    d_val_max_ptr = d_val_max.data();

    const auto [x_last, y_last] = hurricane_eye_track_xy.back();
    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();

    amrex::Real x_eye = x_last;
    amrex::Real y_eye = y_last;

    for (amrex::MFIter mfi(mf_cc_vel); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.validbox();
        const auto& vel_arr = mf_cc_vel.const_array(mfi);

        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];
            amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];
            amrex::Real dist = std::sqrt((x-x_eye)*(x-x_eye) +
                                         (y-y_eye)*(y-y_eye));
            if(k==1 && dist < 200e3) {
                amrex::Real velmag = 0.0;
                for (int comp = 0; comp < ncomp; ++comp) {
                    amrex::Real vel = vel_arr(i, j, k, comp);
                    velmag += vel * vel;
                }
                velmag = std::sqrt(velmag)*3.6; // km/hr
                amrex::Gpu::Atomic::Max(&d_val_max_ptr[0], velmag);
            }
        });
    }

    amrex::Gpu::synchronize();

    amrex::Real h_val_max_local = -1.0e30;
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_val_max.begin(), d_val_max.end(), &h_val_max_local);

    amrex::Real h_val_max_global = -1.0e30;
     #ifdef AMREX_USE_MPI
        MPI_Allreduce(&h_val_max_local, &h_val_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    #else
        h_val_max_global = h_val_max_local;
    #endif

    amrex::Real time_in_hrs = time / 3600.0;
    hurricane_maxvel_vs_time.push_back({time_in_hrs, h_val_max_global});
}

#endif
