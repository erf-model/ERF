#ifndef ERF_HURRICANE_DIAGNOSTICS_H_
#define ERF_HURRICANE_DIAGNOSTICS_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelReduce.H>
#include <limits>

#include "ERF_DataStruct.H"
#include "ERF.H"

/**
 * Routines to compute hurricane diagnostics
 */

#ifndef M_PI
#define M_PI Real(3.14159265358979323846)
#endif

using namespace amrex;

struct {
    Real value;
    int rank;
    } in, out;


void
ERF::ComputeGlobalMinLocation_WRF (const Geometry& lev_geom,
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
    if (rank == owner_rank) {
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

void
ERF::HurricaneTrackerCircle_WRF ()
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

void
ERF::HurricaneEyeTrackerInitial_WRF (const Geometry& lev_geom,
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
                // Atomic min using device pointer from DeviceVector
                Real old = Gpu::Atomic::Min(&d_val_min_ptr[0], dist);
                //Gpu::Atomic::Min(&d_val_min_ptr[0], dist);
                if (dist < old) {
                    // We are the new minimum; record indices
                    d_i_min_ptr[0] = i;
                    d_j_min_ptr[0] = j;
                }
            }
        });
    }

    Real global_val_min;
    int global_i_min, global_j_min;

    ComputeGlobalMinLocation_WRF(lev_geom, S_data,
                                 d_val_min_ptr, d_i_min_ptr, d_j_min_ptr,
                                 global_val_min, global_i_min, global_j_min);
}

void
ERF::HurricaneEyeTrackerNotInitial_WRF (const Geometry& lev_geom,
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
                    Real old = Gpu::Atomic::Min(&d_val_min_ptr[0], pressure);
                    //Gpu::Atomic::Min(&d_val_min_ptr[0], dist);
                    if (old > pressure) {
                        // We are the new minimum; record indices
                        d_i_min_ptr[0] = i;
                        d_j_min_ptr[0] = j;
                    }
                }
            }
        });
    }

    Real global_val_min;
    int global_i_min, global_j_min;

    ComputeGlobalMinLocation_WRF (lev_geom, S_data,
                             d_val_min_ptr, d_i_min_ptr, d_j_min_ptr,
                             global_val_min, global_i_min, global_j_min);
}

void
ERF::HurricaneEyeTracker_WRF (const SolverChoice& sc)
{
    static bool is_start = true;
    int levc=finest_level;

    const MoistureType moisture_type   = sc.moisture_type;
    const Real hurricane_eye_latitude  = sc.hurricane_eye_latitude;
    const Real hurricane_eye_longitude = sc.hurricane_eye_longitude;

    if(is_start){
        HurricaneEyeTrackerInitial_WRF(geom[levc],
                                       vars_new[levc],
                                       hurricane_eye_latitude,
                                       hurricane_eye_longitude);
        is_start = false;
    } else {
        HurricaneEyeTrackerNotInitial_WRF(geom[levc], vars_new[levc],
                                          moisture_type);
    }
    HurricaneTrackerCircle_WRF();
}

#endif
