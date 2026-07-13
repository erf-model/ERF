#include <ERF.H>
#include <ERF_Derive.H>

using namespace amrex;

void
tag_on_distance_from_eye(const Geometry& cgeom, TagBoxArray* tags,
                         const Real eye_x, const Real eye_y, const Real rad_tag);

bool
ERF::FindInitialEye(int levc,
                    const MultiFab& mf_cc_vel,
                    const Real velmag_threshold,
                    Real& eye_x, Real& eye_y)
{
    const auto dx = geom[levc].CellSizeArray();
    const auto prob_lo = geom[levc].ProbLoArray();

    Gpu::DeviceVector<Real> d_coords(2, zero);
    Gpu::DeviceVector<int>  d_found(1,0);

    Real* d_coords_ptr = d_coords.data();
    int*   d_found_ptr = d_found.data();

    for (MFIter mfi(mf_cc_vel); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const Array4<const Real>& vel_arr = mf_cc_vel.const_array(mfi);

        ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Real magnitude = std::sqrt(vel_arr(i,j,k,0) * vel_arr(i,j,k,0) +
                                       vel_arr(i,j,k,1) * vel_arr(i,j,k,1) +
                                       vel_arr(i,j,k,2) * vel_arr(i,j,k,2));

            magnitude *= Real(3.6);

            Real z = prob_lo[2] + (k + myhalf) * dx[2];

            // Check if magnitude exceeds threshold
            if (z < Real(2000.) && magnitude > velmag_threshold) {
                // Use atomic operations to set found flag and store coordinates
                Gpu::Atomic::Add(&d_found_ptr[0], 1); // Mark as found

                Real x = prob_lo[0] + (i + myhalf) * dx[0];
                Real y = prob_lo[1] + (j + myhalf) * dx[1];

                // Store coordinates
                Gpu::Atomic::Add(&d_coords_ptr[0],x); // Store x index
                Gpu::Atomic::Add(&d_coords_ptr[1],y); // Store x index
            }
        });
    }

    // Synchronize to ensure all threads complete their execution
    amrex::Gpu::streamSynchronize(); // Wait for all GPU threads to finish

    Vector<int> h_found(1,0);
    Gpu::copy(Gpu::deviceToHost, d_found.begin(), d_found.end(), h_found.begin());
    ParallelAllReduce::Sum(h_found.data(), h_found.size(), ParallelContext::CommunicatorAll());

    // Broadcast coordinates if found
    if (h_found[0] > 0) {
        Vector<Real> h_coords(2,-bogus_large_value);
        Gpu::copy(Gpu::deviceToHost, d_coords.begin(), d_coords.end(), h_coords.begin());

        ParallelAllReduce::Sum(h_coords.data(), h_coords.size(), ParallelContext::CommunicatorAll());

        eye_x = h_coords[0]/h_found[0];
        eye_y = h_coords[1]/h_found[0];

    } else {
        // Random large negative numbers so we don't trigger refinement in this case
        eye_x = -bogus_large_value;
        eye_y = -bogus_large_value;
    }

    return (h_found[0] > 0);
}

void
tag_on_distance_from_eye(const Geometry& cgeom, TagBoxArray* tags,
                         const Real eye_x, const Real eye_y, const Real rad_tag)
{
    const auto dx      = cgeom.CellSizeArray();
    const auto prob_lo = cgeom.ProbLoArray();

    for (MFIter mfi(*tags); mfi.isValid(); ++mfi) {
        TagBox& tag = (*tags)[mfi];
        auto tag_arr = tag.array();  // Get device-accessible array

        const Box& tile_box = mfi.tilebox(); // The box for this tile

        ParallelFor(tile_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // Compute cell center coordinates
            Real x = prob_lo[0] + (i + myhalf) * dx[0];
            Real y = prob_lo[1] + (j + myhalf) * dx[1];

            Real dist = std::sqrt((x - eye_x)*(x - eye_x) + (y - eye_y)*(y - eye_y));

            if (dist < rad_tag) {
                tag_arr(i,j,k) = TagBox::SET;
            }
        });
    }
}

void
ERF::HurricaneTracker(int levc,
                      Real time,
                      const MultiFab& mf_cc_vel,
                      const Real velmag_threshold,
                      TagBoxArray* tags)
{
    bool is_found;

    Real eye_x, eye_y;

    if (time==zero || hurricane_eye_track_xy.empty()) {
        is_found = FindInitialEye(levc, mf_cc_vel, velmag_threshold, eye_x, eye_y);
    } else {
        is_found = true;
        const auto& last = hurricane_eye_track_xy.back();
        eye_x = last[0];
        eye_y = last[1];
    }

    if (is_found) {
        const int exponent = max_level-1-levc;
        Real rad_tag = std::ldexp(Real(4.e5), exponent);
        tag_on_distance_from_eye(geom[levc], tags, eye_x, eye_y, rad_tag);
    }
}
