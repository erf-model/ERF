#ifndef ERF_ACTUATOR_DEVICE_BUFFERS_H
#define ERF_ACTUATOR_DEVICE_BUFFERS_H

/** Device memory that the actuator spreading and sampling kernels keep from one call to the next. */

#include <algorithm>
#include <cstddef>
#include <vector>

#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>

namespace erf_actuator {

/** Grow v to hold at least n values; never shrink it, so that a buffer keeps its largest size. */
template <class V>
void grow_buffer (V& v, std::size_t n)
{
    if (v.size() < n) { v.resize(n); }
}

/**
 * An array filled on the host and read by the kernels: the host side is a std::vector (clear() keeps its
 * capacity), copied through pinned memory to the device with one asynchronous copy. The pinned and device
 * sides only grow. A call that uploads must synchronise the stream before it returns, so that the next
 * call may refill the host side.
 */
template <class T>
struct StagedBuffer
{
    std::vector<T> host;                    //!< what the next upload() sends
    amrex::Gpu::PinnedVector<T> pinned;     //!< staging, and the host side of download_async()
    amrex::Gpu::DeviceVector<T> device;

    /** Copy host to the device, asynchronously. */
    void upload ()
    {
        grow_buffer(pinned, host.size());
        grow_buffer(device, host.size());
        if (host.empty()) { return; }
        std::copy(host.begin(), host.end(), pinned.begin());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, pinned.begin(), pinned.begin() + host.size(), device.begin());
    }
    /** Copy the first n device values to pinned, asynchronously: they are there after the next stream synchronisation. */
    void download_async (std::size_t n)
    {
        grow_buffer(pinned, n);
        if (n == 0) { return; }
        amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, device.begin(), device.begin() + n, pinned.begin());
    }
    T* data () { return device.data(); }
};

/**
 * The buffers of the spreading kernels (spread_partial_norms, add_scaled_kernels) and of the sampler
 * (sample_velocity_local), kept by their caller so that a step allocates no
 * device memory once the buffers have reached their size. Passing none makes a call use its own.
 */
struct ActuatorDeviceBuffers
{
    // spreading, one entry per periodic image of a point within reach of the domain
    StagedBuffer<amrex::Real> image_pos;     //!< 3 per image: its position (m)
    StagedBuffer<int> image_point;           //!< 1 per image: the point it is an image of
    StagedBuffer<int> image_reach;           //!< 18 per image: per face grid, the faces within its reach (lo x y z, hi x y z)
    // spreading, one entry per point
    StagedBuffer<amrex::Real> point_eps;     //!< 2 per point: 1 / eps^2 and the squared cut-off radius 9 eps^2
    StagedBuffer<amrex::Real> scales;        //!< 3 per point: F_p / S_p (add_scaled_kernels)
    StagedBuffer<amrex::Real> norms;         //!< 3 per point: the accumulators of S_p (spread_partial_norms)
    // spreading, per (box, face grid): which images reach its faces, in point then image order
    StagedBuffer<int> bins;
    // sampling, per point
    StagedBuffer<amrex::Real> sample_pos;    //!< 3 per point (m)
    StagedBuffer<amrex::Real> sample_vel;    //!< 3 per point (m/s)
    StagedBuffer<int> sample_cnt;            //!< 1 per point: how many of this rank's boxes hold it

    long kernel_launches{0};                 //!< kernels launched through these buffers so far (for the tests)
    /** A switch for the unit tests only; a model leaves it false. When true, add_scaled_kernels() runs the
     * face-by-face pass on the host too, which is the only pass a device has, so a CPU build can compare the
     * two. Both give the same sources bit for bit; the face-by-face pass is the slower one on a host. */
    bool force_face_major{false};
};

} // namespace erf_actuator

#endif
