#include <ERF_FireTerrainReader.H>
#include <fstream>
#include <vector>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Gpu.H>

using namespace amrex;

bool read_terrain_onto_fire_grid(
    MultiFab&   z_fire_nd,
    const FireGrid&    fg,
    const std::string& fname)
{
    // Return false if fname is empty
    if (fname.empty()) {
        return false;
    }

    int C = fg.C;
    int nx_terrain = 0;
    int ny_terrain = 0;
    std::vector<Real> x_coords;
    std::vector<Real> y_coords;
    std::vector<Real> z_values;

    // IO rank reads the file
    if (ParallelDescriptor::IOProcessor()) {
        std::ifstream file(fname);
        if (!file.is_open()) {
            amrex::Warning("Could not open terrain file: " + fname);
            return false;
        }

        // Read nx and ny
        file >> nx_terrain >> ny_terrain;

        x_coords.resize(nx_terrain);
        y_coords.resize(ny_terrain);
        z_values.resize(nx_terrain * ny_terrain);

        // Read x-coordinates
        for (int i = 0; i < nx_terrain; ++i) {
            file >> x_coords[i];
        }

        // Read y-coordinates
        for (int j = 0; j < ny_terrain; ++j) {
            file >> y_coords[j];
        }

        // Read z-values (stored contiguous in y: z[ix*ny + iy])
        for (int i = 0; i < nx_terrain * ny_terrain; ++i) {
            file >> z_values[i];
        }

        file.close();
    }

    // Broadcast nx, ny to all ranks
    ParallelDescriptor::Bcast(&nx_terrain, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&ny_terrain, 1, ParallelDescriptor::IOProcessorNumber());

    // Resize on non-IO ranks
    if (!ParallelDescriptor::IOProcessor()) {
        x_coords.resize(nx_terrain);
        y_coords.resize(ny_terrain);
        z_values.resize(nx_terrain * ny_terrain);
    }

    // Broadcast the arrays
    ParallelDescriptor::Bcast(x_coords.data(), nx_terrain, 
                              ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(y_coords.data(), ny_terrain, 
                              ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(z_values.data(), nx_terrain * ny_terrain, 
                              ParallelDescriptor::IOProcessorNumber());

    // Copy to GPU
    //Gpu::DeviceVector<Real> x_device(x_coords.begin(), x_coords.end());
    //Gpu::DeviceVector<Real> y_device(y_coords.begin(), y_coords.end());
    //Gpu::DeviceVector<Real> z_device(z_values.begin(), z_values.end());

    // Copy to GPU
    Gpu::DeviceVector<Real> x_device(nx_terrain);
    Gpu::DeviceVector<Real> y_device(ny_terrain);
    Gpu::DeviceVector<Real> z_device(nx_terrain * ny_terrain);

    Gpu::copy(Gpu::hostToDevice, x_coords.begin(), x_coords.end(), x_device.begin());
    Gpu::copy(Gpu::hostToDevice, y_coords.begin(), y_coords.end(), y_device.begin());
    Gpu::copy(Gpu::hostToDevice, z_values.begin(), z_values.end(), z_device.begin());

    Real* x_ptr = x_device.data();
    Real* y_ptr = y_device.data();
    Real* z_ptr = z_device.data();

    // Problem domain
    Real ProbLo_x = fg.geom.ProbLo(0);
    Real ProbLo_y = fg.geom.ProbLo(1);

    // Fire grid spacing
    auto dx_dy = fg.geom.CellSizeArray();
    Real dx_f = dx_dy[0];
    Real dy_f = dx_dy[1];

    // Get fire domain bounds (for clamping)
    const Box& domain_fire = fg.geom.Domain();
    int i_fire_lo = domain_fire.smallEnd(0);
    int i_fire_hi = domain_fire.bigEnd(0);
    int j_fire_lo = domain_fire.smallEnd(1);
    int j_fire_hi = domain_fire.bigEnd(1);

    // Get terrain domain bounds
    Real x_min = x_coords[0];
    Real x_max = x_coords[nx_terrain - 1];
    Real y_min = y_coords[0];
    Real y_max = y_coords[ny_terrain - 1];

    // Interpolate onto fire grid
    for (MFIter mfi(z_fire_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> z_fire = z_fire_nd.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            int i_f = iv[0];
            int j_f = iv[1];

            // Clamp to domain for ghost cells
            int i_f_clamped = std::max(i_fire_lo, std::min(i_fire_hi, i_f));
            int j_f_clamped = std::max(j_fire_lo, std::min(j_fire_hi, j_f));

            // Compute physical position at fire-grid node
            Real x = ProbLo_x + i_f_clamped * dx_f;
            Real y = ProbLo_y + j_f_clamped * dy_f;

            // Clamp to terrain domain
            x = std::max(x_min, std::min(x_max, x));
            y = std::max(y_min, std::min(y_max, y));

            // Find surrounding terrain grid points via binary search or linear scan
            // For efficiency, use a linear scan (could optimize with binary search)
            int ix_lo = nx_terrain - 2;  // safe fallback for x == x_max
            int iy_lo = ny_terrain - 2;  // safe fallback for y == y_max

            // Find x index
            for (int ix = 0; ix < nx_terrain - 1; ++ix) {
                if (x_ptr[ix] <= x && x <= x_ptr[ix + 1]) {
                    ix_lo = ix;
                    break;
                }
            }

            // Find y index
            for (int iy = 0; iy < ny_terrain - 1; ++iy) {
                if (y_ptr[iy] <= y && y <= y_ptr[iy + 1]) {
                    iy_lo = iy;
                    break;
                }
            }

            int ix_hi = std::min(ix_lo + 1, nx_terrain - 1);
            int iy_hi = std::min(iy_lo + 1, ny_terrain - 1);

            // Get four corner values
            Real z_ll = z_ptr[ix_lo * ny_terrain + iy_lo];
            Real z_lr = z_ptr[ix_hi * ny_terrain + iy_lo];
            Real z_ul = z_ptr[ix_lo * ny_terrain + iy_hi];
            Real z_ur = z_ptr[ix_hi * ny_terrain + iy_hi];

            // Bilinear interpolation weights
            Real dx_x = (x_ptr[ix_hi] > x_ptr[ix_lo]) ? (x - x_ptr[ix_lo]) / (x_ptr[ix_hi] - x_ptr[ix_lo]) : 0.0;
            Real dy_y = (y_ptr[iy_hi] > y_ptr[iy_lo]) ? (y - y_ptr[iy_lo]) / (y_ptr[iy_hi] - y_ptr[iy_lo]) : 0.0;

            dx_x = std::max(0.0, std::min(1.0, dx_x));
            dy_y = std::max(0.0, std::min(1.0, dy_y));

            // Bilinear interpolation
            Real z_lo = z_ll * (1.0 - dx_x) + z_lr * dx_x;
            Real z_hi = z_ul * (1.0 - dx_x) + z_ur * dx_x;
            Real z = z_lo * (1.0 - dy_y) + z_hi * dy_y;

            z_fire(iv, 0) = z;
        });
    }

    return true;
}
