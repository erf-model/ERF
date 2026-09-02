#include <ERF_DustTerrainSlope.H>
#include <AMReX_Gpu.H>
#include <AMReX_ParallelDescriptor.H>
#include <fstream>
#include <vector>


using namespace amrex;

namespace {

bool read_terrain_onto_dust_grid(
    MultiFab& z_dust_nd,
    const DustGrid& dg,
    const std::string& fname)
{
    if (fname.empty()) {
        return false;
    }

    int nx_terrain = 0;
    int ny_terrain = 0;
    std::vector<Real> x_coords;
    std::vector<Real> y_coords;
    std::vector<Real> z_values;

    if (ParallelDescriptor::IOProcessor()) {
        std::ifstream file(fname);
        if (!file.is_open()) {
            amrex::Warning("Could not open dust terrain file: " + fname);
            return false;
        }

        file >> nx_terrain >> ny_terrain;
        x_coords.resize(nx_terrain);
        y_coords.resize(ny_terrain);
        z_values.resize(nx_terrain * ny_terrain);

        for (int i = 0; i < nx_terrain; ++i) {
            file >> x_coords[i];
        }
        for (int j = 0; j < ny_terrain; ++j) {
            file >> y_coords[j];
        }
        for (int n = 0; n < nx_terrain * ny_terrain; ++n) {
            file >> z_values[n];
        }
    }

    ParallelDescriptor::Bcast(&nx_terrain, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&ny_terrain, 1, ParallelDescriptor::IOProcessorNumber());

    if (!ParallelDescriptor::IOProcessor()) {
        x_coords.resize(nx_terrain);
        y_coords.resize(ny_terrain);
        z_values.resize(nx_terrain * ny_terrain);
    }

    ParallelDescriptor::Bcast(
        x_coords.data(), nx_terrain, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(
        y_coords.data(), ny_terrain, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(
        z_values.data(), nx_terrain * ny_terrain, ParallelDescriptor::IOProcessorNumber());

    Gpu::DeviceVector<Real> x_device(nx_terrain);
    Gpu::DeviceVector<Real> y_device(ny_terrain);
    Gpu::DeviceVector<Real> z_device(nx_terrain * ny_terrain);

    Gpu::copy(Gpu::hostToDevice, x_coords.begin(), x_coords.end(), x_device.begin());
    Gpu::copy(Gpu::hostToDevice, y_coords.begin(), y_coords.end(), y_device.begin());
    Gpu::copy(Gpu::hostToDevice, z_values.begin(), z_values.end(), z_device.begin());

    Real* x_ptr = x_device.data();
    Real* y_ptr = y_device.data();
    Real* z_ptr = z_device.data();

    Real prob_lo_x = dg.geom.ProbLo(0);
    Real prob_lo_y = dg.geom.ProbLo(1);
    auto dx_dy = dg.geom.CellSizeArray();
    Real dx_d = dx_dy[0];
    Real dy_d = dx_dy[1];

    const Box& domain_dust = dg.geom.Domain();
    int i_dust_lo = domain_dust.smallEnd(0);
    int i_dust_hi = domain_dust.bigEnd(0);
    int j_dust_lo = domain_dust.smallEnd(1);
    int j_dust_hi = domain_dust.bigEnd(1);

    Real x_min = x_coords.front();
    Real x_max = x_coords.back();
    Real y_min = y_coords.front();
    Real y_max = y_coords.back();

    for (MFIter mfi(z_dust_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> z_dust = z_dust_nd.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            int i_d = iv[0];
            int j_d = iv[1];

            int i_clamped = amrex::max(i_dust_lo, amrex::min(i_dust_hi, i_d));
            int j_clamped = amrex::max(j_dust_lo, amrex::min(j_dust_hi, j_d));

            Real x = prob_lo_x + i_clamped * dx_d;
            Real y = prob_lo_y + j_clamped * dy_d;
            x = amrex::max(x_min, amrex::min(x_max, x));
            y = amrex::max(y_min, amrex::min(y_max, y));

            int ix_lo = nx_terrain - 2;
            int iy_lo = ny_terrain - 2;
            for (int ix = 0; ix < nx_terrain - 1; ++ix) {
                if (x_ptr[ix] <= x && x <= x_ptr[ix + 1]) {
                    ix_lo = ix;
                    break;
                }
            }
            for (int iy = 0; iy < ny_terrain - 1; ++iy) {
                if (y_ptr[iy] <= y && y <= y_ptr[iy + 1]) {
                    iy_lo = iy;
                    break;
                }
            }

            int ix_hi = amrex::min(ix_lo + 1, nx_terrain - 1);
            int iy_hi = amrex::min(iy_lo + 1, ny_terrain - 1);
            Real z_ll = z_ptr[ix_lo * ny_terrain + iy_lo];
            Real z_lr = z_ptr[ix_hi * ny_terrain + iy_lo];
            Real z_ul = z_ptr[ix_lo * ny_terrain + iy_hi];
            Real z_ur = z_ptr[ix_hi * ny_terrain + iy_hi];

            Real fx = (x_ptr[ix_hi] > x_ptr[ix_lo])
                        ? (x - x_ptr[ix_lo]) / (x_ptr[ix_hi] - x_ptr[ix_lo])
                        : 0.0_rt;
            Real fy = (y_ptr[iy_hi] > y_ptr[iy_lo])
                        ? (y - y_ptr[iy_lo]) / (y_ptr[iy_hi] - y_ptr[iy_lo])
                        : 0.0_rt;
            fx = amrex::max(0.0_rt, amrex::min(1.0_rt, fx));
            fy = amrex::max(0.0_rt, amrex::min(1.0_rt, fy));

            Real z_lo = z_ll * (1.0_rt - fx) + z_lr * fx;
            Real z_hi = z_ul * (1.0_rt - fx) + z_ur * fx;
            z_dust(iv, 0) = z_lo * (1.0_rt - fy) + z_hi * fy;
        });
    }

    // ParallelFor is asynchronous; the Gpu::DeviceVectors above free their device
    // allocations when this function returns. Let the kernels finish reading them.
    Gpu::streamSynchronize();

    return true;
}

} // namespace

void compute_dust_terrain_slopes(
    MultiFab& dust_slopes,
    const MultiFab& z_phys_nd,
    const Geometry& geom_atm,
    const DustGrid& dg,
    const std::string& terrain_fname)
{
    // Compute dz/dx and dz/dy on dust grid by interpolating from nodal z_phys_nd
    // or from fine-grid terrain file (if available)
    // Dust grid is C times finer than atmospheric grid

    int C = dg.grid_ratio;
    Real dx_atm = geom_atm.CellSize(0);
    Real dy_atm = geom_atm.CellSize(1);
    Real dx_dust = dg.geom.CellSize(0);
    Real dy_dust = dg.geom.CellSize(1);

    // Atmospheric domain
    const Box& domain_atm = geom_atm.Domain();
    int atm_nlo_x = domain_atm.smallEnd(0);
    int atm_nlo_y = domain_atm.smallEnd(1);

    // Try to use fine-grid terrain if file is available
    bool use_fine_terrain = false;
    std::unique_ptr<MultiFab> z_dust_nd_ptr;

    if (!terrain_fname.empty()) {
        z_dust_nd_ptr = std::make_unique<MultiFab>(dg.ba, dg.dm, 1, 1);
        if (read_terrain_onto_dust_grid(*z_dust_nd_ptr, dg, terrain_fname)) {
            use_fine_terrain = true;
            z_dust_nd_ptr->FillBoundary(dg.geom.periodicity());
            amrex::Print() << "[DUST] Using fine-grid terrain from file: "
                           << terrain_fname << std::endl;
        }
    }

    if (use_fine_terrain) {
        // Compute slopes from fine-grid terrain using centered bilinear averages
        // slopes(i,j,0,0) = 0.5 * ((z(i+1,j,0) - z(i,j,0)) + (z(i+1,j+1,0) - z(i,j+1,0))) / dx_f
        // slopes(i,j,0,1) = 0.5 * ((z(i,j+1,0) - z(i,j,0)) + (z(i+1,j+1,0) - z(i+1,j,0))) / dy_f

        for (MFIter mfi(dust_slopes, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            Array4<Real> slopes = dust_slopes.array(mfi);
            Array4<const Real> z_nd = z_dust_nd_ptr->array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
                int i = iv[0];
                int j = iv[1];

                // Centered bilinear average formula
                Real z_ll = z_nd(i,   j,   0);  // (i,   j)
                Real z_lr = z_nd(i+1, j,   0);  // (i+1, j)
                Real z_ul = z_nd(i,   j+1, 0);  // (i,   j+1)
                Real z_ur = z_nd(i+1, j+1, 0);  // (i+1, j+1)

                // dz/dx: average of differences at j and j+1
                slopes(i, j, 0, 0) = 0.5 * ((z_lr - z_ll) + (z_ur - z_ul)) / dx_dust;

                // dz/dy: average of differences at i and i+1
                slopes(i, j, 0, 1) = 0.5 * ((z_ul - z_ll) + (z_ur - z_lr)) / dy_dust;
            });
        }
    } else {
        // Fallback: use atmospheric terrain with fixed centered bilinear average formula
        // This fixes the dead-code bug in the original implementation

        for (MFIter mfi(dust_slopes, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            Array4<Real> slopes = dust_slopes.array(mfi);
            Array4<const Real> z_nd = z_phys_nd.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
                int i_d = iv[0];
                int j_d = iv[1];

                // Map fire cell to atmospheric cell indices
                // Fire cell (i_f, j_f) is in atmospheric cell (i_a, j_a)
                int i_a = i_d / C;
                int j_a = j_d / C;

                // Nodal indices: z_phys_nd is nodal (one extra in each direction)
                // Node (i_n, j_n) is at corner of cell (i_n-1, j_n-1)
                int i_n_lo = i_a + atm_nlo_x;
                int j_n_lo = j_a + atm_nlo_y;

                // Get the four surrounding nodal values for bilinear interp
                Real z_ll = z_nd(i_n_lo,   j_n_lo,   0);  // (i,   j)
                Real z_lr = z_nd(i_n_lo+1, j_n_lo,   0);  // (i+1, j)
                Real z_ul = z_nd(i_n_lo,   j_n_lo+1, 0);  // (i,   j+1)
                Real z_ur = z_nd(i_n_lo+1, j_n_lo+1, 0);  // (i+1, j+1)

                // Use centered bilinear average formula (fixes dead-code bug)
                // dz/dx: average of differences at j and j+1
                slopes(i_d, j_d, 0, 0) = 0.5 * ((z_lr - z_ll) + (z_ur - z_ul)) / dx_atm;

                // dz/dy: average of differences at i and i+1
                slopes(i_d, j_d, 0, 1) = 0.5 * ((z_ul - z_ll) + (z_ur - z_lr)) / dy_atm;
            });
        }
    }
}
