#include <ERF_TerrainSlope.H>
#include <ERF_FireTerrainReader.H>


using namespace amrex;

void compute_terrain_slopes(
    MultiFab& fire_slopes,
    const MultiFab& z_phys_nd,
    const Geometry& geom_atm,
    const FireGrid& fg,
    const std::string& terrain_fname)
{
    // Compute dz/dx and dz/dy on fire grid by interpolating from nodal z_phys_nd
    // or from fine-grid terrain file (if available)
    // Fire grid is C times finer than atmospheric grid

    int C = fg.C;
    Real dx_atm = geom_atm.CellSize(0);
    Real dy_atm = geom_atm.CellSize(1);
    Real dx_fire = fg.geom.CellSize(0);
    Real dy_fire = fg.geom.CellSize(1);

    // Atmospheric domain
    const Box& domain_atm = geom_atm.Domain();
    int atm_nlo_x = domain_atm.smallEnd(0);
    int atm_nlo_y = domain_atm.smallEnd(1);

    // Try to use fine-grid terrain if file is available
    bool use_fine_terrain = false;
    std::unique_ptr<MultiFab> z_fire_nd_ptr;

    if (!terrain_fname.empty()) {
        z_fire_nd_ptr = std::make_unique<MultiFab>(fg.ba, fg.dm, 1, 1);
        if (read_terrain_onto_fire_grid(*z_fire_nd_ptr, fg, terrain_fname)) {
            use_fine_terrain = true;
            z_fire_nd_ptr->FillBoundary(fg.geom.periodicity());
            amrex::Print() << "[FIRE] Using fine-grid terrain from file: " << terrain_fname << std::endl;
        }
    }

    if (use_fine_terrain) {
        // Compute slopes from fine-grid terrain using centered bilinear averages
        // slopes(i,j,0,0) = 0.5 * ((z(i+1,j,0) - z(i,j,0)) + (z(i+1,j+1,0) - z(i,j+1,0))) / dx_f
        // slopes(i,j,0,1) = 0.5 * ((z(i,j+1,0) - z(i,j,0)) + (z(i+1,j+1,0) - z(i+1,j,0))) / dy_f

        for (MFIter mfi(fire_slopes, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            Array4<Real> slopes = fire_slopes.array(mfi);
            Array4<const Real> z_nd = z_fire_nd_ptr->array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
                int i = iv[0];
                int j = iv[1];

                // Centered bilinear average formula
                Real z_ll = z_nd(i,   j,   0);  // (i,   j)
                Real z_lr = z_nd(i+1, j,   0);  // (i+1, j)
                Real z_ul = z_nd(i,   j+1, 0);  // (i,   j+1)
                Real z_ur = z_nd(i+1, j+1, 0);  // (i+1, j+1)

                // dz/dx: average of differences at j and j+1
                slopes(i, j, 0, 0) = 0.5 * ((z_lr - z_ll) + (z_ur - z_ul)) / dx_fire;

                // dz/dy: average of differences at i and i+1
                slopes(i, j, 0, 1) = 0.5 * ((z_ul - z_ll) + (z_ur - z_lr)) / dy_fire;
            });
        }
    } else {
        // Fallback: use atmospheric terrain with fixed centered bilinear average formula
        // This fixes the dead-code bug in the original implementation

        for (MFIter mfi(fire_slopes, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            Array4<Real> slopes = fire_slopes.array(mfi);
            Array4<const Real> z_nd = z_phys_nd.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
                int i_f = iv[0];
                int j_f = iv[1];

                // Map fire cell to atmospheric cell indices
                // Fire cell (i_f, j_f) is in atmospheric cell (i_a, j_a)
                int i_a = i_f / C;
                int j_a = j_f / C;

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
                slopes(i_f, j_f, 0, 0) = 0.5 * ((z_lr - z_ll) + (z_ur - z_ul)) / dx_atm;

                // dz/dy: average of differences at i and i+1
                slopes(i_f, j_f, 0, 1) = 0.5 * ((z_ul - z_ll) + (z_ur - z_lr)) / dy_atm;
            });
        }
    }
}
